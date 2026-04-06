#!/usr/bin/env python3
"""Tabulate Prokka/RGI annotations into pangenome feature tables."""

from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import multiprocessing as mp
import os
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable


import pandas as pd
from BCBio import GFF


RGI_MUT_SPLIT_RE = re.compile(r"\s*,\s*")
RGI_MUT_RE = re.compile(r"^([A-Za-z\-*])(\d+)([A-Za-z\-*])$")


@dataclass(frozen=True)
class RawAnnotation:
    """Container for a single annotation (RGI or Prokka) from one accession.

    Parameters
    ----------
    accession : str
        Genome accession this annotation belongs to.
    source : str
        Source of the annotation, either "RGI" or "PROKKA".
    contig : str
        Contig identifier where the feature is located.
    start, end : int
        1-based inclusive coordinates of the feature on the contig.
    strand : str
        Strand orientation: '+' or '-'.
    is_rna : bool
        True when the feature is an RNA (tRNA/rRNA/etc.).
    is_plasmid : bool
        True when the contig was predicted as plasmid by Platon.
    function_id : str
        Function identifier (RGI uses A#######, Prokka uses P#######).
    locus_tag : str or None
        Locus tag (from Prokka) when available.
    sequence : str or None
        Sequence associated with the feature (AA for proteins, NT for RNA).
    mutations : tuple[str, ...]
        Normalized mutation strings associated with the feature.
    """

    accession: str
    source: str  # "RGI" | "PROKKA"
    contig: str
    start: int
    end: int
    strand: str  # "+" | "-"
    is_rna: bool
    is_plasmid: bool
    function_id: str  # A####### or P#######
    locus_tag: str | None = None
    sequence: str | None = None
    mutations: tuple[str, ...] = ()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Namespace with attributes:
        - INDIR : input directory containing per-accession annotation folders
        - OUTPREFIX : output prefix for CSV, JSON and cache directory
    """
    parser = argparse.ArgumentParser(
        description="Tabulate genome annotations into feature counts."
    )
    parser.add_argument("INDIR", help="Input directory with per-accession annotations")
    parser.add_argument("OUTPREFIX", help="Output prefix for CSV/JSON/cache")
    return parser.parse_args()


def get_jobs() -> int:
    """Determine number of parallel jobs from environment.

    Reads the JOBS environment variable and returns it as a positive integer.
    If JOBS is unset or cannot be parsed as a positive integer, returns the
    default value 4.

    Returns
    -------
    int
        Number of parallel jobs to use.
    """
    try:
        jobs = int(os.environ.get("JOBS", ""))
        if jobs > 0:
            return jobs
    except ValueError:
        pass
    return 4


def parse_plasmids(platon_tsv: Path) -> set[str]:
    """Parse Platon TSV to get plasmid contig IDs.

    Parameters
    ----------
    platon_tsv : pathlib.Path
        Path to the Platon output TSV file. If the file does not exist or is
        empty, an empty set is returned.

    Returns
    -------
    set[str]
        Set of contig IDs predicted as plasmids (taken from the 'ID' column).
    """
    if not platon_tsv.exists() or platon_tsv.stat().st_size == 0:
        return set()
    with platon_tsv.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["ID"] for row in reader if row.get("ID")}


def parse_gtdb_info(sourmash_csv: Path) -> tuple[str, float]:
    """Read Sourmash CSV and extract GTDB species and ANI.

    The function reads the first non-header row of the sourmash CSV and parses
    the `name` column to extract a two-token species string following the
    accession. If the first token after the accession is 'uncultured' it is
    removed before selecting the two tokens. The `ani` column is parsed as a
    float and may be NaN when empty.

    Parameters
    ----------
    sourmash_csv : pathlib.Path
        Path to the sourmash CSV file produced by the pipeline.

    Returns
    -------
    tuple[str, float]
        (gtdb_species, gtdb_ani). gtdb_ani may be math.nan if not present.
    """
    with sourmash_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        row = next(reader, None)
        if row is None:
            raise ValueError(f"{sourmash_csv} has no rows")
        name_tokens = row.get("name", "").split()
        if len(name_tokens) <= 1:
            raise ValueError(
                f"Could not parse species from sourmash name: {row.get('name')}"
            )
        tokens = name_tokens[1:]
        if tokens and tokens[0].lower() == "uncultured":
            tokens = tokens[1:]
        if len(tokens) < 2:
            raise ValueError(
                f"Could not derive species from sourmash name: {row.get('name')}"
            )
        species = f"{tokens[0]} {tokens[1]}"
        ani_raw = row.get("ani", "")
        ani = float(ani_raw) if ani_raw not in {"", None} else math.nan
        return species, ani


def parse_rgi_mutations(*mut_fields: str) -> tuple[str, ...]:
    """Parse RGI mutation fields into normalized mutation tokens.

    RGI reports mutations in two columns which may contain comma-separated
    values and may include trailing ":<number>" annotations. This function
    splits the provided fields on commas, strips any trailing ":..." suffix,
    normalizes tokens (via :func:`normalize_mutation_token`) and returns a
    tuple of normalized mutation strings.

    Parameters
    ----------
    *mut_fields : str
        One or more raw mutation strings from RGI columns (e.g.,
        'S83I:3329, D87G').

    Returns
    -------
    tuple[str, ...]
        Tuple of normalized mutation strings, zero-length if none were found.
    """
    parsed: list[str] = []
    for field in mut_fields:
        if not field or field in {"n/a", "NA"}:
            continue
        for entry in RGI_MUT_SPLIT_RE.split(field.strip()):
            token = entry.split(":", 1)[0].strip()
            if not token:
                continue
            normalized = normalize_mutation_token(token)
            if normalized is not None:
                parsed.append(normalized)
    return tuple(parsed)


def normalize_mutation_token(token: str) -> str | None:
    """Normalize a single mutation token.

    Accepts simple amino-acid or nucleotide mutation tokens of the form
    'S83I', 'G442E', or with '-' for indels. Returns a canonical string where
    the position is zero-padded to seven digits (e.g., 'S0000083I'). If the
    token does not match the expected pattern, returns None.

    Parameters
    ----------
    token : str
        Raw mutation token to normalize.

    Returns
    -------
    str or None
        Normalized mutation string or None if the token could not be parsed.
    """
    match = RGI_MUT_RE.match(token)
    if not match:
        return None
    ref, pos, alt = match.groups()
    return f"{ref.upper()}{int(pos):07d}{alt.upper()}"


def normalize_aro(aro_value: str) -> str:
    """Normalize an ARO identifier to a seven-digit numeric string.

    Parameters
    ----------
    aro_value : str
        Raw ARO value (may contain non-digit characters). The function
        extracts the trailing digits, trims to seven digits if longer, and
        left-pads with zeros.

    Returns
    -------
    str
        Seven-digit zero-padded ARO identifier.

    Raises
    ------
    ValueError
        If no digits can be extracted from the provided value.
    """
    digits = re.sub(r"\D", "", aro_value or "")
    if not digits:
        raise ValueError(f"Invalid ARO value: {aro_value!r}")
    if len(digits) > 7:
        digits = digits[-7:]
    return digits.zfill(7)


def extract_feature_sequence(record_seq, feature, is_rna: bool) -> str:
    """Retrieve the sequence for a feature from a GFF record.

    For CDS features, prefer the 'translation' qualifier when present (Prokka
    often provides translated peptide sequences). For RNA features, return the
    nucleotide sequence. If no sequence can be obtained, return an empty
    string.

    Parameters
    ----------
    record_seq : Bio.Seq.Seq or similar
        Sequence of the parent record (used to extract nucleotide subsequences).
    feature : BCBio GFF feature
        The feature object containing location and qualifiers.
    is_rna : bool
        If True, treat the feature as RNA and return nucleotide sequence;
        otherwise attempt to return an amino-acid sequence (translation).

    Returns
    -------
    str
        Sequence string (AA for proteins, NT for RNA) or empty string if not
        available.
    """
    if not is_rna:
        translation = feature.qualifiers.get("translation", [])
        if translation:
            return str(translation[0]).strip()
    if record_seq is None or len(record_seq) == 0:
        return ""
    seq = feature.extract(record_seq)
    if is_rna:
        return str(seq).upper()
    return str(seq.translate(to_stop=False)).rstrip("*").upper()


def is_rna_feature_type(feature_type: str) -> bool:
    """Heuristic to decide whether a GFF feature type represents RNA.

    Parameters
    ----------
    feature_type : str
        GFF feature.type string (e.g., 'CDS', 'tRNA').

    Returns
    -------
    bool
        True when the type is considered RNA (any type containing 'rna' but
        not exactly 'CDS').
    """
    t = (feature_type or "").lower()
    return t != "cds" and "rna" in t


def parse_accession_raw(accession_dir: Path) -> dict:
    """Parse all annotation files for one accession directory.

    This function reads Prokka GFF, RGI TSV, sourmash CSV and optional Platon
    TSV for a single accession directory and returns a dict containing:
    - accession
    - gtdb_species
    - gtdb_ani
    - rgi_rows : list of parsed RGI annotation dicts
    - prokka_rows : list of parsed Prokka annotation dicts

    Parameters
    ----------
    accession_dir : pathlib.Path
        Directory containing files named <accession>.* produced by the
        per-accession annotation step.

    Returns
    -------
    dict
        Summary dictionary described above.
    """
    accession = accession_dir.name
    prokka_gff = accession_dir / f"{accession}.prokka.gff"
    rgi_txt = accession_dir / f"{accession}.rgi.txt"
    sourmash_csv = accession_dir / f"{accession}.sourmash.csv"
    platon_tsv = accession_dir / f"{accession}.platon.tsv"

    plasmids = parse_plasmids(platon_tsv)
    gtdb_species, gtdb_ani = parse_gtdb_info(sourmash_csv)

    rgi_df = pd.read_csv(rgi_txt, sep="\t", dtype=str, keep_default_na=False)
    rgi_df = rgi_df[rgi_df["Cut_Off"] != "Loose"]
    rgi_rows = []
    for row in rgi_df.to_dict(orient="records"):  # pyright: ignore[reportCallIssue]
        is_rna = not row.get("Model_type", "").lower().startswith("protein")
        start = int(row["Start"])
        end = int(row["Stop"])
        if start > end:
            start, end = end, start
        muts = parse_rgi_mutations(
            row.get("SNPs_in_Best_Hit_ARO", ""), row.get("Other_SNPs", "")
        )
        rgi_rows.append(
            {
                "contig": row["Contig"],
                "start": start,
                "end": end,
                "strand": row["Orientation"],
                "is_rna": is_rna,
                "is_plasmid": row["Contig"] in plasmids,
                "aro": normalize_aro(row.get("ARO", "")),
                "mutations": muts,
            }
        )

    prokka_rows: list[dict] = []
    with prokka_gff.open("r", encoding="utf-8") as handle:
        for record in GFF.parse(handle):
            for feature in record.features:
                if feature.type != "CDS" and not is_rna_feature_type(feature.type):
                    continue
                start = int(feature.location.start) + 1
                end = int(feature.location.end)
                if start > end:
                    start, end = end, start
                strand = "+" if int(feature.location.strand or 1) >= 0 else "-"
                is_rna = feature.type != "CDS"
                qualifiers = feature.qualifiers
                product = str(
                    qualifiers.get("product", ["hypothetical protein"])[0]
                ).strip()
                locus_tag = str(
                    qualifiers.get("locus_tag", [qualifiers.get("ID", [""])[0]])[0]
                )
                seq = extract_feature_sequence(record.seq, feature, is_rna)
                if not seq:
                    continue
                prokka_rows.append(
                    {
                        "contig": str(record.id),
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "is_rna": is_rna,
                        "is_plasmid": str(record.id) in plasmids,
                        "product": product,
                        "locus_tag": locus_tag,
                        "sequence": seq,
                    }
                )

    return {
        "accession": accession,
        "gtdb_species": gtdb_species,
        "gtdb_ani": gtdb_ani,
        "rgi_rows": rgi_rows,
        "prokka_rows": prokka_rows,
    }


def load_cached_accession_raw(cache_path: Path, accession: str) -> dict | None:
    """Load cached parse result for one accession.

    Parameters
    ----------
    cache_path : pathlib.Path
        Path to the per-accession cache JSON file.
    accession : str
        Expected accession identifier.

    Returns
    -------
    dict or None
        Cached parse result if valid, otherwise None.
    """
    if not cache_path.exists():
        return None
    with cache_path.open("r", encoding="utf-8") as handle:
        cached = json.load(handle)
    if not isinstance(cached, dict):
        return None
    if cached.get("accession") != accession:
        return None
    if "rgi_rows" not in cached or "prokka_rows" not in cached:
        return None
    if not isinstance(cached["rgi_rows"], list) or not isinstance(
        cached["prokka_rows"], list
    ):
        return None
    return cached


def write_cached_accession_raw(cache_path: Path, result: dict) -> None:
    """Persist parsed accession result in cache JSON.

    Parameters
    ----------
    cache_path : pathlib.Path
        Destination cache JSON path.
    result : dict
        Parsed accession result produced by :func:`parse_accession_raw`.
    """
    with cache_path.open("w", encoding="utf-8") as handle:
        json.dump(result, handle, separators=(",", ":"))
        handle.write("\n")


def overlap_ratio_union(a_start: int, a_end: int, b_start: int, b_end: int) -> float:
    """Compute overlap divided by union span between two intervals.

    The function computes the length of the intersection divided by the span of
    the union, i.e. overlap_len / (max(end)-min(start)+1). This is the
    metric used to decide whether a Prokka and RGI annotation overlap >=95%.

    Parameters
    ----------
    a_start, a_end, b_start, b_end : int
        1-based inclusive coordinates for the two intervals.

    Returns
    -------
    float
        Fraction in range [0.0, 1.0] describing overlap/span.
    """
    overlap = max(0, min(a_end, b_end) - max(a_start, b_start) + 1)
    if overlap == 0:
        return 0.0
    span = max(a_end, b_end) - min(a_start, b_start) + 1
    return overlap / span


def load_existing_prokka_mapping(json_path: Path) -> dict[str, str]:
    """Load existing Prokka name->ID mapping from JSON if present.

    The JSON produced by this script maps Prokka IDs (e.g., 'P0000001') to
    product names. This helper loads that JSON and returns a mapping from
    product name to numeric ID string (without the 'P' prefix) for reuse when
    assigning stable IDs across reruns.

    Parameters
    ----------
    json_path : pathlib.Path
        Path to the existing OUTPREFIX.json file.

    Returns
    -------
    dict[str, str]
        Mapping product name -> seven-digit numeric ID string (e.g., 'beta'
        -> '0000001'). Empty dict if json_path does not exist or is invalid.
    """
    if not json_path.exists():
        return {}
    with json_path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        return {}
    mapping = {}
    for k, v in data.items():
        if isinstance(k, str) and isinstance(v, str) and re.match(r"^P\d{7}$", k):
            mapping[v] = k[1:]
    return mapping


def build_prokka_id_map(
    all_products: Iterable[str], existing_name_to_id: dict[str, str]
) -> tuple[dict[str, str], dict[str, str]]:
    """Create a stable mapping from Prokka product names to P####### IDs.

    This function merges existing mappings with new product names discovered
    in the current run. It reserves '0000000' for 'hypothetical protein' and
    assigns incrementing numeric IDs for previously unseen products. It
    returns both name->id and id->product mappings suitable for writing the
    OUTPREFIX.json file.

    Parameters
    ----------
    all_products : iterable of str
        All Prokka product strings observed across accessions.
    existing_name_to_id : dict
        Existing mapping from product name to numeric id string (no 'P' prefix),
        if any (loaded from previous OUTPREFIX.json).

    Returns
    -------
    tuple[dict[str,str], dict[str,str]]
        (name_to_id, id_to_name) where name_to_id maps product->numeric-id
        (e.g., 'beta'->'0000001') and id_to_name maps 'P{numeric}'->product.
    """
    name_to_id = dict(existing_name_to_id)
    name_to_id["hypothetical protein"] = "0000000"
    used = {int(v) for v in name_to_id.values() if v.isdigit()}
    next_id = max([0, *used]) + 1

    for product in sorted(set(all_products)):
        if product == "hypothetical protein":
            continue
        if product in name_to_id:
            continue
        while next_id in used or next_id == 0:
            next_id += 1
        name_to_id[product] = f"{next_id:07d}"
        used.add(next_id)
        next_id += 1

    id_to_name = {f"P{pid}": name for name, pid in name_to_id.items()}
    id_to_name = dict(sorted(id_to_name.items()))
    return name_to_id, id_to_name


def write_source_fastas(
    merged_by_accession: dict[str, list[RawAnnotation]], cache_dir: Path
) -> set[tuple[str, bool]]:
    """Write per-product FASTA files for Prokka annotations.

    For each unique Prokka product identifier (numeric PID) the function
    writes a FASTA file named '<PID>.faa' or '<PID>.fna' into the cache
    directory. Sequence headers are formed as '<accession>__<locus_tag>'. To
    reduce race conditions this function writes files atomically by comparing
    content and only rewriting when necessary.

    Parameters
    ----------
    merged_by_accession : dict
        Mapping accession -> list of RawAnnotation objects (merged RGI+Prokka).
    cache_dir : pathlib.Path
        Directory where per-product FASTA files are stored.

    Returns
    -------
    set[tuple[str,bool]]
        Set of (product_id_without_prefix, is_rna) tuples written to FASTA.
        This is used by later clustering steps.
    """
    by_product_type: dict[tuple[str, bool], list[tuple[str, str]]] = defaultdict(list)
    for anns in merged_by_accession.values():
        for ann in anns:
            if ann.source != "PROKKA":
                continue
            if ann.locus_tag is None or ann.sequence is None:
                raise RuntimeError("Missing locus_tag/sequence for Prokka annotation")
            seq_id = f"{ann.accession}__{ann.locus_tag}"
            by_product_type[(ann.function_id[1:], ann.is_rna)].append(
                (seq_id, ann.sequence)
            )

    for (pid, is_rna), records in by_product_type.items():
        ext = "fna" if is_rna else "faa"
        out_path = cache_dir / f"{pid}.{ext}"
        lines = [f">{seq_id}\n{seq}\n" for seq_id, seq in records]
        new_content = "".join(lines)
        if out_path.exists() and out_path.read_text(encoding="utf-8") == new_content:
            continue
        out_path.write_text(new_content, encoding="utf-8")
    return set(by_product_type.keys())


def run_cmd(
    cmd: list[str], cwd: Path | None = None, stdout_path: Path | None = None
) -> None:
    """Run a subprocess command with optional stdout redirection.

    This wrapper calls subprocess.run with check=True and optionally writes the
    command stdout to the provided file path. Any exception from subprocess is
    propagated to the caller.

    Parameters
    ----------
    cmd : list[str]
        Command and arguments to execute.
    cwd : pathlib.Path or None
        Working directory for the command.
    stdout_path : pathlib.Path or None
        If provided, open this file and redirect the command's stdout into it.
    """
    if stdout_path is None:
        subprocess.run(cmd, check=True, cwd=cwd)
        return
    with stdout_path.open("w", encoding="utf-8") as out_handle:
        subprocess.run(cmd, check=True, cwd=cwd, stdout=out_handle)


def parse_mmseqs_cluster_tsv(cluster_tsv_path: Path) -> list[list[str]]:
    """Parse an mmseqs ``*_cluster.tsv`` into ordered member clusters.

    The cluster TSV emitted by mmseqs uses two tab-delimited columns:
    representative ID and member ID. This helper groups rows by representative
    while preserving first-seen order of both representatives and members.

    Parameters
    ----------
    cluster_tsv_path : pathlib.Path
        Path to the mmseqs ``*_cluster.tsv`` file.

    Returns
    -------
    list[list[str]]
        Ordered clusters where each item is a list of member sequence IDs.
    """
    clusters: dict[str, list[str]] = {}
    with cluster_tsv_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            rep, member = fields[0].strip(), fields[1].strip()
            if not rep or not member:
                continue
            if rep not in clusters:
                clusters[rep] = []
            if member not in clusters[rep]:
                clusters[rep].append(member)
    return list(clusters.values())


def read_fasta_map(fasta_path: Path) -> dict[str, str]:
    """Read FASTA records into a header-to-sequence mapping.

    When a header appears multiple times, a non-empty sequence replaces an
    existing empty sequence. This behavior handles mmseqs flat MSA files where
    representative headers can appear once as a marker and again with aligned
    sequence content.

    Parameters
    ----------
    fasta_path : pathlib.Path
        Path to the FASTA file.

    Returns
    -------
    dict[str, str]
        Mapping from FASTA header (without '>') to sequence string.
    """
    records: dict[str, str] = {}
    header: str | None = None
    seq_lines: list[str] = []
    with fasta_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(seq_lines)
                    if header not in records or (not records[header] and seq):
                        records[header] = seq
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header is not None:
            seq = "".join(seq_lines)
            if header not in records or (not records[header] and seq):
                records[header] = seq
    return records


def split_mmseqs_msa(
    msa_path: Path, out_prefix: str, cluster_tsv_path: Path, source_fasta: Path
) -> list[Path]:
    """Write per-cluster MSA FASTAs using mmseqs cluster membership.

    Cluster membership is sourced from ``*_cluster.tsv`` because inferring
    boundaries from flat MSA text can drop members when mmseqs emits header-only
    records. Aligned sequences are taken from ``msa_path`` when present and
    fall back to the original source FASTA otherwise.

    Parameters
    ----------
    msa_path : pathlib.Path
        Path to the mmseqs-produced flat MSA.
    out_prefix : str
        Prefix used for per-cluster output filenames.
    cluster_tsv_path : pathlib.Path
        Path to ``*_cluster.tsv`` for the same product.
    source_fasta : pathlib.Path
        Source FASTA used for clustering.

    Returns
    -------
    list[pathlib.Path]
        List of written per-cluster MSA FASTA paths.
    """
    ext = "fna" if msa_path.name.endswith(".fna") else "faa"
    old_files = glob.glob(f"{out_prefix}.*.msa.{ext}")
    for file_path in old_files:
        os.unlink(file_path)

    msa_records = read_fasta_map(msa_path) if msa_path.exists() else {}
    source_records = read_fasta_map(source_fasta)
    if cluster_tsv_path.exists():
        clusters = parse_mmseqs_cluster_tsv(cluster_tsv_path)
    else:
        clusters = [[header] for header in source_records]

    out_paths: list[Path] = []
    for idx, members in enumerate(clusters):
        out_path = Path(f"{out_prefix}.{idx:07d}.msa.{ext}")
        with out_path.open("w", encoding="utf-8") as handle:
            for member in members:
                seq = msa_records.get(member) or source_records.get(member)
                if seq is None:
                    raise RuntimeError(
                        f"Missing sequence '{member}' while splitting {cluster_tsv_path.name}"
                    )
                handle.write(f">{member}\n{seq}\n")
        out_paths.append(out_path)
    return out_paths


def msa_has_accession_headers(msa_path: Path) -> bool:
    """Detect whether an MSA file uses accession__locus headers.

    Some older mmseqs result files use numeric identifiers rather than the
    fasta headers embedded in the original sequences. This helper scans the
    first FASTA header and returns True when it detects the '{{accession}}__'
    pattern used by this pipeline.

    Parameters
    ----------
    msa_path : pathlib.Path
        Path to the MSA file to inspect.

    Returns
    -------
    bool
        True when accession__locus headers are present, False otherwise.
    """
    if not msa_path.exists():
        return False
    with msa_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            if "__" in header:
                return True
            # First header is enough to detect old numeric-id files.
            return False
    return False


def fasta_record_count(fasta_path: Path) -> int:
    """Count FASTA records in a file.

    Parameters
    ----------
    fasta_path : pathlib.Path
        Path to the FASTA file.

    Returns
    -------
    int
        Number of records (lines starting with '>') in the file.
    """
    count = 0
    with fasta_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def run_product_pipeline(task: tuple[str, bool, str]) -> tuple[str, bool]:
    """Run clustering, MSA and SNP-calling pipeline for one Prokka product.

    This function executes the sequence of external tools required to cluster
    sequences for a given Prokka product, build per-cluster MSAs and call
    variants with snp-sites. It is intended to be run in parallel across
    products.

    Parameters
    ----------
    task : tuple
        (product_id_without_prefix, is_rna, cache_dir_str). product_id is the
        seven-digit identifier assigned to the Prokka product (without the
        leading 'P'), is_rna indicates nucleotide vs amino-acid mode, and
        cache_dir_str is the path to the OUTPREFIX.cache directory.

    Returns
    -------
    tuple[str, bool]
        The (product_id, is_rna) tuple for the completed task.
    """
    pid, is_rna, cache_dir_str = task
    cache_dir = Path(cache_dir_str)
    ext = "fna" if is_rna else "faa"
    dbtype = "2" if is_rna else "1"
    min_seq = "0.95" if is_rna else "0.9"

    source_fasta = cache_dir / f"{pid}.{ext}"
    prefix = cache_dir / pid
    tmp_dir = cache_dir / f"{pid}_tmp"
    flat_msa = cache_dir / f"{pid}.msa.{ext}"
    clu_msa = tmp_dir / "latest" / "clu_msa"
    tmp_input = tmp_dir / "latest" / "input"
    tmp_clu = tmp_dir / "latest" / "clu"
    clu_seqs = tmp_dir / "latest" / "clu_seqs"
    cluster_tsv = cache_dir / f"{pid}_cluster.tsv"
    seq_count = fasta_record_count(source_fasta)

    try:
        if seq_count <= 1:
            if (
                not flat_msa.exists()
            ) or source_fasta.stat().st_mtime > flat_msa.stat().st_mtime:
                shutil.copyfile(source_fasta, flat_msa)
        else:
            need_cluster = (
                (not flat_msa.exists())
                or source_fasta.stat().st_mtime > flat_msa.stat().st_mtime
                or (not msa_has_accession_headers(flat_msa))
            )
            if need_cluster:
                if tmp_dir.exists():
                    shutil.rmtree(tmp_dir)
                cmd = [
                    "mmseqs",
                    "easy-linclust",
                    "--threads",
                    "1",
                    "--min-seq-id",
                    min_seq,
                    "-c",
                    min_seq,
                    "--cov-mode",
                    "1",
                    "--cluster-mode",
                    "2",
                    "--remove-tmp-files",
                    "0",
                    "--dbtype",
                    dbtype,
                    str(source_fasta),
                    str(prefix),
                    str(tmp_dir),
                ]
                run_cmd(cmd)

                if is_rna:
                    run_cmd(
                        [
                            "mmseqs",
                            "apply",
                            str(clu_seqs),
                            str(clu_msa),
                            "--",
                            "mafft",
                            "--retree",
                            "1",
                            "--maxiterate",
                            "0",
                            "--nofft",
                            "--parttree",
                            "--nuc",
                            "-",
                        ]
                    )
                else:
                    run_cmd(
                        [
                            "mmseqs",
                            "result2msa",
                            str(tmp_input),
                            str(tmp_input),
                            str(tmp_clu),
                            str(clu_msa),
                        ]
                    )
                run_cmd(
                    [
                        "mmseqs",
                        "result2flat",
                        str(tmp_input),
                        str(tmp_input),
                        str(clu_msa),
                        str(flat_msa),
                        "--use-fasta-header",
                        "1",
                    ]
                )

        cluster_paths = split_mmseqs_msa(
            flat_msa, str(cache_dir / pid), cluster_tsv, source_fasta
        )
        old_vcfs = glob.glob(str(cache_dir / f"{pid}.*.vcf"))
        for file_path in old_vcfs:
            os.unlink(file_path)
        for cluster_path in cluster_paths:
            stem = cluster_path.name.replace(f".msa.{ext}", "")
            vcf_path = cache_dir / f"{stem}.vcf"
            with vcf_path.open("w", encoding="utf-8") as out_handle:
                proc = subprocess.run(
                    ["snp-sites", "-v", str(cluster_path)],
                    stdout=out_handle,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=False,
                )
            if proc.returncode != 0:
                stderr = proc.stderr or ""
                if "No SNPs were detected" not in stderr:
                    raise subprocess.CalledProcessError(
                        returncode=proc.returncode,
                        cmd=proc.args,
                        stderr=stderr,
                    )
        return pid, is_rna
    finally:
        if tmp_dir.exists():
            shutil.rmtree(tmp_dir)


def parse_vcf_mutations(vcf_path: Path) -> dict[str, list[str]]:
    """Parse a VCF file produced by snp-sites into per-sample mutation lists.

    Only simple SNPs are returned (single-base ref and alt). The sample
    identifiers are expected to be fasta headers produced by mmseqs 'result2flat'
    (e.g., '<accession>__<locus_tag>'). The return is a dict mapping sample
    header -> list of normalized mutation strings.

    Parameters
    ----------
    vcf_path : pathlib.Path
        Path to the VCF file to parse.

    Returns
    -------
    dict[str, list[str]]
        Mapping of sample header to list of mutation strings.
    """
    sample_mutations: dict[str, list[str]] = defaultdict(list)
    if not vcf_path.exists() or vcf_path.stat().st_size == 0:
        return sample_mutations

    header_samples: list[str] = []
    with vcf_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header_samples = line.rstrip("\n").split("\t")[9:]
                continue
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            pos = int(fields[1])
            ref = fields[3]
            alts = fields[4].split(",")
            sample_fields = fields[9:]
            for sample, gt_field in zip(header_samples, sample_fields):
                gt = gt_field.split(":")[0]
                if gt in {".", "./.", "0", "0/0", "0|0"}:
                    continue
                allele_idx = None
                for token in re.split(r"[\/|]", gt):
                    if token not in {"", "."}:
                        value = int(token)
                        if value > 0:
                            allele_idx = value
                            break
                if allele_idx is None or allele_idx > len(alts):
                    continue
                alt = alts[allele_idx - 1]
                if len(ref) != 1 or len(alt) != 1:
                    continue
                mut = f"{ref.upper()}{pos:07d}{alt.upper()}"
                sample_mutations[sample].append(mut)
    return sample_mutations


def collect_prokka_cluster_map_for_product(
    task: tuple[str, bool, str],
) -> tuple[dict[tuple[str, str], str], dict[tuple[str, str], list[str]]]:
    """Collect family and mutation maps for a single Prokka product."""
    pid, is_rna, cache_dir_str = task
    cache_dir = Path(cache_dir_str)
    seq_to_family: dict[tuple[str, str], str] = {}
    seq_to_mutations: dict[tuple[str, str], list[str]] = defaultdict(list)

    cluster_tsv = cache_dir / f"{pid}_cluster.tsv"
    if cluster_tsv.exists():
        for idx, members in enumerate(parse_mmseqs_cluster_tsv(cluster_tsv)):
            family = f"{'R' if is_rna else 'P'}{idx:07d}"
            for member in members:
                if "__" not in member:
                    continue
                accession, locus_tag = member.split("__", 1)
                seq_to_family[(accession, locus_tag)] = family

    vcf_paths = sorted(Path(p) for p in glob.glob(str(cache_dir / f"{pid}.*.vcf")))
    for vcf_path in vcf_paths:
        sample_muts = parse_vcf_mutations(vcf_path)
        for sample, muts in sample_muts.items():
            if "__" not in sample:
                continue
            accession, locus_tag = sample.split("__", 1)
            seq_to_mutations[(accession, locus_tag)].extend(muts)

    return seq_to_family, dict(seq_to_mutations)


def collect_prokka_cluster_maps(
    cache_dir: Path, product_types: Iterable[tuple[str, bool]]
) -> tuple[dict[tuple[str, str], str], dict[tuple[str, str], list[str]]]:
    """Collect mapping from Prokka sequences to cluster families and mutations.

    This function uses mmseqs ``*_cluster.tsv`` files to define family
    membership and reads per-cluster VCFs to collect mutations, building two
    structures:
      - seq_to_family: (accession, locus_tag) -> family id string (e.g., 'P0000001')
      - seq_to_mutations: (accession, locus_tag) -> list of mutation strings

    Parameters
    ----------
    cache_dir : pathlib.Path
        Path to the OUTPREFIX.cache directory where cluster MSAs and VCFs are
        stored.
    product_types : iterable
        Iterable of (product_id_without_prefix, is_rna) tuples used to seed the
        clustering.

    Returns
    -------
    tuple
        (seq_to_family, seq_to_mutations) as described above.
    """
    seq_to_family: dict[tuple[str, str], str] = {}
    seq_to_mutations: dict[tuple[str, str], list[str]] = defaultdict(list)

    try:
        n_products = len(product_types)  # pyright: ignore[reportArgumentType]
    except TypeError:
        n_products = sum(1 for _ in product_types)

    tasks = [(pid, is_rna, str(cache_dir)) for pid, is_rna in sorted(product_types)]
    if not tasks:
        return seq_to_family, seq_to_mutations

    jobs = get_jobs()
    with mp.Pool(processes=jobs) as pool:
        per_product_results = pool.map(collect_prokka_cluster_map_for_product, tasks)

    for family_map, mut_map in per_product_results:
        seq_to_family.update(family_map)
        for key, muts in mut_map.items():
            seq_to_mutations[key].extend(muts)
    return seq_to_family, seq_to_mutations


def build_outputs(
    merged_by_accession: dict[str, list[RawAnnotation]],
    meta_by_accession: dict[str, tuple[str, float]],
    seq_to_family: dict[tuple[str, str], str],
    seq_to_mutations: dict[tuple[str, str], list[str]],
) -> pd.DataFrame:
    """Aggregate feature counts and build final pandas DataFrame.

    For each accession, iterate through merged annotations (RGI and Prokka),
    compute feature column names according to the pipeline's naming scheme and
    increment counts for genes and specific mutations. RGI annotations take
    precedence and are represented with ARO-based function IDs; Prokka
    annotations use assigned P####### IDs and cluster family identifiers.

    Parameters
    ----------
    merged_by_accession : dict
        Mapping accession -> list of RawAnnotation objects (merged set).
    meta_by_accession : dict
        Mapping accession -> (gtdb_species, gtdb_ani).
    seq_to_family : dict
        Mapping (accession, locus_tag) -> family id string produced by
        clustering (e.g., 'P0000001').
    seq_to_mutations : dict
        Mapping (accession, locus_tag) -> list of mutation strings found in
        the VCFs.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns ['accession','gtdb_species','gtdb_ani', ...features...]
        where feature columns are integers counting occurrences and mutation
        occurrences.
    """
    rows: list[dict] = []
    all_features: set[str] = set()

    for accession in sorted(merged_by_accession):
        species, ani = meta_by_accession[accession]
        feature_counts: Counter[str] = Counter()
        for ann in merged_by_accession[accession]:
            plasmid_prefix = "P" if ann.is_plasmid else "N"
            if ann.source == "RGI":
                family = f"{'R' if ann.is_rna else 'P'}0000000"
                base = f"{plasmid_prefix}_{ann.function_id}_{family}"
                feature_counts[base] += 1
                for mut in ann.mutations:
                    feature_counts[f"{base}_{mut}"] += 1
                continue

            if ann.locus_tag is None:
                raise RuntimeError("Missing locus_tag for Prokka annotation")
            key = (accession, ann.locus_tag)
            family = seq_to_family.get(key)
            if family is None:
                family = f"{'R' if ann.is_rna else 'P'}0000000"
            base = f"{plasmid_prefix}_{ann.function_id}_{family}"
            feature_counts[base] += 1
            for mut in seq_to_mutations.get(key, []):
                feature_counts[f"{base}_{mut}"] += 1

        row = {"accession": accession, "gtdb_species": species, "gtdb_ani": ani}
        row.update(feature_counts)
        all_features.update(feature_counts.keys())
        rows.append(row)

    feature_cols = sorted(all_features)
    df = pd.DataFrame(rows)
    for col in feature_cols:
        if col not in df.columns:
            df[col] = 0
    if feature_cols:
        df[feature_cols] = df[feature_cols].fillna(0).astype(int)
    out_cols = ["accession", "gtdb_species", "gtdb_ani", *feature_cols]
    return df[out_cols]  # pyright: ignore[reportReturnType]


def main() -> None:
    """Main entry point executing the full tabulation workflow.

    The function parses command-line arguments, prepares the cache directory,
    reuses cached per-accession parse results when available, parses only new
    accessions in parallel, persists per-product FASTA files, runs the
    clustering/MSA/SNP pipeline (parallel across products), collects
    cluster-to-sequence mappings, compiles the feature table DataFrame, and
    writes the CSV and JSON outputs.
    """
    args = parse_args()
    indir = Path(args.INDIR)
    outprefix = Path(args.OUTPREFIX)
    csv_path = outprefix.with_suffix(".csv")
    json_path = outprefix.with_suffix(".json")
    cache_dir = Path(f"{outprefix}.cache")
    accession_cache_dir = cache_dir / "accessions"
    jobs = get_jobs()

    cache_dir.mkdir(parents=True, exist_ok=True)
    accession_cache_dir.mkdir(parents=True, exist_ok=True)

    accession_dirs = sorted(
        path
        for path in indir.iterdir()
        if path.is_dir() and not path.name.startswith(".")
    )
    if not accession_dirs:
        raise RuntimeError(f"No accessions found in {indir}")

    results_by_accession: dict[str, dict] = {}
    to_parse_dirs: list[Path] = []
    for accession_dir in accession_dirs:
        accession = accession_dir.name
        cached = load_cached_accession_raw(
            accession_cache_dir / f"{accession}.raw.json", accession
        )
        if cached is not None:
            results_by_accession[accession] = cached
            continue
        to_parse_dirs.append(accession_dir)

    if to_parse_dirs:
        with mp.Pool(processes=jobs) as pool:
            parsed_results = pool.map(parse_accession_raw, to_parse_dirs)
        for res in parsed_results:
            accession = str(res["accession"])
            results_by_accession[accession] = res
            write_cached_accession_raw(
                accession_cache_dir / f"{accession}.raw.json", res
            )

    raw_results = [results_by_accession[path.name] for path in accession_dirs]

    all_products = [row["product"] for res in raw_results for row in res["prokka_rows"]]
    existing_name_to_id = load_existing_prokka_mapping(json_path)
    product_to_id, id_to_product = build_prokka_id_map(
        all_products, existing_name_to_id
    )

    merged_by_accession: dict[str, list[RawAnnotation]] = {}
    meta_by_accession: dict[str, tuple[str, float]] = {}

    for res in raw_results:
        accession = res["accession"]
        meta_by_accession[accession] = (res["gtdb_species"], float(res["gtdb_ani"]))
        rgi_annotations: list[RawAnnotation] = []
        for row in res["rgi_rows"]:
            rgi_annotations.append(
                RawAnnotation(
                    accession=accession,
                    source="RGI",
                    contig=row["contig"],
                    start=row["start"],
                    end=row["end"],
                    strand=row["strand"],
                    is_rna=bool(row["is_rna"]),
                    is_plasmid=bool(row["is_plasmid"]),
                    function_id=f"A{row['aro']}",
                    mutations=tuple(row["mutations"]),
                )
            )

        prokka_annotations: list[RawAnnotation] = []
        for row in res["prokka_rows"]:
            pid = product_to_id[row["product"]]
            prokka_annotations.append(
                RawAnnotation(
                    accession=accession,
                    source="PROKKA",
                    contig=row["contig"],
                    start=row["start"],
                    end=row["end"],
                    strand=row["strand"],
                    is_rna=bool(row["is_rna"]),
                    is_plasmid=bool(row["is_plasmid"]),
                    function_id=f"P{pid}",
                    locus_tag=row["locus_tag"],
                    sequence=row["sequence"],
                )
            )

        merged: list[RawAnnotation] = list(rgi_annotations)
        for p_ann in prokka_annotations:
            should_replace = False
            for r_ann in rgi_annotations:
                if (
                    p_ann.contig == r_ann.contig
                    and p_ann.strand == r_ann.strand
                    and p_ann.is_rna == r_ann.is_rna
                    and overlap_ratio_union(
                        p_ann.start, p_ann.end, r_ann.start, r_ann.end
                    )
                    >= 0.95
                ):
                    should_replace = True
                    break
            if not should_replace:
                merged.append(p_ann)
        merged_by_accession[accession] = merged

    product_types = write_source_fastas(merged_by_accession, cache_dir)
    for accession, anns in merged_by_accession.items():
        merged_by_accession[accession] = [
            replace(ann, sequence=None) if ann.sequence is not None else ann
            for ann in anns
        ]

    raw_results.clear()
    wanted_fastas = {
        cache_dir / f"{pid}.{'fna' if is_rna else 'faa'}"
        for (pid, is_rna) in product_types
    }
    for stale in cache_dir.glob("[0-9][0-9][0-9][0-9][0-9][0-9][0-9].fna"):
        if stale not in wanted_fastas:
            stale.unlink()
    for stale in cache_dir.glob("[0-9][0-9][0-9][0-9][0-9][0-9][0-9].faa"):
        if stale not in wanted_fastas:
            stale.unlink()
    tasks = [(pid, is_rna, str(cache_dir)) for (pid, is_rna) in sorted(product_types)]

    if tasks:
        with mp.Pool(processes=jobs) as pool:
            pool.map(run_product_pipeline, tasks)

    seq_to_family, seq_to_mutations = collect_prokka_cluster_maps(
        cache_dir, product_types
    )

    df = build_outputs(
        merged_by_accession, meta_by_accession, seq_to_family, seq_to_mutations
    )
    df.to_csv(csv_path, index=False)

    with json_path.open("w", encoding="utf-8") as handle:
        json.dump(id_to_product, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
