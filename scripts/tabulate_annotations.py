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
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd
from BCBio import GFF
from Bio import SeqIO


RGI_MUT_SPLIT_RE = re.compile(r"\s*,\s*")
RGI_MUT_RE = re.compile(r"^([A-Za-z\-*])(\d+)([A-Za-z\-*])$")


@dataclass(frozen=True)
class RawAnnotation:
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
    parser = argparse.ArgumentParser(
        description="Tabulate genome annotations into feature counts."
    )
    parser.add_argument("INDIR", help="Input directory with per-accession annotations")
    parser.add_argument("OUTPREFIX", help="Output prefix for CSV/JSON/cache")
    return parser.parse_args()


def get_jobs() -> int:
    try:
        jobs = int(os.environ.get("JOBS", ""))
        if jobs > 0:
            return jobs
    except ValueError:
        pass
    return 4


def parse_plasmids(platon_tsv: Path) -> set[str]:
    if not platon_tsv.exists() or platon_tsv.stat().st_size == 0:
        return set()
    with platon_tsv.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["ID"] for row in reader if row.get("ID")}


def parse_gtdb_info(sourmash_csv: Path) -> tuple[str, float]:
    with sourmash_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        row = next(reader, None)
        if row is None:
            raise ValueError(f"{sourmash_csv} has no rows")
        name_tokens = row.get("name", "").split()
        if len(name_tokens) <= 1:
            raise ValueError(f"Could not parse species from sourmash name: {row.get('name')}")
        tokens = name_tokens[1:]
        if tokens and tokens[0].lower() == "uncultured":
            tokens = tokens[1:]
        if len(tokens) < 2:
            raise ValueError(f"Could not derive species from sourmash name: {row.get('name')}")
        species = f"{tokens[0]} {tokens[1]}"
        ani_raw = row.get("ani", "")
        ani = float(ani_raw) if ani_raw not in {"", None} else math.nan
        return species, ani


def parse_rgi_mutations(*mut_fields: str) -> tuple[str, ...]:
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
    match = RGI_MUT_RE.match(token)
    if not match:
        return None
    ref, pos, alt = match.groups()
    return f"{ref.upper()}{int(pos):07d}{alt.upper()}"


def normalize_aro(aro_value: str) -> str:
    digits = re.sub(r"\D", "", aro_value or "")
    if not digits:
        raise ValueError(f"Invalid ARO value: {aro_value!r}")
    if len(digits) > 7:
        digits = digits[-7:]
    return digits.zfill(7)


def extract_feature_sequence(record_seq, feature, is_rna: bool) -> str:
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
    t = (feature_type or "").lower()
    return t != "cds" and "rna" in t


def parse_accession_raw(accession_dir: Path) -> dict:
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
                product = str(qualifiers.get("product", ["hypothetical protein"])[0]).strip()
                locus_tag = str(qualifiers.get("locus_tag", [qualifiers.get("ID", [""])[0]])[0])
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


def overlap_ratio_union(a_start: int, a_end: int, b_start: int, b_end: int) -> float:
    overlap = max(0, min(a_end, b_end) - max(a_start, b_start) + 1)
    if overlap == 0:
        return 0.0
    span = max(a_end, b_end) - min(a_start, b_start) + 1
    return overlap / span


def load_existing_prokka_mapping(json_path: Path) -> dict[str, str]:
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
) -> dict[tuple[str, bool], list[RawAnnotation]]:
    by_product_type: dict[tuple[str, bool], list[RawAnnotation]] = defaultdict(list)
    for anns in merged_by_accession.values():
        for ann in anns:
            if ann.source != "PROKKA":
                continue
            by_product_type[(ann.function_id[1:], ann.is_rna)].append(ann)

    for (pid, is_rna), anns in by_product_type.items():
        ext = "fna" if is_rna else "faa"
        out_path = cache_dir / f"{pid}.{ext}"
        lines = []
        for ann in anns:
            seq_id = f"{ann.accession}__{ann.locus_tag}"
            lines.append(f">{seq_id}\n{ann.sequence}\n")
        new_content = "".join(lines)
        if out_path.exists() and out_path.read_text(encoding="utf-8") == new_content:
            continue
        out_path.write_text(new_content, encoding="utf-8")
    return by_product_type


def run_cmd(cmd: list[str], cwd: Path | None = None, stdout_path: Path | None = None) -> None:
    if stdout_path is None:
        subprocess.run(cmd, check=True, cwd=cwd)
        return
    with stdout_path.open("w", encoding="utf-8") as out_handle:
        subprocess.run(cmd, check=True, cwd=cwd, stdout=out_handle)


def split_mmseqs_msa(msa_path: Path, out_prefix: str) -> list[Path]:
    ext = "fna" if msa_path.name.endswith(".fna") else "faa"
    old_files = glob.glob(f"{out_prefix}.*.msa.{ext}")
    for file_path in old_files:
        os.unlink(file_path)

    records: list[tuple[str, str]] = []
    header = None
    seq_lines: list[str] = []
    with msa_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_lines)))
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header is not None:
            records.append((header, "".join(seq_lines)))

    clusters: list[list[tuple[str, str]]] = []
    current: list[tuple[str, str]] = []
    seen: set[str] = set()
    for rec_header, seq in records:
        if not seq:
            if current:
                clusters.append(current)
                current = []
                seen = set()
            continue
        if rec_header in seen and current:
            clusters.append(current)
            current = []
            seen = set()
        current.append((rec_header, seq))
        seen.add(rec_header)
    if current:
        clusters.append(current)

    out_paths: list[Path] = []
    for idx, cluster in enumerate(clusters):
        out_path = Path(f"{out_prefix}.{idx:07d}.msa.{ext}")
        with out_path.open("w", encoding="utf-8") as handle:
            for rec_header, seq in cluster:
                handle.write(f">{rec_header}\n{seq}\n")
        out_paths.append(out_path)
    return out_paths


def msa_has_accession_headers(msa_path: Path) -> bool:
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
    count = 0
    with fasta_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def run_product_pipeline(task: tuple[str, bool, str]) -> tuple[str, bool]:
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
    seq_count = fasta_record_count(source_fasta)

    if seq_count <= 1:
        if (not flat_msa.exists()) or source_fasta.stat().st_mtime > flat_msa.stat().st_mtime:
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

    cluster_paths = split_mmseqs_msa(flat_msa, str(cache_dir / pid))
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


def parse_vcf_mutations(vcf_path: Path) -> dict[str, list[str]]:
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


def collect_prokka_cluster_maps(
    cache_dir: Path, by_product_type: dict[tuple[str, bool], list[RawAnnotation]]
) -> tuple[dict[tuple[str, str], str], dict[tuple[str, str], list[str]]]:
    seq_to_family: dict[tuple[str, str], str] = {}
    seq_to_mutations: dict[tuple[str, str], list[str]] = defaultdict(list)

    for pid, is_rna in sorted(by_product_type):
        ext = "fna" if is_rna else "faa"
        cluster_paths = sorted(Path(p) for p in glob.glob(str(cache_dir / f"{pid}.*.msa.{ext}")))
        for cluster_path in cluster_paths:
            family_idx = cluster_path.name.split(".")[1]
            family = f"{'R' if is_rna else 'P'}{family_idx}"
            for rec in SeqIO.parse(str(cluster_path), "fasta"):
                if "__" not in rec.id:
                    continue
                accession, locus_tag = rec.id.split("__", 1)
                seq_to_family[(accession, locus_tag)] = family

            vcf_path = cache_dir / cluster_path.name.replace(f".msa.{ext}", ".vcf")
            sample_muts = parse_vcf_mutations(vcf_path)
            for sample, muts in sample_muts.items():
                if "__" not in sample:
                    continue
                accession, locus_tag = sample.split("__", 1)
                seq_to_mutations[(accession, locus_tag)].extend(muts)
    return seq_to_family, seq_to_mutations


def build_outputs(
    merged_by_accession: dict[str, list[RawAnnotation]],
    meta_by_accession: dict[str, tuple[str, float]],
    seq_to_family: dict[tuple[str, str], str],
    seq_to_mutations: dict[tuple[str, str], list[str]],
) -> pd.DataFrame:
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
    args = parse_args()
    indir = Path(args.INDIR)
    outprefix = Path(args.OUTPREFIX)
    csv_path = outprefix.with_suffix(".csv")
    json_path = outprefix.with_suffix(".json")
    cache_dir = Path(f"{outprefix}.cache")
    jobs = get_jobs()

    cache_dir.mkdir(parents=True, exist_ok=True)

    accession_dirs = sorted(
        path for path in indir.iterdir() if path.is_dir() and not path.name.startswith(".")
    )
    if not accession_dirs:
        raise RuntimeError(f"No accessions found in {indir}")

    with mp.Pool(processes=jobs) as pool:
        raw_results = pool.map(parse_accession_raw, accession_dirs)

    all_products = [row["product"] for res in raw_results for row in res["prokka_rows"]]
    existing_name_to_id = load_existing_prokka_mapping(json_path)
    product_to_id, id_to_product = build_prokka_id_map(all_products, existing_name_to_id)

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
            replace = False
            for r_ann in rgi_annotations:
                if (
                    p_ann.contig == r_ann.contig
                    and p_ann.strand == r_ann.strand
                    and p_ann.is_rna == r_ann.is_rna
                    and overlap_ratio_union(p_ann.start, p_ann.end, r_ann.start, r_ann.end) >= 0.95
                ):
                    replace = True
                    break
            if not replace:
                merged.append(p_ann)
        merged_by_accession[accession] = merged

    by_product_type = write_source_fastas(merged_by_accession, cache_dir)
    wanted_fastas = {
        cache_dir / f"{pid}.{'fna' if is_rna else 'faa'}" for (pid, is_rna) in by_product_type
    }
    for stale in cache_dir.glob("[0-9][0-9][0-9][0-9][0-9][0-9][0-9].fna"):
        if stale not in wanted_fastas:
            stale.unlink()
    for stale in cache_dir.glob("[0-9][0-9][0-9][0-9][0-9][0-9][0-9].faa"):
        if stale not in wanted_fastas:
            stale.unlink()
    tasks = [(pid, is_rna, str(cache_dir)) for (pid, is_rna) in sorted(by_product_type.keys())]

    if tasks:
        with mp.Pool(processes=jobs) as pool:
            pool.map(run_product_pipeline, tasks)

    seq_to_family, seq_to_mutations = collect_prokka_cluster_maps(cache_dir, by_product_type)

    df = build_outputs(merged_by_accession, meta_by_accession, seq_to_family, seq_to_mutations)
    df.to_csv(csv_path, index=False)

    with json_path.open("w", encoding="utf-8") as handle:
        json.dump(id_to_product, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
