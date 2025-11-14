#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Create a pangenome-like table from the annotations created using the
annotate_genome.sh script in INDIR and save results into OUTDIR. If the JOBS
environment variable is set, will use JOBS processes, and four otherwise.
Requires mmseqs, proteinortho, mafft and snp-sites to be available in PATH."""

# Using Python 3.13.5
# Standard library imports
import atexit
import collections as cls
import csv
import functools as ft
import glob
import multiprocessing as mp
import os
import re
import shutil
import subprocess as sp
import sys
import tempfile

# Third-party library imports
import pandas as pd     # 2.2.3
pd.set_option('future.no_silent_downcasting', True)

# Globals
prokka_managed_dict = mp.Manager().dict()
prokka_managed_dict["hypothetical protein"] = 0


def check_dependencies() -> bool:
    """Checks for executable dependencies, returning False if anything is
    missing, and True otherwise."""

    deps = ["mmseqs", "proteinortho", "mafft", "snp-sites"]
    missing = []

    for dep in deps:
        if shutil.which(dep) is None:
            missing.append(dep)

    if missing:
        print(
            f"error: missing executable dependencies: {', '.join(missing)}",
            file=sys.stderr
        )
        return False
    return True


def mmseqs_gpu_available() -> bool:
    """Returns True if nvidia-smi exists and mmseqs has GPU support."""

    cmd = ["mmseqs", "easy-search", "--help"]
    nvidia_smi = False if shutil.which("nvidia-smi") is None else True
    return "--gpu" in sp.run(cmd, stdout=sp.PIPE).stdout.decode() and nvidia_smi


def get_species(genome: str, indir: str):
    """Returns a dictionary with best GTDB species match of a genome with its
    associated ANI."""

    file = os.path.join(indir, genome, f"{genome}.sourmash.csv")
    
    # Best match is listed as first entry in Sourmash CSV output
    with open(file) as handle:
        reader = csv.DictReader(handle)
        best_match = next(reader)

    # Remove GenBank identifier and the "uncultured" keyword from species name
    species = best_match["name"].split(" ", 1)[1].removeprefix("uncultured ")
    ani = float(best_match["ani"])

    return {"gtdb_species": species, "gtdb_ani": ani}


def get_overlap(query: dict, target: dict) -> float:
    """Given two dictionaries with 'contig', 'start', 'stop', 'strand' and
    'type' keys describing two sequence features, determine their overlap."""

    # Skip obvious non-overlapping features
    if (
        query["contig"] != target["contig"] or
        query["strand"] != target["strand"] or
        query["type"] != target["type"]
    ):
        return 0.0

    # Calculate overlap for features on the same contig and strand
    query_length = query["stop"] - query["start"]
    target_length = target["stop"] - target["start"]
    start_overlap = max(query["start"], target["start"])
    stop_overlap = min(query["stop"], target["stop"])
    overlap_length = max(0, stop_overlap - start_overlap)
    total_length = query_length + target_length - overlap_length

    return overlap_length / total_length


def get_resolved_annotations(genome: str, indir: str):
    """Returns a list of annotations by resolving annotations from Platon,
    Prokka, and RGI."""

    # Create a set with plasmid contig names found in genome
    platon_file = os.path.join(indir, genome, f"{genome}.platon.tsv")

    with open(platon_file) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        plasmids = {plasmid["ID"] for plasmid in reader}

    # Retrieve RGI non-loose predictions with coordinates
    rgi_file = os.path.join(indir, genome, f"{genome}.rgi.txt")
    rgi_preds = {
        "contig": [], "start": [], "stop": [], "features": [], "strand": [],
        "type": []
    }

    with open(rgi_file) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            if row["Cut_Off"] == "Loose": continue
            plasmid = "P" if row["Contig"] in plasmids else "N"
            gene_type = "R" if row["Model_type"].startswith("rRNA") else "P"
            features = [f"{plasmid}_A{row['ARO']}_{gene_type}0"]

            # List all SNPs, remove n/a as it represents nothing
            snps = (
                set(row["SNPs_in_Best_Hit_ARO"].split(", ")) |
                set(row["Other_SNPs"].split(", "))
            ) - {"n/a"}

            # If SNPs are found, add SNP features to feature list
            for snp in snps:
                features.append(f"{features[0]}_A{snp}")

            rgi_preds["contig"].append(row["Contig"])
            rgi_preds["start"].append(int(row["Start"]))
            rgi_preds["stop"].append(int(row["Stop"]))
            rgi_preds["features"].append(features)
            rgi_preds["strand"].append(row["Orientation"])
            rgi_preds["type"].append(gene_type)

    # Retrieve Prokka predictions and merge with RGI predictions
    prokka_file = os.path.join(indir, genome, f"{genome}.prokka.gff")
    final_preds = cls.Counter()
    # If overlap is greater than this, it refers to the same location
    overlap_threshold = 0.95

    with open(prokka_file) as handle:

        # Ignore commented lines in file
        reader = csv.reader(
            (line for line in handle if not line.startswith("#")),
            delimiter="\t"
        )

        for row in reader:

            # Stop on fasta header
            if row[0].startswith(">"):
                break

            query = {
                "contig": row[0],
                "start": int(row[3]),
                "stop": int(row[4]),
                "strand": row[6],
                "type": "P" if row[2] == "CDS" else "R"
            }

            # Check if RGI has a prediction on the same location
            best_overlap = 0.0
            best_index = -1
            for index in range(len(rgi_preds["contig"])):
                target = {
                    "contig": rgi_preds["contig"][index],
                    "start": rgi_preds["start"][index],
                    "stop": rgi_preds["stop"][index],
                    "strand": rgi_preds["strand"][index],
                    "type": rgi_preds["type"][index],
                }
                overlap = get_overlap(query, target)
                if overlap > best_overlap:
                    best_overlap = overlap
                    best_index = index

            if best_overlap > overlap_threshold:
                for pred in rgi_preds["features"][best_index]:
                    final_preds[pred] += 1
                continue

            # Update global Prokka product dictionary on new products
            product = dict(el.split("=") for el in row[8].split(";"))["product"]
            if product not in prokka_managed_dict:
                prokka_managed_dict[product] = len(prokka_managed_dict)
            product_id = prokka_managed_dict[product]

            # Create feature name
            plasmid = "P" if query["contig"] in plasmids else "N"
            feature = f"{plasmid}_P{product_id}_{query['type']}0"
            final_preds[feature] += 1

    return final_preds


def main() -> int:
    """Driver code."""

    prog = os.path.basename(__file__)
    jobs_env = os.environ.get("JOBS")
    if jobs_env == "" or jobs_env is None:
        jobs = os.cpu_count()
    else:
        jobs = int(jobs_env)

    # Print help message if incorrect number of parameters or missing deps
    if len(sys.argv) != 3 or not check_dependencies():
        print(
            f"usage: {prog} INDIR OUTDIR",
            __doc__,
            f"example: JOBS=10 {prog} annotations/ pangenome/",
            file=sys.stderr,
            sep="\n\n"
        )
        return 1

    # Handle input, output and temporary directories
    indir = os.path.abspath(sys.argv[1])
    outdir = os.path.abspath(sys.argv[2])
    tmp = tempfile.mkdtemp(prefix=".tmp.", dir=os.path.dirname(outdir))
    rmtree = ft.partial(shutil.rmtree, ignore_errors=True)
    atexit.register(rmtree, tmp)
    os.chdir(tmp)

    if os.path.exists(outdir):
        print(f"error: '{sys.argv[2]}' exists, not overwriting")
        return 1

    # List annotations and prepare data for output table
    genomes = glob.glob("*", root_dir=indir)
    data = {genome: {} for genome in genomes}

    with mp.Pool(jobs) as pool:

        # Get GTDB species match and ANI
        func = ft.partial(get_species, indir=indir)
        species = pool.map(func, genomes)
        data.update(dict(zip(genomes, species)))
        
        # Resolve annotations
        func = ft.partial(get_resolved_annotations, indir=indir)
        annotations = pool.map(func, genomes)
        for genome, annotation in zip(genomes, annotations):
            data[genome].update(annotation)

    # Save outputs
    df = pd.DataFrame(data).T.fillna(0)
    df.index.name = "accession"
    df.to_csv(os.path.join(tmp, "table.csv"))
    with open(os.path.join(tmp, "prokka.txt"), "w") as handle:
        prokka = sorted(prokka_managed_dict, key=prokka_managed_dict.get) # type: ignore
        print(*prokka, sep="\n", file=handle)
    shutil.move(tmp, outdir)

    return 0


if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw|\.exe)?$', '', sys.argv[0])
    sys.exit(main())
