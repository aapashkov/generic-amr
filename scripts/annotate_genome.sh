#!/bin/bash
{
# Set Bash strict mode
set -euo pipefail
IFS=$'\t\n'

help="usage: ${0##*/} INPUT OUTDIR SBT PLATON

Annotate an INPUT fasta genome and save results to OUTDIR. Must provide a
Sourmash GTDB SBT zip file for taxonomic classification, which you can find
here: https://sourmash.readthedocs.io/en/latest/databases.html. The PLATON
database must also be provided for plasmid detection. The RGI database must be
loaded for global usage. Requires 'prokka', 'rgi', 'barrnap', 'cmscan',
'platon', and 'sourmash' to be available in PATH. 

example: ${0##*/} genome.fna annotation/ gtdb.sbt.zip"

# Print help message if insufficient number of positional arguments
if [ "$#" != 4 ]; then
  printf '%s\n' "$help"
  exit 1
fi

# Make sure output directory does not exist
if [ -e "$2" ]; then
  printf '%s\n' "error: '$2' exists, not overwriting"
  exit 1
fi

input=$(readlink -f "$1")
outdir=$(readlink -f "$2")
base=$(basename "$outdir")
sbt=$(readlink -f "$3")
platon=$(readlink -f "$4")

# Check dependencies
for prog in prokka rgi barrnap cmscan platon sourmash; do
  if ! command -v "$prog" > /dev/null; then
    printf 'error: %s: command not found\n' "$prog"
    exit 1
  fi
done

# Create temporary directory
tmp=$(mktemp -d "$(dirname "$outdir")/.tmp.XXXXXX")
trap 'rm -rf "$tmp"' EXIT
cd "$tmp"

# Classify genome by sketching it with Sourmash
sourmash sketch dna -o "${base}.sourmash.sig" --name "$base" "$input"
sourmash search -o "${base}.sourmash.csv" --best-only --threshold 0.0 \
  "${base}.sourmash.sig" "$sbt"
species=$(sed -n '2p' "${base}.sourmash.csv" | cut -f4 -d,)

# Remove GenBank identifier and 'uncultured' keyword and split species name
species="${species#* }"
species="${species#uncultured }"
genus="${species%% *}"
species="${species#* }"
species="${species%% *}"

# Run annotation with Prokka
prokka --outdir . --prefix "${base}.prokka" --locustag GAMR --force --rfam \
  --genus "$genus" --species "$species" --strain "$base" --cpus 1 "$input"

# Run annotation with RGI, safely removing the '.fai' file it sometimes creates
if [ -e "${input}.fai" ]; then
  index_exists=true
fi
rgi main --include_loose --debug -i "$input" -o "${base}.rgi" -n 1 --clean \
  2>&1 | tee "${base}.rgi.log"
if [ -z ${index_exists+x} ]; then
  rm -f "${input}.fai"
fi

# Run annotation with Platon
platon -vp "${base}.platon" -t 1 "$input" -d "$platon"

# Move results out of temporary directory
mv -v "$tmp" "$outdir"

exit 0
} >&2
