#!/bin/bash
{
# Set Bash strict mode
set -euo pipefail
IFS=$'\t\n'

# Prints a message to stderr
log () {
  printf '%s\n' "${1}" >&2
}

# Print help message if incorrect number of arguments and exit
help="usage: $(basename "$0") FASTA OUTDIR SBT

Annotate a bacterial FASTA genome and save results to OUTDIR. Requires passing
the Sourmash SBT zip file built from GTDB RS266, which you can download here:
https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db.new/gtdb-rs226/gtdb-reps-rs226-k31.dna.zip.
Requires the following commands to be available in PATH: prokka, rgi, barrnap, cmscan, sourmash

example: $(basename "$0") genome.fna result gtdb.sbt.zip"
if [ "$#" != 3 ]; then
  log "$help"
  exit 1
fi

input=$(readlink -f "$1")
outdir=$(readlink -f "$2")
outbase=$(basename "$outdir")
sbt=$(readlink -f "$3")

# Check dependencies
for prog in prokka rgi barrnap cmscan sourmash; do
  if ! command -v "$prog" > /dev/null; then
    log 'error: missing required dependency '\'"$prog"\'
    exit 1
  fi
done

# Make sure output directory does not exist
if [ -e "$outdir" ]; then
  log 'error: output '\'"$outdir"\'' already exists, not overwriting'
  exit 1
fi

# Validate SBT file
# if ! md5sum -c <<< "179760fcd05ff75b7515a16e46ad537c ${sbt}" 2>&1; then
#   log 'error: the file '\'"$sbt"\'' is not a GTDB RS226 SBT file with representatives of kmer size 31'
#   exit 1
# fi

# Create temporary directory and set its removal on exit
tmp=$(mktemp -d "$(dirname "$outdir")/.tmp.XXXXXX")
# trap 'rm -rf "$tmp"' EXIT
cd "$tmp"

# Classify genome by sketching it
IFS=',' read -r _ _ _ species _ _ _ _ <<< "$(
  sourmash sketch dna -o - "$input" | \
  sourmash search -o - --best-only - "$sbt" | \
  tail -1
)"

# Remove GenBank identifier and 'uncultured' keyword and split species name
species="${species#* }"
species="${species#uncultured }"
genus="${species%% *}"
species="${species#* }"
species="${species%% *}"

# Run annotation with Prokka
prokka --outdir "$outbase" --prefix "$outbase" --locustag GAMR \
  --genus "$genus" --species "$species" --strain "$outbase" --cpus 10 --rfam \
  "$input"

exit 0
}
