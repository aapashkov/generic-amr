#!/bin/bash
{
set -euo pipefail
IFS=$'\t\n'

help="usage: ${0##*/} [-o OUTDIR] ACCESSION...

Download SRA ACCESSIONs into OUTDIR, which defaults to working directory.
Requires 'fasterq-dump' and 'gzip' to be available in PATH.

example: ${0##*/} -o reads/ SRR33085175 ERR12520663"

# Parse output directory, removing trailing '/'
outdir='.'
while getopts ':o:' opt; do
  case ${opt} in
    o)
      outdir="$OPTARG"
      ;;
    \?)
      printf 'error: invalid option: -%s\n' "$OPTARG"
      exit 1
      ;;
  esac
done
outdir="${outdir%/}"

# Print help message if insufficient number of positional arguments
shift "$((OPTIND - 1))"
if [ "$#" -lt 1 ]; then
  printf '%s\n' "$help"
  exit 1
fi

# Check that output directory is writeable
if [[ -d "$outdir" && -w "$outdir" ]]; then
  :
else
  printf '%s\n' "error: '$outdir' is not a writeable directory"
  exit 1
fi

# Check for dependencies
for cmd in fasterq-dump gzip; do
  if ! command -v $cmd > /dev/null; then
    printf '%s\n' "error: $cmd: command not found"
    exit 1
  fi
done

for accession in "$@"; do
  # Skip if accession exists
  output="${outdir}/${accession}"
  if [[ -e "$output" ]]; then
    printf '%s\n' "warning: '$output' exists, not overwriting"
    continue
  fi

  # First download to temp location, and move to final location on success
  tmp=$(mktemp -d "${outdir}/.tmp.XXXXXX")
  trap 'rm -rf '"$tmp" EXIT

  # Download reads
  success=0
  for ((try=0;try<10;try++)); do
    # shellcheck disable=SC2016
    if prefetch --progress --output-directory "$tmp" "$accession" && \
      fasterq-dump --seq-defline '@$ac.$si' \
        --qual-defline '+' \
        --details \
        --outdir "$tmp/$accession" \
        --temp "$tmp" \
        "$tmp/$accession/${accession}.sra"; then
      success=1
      break
    fi
    printf '%s\n' "warning: retrying to fetch reads for '$accession'"
    sleep 1
  done

  if [[ $success == 0 ]]; then
    printf '%s\n' "error: could not fetch reads for '$accession'"
    exit 1
  fi

  # Append '.1' and '.2' to read IDs if files include paired-end data.
  # For example: @ERR12520663.5 → @ERR12520663.5.1 or @ERR12520663.5.2
  paired=false
  for pair in 1 2; do
    file="$tmp/$accession/${accession}_${pair}.fastq"
    if [[ -e "$file" ]]; then
      paired=true
      sed -i '/^@'"$accession"'/s/$/.'"$pair"'/' "$file"
    fi
  done

  # If both paired and unpaired files exist, append '.3' to the unpaired IDs.
  # This is a requirement for the MIRA assembler.
  file="$tmp/$accession/${accession}.fastq"
  if [[ "$paired" == "true" && -e "$file" ]]; then
    sed -i '/^@'"$accession"'/s/$/.3/' "$file"
  fi

  # Compress and move outputs out of temp directory
  gzip -v "$tmp/$accession"/*.fastq
  rm -v "$tmp/$accession/${accession}.sra"
  mv -v "$tmp/$accession" "$output"
  rm -rf "$tmp"

done

exit 0
} >&2
