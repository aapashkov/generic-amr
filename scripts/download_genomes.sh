#!/bin/bash
{
# Set Bash strict mode
set -euo pipefail
IFS=$'\t\n'

help="usage: ${0##*/} [-o OUTDIR] ACCESSION...

Download genome ACCESSIONs into OUTDIR, which defaults to working directory.
Valid ACCESSION strings include NCBI Genome/RefSeq identifiers (starting with
GCA_ or GCF_), ENA sequence assembly analysis identifiers (starting with ERZ),
or BV-BRC identifiers (two integers separated by a dot). Requires 'wget' or
'curl', and 'jq' to be available in PATH.

example: ${0##*/} -o genomes/ GCA_000005845.2 ERZ3086155 170673.13"

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

# Select a program to download files with
if command -v wget > /dev/null; then
  downloader="wget"
  params="-qO-"
elif command -v curl > /dev/null; then
  downloader="curl"
  params="-sL"
else
  printf '%s\n' 'error: neither wget nor curl is installed'
  exit 1
fi

# Check for jq
if ! command -v jq > /dev/null; then
  printf '%s\n' 'error: jq: command not found'
  exit 1
fi

for accession in "$@"; do

  # Skip if accession exists
  output="${outdir}/${accession}.fna"
  if [[ -e "$output" ]]; then
    printf '%s\n' "warning: '$output' exists, not overwriting"
    continue
  fi

  # First download to temp location, and move to final location on success
  tmp=$(mktemp "${outdir}/.tmp.XXXXXX")
  trap 'rm -f '"$tmp" EXIT

  # Try to download genome five times before exiting on error
  success=0
  for ((try=0;try<5;try++)); do

    # NCBI Genome/RefSeq identifiers
    if [[ $accession =~ ^(GCA|GCF)_[0-9]+\.[0-9]+$ ]]; then

      # Request FTP URL of genome
      url=$($downloader $params "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$accession/links" | \
        jq -r '.assembly_links[] | select(.assembly_link_type == "FTP_LINK") | .resource_link' 2> /dev/null || printf '')

      # Retry on error
      if [ -z "$url" ]; then
        printf '%s\n' "warning: retrying to download '$accession'"
        sleep 1
        continue
      fi

      # Extract URL from response and download genome
      url="${url/https:/ftp:}/${url##*/}_genomic.fna.gz"
      if $downloader $params "$url" | gunzip -c > "$tmp"; then
        success=1
        break

      # Retry on error
      else
        printf '%s\n' "warning: retrying to download '$accession'"
        sleep 1
        continue
      fi

    # ENA sequence assembly analysis identifiers
    elif [[ $accession =~ ^ERZ[0-9]+$ ]]; then

      # Request FTP URL of genome
      url=$($downloader $params "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${accession}&result=analysis&fields=analysis_type,submitted_ftp" | \
        grep -m 1 SEQUENCE_ASSEMBLY || printf '')

      # Retry on error
      if [ -z "$url" ]; then
        printf '%s\n' "warning: retrying to download '$accession'"
        sleep 1
        continue
      fi

      # Extract URL from response and download genome
      url="ftp://ftp${url#*ftp}"
      url="${url%%	*}"
      if $downloader $params "$url" | gunzip -c > "$tmp"; then
        success=1
        break

      # Retry on error
      else
        printf '%s\n' "warning: retrying to download '$accession'"
        sleep 1
        continue
      fi

    # BV-BRC identifiers
    elif [[ $accession =~ ^[0-9]+\.[0-9]+$ ]]; then
      url="https://www.bv-brc.org/api/genome_sequence/?eq(genome_id,${accession})"

      if $downloader $params --header 'Accept: text/csv' "$url" | grep -qFm1 "$accession"; then
        $downloader $params "$url" | \
          jq -r '.[] | ">" + .sequence_id + "\n" + (.sequence | ascii_upcase)' \
          > "$tmp"
        success=1
        break

      # Retry on error
      else
        printf '%s\n' "warning: retrying to download '$accession'"
        sleep 1
        continue
      fi

    # Invalid identifiers
    else
      printf '%s\n' "error: '$accession' is not a valid identifier"
      exit 1
    fi
  done

  if [[ $success == 0 ]]; then
    printf '%s\n' "error: could not fetch assembly for '$accession'"
    exit 1
  fi

  sleep 1
  mv "$tmp" "$output"
  printf '%s\n' "info: downloaded '$accession' into '$output'"

done

exit 0
} >&2
