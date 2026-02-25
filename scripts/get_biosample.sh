#!/bin/bash
{
# Set Bash strict mode
set -euo pipefail
IFS=$'\t\n'

help="usage: ${0##*/} ACCESSION...

Retrieve the BioSamples of genome ACCESSIONs and print them to stdout in TSV.
Valid ACCESSION strings include NCBI Genome/RefSeq identifiers (starting with
GCA_ or GCF_), ENA sequence assembly analysis identifiers (starting with ERZ),
or BV-BRC identifiers (starting with BVBRC_). Requires 'wget' or
'curl', and 'jq' to be available in PATH.

example: ${0##*/} GCA_000005845.2 ERZ3086155 BVBRC_170673.13"

# Print help message if insufficient number of positional arguments
if [ "$#" -lt 1 ]; then
  printf '%s\n' "$help" >&2
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
  printf '%s\n' 'error: neither wget nor curl is installed' >&2
  exit 1
fi

# Check for jq
if ! command -v jq > /dev/null; then
  printf '%s\n' 'error: jq: command not found' >&2
  exit 1
fi

for accession in "$@"; do

  # Try to retrieve metadata five times before exiting on error
  success=0
  biosample=''
  for ((try=0;try<5;try++)); do

    # NCBI Genome/RefSeq identifiers
    if [[ $accession =~ ^(GCA|GCF)_[0-9]+\.[0-9]+$ ]]; then
      url="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${accession}/dataset_report"
      selector='.reports[0].assembly_info.biosample.accession'

    # ENA sequence assembly analysis identifiers
    elif [[ $accession =~ ^ERZ[0-9]+$ ]]; then
      url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${accession}&result=analysis&fields=sample_accession&format=json"
      selector='.[0].sample_accession'

    # BV-BRC identifiers
    elif [[ $accession =~ ^BVBRC_[0-9]+\.[0-9]+$ ]]; then
      acc=${accession#*_}
      url="https://www.bv-brc.org/api/genome/${acc}"
      selector='.biosample_accession'

    # Invalid identifiers
    else
      printf '%s\n' "error: '$accession' is not a valid identifier" >&2
      exit 1
    fi

    # Run retrieval command
    if biosample=$($downloader $params "$url" | jq -r "$selector"); then
      success=1
      break
    else
      printf '%s\n' "warning: retrying to retrieve '$accession'" >&2
      sleep 1
      continue
    fi
  done

  if [[ $success == 0 ]]; then
    printf '%s\n' "error: could not fetch metadata for '$accession'" >&2
    exit 1
  fi

  printf '%s\t%s\n' "$accession" "$biosample"
  sleep 1
done
}
