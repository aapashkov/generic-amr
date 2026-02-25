#!/bin/bash
{
set -euo pipefail
IFS=$'\t\n'

help="usage: ${0##*/} INDIR OUTDIR KRAKENDB [PLATFORM]

Trim reads from INDIR and save results to OUTDIR. Files in INDIR must comprise
one to three fastq.gz files and should start with the BASEname of INDIR, with
the following structure:

  BASE/BASE.fastq.gz   → single/unpaired reads
  BASE/BASE_1.fastq.gz → forward reads
  BASE/BASE_2.fastq.gz → reverse reads

Read identifiers must match '@BASE.X[.Y]', where X is the read ID, and Y must be
1, 2, or 3, respectively corresponding to forward, reverse and unpaired reads.
Single reads (without paired-end data) must not contain '.Y'. A KRAKENDB with at
least human sequences is required to remove human-mapping reads. When BASE is an
accession from SRA, PLATFORM is not necessary as it is fetched directly from
SRA. Otherwise, or if you want to overwrite it, PLATFORM must be provided and be
one of: ILLUMINA, ION_TORRENT, LS454, OXFORD_NANOPORE, PACBIO_SMRT. Requires
'efetch' and 'xtract' (from NCBI Entrez Direct), 'trim_galore', 'pear',
'porechop', 'seqkit', 'fastp', 'kraken2' to be available in PATH.

example: ${0##*/} raw/SRR33085175/ trimmed/ krakendb/ OXFORD_NANOPORE"

# Print help message if insufficient number of positional arguments
if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
  printf '%s\n' "$help"
  exit 1
fi
indir=$(readlink -f "$1")
outdir=$(readlink -f "$2")
krakendb=$(readlink -f "$3")
user_platform=${4:-UNDEFINED}
base=$(basename "$indir")

# Check for dependencies
for cmd in efetch xtract trim_galore pear porechop seqkit fastp kraken2; do
  if ! command -v $cmd > /dev/null; then
    printf '%s\n' "error: $cmd: command not found"
    exit 1
  fi
done

# Check that output directory is writeable
if [[ -d "$outdir" && -w "$outdir" ]]; then
  :
else
  printf '%s\n' "error: '$2' is not a writeable directory"
  exit 1
fi

# Check if output exists
if [[ -e "$outdir/$base" ]]; then
  printf '%s\n' "error: '$(basename "$2")/$base' exists, not overwriting"
  exit 1
fi

# Make sure "Homo sapiens" exists in provided Kraken2 database
if [[ -e "$krakendb/.gamr-validated" ]] || \
  kraken2-inspect --db "$krakendb" | grep -q 'Homo sapiens' && \
  touch "$krakendb/.gamr-validated"; then

  printf '%s\n' "info: '$3' is a valid database"
else
  printf '%s\n' "error: '$3' is missing human data"
  exit 1
fi

# Try to fetch platform from SRA
platform=$(
  efetch -db sra -id "$base" -format xml 2> /dev/null |
  xtract -pattern PLATFORM -element '$'
)

# If user has set a platform, overwrite it
if [[ "$user_platform" != "UNDEFINED" ]]; then
  printf '%s\n' "info: overwriting platform '$platform' to '$user_platform'"
  platform=$user_platform
else
  printf '%s\n' "info: fetched platform '$platform'"
fi

# Exit if an incompatible platform is provided
case "$platform" in
  ILLUMINA|ION_TORRENT|LS454|OXFORD_NANOPORE|PACBIO_SMRT)
    :
    ;;
  *)
    printf '%s' "error: '$platform' is not a valid platform, must be one of: "
    printf '%s\n' "ILLUMINA, ION_TORRENT, LS454, OXFORD_NANOPORE, PACBIO_SMRT."
    exit 1
    ;;
esac

# Create temporary directory
tmp=$(mktemp -d "$(dirname "$outdir")/.tmp.XXXXXX")
trap 'rm -rf "$tmp"' EXIT
cd "$tmp"

if [[ "$platform" == "ILLUMINA" ]]; then

  # Trim-galore already adapter and quality trims Illumina data with sensible
  # defaults. We simply change to paired configuration when needed
  if [[ -e "$indir/${base}.fastq.gz" ]]; then
    trim_galore \
      --length 40 \
      --no_report_file \
      --basename "$base" \
      --cores 1 \
      "$indir/${base}.fastq.gz"
    mv -v "${base}_trimmed.fq.gz" "${base}.fastq.gz"
  fi
  if [[ -e "$indir/${base}_1.fastq.gz" && -e "$indir/${base}_2.fastq.gz" ]]; then
    trim_galore \
      --length 40 \
      --no_report_file \
      --basename "$base" \
      --cores 1 \
      --paired \
      "$indir/${base}_1.fastq.gz" "$indir/${base}_2.fastq.gz"
    mv -v "${base}_val_1.fq.gz" "${base}_1.fastq.gz"
    mv -v "${base}_val_2.fq.gz" "${base}_2.fastq.gz"
  fi
elif [[ "$platform" == "OXFORD_NANOPORE" ]]; then

  # Some Oxford Nanopore runs may have paired-end reads, we merge those by pair
  # and concatenate them with any unpaired data we might have
  input="$indir/${base}.fastq.gz"
  if [[ -e "$indir/${base}_1.fastq.gz" && -e "$indir/${base}_2.fastq.gz" ]]; then
    pear -f "$indir/${base}_1.fastq.gz" -r "$indir/${base}_2.fastq.gz" -o tmp

    # Reverse "reversed" reads before concatenation
    seqkit seq \
      --validate-seq \
      --reverse \
      --complement \
      --seq-type DNA \
      --threads 1 \
      tmp.unassembled.reverse.fastq > tmp.unassembled.reverse.fastq.bak

    mv tmp.unassembled.reverse.fastq.bak tmp.unassembled.reverse.fastq
    cat tmp.*assembled*.fastq > tmp.fastq
    if [[ -e "$input" ]]; then
      gunzip -c "$input" >> tmp.fastq
    fi
    input=tmp.fastq
  fi

  # Remove adapters, and then remove reads <200 or with phred quality <7
  if [[ -e "$input" ]]; then
    porechop --input "$input" \
      --threads 1 \
      --format fastq | \
        fastp \
          --stdin \
          --stdout \
          --disable_adapter_trimming \
          --thread 1 \
          --length_required 200 \
          --qualified_quality_phred 7 \
          --html tmp.html \
          --json tmp.json | gzip -c > "${base}.fastq.gz"
  fi
  rm -f tmp.*

elif [[ "$platform" == "ION_TORRENT" ]]; then

  # Treat paired-end Ion Torrent data as treated in Oxford Nanopore above
  input="$indir/${base}.fastq.gz"
  if [[ -e "$indir/${base}_1.fastq.gz" && -e "$indir/${base}_2.fastq.gz" ]]; then
    pear -f "$indir/${base}_1.fastq.gz" -r "$indir/${base}_2.fastq.gz" -o tmp
    
    # Reverse "reversed" reads before concatenation
    seqkit seq \
      --validate-seq \
      --reverse \
      --complement \
      --seq-type DNA \
      --threads 1 \
      tmp.unassembled.reverse.fastq > tmp.unassembled.reverse.fastq.bak

    mv tmp.unassembled.reverse.fastq.bak tmp.unassembled.reverse.fastq
    cat tmp.*assembled*.fastq > tmp.fastq
    if [[ -e "$input" ]]; then
      gunzip -c "$input" >> tmp.fastq
    fi
    input=tmp.fastq
  fi

  # Adapter sequences taken from the following Roche url:
  # https://web.archive.org/web/20260204215727/https://rochesequencingstore.com/wp-content/uploads/2017/10/KAPA-Adapter-Kits-for-Ion-Torrent-Systems.pdf
  adapters=">Ion_Torrent_A
CCATCTCATCCCTGCGTGTCTCCGACTCAG
>Ion_Torrent_P1
CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT"
  if [[ -e "$input" ]]; then
    fastp \
      -i "$input" \
      -o "${base}.fastq.gz" \
      --adapter_fasta <(printf '%s\n' "$adapters") \
      --thread 1 \
      --cut_tail \
      --length_required 50 \
      --html tmp.html \
      --json tmp.json
  fi
  rm -f tmp.*

elif [[ "$platform" == "LS454" ]]; then

  # Adapter sequences taken from the following Roche url:
  # https://web.archive.org/web/20240613015329/https://resources.qiagenbioinformatics.com/manuals/biomedicalgenomicsworkbench/500/index.php?manual=Import_Roche_454.html
  adapters=">LS454_FLX
GTTGGAACCGAAAGGGTTTGAATTCAAACCCTTTCGGTTCCAAC
>LS454_Titanium
TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG
>LS454_Linker
CTGCTGTACCGTACATCCGCCTTGGCCGTACAGCAG"

  # Roughly based on the steps outlined in PRINSEQ docs (see 'less stringent'):
  # https://web.archive.org/web/20230601000000*/https://prinseq.sourceforge.net/Preprocessing_454_SFF_chart.pdf
  if [[ -e "$indir/${base}_1.fastq.gz" && -e "$indir/${base}_2.fastq.gz" ]]; then
    fastp \
      --in1 "$indir/${base}_1.fastq.gz" \
      --in2 "$indir/${base}_2.fastq.gz" \
      --out1 "${base}_1.fastq.gz" \
      --out2 "${base}_2.fastq.gz" \
      --thread 1 \
      --adapter_fasta <(printf '%s\n' "$adapters") \
      --cut_front --cut_tail --cut_mean_quality 15 \
      --length_required 60 \
      --html tmp.html \
      --json tmp.json
  fi

  if [[ -e "$indir/${base}.fastq.gz" ]]; then
    fastp \
      -i "$indir/${base}.fastq.gz" \
      -o "${base}.fastq.gz" \
      --adapter_fasta <(printf '%s\n' "$adapters") \
      --cut_front --cut_tail --cut_mean_quality 15 \
      --length_required 60 \
      --html tmp.html \
      --json tmp.json
  fi
  rm -f tmp.*

elif [[ "$platform" == "PACBIO_SMRT" ]]; then

  # Treat paired-end PacBio data as treated in Oxford Nanopore and Ion Torrent
  input="$indir/${base}.fastq.gz"
  if [[ -e "$indir/${base}_1.fastq.gz" && -e "$indir/${base}_2.fastq.gz" ]]; then
    pear -f "$indir/${base}_1.fastq.gz" -r "$indir/${base}_2.fastq.gz" -o tmp

    # Reverse "reversed" reads before concatenation
    seqkit seq \
      --validate-seq \
      --reverse \
      --complement \
      --seq-type DNA \
      --threads 1 \
      tmp.unassembled.reverse.fastq > tmp.unassembled.reverse.fastq.bak

    mv tmp.unassembled.reverse.fastq.bak tmp.unassembled.reverse.fastq
    cat tmp.*assembled*.fastq > tmp.fastq
    if [[ -e "$input" ]]; then
      gunzip -c "$input" >> tmp.fastq
    fi
    input=tmp.fastq
  fi

  # Adapter sequences taken from HiFiAdapterFilt source code:
  # https://github.com/sheinasim-USDA/HiFiAdapterFilt/blob/v3.0.0/DB/pacbio_vectors_db
  adapters=">PacBio_Blunt_Adapter
ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT
>PacBio_C2_Primer
AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA"

  # Length and quality thresholds taken from the following BioStars thread:
  # https://web.archive.org/web/20260206232328/https://www.biostars.org/p/361309/#9616293
  fastp \
    -i "$input" \
    -o "${base}.fastq.gz" \
    --adapter_fasta <(printf '%s\n' "$adapters") \
    --thread 1 \
    --length_required 500 \
    --qualified_quality_phred 10 \
    --html tmp.html \
    --json tmp.json

  rm -f tmp.*

fi

if [ "$(find . -prune -empty)" = "." ]; then
  printf '%s\n' "error: nothing was output, check your input filenames"
  exit 1
fi

# Remove reads mapping to human:
# 1. Kraken2 classifies each read (or pair of reads)
# 2. Grep extracts read IDs matching to Homo sapiens
#   2.5. If paired, separate files for forward and reverse human-matching reads
# 3. Seqkit removes reads matched by grep (both pairs if paired)
if [[ -e "${base}.fastq.gz" ]]; then

  # We use '|| :;' for grep if nothing is matched and ensuring zero exit status
  kraken2 \
    --db "$krakendb" \
    --gzip-compressed \
    --use-names \
    "${base}.fastq.gz" | \
    { grep -F 'Homo sapiens (taxid 9606)' || :; } | \
    cut -f 2 | \
    seqkit grep \
      --threads 1 \
      --invert-match \
      --out-file "${base}.bak.fastq.gz" \
      --pattern-file - \
      "${base}.fastq.gz"
  mv "${base}.bak.fastq.gz" "${base}.fastq.gz"
fi

# At this point, if forward reads exist, so do the reverse ones.
if [[ -e "${base}_1.fastq.gz" ]]; then
  kraken2 \
    --db "$krakendb" \
    --gzip-compressed \
    --use-names \
    --paired \
    "${base}_1.fastq.gz" "${base}_2.fastq.gz" | \
    { grep -F 'Homo sapiens (taxid 9606)' || :; } | \
    cut -f 2 > tmp.human.txt

  # Kraken2 reports forward reads by default, start by removing those
  seqkit grep \
    --threads 1 \
    --invert-match \
    --out-file "${base}_1.bak.fastq.gz" \
    --pattern-file tmp.human.txt \
    "${base}_1.fastq.gz"
  mv "${base}_1.bak.fastq.gz" "${base}_1.fastq.gz"

  # Then rename '.1' to '.2' in pattern file and repeat with reverse file
  sed -i 's/\.1$/.2/' tmp.human.txt
  seqkit grep \
    --threads 1 \
    --invert-match \
    --out-file "${base}_2.bak.fastq.gz" \
    --pattern-file tmp.human.txt \
    "${base}_2.fastq.gz"
  mv "${base}_2.bak.fastq.gz" "${base}_2.fastq.gz"

  rm -f tmp.human.txt
fi

mv -v "$tmp" "${outdir}/${base}"

exit 0
} >&2
