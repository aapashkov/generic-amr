# Generic AMR Script Documentation

### annotate_genome.sh

```text
usage: annotate_genome.sh INPUT OUTDIR SBT PLATON

Annotate an INPUT fasta genome and save results to OUTDIR. Must provide a
Sourmash GTDB SBT zip file for taxonomic classification, which you can find
here: https://sourmash.readthedocs.io/en/latest/databases.html. The PLATON
database must also be provided for plasmid detection. The RGI database must be
loaded for global usage. Requires 'prokka', 'rgi', 'barrnap', 'cmscan',
'platon', and 'sourmash' to be available in PATH. 

example: annotate_genome.sh genome.fna annotation/ gtdb.sbt.zip
```

### download_genomes.sh

```text
usage: download_genomes.sh [-o OUTDIR] ACCESSION...

Download genome ACCESSIONs into OUTDIR, which defaults to working directory.
Valid ACCESSION strings include NCBI Genome/RefSeq identifiers (starting with
GCA_ or GCF_), ENA sequence assembly analysis identifiers (starting with ERZ),
or BV-BRC identifiers (starting with BVBRC_). Requires 'wget' or 'curl', and
'jq' to be available in PATH.

example: download_genomes.sh -o genomes/ GCA_000005845.2 ERZ3086155 BVBRC_170673.13
```

### download_reads.sh

```text
usage: download_reads.sh [-o OUTDIR] ACCESSION...

Download SRA ACCESSIONs into OUTDIR, which defaults to working directory.
Requires 'fasterq-dump' and 'gzip' to be available in PATH.

example: download_reads.sh -o reads/ SRR33085175 ERR12520663
```

### get_biosample.sh

```text
usage: get_biosample.sh ACCESSION...

Retrieve the BioSamples of genome ACCESSIONs and print them to stdout in TSV.
Valid ACCESSION strings include NCBI Genome/RefSeq identifiers (starting with
GCA_ or GCF_), ENA sequence assembly analysis identifiers (starting with ERZ),
or BV-BRC identifiers (two integers separated by a dot). Requires 'wget' or
'curl', and 'jq' to be available in PATH.

example: get_biosample.sh GCA_000005845.2 ERZ3086155 170673.13
```

### tabulate_annotations.py

```text
usage: tabulate_annotations.py INDIR OUTDIR

Create a pangenome-like table from the annotations created using the
annotate_genome.sh script in INDIR and save results into OUTDIR. If the JOBS
environment variable is set, will use JOBS processes, and four otherwise.
Requires mmseqs, proteinortho, mafft and snp-sites to be available in PATH.

example: JOBS=10 tabulate_annotations.py annotations/ pangenome/
```

### trim_reads.sh

```text
usage: trim_reads.sh INDIR OUTDIR KRAKENDB [PLATFORM]

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

example: trim_reads.sh raw/SRR33085175/ trimmed/ OXFORD_NANOPORE
```
