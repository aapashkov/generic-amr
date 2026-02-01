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
or BV-BRC identifiers (two integers separated by a dot). Requires 'wget' or
'curl', and 'jq' to be available in PATH.

example: download_genomes.sh -o genomes/ GCA_000005845.2 ERZ3086155 170673.13
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

### download_reads.sh

```text
usage: download_reads.sh [-o OUTDIR] ACCESSION...

Download SRA ACCESSIONs into OUTDIR, which defaults to working directory.
Requires 'fasterq-dump' and 'gzip' to be available in PATH.

example: download_reads.sh -o reads/ SRR33085175 ERR12520663
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
