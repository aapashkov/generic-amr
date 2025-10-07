# Generic AMR: A general purpose pipeline for AMR analyses

> [!CAUTION]
> **This pipeline is still experimental**. Use at your own risk.

**Generic AMR** is a pipeline that extracts features from bacterial genomes
useful for antimicrobial resistance (AMR) prediction.

### Quickstart

1. Install [Docker Engine or Desktop](https://docs.docker.com/engine/install/)
and [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) on your
system.
2. Open a terminal, and run the following commands:

```shell
: # 🌎 Clone this repository
git clone https://github.com/aapashkov/generic-amr

: # 📌 Change into base directory
cd generic-amr

: # 📄 Place your genome accessions in accessions.txt

: # ⚡ Run pipeline (you might need to run it as docker-compose on some systems)
docker compose run --rm pipeline
```

### Configuration

The `accessions.txt` file lists, one per line, the genome accessions to process
and includes some example genome accessions to get started. These genome
accessions may be any of the following:

- **[NCBI Genome/RefSeq](https://ncbi.nlm.nih.gov/datasets/genome) identifiers**. They
start with GCA_ or GCF_. For example,
[GCA_000005845.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000005845.2).
- **[ENA](https://www.ebi.ac.uk/ena/browser/home) sequence assembly analysis
identifiers**. They start with ERZ. For example,
[ERZ3086155](https://www.ebi.ac.uk/ena/browser/view/ERZ3086155).
- **[BV-BRC](https://www.bv-brc.org/searches/GenomeSearch) identifiers**. They
consist of two integers separated by a dot. For example,
[170673.13](https://www.bv-brc.org/view/Genome/170673.13).
- **Private genomes**. When including these, place them in the `data/genomes`
directory, decompressed and ending with `.fna`. Then, add their basenames (i.e.,
without the `.fna` extension) to `accessions.txt`.

Create a `.env` file in the root directory of this repository storing pipeline
configuration details. Currently, the following environment variables are
supported:

- `UID` (int) - your user ID, defaults to root if not specified.
- `GID` (int) - your group ID, defaults to root if not specified.
- `JOBS` (int) - number of CPUs to use for parallelization, defaults to all
available CPUs.

If you wish to connect to the Docker container to perform manual operations
instead of running the pipeline, use the following command:

```shell
: # 📦 Connect to Docker container without running the pipeline
docker compose run --rm pipeline bash
```

### Output description

All pipeline outputs are stored in the `data` directory. You should expect a
file structure like the following (where `{accession}` is a single entry in
`accessions.txt`):

```shell
├── 📁 data/                # Pipeline outputs
│   ├── 📁 annotations/     # Genome annotations
│   │   └── 📁 {accession}/
│   │       ├── 📄 {accession}.platon.*     # Platon plasmid predictions
│   │       ├── 📄 {accession}.prokka.*     # Prokka annotations
│   │       ├── 📄 {accession}.rgi.*        # RGI annotations
│   │       └── 📄 {accession}.sourmash.*   # Sourmash signatures and taxonomies
│   ├── 📁 databases/       # Reference databases
│   │   ├── 📁 platon/      # Platon database
│   │   └── 📄 gtdb.sbt.zip # Sourmash database with genomes from GTDB RS226
│   ├── 📁 genomes/         # Decompressed genomes with .fna extension
│   │   └── 📄 {accession}.fna
│   └── 📁 logs/            # Pipeline execution logs, one per accession
├── 📁 env/                 # Environment definitions
│   ├── 📄 apt.txt          # Debian's apt package manager dependencies
│   └── 📄 pip.txt          # Python's pip requirements file
├── 📄 compose.yml          # Docker Compose execution parameters
├── 📄 Dockerfile           # Docker build definition
├── 📄 LICENCE              # MIT license
├── 📄 Makefile             # Pipeline steps
└── 📄 README.md            # Documentation
```
