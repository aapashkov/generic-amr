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

: # 📄 Place your genomes with .fna extension in data/genomes directory

: # ⚡ Run pipeline
docker compose run --rm pipeline
```

### Configuration

Create a `.env` file in the root directory of this repository storing pipeline
configuration details. Currently, the following environment variables are
supported:

- `UID` (int) - your user ID, defaults to root if not specified.
- `GID` (int) - your group ID, defaults to root if not specified.
- `JOBS` (int) - number of CPUs to use for parallelization, defaults to all
available CPUs.

Genomes are read from `data/genomes`. They must be decompressed and their
filenames must end with `.ext`.

If you wish to connect to the Docker container to perform manual operations
instead of running the pipeline, use the following command:

```shell
: # 📦 Connect to Docker container without running the pipeline
docker compose run --rm pipeline bash
```

### Output description

All pipeline outputs are stored in the `data` directory. You should expect a
file structure like the following:

```shell
├── 📁 data/                # Pipeline outputs
│   ├── 📁 annotations/     # Prokka annotations
│   ├── 📁 databases/       # Reference databases
│   │   └── 📄 gtdb.sbt.zip # Sourmash database with bacterial and archaeal genomes from GTDB RS226
│   ├── 📁 genomes/         # Input decompressed genomes with .fna
│   ├── 📁 logs/            # Pipeline execution logs
│   ├── 📁 plasmids/        # Plasmid detection results
│   ├── 📁 resistance/      # RGI predictions from CARD reference
│   ├── 📁 signatures/      # Sourmash signature files and taxonomical predictions
│   │   ├── 📄 *.sig        # Sourmash signature files
│   │   └── 📄 *.tax        # Predicted species name from GTDB RS226 using Sourmash
│   └── 📁 tmp/             # Temporary directory location
├── 📁 env/                 # Environment definitions
│   ├── 📄 apt.txt          # Debian's apt package manager dependencies
│   └── 📄 pip.txt          # Python's pip requirements file
├── 📄 compose.yml          # Docker execution parameters
├── 📄 Dockerfile           # Docker build definition
├── 📄 LICENCE              # MIT license
├── 📄 Makefile             # Pipeline steps
└── 📄 README.md            # Documentation
```
