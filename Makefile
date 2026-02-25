SHELL := bash
.ONESHELL:
.SHELLFLAGS := -euo pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

ifeq ($(origin .RECIPEPREFIX), undefined)
  $(error This Make does not support .RECIPEPREFIX. Please use GNU Make 4.0 or later)
endif
.RECIPEPREFIX = >

# URLs for external resources
c = :
gtdb = https$(c)//farm.cse.ucdavis.edu/~ctbrown/sourmash-db.new/gtdb-rs226/gtdb-reps-rs226-k31.dna.zip
platon = https$(c)//zenodo.org/record/4066768/files/db.tar.gz

# Retrieve accessions, removing duplicate entries
accessions = $(shell cut -f1 accessions.txt | awk '!seen[$$0]++')

# Each accession will require a genome and an annotation
genomes = $(accessions:%=data/genomes/%.fna)
annotations = $(accessions:%=data/annotations/%)

# Explicitly ask make to keep read files (because they are intermediate)
raw_reads = $(accessions:%=data/reads/raw/%)
trimmed_reads = $(accessions:%=data/reads/trimmed/%)
.SECONDARY: $(raw_reads) $(trimmed_reads)

# Genomes with specific prefixes will always be downloaded and not assembled
prefixes = GCA_ GCF_ ERZ BVBRC_

# Function that looks up sequencing technology of an accession (if provided)
# Usage: $(call get_tech,accession)
define get_tech
$(strip $(shell awk -F'\t' -v key='$(1)' '$$1 == key { print ($$2 ? $$2 : ""); exit }' 'accessions.txt'))
endef

# Run the entire pipeline (default when running simply 'make')
all: $(genomes) $(annotations) data/pangenome
> @date +'[%F %T] finished' 1>&2
.PHONY: all

# Downloads platon database
data/databases/platon:
> @date +'[%F %T] downloading platon database' 1>&2 && \
  mkdir -p data/databases && \
  tmp=$$(mktemp -d data/databases/.tmp.XXXXXX) && \
  trap 'rm -rf "$$tmp"' EXIT && \
  wget -qO- "$(platon)" | tar -zxC "$$tmp" && \
  mv "$$tmp/db" "$@" && \
  date +'[%F %T] finished downloading platon database' 1>&2

# Downloads Sourmash GTDB reference
data/databases/gtdb.sbt.zip:
> @date +'[%F %T] downloading gtdb' 1>&2 && \
  mkdir -p data/databases && \
  wget -qO "$@" "$(gtdb)" && \
  date +'[%F %T] finished downloading gtdb' 1>&2

# Downloads human kraken database
data/databases/kraken:
> @date +'[%F %T] building kraken human database' 1>&2 && \
  mkdir -p data/{databases,logs} && \
  tmp=$$(mktemp -d data/databases/.tmp.XXXXXX) && \
  trap 'rm -rf "$$tmp"' EXIT && \
  ( kraken2-build --download-taxonomy --use-ftp --db "$$tmp" && \
  kraken2-build --download-library human --use-ftp --no-masking --db "$$tmp" && \
  kraken2-build --build --db "$$tmp" && \
  kraken2-build --clean --db "$$tmp" ) &>> data/logs/kraken-build.log && \
  touch "$$tmp/.gamr-validated" && \
  mv "$$tmp" "$@" && \
  date +'[%F %T] finished building kraken human database' 1>&2

# Define a rule for each prefix referencing genomes that have to be downloaded
define download_template
data/genomes/$(1)%.fna:
> @date +'[%F %T] downloading $(1)$$*' 1>&2 && \
  mkdir -p data/{genomes,logs} && \
  ./scripts/download_genomes.sh -o data/genomes '$(1)$$*' \
    &>> 'data/logs/$(1)$$*.log' && \
  date +'[%F %T] finished downloading $(1)$$*' 1>&2
endef
$(foreach p,$(prefixes),$(eval $(call download_template,$(p))))

# Genomes that are not downloaded must be assembled from reads
data/genomes/%.fna: data/reads/trimmed/%
> @date +'[%F %T] assembling $* (NOT IMPLEMENTED)' 1>&2 && \
  date +'[%F %T] finished assembling $* (NOT IMPLEMENTED)' 1>&2

# Downloads raw reads from SRA
data/reads/raw/%:
> @date +'[%F %T] downloading $*' 1>&2 && \
  mkdir -p data/{reads/raw,logs} && \
  ./scripts/download_reads.sh -o data/reads/raw '$*' &>> 'data/logs/$*.log' && \
  date +'[%F %T] finished downloading $*' 1>&2

# Trims reads
data/reads/trimmed/%: data/reads/raw/% data/databases/kraken
> @date +'[%F %T] trimming $*' 1>&2 && \
  mkdir -p data/{reads/trimmed,logs} && \
  ./scripts/trim_reads.sh '$<' data/reads/trimmed data/databases/kraken \
    $(call get_tech,$*) &>> 'data/logs/$*.log' && \
  date +'[%F %T] finished trimming $*' 1>&2

# Annotates genomes
data/annotations/%: data/genomes/%.fna data/databases/platon data/databases/gtdb.sbt.zip
> @date +'[%F %T] annotating $*' 1>&2 && \
  mkdir -p data/{annotations,logs} && \
  ./scripts/annotate_genome.sh '$<' '$@' data/databases/gtdb.sbt.zip \
    data/databases/platon &>> 'data/logs/$*.log' && \
  date +'[%F %T] finished annotating $*' 1>&2

# Builds final output table
data/pangenome: $(annotations)
> @date +'[%F %T] building pangenome' 1>&2 && \
  mkdir -p data/logs && \
  ./scripts/tabulate_annotations.py data/annotations data/pangenome &>> data/logs/pangenome.log && \
  date +'[%F %T] finished building pangenome' 1>&2
