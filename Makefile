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

c = :
gtdb = https$(c)//farm.cse.ucdavis.edu/~ctbrown/sourmash-db.new/gtdb-rs226/gtdb-reps-rs226-k31.dna.zip
platon = https$(c)//zenodo.org/record/4066768/files/db.tar.gz

accessions = $(shell cat accessions.txt)
genomes = $(accessions:%=data/genomes/%.fna)
annotations = $(accessions:%=data/annotations/%)

all: $(genomes) $(annotations) data/pangenome
> @date +'[%F %T] finished' 1>&2
.PHONY: all

data/databases/platon:
> @date +'[%F %T] downloading platon database' 1>&2 && \
  mkdir -p data/databases && \
  tmp=$$(mktemp -d data/databases/.tmp.XXXXXX) && \
  trap 'rm -rf "$$tmp"' EXIT && \
  wget -qO- "$(platon)" | tar -zxC "$$tmp" && \
  mv "$$tmp/db" "$@" && \
  date +'[%F %T] finished downloading platon database' 1>&2

data/databases/gtdb.sbt.zip:
> @date +'[%F %T] downloading gtdb' 1>&2 && \
  mkdir -p data/databases && \
  wget -qO "$@" "$(gtdb)" && \
  date +'[%F %T] finished downloading gtdb' 1>&2

data/genomes/%.fna: data/databases/platon data/databases/gtdb.sbt.zip
> @sleep $${$$$(c) -1} && \
  date +'[%F %T] downloading $*' 1>&2 && \
  mkdir -p data/{genomes,logs} && \
  ./scripts/download_genomes.sh -o data/genomes '$*' &>> 'data/logs/$*.log' && \
  date +'[%F %T] finished downloading $*' 1>&2

data/annotations/%: data/genomes/%.fna
> @date +'[%F %T] annotating $*' 1>&2 && \
  mkdir -p data/{annotations,logs} && \
  ./scripts/annotate_genome.sh '$<' '$@' data/databases/gtdb.sbt.zip \
    data/databases/platon &>> 'data/logs/$*.log' && \
  date +'[%F %T] finished annotating $*' 1>&2

data/pangenome: $(annotations)
> @date +'[%F %T] building pangenome' 1>&2 && \
  mkdir -p data/logs && \
  ./scripts/tabulate_annotations.py data/annotations data/pangenome &>> data/logs/pangenome.log && \
  date +'[%F %T] finished building pangenome' 1>&2
