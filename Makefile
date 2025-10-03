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

accessions = $(shell basename -as .fna data/genomes/*.fna)
signatures = $(accessions:%=data/signatures/%.sig)
taxonomies = $(signatures:%.sig=%.tax)
annotations = $(accessions:%=data/annotations/%)
resistance = $(accessions:%=data/resistance/%)
plasmids =$(accessions:%=data/plasmids/%)

all: $(signatures) $(plasmids) $(taxonomies) $(resistance) $(annotations)
> @date +'[%F %T] finished' 1>&2

data/databases/platon:
> @date +'[%F %T] downloading platon database' 1>&2 && \
  mkdir -p data/{databases,tmp} && \
  tmp=$$(mktemp -dp data/tmp) && \
  trap 'rm -rf "$$tmp"' EXIT && \
  wget -qO- "$(platon)" | tar -zxC "$$tmp" && \
  mv "$$tmp/db" "$@" && \
  date +'[%F %T] finished downloading platon database' 1>&2

data/databases/gtdb.sbt.zip:
> @date +'[%F %T] downloading gtdb' 1>&2 && \
  mkdir -p data/databases && \
  wget -qO "$@" "$(gtdb)" && \
  date +'[%F %T] finished downloading gtdb' 1>&2

data/signatures/%.sig: data/databases/gtdb.sbt.zip
> @date +'[%F %T] sketching $*' 1>&2 && \
  mkdir -p data/{signatures,logs} && \
  sourmash sketch dna --name "$*" -o "$@" "data/genomes/$*.fna" \
    &> "data/logs/sketch_$*.log" && \
  date +'[%F %T] finished sketching $*' 1>&2

data/signatures/%.tax: data/signatures/%.sig
> @date +'[%F %T] classifying $*' 1>&2 && \
  mkdir -p data/logs && \
  sourmash search --best-only -o - "$<" data/databases/gtdb.sbt.zip \
    2> "data/logs/tax_$*.log" | tail -1 | cut -f4 -d, | cut -f2,3 -d ' ' \
    > "$@" && \
  date +'[%F %T] finished classifying $*' 1>&2

data/resistance/%: data/genomes/%.fna
> @date +'[%F %T] detecting resistance of $*' 1>&2 && \
  mkdir -p data/{tmp,resistance,logs} && \
  tmp=$$(mktemp -dp data/tmp) && \
  trap 'rm -rf "$$tmp"' EXIT && \
  rgi main -i "$<" -o "$$tmp/$*" --include_loose --clean -n 1 \
    &> "data/logs/res_$*.log" && \
  ! grep -q Error "data/logs/res_$*.log" && \
  [[ -f "$$tmp/$*.txt" ]] && [[ -f "$$tmp/$*.json" ]] && \
  mv "$$tmp" "$@" && \
  date +'[%F %T] finished detecting resistance of $*' 1>&2

data/annotations/%: data/signatures/%.tax
> @date +'[%F %T] annotating $*' 1>&2 && \
  mkdir -p data/{tmp,annotations,logs} && \
  tmp=$$(mktemp -dp data/tmp) && \
  trap 'rm -rf "$$tmp"' EXIT && \
  IFS=' ' read -r genus species < "$<" && \
  prokka --outdir "$$tmp" --force --prefix "$*" --locustag GAMR \
    --genus "$$genus" --species "$$species" --strain "$*" --cpus 1 \
    --rfam "data/genomes/$*.fna" \
    &> "data/logs/annotation_$*.log" && \
  mv "$$tmp" "$@" && \
  date +'[%F %T] finished annotating $*' 1>&2

data/plasmids/%: data/databases/platon
> @date +'[%F %T] identifying plasmids in $*' 1>&2 && \
  mkdir -p data/{plasmids,tmp} && \
  tmp=$$(mktemp -dp data/tmp) && \
  trap 'rm -rf "$$tmp"' EXIT && \
  platon -t 1 -d data/databases/platon -p "$*" -o "$$tmp" \
    "data/genomes/$*.fna" > "$$tmp/$*.tsv" && \
  mv "$$tmp" "$@" && \
  date +'[%F %T] finished identifying plasmids in $*' 1>&2
