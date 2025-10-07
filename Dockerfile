FROM debian:trixie-20250929

ARG DEBIAN_FRONTEND="noninteractive"
ARG OLIGOARRAYAUX_URL="https://anaconda.org/bioconda/oligoarrayaux/3.8.1/download/linux-64/oligoarrayaux-3.8.1-pl5321h9948957_0.tar.bz2"
ARG OLIGOARRAYAUX_MD5SUM="97c83c978d07d0a01ec2f96cfc956688"
ARG CARD_URL="https://card.mcmaster.ca/download/0/broadstreet-v4.0.1.tar.bz2"
ARG CARD_MD5SUM="3471be49d19cfee551e52a9a0d82e548"
ARG MUMMER_URL="https://anaconda.org/bioconda/mummer4/4.0.1/download/linux-64/mummer4-4.0.1-pl5321h9948957_0.tar.bz2"
ARG MUMMER_MD5SUM="69c9b49ed87fcf4142b6ca59ab8c4ba0"
ENV PYTHONWARNINGS="ignore::SyntaxWarning"
ENV MPLCONFIGDIR="/tmp/matplotlib"
COPY /env/apt.txt /env/pip.txt /tmp/

RUN apt-get update && \
    # Install dependencies from package managers
    xargs -a /tmp/apt.txt apt-get install --no-install-recommends --yes && \
    pip install --no-cache-dir --break-system-packages --no-deps --requirement /tmp/pip.txt && \
    # Setup OligoArrayAux (an RGI dependency)
    wget -O /tmp/oligoarrayaux.tar.bz2 --no-hsts "${OLIGOARRAYAUX_URL}" && \
    printf '%s\t/tmp/oligoarrayaux.tar.bz2\n' "${OLIGOARRAYAUX_MD5SUM}" | md5sum -c - && \
    tar -C /usr/local -jxf /tmp/oligoarrayaux.tar.bz2 && \
    # Setup MUMmer (a Platon dependency)
    wget -O /tmp/mummer.tar.bz2 --no-hsts "${MUMMER_URL}" && \
    printf '%s\t/tmp/mummer.tar.bz2\n' "${MUMMER_MD5SUM}" | md5sum -c - && \
    tar -C /usr/local -jxf /tmp/mummer.tar.bz2 && \
    # Setup Prokka database (symlink required for nonroot users)
    prokka --dbdir /usr/share/prokka/db --setupdb && \
    mkdir -p /.local/lib/prokka && \
    ln -s /usr/share/prokka/db /.local/lib/prokka/db && \
    # Setup RGI database (sample workload required to build full database)
    wget -O /tmp/card.tar.bz2 --no-hsts "${CARD_URL}" && \
    printf '%s\t/tmp/card.tar.bz2\n' "${CARD_MD5SUM}" | md5sum -c - && \
    tar -C /tmp -jxf /tmp/card.tar.bz2 ./card.json && \
    rgi load -i /tmp/card.json && \
    rgi main -t protein -n 1 -i /usr/share/doc/proteinortho/examples/E.faa -o /tmp/result --clean && \
    # Clean up
    apt-get clean && \
    rm -rf /tmp/*

WORKDIR /ext
