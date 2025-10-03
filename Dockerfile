FROM debian:13

ARG DEBIAN_FRONTEND="noninteractive"
ARG OLIGOARRAYAUX_URL="https://www.unafold.org/download/oligoarrayaux-3.8-1.x86_64.rpm"
ARG OLIGOARRAYAUX_MD5SUM="6bd4817e1e75e1c6f66ef23ff31830b4"
ARG CARD_URL="https://card.mcmaster.ca/download/0/broadstreet-v4.0.1.tar.bz2"
ARG CARD_MD5SUM="3471be49d19cfee551e52a9a0d82e548"
ENV PYTHONWARNINGS="ignore::SyntaxWarning"
ENV MPLCONFIGDIR="/tmp/matplotlib"
COPY /env/apt.txt /env/pip.txt /tmp/

RUN apt-get update && \
    # Install dependencies from package managers
    xargs -a /tmp/apt.txt apt-get install --no-install-recommends --yes && \
    pip install --no-cache-dir --break-system-packages --no-deps --requirement /tmp/pip.txt && \
    # Setup OligoArrayAux (an RGI dependency)
    wget -O /tmp/oligoarrayaux.rpm --no-hsts "${OLIGOARRAYAUX_URL}" && \
    printf '%s\t/tmp/oligoarrayaux.rpm\n' "${OLIGOARRAYAUX_MD5SUM}" | md5sum -c - && \
    rpm2cpio /tmp/oligoarrayaux.rpm | cpio -idmv -D / && \
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
    # Platon works just fine with MUMmer 3.0.0
    sed -i 's/4,0,0/3,0,0/' /usr/local/lib/python3.13/dist-packages/platon/utils.py && \
    # Clean up
    apt-get clean && \
    rm -rf /tmp/*

WORKDIR /ext
