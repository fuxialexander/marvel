FROM continuumio/miniconda3:latest
LABEL authors="Alexander Fu Xi" \
      description="Docker image containing all requirements for MARVEL pipeline"

RUN apt-get update && apt-get install -y procps && apt-get clean -y

COPY environment.yml /
COPY glmpath_0.98.tar.gz /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-marvel-1.1dev/bin:$PATH
RUN R -e "install.packages(c('glmpath'), dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R CMD INSTALL glmpath_0.98.tar.gz


RUN set -e \
      && ln -sf /bin/bash /bin/sh

RUN set -e \
      && apt-get -y update \
      && apt-get install --yes build-essential \
      && apt-get -y dist-upgrade \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        autoconf ca-certificates curl gcc libbz2-dev libcurl4-gnutls-dev \
        libgsl-dev libperl-dev liblzma-dev libssl-dev libz-dev make \
        libncurses5-dev zlib1g-dev \
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*

RUN set -e \
      && wget https://raw.githubusercontent.com/dceoy/print-github-tags/master/print-github-tags \
      && mv print-github-tags /usr/local/bin/print-github-tags
RUN set -eo pipefail \
      && chmod +x /usr/local/bin/print-github-tags \
      && print-github-tags --release --latest --tar samtools/htslib \
        | xargs -i curl -SL {} -o /tmp/htslib.tar.gz \
      && tar xvf /tmp/htslib.tar.gz -C /usr/local/src --remove-files \
      && mv /usr/local/src/htslib-* /usr/local/src/htslib \
      && cd /usr/local/src/htslib \
      && autoheader \
      && autoconf \
      && ./configure \
      && make \
      && make install

RUN set -eo pipefail \
      && print-github-tags --release --latest --tar samtools/bcftools \
        | xargs -i curl -SL {} -o /tmp/bcftools.tar.gz \
      && tar xvf /tmp/bcftools.tar.gz -C /usr/local/src --remove-files \
      && mv /usr/local/src/bcftools-* /usr/local/src/bcftools \
      && cd /usr/local/src/bcftools \
      && autoheader \
      && autoconf \
      && ./configure --enable-libgsl --enable-perl-filters \
      && make \
      && make install

RUN set -eo pipefail \
      && print-github-tags --release --latest --tar samtools/samtools \
        | xargs -i curl -SL {} -o /tmp/samtools.tar.gz \
      && tar xvf /tmp/samtools.tar.gz -C /usr/local/src --remove-files \
      && mv /usr/local/src/samtools-* /usr/local/src/samtools \
      && cd /usr/local/src/samtools \
      && autoheader \
      && autoconf \
      && ./configure \
      && make \
      && make install