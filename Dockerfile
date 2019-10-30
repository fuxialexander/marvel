FROM nfcore/base:1.7
LABEL authors="Alexander Fu Xi" \
      description="Docker image containing all requirements for MARVEL pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-marvel-1.0dev/bin:$PATH
