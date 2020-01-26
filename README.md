# ![marvel](docs/images/marvel_banner.png)

**Multigranular Analysis of Regulatory Variants on the Epigenomic Landscape**.

[![Build Status](https://travis-ci.com/fuxialexander/marvel.svg?branch=master)](https://travis-ci.com/fuxialexander/marvel)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/marvel.svg)](https://hub.docker.com/r/nfcore/marvel)

## Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

i. Install [`nextflow`](https://nf-co.re/usage/installation)

ii. Install one of [`docker`](https://docs.docker.com/engine/installation/), [`singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`conda`](https://conda.io/miniconda.html)

iii. Download the pipeline and test it on a minimal dataset with a single command

```bash
nextflow run marvel -profile test,<docker/singularity/conda>
```

iv. Start running your own analysis!

<!-- TODO nf-core: Update the default command above used to run the pipeline -->
```bash
nextflow run marvel -profile <docker/singularity/conda> --enhancer_bed "enhancer.bed" --enhancer_bed "promoter.bed" --genome hg19
```

See [usage docs](docs/usage.md) for all of the available options when running the pipeline.

## Documentation

The marvel pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
    * [Reference genomes](https://nf-co.re/usage/reference_genomes)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)


## Credits

MARVEL is implemented using a boilplate created by the nf-core team (https://nf-co.re/).

## Citation

<!-- TODO: to be updated -->
If you use MARVEL for your analysis, please cite it as: 
Alexander Xi Fu, Kathy Nga-Chu Lui, Clara Sze-Man Tang, Ray Kit Ng, Frank Pui-Ling Lai, Sin-Ting
Lau, Zhixin Li, Maria-Mercè Gracia-Barcelo, Pak-Chung Sham, Paul Kwang-Hang Tam, Elly Sau-Wai
Ngan and Kevin Y. Yip. Whole-genome analysis of noncoding genetic variations identifies
multigranular regulatory element perturbations associated with Hirschsprung disease. (2020)
(DOI AND JOURNAL TO BE ADDED).