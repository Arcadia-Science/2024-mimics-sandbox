# 2024-mimics-sandbox

[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)
[![Snakemake](https://img.shields.io/badge/snakemake--green)](https://snakemake.readthedocs.io/en/stable/)

## Purpose

This repository looks for structural protein mimicry between parasites/viruses that infect specific hosts.
This repository is an idea sandbox.
The scripts aren't guaranteed to work, to have complete documentation, or to achieve their intended purpose.
For a complete and functional code base, please see the repository (Arcadia-Science/2024-mimic-benchmarking)[https://github.com/Arcadia-Science/2024-mimic-benchmarking/tree/v1.1] associated with the pub "[A method for computational discovery of viral structural mimics](https://doi.org/10.57844/arcadia-1eu9-gcsx)."

## Installation and Setup

This repository uses Snakemake to run the pipeline and conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/projects/miniconda/en/latest/). After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the pipeline run environment.

```{bash}
mamba env create -n mimicks --file envs/dev.yml
conda activate mimicks
```

Snakemake manages rule-specific environments via the `conda` directive and using environment files in the [envs/](./envs/) directory. Snakemake itself is installed in the main development conda environment as specified in the [dev.yml](./envs/dev.yml) file.

To start the pipeline, run:

```{bash}
snakemake -s <snakefile> --software-deployment-method conda -j 8
```

### Compute Specifications

#### `eukaryotic.snakefile`

We ran this snakefile on a MacBook Pro, 2021 (Apple M1 Max chip, 64 GB RAM, 10 cores).
We used up to 6 cores to run the snakefile.
The complete snakefile ran in ~4 hours.
The input and output files take up approximately 60GB of space.

#### `viral.snakefile`

We ran this snakefile on a MacBook Pro, 2021 (Apple M1 Max chip, 64 GB RAM, 10 cores).
We used up to 8 cores to run the snakefile.
The complete snakefile ran in 50 minutes with 8 cores (however runtime will depend on download speed as there are a few large downloads).
The input and output files take up approximately 95GB of space.

## Contributing

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide-credit-for-contributions.md).
