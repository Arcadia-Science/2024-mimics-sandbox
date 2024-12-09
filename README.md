# 2024-parasite-human-mimicks

[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)
[![Snakemake](https://img.shields.io/badge/snakemake--green)](https://snakemake.readthedocs.io/en/stable/)

## Purpose

This repository looks for structural protein mimicry between parasites/viruses that infect humans and human proteins.

Molecular mimicry was defined in 1964 as the sharing of antigens between parasite and host (reference: 10.1086/282313).
Pathogens and parasites produce a range of mimics that resemble host components in both form and function:
* The frequency of dinucleotide CpG islands in Influenza virus genomes decreased to levels consistent with the human genome as virus adapted to humans. This makes it less likely that the presence of Influenza will activate host immune responses (reference: 10.1371/journal.ppat.1000079).
* Poxvirus particle packages are covered in host phospholipids so they take on the appearance of apoptotic cellular fragments that are engulfed by neighboring cells, facilitating viral spread (reference: 10.1126/science.1155164).
* Parasitic nematode Ascaris lumbricoides possesses A- and B-like blood group antigens in its polysaccharides (reference: 10.1093/infdis/74.2.81). 

This repository focuses on detecting **structural protein mimics**.
This refers to a parasite structurally mimicking a host protein such that the mimicry confers a fitness advantage to the parasite, by either co-opting or disrupting the function of the mimicked host protein (reference: 10.1038/nrmicro2222; 10.3389/fpara.2023.1162697).
We aim to detect both **perfect** and **imperfect** mimics (reference: 10.1038/nrmicro2222).
* Perfect mimics co-opt host functions to favor pathogen fitness
* Imperfect mimics resemble host components but perform functions that are distinct from those of the host factors they model and are for the benefit of the pathogen.

We outline our approaches for identifying structural mimicry below.

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

## Data

TODO: Add details about the description of input / output data and links to Zenodo depositions, if applicable.

## Overview

### Description of the folder structure

### Description of how the tool works

TODO: add a thorough description of how the tool works and should be used. Consider adding a quickstart guide for users who want to run the pipeline, and/or a demo dataset that they can use to test the pipeline.  

### Compute Specifications

#### `eukaryotic.snakefile`

We ran this snakefile on a MacBook Pro, 2021 (Apple M1 Max chip, 64 GB RAM, 10 cores).
We used up to 6 cores to run the snakefile.
The complete snakefile ran in ~4 hours.
The input and output files take up approximately 60GB of space.

#### `viral.snakefile`

We ran this snakefile on a MacBook Pro, 2021 (Apple M1 Max chip, 64 GB RAM, 10 cores).
We used up to 8 cores to run the snakefile.
The complete snakefile ran in 5 hours with 8 cores (however runtime will depend on download speed as there are a few large downloads).
The input and output files take up approximately 95GB of space.

## Contributing

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide-credit-for-contributions.md).
