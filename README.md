# iChIPv2

!['anyImage?'](./img/logo.png)

This pipeline shows the steps followed in (--paper--) to process iChIPv2 data from Fastq files to peak coordinates and further analysis.
It takes as input Fastq files and reference genomes and organises the data by the use of 3 Samplesheets: i) reference genomes, ii) fastq files info, and iii) peak calling info.

# Installation  
If you don't have it yet, first start by installing Conda [Miniforge](https://github.com/conda-forge/miniforge#miniforge).

We can then move to install tool dependencies by the yaml file found [here](iChIPv2_environment.yaml):

```bash
# create environment
conda env create -n iChIPv2 -f iChIPv2_environment.yaml

# activate environment
conda activate iChIPv2
```

# Test run
Once dependencies are installed, make sure the pipeline works by running `pending` with sample data:

```bash
Code to run test
```
