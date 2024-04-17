# iChIPv2

!['anyImage?'](./img/logo.png)

This pipeline shows the steps followed in (--paper--) to process iChIPv2 data from Fastq files to peak coordinates and further analysis.
It takes as input Fastq files and reference genomes and organises the data using 2 Sample tables and a config file defined below.

# Installation  
If you don't have it yet, first start by installing Conda [Miniforge](https://github.com/conda-forge/miniforge#miniforge).

We can then move to install tool dependencies with the yaml file found [here](iChIPv2_environment.yml):

```bash
# create environment
conda env create -n iChIPv2 -f iChIPv2_environment.yml

# activate environment
conda activate iChIPv2
```

# Sample tables format
File examples [here](examples/sampleTables/)

## 1. Genome table
Tab (\t) separated file with 2 mandatory columns:  
| Column ID        | Description          | 
| ------------- |:-------------:|
| **sps_id**      | species identifier to be linked to samples as indicated in "Sample table" (e.g. Scer) | 
| **genome_file**     | full path to genome fasta file (NOTE: should be the same file basename used to generate the bowtie2 and fasta indexes, see instructions)   |  

This file is made once you have located and indexed the reference genomes you will use.
  
**NOTE:** All fasta files need to have Bowtie2 and Samtools indexes. How to get them:
```bash
# Bowtie2 index
bowtie2-build SPECIES-REF-GENOME.fasta
# Samtools index (from same folder as file)
samtools faidx SPECIES-REF-GENOME.fasta
```

## 2. Sample table
Tab (\t) separated file with 6 mandatory columns:  
| Column ID        | Description          | 
| ------------- |:-------------:|
| **sample_id**      | unique identifier for the ChIP-seq sample (e.g. Scer_H3K36me3_rep1) | 
| **fastq_1**      | comma-separated list of gzipped fastq files corresponding to read_1 | 
| **fastq_2**     | comma-separated list of gzipped fastq files corresponding to read_2     |  
| **sps_id** | species identifier to be linked to reference genome as indicated in "Genome table" (e.g. Scer)    |  
| **is.control** | whether the sample is to be used as control (e.g. H3 or sequenced input DNA). Possible values: Y/N    |  
| **control** | sample to be used as background control for peak calling. Possible values: sample_id within the same experiment, full path to an existing bam file, empty (in this case, peak calling will be done without background sample option)    |  
  
This file is made once you have located the FASTQ.gz files to analyze.

# Config file format
File example [here](examples/configFile.sh)
This file is formatted as a bash script with all the needed input variables for the pipeline. It will be sourced in all the script to load the variable values, hence, i) don't remove the `#!/usr/bin/env bash` string from the first line, ii) add your own comments (if you ever do it) on lines starting with a `#`, iii) give a proper value to all the defined variables, and always follow [bash syntax](http://www.compciv.org/topics/bash/variables-and-substitution/)

VARIABLE DESCRIPTION TO BE ADDED  


# Running the scripts
## Test run
Once dependencies are installed, make sure the pipeline works by running `pending` with sample data:

```bash
Code to run test
```

## Run it with your own data
### 1. Create the Sample table and config files
Create the Sample table and config files and move them (or not) to a preferred location.

### 2. Run alignment
Locate the config file and execute 00_alignmentLauncher.sh script including the path as an input
```bash
bash 00_alignmentLauncher.sh /PATH/TO/configFile.sh
```
If you chose "noQueue" mode in the config file, you can execute as a job for your clusters queuing system. In this mode, all the files in the Sample table will be aligned one after the other without parallelization. This is the only approach to use if you run the analysis on a computer without a queueing system and the best one to use if you are not very familiar with Linux and/or queuing systems. 

Output files:  
-BAM files --> ${out_dir}/bam_files  
-Bigwig files --> ${out_dir}/bw_files  
-Alignment information files --> ${out_dir}/QC  
-Alignment step files (to restart broken runs) and temporal files -->  ${out_dir}/temp  
  
Note that we store a file named "pipelineStep_${fileLabel}.txt" in "${out_dir}/temp". This file contains one ID for each successful pipeline step and must be deleted if, for any reason, we want to re-start from zero the alignment of a given file or align files with the same sample_id as previously analyzed ones. 

### 3. Run peak calling with MACS2
Locate the config file and execute 01_peakCallingLauncher.sh script including the path as an input
```bash
bash 01_peakCallingLauncher.sh /PATH/TO/configFile.sh
```
This script should only be launched after 00_alignmentLauncher.sh has finished running.

Output files:  
-MACS2 peak files are stored in ${out_dir}/peakCalling/${sps_id}  
-MACS2 Log files --> ${out_dir}/peakCalling/${sps_id}/logs  
-Per file FRiP values --> ${out_dir}/QC/${sps_id}/summaryPeak_*.txt  

### 4. Run the QC
Locate the config file and execute 02_generate_QC.sh script including the path as an input
```bash
bash 02_generate_QC.sh /PATH/TO/configFile.sh
```
This script should only be launched after 01_peakCallingLauncher.sh has finished running.

Output files:  
-Table with all Alignment and Peak Calling stats --> ${out_dir}/QC/gathered_QC.tsv  
-Plots showing CPM density distributions of samples and their respective controls, for each ${sps_id} --> ${out_dir}/QC/${sps_id}_cpmDistrib_1kb.pdf

