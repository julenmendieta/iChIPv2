#!/usr/bin/env bash

# Path to Sample table
sample_table="/nfs/users/asebe/cnavarrete/proj/iChIP/fastq_sel/Scer_reads/alignment_samplesheet_v2.tsv"
# Path to genome table
genome_table="/nfs/users/asebe/cnavarrete/proj/iChIP/fastq_sel/Scer_reads/genome_samplesheet.tsv"
# Path to Input FastQ folder. It MUST contain all the FastQ files to be
## aligned (controls can be elsewhere if indicated by path in alignment samplesheet) 
in_dir="/nfs/users/asebe/cnavarrete/proj/iChIP/fastq_sel/Scer_reads/"
# Path to Output folder
out_dir="/nfs/users/asebe/cnavarrete/proj/iChIP/new_pipeline/Scer_output/"
# Path to Git repository
gitP="/nfs/users/asebe/cnavarrete/git_repo/iChIPv2/"
# Path to iChIPv2 Conda env bin folder
## You can find it by loading the environmnet "conda activate iChIPv2"
## and then executing:
## which python | sed "s/\/python//g"
condaEnv="/nfs/users/asebe/cnavarrete/miniconda3/envs/iChIPv2/bin/"
# N CPU for alignment
nCPU=8
# Variable to define how to run the alignment, either by queuing system
## or each sample at a time
## Set it to "queue" (you must have adapted the scripts for it)  or "noQueue"
jobMode="queue"
