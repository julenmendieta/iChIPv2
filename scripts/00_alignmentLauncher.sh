#!/bin/bash
# -*- ENCODING: UTF-8 -*-

##==============================================================================
## SLURM VARIABLES
#SBATCH --job-name=Chip_fqToBw
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=00-08:00:00
#SBATCH -p short
#SBATCH -o /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.out  
#SBATCH -e /home/jmendietaes/jobsSlurm/outErr/%x_%A_%a.err 
##SBATCH --dependency=afterany:624523

#########################   USER-DEFINED INPUT   ###############################

# Input Paths
## Path to alignment samplesheet
alignS="/home/jmendietaes/data/2021/Arnau/pipelineTest/sampleSheets/alignment_samplesheet.csv"
## Path to genome samplehseet
genomeS="/home/jmendietaes/data/2021/Arnau/pipelineTest/sampleSheets/genome_samplesheet.csv"
## Path to Output folder
out_dir="/home/jmendietaes/data/2021/Arnau/pipelineTest/output" 
## Path to Git repository
gitP="/home/jmendietaes/programas/iChIPv2"

# Input parameters
## Number of threads for each job
nCPU=8  
## Conda environmnet name to be loaded
condaEnv="iChIPv2"
## Variable to define how to run the alignment, either by queuing sustem
## or each sample at a time
## Set it to "Slurm", "Torque", or "noQueue"
jobMode="Slurm"

##################################   CODE   ####################################
# Get number of files to align (header is index zero, that wont check)
nJobs=$(wc -l ${alignS} | awk '{print $1'})

# Purge current modules before starting job (comment it if no module system is
# installed)
module purge

# Load conda environment
conda activate ${condaEnv}

# Remove from samplesheets newline characters that might come from working 
# with excel
sed -i 's/\r$//' ${alignS}

if [[ ${jobMode} == "Slurm" ]]; then 
    
    sbatch --array=1-${N} --job-name=iChIP-align --cpus-per-task=${nCPU} \
        --mem=30G --time=00-08:00:00 \
        ${gitP}/scripts/sub-scripts/00a_alignment.sh \
        ${alignS} ${genomeS} ${out_dir} ${nCPU} ${jobMode}

elif [[ ${jobMode} == "Torque" ]]; then
    echo "Pending to written"

elif [[ ${jobMode} == "noQueue" ]]; then 
    for nj in $(seq 1 ${nJobs}); do
        bash ${gitP}/scripts/sub-scripts/00a_alignment.sh \
            ${alignS} ${genomeS} ${out_dir} ${nCPU} ${nj}
    done

else
    echo "Wrong jobMode!"
    # "command not found"
    exit 127
fi


exit 0