#!/bin/bash
# -*- ENCODING: UTF-8 -*-

#########################   USER-DEFINED INPUT   ###############################

# Input Paths
## Path to Peak Calling samplesheet
peakS="/home/jmendietaes/data/2021/Arnau/pipelineTest/sampleSheets/peakCalling_samplesheet.csv"
## Path to genome samplehseet
genomeS="/home/jmendietaes/data/2021/Arnau/pipelineTest/sampleSheets/genome_samplesheet.csv"
## Path to Output folder
out_dir="/home/jmendietaes/data/2021/Arnau/pipelineTest/output2" 
## Path to Git repository
gitP="/home/jmendietaes/programas/iChIPv2"

# Input parameters
## No multithread in peak calling step with MACS2
#nCPU=1  
## Conda environmnet name to be loaded. Path to env bin folder
# You can find it by loading the environmnet "conda activate iChIPv2"
# and then executing:
# which python | sed "s/\/python//g"
condaEnv="/home/jmendietaes/programas/miniconda3/envs/iChIPv2/bin"
## Variable to define how to run the alignment, either by queuing sustem
## or each sample at a time
## Set it to "Slurm", "Torque", or "noQueue"
jobMode="Slurm"

#################################  HOW TO RUN ME  ##############################
#bash /home/jmendietaes/programas/iChIPv2/scripts/01_peakCallingLauncher.sh

##################################   CODE   ####################################
# Get number of files to align (header row has to be removed)
nJobs=$(wc -l ${peakS} | awk '{print $1'})
nJobs=$((${nJobs} - 1))

# Purge current modules before starting job (comment it if no module system is
# installed)
module purge

# Load conda environment
#conda activate ${condaEnv}
export PATH="${condaEnv}:$PATH"


# Remove from samplesheets newline characters that might come from working 
# with excel
sed -i 's/\r$//' ${peakS}

if [[ ${jobMode} == "Slurm" ]]; then 
    
    sbatch --array=1-${nJobs} --job-name=iChIP-peakCall --cpus-per-task=1 \
        --mem=30G --time=00-08:00:00 -p short \
        ${gitP}/scripts/sub-scripts/01a_peakCalling.sh \
        ${peakS} ${genomeS} ${out_dir} ${jobMode}

elif [[ ${jobMode} == "Torque" ]]; then
    echo "Pending to written"

elif [[ ${jobMode} == "noQueue" ]]; then 
    for nj in $(seq 1 ${nJobs}); do
        bash ${gitP}/scripts/sub-scripts/01a_peakCalling.sh \
            ${peakS} ${genomeS} ${out_dir} ${nj}
    done

else
    echo "Wrong jobMode!"
    # "command not found"
    exit 127
fi


exit 0