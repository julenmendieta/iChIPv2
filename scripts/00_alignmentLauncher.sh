#!/bin/bash
# -*- ENCODING: UTF-8 -*-

#########################   USER-DEFINED INPUT   ###############################

# Input Paths
## Path to alignment samplesheet
alignS="/nfs/users/asebe/cnavarrete/proj/iChIP/fastq_sel/Nvec_reads/alignment_samplesheet_test.csv"
## Path to genome samplehseet
genomeS="/nfs/users/asebe/cnavarrete/proj/iChIP/new_pipeline/genome_samplesheet.csv"
## Path to Output folder
out_dir="/nfs/users/asebe/cnavarrete/proj/iChIP/new_pipeline/Nvec_output/"
## Path to Git repository
gitP="/nfs/users/asebe/cnavarrete/git_repo/iChIPv2/"

# Input parameters
## Number of threads for each job
nCPU=8  
## Conda environmnet name to be loaded. Path to env bin folder
# You can find it by loading the environmnet "conda activate iChIPv2"
# and then executing:
# which python | sed "s/\/python//g"
condaEnv="/nfs/users/asebe/cnavarrete/miniconda3/envs/iChIPv2/bin/"
## Variable to define how to run the alignment, either by queuing sustem
## or each sample at a time
## Set it to "Slurm", "SGE", or "noQueue"
jobMode="SGE"

#################################  HOW TO RUN ME  ##############################
#bash /home/jmendietaes/programas/iChIPv2/scripts/00_alignmentLauncher.sh

##################################   CODE   ####################################
# Get number of files to align (header row has to be removed)
nJobs=$(wc -l ${alignS} | awk '{print $1'})
nJobs=$((${nJobs} - 1))

# Purge current modules before starting job (comment it if no module system is
# installed)
module purge

# Load conda environment
#conda activate ${condaEnv}
export PATH="${condaEnv}:$PATH"


# Remove from samplesheets newline characters that might come from working 
# with excel
sed -i 's/\r$//' ${alignS}

if [[ ${jobMode} == "Slurm" ]]; then 
    
    sbatch --array=1-${nJobs} --job-name=iChIP-align --cpus-per-task=${nCPU} \
        --mem=30G --time=00-08:00:00 -p short \
        ${gitP}/scripts/sub-scripts/00a_alignment.sh \
        ${alignS} ${genomeS} ${out_dir} ${nCPU} ${jobMode} ${condaEnv}

elif [[ ${jobMode} == "SGE" ]]; then

    qsub -t 1-${nJobs} -N iChIP-align \
         -q long-centos79,mem_512 -pe smp ${nCPU}\
         ${gitP}/scripts/sub-scripts/00a_alignment.sh \
         ${alignS} ${genomeS} ${out_dir} ${nCPU} ${jobMode} ${condaEnv}
         
         
elif [[ ${jobMode} == "noQueue" ]]; then 
    for nj in $(seq 1 ${nJobs}); do
        bash ${gitP}/scripts/sub-scripts/00a_alignment.sh \
            ${alignS} ${genomeS} ${out_dir} ${nCPU} ${nj} ${condaEnv}
    done

else
    echo "Wrong jobMode!"
    # "command not found"
    exit 127
fi


exit 0