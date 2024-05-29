#!/bin/bash
# -*- ENCODING: UTF-8 -*-


#########################   USER-DEFINED INPUT   ###############################
# Path to config file
configFile=$1
#configFile=/nfs/users/asebe/cnavarrete/proj/iChIP/fastq_sel/Scer_reads/configFile.sh

#################################  HOW TO RUN ME  ##############################
#bash /nfs/users/asebe/cnavarrete/git_repo/iChIPv2/scripts/02_QC-launcher.sh /nfs/users/asebe/cnavarrete/proj/iChIP/fastq_sel/Scer_reads/configFile.txt




#############   LOAD USER-DEFINED VARIABLES IN CONFIG   ########################
# make sure the user didn't add a new non-commented line
nVariable=$(grep -vE '^(\s*$|#)' ${configFile} | wc -l)
if [[ ${nVariable} == 8 ]] ; then 
    echo "Config File has correct number of variables/non-commented lines"
else
    echo -e "\nERROR: Config File format is WRONG";
    echo -e "Maximum number of allowed variables/non-commented lines is ${nVariable}, when should be 8\n";
    exit 1;
fi

source ${configFile}
##################################   CODE   ####################################

# Purge current modules before starting job (comment it if no module system is
# installed)
module purge

# Load conda environment
#conda activate ${condaEnv}
export PATH="${condaEnv}:$PATH"

# Remove from Sample table newline characters that might come from working 
# with excel
sed -i 's/\r$//' ${sample_table}

if [[ ${jobMode} == "queue" ]]; then 
    
    # Slurm example
    # sbatch --array=1 --job-name=iChIP-QC --cpus-per-task=${nCPU} \
    #     --mem=30G --time=00-08:00:00 -p short \
    #     ${gitP}/scripts/sub-scripts/02a_generate_QC.sh \
    #     ${configFile}

    # SGE example
    qsub -t ${nCPU} -N iChIP-QC \
         -q long-centos79,mem_512 \
         ${gitP}/scripts/sub-scripts/02a_generate_QC.sh \
         ${configFile}

elif [[ ${jobMode} == "noQueue" ]]; then 
    for nj in $(seq 1 ${nJobs}); do
        bash ${gitP}/scripts/sub-scripts/02a_generate_QC.sh \
            ${configFile}
    done

else
    echo "Wrong jobMode!"
    # "command not found"
    exit 127
fi


exit 0