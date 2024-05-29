#!/bin/bash
# -*- ENCODING: UTF-8 -*-


#########################   USER-DEFINED INPUT   ###############################
# Path to config file
configFile=$1
#configFile=/nfs/users/asebe/cnavarrete/proj/iChIP/fastq_sel/Scer_reads/configFile.sh

#################################  HOW TO RUN ME  ##############################
#bash /nfs/users/asebe/cnavarrete/git_repo/iChIPv2/scripts/02_generate_QC.sh /nfs/users/asebe/cnavarrete/proj/iChIP/fastq_sel/Scer_reads/configFile.txt




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

# Get gathered table
outGather=${out_dir}/QC/gathered_QC.tsv
echo -e "sampleName\tnReadPairs\talignedReadPairs\talignedRate\tQ>=30\t\
            dupFiltered\tdupRate\tN_narrowPeak\tnarrowFRiP\t\
            N_broadPeak\tbroadFRiP" > ${outGather}
for sps_id in $(find ${out_dir}/QC/* -maxdepth 0 -printf "%f\n"); do
    if [[ ${sps_id} != "gathered_QC.tsv" ]]; then 
        files=$(find ${out_dir}/QC/${sps_id}/* -maxdepth 0 -printf "%f\n")
        IDs=$(for f in "${files}"; do echo ${f} | sed "s/summaryAlign_//g" | \
            sed "s/summaryPeak_//g" | sed "s/.txt//g"; done | \
            tr ' ' '\n' | sort | uniq)
        for id_ in ${IDs}; do
            alignFile=${out_dir}/QC/${sps_id}/summaryAlign_${id_}.txt
            peakFile=${out_dir}/QC/${sps_id}/summaryPeak_${id_}.txt
            
            if [ -f "${alignFile}" ]; then
                toAdd=$(tail -n 1 ${alignFile})
            else
                toAdd=$(echo -e "${id_}\t \t \t \t \t \t ")
            fi
            if [ -f "${peakFile}" ]; then
                toAdd_=$(grep -A 1 "narrowFRiP$" ${peakFile} | tail -n 1 | \
                    awk '{print $2"\t"$3}')
                toAdd=$(echo -e "${toAdd}\t${toAdd_}")
                toAdd_=$(grep -A 1 "broadFRiP$" ${peakFile} | tail -n 1 | \
                    awk '{print $2"\t"$3}')
                toAdd=$(echo -e "${toAdd}\t${toAdd_}")
            fi
            echo -e "${toAdd}" >> ${outGather}
        done
    fi
done


# Overwrite nCPU, lets leave it as single thread
#nCPU=1
# Call distribution QC script
python ${gitP}/scripts/sub-scripts/02b_py3_wholeGenome_binDistrib.py \
        ${sample_table} ${outGather} ${out_dir}/bam_files \
        ${out_dir} ${nCPU} ${condaEnv}


exit 0