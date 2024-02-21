#!/bin/bash
 
#########################   INPUT PARAMETERS  #############################

# Path to alignment samplesheet
alignS=$1
# Path to genome samplehseet
genomeS=$2
# Path to Output folder
out_dir=$3
# Number of threads for each job
nCPU=$4 
# Alignment file index to select fastq files
alignLine=$5


#################################   EXTRA   ####################################
##===============================================================================

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command failed with exit code $?."' EXIT

##===============================================================================


##################################   CODE   ####################################

# Get alignment file index in case this is a job from a queueing system
# If not Slurm or Torque we asume its a number with the line from ${alignS}
if [[ ${alignLine} == "Slurm" ]]; then 
    alignLine=$SLURM_ARRAY_TASK_ID
elif [[ ${jobMode} == "Torque" ]]; then
    echo "Pending to be done"
    exit 1
fi

# Create alignment file folders
if [ ! -e ${out_dir}/bam_files/valid/all ]; then
    mkdir -p ${out_dir}/bam_files/valid/all
fi
if [ ! -e ${out_dir}/bw_files ]; then
	mkdir -p ${out_dir}/bw_files
fi
if [ ! -e ${out_dir}/temp ]; then
	mkdir -p ${out_dir}/temp
fi
if [ ! -e ${out_dir}/QC/ ] ; then
    mkdir -p ${out_dir}/QC/
fi


# Get alignment input info
## Get files to align (this is index 1, so +1 to skip header)
alignContent=`sed "$((${alignLine} + 1))q;d" ${alignS}`
content=(${alignContent//,/ })
fastq_1=${content[0]}
fastq_2=${content[1]}
refID=${content[2]}
newName=${content[3]}

echo "${fastq_1} ${fastq_2} ${refID} ${newName}"
## Define file label
if [ ! -z "${newName}" ]; then 
    fileLabel=${newName};
else 
    f1=$(basename ${fastq_1})
    f2=$(basename ${fastq_2})
    # Get common text in read1 and read2 names
    p=0
    while [[ ${f1:0:p} == ${f2:0:p} ]] ; do 
        (( ++p ))
    done
    # Store common text and delete last character if its '_'
    fileLabel=$(echo ${f1:0:p-2} | sed 's/_$//')
fi

# Get reference genome path
for i in $(cat ${genomeS}); do
    content=(${i//,/ })
    refID_=${content[0]}
    refGenome_=${content[1]}
    if [[ ${refID} == ${refID_} ]]; then 
        refGenome=$(echo ${refGenome_} | tr -d '\r')

    fi
done


# Create step control and ouput info files
stepControl="${out_dir}/temp/pipelineStep_${fileLabel}.txt"
summaryFile="${out_dir}/QC/summary_${fileLabel}.txt"

echo -e "${fileLabel}\n" > ${summaryFile}


echo -e "Aligning ${fileLabel}..."


##########################################################
# Alignment
##########################################################

# We shouldn't create the index withing the jobs, better if its done before
# Mapping bowtie2
echo -e "Bowtie2 Alignment"

bamfile=${out_dir}/temp/${fileLabel}.bam
bowtie2 -X2000 --mm --threads ${nCPU} \
        -x "${refGenome}" -1 ${fastq_1} -2 ${fastq_2} \
        2>${out_dir}/temp/${fileLabel}.log | \
    sambamba view -t ${nCPU} -f bam -S /dev/stdin | \
    sambamba sort -t ${nCPU} /dev/stdin -o ${bamfile}

echo -e "${fileLabel}\n" > ${summaryFile}
echo -e "BOWTIE2 ALIGNMENT" >> ${summaryFile}
cat ${out_dir}/temp/${fileLabel}.log >> ${summaryFile}
# Get flagstat results
echo -e "\nFlagstat" >> ${summaryFile}
sambamba flagstat ${bamfile} -t ${nCPU} >> ${summaryFile}
# get number of aligned reads (both mates)
nAligned=$(sambamba view -c -F 'first_of_pair' --num-filter 3 ${bamfile})


##########################################################
# Quality filtering with Sambamba
##########################################################

bamfile_q30=${out_dir}/temp/${fileLabel}.q30.bam
sambamba view -f bam -F "mapping_quality >= 30" ${bamfile} | \
    sambamba sort -t ${nCPU} /dev/stdin -o ${bamfile_q30}

count_q30=$(sambamba view -c -F 'first_of_pair' --num-filter 3 ${bamfile_q30} \
            -t ${nCPU})
echo -e "\n\nQUALITY FILTERING" >> ${summaryFile}
echo -e "Reads after q30 filtering:\t${count_q30}" >> ${summaryFile}


##########################################################
# Remove duplicates with Picard 
##########################################################

bamfile_q30rm=${out_dir}/bam_files/${fileLabel}.q30.rmdup.bam
sambamba markdup -r -t ${nCPU} ${bamfile_q30} ${bamfile_q30rm}

count_dup=$(sambamba view -c -F 'first_of_pair' --num-filter 3 \
            ${bamfile_q30rm} -t ${nCPU})
echo -e "\n\nDUPLICATE REMOVAL" >> ${summaryFile}
echo -e "Reads after duplicate removal:\t${count_dup}" >> ${summaryFile}
# Get flagstat results
echo -e "\nFlagstat" >> ${summaryFile}
sambamba flagstat ${bamfile_q30rm} -t ${nCPU} >> ${summaryFile}


##########################################################
# Store final summary
##########################################################
# Add alignment summary
#alignRate=$(cat $summaryFile | grep "overall alignment rate" | \
#    sed 's/ overall.*//g') # Rate counting separately R1 & R2
nReads=$(cat $summaryFile | grep "reads; of these:" | sed 's/ reads.*//g')
alignRate=$(echo "print(f'{round(${nAligned}/${nReads} * 100, 2)}%')" | python)
qFiltered=$((${nAligned} - ${count_q30}))
dupRate=$(echo "print(f'{100 - round(${count_dup}/${count_q30} * 100, 2)}%')" \
            | python)
echo -e "\n\nFINAL ALIGNMENT SUMMARY" >> ${summaryFile}
echo -e "sampleName\tnReadPairs\talignedReadPairs\talignedRate\tQ>=30\t\
        dupFiltered\tdupRate" >> ${summaryFile}
echo -e "${fileLabel}\t${nReads}\t${nAligned}\t${alignRate}\t${count_q30}\t\
        ${count_dup}\t${dupRate}" >> ${summaryFile}


##########################################################
# Create Bigwigs
##########################################################


##########################################################
# End
##########################################################
echo -e "Finished -----------------------------------------------------------\n"

exit 0