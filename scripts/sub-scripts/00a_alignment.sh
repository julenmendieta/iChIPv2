#!/bin/bash
 
#########################   INPUT PARAMETERS  #############################

# Path to alignment samplesheet
sample_table=$1
# Path to genome samplehseet
genome_table=$2
# Path to input FastQ folder
in_dir=$3
# Path to Output folder
out_dir=$4
# Number of threads for each job
nCPU=$5 
# Alignment file index to select fastq files
alignLine=$6
# Path to conda environment folder
condaEnv=$7

# Set to "Yes" to remove temporal files
cleanUp="Yes"

#################################   EXTRA   ####################################
##===============================================================================

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command failed with exit code $?."' EXIT

##===============================================================================

##################################  CHECK   ####################################

# Make sure the input Sample table file format is correct
maxSection=$(awk -F '\t' 'FNR>1 && NF > maxNF {maxNF = NF} END {print maxNF+0}' ${sample_table})
if [[ ${maxSection} == 6 ]]; then 
    echo "Sample table has correct number of TAB separated columns"
else
    echo -e "\nERROR: Sample table format might be WRONG";
    echo -e "Maximum number of TAB separated sections is ${maxSection}, when should be 6\n";
    exit 1;
fi
# Make sure the input Genome table file format is correct
maxSection=$(awk -F '\t' 'FNR>1 && NF > maxNF {maxNF = NF} END {print maxNF+0}' ${genome_table})
if [[ ${maxSection} == 2 ]]; then 
    echo "Genome table has correct number of TAB separated columns"
else
    echo -e "\nERROR: Genome table format might be WRONG":
    echo -e "Maximum number of TAB separated sections is ${maxSection}, when should be 2\n";
    exit 1
fi


##############################   GETING READY   ################################


# Load conda environment
#conda activate ${condaEnv}
export PATH="${condaEnv}:$PATH"

# Get alignment file index in case this is a job from a queueing system
# If not Slurm or SGE we asume its a number with the line from ${sample_table}
if [[ ${alignLine} == "queue" ]]; then 
    # Slurm queue example
    #alignLine=$SLURM_ARRAY_TASK_ID
    # SGE queue example
    alignLine=$SGE_TASK_ID

    # Make sure user set up this correctly
    ## Define file label
    if [ -z "${alignLine}" ]; then
        echo -e "\nERROR: You decided to call the script as jobMode='queue' but \
the variable to retrieve TASK ID from the job array is NULL";
        echo -e "Fix it and try again: Lines 71 (alignLine) and 82 (chipLine) \
in 00a_alignment.sh and 01a_peakCaling.sh, repectively.\n";
        exit 1;
    fi

fi

# Get alignment input info
## Get files to align (this is index 1, so +1 to skip header)
alignContent=`sed "$((${alignLine} + 1))q;d" ${sample_table}`
content=(${alignContent//\\t/ })
sample_id=${content[0]}
fastq_1=${content[1]}
fastq_2=${content[2]}
sps_id=${content[3]}

echo -e "in_dir: ${in_dir}\n"
echo -e "\nfastq_1: ${fastq_1}\nfastq_2: ${fastq_2}"
echo -e "sps_id: ${sps_id}\nsample_id: ${sample_id}\n"

## Step control to make sure we have input info
if [ -z "${sample_id}" ]; then 
    echo "ERROR: Sample file issue with sample_id"
    exit 1
fi
if [ -z "${fastq_1}" ]; then 
    echo "ERROR: Sample file issue with fastq_1"
    exit 1
fi
if [ -z "${fastq_2}" ]; then 
    echo "ERROR: Sample file issue with fastq_2"
    exit 1
fi
if [ -z "${sps_id}" ]; then 
    echo "ERROR: Sample file issue with sps_id"
    exit 1
fi

# Create alignment file folders
if [ ! -e ${out_dir}/bam_files/${sps_id} ]; then
    mkdir -p ${out_dir}/bam_files/${sps_id}
fi
if [ ! -e ${out_dir}/bw_files/${sps_id} ]; then
	mkdir -p ${out_dir}/bw_files/${sps_id}
fi
if [ ! -e ${out_dir}/temp ]; then
	mkdir -p ${out_dir}/temp
fi
if [ ! -e ${out_dir}/QC/${sps_id} ] ; then
    mkdir -p ${out_dir}/QC/${sps_id}
fi

# ## Define file label
# maxSection=$(echo $alignContent | awk -F ',' 'NF > maxNF {maxNF = NF} END {print maxNF+0}')
# if [[ ${maxSection} == 4 ]]; then 
#     sample_id=${newName};
# else 
#     f1=$(basename ${fastq_1})
#     f2=$(basename ${fastq_2})
#     # Get common text in read1 and read2 names
#     p=0
#     while [[ ${f1:0:p} == ${f2:0:p} ]] ; do 
#         (( ++p ))
#     done
#     # Store common text and delete last character if its '_'
#     sample_id=$(echo ${f1:0:p-2} | sed 's/_$//')
# fi

# Get reference genome path
nLines=$(wc -l ${genome_table}| awk '{print $1}')
for i in $(seq 1 ${nLines}); do
    content=`sed "$((${i}))q;d" ${genome_table}`
    content=(${content//\\t/ })
    sps_id_=${content[0]}
    refGenome_=${content[1]}
    if [[ ${sps_id} == ${sps_id_} ]]; then 
        refGenome=$(echo ${refGenome_} | tr -d '\r')

    fi
done
# Check if it has bowtie two indexes
if [ -f "${refGenome}.1.bt2" ]; then
    echo "Bowtie2 Index file exists"
else 
    echo -e "\nERROR. There is no Bowtie 2 Index for ${refGenome}"
    echo "We would spect to find a file called '${refGenome}.1.bt2' among others"
    exit 1
fi



# Create step control and ouput info files
stepControl="${out_dir}/temp/pipelineStep_${sample_id}.txt"
if [ ! -e ${stepControl} ] ; then
    touch ${stepControl}
fi
summaryFile="${out_dir}/QC/${sps_id}/summaryAlign_${sample_id}.txt"


echo -e "Aligning ${sample_id}..."


################################################################################
# Alignment
################################################################################

# We shouldn't create the index withing the jobs, better if its done before
# Mapping bowtie2
echo -e "Starting Bowtie2 Alignment -------------------------------------- \n"
bamfile=${out_dir}/temp/${sample_id}.bam

# Add full path to input fastq files
fastq_1=(${fastq_1//,/ })
fastq_1=("${fastq_1[@]/#/${in_dir}}")
for file in "${fastq_1[@]}"; do
    if ! [ -f "${file}" ]; then
        echo "ERROR: File fastq_1: ${file} doesn't exist"
        exit 1
    fi
done
fastq_1=$(IFS=, ; echo "${fastq_1[*]}")

fastq_2=(${fastq_2//,/ })
fastq_2=("${fastq_2[@]/#/${in_dir}}")
for file in "${fastq_2[@]}"; do
    if ! [ -f "${file}" ]; then
        echo "ERROR: File fastq_2: ${file} doesn't exist"
        exit 1
    fi
done
fastq_2=$(IFS=, ; echo "${fastq_2[*]}")

# check content of first line of step control file
linec=`sed "1q;d" ${stepControl}`
if [[ ${linec} != "Alignment" ]]; then 
    echo -e "${sample_id}\n" > ${summaryFile}
    echo -e "in_dir: ${in_dir}\n" > ${summaryFile}
    echo -e "\nfastq_1: ${fastq_1}\nfastq_2: ${fastq_2}" > ${summaryFile}
    echo -e "sps_id: ${sps_id}\nsample_id: ${sample_id}\n" > ${summaryFile}

    # bowtie2 align
    bowtie2 -X2000 --mm --threads ${nCPU} \
            -x "${refGenome}" -1 ${fastq_1} -2 ${fastq_2} \
            2>${out_dir}/temp/${sample_id}.log | \
        sambamba view -t ${nCPU} -f bam -S /dev/stdin | \
        sambamba sort -t ${nCPU} /dev/stdin -o ${bamfile}

    # The issue of storing bowtie2 output 

    # Alignment info storage
    echo -e "${sample_id}\n" > ${summaryFile}
    echo -e "BOWTIE2 ALIGNMENT" >> ${summaryFile}
    cat ${out_dir}/temp/${sample_id}.log >> ${summaryFile}
    # Get flagstat results
    echo -e "\nFlagstat" >> ${summaryFile}
    sambamba flagstat ${bamfile} -t ${nCPU} >> ${summaryFile}

    echo -e "Bowtie2 Alignment - done ------------------------------------ \n"
    # store stage control info
    echo "Alignment" >> ${stepControl}
else
    echo -e "Bowtie2 Alignment - already done before --------------------- \n"
fi
# get number of aligned reads (both mates)
nAligned=$(sambamba view -c -F 'first_of_pair' --num-filter 3 ${bamfile})


################################################################################
# Quality filtering with Sambamba
################################################################################

echo -e "Starting Quality filtering -------------------------------------- \n"
bamfile_q30=${out_dir}/temp/${sample_id}.q30.bam

# check content of second line of step control file
linec=`sed "2q;d" ${stepControl}`
if [[ ${linec} != "Filtering" ]]; then 

    sambamba view -f bam -F "mapping_quality >= 30" ${bamfile} | \
        sambamba sort -t ${nCPU} /dev/stdin -o ${bamfile_q30}

    count_q30=$(sambamba view -c -F 'first_of_pair' \
                    --num-filter 3 ${bamfile_q30} -t ${nCPU})

    # Filtering info storage
    echo -e "\n\nQUALITY FILTERING" >> ${summaryFile}
    echo -e "Reads after q30 filtering:\t${count_q30}" >> ${summaryFile}

    echo -e "Quality filtering - done ------------------------------------ \n"
    # store stage control info
    echo "Filtering" >> ${stepControl}
else
    count_q30=$(sambamba view -c -F 'first_of_pair' \
                    --num-filter 3 ${bamfile_q30} -t ${nCPU})
    echo -e "Quality filtering - already done before --------------------- \n"
fi


################################################################################
# Remove duplicates with Picard 
################################################################################

echo -e "Starting Duplicate removal -------------------------------------- \n"
bamfile_q30rm=${out_dir}/bam_files/${sps_id}/${sample_id}.q30.rmdup.bam

# check content of third line of step control file
linec=`sed "3q;d" ${stepControl}`
if [[ ${linec} != "Duplicates" ]]; then 

    # mark suplicates with sambamba
    sambamba markdup -r -t ${nCPU} ${bamfile_q30} ${bamfile_q30rm}
    count_dup=$(sambamba view -c -F 'first_of_pair' --num-filter 3 \
                ${bamfile_q30rm} -t ${nCPU})

    # Duplicate removal info storage
    echo -e "\n\nDUPLICATE REMOVAL" >> ${summaryFile}
    echo -e "Reads after duplicate removal:\t${count_dup}" >> ${summaryFile}
    # Get flagstat results
    echo -e "\nFlagstat" >> ${summaryFile}
    sambamba flagstat ${bamfile_q30rm} -t ${nCPU} >> ${summaryFile}

    echo -e "Duplicate removal - done ------------------------------------ \n"
    # store stage control info
    echo "Duplicates" >> ${stepControl}
else
    count_dup=$(sambamba view -c -F 'first_of_pair' --num-filter 3 \
                ${bamfile_q30rm} -t ${nCPU})
    echo -e "Duplicate removal - already done before --------------------- \n"
fi


################################################################################
# Store final summary
################################################################################

echo -e "Starting Alignment summary -------------------------------------- \n"

# check content of forth line of step control file
linec=`sed "4q;d" ${stepControl}`
if [[ ${linec} != "Summary" ]]; then 

    # Get alinment info
    #alignRate=$(cat $summaryFile | grep "overall alignment rate" | \
    #    sed 's/ overall.*//g') # Rate counting separately R1 & R2
    nReads=$(cat $summaryFile | grep "reads; of these:" | sed 's/ reads.*//g')
    alignRate=$(echo "print(f'{round(${nAligned}/${nReads} * 100, 2)}%')" \
                | python)
    qFiltered=$((${nAligned} - ${count_q30}))
    dupRate=$(echo "print(f'{100 - round(${count_dup}/${count_q30} * 100, 2)}%')" \
                | python)

    # Store alignment info
    echo -e "\n\nFINAL ALIGNMENT SUMMARY" >> ${summaryFile}
    echo -e "sampleName\tnReadPairs\talignedReadPairs\talignedRate\tQ>=30\t\
            dupFiltered\tdupRate" >> ${summaryFile}
    echo -e "${sample_id}\t${nReads}\t${nAligned}\t${alignRate}\t${count_q30}\t\
            ${count_dup}\t${dupRate}" >> ${summaryFile}

    echo -e "Alignment summary - done ------------------------------------ \n"
    # store stage control info
    echo "Summary" >> ${stepControl}
else
    echo -e "Alignment summary - already done before --------------------- \n"
fi

################################################################################
# Create Bigwigs
################################################################################

echo -e "Starting Bigwig creation --------------------------------------- \n"
bigWigOut=${out_dir}/bw_files/${sps_id}/${sample_id}.q30.rmdup.CPM.bw
# check content of fifth line of step control file
linec=`sed "5q;d" ${stepControl}`
if [[ ${linec} != "Bigwig" ]]; then 
    
    bamCoverage -b ${bamfile_q30rm} -o ${bigWigOut} \
                -e --scaleFactor 1 --normalizeUsing CPM --outFileFormat bigwig \
                --binSize 10 --numberOfProcessors ${nCPU}

    echo -e "Bigwig creation - done ------------------------------------- \n"
    # store stage control info
    echo "Bigwig" >> ${stepControl}
else
    echo -e "Bigwig creation - already done before ---------------------- \n"
fi

################################################################################
# Delete temporal files
################################################################################

if [[ ${cleanUp} == "Yes" ]]; then 
    echo -e "Deleting temporal files ------------------------------------ \n"
    rm ${out_dir}/temp/${sample_id}.*bam
    rm ${out_dir}/temp/${sample_id}.*bai
    rm ${out_dir}/temp/${sample_id}.*log
fi


################################################################################
# End
################################################################################
echo -e "Finished --------------------------------------------------------\n"

exit 0
