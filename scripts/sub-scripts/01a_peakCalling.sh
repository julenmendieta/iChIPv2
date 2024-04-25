#!/bin/bash
 
#########################   INPUT PARAMETERS  #############################

# Path to alignment samplesheet
sample_table=$1
# Path to genome samplehseet
genome_table=$2
# Path to Output folder
out_dir=$3
# Alignment file index to select fastq files
chipLine=$4
# Path to conda environment folder
condaEnv=$5


# Set to "Yes" to delete MACS2 log files (${out_dir}/peakCalling/${sps_id}/logs/)
# Were ${sps_id} is the ID of the species
cleanLogs="No"
# Set to "Yes" to delete MACS2 gappedPeak (only keep MACS2 xls and BED6+4 files)
# (${out_dir}/peakCalling/${sps_id})
cleanGapped="Yes"

#################################   EXTRA   ####################################
##===============================================================================

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command failed with exit code $?."' EXIT

##==============================================================================

# function to check if the given first file doesnt exist or is older than 
# the second input file
fileNotExistOrOlder () {
    # check if the file exists of it was created with a previous bam version 
    analyse="no"
    if [ ! -e $1 ]; then
        analyse="yes"
    # only proceed if the output file is older than the bam file
    # in this way if we resequenced and kept the name the analysis 
    # will be repeated
    else
        for tfile in $2; do
            if [[ $1 -ot ${tfile} ]] ; then
                analyse="yes"
                echo $1" older than "${tfile}
            fi
        done
    fi
}

##==============================================================================

##################################  Checks   ###################################
# Make sure the input Sample table file format is correct
maxSection=$(awk -F '\t' 'FNR>1 && NF > maxNF {maxNF = NF} END {print maxNF+0}' ${sample_table})
if [[ ${maxSection} == 6 ]]; then 
    echo "Sample table has correct number of TAB separated columns"
else
    echo -e "\nERROR: Sample table format might be WRONG";
    echo -e "Maximum number of TAB separated sections is ${maxSection}, when should be 6\n";
    exit 1;
fi

##############################   GETING READY   ################################

# Load conda environment
#conda activate ${condaEnv}
export PATH="${condaEnv}:$PATH"

# Get alignment file index in case this is a job from a queueing system
# If not Slurm or SGE we asume its a number with the line from ${sample_table}
if [[ ${chipLine} == "queue" ]]; then 
    # Slurm queue example
    #chipLine=$SLURM_ARRAY_TASK_ID
    # SGE queue example
    chipLine=$SGE_TASK_ID

    # Make sure user set up this correctly
    ## Define file label
    if [ -z "${chipLine}" ]; then
        echo -e "\nERROR: You decided to call the script as jobMode='queue' but \
                     the variable to retrieve TASK ID from the job array is NULL";
        echo -e "Fix it and try again: Line 82 (chipLine) in 01a_peakCaling.sh.\n";
        exit 1;
    fi
fi

# Get alignment input info
## Get files to align (this is index 1, so +1 to skip header)
peakContent=`sed "$((${chipLine} + 1))q;d" ${sample_table}`
content=(${peakContent//\\t/ })
sample_id=${content[0]}
sps_id=${content[3]}
is_control=${content[4]}
controlFile=${content[5]}

## Step control to make sure we have input info
if [ -z "${sample_id}" ]; then 
    echo "ERROR: Sample file issue with sample_id, empty/null variable";
    exit 1;
fi
if [ -z "${sps_id}" ]; then 
    echo "ERROR: Sample file issue with sps_id, empty/null variable";
    exit 1;
fi
if [ -z "${is_control}" ]; then 
    echo "ERROR: Sample file issue with is.control, empty/null variable";
    exit 1;
fi
if [ -z "${controlFile}" ]; then 
    echo "WARNING: ControlFile variable is empty/null. We will run peak calling\
 without background";
fi

# End job if this is a control
if [[ ${is_control} == "Y" ]]; then 
    echo "Ops, we ran a peak calling job for a control. Stoping here...";
    exit 0;
fi

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
if [ -f "${refGenome}.fai" ]; then
    echo "Fasta Index file exists"
else 
    echo -e "\nERROR. There is no Fasta .fai Index for ${refGenome}"
    echo "We would spect to find a file called '${refGenome}.fai'"
    exit 1
fi

# Get selected bam file for peak calling
## Main sample
bamfile="${out_dir}/bam_files/${sps_id}/${sample_id}.q30.rmdup.bam"
if ! [ -f "${bamfile}" ]; then
    echo -e "\nERROR: BAM file for ${sample_id}: ${bamfile} doesn't exist"
    exit 1
fi

## Control
## Store path to control file and deal with cases in which user provides path 
## instead of name
if ! [ -z "${controlFile}" ]; then 
    if [ "${controlFile:0:1}" = "/" ]; then
        echo -e "\nWARNING: Control file for peak calling will come from a folder \
    different from\n${out_dir}/bam_files"
        echo ${controlFile}
        controlbam="${controlFile}"
        controlFile=$(basename $controlFile)
    else
        controlbam="${out_dir}/bam_files/${sps_id}/${controlFile}.q30.rmdup.bam"
    fi
    if ! [ -f "${controlbam}" ]; then
        echo -e "\nERROR: BAM file for ${controlFile}: ${controlbam} doesn't exist"
        exit 1
    fi
fi


# Create peak calling file folders
if [ ! -e ${out_dir}/peakCalling/${sps_id}/logs ]; then
    mkdir -p ${out_dir}/peakCalling/${sps_id}/logs
fi
if [ ! -e ${out_dir}/QC/${sps_id} ] ; then
    mkdir -p ${out_dir}/QC/${sps_id}
fi

################################################################################
# PEAK CALLING: MACS2
################################################################################

# We shouldn't create the index withing the jobs, better if its done before
# Mapping bowtie2
echo -e "Starting Peak calling ------------------------------------------- \n"

cd ${out_dir}/peakCalling/${sps_id}

# proceed to call peaks for no controls
summaryFile="${out_dir}/QC/${sps_id}/summaryPeak_${sample_id}.txt"

# This variable will help you to add code for narrow peak calling as 
# well while avoiding to get the number of reads twice
total_reads="empty"

# Get reference genome size from Samtools index
genome_size=$(awk '{SUM+=$2}END{print SUM}' ${refGenome}.fai)


## MACS2 narrow peaks
peaktype="narrowPeak"
# check if the file exists of it was created with a previous bam version 
fileNotExistOrOlder "${out_dir}/peakCalling/${sps_id}/${sample_id}_peaks_${peaktype}.xls" \
                    "${bamfile} ${controlbam}"
# this outputs analyse as yes or no in lowercase
if [[ ${analyse} == "yes" ]]; then
    # echo file and control
    echo "Bams used for narrow peak calling:"
    echo ${bamfile}
    echo ${controlbam}
    echo

    if ! [ -z "${controlFile}" ]; then 
        echo -e "sps_id: ${sps_id}\nout_dir: ${out_dir}" > ${summaryFile}
        echo -e "bamFile: ${bamfile}\ncontrolBam:${controlFile}\n" >> ${summaryFile}
    else
        echo -e "sps_id: ${sps_id}\nout_dir: ${out_dir}" > ${summaryFile}
        echo -e "bamFile: ${bamfile}\ncontrolBam: None\n" >> ${summaryFile}
    fi

    if [[ ${total_reads} == "empty" ]]; then
        # This also gets orphan reads, so the number might be different to
        # the values from the previous scripts, that looks for aligned pairs
        total_reads=$(sambamba view -c ${bamfile})
    fi

    macs2 callpeak \
            -t ${bamfile} \
            -c ${controlbam} \
            --max-gap 300 \
            -f BAMPE \
            -g ${genome_size} \
            -n ${sample_id} \
            --outdir ${out_dir}/peakCalling/${sps_id} 2> \
            ${out_dir}/peakCalling/${sps_id}/logs/${sample_id}_macs2_${peaktype}.log
    
    # -q 0.05 as default (at least for broad)
    mv ${out_dir}/peakCalling/${sps_id}/${sample_id}_peaks.xls \
        ${out_dir}/peakCalling/${sps_id}/${sample_id}_peaks_${peaktype}.xls

    npeaksNarrow=$(cat \
        ${out_dir}/peakCalling/${sps_id}/${sample_id}_peaks.${peaktype} | \
        wc -l)
    reads_in_peaksNarrow=$(bedtools sort -i \
        ${out_dir}/peakCalling/${sps_id}/${sample_id}_peaks.${peaktype} \
        | bedtools merge -i stdin | bedtools intersect -u -nonamecheck \
        -a ${bamfile} -b stdin -ubam | sambamba view -c /dev/stdin)
    FRiP_narrow=$(awk "BEGIN {print "${reads_in_peaksNarrow}"/"${total_reads}"}")
    # report
    echo -e "sampleName\tN_narrowPeak\tnarrowFRiP" >> ${summaryFile}
    echo -e "${sample_id}\t${npeaksNarrow}\t${FRiP_narrow}" >> ${summaryFile}
 
fi


## MACS2 broad peaks
peaktype="broadPeak"
# check if the file exists of it was created with a previous bam version 
fileNotExistOrOlder "${out_dir}/peakCalling/${sps_id}/${sample_id}_peaks_${peaktype}.xls" \
                    "${bamfile} ${controlbam}"
# this outputs analyse as yes or no in lowercase
if [[ ${analyse} == "yes" ]]; then
    # echo file and control
    echo "Bams used for broad peak calling:"
    echo ${bamfile}
    echo ${controlbam}
    echo

    if [[ ${total_reads} == "empty" ]]; then
        total_reads=$(sambamba view -c ${bamfile})
    fi

    macs2 callpeak \
            -t ${bamfile} \
            -c ${controlbam} \
            --broad \
            --max-gap 300 \
            -f BAMPE \
            -g ${genome_size} \
            -n ${sample_id} \
            --outdir ${out_dir}/peakCalling/${sps_id} 2> \
            ${out_dir}/peakCalling/${sps_id}/logs/${sample_id}_macs2_${peaktype}.log
    
    # -q 0.05 as default (at least for broad)
    mv ${out_dir}/peakCalling/${sps_id}/${sample_id}_peaks.xls \
        ${out_dir}/peakCalling/${sps_id}/${sample_id}_peaks_${peaktype}.xls

    npeaksBroad=$(cat \
        ${out_dir}/peakCalling/${sps_id}/${sample_id}_peaks.${peaktype} | \
        wc -l)
    reads_in_peaksBroad=$(bedtools sort -i \
        ${out_dir}/peakCalling/${sps_id}/${sample_id}_peaks.${peaktype} \
        | bedtools merge -i stdin | bedtools intersect -u -nonamecheck \
        -a ${bamfile} -b stdin -ubam | sambamba view -c /dev/stdin)
    FRiP_broad=$(awk "BEGIN {print "${reads_in_peaksBroad}"/"${total_reads}"}")
    # report
    echo -e "sampleName\tN_broadPeak\tbroadFRiP" >> ${summaryFile}
    echo -e "${sample_id}\t${npeaksBroad}\t${FRiP_broad}" >> ${summaryFile}
    
fi

# Get maximum FRiP value
FRiP_narrow=$(grep -A 1 "narrowFRiP" ${summaryFile} | tail -n 1 | \
                awk '{print $3}')
FRiP_broad=$(grep -A 1 "broadFRiP" ${summaryFile} | tail -n 1 | \
                awk '{print $3}')
if (( $(echo "$FRiP_narrow >= $FRiP_broad" |bc -l) )); then
    echo -e "\nmaxFRiP: narrowFRiP = ${FRiP_narrow}" >> ${summaryFile}
else
    echo -e "\nmaxFRiP: broadFRiP = ${FRiP_broad}" >> ${summaryFile}
fi

echo -e "Peak calling - finished ------------------------------------------\n"


################################################################################
# Delete temporal files
################################################################################

if [[ ${cleanLogs} == "Yes" ]]; then 
    echo -e "Deleting Log files ----------------------------------------- \n"
    if [ -f "${out_dir}/peakCalling/${sps_id}/logs/${sample_id}_macs2_${peaktype}.log" ]; then
        rm ${out_dir}/peakCalling/${sps_id}/logs/${sample_id}_macs2_${peaktype}.log
    fi
fi
if [[ ${cleanGapped} == "Yes" ]]; then 
    echo -e "Deleting Gapped files -------------------------------------- \n"
    if [ -f "${out_dir}/peakCalling/${sps_id}/*.gappedPeak" ]; then
        rm ${out_dir}/peakCalling/${sps_id}/${sample_id}_peaks.gappedPeak
    fi
fi

################################################################################
# End
################################################################################
echo -e "Finished --------------------------------------------------------\n"

exit 0