#!/bin/bash
 
#########################   INPUT PARAMETERS  #############################

# Path to alignment samplesheet
peakS=$1
# Path to genome samplehseet
genomeS=$2
# Path to Output folder
out_dir=$3
# Alignment file index to select fastq files
chipLine=$4
# Path to conda environment folder
condaEnv=$5


# Set to "Yes" to delete MACS2 log files (${out_dir}/peakCalling/${mainLabel}/logs/)
# Were ${mainLabel} is the name of the folder containing the BAM files
cleanLogs="No"
# Set to "Yes" to delete MACS2 gappedPeak (only keep MACS2 xls and BED6+4 files)
# (${out_dir}/peakCalling/${mainLabel}/peaks)
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
# Make sure the input shamplesheet files format is correct
maxSection=$(awk -F ',' 'NF > maxNF {maxNF = NF} END {print maxNF+0}' ${peakS})
if [[ ${maxSection} == 3 ]] ; then 
    echo -e "Peak samplesheet has correct number of comma separated columns"
else
    echo -e "\nERROR: Peak samplesheet format might be WRONG";
    echo -e "Maximum number of comma separated sections is ${maxSection}, when should be 3\n";
    exit 1;
fi

##############################   GETING READY   ################################

# Load conda environment
#conda activate ${condaEnv}
export PATH="${condaEnv}:$PATH"

# Get alignment file index in case this is a job from a queueing system
# If not Slurm or SGE we asume its a number with the line from ${peakS}
if [[ ${chipLine} == "Slurm" ]]; then 
    chipLine=$SLURM_ARRAY_TASK_ID
elif [[ ${chipLine} == "SGE" ]]; then
    chipLine=$SGE_TASK_ID
    echo "TASK ID"
    echo $SGE_TASK_ID
fi

# Get alignment input info
## Get files to align (this is index 1, so +1 to skip header)
peakContent=`sed "$((${chipLine} + 1))q;d" ${peakS}`
content=(${peakContent//,/ })
inPath=${content[0]}
controlFile=${content[1]}
refID=${content[2]}

echo -e "inPath: ${inPath}\ncontrolFile: ${controlFile}"
echo -e "refID: ${refID}"

## Step control to make sure we have input info
if [ -z "${inPath}" ]; then 
    echo "ERROR: Sample file issue with inPath, empty/null variable";
    exit 1;
fi
if [ -z "${controlFile}" ]; then 
    echo "ERROR: Sample file issue with controlFile, empty/null variable";
    exit 1;
fi
if [ -z "${refID}" ]; then 
    echo "ERROR: Sample file issue with refID, empty/null variable";
    exit 1;
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

# Get a list of BAM files and labels
allbams=$(find ${inPath}/*bam -printf "${inPath}/%f\n" | \
            tr '\n' ' ')
allLabels=`for i in $allbams; do basename ${i} | cut -d '.' -f 1 \
                ; done | tr '\n' ' '`

# Store path to control file and deal with cases in which user provides path instead of name
if [ "${controlFile:0:1}" = "/" ]; then
     echo -e "\nWARNING: 'controlFile' should be a file name within inPath"
     echo -e "Anyways, we will use the provided file path to look for the control BAM file\n"
     controlbam="${controlFile}"
     controlFile=$(basename $controlFile)
else
    controlbam="${inPath}/${controlFile}"
fi
controlLabel=$(echo ${controlFile} | cut -d '.' -f 1)

# Get the name of the folder containing the BAM files
mainLabel=$(basename $inPath)

# Create alignment file folders
if [ ! -e ${out_dir}/peakCalling/${mainLabel}/peaks ]; then
    mkdir -p ${out_dir}/peakCalling/${mainLabel}/logs
    mkdir -p ${out_dir}/peakCalling/${mainLabel}/peaks
fi


##################################   CODE   ####################################


################################################################################
# PEAK CALLING: MACS2
################################################################################

# We shouldn't create the index withing the jobs, better if its done before
# Mapping bowtie2
echo -e "Starting Peak calling ------------------------------------------- \n"

cd ${out_dir}/peakCalling/${mainLabel}/peaks

## we go for each species
for bam in ${allbams}; do

    label=$(basename $bam | cut -d '.' -f 1)

    # proceed to call peaks for no controls
    if [[ ${controlLabel} != ${label} ]] ; then
        
        
        summaryFile="${out_dir}/QC/peakSummary_${label}.txt"

        # This variable will help you to add code for narrow peak calling as 
        # well while avoiding to get the number of reads twice
        total_reads="empty"

        # Get reference genome size from Samtools index
        genome_size=$(awk '{SUM+=$2}END{print SUM}' ${refGenome}.fai)

        ## MACS2 broad peaks
        peaktype="broadPeak"
        # check if the file exists of it was created with a previous bam version 
        fileNotExistOrOlder "${out_dir}/peakCalling/${mainLabel}/peaks/${label}_peaks_${peaktype}.xls" \
                            "${bam} ${controlbam}"
        # this outputs analyse as yes or no in lowercase
        if [[ ${analyse} == "yes" ]]; then
            # echo file and control
            echo "Bams used for peak calling:"
            echo $bam
            echo $controlbam
            echo

            if [[ $total_reads == "empty" ]]; then
                total_reads=$(sambamba view -c ${bam})
            fi

            macs2 callpeak \
                    -t ${bam} \
                    -c ${controlbam} \
                    --broad \
                    --max-gap 300 \
                    -f BAMPE \
                    -g ${genome_size} \
                    -n $label \
                    --outdir ${out_dir}/peakCalling/${mainLabel}/peaks 2> \
                    ${out_dir}/peakCalling/${mainLabel}/logs/${label}_macs2_${peaktype}.log
            
            # -q 0.05 as default (at least for broad)
            mv ${out_dir}/peakCalling/${mainLabel}/peaks/${label}_peaks.xls \
                ${out_dir}/peakCalling/${mainLabel}/peaks/${label}_peaks_${peaktype}.xls

            npeaks=$(cat \
                ${out_dir}/peakCalling/${mainLabel}/peaks/${label}_peaks.${peaktype} | \
                wc -l)
            reads_in_peaks=$(bedtools sort -i \
                ${out_dir}/peakCalling/${mainLabel}/peaks/${label}_peaks.${peaktype} \
                | bedtools merge -i stdin | bedtools intersect -u -nonamecheck \
                -a ${bam} -b stdin -ubam | sambamba view -c /dev/stdin)
            FRiP=$(awk "BEGIN {print "${reads_in_peaks}"/"${total_reads}"}")
            # report
            echo -e "NUMBER OF BROAD PEAKS\t${npeaks}" >> ${summaryFile}
            echo -e "total_reads\treads_in_peaks\tFRIP" >> ${summaryFile}
            echo -e "${total_reads}\t${reads_in_peaks}\t${FRiP}" >> ${summaryFile}
            echo -e "\n${label}: ${FRiP}" >> ${summaryFile}
        fi

        
    fi
done

echo -e "Peak calling - finished ------------------------------------------\n"


################################################################################
# Delete temporal files
################################################################################

if [[ ${cleanLogs} == "Yes" ]]; then 
    echo -e "Deleting Log files ----------------------------------------- \n"
    rm ${out_dir}/peakCalling/${mainLabel}/logs/*.log
fi
if [[ ${cleanGapped} == "Yes" ]]; then 
    echo -e "Deleting Gapped files -------------------------------------- \n"
    rm ${out_dir}/peakCalling/${mainLabel}/peaks/*.gappedPeak
fi

################################################################################
# End
################################################################################
echo -e "Finished --------------------------------------------------------\n"

exit 0