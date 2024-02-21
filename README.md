# iChIPv2

!['anyImage?'](./img/logo.png)

This pipeline shows the steps followed in (--paper--) to process iChIPv2 data from Fastq files to peak coordinates and further analysis.
It takes as input Fastq files and reference genomes and organises the data by the use of 3 Samplesheets: i) reference genomes, ii) fastq files info, and iii) peak calling info.

# Installation  
If you don't have it yet, first start by installing Conda [Miniforge](https://github.com/conda-forge/miniforge#miniforge).

We can then move to install tool dependencies with the yaml file found [here](iChIPv2_environment.yml):

```bash
# create environment
conda env create -n iChIPv2 -f iChIPv2_environment.yml

# activate environment
conda activate iChIPv2
```

# Sample sheets format
File examples [here](examples/sampleSheets)
## 1. Alignment Samplesheet
Comma separated file with 3 mandatory + 1 optional columns:  
| Column ID        | Description          | 
| ------------- |:-------------:|
| **fastq_1**      | FastQ file name for Illumina short reads 1. File has to be gzipped | 
| **fastq_2**     | FastQ file name for Illumina short reads 2. File has to be gzipped     |  
| **refID** | ID to be linked with a reference genome in the Genome Samplesheet    |
| newName (optional) |  In case you want to name output files using this label instead of file name    |

 ```bash
# Trick to create a sample sheet if your files are named as cellID_extra?_S[1-9].R[12]_001.fastq.gz
# We will link cellID to a reference genome in the "Genome Samplesheet"
# Command to be run inside the folder containing the Fastq.gz files
# Examples: Nvec_H3K36me3_301121_S27_R1_001.fastq.gz & Nvec_H3K36me3_301121_S27_R2_001.fastq.gz
for i in *fastq.gz; do
    echo $i | sed 's/_R._001.fastq.gz//g' ;
done | sort | uniq > samplesNames.txt

echo "fastq_1,fastq_2,refID" > alignment_samplesheet.csv ;
for i in $(cat samplesNames.txt); do 
    f1=$(find $PWD -name "${i}_R1_001.fastq.gz");
    f2=$(find $PWD -name "${i}_R2_001.fastq.gz");
    refID=(${i/_/ }); refID=${refID[0]};
    echo "${f1},${f2},${refID}";
done >> alignment_samplesheet.csv
```

## 2. Genome Samplesheet
Comma separated file with 2 mandatory columns:  
| Column ID        | Description          | 
| ------------- |:-------------:|
| **refID**      | Cell or species ID | 
| **refGenome**     | Path to reference genome fasta file   |  

**NOTE:** All fasta files need to have Bowtie2 and Samtools indexes. How to get them:
```bash
# Bowtie2 index
bowtie2-build SPECIES-REF-GENOME.fasta
# Samtools index (from same folder as file)
samtools faidx SPECIES-REF-GENOME.fasta
```

## 3. Peak calling Samplesheet
Comma-separated file with 3 mandatory columns. Each "inPath" folder must only contain files for the same species:  
| Column ID        | Description          | 
| ------------- |:-------------:|
| **inPath**      | Path to folder containing BAM files to be analysed | 
| **controlFile**     | Full name of the control bamfile (input, IgG, H3, etc.) for all the samples in **inPath** folder   |  
| **refID** | ID to be linked with a reference genome in the Genome Samplesheet  |


# Running the scripts
## Test run
Once dependencies are installed, make sure the pipeline works by running `pending` with sample data:

```bash
Code to run test
```

## Run it with your own data
### 1. Create the Samplesheet files
Create the Samplesheet files and move them (or not) to a preferred location.

### 2. Run alignment
Open 00_alignmentLauncher.sh script and modify the variables in USER-DEFINED INPUT and the queuing system parameters if needed. Then execute it with bash or within a cluster job (only valid for "noQueue" mode).  
The script has code that works in Slurm and (Torque FUTURE) queueing systems. However, if it is set to "noQueue", all the files in the Alignment Samplesheet will be aligned one after the other without parallelization. This is the only approach to use if you run the analysis on a computer without a queueing system and the best one to use if our code to run jobs is not working in your system. If so, you can always create multiple Alignment Samplesheets and run them through your queuing system in "noQueue" mode.

Output files:  
-BAM files --> ${out_dir}/bam_files  
-Bigwig files --> ${out_dir}/bw_files  
-Alignment information files --> ${out_dir}/QC  
-Alignment step files (to restart broken runs) and temporal files -->  ${out_dir}/temp  
  
Note that we store a file named "pipelineStep_${fileLabel}.txt" in "${out_dir}/temp". This file contains one ID for each successful pipeline step and must be deleted if, for any reason, we want to re-start from zero the alignment of a given file or align files with the same name as previously analyzed ones. 

### 3. Check the alignment output
Take a look at your data and discard failed experiments. Once you're sure to continue ahead with the peak calling, we recommend moving the bam files to "${out_dir}/bam_files/valid/all". Then, create a species-specific folder in "${out_dir}/bam_files/valid" and link the selected bam files. Example:
```bash
out_dir="/PATH/TO/OUTPUT"
# Create a folder for this set of BAM files. In this case, we call it Scer for Saccharomyces cerevisiae
mkdir ${out_dir}/bam_files/valid/Scer
# Link Scer bam files to that folder (will be accessible from there but not duplicated)
# The link has to be than from inside the destination folder
cd ${out_dir}/bam_files/valid/Scer
ln -s ${out_dir}/bam_files/valid/all/Scer* .
```

### 4. Run broad peak calling with MACS2
Open 01_peakCallingLauncher.sh script and modify the variables in USER-DEFINED INPUT and the queueing system parameters if needed. Then execute it with bash or within a cluster job (only valid for "noQueue" mode). The details regarding the queueing system are the same as in the Alignment step.

Output files:  
-MACS2 output files are stored in ${out_dir}/peakCalling/${mainLabel}, were ${mainLabel} is the same as the name of the folder containing the analysed BAM files (Scer in our previous example)  
-MACS2 peak files --> ${out_dir}/peakCalling/${mainLabel}/peaks  
-MACS2 Log files --> ${out_dir}/peakCalling/${mainLabel}/logs  
-Per file FRiP values --> ${out_dir}/QC/peakSummary_*.txt  

