# iChIPv2

!['anyImage?'](./img/logo.png)

This pipeline shows the steps followed in (--paper--) to process iChIPv2 data from Fastq files to peak coordinates and further analysis.
It takes as input Fastq files and reference genomes and organises the data by the use of 3 Samplesheets: i) reference genomes, ii) fastq files info, and iii) peak calling info.

# Installation  
If you don't have it yet, first start by installing Conda [Miniforge](https://github.com/conda-forge/miniforge#miniforge).

We can then move to install tool dependencies with the yaml file found [here](iChIPv2_environment.yaml):

```bash
# create environment
conda env create -n iChIPv2 -f iChIPv2_environment.yaml

# activate environment
conda activate iChIPv2
```

# Test run
Once dependencies are installed, make sure the pipeline works by running `pending` with sample data:

```bash
Code to run test
```
# Sample sheets format
## 1. Alignment Samplesheet
Comma separated file with 3 mandatory + 1 optional columns:  
| Column ID        | Description          | 
| ------------- |:-------------:|
| **fastq_1**      | FastQ file name for Illumina short reads 1. File has to be gzipped | 
| **fastq_2**     | FastQ file name for Illumina short reads 2. File has to be gzipped     |  
| **refID** | ID to be linked with a reference genome in the Genome Sampleshee    |
| newName (optional) |  In case you want to name output files using this label instead of file name    |

 ```bash
# Trick to create sample sheet if your files are named as cellID_extra?_S[1-9].R[12]_001.fastq.gz
# We will link cellID to a reference genome in the "Genome Samplesheet"
# Examples: Nvec_H3K36me3_301121_S27_R1_001.fastq.gz & Nvec_H3K36me3_301121_S27_R2_001.fastq.gz
for i in *fastq.gz; do
    echo $i | sed 's/_R._001.fastq.gz//g' ;
done | sort | uniq > samplesNames.txt

echo "fastq_1,fastq_2,refID" > alignment_samplesheet.csv ;
for i in $(cat samplesNames.txt); do 
    f1=$(find ${i}_R1_001.fastq.gz);
    f2=$(find ${i}_R2_001.fastq.gz);
    refID=(${i/_/ }); refID=${refID[0]};
    echo "${f1},${f2},${refID}";
done >> alignment_samplesheet.csv
```

## 2. Genome Samplesheet
Comma separated file with 2 mandatory columns:  
| Column ID        | Description          | 
| ------------- |:-------------:|
| **refID**      | Cell or species ID | 
| **refGenome**     | Path to reference genome fasta file. Bowtie 2 index llslasdlasldalsd    |  


## 3. Peak calling Samplesheet
Comma separated file with 3 mandatory columns:  
| Column ID        | Description          | 
| ------------- |:-------------:|
| **inPath**      | Path to folder containing BAM files to be analysed | 
| **controlPath**     | Path to control bamfile (input, IgG, H2, etc.) for all samples in folder   |  
| **refID** | ID to be linked with a reference genome in the Genome Samplesheet  |

