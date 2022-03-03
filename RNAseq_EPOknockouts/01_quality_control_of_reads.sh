#!/bin/bash

## Quality Control using FASTQC & MULTIQC
## CHARLI HARLOW

# create directory to store the fastqc output in - one for each sample and each read

for ID in Empty1 Empty2 Empty3 Empty4 KOA1 KOA2 KOA3 KOA5 KOB1 KOB211 KOB3 KOB5
do
mkdir -p /gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/RNA_sequencing_analysis/ftp1.sequencing.exeter.ac.uk/H0253/11_trimmed/additional_trimming/additional_trimming_fastqc/${ID}_r1/
mkdir -p /gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/RNA_sequencing_analysis/ftp1.sequencing.exeter.ac.uk/H0253/11_trimmed/additional_trimming/additional_trimming_fastqc/${ID}_r2/

# Load in module to run FastQC
module load FastQC/0.11.7-Java-1.8.0_162

# Run FastQC
# Read1 
fastqc /gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/RNA_sequencing_analysis/ftp1.sequencing.exeter.ac.uk/H0253/11_trimmed/additional_trimming/${ID}_additional_trimming_r1.fq.gz --extract -o ./additional_trimming_fastqc/${ID}_r1/
# Read2
fastqc /gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/RNA_sequencing_analysis/ftp1.sequencing.exeter.ac.uk/H0253/11_trimmed/additional_trimming/${ID}_additional_trimming_r2.fq.gz --extract -o ./additional_trimming_fastqc/${ID}_r2/

done

# Run MultiQC
# Change directory to faatQC directory
cd ./additional_trimming_fastqc/


# Load in module
module load MultiQC/1.2-intel-2017b-Python-2.7.14

# Run multiQC
multiqc . --dirs --interactive -o ./multiQC/


# Combine all picture files together to make one pdf

# For all the ones which have warnings or have failed
for ID in Empty1 Empty2 Empty3 Empty4 KOA1 KOA2 KOA3 KOA5 KOB1 KOB211 KOB3 KOB5
do
grep -v PASS ${ID}/3064_${ID}_trimmed_r1_fastqc/summary.txt |montage txt:-${ID}/3064_${ID}_trimmed_r1_fastqc/Images/*png -tile x3 -geometry +0.1+0.1 -title ${ID} ${ID}.png done < file_names.txt

convert *png fastqc_summary_warnings.pdf

# For all those which have PASSED
for ID in Empty1 Empty2 Empty3 Empty4 KOA1 KOA2 KOA3 KOA5 KOB1 KOB211 KOB3 KOB5
do
grep PASS ${ID}/3064_${ID}_trimmed_r1_fastqc/summary.txt |montage txt:-${ID}/3064_${ID}_trimmed_r1_fastqc/Images/*png -tile x3 -geometry +0.1+0.1 -title ${ID} ${ID}.png done < file_names.txt

convert *png fastqc_summary_passed.pdf
