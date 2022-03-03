#!/bin/bash

## Trimming using Trimmomatic
## Use this script to carry out additional trimming using trimmomatic if required. 

## CHARLI HARLOW

# Download Trimmomatic

wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip

# Unzip Trimmomatic
unzip Trimmomatic-0.39.zip

# Load in java module needed for Trimmomatic
module load java/1.8.0_92

# may need to edit depending on settings required
# For more details about the trimmomatic package please see http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf

# Paired-End
for ID in Empty1 Empty2 Empty3 Empty4 KOA1 KOA2 KOA3 KOA5 KOB1 KOB211 KOB3 KOB5
do
java -jar ../Trimmomatic-0.39/trimmomatic-0.39.jar 
PE \ #specify paired-end or single-end
-phred33 #specifies the base quality encoding. Could be changed to phred64 if using older sequencing machines
-trimlog trim_empty1_101219.log \ # creates a log file of all the trimming
-basein /gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/EPO_RNA_Seq/ftp1.sequencing.exeter.ac.uk/H0253/11_trimmed/3064_Empty1_trimmed_r1.fq.gz \ # specifies the input file name. Use this option if all input names have the same common naming pattern so the reverse read file can automatically be detected. If not then remove this option and just specify the two files e.g input_filename_r1.fq.gz input_filename_r2.fq.gz
-baseout 3064_Empty1_trimmed_r1_filtered.fq.gz \ # specifies the output file names. Four files will be produced 1-Paired forward reads, 2-Unpaired forward reads 3-Paired reverse reads 4-Unpaired reverse reads
ILLUMINACLIP:/gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/EPO_RNA_Seq/ftp1.sequencing.exeter.ac.uk/H0253/11_trimmed/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads \
HEADCROP:10 \ #Cut the specified number of bases from the start of the read. Particular important to add this step in if looking at the QC report and realising that the first 10 bases for example have a phred score lower than
LEADING:3 \ #Cut bases off the start of a read, if below a threshold quality
TRAILING:3 \ # Cut bases off the end of a read, if below a threshold quality  
SLIDINGWINDOW:4:15 \ #Performs a sliding window trimming approach. It starts scanning at the 5’ end and clips the read once the average quality within the window falls below a threshold. 
MINLEN:36 \ # Drop the read if it is below a specified length

# Single-End
# Remove the hashtags if wanting to run for single ended
# for ID in Empty1 Empty2 Empty3 Empty4 KOA1 KOA2 KOA3 KOA5 KOB1 KOB211 KOB3 KOB5
# do
# java -jar /gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/EPO_RNA_Seq/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 -trimlog trim_{$ID}_101219.log /gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/EPO_RNA_Seq/ftp1.sequencing.exeter.ac.uk/H0253/11_trimmed/3064_{$ID}_trimmed_r1.fq.gz 3064_{$ID}_trimmed_r1_extratrimming.fq.gz ILLUMINACLIP:/gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/EPO_RNA_Seq/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:10 SLIDINGWINDOW:4:15 MINLEN:36
# done


