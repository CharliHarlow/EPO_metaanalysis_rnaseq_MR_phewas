#!/bin/bash
#MOAB -V
# set verbose output

#MOAB -d .
# set working directory to .

#MOAB -q mrcq
# submit to the mrcq

#MOAB -l walltime=7:00:00:00
# set the maximum wallclock time for the job. Don't just set this to a very large number or you risk your job being terminated!

#MOAB -A Research_Project-MRC158833
# research project to submit under. I don't know exactly how this flag works at the moment

#MOAB -l procs=16
# either this to specify number of processors. This will distribute over nodes as MOAB sees fit. Not needed if only 1 slot required

# #MOAB -l nodes=1:ppn=13
# or this, nodes=number of nodes required. ppn=number of processors per node. DO NOT USE UNLESS YOU WANT TO RESERVE A WHOLE NODE

#MOAB -m e -M cs660@exeter.ac.uk
# email me at job completion with stats

#MOAB -l naccesspolicy=SHARED
# allow other jobs to use the same node

# create alias for the directory with reference genome
Ref_dir="/gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/RNA_sequencing_analysis/Ref_genome"

# load in STAR module
module load STAR/2.7.1a-foss-2018b

# create directory for Star Index files
# mkdir /gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/RNA_sequencing_analysis/StarIndex

# Generate STAR index files

# STAR --runMode genomeGenerate --genomeDir /gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/RNA_sequencing_analysis/StarIndex/ --genomeFastaFiles ${Ref_dir}/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile ${Ref_dir}/ensembl/Homo_sapiens.GRCh38.98.gtf --sjdbOverhang 74

# Run alignment for each sample. If Single End then need to change some options and input names
mkdir /gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/RNA_sequencing_analysis/Star_Alignment_additional_trimming/

index_dir="/gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/RNA_sequencing_analysis/StarIndex"

for SAMPLE in Empty1 Empty2 Empty3 Empty4 KOA1 KOA2 KOA3 KOA5 KOB1 KOB211 KOB3 KOB5
do
STAR --runMode alignReads --genomeDir ${index_dir}/ --readFilesIn /gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/RNA_sequencing_analysis/ftp1.sequencing.exeter.ac.uk/H0253/11_trimmed/3064_${SAMPLE}_trimmed_r1.fq.gz /gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/RNA_sequencing_analysis/ftp1.sequencing.exeter.ac.uk/H0253/11_trimmed/3064_${SAMPLE}_trimmed_r2.fq.gz --readFilesCommand zcat --outFileNamePrefix /gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/RNA_sequencing_analysis/Star_Alignment_default/test_2/${SAMPLE}_alignment_default --outSAMtype BAM SortedByCoordinate --runThreadN 4 
done
