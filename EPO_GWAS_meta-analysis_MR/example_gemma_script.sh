#!/bin/bash

# Charli Harlow
Î# Example Gemma script used for GWAS of EPO 


# Use the following script template to submit scripts running only a single thread (i.e. using a single core on a single node)
#

#The -V switch means that the same environment variables as currently set will be used on the nodes
#The -cwd switch means run in the current directory
#The -q serial.q switch means run in the serial queue

#$ -V -cwd -q serial.q
. /etc/profile.d/modules.sh    # Do not delete this line if you have difficulty loading modules

#The following line will kill the job automatically after 72 hrs and 10 mins - please contact the sysadmin before changing
#$ -l h_rt=72:10:00


#Uncomment the following line and enter a valid email address to receive email notification of job completion

#$ -m e -M cs660@exeter.ac.uk
echo Running on ; hostname

#Place your command lines below this
#$ -t 1-22

/users/cs660/programs/gemma -g /scratch/cs660/Inchianti/bimbam/chr$SGE_TASK_ID.bimbam -p /scratch/cs660/Inchianti/Phenotype_and_GRM_files/Phenotype_file_for_GWAS_unix.txt -n 2 -k /scratch/cs660/Inchianti/Phenotype_and_GRM_files/GRM_autosomes.cXX.txt -lmm 4 -o epo_raw_gwas_results_chr$SGE_TASK_ID

