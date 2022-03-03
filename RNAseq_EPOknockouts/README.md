# RNA Sequencing Analysis Pipeline for the analysis of _EPO_ knock-outs in HEK-293 cells #
### Charli (Stoneman) Harlow ###
### charliharlow@outlook.com
### April 2020 ###
### University of Exeter ###

## RNA Sequencing Analysis
Follow this document to carry out basic RNA sequencing Analysis.

In this experiment, I performed RNA sequencing on 12 samples; 4 X WT *(WT1-4)* (cells treated with empty vector), 4 X KOA *(KOA1-4)* (each KOA sample has come from the same single cell which was clonally expanded), 4 X KOB  *(KOB1-4)* (each KOB sample has come from the same single cell which was clonally expanded)

## Download sequencing data from Exeter Sequencing Service
Once you have received email from Exeter Sequencing Service, you can download the results onto server using the following command
```
wget -r -p ftp://Username:Password@ftp1.sequencing.exeter.ac.uk 
```
•	To run this in the background so if you log off it still continues to download add in the -b option. This will produce a wget-log file with information on where it is in the downloading. 

•	To look into the wget-log file use: less wget-log

## Quality Control of the raw reads

1. QC analysis was performed on the RNA raw data using FastQC and then MultiQC to compile all quality control assessments on each sample into a single html document. 

FastQC 
```
#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=24:00:00 # maximum walltime for the job e.g 01:00:00
#SBATCH -A Research_Project-MRC158833 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mem=100G # specify bytes of memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=cs660@exeter.ac.uk # email address

# Make directory for fastqc results to go into

for ID in Empty1 Empty2 Empty3 Empty4 KOA1 KOA2 KOA3 KOA5 KOB1 KOB211 KOB3 KOB5
do
mkdir -p /gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/RNA_sequencing_analysis/ftp1.sequencing.exeter.ac.uk/H0253/01_raw_reads/fastqc/${ID}_raw_reads_r1/
mkdir -p /gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/RNA_sequencing_analysis/ftp1.sequencing.exeter.ac.uk/H0253/01_raw_reads/fastqc/${ID}_raw_reads_r2/

# Change directory and name of input files and output files

#load in fast qc module
module load FastQC/0.11.7-Java-1.8.0_162
# run fast qc
# change input path to directory

fastqc /gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/RNA_sequencing_analysis/ftp1.sequencing.exeter.ac.uk/H0253/01_raw_reads/3064_${ID}_r1.fq.gz --extract -o ./fastqc/${ID}_raw_reads_r1/
fastqc /gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/RNA_sequencing_analysis/ftp1.sequencing.exeter.ac.uk/H0253/01_raw_reads/3064_${ID}_r2.fq.gz --extract -o ./fastqc/${ID}_raw_reads_r2/
done


# move into fast qc directory
cd /fastqc/

# run MultiQC to combine alll fast qc files

module load MultiQC/1.2-intel-2017b-Python-2.7.14

multiqc . --dirs --interactive -o ./multiQC_rawreads/
```

2. Open up the html file (e.g. `raw-reads_multiqc_report.html`) produced and view your QC results.

- In the QC file: under the General Statistics section, check that % Dups is around 50%, %GC is around 50% and length is similar across samples. The total number of sequences for each sample can also be seen in the final column. 

- Under FastQC: Sequence Quality Histograms check that all bases have a Phred score > 28 so are in the green section of the graph. If the first 5 bases are not some additional trimming may be needed.

- Per Sequence GC content should have a normal distribution with a peak just under 50%

- Per Base N Content should be very close to 0

- It is fine if you have a sudden rise at the end of the sequence length distribution plot as samples have been trimmed by sequencing service to be the same length.

- A little peak in sequence duplications peak is seen commonly so do not be concerned about this

- There should be no over-represented sequences

- There should be no adapter content as these should have already been removed by Cutadapt by the sequencing team.

## Trimming using CutAdapt

Trimming was performed to remove low quality bases (<Q22) from the 3' end and any adapter sequences from the reads. Adapter sequences were defined here: <http://emea.support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html>. Reads shorter than 25 were also discarded. 

Cutadapt version 1.13 was used

``` 
cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --minimum-length 25 \
    -q 22 \ 
    -o trimmed.R1.fastq.gz -p trimmed.R2.fastq.gz \
    reads.R1.fastq.gz reads.R2.fastq.gz
``` 

## 1. Download the raw reads from Exeter Sequencing service onto server ##
Once you have received email from Exeter Sequencing Service, you can download the results onto server using the following command
```
wget -r -p ftp://Username:Password@ftp1.sequencing.exeter.ac.uk 
```
•	To run this in the background so if you log off it still continues to download add in the -b option. This will produce a wget-log file with information on where it is in the downloading. 

•	To look into the wget-log file use: less wget-log

## Repeat FastQC and MultiQC on the trimmed reads to check trimming

Repeat FastQC and MultiQC as above changing the names of the input files to the trimmed fastq files. 

## Additional Trimming Using Trimmomatic 
If additional trimming is required, trimmomatic can be used, for examples sometimes the first 5 bases may need to be removed if they have a poor Phred score (<28). *More details on the software can be found here: <http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf>*

1. Download Trimmomatic onto the server and set up for running Trimmomatic

```
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
# load in the correct java module
module load java/1.8.0_92
```

2. Run Trimmomatic
For details on the additional options, read the Trimmomatic Manual: <http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf>

```
java -jar ../Trimmomatic-0.39/trimmomatic-0.39.jar PE \ #specify paired-end or single-end
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
```

## Quality control check of the additional trimming
Using FastQC or MultiQC again as described above. See script: _**02_quality_control_of_reads.sh**_

Picture files produced by FastQC can be combined into on pdf file using the following commands: 

For all the ones that have **PASSED** the QC
```
cd /fastqc/

for ID in Empty1 Empty2 Empty3 Empty4 KOA1 KOA2 KOA3 KOA4 KOB1 KOB2 KOB3 KOB4
do
grep PASS ${ID}_fastqc/summary.txt | montage txt:- ${ID}_fastqc/Images/*png -tile x3 -geometry +0.1+0.1 -title ${ID} ${ID}_pass.png
done
#combine images
convert *_pass.png fastqc_summary_pass.pdf
```

For all the ones that have **NOT PASSED** *(either warning or failed)* the QC
```
cd /fastqc/

for ID in Empty1 Empty2 Empty3 Empty4 KOA1 KOA2 KOA3 KOA4 KOB1 KOB2 KOB3 KOB4
do
grep -v PASS ${ID}_fastqc/summary.txt | montage txt:- ${ID}_fastqc/Images/*png -tile x3 -geometry +0.1+0.1 -title ${ID} ${ID}_warnings.png
done
#combine images
convert *_warnings.png fastqc_summary_warning.pdf
```
## Alignment of the reads to the reference genome
Alignment was carried out using the STAR software: <https://github.com/alexdobin/STAR/raw/master/doc/STARmanual.pdf> 

1. Downloading the reference genome
- Create directory for reference genome to be stored
```
mkdir ~/Ref_genome
```

- Download the fasta files and gzip files for the reference genome
```
wget ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz

# unzip both for use in STAR
gunzip Homo_sapiens.GRCh38.98.gtf.gz
gunzip  Homo_sapiens.GRCh38.98.gtf.gz
```

2. Generate the Index file
Only one index file needs to be created per reference genome. The index file will contain all the information from the reference genome in a compressed format that is optimized for efficient access and comparison with the query read sequences. The main input files for this step therefore encompass the reference genome sequence and an annotation file.

``` 
mkdir StarIndex

Ref_dir="/gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/RNA_sequencing_analysis/Ref_genome/ensembl/"

module load STAR/2.7.1a-foss-2018bSTAR --runMode genomeGenerate \ # tells STAR to generate index
--genomeDir /gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/RNA_sequencing_analysis/StarIndex/ \ #where to store the results
--genomeFastaFiles ${Ref_dir}/GRCh38.p13.genome.fa \ # make sure this file is unzipped
--sjdbGTFfile ${Ref_dir}/gencode.v32.chr_patch_hapl_scaff.annotation.gtf \ # make sure this file is unzipped
--sjdbOverhang 100 \ # specifies the length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database. Ideally, this length should be equal to the ReadLength-1, where ReadLength is the length of the reads.

```

3. Run the mapping of the fastq files to the reference genome
*Make sure when reading files in if you are mapping paired end data then read1 and read2 and space separated NOT separated by comma.* 


```
mkdir alignment_STAR
ref_dir="/gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/RNA_sequencing_analysis/StarIndex"
module load STAR/2.7.1a-foss-2018b
for SAMPLE in Empty1 Empty2 Empty3 Empty4 KOA1 KOA2 KOA3 KOA4 KOB1 KOB2 KOB3 KOB4
do
STAR --runMode alignReads \
--genomeDir ${ref_dir}/ \ #directory containing the index files
--readFilesIn ${SAMPLE}_1P.fq.gz ${SAMPLE}_2P.fq.gz \ #read in the fast q files for aligning
--readFilesCommand zcat \ #necessary because of gzipped fastq files
--outFileNamePrefix /gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/RNA_sequencing_analysis/Alignment/${SAMPLE}_ \ #output files name prefix (including full or relative path). Can only be defined on the command line.
--outSAMtype BAM SortedByCoordinate \ # output in BAM format, sorted by coordinate. This option will allocate extra memory for sorting which can be specified by –limitBAMsortRAM
--outReadsUnmapped Fastx \ # default is none, fastx will output unmapped and partially mapped reads in separate files
--runThreadN 4 \
--twopassMode Basic \ # STAR will perform mapping, then extract novel junctions which will be inserted into the genome index which will then be used to re-map all reads
--outFilterMultimapNmax 1 \ # only reads with 1 match in the reference will be returned as aligned(
done
```


Check alignment scores in the .log.final.out files. Looking for:
- Uniquely mapped reads > 80%
- Reads mapped to multiple loci 
- Reads unmapped
- Average mapping length should be the length of the sequence

Can also run MultiQC on the alignment files
``` 
module load MultiQC/1.2-intel-2017b-Python-2.7.14	

multiqc ./alignment_STAR/ --dirs --interactive -o ./multiQC/

```

##  Visualising the STAR Alignments

This was carried out in R Studio using the script _** Visualising_alignments_script.R** _ to produce bar plots summarising the Star Alignment results. 

## Creating index files for the bam files
Most downstream tools require an index file (.bam.bai) of the bam files produced after alignment.

SAMtools was used to create this file

```
module load SAMtools 
for ID in Empty1 Empty2 Empty3 Empty4 KOA1 KOA2 KOA3 KOA4 KOB1 KOB2 KOB3 KOB4
do 
samtools index ${ID}_Aligned.sortedByCoord.out.bam
done
```

## Viewing alignment files in IGV
To do this, can either load each individual bam file into IGV and view the region you desire or you can merge together the bam files for each cell line. To view in IGV. an index file is required

1. Merge .bam files into one
```
module load SAMtools
samtools merge -o WT_merged.bam WT1_Aligned.sortedByCoord.out.bam WT2_Aligned.sortedByCoord.out.bam WT3_Aligned.sortedByCoord.out.bam WT4_Aligned.sortedByCoord.out.bam
samtools index WT_merged.bam
```

2. View in IGV. 
- For more information on IGV and how to download,look here:
<http://software.broadinstitute.org/software/igv/download>. There are really good online tutorials into how to best use IGV.

a. Load the corresponding genome e.g. hg19 from the left hand side drop down list on the upper left of the IGV window
b.	Load in the bam file which you want to view
  i. Load from file
  ii.	Select bam file. IGV will then automatically detect and look for the index file of this bam file which we created earlier. 
  iii.	Loading a BAM file creates up to 3 associated tracks:  
          1.	Alignment Track to view individual aligned reads
          2.	Coverage Track to view depth of coverage
          3.	Splice Junction Track which provides an alternative view of reads spanning splice junctions
  iv.	You can then zoom into the gene of interest by either typing the gene name into the search bar or by entering the coordinates of the gene
  v.	Investigate whether the region expected has been deleted i.e no reads are mapping to this area
  vi.	It is also important to have a look at a few housekeeping genes to ensure alignment has worked and that there is good coverage and read depth across the genes.

## Quality Control of the Alignment files using QoRTs
This step is often not carried out if the alignment stats look good but can be done for sanity check and to get more figures that can be used to show alignment looks good. 

A suitable R environment is required to generate plots using QoRTs:
1.	Navigate to your packages directory and input the following commands to download the required files (using Anaconda here is useful as administrator privileges are absent on the HPC UNIX environment)

```
wget https://repo.continuum.io/archive/Anaconda2-5.3.1-Linux-x86_64.sh
wget https://repo.anaconda.com/archive/Anaconda3-2019.07-Linux-x86_64.sh
# Unpack the files:
bash Anaconda2-5.3.1-Linux-x86_64.sh
bash Anaconda3-2019.07-Linux-x86_64.sh
# Declare the variable:
export PATH=/home/ubuntu/software/anaconda2bin:$PATH
# Activate the Conda environment
Conda activate
# Install and create a new R environment via the installed Conda environment:
conda create -n r_env r-essentials r-base
conda activate r_env
# Install the Java dev-kit JRE in conda:
conda install -c cyclus java-jdk
# Install R-QoRTs in Conda:
conda install -c bioconda r-qorts
# Install Bioconductor-deseq2 in Conda:
conda install -c bioconda bioconductor-deseq2
# Check that the packages installed:
Conda list
# Navigate back to the packages directory and download the QoRTs package:
wget -O qorts.jar https://github.com/hartleys/QoRTs/archive/v1.3.6.tar.gz
```
2. Once suitable environment has been made, QoRTs can be run
- Load in Anaconda modules
``` 
module load Anaconda3/5.2.0
```

-	Activate the R environment

```
source activate r_env
```

- Load in required Java module for QoRTs
```
module load Java/1.8.0_92
```
- Run QoRTs command
```
for ID in Empty1 Empty2 Empty3 Empty4 KOA1 KOA2 KOA3 KOA4 KOB1 KOB2 KOB3 KOB4
do
java -Xmx4g -jar ./programme/QoRTs.jar QC --stranded --maxReadLength 75 --generatePdfReport ./Alignment/${ID}_Aligned.sortedByCoord.out.bam [path/to/directory/ref_genome]/Homo_sapiens.GRCh38.98.gtf [path/to/output/directory/]
```

3. Review the quality of the aligned reads in the PDF file found in the QoRTs output directory

## Gene Quantification

There are two ways to do this;

1.	Assign all the reads to a given gene
2.	Infer the quantity of individual transcripts 

There are also several different programs which can be used to carry out gene quantification. For this purpose, to count the number of reads for each gene, I used the featureCounts package. 

1.	Download the feature counts package which is present in the subread form: <https://sourceforge.net>
2.	
```
wget -r https://sourceforge.net/projects/subread/files/subread-2.0.0/subread-2.0.0-Linux-x86_64.tar.gz/download
```
2.	Uncompress and unpackage the .tar.gz file
```
tar -zxvf download
```
3.	Create directory to put feature count results in
```
mkdir featureCounts
```
4.	If carrying out featureCounts on paired-end data, ensure that the .bam files are sorted by read name NOT BY coordinate like STAR output files
*If not sorted by read name, feature counts will assume that almost all reads are not properly paired*
```
module load SAMtools
samtools sort -n -o /path/to/directory/output_name_sortByReadName.bam /path/to/directory/input_name.out.bam
```

5.	Run feature counts to count reads per gene
```
# Set up alias for feature counts directory
featureCounts=~/path/to/package/directory/subread-2.0.0-Linux-x86_64/bin/featureCounts
# Run feature counts (if you have not set up alias as above command states,make sure to supply full path-to-directory/featureCounts at start of below command to run featureCounts
$featureCounts -a /path/to/directory/to/reference-genome/Homo_sapiens.GRCh38.98.gtf -T 8 -p -g gene_id -F GTF -o /path/to/directory/output/featureCounts/feature_count_results.txt /path/to/directory/input_file_name.bam 2> /path/to/directory/for/log/file/featurecounts_screen_output.log
# If you want to run for all bam files then use * for wildcard e.g /path/to/directory/input_files/*.bam
```

Options:
- -T indicates number of threads to use
-	-p indicates paired-end
-	-F indicates type of reference file
-	-f indicates what level to perform the assignment at – default is to perform assignment at gene-level (meta-feature). If you specify -f then it will perform quantification at exon-level (feature level)
-	-g indicates how to name the genes. Default is gene_id. Can change this to gene_name if you want output to contain the gene names not the accession numbers. Often it is better to use accession numbers as genes can have more than one gene-name
-	2> indicates to make a log file showing the screen output as featurecounts is running


featureCounts outputs several different files;
      1.	Count table
*This includes annotation columns (`Geneid', `Chr', `Start', `End', `Strand' and `Length') and data columns (read counts for genes for each library)* 
      2.	Counting summary file 
*This file includes information on the total number of alignments that were successfully assigned to genes and the number of alignments which failed to be assigned due to multiple reasons.*
*Note that the counting summary includes the number of alignments, not the number of reads. Number of alignments will be higher than the number of reads when multi-mapping reads are included since each multi-mapping read contains more than one alignment.*
      3.	Screen ouput
*Contains information similar to the summary file with number and percentages of successfully assigned alignments*

## Visualising featureCounts results

It is good practice to plot the statistics produced from featureCounts so you can assess if the quantification worked and can compare the number of alignments that have been assigned vs those that have not. To do this see the script **_Visualising_featureCounts.R_**

## Differential Gene Analysis

This was carried out in R using the DeSeq2 package. For run through example and R script, see Differential_Gene_Expression_Analysis.html

## Gene Ontology Analysis

Once you have a list of differentially expressed genes, gene ontology enrichment analysis can be carried out using DAVID (https://david.ncifcrf.gov ) or Enrichr (https://amp.pharm.mssm.edu/Enrichr/) – both give similar results. Just plug in a list of gene names to both online tools.

The main categories to look for enriched terms if KEGG pathways, GO biological pathways, GO molecular function, GO Cellular Component, OMIM Disease. 

I also performed Gene Ontology Analysis using GoSeq2 for comparison purposes. How this was run can be found in the script _**GOseq2.R**_


