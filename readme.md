## Identification and single-base gene-editing functional validation of a cis-EPO variant for use as a genetic proxy for EPO-increasing therapies. ##

### Charli Harlow ###
### charliharlow@outlook.com ###
### March 2022 ###

All data, code and scripts related to the manuscript 'Identification and single-base gene-editing functional validation of a cis-EPO variant for use as a genetic proxy for EPO-increasing therapies'. Scripts include GWAS meta-analysis, RNA-seq analysis of the EPO knock-outs, MR and PheWAS.

## GWAS meta-analysis

GWAS was performed on EPO (after inverse normalisation) in each of the four contributing study cohorts: Health ABC, PREVEND, InCHANTI and BLSA. GWAS was performed using GEMMA. An example script can be seen in EPO_GWAS_meta-analysis_MR/example_gemma_script.sh

```
/programs/gemma -g /scratch/cs660/Inchianti/bimbam/chr$SGE_TASK_ID.bimbam -p /scratch/cs660/Inchianti/Phenotype_and_GRM_files/Phenotype_file_for_GWAS_unix.txt -n 1 -k /scratch/cs660/Inchianti/Phenotype_and_GRM_files/GRM_autosomes.cXX.txt -lmm 4 -o gwas_results_chr22
```

-n option indicates which phenotype column to run. 

Prior to performing GWAS using gemma the vcf files were converted to bimbam files (the input for Gemma) using the perl script EPO_GWAS_meta-analysis_MR/vcf_to_gemma.pl
```
perl vcf_to_gemma.pl --vcf /mnt/Data9/cs660/Inchianti_imputation/chr1.dose.vcf.gz —out /mnt/Data9/cs660/Inchianti_imputation/bimbam/chr1.bimbam >& /mnt/Data9/cs660/Inchianti_imputation/bimbam/chr1.csl 
```
## Meta-analysis
Meta-analysis was performed after quality control checks had been carried out on the individual GWAS' using METAL. The filters applied and script used can be seen in EPO_GWAS_meta-analysis_MR/metal_script.txt

```
metal metal_script.txt
```
Example scripts used for the Manhattan Plot and the QQ plot can be seen in EPO_GWAS_meta-analysis_MR.

## Validation of the _cis-EPO_ SNP
Colocalisation was used to assess whether the _cis-EPO_ SNP was likely the most causal variant in both the hepatic eQTL analysis and the EPO meta-analysis. Colocalisation was performed using the coloc package in R: EPO_GWAS_meta-analysis_MR/coloc.R

A Miami plot was produced for the 500 kb region surrounding the _cis-EPO_ SNP using the script: EPO_GWAS_meta-analysis_MR/miami_plot.R

## RNA-Seq analysis of the whole _EPO_ gene knock-outs

Whole _EPO_ gene knock-outs were generated using CRISPR-Cas9 gene-editing. Whole transcriptomic analysis was then performed on 12 samples (4 X KOA, 4 X KOB and 4 X wild-type controls). The pipeline followed is described in RNAseq_EPOknockouts/README.md and all scripts generated can be found in RNAseq_EPOknockouts. 

## MR analysis 

Single SNP MR analysis was performed to assess the causal association between circulating EPO levels and outcomes of interest using the EPO_GWAS_meta-analysis_MR/single_snp_MR.R script with the _cis-EPO_ SNP as the instrument.

The SNP-exposure association statistics were obtained from EPO meta-analysis. 

The SNP-outcome association statistics were obtained by perfromed a meta-analysis of recently published publicly avaliable GWAS and UK Biobank. 
