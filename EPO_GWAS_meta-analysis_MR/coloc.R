## Colocalisation of eQTL data and EPO meta-analysis
## Charli Harlow

## Carried out using Coloc package ##

library(coloc)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("snpStats")

if(requireNamespace("snpStats")) {
  library(snpStats)}

data <- read.table("~/Documents/EPO Project/Colocalisation/eQTL_and_GWAS_500kb_rs1617640.txt", header=T, stringsAsFactors=F)

# To run the command below, you need the two betas, two ses (which are squared) and the N.
# For one dataset I provide the SD of the phenotype tested as I know it – but for other I do not - so coloc tries to estimate this from the data provided.
# Note also it doesn’t matter whether you use maf or effect allele frequency but given the documentation I tried to be consistent with it when using the command

my.res <- coloc.abf(dataset1=list(beta=data$beta_acr, varbeta=data$se_acr^2, sdY=1, N=6127, type="quant"),
                    dataset2=list(beta=data$beta_fg,  varbeta=data$se_acr^2,        N=861,  type="quant"),
                    MAF=data$maf)
## Note:
## H0: neither trait has a genetic association in the region
## H1: only trait 1 has a genetic association in the region
## H2: only trait 2 has a genetic association in the region
## H3: both traits are associated, but with different causal variants
## H4: both traits are associated and share a single causal variant

# If you did a p-value analyses you should get similar results (unless N is very small). Above methods described as being slightly better and should be used if you have the betas and SES.

my.res_pval <- coloc.abf(dataset1=list(pvalues=data$p_gwas, sdY=1, N=6127, type="quant"),
                         dataset2=list(pvalues=data$p_eqtl, N=861,  type="quant"),
                         MAF=data$maf)
