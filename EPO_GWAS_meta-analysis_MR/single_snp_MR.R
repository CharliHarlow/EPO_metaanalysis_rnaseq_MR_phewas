#ÂÂ# Single SNMR using the cis-EPO variant
## Charli Harlow

## Read in package
library(TwoSampleMR)

## read in exposure data (genotype-exposure)
exp <- read.csv("~/Documents/EPO_Consortium/Analysis of results/MRBase/epo.csv", header=T)
random_exp_dat <- format_data(exp, type="exposure")

## read in outcome data (Genotype-outcome)
## For binary traits make sure it is LOG(ODDS) not OR
outcome <- read_outcome_data(
  snps = exp$SNP,
  filename = "~/Documents/EPO_Consortium/Analysis of results/MRBase/Stroke.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "ln(OR)",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  phenotype_col = "Phenotype",
  ncase_col = "ncase",
  ncontrol_col = "ncontrol")

## make sure alleles are flipped the same way

dat <- harmonise_data(
  exposure_dat = random_exp_dat, 
  outcome_dat = outcome
)  

## Perform MR using single snp instrument - Wald ratio estimator

res <- mr(dat)
res1 <- generate_odds_ratios(res)
write.table(res1, "~/Documents/EPO_Consortium/Analysis of results/MRBase/epo-stroke_mr_results.csv", sep=",", row.names=F)
res1 <- mr_singlesnp(dat)

## For continuous traits

outcome <- read_outcome_data(
  snps = exp$SNP,
  filename = "~/Documents/EPO_Consortium/Analysis of results/MRBase/SBP.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "effect",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  phenotype_col = "Phenotype",
  samplesize_col = "n")

## make sure alleles are flipped the same way

dat <- harmonise_data(
  exposure_dat = random_exp_dat, 
  outcome_dat = outcome
)  

## Perform MR using single snp instrument - Wald ratio estimator

res <- mr(dat)

res1 <- mr_singlesnp(dat)
