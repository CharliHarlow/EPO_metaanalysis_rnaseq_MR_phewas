### Gene Ontology Analysis using GoSeq2 ###
### Charli (Stoneman) Harlow ###


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("goseq")
BiocManager::install("org.Hs.eg.db")
BiocManager::install(c("rlang","TxDb.Hsapiens.UCSC.hg38.knownGene"))

library(goseq)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
## Goseq2 analysis
## https://bioconductor.org/packages/devel/bioc/vignettes/goseq/inst/doc/goseq.pdf
setwd("~/Documents/EPO Project/CRISPR/Whole gene knock-out/Confirming EPO KO/RNA Sequencing/DeSeq2/Overlapping DEG analysis/GEO analysis/GO")
KOA_KOB_values <- read.table("~/Documents/EPO Project/CRISPR/Whole gene knock-out/Confirming EPO KO/RNA Sequencing/DeSeq2/Overlapping DEG analysis/overlapping_consistent_DEGs_p05.txt", sep="\t", header=T)
### Using the overlapping genes with FDR P<0.05
KOA_KOB_values_2 <- read.table("~/Documents/EPO Project/CRISPR/Whole gene knock-out/Confirming EPO KO/RNA Sequencing/DeSeq2/Overlapping DEG analysis/Overlapping_genes_DEGs_p05_KOA_KOB_analysis.txt", sep="\t", header=T)
KOA_KOB_values$DE <- 1
de.genes <- KOA_KOB_values$Geneid
de.genes
## Read in all genes assayed
gtf <- read.table("~/Downloads/Homo_sapiens.GRCh38.98_gene_annotation_table.txt", sep="\t", header=T)
gtf$Geneid <- gtf$gene_id

library(biomaRt)
listMarts()    # to see which database options are present
ensembl=useMart("ensembl")  # using ensembl database data
listDatasets(ensembl)     # function to see which datasets are present in ensembl
ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)   # from ensembl using homosapien gene data
listFilters(ensembl)  # check which filters are available
listAttributes(ensembl) # check attributes are available to select.More information on ensembl data base
genes.with.id=getBM(attributes=c("ensembl_gene_id", "entrezgene_id"), mart= ensembl) # fuction to get  gene id's and gene name from data base

genes.with.id$Geneid <- genes.with.id$ensembl_gene_id

gtf_entrez <- merge(gtf, genes.with.id, by="Geneid")

## Load in the count data/file with list of gene names
read.counts <- read.table("~/Documents/EPO Project/CRISPR/Whole gene knock-out/Confirming EPO KO/RNA Sequencing/FeatureCounts/feature_count_results_geneid.txt", header=T)
read.counts.goseq <- read.table("iCloud Drive (Archive)/Documents/EPO Project/CRISPR/Whole gene knock-out/Confirming EPO KO/RNA Sequencing/FeatureCounts/feature_count_results_geneid.txt", header = TRUE)
read.counts.goseq <- read.counts

#replace all row names with the names of genes
row.names(read.counts.goseq) <- read.counts.goseq$Geneid
#remove the columns which contain no count data
read.counts.goseq <- read.counts.goseq[,-c(2:6)]

# give meaningful sample names - this can be achieved via numerous approaches
names(read.counts.goseq) <- c("Geneid", "WT1", "WT2","WT3","WT4","KOA1", "KOA2", "KOA3", "KOA4", "KOB1", "KOB2", "KOB3", "KOB4")


# Check data is what we expect
str(read.counts.goseq)
head(read.counts.goseq, n = 3)

## Merge the gtf annotation data frame with the read counts file by geneid to obtain gene_names as well
read.counts.goseq.genename <- merge(read.counts.goseq, gtf, by="Geneid")

## make a list of all gene assaysed
all_genes <- read.counts.goseq.genename$gene_id

## create a gene list for deseq2, 1=DE, 0=no DE
gene.vector=as.integer(all_genes%in%de.genes)
names(gene.vector)=all_genes
head(gene.vector)
table(gene.vector)

## see which organisms are supported 
supportedOrganisms()
supportedOrganisms()[supportedOrganisms()$Genome=="hg19",]

## obtain data for hg38
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicFeatures")

install.packages("RMariaDB")
library(RMariaDB)

library(GenomicFeatures)

## Get from USCS
## The function mUCSC downloads UCSC Genome Bioinformatics transcript tables (e.g.   Gene", "refGene", "ensGene") for a genome build (e.g. "mm9", "hg19"). Use the supportedUCSCtables utility function to get the list of tables known to work with makeTxDbFromUCSC.
supportedUCSCtables(genome="hg38")

hg38_txdb <- makeTxDbFromUCSC(genome="hg38", tablename="knownGeneOld8")

#hg38
# ensGene

#Fitting the Probability Weighting Function (PWF)
pwf=nullp(gene.vector,"hg19", "ensGene")
head(pwf)


##Using the Wallenius approximation
GO.wall=goseq(pwf,"hg19", "ensGene")
head(GO.wall)

##Using random sampling
## It may sometimes be desirable to use random sampling to generate the null distribution for category membership. For example, to check consistency against results from the Wallenius approximation. This is easily accomplished by using the method option to specify sampling and the repcnt option to specify the number of samples to generate:
GO.samp=goseq(pwf,"hg19","ensGene",method="Sampling",repcnt=1000)
head(GO.samp)

# plot p-values against each other
library(ggplot2)
plot(log10(GO.wall[,2]), log10(GO.samp[match(GO.wall[,1],GO.samp[,1]),2]),xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)",xlim=c(-3,0))
abline(0,1,col=3,lty=2)

# By default, goseq tests all three major Gene Ontology branches; Cellular Components, Biological Processes and Molecular Functions. 
# limit what is analysed 
# Molecular Functions
GO.MF=goseq(pwf,"hg19","ensGene",test.cats=c("GO:MF"))
head(GO.MF)
enriched_MF <- GO.MF[p.adjust(GO.MF$over_represented_pvalue, method="BH")<.05, ]
library(ggplot2)

ggplot(enriched_MF, aes(y=term, x=-log10(over_represented_pvalue))) +
  geom_bar(stat="identity") + theme_classic() + xlab("-log10 P-value") + ylab("Molecular function")


#Biological processes 
GO.BP=goseq(pwf,"hg19","ensGene",test.cats=c("GO:BP"))
head(GO.BP)
enriched_BP <- GO.BP[p.adjust(GO.BP$over_represented_pvalue, method="BH")<.05, ]

top101 <- GO.BP[1:15, ]

library(ggplot2)
library(dplyr)

top101 <- mutate(top101, cat1 = paste(term,"(",category,")")) 
#top101$cat <- paste(top101$term,"(",top101$category,")")
top101$cat <- factor(top101$cat, levels = top101$cat[order(top101$over_represented_pvalue, decreasing=TRUE)])

ggplot(top101, aes(y=cat, x=-log10(over_represented_pvalue))) +
  geom_bar(stat="identity", aes(fill = over_represented_pvalue), show.legend = F) + 
  theme_classic() + xlab("-log10 P-value") + ylab("Biological Pathway") + 
  geom_text(aes(label=cat),size=3.5, hjust=1, color="white") + theme(
                                                                       axis.text.y = element_blank(),
                                                                       axis.ticks.y = element_blank())
top101_MF <- GO.MF[1:15, ]

top101_MF <- mutate(top101_MF, cat = paste(term,"(",category,")")) 
top101_MF$cat <- factor(top101_MF$cat, levels = top101_MF$cat[order(top101_MF$over_represented_pvalue, decreasing=TRUE)])

ggplot(top101_MF, aes(y=cat, x=-log10(over_represented_pvalue))) +
  geom_bar(stat="identity", aes(fill = over_represented_pvalue), show.legend = F) + 
  theme_classic() + xlab("-log10 P-value") + ylab("Molecular Function") + 
  geom_text(aes(label=cat),size=3.5, hjust=1, color="white") + theme(
                                                                     axis.text.y = element_blank(),
                                                                     axis.ticks.y = element_blank())
#Making sense of the results
#Having performed the GO analysis, you may now wish to interpret the results. If you wish to identify categories significantly enriched/unenriched below some p-value cutoff, it is necessary to first apply some kind of multiple hypothesis testing correction. For example, GO categories over enriched using a .05 FDR cutoff [Benjamini and Hochberg, 1995] are:
enriched <- GO.wall[p.adjust(GO.wall$over_represented_pvalue, method="BH")<.05, ]

write.table(enriched, "3501_overlapping_DEGs_GO_terms.txt", sep="/t", row.names = F, quote=F)

enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue,method="BH")<.05]
head(enriched.GO)

library(GO.db)

for(go in enriched.GO[1:263]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}



## KEGG Pathway

# Get the mapping from ENSEMBL 2 Entrez
en2eg=as.list(org.Hs.egENSEMBL2EG)
# Get the mapping from Entrez 2 KEGG
eg2kegg=as.list(org.Hs.egPATH)
# Define a function which gets all unique KEGG IDs
# associated with a set of Entrez IDs
grepKEGG=function(id,mapkeys){unique(unlist(mapkeys[id],use.names=FALSE))} # Apply this function to every entry in the mapping from
# ENSEMBL 2 Entrez to combine the two maps
kegg=lapply(en2eg,grepKEGG,eg2kegg)
head(kegg)

pwf=nullp(gene.vector,"hg19","ensGene")
KEGG=goseq(pwf,gene2cat=kegg)
head(KEGG)
KEGG1=goseq(pwf,'hg19','ensGene',test.cats="KEGG")
head(KEGG1)

## COrrecting for count bias too
counts <- read.counts.goseq[,-c(1)]
counts <- as.data.frame(lapply(counts, as.numeric))
countbias=rowSums(counts)
length(countbias)
length(gene.vector)

pwf.counts=nullp(gene.vector,bias.data=countbias)
GO.counts=goseq(pwf.counts,"hg19","ensGene")
head(GO.counts)


### Repeating analysis using the consistent DEGs
consistent$DE <- 1
de.genes <- consistent$Geneid

## create a gene list for deseq2, 1=DE, 0=no DE
gene.vector=as.integer(all_genes%in%de.genes)
names(gene.vector)=all_genes
head(gene.vector)
table(gene.vector)

#Fitting the Probability Weighting Function (PWF)
pwf=nullp(gene.vector,"hg19","ensGene")
head(pwf)


##Using the Wallenius approximation
GO.wall.consistent=goseq(pwf,"hg19","ensGene")
head(GO.samp.consistent)

##Using random sampling
GO.samp.consistent=goseq(pwf,"hg19","ensGene",method="Sampling",repcnt=1000)
head(GO.samp.consistent)

# plot p-vaoues against each other
library(ggplot2)
plot(log10(GO.wall.consistent[,2]), log10(GO.samp.consistent[match(GO.wall.consistent[,1],GO.samp.consistent[,1]),2]),xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)",xlim=c(-3,0))
abline(0,1,col=3,lty=2)

#Making sense of the results
#Having performed the GO analysis, you may now wish to interpret the results. If you wish to identify categories significantly enriched/unenriched below some p-value cutoff, it is necessary to first apply some kind of multiple hypothesis testing correction. For example, GO categories over enriched using a .05 FDR cutoff [Benjamini and Hochberg, 1995] are:
enriched.consistent <- GO.wall.consistent[p.adjust(GO.wall.consistent$over_represented_pvalue,method="BH")<.05, ]
enriched.GO.consistent=GO.wall.consistent$category[p.adjust(GO.wall.consistent$over_represented_pvalue,method="BH")<.05]
head(enriched.consistent)

library(GO.db)

for(go in enriched.GO[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}

write.table(enriched.314, "Overlapping DEG analysis/GEO analysis/GO/314_overlapping_DEGs_GO_terms.txt", sep="/t", row.names = F, quote=F)




###### For overlapping DEGs with P<0.05 & log2FC>2
overlapping.314 <- read.table("~/Documents/EPO Project/CRISPR/Whole gene knock-out/Confirming EPO KO/RNA Sequencing/DeSeq2/Overlapping DEG analysis/Overlapping_genes_KOA_KOB_analysis.txt", sep="\t", header=T)
overlapping.314$DE <- 1
de.genes <- overlapping.314$Geneid

## Read in all genes assayed
gtf <- read.table("~/Downloads/Homo_sapiens.GRCh38.98_gene_annotation_table.txt", sep="\t", header=T)

gtf$Geneid <- gtf$gene_id

## Load in the count data/file with list of gene names
read.counts <- read.table("iCloud Drive (Archive)/Documents/EPO Project/CRISPR/Whole gene knock-out/Confirming EPO KO/RNA Sequencing/FeatureCounts/feature_count_results_geneid.txt", header = TRUE)

#replace all row names with the names of genes
row.names(read.counts) <- read.counts$Geneid
#remove the columns which contain no count data
read.counts <- read.counts[,-c(2:6)]

# give meaningful sample names - this can be achieved via numerous approaches
names(read.counts) <- c("Geneid","WT1", "WT2","WT3","WT4","KOA1", "KOA2", "KOA3", "KOA4", "KOB1", "KOB2", "KOB3", "KOB4")

# Check data is what we expect
str(read.counts)
head (read.counts, n = 3)

## Merge the gtf annotation data frame with the read counts file by geneid to obtain gene_names as well
read.counts.geneid <- merge(read.counts, gtf, by="Geneid")

## make a list of all gene assaysed
all_genes <- read.counts.geneid$gene_id

## create a gene list for deseq2, 1=DE, 0=no DE
gene.vector=as.integer(all_genes%in%de.genes)
names(gene.vector)=all_genes
head(gene.vector)
table(gene.vector)

## see which organisms are supported 
supportedOrganisms()
supportedOrganisms()[supportedOrganisms()$Genome=="hg19",]
# ensGene

#Fitting the Probability Weighting Function (PWF)
pwf=nullp(gene.vector,"hg19","ensGene")
head(pwf)


##Using the Wallenius approximation
GO.wall.314=goseq(pwf,"hg19","ensGene")
head(GO.wall.314)

##Using random sampling
GO.samp.314=goseq(pwf,"hg19","ensGene",method="Sampling",repcnt=1000)
head(GO.samp.314)

# plot p-vaoues against each other
library(ggplot2)
plot(log10(GO.wall[,2]), log10(GO.samp[match(GO.wall[,1],GO.samp[,1]),2]),xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)",xlim=c(-3,0))
abline(0,1,col=3,lty=2)

# By default, goseq tests all three major Gene Ontology branches; Cellular Components, Biological Processes and Molecular Functions. 
# limit what is analysed 
GO.MF.314=goseq(pwf,"hg19","ensGene",test.cats=c("GO:MF"))
GO.BP.314=goseq(pwf,"hg19","ensGene",test.cats=c("GO:BP"))

enriched.GO.314=GO.wall.314$category[p.adjust(GO.wall.314$over_represented_pvalue,method="BH")<.05]
enriched.MF.314=GO.MF.314$category[p.adjust(GO.MF.314$over_represented_pvalue,method="BH")<.05]
enriched.BP.314=GO.BP.314$category[p.adjust(GO.BP.314$over_represented_pvalue,method="BH")<.05]

enriched.314 <- GO.wall.314[p.adjust(GO.wall.314$over_represented_pvalue,method="BH")<.05, ]
enriched.314.MF <- GO.MF.314[p.adjust(GO.MF.314$over_represented_pvalue,method="BH")<.05, ]
enriched.314.BP <- GO.BP.314[p.adjust(GO.BP.314$over_represented_pvalue,method="BH")<.05, ]


#Making sense of the results
#Having performed the GO analysis, you may now wish to interpret the results. If you wish to identify categories significantly enriched/unenriched below some p-value cutoff, it is necessary to first apply some kind of multiple hypothesis testing correction. For example, GO categories over enriched using a .05 FDR cutoff [Benjamini and Hochberg, 1995] are:
enriched.314 <- GO.wall.314[p.adjust(GO.wall.314$over_represented_pvalue,method="BH")<.05, ]
enriched.GO.314=GO.wall.314$category[p.adjust(GO.wall.314$over_represented_pvalue,method="BH")<.05]
head(enriched.GO.314)

library(GO.db)

for(go in enriched.GO[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}

write.table(enriched.314, "Overlapping DEG analysis/GEO analysis/GO/314_overlapping_DEGs_GO_terms.txt", sep="/t", row.names = F, quote=F)



### Enrichr results
enrichr.bp <- read.delim("~/Documents/EPO Project/CRISPR/Whole gene knock-out/Confirming EPO KO/RNA Sequencing/DeSeq2/Overlapping DEG analysis/GEO analysis/Enrichr/Biological_pathways_2021.txt")
enrichr.bp.20 <- enrichr.bp[1:20,]
enrichr.bp.20$Term <- factor(enrichr.bp.20$Term, levels = enrichr.bp.20$Term[order(enrichr.bp.20$P.value, decreasing=TRUE)])

ggplot(enrichr.bp.20, aes(y=Term, x=-log10(P.value))) +
  geom_bar(stat="identity", aes(fill = P.value), show.legend = F) + 
  theme_classic() + xlab("-log10 P-value") + ylab("Biological Pathway 2021") + 
  geom_text(aes(label=Term),size=5, hjust=1, color="white") + theme(
    axis.text.y = element_blank(),
    axis.title.x = element_text(size=18),
    axis.text.x = element_text(size=18),
    axis.title.y = element_text(size=20),
    axis.ticks.y = element_blank())


enrichr.mf <- read.delim("~/Documents/EPO Project/CRISPR/Whole gene knock-out/Confirming EPO KO/RNA Sequencing/DeSeq2/Overlapping DEG analysis/GEO analysis/Enrichr/Molecular_functions_2021.txt")
enrichr.mf.20 <- enrichr.mf[1:20,]
enrichr.mf.20$Term <- factor(enrichr.mf.20$Term, levels = enrichr.mf.20$Term[order(enrichr.mf.20$P.value, decreasing=TRUE)])

ggplot(enrichr.mf.20, aes(y=Term, x=-log10(P.value))) +
  geom_bar(stat="identity", aes(fill = P.value), show.legend = F) + 
  theme_classic() + xlab("-log10 P-value") + ylab("Molecular Function 2021") + 
  geom_text(aes(label=Term),size=5, hjust=1, color="white") + theme(
    axis.text.y = element_blank(),
    axis.title.x = element_text(size=18),
    axis.text.x = element_text(size=18),
    axis.title.y = element_text(size=20),
    axis.ticks.y = element_blank())



enrichr.kegg <- read.delim("~/Downloads/KEGG_2021_Human_table.txt")
enrichr.kegg.20 <- enrichr.kegg[1:10,]
enrichr.kegg.20$Term <- factor(enrichr.kegg.20$Term, levels = enrichr.kegg.20$Term[order(enrichr.kegg.20$P.value, decreasing=TRUE)])

ggplot(enrichr.kegg.20, aes(y=Term, x=-log10(P.value))) +
  geom_bar(stat="identity", aes(fill = P.value), show.legend = F) + 
  theme_classic() + xlab("-log10 P-value") + ylab("KEGG Pathway 2021") + 
  scale_fill_gradient(low = "lightgrey", high = "darkgrey", na.value = NA) +
  geom_text(aes(label=Term),size=7, hjust=0.87, color="black") + 
  theme(
    axis.text.y = element_blank(),
    axis.title.x = element_text(size=22),
    axis.text.x = element_text(size=22),
    axis.title.y = element_text(size=22),
    axis.ticks.y = element_blank())


