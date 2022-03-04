### PheWAS plot of rs1617640 ###
## Charli Harlow ##

## Read in phewas results
phewas <- read.table("~/Documents/EPO_Consortium/PheWAS of EPO SNP/rs1617640_PheWAS_290919.txt", header = T, sep = "\t")

## Remove the . at the end of each Trait name
phewas$Trait <- as.character(phewas$Trait)
phewas$Trait = substr(phewas$Trait,1,nchar(phewas$Trait)-1)

## Read in refined phewas list
refined_list <- read.table("~/Documents/Scripts/R/phewasTraits_Rachel_refined_long_list.txt",  header = T, sep = "\t")

## Merge the refined phewas list and the results file based on Trait name
results <- merge(phewas, refined_list, by.x="Trait", by.y="Trait")

## Make new file with all results p<0.05

filtered_results <- subset(results, P_BOLT_LMM<0.05)
write.table(filtered_results, file="~/Documents/EPO_Consortium/PheWAS of EPO SNP/rs1617640_phewas_0.05_290919.txt", sep="\t", row.names = F)

## Plot the results
################### phewas results divided by category ###########################
#for this plot I made a new excel file with only p<0.05 results
#and then I manually added a category for each trait
library(ggplot2)

# read in file
phewas <- read.table("~/Documents/EPO Project/PheWAS of EPO SNP/rs1617640_phewas_0.05_290919.txt", sep="\t", header=T)
plot <- read.table("~/Documents/EPO_Consortium/PheWAS of EPO SNP/rs1617640_phewas_0.05_290919.txt", sep="\t", header=T)

# subet the p-value thresholds

p_value <- data.frame(label=c("Bonferroni Corrected P-value (P<5.75e-05)","Genome-wide significance P-value (P<5e-08)"), value=c(-log10(5.75e-05), -log10(5e-08)))

# create plot of points and categories
main_plot <- ggplot(data=phewas, aes(x=Category, y=-log10(P_BOLT_LMM), col=Gender)) +
  geom_jitter(width=0.3, size=3) +
  #geom_hline(data=p_value, aes(yintercept=value), linetype=c("dotted", "dashed")) + #both lines
  #geom_hline(aes(yintercept=-log10(1.5e-05)), linetype="dotted", colour="grey3")+  #bonferroni correction
  #geom_hline(aes(yintercept=-log10(5e-08)), linetype="dashed", colour="grey3")+ #genome significance
  scale_colour_manual(values = c("green1", "dodgerblue3", "deeppink2"), labels=c("Both", "Female Only", "Male Only"))+
  theme_classic()+
  ylab("-log10(p)")+xlab("Category") +
  #geom_label(aes(label=P_BOLT_LMM<5e-08),hjust=1,vjust=0, size=3)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = unit(c(1,3,1,1), "lines")) 

# add on the thresholds and label them
## create threshold line label names
labels <- p_value$label

labelled_threshold <- main_plot + geom_hline(data=p_value, aes(yintercept=value, linetype=factor(label))) + 
  scale_linetype_manual(name = "P-value threshold", values = c("dotted","dashed"), labels = labels)

labelled_threshold
ggsave("~/Documents/EPO Project/PheWAS of EPO SNP/PheWAS_categories.png", plot = labelled_threshold, width = 30, height = 20, units = "cm")

