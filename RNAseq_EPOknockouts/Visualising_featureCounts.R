### Visualising FeatureCounts results #####
### Plotting the results ###
### Charli (Stoneman) Harlow ###


featurecounts <- read.csv("~/Documents/EPO Project/CRISPR/Whole gene knock-out/Confirming EPO KO/RNA Sequencing/FeatureCounts/featureCounts_summary1.csv")
library(ggplot2)
## Plot the counts 
subset <- subset(featurecounts, featurecounts$Type=="Assigned" | featurecounts$Type=="Unassigned_MultiMapping"|featurecounts$Type=="Unassigned_NoFeatures"|featurecounts$Type=="Unassigned_Ambiguity")

png("~/Documents/EPO Project/CRISPR/Whole gene knock-out/Confirming EPO KO/RNA Sequencing/FeatureCounts/feature_counts_summary.png", height=1000, width=1500)
ggplot(data=subset, aes(y=Counts, x=Name1)) +
  geom_bar(aes(fill=Type), stat="identity",position="stack") +
  theme_bw(base_size = 16) +
  ylab("Count") + 
  xlab("Sample") +
  #scale_fill_manual(values=mycolors) +
  theme(legend.position="bottom",
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.4, "cm"),
        legend.title=element_blank(),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))
dev.off()

featurecounts <- read.table("~/Downloads/trimmed_feature_count_results_geneid.txt.summary", sep="\t", head=T)
library(ggplot2)
## Plot the counts 
subset <- subset(featurecounts, featurecounts$Type=="Assigned" | featurecounts$Type=="Unassigned_MultiMapping"|featurecounts$Type=="Unassigned_NoFeatures"|featurecounts$Type=="Unassigned_Ambiguity")

subset <- subset(featurecounts, featurecounts$Type=="Assigned" | featurecounts$Type=="Unassigned_MultiMapping"|featurecounts$Type=="Unassigned_NoFeatures"|featurecounts$Type=="Unassigned_Ambiguity")


png("~/Documents/EPO Project/CRISPR/Whole gene knock-out/Confirming EPO KO/RNA Sequencing/FeatureCounts/feature_counts_summary.png", height=1000, width=1500)
ggplot(data=subset, aes(y=Counts, x=Name1)) +
  geom_bar(aes(fill=Type), stat="identity",position="stack") +
  theme_bw(base_size = 16) +
  ylab("Count") + 
  xlab("Sample") +
  #scale_fill_manual(values=mycolors) +
  theme(legend.position="bottom",
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.4, "cm"),
        legend.title=element_blank(),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))
dev.off()



