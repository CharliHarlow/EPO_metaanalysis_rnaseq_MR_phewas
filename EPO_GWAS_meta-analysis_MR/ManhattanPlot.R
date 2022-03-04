#! /usr/bin/env Rscript

# Charli (Stoneman) Harlow
# Manhattan Plot for meta-analysis of EPO 

args <- commandArgs(TRUE)

results <- args[1]
output  <- args[2]


library(ggplot2)
library(plyr)

setwd("~/Downloads/")
#print(paste("Loading data in ", results, sep=""))
association_data<-read.table("epo_invnt_metaanalysis_results_manhattan_cols_head.txt", header=TRUE)
names(association_data) = tolower(names(association_data))

print("Ensuring chromosome X consistency")
# relabel chr X and XY = 23 if necessary
if ("X" %in% association_data$chr) {
  association_data$chr = revalue(association_data$chr, c("X" = "23"))
}
if ("XY" %in% association_data$chr) {
  association_data$chr = revalue(association_data$chr, c("XY" = "23"))
}

## Check it has worked

sum(association_data$chr == 23)

print("Setting up plot")

sorted_labels <- paste(sort(as.integer(levels(association_data$chr))))
association_data$chr <- factor(association_data$chr, levels = sorted_labels)


# create X-values relative to the last bp of the last chromosome for plot
association_data$x_pos = association_data$pos
last_bp = 0

for (i in 1:23) {
  association_data$x_pos[association_data$chr==i] <- (last_bp + association_data$pos[association_data$chr==i])
  last_bp = last_bp + max(association_data$pos[association_data$chr==i])}

# sort the association data     by position along X-axis
association_data<-association_data[with(association_data, order(x_pos)), ]

sig_snp <- association_data[association_data$chr == "7" & association_data$pos == "100317298", ]

library(devtools)
# generate the plot
#print(paste("Generating manhattan in ", output, sep=""))
tiff(filename="epoadjhb_none_EA_manhattan_label.png", compression="lzw", width = 1200, height = 600, units = "px")
#mh = ggplot(association_data, aes(x=x_pos, y=-log10(p), color=chr)) + geom_point(shape=20, size = 1)
#mh + xlab("") + theme(axis.text.x = element_blank(), legend.position="bottom") +  scale_color_manual(values=c("deeppink2", "orange", "green2", "blue", "yellow", "purple",  "magenta", "darkgreen", "gold", "firebrick", "yellowgreen", "red", "black", "turquoise3", "tomato", "darkblue", "chocolate", "violet", "slategray4", "OrangeRed", "darkblue", "deeppink", "aquamarine3")) + geom_hline(aes(yintercept=-log10(5e-08)))

y <- ggplot(association_data, aes(x=x_pos, y=-log10(p), color=chr)) + 
        geom_point(shape=20, size = 1) + xlab("") + theme(axis.text.x = element_blank(), legend.position="bottom") +  
        scale_color_manual(values=c("deeppink2", "orange", "green2", "blue", "yellow", "purple",  "magenta", "darkgreen", "gold", "firebrick", "yellowgreen", "red", "black", "turquoise3", "tomato", "darkblue", "chocolate", "violet", "slategray4", "OrangeRed", "darkblue", "deeppink", "aquamarine3")) +
        geom_hline(aes(yintercept=-log10(5e-08)))

print(y + annotate("point", x = 1333430944, y = 2.123609, label = "rs1617640", colour="blue"))

dev.off()




#print("Done")
rm(list=ls())

