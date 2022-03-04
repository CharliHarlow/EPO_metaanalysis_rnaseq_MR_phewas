## Miami Plot for the eQTL data and EPO meta-analysis data ##
## Charli Harlow ##

library(Cairo)
library(ggplot2)
library(plyr)

## For one manhattan plot
# Read in input data
data <- read.table("~/Documents/EPO Project/Colocalisation/eQTL_and_GWAS_500kb_rs1617640.txt", header=T, stringsAsFactors=F)
names(data) = tolower(names(data))

# relabel chr X, Y, XY and MT with their numeric values if necessary
if ("X" %in% data$chr) {
  data$chr = revalue(data$chr, c("X" = "23"))
}
if ("Y" %in% data$chr) {
  data$chr = revalue(data$chr, c("Y" = "24"))
}
if ("XY" %in% data$chr) {
  data$chr = revalue(data$chr, c("XY" = "25"))
}
if ("MT" %in% data$chr) {
  data$chr = revalue(data$chr, c("MT" = "26"))
}

## Force chr and bp as numeric not int
data[12:13] <- lapply(data[12:13], as.numeric)

## make chr a factor
if(is.numeric(data$chr)) data$CHR <- as.factor(data$chr)
## sort the chr into levels
sorted_labels <- paste(sort(as.integer(levels(data$chr))))
data$CHR <- factor(data$chr, levels = sorted_labels)

## x position = base pair, set first bp position and last bp position then order data by this
data$x_pos = data$bp
last_bp = 0
last_bp_vec = rep(0,22)
for (i in 1:22) {
  if(length(data$BP[data$CHR==i])!=0) {
    data$x_pos[data$CHR==i] <- (last_bp + data$BP[data$CHR==i])
    last_bp = last_bp + max(data$BP[data$CHR==i])
    last_bp_vec[i] <- max(data$x_pos[data$CHR==i])
  } else {
    last_bp = last_bp + 2E7
    last_bp_vec[i] <- last_bp
  }
}
data<-data[with(data, order(x_pos)), ]

## Make MIAMI plot
## Add extra geom_point layer - one with -log10(P), other with +log10(P)
## For example - transform GWAS data -log10(P) so will go on top, transform eQTL data +log10(P) so will go on bottom

CairoPNG(filename="miami_plot.png", width = 1200, height = 600, units = "px")

print((ggplot(data, aes(x=x_pos), color=chr)) + 
        geom_point(y=(-log10(data$p_gwas)), colour= "red", shape=20, size = 1) + 
        geom_point(y=(+log10(data$p_eqtl)), colour = "blue", shape=20, size=1) +
        xlab("") + theme(axis.text.x = element_blank(), legend.position="bottom") + 
        scale_color_manual(values=c("deeppink2", "orange", "green2", "blue", "yellow", "purple",  "magenta", "darkgreen", "gold", "firebrick", "yellowgreen", "red", "black", "turquoise3", "tomato", "darkblue", "chocolate", "violet", "slategray4", "OrangeRed", "darkblue", "deeppink", "aquamarine3","steelblue1","mediumorchid","yellowgreen")) +
        geom_hline(aes(yintercept=-log10(5e-08))) + scale_x_continuous(breaks = last_bp_vec,minor_breaks=NULL) + ggtitle("Colocalisation plot"))
dev.off()


## Colocalisation plot for one region (zoomed into 500bp either side of rs1617640)
data <- data[data$BP>99938955 && data$BP<100816347 && data$Chr==7, ]

## Create label for rs1617640
epo_snp <- data[data$rsid=="rs1617640", ]

## Create plot with snp labelled
CairoPNG(filename="~/Documents/EPO Project/Colocalisation/eqtl_gwas_colocalisation_plot_2.png", width = 600, height = 600, units = "px")

CairoSVG(filename="~/Documents/EPO Project/Colocalisation/eqtl_gwas_colocalisation_plot.svg", width = 600, height = 600)
par(mar=c(1,1,1,1), mfrow=c(2,1), cex=1.0, cex.main=0.8, cex.axis=0.8)
ggplot(data, aes(x=BP), colour=chr) + geom_point(aes(x=BP,y=-log10(p_gwas)), shape=20, size=3, colour="magenta") + 
  geom_point(aes(x=BP,y=+log10(p_eqtl)), shape=20, size=3, colour="dodgerblue3") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=90), legend.position="topright", panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
  geom_point(data=epo_snp, aes(x=BP, y=-log10(p_gwas)), colour="black", size=2) +
  geom_point(data=epo_snp, aes(x=BP, y=+log10(p_eqtl)), colour="black", size=2) +
  geom_label(aes(label=rsid, y=-log10(p_gwas)), data=epo_snp, size = 3, nudge_y = 0.18, nudge_x = 60000) +
  geom_label(aes(label=rsid, y=+log10(p_eqtl)), data=epo_snp, size=3, nudge_y = 0.18, nudge_x=60000) +
  labs(color=" ") +
  ylab('log10(P value)                    -log10(P value)')+ xlab('Position on Chr 7') +
  scale_color_manual(name="Study", 
                     labels = c("EPO Meta-analysis", "Liver eQTL"), 
                     values = c("EPO Meta-analysis"="magenta", 
                                "Liver eQTL"="dodgerblue3")) +
  geom_hline(aes(yintercept=-log10(5e-08)), linetype="dashed", colour="darkgrey") + 
  geom_hline(aes(yintercept=+log10(5e-08)), linetype="dashed", colour="darkgrey") +
  geom_hline(aes(yintercept=0), linetype=c(1), colour="darkgrey") +
  scale_x_continuous(breaks = seq(99938955, 100900000, 100000))

par(xpd=TRUE)
plot.new()
legend("right", bty="n", cex=0.8, title="", legend=c("EPO Meta-analysis", "Liver eQTL"), col=c("magenta", "dodgerblue3"), pch=c(16, 16))

dev.off()  

