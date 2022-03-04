## Forest Plot of MR Results ##
## Charli Harlow ##

library(meta)

cont <- read.table("~/Documents/EPO Project/MR rs1617640 analysis/MRBase/Continuous_MR.txt", sep='/t')

cont$Outcome <- as.character(cont$Outcome)
cont$E95 <- as.character(cont$E95)


tabletext <- cbind(c("Risk Factor",cont$Outcome), 
                   c("Effect Estimate (95% CI)",cont$E95),
                   c("P-Value",cont$P.value))

library(forestplot)

png("~/Documents/EPO Project/MR rs1617640 analysis/MRBase/Forestplot_continuousTraits_Waldratio_2.png", width=1000, height=800)
forestplot(labeltext=tabletext, graph.pos=2, align=T,
           mean=c(NA,cont$Effect.estimate), 
           lower=c(NA,cont$Effect.estimate-(cont$SE*1.96)), upper=c(NA,cont$Effect.estimate+(cont$SE*1.96)),
           #title="Hazard Ratio",
           xlab="Decreased levels                                                                                 Increased levels \n Effect estimate for risk factor per 1SD increase in EPO levels",
           hrzl_lines=list("2" = gpar(lwd=1, col="grey50")), 
           #"3" = gpar(lwd=60, lineend="butt", columns=c(2:4), col="#99999922")),
           #"15" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922"),
           #"23" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922"),
           #"31" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922")),
           txt_gp=fpTxtGp(label=gpar(cex=1.0),
                          ticks=gpar(cex=1.0),
                          xlab=gpar(cex = 1.0),
                          title=gpar(cex = 1.2)),
           xticks=c(-2.0,-1.6,-1.2,-0.8,-0.4,0.0,0.4, 0.8,1.2,1.6, 2.0),
           col=fpColors(box="black", lines="black", zero = "grey50"),
           zero=1, cex=1.0, boxsize=0.08, colgap=unit(5,"mm"),
           lwd.zero=1, lwd.ci=1, ci.vertices=TRUE, ci.vertices.height = 0.1)
dev.off()

### Combining meta-analysis forrest plots
# set up the correct text for the forest plot
cat <- read.table("~/Documents/EPO Project/MR rs1617640 analysis/EPO SNPs - Outcome/Meta-analyses/Meta_analysis_results_categorical_for_combined_forest.txt", header=T, sep='\t')

cat$Outcome <- as.character(cat$Outcome)
cat$E95 <- as.character(cat$E95)
cat$Study <- as.character(cat$Study)

## Labels defining subgroups are a little indented!
subgps <- c(2,3,4,6,7,8,10,11,12)
cat$Study[subgps] <- paste("  ",cat$Study[subgps])


tabletext <- cbind(c("Study",cat$Study), 
                   c("Odds Ratio (95% CI)",cat$E95),
                   c("P-Value",cat$P.value))


png("~/Documents/EPO Project/MR rs1617640 analysis/EPO SNPs - Outcome/Meta-analyses/Combined_forest_plot_Categorical.png", width=1000, height=800)
forestplot(labeltext=tabletext, graph.pos=2, align=T,
           mean=c(NA,cat$OR), 
           lower=c(NA,cat$L95), upper=c(NA,cat$U95),
           #title="Hazard Ratio",
           xlab="<-- EPO-raising allele decreases risk of disease --                 -- EPO-raising allele increases risk of disease -->",
           hrzl_lines=list("2" = gpar(lwd=1, col="grey50"), 
           #"4" = gpar(lwd=40, lineend="butt", columns=c(1:4), col="#99999922"),
           "4" = gpar(lwd=300, lineend="butt", columns=c(1:4), col="#99999922"),
           #"3" = gpar(lwd=1, lty=2, col="grey50"),
           "5" = gpar(lwd=1, lty=2, col="grey50"),
           "6" = gpar(lwd=1, lty=2, col="grey50"),
           "9" = gpar(lwd=1, lty=2, col="grey50"),
           "10" = gpar(lwd=1, lty=2, col="grey50"),
           #"12" = gpar(lwd=1, lty=2, col="grey50")),
           "12" = gpar(lwd=300, lineend="butt", columns=c(1:4), col="#99999922"),
           "13" = gpar(lwd=1, lty=2, col="grey50")),
           #"14" = gpar(lwd=1, lty=2, col="grey50")),
           fn.ci_sum=function(col, size, ...) {
           fpDrawDiamondCI(clr.line = "blue", clr.marker = "blue", size=0.18, vertices.height=0.1, lwd=1,...)
           },
           
           #"15" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922"),
           #"23" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922"),
           #"31" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922")),
           is.summary=c(FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE),
           txt_gp = fpTxtGp(summary = list(
             gpar(fontface="bold"), 
             gpar(fontface="bold")),
             label=gpar(cex=1.0),
             ticks=gpar(cex=1.0),
             xlab=gpar(cex = 0.8),
             title=gpar(cex = 1.2)
           
           ),
           xticks=c(0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4),
           col=fpColors(box=c(rep("black", times=2),"red"), lines="black", zero = "grey50"),
           fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI,fpDrawNormalCI,fpDrawNormalCI, fpDrawDiamondCI,fpDrawNormalCI, fpDrawNormalCI,fpDrawNormalCI,fpDrawDiamondCI, fpDrawNormalCI, fpDrawNormalCI,fpDrawNormalCI, fpDrawDiamondCI),
           zero=1, cex=1.0, boxsize=c(0.08,0.08, 0.08, 0.08, 0.15, 0.08,0.08, 0.08, 0.15,0.08,0.08, 0.08, 0.15), colgap=unit(5,"mm"),
           lwd.zero=1, lty.zero=3, lwd.ci=1, ci.vertices=TRUE, ci.vertices.height = 0.1)

dev.off()

### Combining meta-analysis forrest plots
# set up the correct text for the forest plot
cont <- read.table("~/Documents/EPO Project/MR rs1617640 analysis/EPO SNPs - Outcome/Meta-analyses/Meta_analysis_results_continuous_for_combined_forest.txt", header=T, sep='\t')

cont$Outcome <- as.character(cont$Outcome)
cont$E95 <- as.character(cont$E95)
cont$Study <- as.character(cont$Study)

## Labels defining subgroups are a little indented!
subgps <- c(2,3,4,6,7,8,10,11,12)
cont$Study[subgps] <- paste("  ",cont$Study[subgps])


tabletext <- cbind(c("Study",cont$Study), 
                   c("Effect estimate (95% CI)",cont$E95),
                   c("P-Value",cont$P.value))


png("~/Documents/EPO Project/MR rs1617640 analysis/EPO SNPs - Outcome/Meta-analyses/Combined_forest_plot_Continuous.png", width=1000, height=800)
forestplot(labeltext=tabletext, graph.pos=2, align=T,
           mean=c(NA,cont$OR), 
           lower=c(NA,cont$L95), upper=c(NA,cont$U95),
           #title="Hazard Ratio",
           xlab="<-- EPO-raising allele decreases levels --                 -- EPO-raising allele increases level  -->",
           hrzl_lines=list("2" = gpar(lwd=1, col="grey50"), 
                           #"4" = gpar(lwd=40, lineend="butt", columns=c(1:4), col="#99999922"),
                           "4" = gpar(lwd=300, lineend="butt", columns=c(1:4), col="#99999922"),
                           #"3" = gpar(lwd=1, lty=2, col="grey50"),
                           "5" = gpar(lwd=1, lty=2, col="grey50"),
                           "6" = gpar(lwd=1, lty=2, col="grey50"),
                           "9" = gpar(lwd=1, lty=2, col="grey50"),
                           "10" = gpar(lwd=1, lty=2, col="grey50"),
                           #"12" = gpar(lwd=1, lty=2, col="grey50")),
                           "12" = gpar(lwd=300, lineend="butt", columns=c(1:4), col="#99999922"),
                           "13" = gpar(lwd=1, lty=2, col="grey50")),
           #"14" = gpar(lwd=1, lty=2, col="grey50")),
           fn.ci_sum=function(col, size, ...) {
             fpDrawDiamondCI(clr.line = "blue", clr.marker = "blue", size=0.18, vertices=T,vertices.height=0.1, lwd=1,...)
           },
           is.summary=c(FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE),
           txt_gp = fpTxtGp(summary = list(
             gpar(fontface="bold"), 
             gpar(fontface="bold")),
             label=gpar(cex=1.0),
             ticks=gpar(cex=1.0),
             xlab=gpar(cex = 0.8),
             title=gpar(cex = 1.2)
             
           ),
           xticks=c(-0.16, -0.12, -0.08, -0.04, 0, 0.04, 0.08, 0.12, 0.16),
           col=fpColors(box=c(rep("black", times=2),"red"), lines="black", zero = "grey50"),
           fn.ci_norm = c(fpDrawNormalCI, fpDrawNormalCI,fpDrawNormalCI,fpDrawNormalCI, fpDrawDiamondCI,fpDrawNormalCI, fpDrawNormalCI,fpDrawNormalCI,fpDrawDiamondCI, fpDrawNormalCI, fpDrawNormalCI,fpDrawNormalCI, fpDrawDiamondCI),
           zero=0, cex=1.0, boxsize=c(0.08,0.08, 0.08, 0.08, 0.15, 0.08,0.08, 0.08, 0.15,0.08,0.08, 0.08, 0.15), colgap=unit(6,"mm"),
           lwd.zero=1, lty.zero=3, lwd.ci=1, ci.vertices=TRUE, ci.vertices.height = 0.1)

dev.off()

