## ***** Visualising STAR alignments *****
## *****  Charli (Stoneman) Harlow  *****
## *****  16th January 2020 *****

## Set working directory
setwd("/gpfs/mrc0/projects/Research_Project-MRC158833/cs660/EPO_project/RNA_sequencing_analysis/Star_Alignment_default/test_2")

## Setting up the functions needed
# ################################################################################
# The functions of this script are not very generalized.
# Most of them expect a long data frame with exactly the following format:
# head(align.results.df)
#                                                                                           V1          V2 sample replicate
#/Users/cs660/Downloads//Empty1_10bp_Log.final.out.1 Mapping speed, Million of reads per hour      264.01  Empty         1
#/Users/cs660/Downloads//Empty1_10bp_Log.final.out.2                    Number of input reads 13493936.00  Empty         1
#/Users/cs660/Downloads//Empty1_10bp_Log.final.out.3                Average input read length       89.00  Empty         1
#/Users/cs660/Downloads//Empty1_10bp_Log.final.out.5             Uniquely mapped reads number 11782768.00  Empty         1
#/Users/cs660/Downloads//Empty1_10bp_Log.final.out.6                  Uniquely mapped reads %       87.32  Empty         1
#/Users/cs660/Downloads//Empty1_10bp_Log.final.out.7                    Average mapped length       89.56  Empty         1
##################################################################################


PlottingCorrelation <- function(DF, Var1, Var2, Var1.Label, Var2.Label){
  # Convenience function for the simple plot() function that allows for separate
  # definition of labels and columns that should be compared against each other
  # usage: PlottingCorrelation(DF=aligned.reads.df,
  #               Var1="Number of input reads", Var2="Uniquely mapped reads %",
  #               Var1.Label = "InputReads", Var2.Label="UniquelyMappedFraction")
  m <- matrix(data = c(DF$V2[which(DF$V1 == Var1)],
                       DF$V2[which(DF$V1 == Var2)]),
              ncol = 2)
  colnames(m) <- c(Var1.Label,Var2.Label)
  plot(m)
}

# Setting up colour for plot
library(RColorBrewer)
nb.cols <- 5 # change this to number of coloumns
mycolors <- colorRampPalette(brewer.pal(5, "YlOrRd"))(nb.cols)

# Setting up the function to produce the plot of 4 side by side
PlottingAlignmentResults <- function(Filter, DF, Legend=TRUE, PlotMedian = TRUE){
  # this function extracts those lines that correspond to the value stored in Filter and generates a bar plot where each sample is shown with a different color
  
  library(ggplot2)
  library(grid) # for unit() function
  filtered.df <- DF[which(DF$V1 == Filter),]
  medians <- as.data.frame(aggregate(V2~sample, data=filtered.df, FUN=median))
  filtered.df <- merge(filtered.df, medians, by.x = "sample", by.y = "sample", all.x=TRUE)
  
  p <- ggplot(data=filtered.df, aes(fill=replicate, y=V2.x, x=sample)) +
    geom_bar(stat="identity",position=position_dodge()) +
    theme_bw(base_size = 8) +
    scale_fill_manual(values=mycolors) +
    theme(legend.position="bottom",
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.5, "cm"),
          legend.title=element_blank()) +
    coord_flip() + ylab("") + ggtitle(Filter)
  
  if(PlotMedian){
    p <- p + geom_errorbar(aes(y=V2.y, ymax=V2.y, ymin=V2.y), linetype="dashed")
  }
  
  if(!Legend){
    p <- p +  theme(legend.position="none")
  }
  
  return(p)
}

# Set up the Legend function
ExtractLegend <- function(Plot){
  library(ggplot2)
  G <- ggplotGrob(Plot)$grobs
  Legend <- G[[which(sapply(G, function(x) x$name) == "guide-box")]]
  Lheight <- sum(Legend$height)
  return(list(legend = Legend, lheight=Lheight))
}

Extract.Histo.Info <- function(InList, # output from hist(...,plot=FALSE)
                               Percentage = TRUE
){
  # this function uses the information from hist() that are stored in lists
  # to make a data frame that's suitable for bar plots
  ll <- length(InList$breaks)
  out.df <- data.frame(breaks1 = InList$breaks[c(1:ll-1)],
                       breaks2 = InList$breaks[c(2:ll)],
                       counts = InList$counts)
  if(Percentage){
    out.df <- transform(out.df,
                        breaks = paste(out.df$breaks1*100, "-", out.df$breaks2*100, sep = ""),
                        breaks1 = NULL, breaks2 = NULL)
  }else{
    out.df <- transform(out.df,
                        breaks = paste(out.df$breaks1, "-", out.df$breaks2, sep = ""),
                        breaks1 = NULL, breaks2 = NULL)
  }
  out.df$breaks <- factor(out.df$breaks, levels = unique(as.character(out.df$breaks)), ordered = TRUE)
  return(out.df)
}


# Load in libraries
# library(Cairo)
library(ggplot2)
library(gridExtra) # for composite plotting of ggplots

# Read in the final.log.out files
infiles <- list.files(pattern="Log.final.out", full.names = TRUE) # listing the files to be read in
head(infiles)
list(infiles)
# Should look like:
## [[1]]
## [1] "./Empty1Log.final.out" "./Empty2Log.final.out" "./Empty3Log.final.out" "./Empty4Log.final.out"
## [5] "./KOA1Log.final.out"   "./KOA2Log.final.out"   "./KOA3Log.final.out"   "./KOA5Log.final.out"  
## [9] "./KOB1Log.final.out"   "./KOB211Log.final.out" "./KOB3Log.final.out"   "./KOB5Log.final.out" 

align.results <- lapply(infiles, function(x) read.table(x, sep="|", strip.white=TRUE, stringsAsFactor=FALSE, skip=3, fill = TRUE, header = FALSE)) #iterating over the file list to generate a list of data frames
typeof(align.results) #check its a list
## [1] "list"

head(align.results[[1]])
## V1       V2
## 1 Mapping speed, Million of reads per hour    81.46
## 2                    Number of input reads 13463454
## 3                Average input read length      149
## 4                            UNIQUE READS:         
##  5             Uniquely mapped reads number 12551309
## 6                  Uniquely mapped reads %   93.23%
  
align.results <- lapply(align.results, function(x)transform(x, V2 = as.numeric(gsub("%", "", x$V2) ))) #remove the % from some of numbers so just a number
names(align.results) <- gsub("(Empty|KOA|KOB)*(\\_[0-12]*)*", "\\1\\2", infiles) # some cosmetics of each data frames name, specific for the sample names of the files used here. Instead of 0-12, it will name them by the names in the infiles)
names(align.results)
align.results.df <- as.data.frame(do.call(rbind, align.results)) # catenating all data frames of align.results together
align.results.df <- align.results.df[complete.cases(align.results.df),] # remove lines without any values
head(align.results.df) 


# adding an additional column of sample names
# Sample name looks like this before:./Empty1Log.final.out.1	
# Want it to look like: Empty, KOA, KOB

align.results.df$sample <- gsub("(//*)\\_.*.", "\\1", row.names(align.results.df)) #create a new column with sample name that is the same as the row names
align.results.df$sample <- sub("/*./", "", align.results.df$sample) #remove the ./ from each sample name
align.results.df$sample <- sub("_alignment.*", "", align.results.df$sample) #remove the Log.final.out.# from each sample name

align.results.df$sample <- gsub('[[:digit:]]+', '', align.results.df$sample) # remove the digits from the end of the sample names

head(align.results.df)

# V1          V2 sample
# ./Empty1Log.final.out.1 Mapping speed, Million of reads per hour       81.46  Empty
# ./Empty1Log.final.out.2                    Number of input reads 13463454.00  Empty
# ./Empty1Log.final.out.3                Average input read length      149.00  Empty
# ./Empty1Log.final.out.5             Uniquely mapped reads number 12551309.00  Empty
# ./Empty1Log.final.out.6                  Uniquely mapped reads %       93.23  Empty
# ./Empty1Log.final.out.7                    Average mapped length      149.45  Empty


# adding an additional column of replicate id
# we want each replicate to be number 1-4
# So Empty1 = replicate 1, Empty2 = replicate 2, Empty 3 = replicate 3....
# KOA1 = replicate 1, KOA2 = replicate 2, KOA3 = replicate 3, KOA5 = replicate 4
# KOB1 = replicate 1, KOB211 = replicate 2, KOB3 = replicate 3, KOB5 = replicate 4


align.results.df$replicate <- gsub("(//*)\\_.*.", "\\1", row.names(align.results.df)) #create a new column with sample name
align.results.df$replicate <- sub("/*./", "", align.results.df$replicate) #remove the ./ from each sample name
align.results.df$replicate <- sub("_.*", "", align.results.df$replicate) #remove the bits after the _ from each sample name
align.results.df$replicate <- gsub("[^0-9.-]+", "", align.results.df$replicate) #remove chracters apart from the number after the sample name e.g. Empty1 = 1
align.results.df$replicate <- ifelse(align.results.df$replicate == 5, 4, 
                                     ifelse(align.results.df$replicate == 211, 2, 
                                            align.results.df$replicate)) #replace some of the funny sample names to be replicates 1,2,3&4
align.results.df$replicate <- as.factor(as.numeric(align.results.df$replicate)) #change column to a factor variable


# find the unique fields
unique(align.results.df$V1)


# define those entries that we are interested in:

filters = c("Number of input reads", "Uniquely mapped reads %", "% of reads mapped to multiple loci",
            "% of reads mapped to too many loci", "% of reads unmapped: too short", "% of reads unmapped: other")

filters_N = c("Number of input reads", "Uniquely mapped reads number",
            "Number of splices: Total", "Number of splices: Non-canonical", "Number of reads unmapped: other")

# Plot the graphs for these filters
# CairoPNG(filename="./Alignment_mapping_percentages_bar_chart.png", height=600, width=1200)

plots <- lapply(filters, function(x) 
  PlottingAlignmentResults(x, align.results.df, Legend = FALSE))

# Add legend 
my.legend <- ExtractLegend(PlottingAlignmentResults(align.results.df, 
                                                    Filter = filters[1], 
                                                    Legend=TRUE))
# combining plots and legend
grid.arrange(arrangeGrob(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]],plots[[6]], nrow=2),
             my.legend$legend, nrow=2,
             heights= unit.c(unit(1, "npc") - my.legend$lheight, my.legend$lheight)
)

## Plotting the numbers

plots <- lapply(filters_N, function(x) 
  PlottingAlignmentResults(x, align.results.df, Legend = FALSE))

# Add legend 
my.legend <- ExtractLegend(PlottingAlignmentResults(align.results.df, 
                                                    Filter = filters_N[1], 
                                                    Legend=TRUE))
# combining plots and legend
grid.arrange(arrangeGrob(plots[[1]], plots[[2]], plots[[3]], plots[[4]],plots[[5]], nrow=2),
             my.legend$legend, nrow=2,
             heights= unit.c(unit(1, "npc") - my.legend$lheight, my.legend$lheight)
)

## Stacked charts
filtered.df <- subset(align.results.df, align.results.df$V1 == "% of reads mapped to multiple loci" | align.results.df$V1 == "% of reads unmapped: too short" | align.results.df$V1 == "% of reads unmapped: other" | align.results.df$V1 ==  "Uniquely mapped reads %")

filtered.df$samplename <- gsub("(//*)\\_.*.", "\\1", row.names(filtered.df)) #create a new column with sample name
filtered.df$samplename <- sub("./*/", "", filtered.df$samplename) #remove the ./ from each sample name
filtered.df$samplename <- sub("_.*", "", filtered.df$samplename) #remove the bits after the _ from each sample name

filtered.df$V1 <- as.factor(filtered.df$V1)
ggplot(data=filtered.df, aes(y=V2, x=samplename)) +
  geom_bar(aes(fill=factor(filtered.df$V1, levels=c("Uniquely mapped reads %","% of reads unmapped: too short", "% of reads unmapped: other", "% of reads mapped to multiple loci"))), stat="identity",position="stack") +
  theme_bw(base_size = 12) +
  scale_fill_manual(values = c("maroon2","yellow", "black", "dodgerblue1"))+
  ylab("Percentage of reads (%)") + 
  xlab("Sample") +
  #scale_fill_manual(values=mycolors) +
  theme(legend.position="bottom",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"),
        legend.title=element_blank())

# Create a table with results we want
list(unique(align.results.df$V1))

table.df <- subset(align.results.df, align.results.df$V1 == "% of reads mapped to multiple loci" | align.results.df$V1 == "% of reads unmapped: too short" | align.results.df$V1 == "% of reads unmapped: other" | align.results.df$V1 ==  "Uniquely mapped reads %" | align.results.df$V1 == "Number of input reads" | align.results.df$V1 ==  "Number of reads mapped to multiple loci" | align.results.df$V1 == "Uniquely mapped reads number"  | align.results.df$V1 == "Number of reads mapped to too many loci" | align.results.df$V1 == "Number of reads unmapped: too short" )

table.df$samplename <- gsub("(//*)\\_.*.", "\\1", row.names(table.df)) #create a new column with sample name
table.df$samplename <- sub("./*/", "", table.df$samplename) #remove the ./ from each sample name
table.df$samplename <- sub("_.*.", "", table.df$samplename) #remove the bits after the _ from each sample name

write.table(table.df, file="./alignment_summary_table_test.txt", row.name=F, quote=F, sep="\t")

