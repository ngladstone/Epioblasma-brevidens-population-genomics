#assessing populations data for Epioblasma (brevidens only)

#for helpful information on filterig VCF files, follow this link:
#https://grunwaldlab.github.io/Population_Genetics_in_R/qc.html#missing-data

#loading libraries
library(vcfR)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(cowplot)

#importing VCF
vcf <- read.vcfR("populations.snps.vcf.gz",verbose=FALSE)

#extracting a matrix of variant depths from the VCF file
dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)

#quantifying missingess for all samples
myMiss <- apply(dp, MARGIN = 2, function(x){ sum(is.na(x)) })
myMiss <- myMiss/nrow(vcf)

palette(brewer.pal(n=12, name = 'Set3'))

par(mar = c(12,4,4,2))
barplot(myMiss, las = 2, col = 1:12)
title(ylab = "Missingness (%)")

#creating boxplots to assess seq depth across all samples
par(mar=c(12,4,4,2))
boxplot(dp, col=2:8, las=3)
title(ylab = "Depth (DP)")

#creating violin plot to assess seq depth across all samples
par(mar=c(5,4,4,2))

# Melt our matrix into a long form data.frame.
dpf <- melt(dp, varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)
dpf <- dpf[ dpf$Depth > 0,]

# Create a row designator.
# You may want to adjust this
#samps_per_row <- 20
samps_per_row <- 18
myRows <- ceiling(length(levels(dpf$Sample))/samps_per_row)
myList <- vector(mode = "list", length = myRows)

for(i in 1:myRows){
  myIndex <- c(i*samps_per_row - samps_per_row + 1):c(i*samps_per_row)
  myIndex <- myIndex[myIndex <= length(levels(dpf$Sample))]
  myLevels <- levels(dpf$Sample)[myIndex]
  myRegex <- paste(myLevels, collapse = "$|^")
  myRegex <- paste("^", myRegex, "$", sep = "")
  myList[[i]] <- dpf[grep(myRegex, dpf$Sample),]
  myList[[i]]$Sample <- factor(myList[[i]]$Sample)
}

# Create the plot.
myPlots <- vector(mode = "list", length = myRows)
for(i in 1:myRows){
  myPlots[[i]] <- ggplot(myList[[i]], aes(x=Sample, y=Depth)) + 
    geom_violin(fill="#8dd3c7", adjust=1.0, scale = "count", trim=TRUE)
  
  myPlots[[i]] <- myPlots[[i]] + theme_bw()
  myPlots[[i]] <- myPlots[[i]] + theme(axis.title.x = element_blank(), 
                                       axis.text.x = element_text(angle = 60, hjust = 1))
  myPlots[[i]] <- myPlots[[i]] + scale_y_continuous(trans=scales::log2_trans(), 
                                                    breaks=c(1, 10, 100, 800),
                                                    minor_breaks=c(1:10, 2:10*10, 2:8*100))
  myPlots[[i]] <- myPlots[[i]] + theme( panel.grid.major.y=element_line(color = "#A9A9A9", size=0.6) )
  myPlots[[i]] <- myPlots[[i]] + theme( panel.grid.minor.y=element_line(color = "#C0C0C0", size=0.2) )
}

# Plot the plot.
plot_grid(plotlist = myPlots, nrow = myRows)