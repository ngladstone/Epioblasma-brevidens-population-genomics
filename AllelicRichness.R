library(diveRsity)
library(devtools)

##READ in GENEPOP FILE. I used from Stacks
test_results<-basicStats(infile="populations.r90S.red.snps.gen",
                         fis_ci = TRUE, ar_ci = TRUE, fis_boots = 999, 
                         ar_boots = 999, mc_reps = 9999, 
                         rarefaction = FALSE, ar_alpha = 0.05, 
                         fis_alpha = 0.05)

##Split out Allelic Richness, make into table.  NOTE: Could modify here (and later with variable name) to split out other statistics
allelicRichness <- test_results$ar

##Data check. Note overall values.
#allelicRichness
#allelicRichness["overall",]

##seperate by population and remove loci that were not genotyped for that population (i.e., Allelic Richness of 0)
pop1_allelicR<-allelicRichness[,2]
pop1_allelicR<-as.data.frame(pop1_allelicR)
pop1_allelicR_zeroesRemoved<-pop1_allelicR[apply(pop1_allelicR!=0,1,all),]

pop2_allelicR<-allelicRichness[,3]
pop2_allelicR<-as.data.frame(pop2_allelicR)
pop2_allelicR_zeroesRemoved<-pop2_allelicR[apply(pop2_allelicR!=0,1,all),]

pop3_allelicR<-allelicRichness[,4]
pop3_allelicR<-as.data.frame(pop3_allelicR)
pop3_allelicR_zeroesRemoved<-pop3_allelicR[apply(pop3_allelicR!=0,1,all),]

#Calculate Mean
pop1mean<-mean(pop1_allelicR_zeroesRemoved)
pop2mean<-mean(pop2_allelicR_zeroesRemoved)
pop3mean<-mean(pop3_allelicR_zeroesRemoved)


##Calculate Sstandard Deviation
pop1sd<-sd(pop1_allelicR_zeroesRemoved)
pop2sd<-sd(pop2_allelicR_zeroesRemoved)
pop3sd<-sd(pop3_allelicR_zeroesRemoved)
