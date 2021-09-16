library(adegenet)
library(poppr)
library(ape)


ebrev.geneid.nout <- read.genepop("populations.r90S.red.snps.gen")
ebrev_strata <- read.csv(file = "ebrev_strata.csv",header = TRUE)
ebrev_strata
strata(ebrev.geneid.nout)=ebrev_strata
ebrev.geneid.nout



##AMOVA
ebrev.Region.amova = poppr.amova(ebrev.geneid.nout, ~riv/pop, cutoff = 0.9)
ebrev.pop.amova=poppr.amova(ebrev.geneid.nout,~pop,cutoff=0.9)

####Print Results
ebrev.Region.amova
ebrev.pop.amova
##Randomization Test

ebrev.Region.rtest<-randtest(ebrev.Region.amova,nrepet = 999)
ebrev.Region.rtest

ebrev.site.amova.rtest<-randtest(ebrev.pop.amova,nrepet=999)
ebrev.site.amova.rtest
plot(ebrev.site.amova.rtest)


##Flinger tests of homogeneity
##Fligner tests apparently handle non-normal datasets well.


#Seperate out collection sites and remove loci with only missing data within collection site.
seperated<-seppop(ebrev.geneid.nout)
ebrev.jordan.noFullyMissing<-missingno(seperated$Jordan9_2,cutoff = 0.99)
ebrev.loving.noFullyMissing<-missingno(seperated$LovingH9,cutoff = 0.99)
ebrev.black.noFullyMissing<-missingno(seperated$Black9_2,cutoff = 0.99)
ebrev.brown.noFullyMissing<-missingno(seperated$Brown9_2,cutoff = 0.99)

#Calculate population summary statistics
div<-summary(ebrev.geneid.nout)
jordan_div<-summary(ebrev.jordan.noFullyMissing)
loving_div<-summary(ebrev.loving.noFullyMissing)
black_div<-summary(ebrev.black.noFullyMissing)
brown_div<-summary(ebrev.brown.noFullyMissing)

##Perform fligner test
fligner.test(list(div$Hexp, div$Hobs))
fligner.test(list(jordan_div$Hexp, jordan_div$Hobs))
fligner.test(list(loving_div$Hexp, loving_div$Hobs))
fligner.test(list(black_div$Hexp, black_div$Hobs))
fligner.test(list(brown_div$Hexp, brown_div$Hobs))

mean(brown_div$Hexp)
mean(brown_div$Hobs)