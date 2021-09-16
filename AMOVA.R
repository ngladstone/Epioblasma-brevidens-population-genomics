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
