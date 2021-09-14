library(hierfstat)
library(ade4)
library(ecodist)

##Generate Pairwise Fst
data1<-read.fstat("Epio_r90M_fstat.dat")

##weir-cokerham Fst calculation
pairwise.WCfst(data1,diploid = TRUE)




####MANTEL TEST#####

##Read distances
Distances<-read.csv("GeoDist.csv", header=FALSE)
distances<-as.matrix(Distances)
Distances2<-as.dist(Distances)
Distances2

##Fst values inferred with STACKS
FST_pairwise<-read.csv("STACKS_Fst.csv", header=FALSE)
FST_pairwise2<-as.dist(FST_pairwise)
FST_pairwise2

##Fst_Values inferred above
FST_weirC<-read.csv("weir_cockerham_Fst.csv", header=FALSE)
FST_weirC2<-as.dist(FST_weirC)
FST_weirC2

##Run mantel analysis with 1000 permutations for randomization tests
ibd<-mantel.randtest(FST_pairwise2,Distances2, nrepet = 1000)
plot(ibd)
ibd

#MRM analyses
MRM.results<-MRM(FST_pairwise2 ~ Distances2,nperm=10000)
MRM.results
