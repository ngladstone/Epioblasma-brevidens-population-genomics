library(adegenet)
library(ape)
library(poppr)
library(na.tools)

data <- read.genepop(file="populations.r90S.red.snps.gen")

grp <- find.clusters(data, max.n.clust=20)

#cross validation to identify number of PCs to retain for DAPC

temp <- optim.a.score(dapc)

#dapc

dapc <- dapc(data,grp$grp)
#1 PCs retained

pdf("Ebrev_DAPC.pdf")
myCol <- c("green","red","blue")
rivCol <- c("blue", "red")
scatter(dapc, ratio.pca=0.3,solid=0.4,cex=3,clab=0,cell=0,posi.da="bottomright",
        bg="white",pch=17:22, cstar=0, col=myCol, scree.pca=TRUE, posi.pca="topleft")

dev.off()

#also safe as png for later edits
#Cluster 1 = Clinch (Green)
#Cluster 2 = Bear (Red)
#Cluster 3 = Troublesome (Blue)

#ascore