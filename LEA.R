library(BiocManager)
library(LEA)

#set working directory to wherever vcf file is first
#load data
input.file <- "populations.r90S.red.snps.vcf"

#convert vcf to geno format for snmf
ebrev <- vcf2geno(input.file, output.file = "ebrev.geno",
         force = TRUE)

#Trying K clusters 1-10, and running 10 iterations per K
obj.snmf = snmf(ebrev, K = 1:10, project = "new",
                repetitions = 1, tolerance = 0.00001, entropy=TRUE, ploidy = 2)

# plot cross-entropy criterion of all runs of the project
#This plot helps us identify the number of clusters to identify
#as the most likely for our population (looks like K=3)
plot(obj.snmf, cex = 1.2, col = "lightblue", pch = 19)

# get the cross-entropy value for each run at K = ?
#I only did one run here for K = 3
ce <- cross.entropy(obj.snmf, K = 2)

#best run is K = 3
best <- which.min(ce)

#qmatrix object contains the matrix of ancestry coefficients for each
#individual and for K = 3 clusters.
qmatrix = Q(obj.snmf, K = 2, run = best)

#generating barplot of ancestry coefficients
barplot(t(qmatrix), col = c("orange", "violet", "lightgreen"),
        border = NA, space = 0, xlab = "individuals",
        ylab = "Admixture coefficients")

#now calling the ancestral genotype freq with K = 3
gmatrix = G(obj.snmf, K = 2, run = best)
head(gmatrix)

#and plotting
barplot(gmatrix,  border = NA, space = 0.2,  col = c("orange","violet","lightgreen"),
        xlab = "Populations", ylab = "Ancestry proportions",  main = "Ancestry matrix") -> gp

show(obj.snmf)
