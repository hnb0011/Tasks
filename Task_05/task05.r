setwd('C:\\Users\\haile\\Desktop\\Evolution\\Tasks\\Task_05')
install.packages("learnPopGen")
library (learnPopGen)
?coalescent.plot()
coalescent.plot()

		coalescent.plot(n=20,ngen=30,col.order="alternating")
		object<-coalescent.plot()
		print(object)
		plot(object)

Q1 20 These are modified by changing the n value
Q2 4Ne generations
Q3 It looks like about 2 offspring is average. It varies by one--some have one offspring and some have 3.
Q4 The fitness of an allele is shown in this plot, because the colors separate out the number of offspring per allele.
Q5 No
install.packages ("coala")
install.packages ("phytools")
library (coala)
library (phytools)
?coala
?phytools
model <- coal_model(sample_size = 5, loci_number =10, loci_length = 500, ploidy = 2) +
feat_mutation (10) +
feat_recombination (10) +
sumstat_trees() +
sumstat_nucleotide_div()
stats <- simulate (model, nsim = 1)
Diversity <- stats$pi
Nloci <- length (stats$trees)
t1 <- read.tree (text=stats$trees [[1]][1])
plot (t1)
axisPhylo()
Q6: This phylogenetic tree does not examine individuals, but the allele.
Agel <- max(nodeHeights (t1))
t2 <- read.tree(text=stats$trees[[2]][1])
plot(t2)
axisPhylo()
Q7: No
par(mfrow=c(1,2))
plot(t1)
axisPhylo()
plot(t2)
axisPhylo()
compare.chronograms(t1,t2)
t1_1 <- read.tree(text=stats$tree [[1]][1])
t1_2 <- read.tree (text=stats$trees [[1]][2])
compare.chronograms(t1_1, t1_2)
for (locus in 1:Nloci) {
ntrees <- length (stats$trees[[locus]])
for (n in 1:ntrees) {
if (locus == 1&& n --1){
outPhy <- read.tree (text=stats$trees [[locus]][n])
} 
else {
outPhy <- ape:::c.phylo(outPhy, read.tree(text=stats$trees [[locus]][n]))
}
}
}
par (mfrow=c(1,1))
densityTree(outPhy)
t1_1 <- read.tree(text=stats$tree [[1]][2])
t1_2 <- read.tree (text=stats$trees [[1]][2])
compare.chronograms(t1_1, t1_2)
for (locus in 1:Nloci) {
ntrees <- length (stats$trees[[locus]])
for (n in 1:ntrees) {
if (locus == 1&& n --1){
outPhy <- read.tree (text=stats$trees [[locus]][n])
} 
else {
outPhy <- ape:::c.phylo(outPhy, read.tree(text=stats$trees [[locus]][n]))
}
}
}
par (mfrow=c(1,1))
densityTree(outPhy)
model3 <- coal_model(10, 50) +
feat_mutation (par_prior ("theta", sample.int(100,1))) +
sumstat_nucleotide_div()
stats <- simulate (model3, nsim = 40)
mean_pi <- sapply (stats, function (x) mean (x$pi))
theta <- sapply (stats, function (x) x$pars [["theta"]])
plot(mean_pi, theta)
?abline
abline(lsfit(mean_pi, theta))
model4 <- coal_model(sample_size = 25, loci_number =5, loci_length = 250, ploidy = 2) +
feat_mutation (5) +
feat_recombination (5) +
sumstat_trees() +
sumstat_nucleotide_div()
stats <- simulate (model, nsim = 1)
Diversity <- stats$pi
Nloci <- length (stats$trees)
t1 <- read.tree (text=stats$trees [[1]][1])
plot (t1)
axisPhylo()
for (locus in 1:Nloci) {
ntrees <- length (stats$trees[[locus]])
for (n in 1:ntrees) {
if (locus == 1&& n --1){
outPhy <- read.tree (text=stats$trees [[locus]][n])
} 
else {
outPhy <- ape:::c.phylo(outPhy, read.tree(text=stats$trees [[locus]][n]))
}
}
}
par (mfrow=c(1,1))
densityTree(outPhy) 
model5 <- coal_model(sample_size = 15, loci_number =15, loci_length = 750, ploidy = 3) +
feat_mutation (15) +
feat_recombination (15) +
sumstat_trees() +
sumstat_nucleotide_div()
stats <- simulate (model, nsim = 1)
Diversity <- stats$pi
Nloci <- length (stats$trees)
t1 <- read.tree (text=stats$trees [[1]][1])
plot (t1)
axisPhylo()
for (locus in 1:Nloci) {
ntrees <- length (stats$trees[[locus]])
for (n in 1:ntrees) {
if (locus == 1&& n --1){
outPhy <- read.tree (text=stats$trees [[locus]][n])
} 
else {
outPhy <- ape:::c.phylo(outPhy, read.tree(text=stats$trees [[locus]][n]))
}
}
}
par (mfrow=c(1,1))
densityTree(outPhy)
par (mfrow=c(1,1))
densityTree(model5)