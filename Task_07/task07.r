setwd('C:\\Users\\haile\\Desktop\\Evolution\\Tasks\\Task_07')
install.packages("phytools")
install.packages("ape")
library("phytools")
library("ape")
text.string <-
	"(((((((cow, pig), whale),(bat,(lemur, human))), (robin, iguana)), coelacanth
	), (gold_fish, trout)), shark);"
vert.tree <- read.tree(text=text.string)
plot(vert.tree, edge.width=2)
nodelabels(frame='circle', bg='white', cex=1)
Q1: The shark
vert.tree
Q2: No. It is "rooted"
str(vert.tree)
tree <- read.tree(text='(((A,B), (C,D)), E);')
plotTree(tree, offset=1)
tiplabels(frame='circle', bg='lightblue', cex=1)
nodelabels(frame='circle', bg='white', cex=1)
tree$tip.label
tree$edge
AnolisTree <- force.ultrametric(read.tree("https://jonsmitchell.com/data/anolis.tre"))
par(las=1)
hist(AnolisTree$edge.length, col='black', border='white', main='', xlab='edge lengths for Anolis tree', ylim=c(0,50), xlim=c(0,6))
tipEdges <- which(AnolisTree$edge[,2] <= Ntip(AnolisTree))
tipEdges
Lengths <- AnolisTree$edge.length
names(Lengths) <- AnolisTree$tip.label
names(Lengths)[which(Lengths == min(Lengths))]
plot(AnolisTree, cex=0.25)
Labs <- sapply(AnolisTree$edge.length, round, digits=2)
edgelabels(text=Labs, cex=0.25)
?plot.phylo
Q3:
tree <- read.tree(text='(((A,B), (C,D)), E);')
plot.phylo(tree, type='phylogram', show.tip.label=FALSE, edge.color='green')

Q4
plot.phylo(tree, type='radial')

Q5
plot.phylo(tree, tip.color = 'red')

Q6-8:
plot(AnolisTree, cex=0.25)
Labs <- sapply(AnolisTree$edge.length, round, digits=2)
edgelabels(text=Labs, cex=0.25)
which(Lengths == min(Lengths))
Q: Anolis_occultus 
names(Lengths)
AnolisTree2 <- drop.tip(AnolisTree, 'Anolis_occultus')
plot(AnolisTree2, cex=0.25)

ltt(AnolisTree)
abline(0, 1, lwd=2, col='red', lty=2)
Q: The red dotted line is a positive slope. The other line shows a curve that seems to begin flattening out by the end of the graph. This tells me that the lizards' population is also leveling out.

Q10:
?fit.bd()
fit.bd(AnolisTree, rho = 0.2)

