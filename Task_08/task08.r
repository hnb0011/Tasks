setwd('C:\\Users\\haile\\Desktop\\Evolution\\Tasks\\Task_08')
library("phytools")
library("ape")
tree <- read.tree("https://jonsmitchell.com/data/anolis.tre")
plot(tree, type="fan")
tree$tip.label
Q1: There are 82 branch tips and yes, the lengths are present.
Data <- read.csv ("https://jonsmitchell.com/data/svl.csv", stringsAsFactors=F, row.names=1)
head(Data)
dim(Data)
Q2: This object contains various Anolis species with their snout-vent lengths. The dimension is 82 x 1
sv1<- setNames(Data$svl, rownames(Data))
sv1
Ancestors <- fastAnc(tree, sv1,vars=TRUE,CI=TRUE)
Ancestors
?fastAnc
Q3: CI95 is a 95% confidence interval, and the estimated values given are a result of re-rooting and analyzing the tree. So the data is stored in the tree.
Q4: One is the maximum liklihood of states arising in the model and the other is assuming the model is 95% correct with the addition of the confidence interal. There is still a 5% chance of being wrong, though.
par (mar=c(0.1, 0.1, 0.1, 0.1))
plot(tree, type= "fan", lwd=2, show.tip.label=F)
tiplabels (pch=16, cex= 0.25*sv1[tree$tip.label])
nodelabels (pch=16, cex=0.25*Ancestors$ace)
obj <- contMap(tree, sv1, plot=F)
plot (obj, type="fan", legend=0.7*max(nodeHeights(tree)), sig=2, fsize=c(0.7, 0.9))
fossilData <- data.frame(sv1=log(c(25.4, 23.2, 17.7, 19.7, 24, 31)), tip1=c("Anolis_aliniger", "Anolis_aliniger", "Anolis_occultus", "Anolis_ricordii", "Anolis_cristatellus", "Anolis_occultus"), tip2=c("Anolis_chlorocyanus", "Anolis_coelestinus", "Anolis_hendersoni", "Anolis_cybotes", "Anolis_angusticeps", "Anolis_angusticeps"))
Q5: 
fossilNodes <- c()
nodeN <- c()
{
for(i in 1:nrow(fossilData)) 
if (i == 1) {
i <- 1 
print(Ancestors)
}
}
Node <- fastMRCA(tree, fossilData[i, "tip1"], fossilData[i, "tip2"])
fossilNodes[i] <- fossilData[i, "sv1"]
nodeN [i] <- Node
names(fossilNodes) <- nodeN
AncestorswithFossils <- fastAnc(tree, sv1, anc.states=fossilNodes[i], CI=TRUE,var=TRUE)
AncestorswithFossils
Q7: Because fossils allow scientists to fill in some gaps under certain circumstances, i.e. with a reference fossil, the estimated ancestral size would increase.
Q8-Q10:
install.packages("geiger")
library("geiger")
?fitContinuous()
?fastAnc
BM <-fitContinuous(tree, sv1, model="BM")
BM
OU <- fitContinuous(tree, sv1, model="OU")
OU
EB <- fitContinuous(tree, sv1, model="EB")
Q: AIC_BM= -6.5 // AIC_OU= -4.5 // AIC_EB= -6.9 
Q: EB is the best model with the lowest AIC value. 
?fastAnc
Q: fastAnc relies on more assumptions than fitContinuous, and so I would think the slower method of fitContinuous is the better option as it takes into account continual data. 