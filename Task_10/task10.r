setwd('C:\\Users\\haile\\Desktop\\Evolution\\Tasks\\Task_10')
install.packages("diversitree")
library("diversitree")
transition_0to1 <- 0.1
transition_1to0 <- 0.1
speciation_0 <- 0.2
extinction_0 <- 0.1
speciation_1 <- 0.2
extinction_1 <- 0.1
maxN <- 1e3
maxT <- 50
Pars <- c(speciation_0, speciation_1, extinction_0, extinction_1, transition_0to1, transition_1to0)
simTree <- tree.bisse(Pars, max.taxa = maxN, max.t = maxT)
plot(simTree, cex=0.01)
str(simTree)
?tree.bisse
stateTable <- table(simTree$tip.state)
stateTable / sum(stateTable)
simTree_2 <- tree.bisse(Pars, max.taxa = maxN, max.t = maxT)
plot(simTree_2, cex=0.01)
str(simTree_2)
stateTable_2 <- table(simTree_2$tip.state)
stateTable_2 / sum(stateTable_2)
simTree_3 <- tree.bisse(Pars, max.taxa = maxN, max.t = maxT)
plot(simTree_3, cex=0.01)
str(simTree_3)
stateTable_3 <- table(simTree_3$tip.state)
stateTable_3 / sum(stateTable_3)
simTree_4 <- tree.bisse(Pars, max.taxa = maxN, max.t = maxT)
plot(simTree_4, cex=0.01)
str(simTree_4)
stateTable_4 <- table(simTree_4$tip.state)
stateTable_4 / sum(stateTable_4)
simTree_5 <- tree.bisse(Pars, max.taxa = maxN, max.t = maxT)
plot(simTree_5, cex=0.01)
str(simTree_5)
stateTable_5 <- table(simTree_5$tip.state)
stateTable_5 / sum(stateTable_5)
a <- stateTable / sum(stateTable)
b <- stateTable_2 / sum(stateTable_2)
c <- stateTable_3 / sum(stateTable_3)
d <- stateTable_4 / sum(stateTable_4)
e <- stateTable_5 / sum(stateTable_5)
scatterP <- plot(c(a, b, c,d, e))


