setwd('C:\\Users\\haile\\Desktop\\Evolution\\Tasks\\Task_08')
library("phytools")
Q1-3:
trees <- list()
births <- c()
Fractions <- c()
for(i in 1:100) {
	births[i] <- runif(1)
	Fractions[i] <- runif(1)
	trees[[i]] <- pbtree(b = births[i], d = (births[i] * Fractions[i]), n = 100, nsim = 1)
}
trees[[i]]
plot(trees[[i]])
library("geiger")
?geiger
install.packages("TreeTools")
library("TreeTools")
y
?TreeTools
?sapply
tips <- sapply(trees, NTip)
tips2 <- log(as.numeric(tips[100]))
tips2
div <- bd.ms(trees[[i]])
div
trees2 <- list()
births2 <- c ()
Fractions2 <- c()
for (var in 1:100)
{
births2[i] <- runif(1)
Fractions2[i] <- runif(1)
trees2[[i]] <- pbtree(b = births2[i], d = (births2[i] * Fractions2[i]), n = 100, nsim = 1)
}
trees2[[i]]
plot(trees2[[i]])
tips2 <- sapply(trees2, NTip)
tips22 <- log(as.numeric(tips2[100]))
tips22
div2 <- bd.ms(trees2[[i]])
div2
trees3 <- list()
births3 <- c ()
Fractions3 <- c()
for (var in 1:100)
{
births3[i] <- runif(1)
Fractions3[i] <- runif(1)
trees3[[i]] <- pbtree(b = births3[i], d = (births3[i] * Fractions3[i]), n = 100, nsim = 1)
}
trees3[[i]]
plot(trees3[[i]])
tips3 <- sapply(trees3, NTip)
tips32 <- log(as.numeric(tips3[100]))
tips32
div3 <- bd.ms(trees3[[i]])
div3
trees4 <- list()
births4 <- c ()
Fractions4 <- c()
for (var in 1:100)
{
births4[i] <- runif(1)
Fractions4[i] <- runif(1)
trees4[[i]] <- pbtree(b = births4[i], d = (births4[i] * Fractions4[i]), n = 100, nsim = 1)
}
trees4[[i]]
plot(trees4[[i]])
tips4 <- sapply(trees4, NTip)
tips42 <- log(as.numeric(tips4[100]))
tips42
div4 <- bd.ms(trees4[[i]])
div4
ftips<- c(tips2, tips22, tips32, tips42)
fdiv <- c(div, div2, div3, div4)
plot(fdiv, ftips, xlab= 'Net Diversification', ylab= 'Tips')
Q4: No correlation. A loose positive correlation at best.
speciation <- (bd.km(trees[[i]]))
speciation 
plot(speciation)
BranchL[[i]] <- mean(y*as.numeric(trees[100]))
Q5: Speciation is around .3 lower than average branch length
Q6: cor(speciation, BranchL)
LTREE <- trees[[i]]
rates <- c()
traits <- list()
for (i in 1:100)
{
rates[i] <- runif(1)
traits[[i]] <- fastBM(tree = LTREE, sig2 = rates[i])
}
AVGT <- sapply(traits, mean)
AVGT
AVGR <- sapply(rates, mean)
Q8: COR <- cor(AVGT, AVGR)
print(COR)
plot(AVGR~AVGT)
abline(lm(AVGR~AVGT), col="blue")
Q9: VART <- sapply(traits, var)
cor(VART, rates)
Q10: T1 <- traits[1]
T1
T2 <- traits[2]
T3 <- cbind(traits[[1]], traits[[2]])
T3
var(T3)
cor(T3[,1], T3[,2])
plot(T3[,1], T3[,2])
abline (lm(T3[,1]~T3[,2]), col="blue")