source("http://jonsmitchell.com/code/fxn05.R")
Pop1 <- simPop(Popsize = 50, nGenerations = 100, initial_p = 0.5, h= 1, s = 0)
plot (1:nrow(Pop1), Pop1 [,1], ylim=c(0,1), type = "l", xlab="generation", ylab="allele freq.", lwd=2)
lines(1:nrow(Pop1), Pop1[,2], lwd=2, col='red')
legend("topleft", legend = c("a", "b"), col = c("black", "red"), lwd = 2, bty="n")
plotFit( nruns = 10, n = 50, ngens = 100, init_p = 0.5, h = 1, s = 0)
Expectation <- c(10, 10, 10, 10)
Observed <- c(15, 15, 5, 5)
Chisq <- sum(((Expectation - Observed) ^2)/Expectation)
barplot(rbind(Expectation, Observed), beside =T, main = bquote(chi^2~ "=" ~. (Chisq)), legend.text=c("expected" , "observed"))
Expectation <- c(10, 10, 10, 10)
Observed <- c(5, 0, 0, 35)
Chisq <- sum(((Expectation - Observed) ^2)/Expectation)
barplot(rbind(Expectation, Observed), beside =T, main = bquote(chi^2~ "=" ~. (Chisq)), legend.text=c("expected" , "observed"))
Expectation <- c(10, 10, 10, 10)
Observed <- c(2, 3, 10, 30)
Chisq <- sum(((Expectation - Observed) ^2)/Expectation)
barplot(rbind(Expectation, Observed), beside =T, main = bquote(chi^2~ "=" ~. (Chisq)), legend.text=c("expected" , "observed"))
Expectation <- c(10, 10, 10, 10)
Observed <- c(5, 5, 5, 5)
Chisq <- sum(((Expectation - Observed) ^2)/Expectation)
barplot(rbind(Expectation, Observed), beside =T, main = bquote(chi^2~ "=" ~. (Chisq)), legend.text=c("expected" , "observed"))
results <- read.csv("http://jonsmitchell.com/data/biol112labresults.csv", stringsAsFactors=F)
counts <- results[,c("yellow", "red", "green", "blue", "black", "tan")]
backgrounds <- c("White", "Red", "Yellow", "Green", "Blue", "Black")
backgroungCol <- c("white" , "#d53e4f", "#fee08b", "#abdda4", "#3288bd", "black")
calcChi (counts[1,])
Chisqs <- apply(counts, 1, calcChi)
setwd('C:\\Users\\haile\\Desktop\\Evolution\\Tasks\\Task_04')
results <- read.csv("http://jonsmitchell.com/data/biol112labresults.csv", stringsAsFactors=F)
counts <- results[,c("yellow", "red", "green", "blue", "black", "tan")]
backgrounds <- c("White", "Red", "Yellow", "Green", "Blue", "Black")
backgroungCol <- c("white" , "#d53e4f", "#fee08b", "#abdda4", "#3288bd", "black")
calcChi(counts[1,])
Chisqs <- apply(counts, 1, calcChi)
source("http://jonsmitchell.com/code/fxn05.R")
results <- read.csv("http://jonsmitchell.com/data/biol112labresults.csv", stringsAsFactors=F)
counts <- results[,c("yellow", "red", "green", "blue", "black", "tan")]
backgrounds <- c("White", "Red", "Yellow", "Green", "Blue", "Black")
backgroundCol <- c("white" , "#d53e4f", "#fee08b", "#abdda4", "#3288bd", "black")
calcChi(counts[1,])
Chisqs <- apply(counts, 1, calcChi)
plotChis(counts)
plotChis(counts)
plotChis(counts)
plotChis(counts)
Avg <- mean(Chisqs)
backgroundAvgs <- tapply(Chisqs, results[,3], mean)
propSig <- length(which (Chisqs > 11.70) / length(Chisqs))
percSig <- round(100 * propSig)
par(las =1, mar= c(4, 4, 1, 1), mgp = c(2, 0.5, 0), tck =-0.01, cex.axis=1)
hist (Chisqs, main="", xlab="chi-squared values", ylab="frequency")
par (las = 1, mar = c(4, 4, 1, 1), mgp = c(2, 0.5, 0), tck = -0.01, cex.axis=1)
plot(1, 1, xlim=c(0, 400), ylim=c(1, 8.5), xlab="", ylab="", type="n", yaxt="n")
axis(2, at = 1:length(backgrounds), labels = backgrounds)
mtext(side=1, expression(chi^2), cex=1.75, line=2.5)
counter <- 1 
for (i in backgrounds) {
Data <- Chisqs [which(results [,3] == i)]
addHist(Y=counter, Dat=Data, Color= backgroundCol[counter])
counter <- counter + 1
}
abline (v= 11.70, lty=2, lwd = 2, col= 'black')
Simulation <- simDraws(10000)
addHist(Y=7, Dat=Simulation, Color="lightgray")
mtext(side=2, at=7, line=0, lwd=2)
Fit <- c(1, 1, 1, 1, 1, 1)
names(Fit) <- 1:6
Simulation2 <- simDraws(1e4, w = Fit)
addHist(Y=8, Dat=Simulation2, Color=rgb(0, 0, 0, 0.25))
Fit <- c(0.1, 1, 1, 1, 1, 1)
names(Fit) <- 1:6
Simulation3 <- simDraws(1e4, w = Fit)
addHist(Y=8, Dat=Simulation3, Color=rgb(0,0,0,0.25))
Fit <- c(0.5, 0.6, 0.7, 1, 1, 1)
names(Fit) <- 1:6
Simulation4 <- simDraws(1e4, w = Fit)
addHist (Y=8, Dat=Simulation4, Color=rgb(0, 0, 0, 0.25))
Fit <- c(0.1, 0.2, 0.3, 0.4, 0.5, 1)
names(Fit) <- 1:6
Simulation5 <- simDraws(1e4, w = Fit)
addHist (Y=8, Dat=Simulation4, Color=rgb(0, 0, 0, 0.25))
Fit <- c(0.1, 0.1, 0.1, 0.1, 0.1, 1)
names(Fit) <- 1:6
Simulation6 <- simDraws(1e4, w = Fit)
addHist (Y=8, Dat=Simulation4, Color=rgb(0, 0, 0, 0.25))
mtext(side=2, at=8, line=0, sel. sim."
Simulation7 <- c(Simulation2, Simulation3, Simulation4, Simulation5, Simulation6)
addHist(Y=8, Dat=Simulation7, Color=rgb(0, 0, 1, 0.25))
simDraws
simPop
simPop1 <-function(Popsize=100, nGenerations=100, h=1, s=0, initial_p=0.5, mu = 0, twoway = TRUE, w = NULL)       {
        if (is.null(w)) {
                # make fitness vector using h & s values
                w <- c(aa = 1, ab = 1 - h * s, bb = 1 - s)
        }
        else if (!is.null(w))   {
                if (length(w) != 3)     {
                        stop(cat("Not enough fitness values supplied. Set w to NULL or provide 3 fitnesses"))
                }
                if (is.null(names(w)))  {
                        names(w) <- c("aa", "ab", "bb")
                        cat("Setting names of w to aa, ab, and bb")
                }
        }
        # initialize the simulation with a starting population
        All_Individuals <- makePop(Popsize, intP = initial_p)

        # set up storage objects for the simulation output
        freqa <- c()
        freqa[1] <- sum(All_Individuals == "a") / length(All_Individuals)

        freqb <- c()
        freqb[1] <- sum(All_Individuals == "b") / length(All_Individuals)
        
        Gens <- apply(All_Individuals, 1, function(x) paste(sort(x), collapse=""))
        Genotypes <- matrix(0, nrow=nGenerations, ncol=3)
        colnames(Genotypes) <- names(w)
        GenCounts <- tapply(Gens, Gens, length)
        Genotypes[1,names(GenCounts)] <- GenCounts
        
        for (generation in 2:nGenerations)      {
                # Last generation's kids are this generation's parents!
                # What are their fitnesses?
                Fitnesses <- apply(All_Individuals, 1, checkFitness, w=w)

                # Half of the kids are female...
                Women <- sample(1:Popsize, size=(Popsize / 2))
                Mothers <- All_Individuals[Women,]

                # ...the other half are male
                Men <- setdiff(1:Popsize, Women)
                Fathers <- All_Individuals[Men,]

                # Now we do what we did before!
                nIndividuals <- Popsize
                for (child in 1:nIndividuals)   {
                        # When using sample(), you can actually SET the probability that
                        #       a given individual is sampled. 
                        father <- sample(1:nrow(Fathers), size=1, prob=Fitnesses[Men])
                        mother <- sample(1:nrow(Mothers), size=1, prob=Fitnesses[Women])

                        Child <- mating(Mothers[mother,], Fathers[father,], mu = mu, twoway = twoway)

                        cGenotype <- paste(Child[1], Child[2], sep="")
                        Genotypes[generation, cGenotype] <- Genotypes[generation, cGenotype] + 1
                        All_Individuals[child,] <- Child
                }

                nA_Individuals <- length(which(All_Individuals == "a"))
                freqa[generation] <- nA_Individuals / length(All_Individuals)

                na_Individuals <- length(which(All_Individuals == "b"))
                freqb[generation] <- na_Individuals / length(All_Individuals)

                if (max(c(freqa[generation],freqb[generation])) == 1 && mu == 0)        {
                        Diff <- nGenerations - generation
                        freqa <- c(freqa, rep(freqa[generation], Diff))
                        freqb <- c(freqb, rep(freqb[generation], Diff))
                        Genotypes <- rbind(Genotypes, matrix(Genotypes[generation,], ncol=ncol(Genotypes), nrow=Diff, byrow= TRUE))
                        break;
                }
        }
        output <- data.frame(freqa=freqa, freqb=freqb, Genotypes=Genotypes[1:length(freqa),])
        return(output)
}
simPop0.5 <-function(Popsize=100, nGenerations=100, h=0.5, s=0, initial_p=0.5, mu = 0, twoway = TRUE, w = NULL)       {
        if (is.null(w)) {
                # make fitness vector using h & s values
                w <- c(aa = 1, ab = 1 - h * s, bb = 1 - s)
        }
        else if (!is.null(w))   {
                if (length(w) != 3)     {
                        stop(cat("Not enough fitness values supplied. Set w to NULL or provide 3 fitnesses"))
                }
                if (is.null(names(w)))  {
                        names(w) <- c("aa", "ab", "bb")
                        cat("Setting names of w to aa, ab, and bb")
                }
        }
        # initialize the simulation with a starting population
        All_Individuals <- makePop(Popsize, intP = initial_p)

        # set up storage objects for the simulation output
        freqa <- c()
        freqa[1] <- sum(All_Individuals == "a") / length(All_Individuals)

        freqb <- c()
        freqb[1] <- sum(All_Individuals == "b") / length(All_Individuals)
        
        Gens <- apply(All_Individuals, 1, function(x) paste(sort(x), collapse=""))
        Genotypes <- matrix(0, nrow=nGenerations, ncol=3)
        colnames(Genotypes) <- names(w)
        GenCounts <- tapply(Gens, Gens, length)
        Genotypes[1,names(GenCounts)] <- GenCounts
        
        for (generation in 2:nGenerations)      {
                # Last generation's kids are this generation's parents!
                # What are their fitnesses?
                Fitnesses <- apply(All_Individuals, 1, checkFitness, w=w)

                # Half of the kids are female...
                Women <- sample(1:Popsize, size=(Popsize / 2))
                Mothers <- All_Individuals[Women,]

                # ...the other half are male
                Men <- setdiff(1:Popsize, Women)
                Fathers <- All_Individuals[Men,]

                # Now we do what we did before!
                nIndividuals <- Popsize
                for (child in 1:nIndividuals)   {
                        # When using sample(), you can actually SET the probability that
                        #       a given individual is sampled. 
                        father <- sample(1:nrow(Fathers), size=1, prob=Fitnesses[Men])
                        mother <- sample(1:nrow(Mothers), size=1, prob=Fitnesses[Women])

                        Child <- mating(Mothers[mother,], Fathers[father,], mu = mu, twoway = twoway)

                        cGenotype <- paste(Child[1], Child[2], sep="")
                        Genotypes[generation, cGenotype] <- Genotypes[generation, cGenotype] + 1
                        All_Individuals[child,] <- Child
                }

                nA_Individuals <- length(which(All_Individuals == "a"))
                freqa[generation] <- nA_Individuals / length(All_Individuals)

                na_Individuals <- length(which(All_Individuals == "b"))
                freqb[generation] <- na_Individuals / length(All_Individuals)

                if (max(c(freqa[generation],freqb[generation])) == 1 && mu == 0)        {
                        Diff <- nGenerations - generation
                        freqa <- c(freqa, rep(freqa[generation], Diff))
                        freqb <- c(freqb, rep(freqb[generation], Diff))
                        Genotypes <- rbind(Genotypes, matrix(Genotypes[generation,], ncol=ncol(Genotypes), nrow=Diff, byrow= TRUE))
                        break;
                }
        }
        output <- data.frame(freqa=freqa, freqb=freqb, Genotypes=Genotypes[1:length(freqa),])
        return(output)
}
simPop2
simPop0.5
simDraws1 <-
function(nruns, ncols=6, nstart=10, nrounds=3, mu=1, w=NULL)  {
        Chiout <- c()
        for (j in 1:nruns)      {
                Start <- rep(1:ncols, nstart)
                Pop <- Start
                for (i in 1:nrounds)    {
                        if (is.null(w)) {
                                Draws <- sample(Pop, 20, replace = F)
                        }
                        else if (!is.null(w))   {
                                if (length(setdiff(unique(Pop), names(w))) == 0)        {
                                        Draws <- sample(Pop, 20, replace=F, prob=w[Pop])
                                }
                                else if (length(setdiff(unique(Pop), names(w))) != 0)   {
                                        cat("Not enough fitness values! ", setdiff(unique(Pop), names(w)))
                                }
                        }
                        Pop <- sort(c(Draws,Draws,Draws))
                }
                Summary <- c()
                for (k in 1:ncols)      {
                        Summary[k] <- length(which(Pop == k))
                }

                Chiout[j] <- sum(((Summary - nstart)^2) / nstart)
        }
        return(Chiout)
}
<bytecode: 0x0ae88068
simDraws0.5 <-
function(nruns, ncols=6, nstart=10, nrounds=3, mu=0.5, w=NULL)  {
        Chiout <- c()
        for (j in 1:nruns)      {
                Start <- rep(1:ncols, nstart)
                Pop <- Start
                for (i in 1:nrounds)    {
                        if (is.null(w)) {
                                Draws <- sample(Pop, 20, replace = F)
                        }
                        else if (!is.null(w))   {
                                if (length(setdiff(unique(Pop), names(w))) == 0)        {
                                        Draws <- sample(Pop, 20, replace=F, prob=w[Pop])
                                }
                                else if (length(setdiff(unique(Pop), names(w))) != 0)   {
                                        cat("Not enough fitness values! ", setdiff(unique(Pop), names(w)))
                                }
                        }
                        Pop <- sort(c(Draws,Draws,Draws))
                }
                Summary <- c()
                for (k in 1:ncols)      {
                        Summary[k] <- length(which(Pop == k))
                }

                Chiout[j] <- sum(((Summary - nstart)^2) / nstart)
        }
        return(Chiout)
}
<bytecode: 0x0ae88068
source("http://jonmitchell.com/code/fxn05.R")
results <- read.csv("http://jonsmitchell.com/data/biol112labresults.csv", stringsAsFactors=F)
counts <- results[,c("yellow", "red", "green", "blue", "black", "tan")]
backgrounds <- c("White", "Red", "Yellow", "Green", "Blue", "Black")
backgroungCol <- c("white" , "#d53e4f", "#fee08b", "#abdda4", "#3288bd", "black")
calcChi (counts[1,])
Chisqs <- apply(counts, 1, calcChi)
plotChis(counts)
Avg <- mean(Chisqs)
backgroundAvgs <- tapply(Chisqs, results[,3], mean)
propSig <- length(which (Chisqs > 11.70) / length(Chisqs))
percSig <- round(100 * propSig)
par(las =1, mar= c(4, 4, 1, 1), mgp = c(2, 0.5, 0), tck =-0.01, cex.axis=1)
hist (Chisqs, main="", xlab="chi-squared values", ylab="frequency")
par (las = 1, mar = c(4, 4, 1, 1), mgp = c(2, 0.5, 0), tck = -0.01, cex.axis=1)
plot(1, 1, xlim=c(0, 400), ylim=c(1, 8.5), xlab="", ylab="", type="n", yaxt="n")
axis(2, at = 1:length(backgrounds), labels = backgrounds)
mtext(side=1, expression(chi^2), cex=1.75, line=2.5)
counter <- 1 
for (i in backgrounds) {
Data <- Chisqs [which(results [,3] == i)]
addHist(Y=counter, Dat=Data, Color= backgroundCol[counter])
counter <- counter + 1
}
abline (v= 11.70, lty=2, lwd = 2, col= 'black')
SimulationN <- simDraws1(10000)
addHist(Y=7, Dat=Simulation, Color="lightgray")
mtext(side=2, at=7, line=0, lwd=2)
Fit <- c(1, 1, 1, 1, 1, 1)
names(Fit) <- 1:6
Simulation8 <- simDraws1(1e4, w = Fit)
addHist(Y=8, Dat=Simulation8, Color=rgb(0, 0, 0, 0.25))
Fit <- c(0.1, 1, 1, 1, 1, 1)
names(Fit) <- 1:6
Simulation9 <- simDraws1(1e4, w = Fit)
addHist(Y=8, Dat=Simulation9, Color=rgb(0,0,0,0.25))
Fit <- c(0.5, 0.6, 0.7, 1, 1, 1)
names(Fit) <- 1:6
Simulation10 <- simDraws1(1e4, w = Fit)
addHist (Y=8, Dat=Simulation10, Color=rgb(0, 0, 0, 0.25))
Fit <- c(0.1, 0.2, 0.3, 0.4, 0.5, 1)
names(Fit) <- 1:6
Simulation11 <- simDraws1(1e4, w = Fit)
addHist (Y=8, Dat=Simulation11, Color=rgb(0, 0, 0, 0.25))
Fit <- c(0.1, 0.1, 0.1, 0.1, 0.1, 1)
names(Fit) <- 1:6
Simulation12 <- simDraws1(1e4, w = Fit)
addHist (Y=8, Dat=Simulation12, Color=rgb(0, 0, 0, 0.25))
mtext(side=2, at=8, line=0, sel. sim.") 
Simulation13 <- c(Simulation8, Simulation9, Simulation10, Simulation11, Simulation12)
addHist(Y=8, Dat=Simulation7, Color=rgb(0, 0, 1, 0.25)) 
plot( Simulation2)
plot( Simulation8)
plot(Simulation7)
plot(Simulation13)