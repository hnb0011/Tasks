setwd('C:\\Users\\haile\\Desktop\\Evolution\\Tasks\\Task_06')
source ("http://jonsmitchell.com/code/reformatData07.R")
source ("http://jonsmitchell.com/code/simFxn.R")
dim(overallFreq)
nrow(overallFreq)
plot (1, 1, type="n", xlim=c(1998, 2013), ylim=c(0, 1))
s <- apply(overallFreq, 2, function(x) lines(overallFreq[,1], x, col= rgb(0,0,0,0.01)))
rescaleFreq <- apply(overallFreq [,3:ncol(overallFreq)], 2, function(x) x - x[1])
plot(1, 1, type ="n", xlim=c (1998, 2013), ylim=c(-0.25, 0.25))
s <- apply (rescaleFreq, 2, function(x) lines (overallFreq[,1] ,x, col=rgb(0,0,0,0.01)))
dYear <- c()
dAlleles <- c()
for (i in 3:ncol (overallFreq)) {
dYear <- c(dYear, overallFreq [,1])
Vec <- overallFreq [,i]
Init <- overallFreq [1,i]
dAlleles <- c(dAlleles, Vec - Init)
}
smoothScatter (dYear, dAlleles, colramp = Pal, nbin = 100)
smoothScatter (dYear, dAlleles, colramp = Pal, nbin = 500)
smoothScatter (dYear, dAlleles, colramp = Pal, nbin = 10)
smoothScatter (dYear dAlleles, colramp = Pal, nbin = 100, xlab="year", 
ylab="change in allele freq. since 1998"
addFit (nruns = 50, n = 100, ngens = 18, startT = 1997, simCol = "gray40", rescale =TRUE)
addFit (nruns = 50, n = 10, ngens = 18, startT = 1997, simCol = "red", rescale =TRUE)
addFit (nruns = 50, n = 1000, ngens = 18, startT = 1997, simCol = "green", rescale =TRUE)
addFit (nruns = 50, n = 926, ngens = 18, startT = 1997, simCol = "green", rescale =TRUE)
addFit (nruns = 50, n = 926, h=0.7, s=0.05, ngens = 18, startT = 1997, simCol = "green", rescale =TRUE)
addFit (nruns = 50, n = 926, h=0.7, s=0.025, ngens = 18, startT = 1997, simCol = "black", rescale =TRUE)
Q: my chosen values are h=0.7, s=0.025, and n=926
plot (alleleFreqs$d_freq, alleleFreqs$d_imm, xlim=c(-0.15, 0.15), xlab="
overall freq. change", ylab="freq. change in subset")
points(alleleFreqs$d_freq, alleleFreqs$d_birth, col= 'blue')
points(alleleFreqs$d_freq, alleleFreqs$d_surv, col ='red')
setwd("C:\\Users\\haile\\Desktop\\Evolution\\Tasks\\Task_06")
today<-format(Sys.Date(),format="%d%b%Y")
install.packages ("plyr")
install.packages ("doParallel")
library(plyr)
library(doParallel)
snplist<-read.table('SNPlist.r',header=TRUE,stringsAsFactors=FALSE)
download.packages ('foreach')
download.packages ('doParallel')
library(foreach)
library(doParallel)
nloci<-100000
indivlist<-read.table('IndivList.r',header=TRUE,sep="\t",stringsAsFactors=FALSE)
plot (snplist$r_freq, indivlist$r_freq, type ="n", xlim=c (1998, 2013), ylim=c(-0.25, 0.25))
s <- apply (rescaleFreq, 2, function(x) lines (overallFreq[,1] ,x, col=rgb(0,0,0,0.01)))