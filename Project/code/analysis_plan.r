setwd('C:\\Users\\haile\\Desktop\\Evolution\\Tasks\\Project')
Data <- read.csv('C:\\Users\\haile\\Desktop\\Evolution\\Tasks\\Hypothesis_02\\evohyp.csv', stringsAsFactors=F)
write.csv(Data,'rawdata.csv',quote=F)
head(Data)
Data[3:131, 15]
Expected <- Data[3:131, 15]
Data[3:131, 14]
Observed <- Data[3:131, 14]
lm(Expected~Observed)
lm <- lm(Expected~Observed)
plot(lm)
plot(Expected, Observed)
Plot <- plot(Expected, Observed)


setwd('C:\\Users\\haile\\Desktop\\Evolution\\Tasks\\Project')
Data <- read.csv('C:\\Users\\haile\\Desktop\\Evolution\\Tasks\\Project\\data\\evohyp.csv', stringsAsFactors=F)
HPVStrain <- Data[4:12, 3]
NumOccurences <- Data[4:12, 4]
plot (Data$HPVStrain, Data$NumOccurences, xlim=c (0,60), ylim=c(0,2600))
barplot(NumOccurences~HPVStrain, data=Data)

setwd('C:\\Users\\haile\\Desktop\\Evolution\\Tasks\\Project\\Data')
Data <- read.csv('evohyp2_toread.csv', stringsAsFactors=F, row.names=1)

head(Data)
Vectored <- as.vector(Data)
muA <- length(grep("a", Vectored))
muT <- length(grep("t", Vectored))
muC <- length(grep("c", Vectored))
muG <- length(grep("g", Vectored))
muSH <- length(grep("sh", Vectored))

par(las=1, mgp=c(2.5, 0.5, 0), tck=-0.01, mar=c(4,5,1,1))
barplot(c(muA, muT, muC, muG, muSH), names.arg=c("a", "t", "c", "g", "del."), col='white', ylim=c(0, 300), ylab="# of instances")

Amuts <- apply(Data, 1, function(x) length(grep("a", x)))
Tmuts <- apply(Data, 1, function(x) length(grep("t", x)))
Cmuts <- apply(Data, 1, function(x) length(grep("c", x)))
Gmuts <- apply(Data, 1, function(x) length(grep("g", x)))
SHmuts <- apply(Data, 1, function(x) length(grep("sh", x)))

Pch <- 16
Col <- "red"
Type <- "p"
par(mfrow=c(3,2), las=1, mgp=c(2.5, 0.5, 0), tck=-0.01, mar=c(4,5,1,1))
plot(as.numeric(names(Amuts)), Amuts, xlab="", ylab="a", type=Type, pch=Pch, col=Col)
plot(as.numeric(names(Tmuts)), Tmuts, xlab="", ylab="t", type=Type, pch=Pch, col=Col)
plot(as.numeric(names(Cmuts)), Cmuts, xlab="", ylab="c", type=Type, pch=Pch, col=Col)
plot(as.numeric(names(Gmuts)), Gmuts, xlab="", ylab="g", type=Type, pch=Pch, col=Col)
plot(as.numeric(names(SHmuts)), SHmuts, xlab="", ylab="del.", type=Type, pch=Pch, col=Col)

