setwd('C:\\Users\\haile\\Desktop\\Evolution\\Tasks\\Task_02')
Data <- read.csv ('http://jonsmitchell.com/data/beren.csv', stringsAsFactors=F)
write.csv(Data, 'rawdata.csv', quote=F)
length(Data)
nrow(Data)
ncol(Data)
colnames(Data)
head(Data)
Data [1,]
Data [2,]
Data [1:3, 4]
Data[1:5, 1:3]
Data [257,3]
Feeds <- which(Data [,9] == 'bottle')
berenMilk <- Data[Feeds,]
head(berenMilk)
nrow(berenMilk)
Feeds <- which(Data[,'event'] == 'bottle')
nrow(berenMilk)
Feeds <- which(Data$event == 'bottle')
nrow(berenMilk)
Feeds == which(Data$event == 'bottle')
which(Data$event == 'bottle') == which(Data[,'event'] == 'bottle')
Feeds == which(Data[,'event'] == 'bottle')
dayID <- apply(Data, 1, function(x) paste (x[1:3], collapse='-'))
dateID <- sapply(dayID, as.Date, format = "&Y-%m-%d", origin="2019-04-18")
Data$age <- dateID - dateID [which(Data$event == 'birth')]
head(Data)
beren2 <- Data
beren3 <- beren2[order(beren2$age) ,]
head(beren2)
head(beren3)
write.csv(beren3, 'beren_new.csv', quote=F, row.names=FALSE)
head(beren3)
Q1: HI does not work because no weight data is recorded
Q1: HII does not work because a nap and a bottle are both events not independent of each other in the data
Feeds <- which(beren3$event == "bottle")
avgMilk <- mean(beren3$value [Feeds])
avgMilk
avgFeed <- tapply(beren3$value[Feeds], beren3$age[Feeds], mean)
head(avgMilk)
head(avgFeed)
length(avgFeed)
nrow(avgFeed)
ncol(avgFeed)
beren <- read.csv("http://jonsmitchell.com/data/beren.csv", stringsAsFactors=F)
dayID <- apply(beren, 1, function(x) paste(x[1:3], collapse="-"))
dateID <- sapply(dayID, as.Date, format = "%Y-%m-%d", origin = "2019-04-18")
beren$age <- dateID - dateID[which(beren$event == "birth")]
beren2 <- beren[order(beren$age),]
beren2$value <- as.numeric(beren2$value)
head(Feeds)
varFeed <- tapply(beren3$value [Feeds], beren3$age[Feeds], var)
totalFeed <- tapply(beren3$value[Feeds], beren3$age[Feeds], sum)
numFeeds <- tapply(beren3$value[Feeds], beren3$age[Feeds], length)
cor(beren3$value[Feeds], beren3$age[Feeds])
cor.test(beren3$value[Feeds], beren3$age[Feeds])
head(Feeds)
setwd('C:\\Users\\haile\\Desktop\\Evolution\\Tasks\\Task_02')
beren <- read.csv("http://jonsmitchell.com/data/beren.csv", stringsAsFactors=F)
write.csv(beren, 'rawdata.csv', quote=F)
length(beren)
nrow(beren)
ncol(beren)
colnames(beren)
head(beren)
beren [1,]
beren [2,]
beren [1:3, 4]
beren[1:5, 1:3]
beren [257,3]
Feeds <- which(beren[,9] == 'bottle')
berenMilk <- beren[Feeds,]
head(berenMilk)
Feeds <- which(beren[, 'event'] == 'bottle')
Feeds <- which(beren$event == 'bottle')
dayID <- apply(beren, 1, function(x) paste(x[1:3], collapse="-"))
dateID <- sapply(dayID, as.Date, format = "%Y-%m-%d", origin = "2019-04-18")
beren$age <- dateID - dateID[which(beren$event == "birth")]
head(beren)
beren2 <- beren
beren2 <- beren[order(beren$age),]
beren2$value <- as.numeric(beren2$value)
write.csv(beren2, 'beren_new.csv', quote=F, row.names=FALSE)
Feeds <- which(beren2$event == "bottle")
avgMilk <- mean(beren2$value[Feeds])
avgFeed <- tapply(beren2$value[Feeds], beren2$age[Feeds], mean)
head(avgFeed)
varFeed <- tapply(beren2$value [Feeds], beren2$age[Feeds], var)
totalFeed <- tapply(beren2$value[Feeds], beren2$age[Feeds], sum)
numFeeds <- tapply(beren2$value[Feeds], beren2$age[Feeds], length)
cor(beren2$value [Feeds], beren2$age[Feeds])
cor.test(beren2$value[Feeds], beren2$age[Feeds])
berenCor <- cor.test (beren2$value[Feeds], beren2$age[Feeds])
summary(berenCor)
berenANOVA <- aov(beren2$value[Feeds] ~ beren2$caregiver[Feeds])
boxplot(beren2$value[Feeds] ~ beren2$caregiver[Feeds], xlab= "who gave the bottle", ylab= "amount of milk consumed (oz)")
plot(as.numeric(names(totalFeed)), totalFeed, type="b", pch=16, xlab="age in days", ylab="ounces of milk")
abline(h=mean(totalFeed), lty=2, col='red')
pdf("r02b-totalMilkByDay.pdf", height = 4, width = 4)
par(las=1, mar=c(5,5,1,1), mgp=c(2,0.5, 0), tck=-0.01)
plot(as.numeric(names(totalFeed)), totalFeed, type="b", pch=16, xlab="age in days", ylab="ounches of milk")
abline(h=mean(totalFeed), lty=2, col='red')
dev.off()
Q2: It is hard to interpret because gridlines are lacking and the x and y axes' labels are far apart. 
source("http://jonsmitchell.com/code/plotFxn02b.R")
unique(beren2$event)
beren3 <- beren2
pdf("r02b-cumulativeMilkByDay.pdf", height = 4, width = 4)
Naps <- which(beren3$event == "nap")
Naps
beren4 <- beren3[Naps, ]
head(beren4)
startID <- apply(beren3, 1, function(x) paste(x[5:6], collapse='-'))
endID <- apply(beren3, 1, function(x) paste(x[7:8], collapse='-'))
beren3$time <- endID - startID
startID2 <- na.omit(startID)
endID2 <- na.omit(endID)
startID <- apply(beren4, 1, function(x) paste(x[5:6], collapse=':'))
endID <- apply(beren4, 1, function(x) paste(x[7:8], collapse=':'))
beren4$time <- endID - startID[which(beren4$event == 'nap')]
timeID <- sapply(endID - startID, format = "%H:%M", origin = "12:30")
timeID <- sapply(endID, as.numeric, format = "%H:%M", origin = "12:30")
beren4$time <- as.numeric(timeID2) - as.numeric(startID2)[which(beren4$event == 'nap')]
timeID2 <- na.omit(timeID, na.action = "omit", fill = 0)
startID
endID
startID <- as.POSIXlt(startID, format='%H:%M')
endID <- as.POSIXlt(endID, format='%H:%M')
endID-startID
totalTime <- tapply(beren4$value[Naps], beren4$time[Naps], sum)
timeline <- c(0, 10)
plot.window(xlim=timeline, ylim=timeline)
plot(as.numeric(names(totalTime)), totalTime, type="b", pch=16, xlab="Day", ylab="Total-Time-Slept")
plot(totalTime, type="b", pch=16, xlab="Day", ylab="Total-Time-Slept", xlim=timeline, ylim=timeline)
timeline <- c(0, 100)
plot(totalTime, type="b", pch=16, xlab="Day", ylab="Total-Time-Slept", xlim=timeline, ylim=timeline)
timeline <- c(0, 100)
totalT <- endID-startID
cor.test(beren3$time[Naps,], beren3$start[Naps,])
beren3$start <- beren3[5:6]
cor.test(totalT[Naps,], startID[Naps,])
cor.test(totalT[Naps,], startID)
cor.test(totalT, startID)



