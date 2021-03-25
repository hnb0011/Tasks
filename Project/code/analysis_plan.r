setwd('C:\\Users\\haile\\Desktop\\Evolution\\Tasks\\Hypothesis_02')
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


setwd('C:\\Users\\haile\\Desktop\\Evolution\\Tasks\\Hypothesis_02')
Data <- read.csv('C:\\Users\\haile\\Desktop\\Evolution\\Tasks\\Hypothesis_02\\evohyp.csv', stringsAsFactors=F)
HPVStrain <- Data[4:12, 3]
NumOccurences <- Data[4:12, 4]
plot (Data$HPVStrain, Data$NumOccurences, xlim=c (0,60), ylim=c(0,2600))
barplot(NumOccurences~HPVStrain, data=Data)