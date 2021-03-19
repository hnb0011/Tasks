#script to model variance in allele frequency change over time
#run sims for 100,000 SNPs to get empirical estimates of covariances and errors
#Nancy Chen & Graham Coop
#Last updated: 5 July 2018

library(plyr)
library(foreach)
library(doParallel)

#number of SNPs to simulate
nloci<-100000

#get date script is run
today<-format(Sys.Date(),format="%d%b%Y")

#get input files: fixed list of indiv in each Category each year
indivlist<-read.table('IndivList.r',header=TRUE,sep="\t",stringsAsFactors=FALSE)
indivlist$Indiv<-as.character(indivlist$Indiv)

#get sample allele freq for simulations
load('indivlistgeno.rdata')
indivlistgeno$Indiv<-as.character(indivlistgeno$Indiv)
indivlistgeno$Dad<-as.character(indivlistgeno$Dad)
indivlistgeno$Mom<-as.character(indivlistgeno$Mom)

datafreq1990<-laply(names(indivlistgeno)[7:10737],function(x) 
	sum(indivlistgeno[indivlistgeno$Year==1990,x],na.rm=TRUE)/
	(2*sum(!is.na(indivlistgeno[indivlistgeno$Year==1990,x]))))

#sample allele freq
simfreq<-sample(datafreq1990,nloci,replace=TRUE)

#add parents
ped<-read.table('FSJpedgeno2018.ped',header=FALSE,sep=' ',stringsAsFactors=FALSE)
indivlist1<-merge(indivlist,ped[,c(2:4)],by.x='Indiv',by.y='V2')
indivlist1<-indivlist1[order(indivlist1$Year,as.numeric(indivlist1$Indiv)),]
names(indivlist1)<-c('Indiv','Year','Category','Genotyped','Dad','Mom')
indivlist1$Dad<-as.character(indivlist1$Dad)
indivlist1$Mom<-as.character(indivlist1$Mom)

#get unique indiv
simindivgeno<-indivlist1[!duplicated(indivlist1$Indiv),]

#separate into adults vs nestlings
simindivgenoAdults<-
	simindivgeno[simindivgeno$Category!='nestling' | simindivgeno$Mom=='0',]
	
simindivgenoNestlings<-
	simindivgeno[simindivgeno$Category=='nestling' & simindivgeno$Mom!='0',]

#simulate genotypes for adults
num.adults<-nrow(simindivgenoAdults)
adult.genos<-sapply(1:nloci,function(loc){
	freq<-simfreq[loc]
	HWE<-c(freq^2,2*freq*(1-freq),(1-freq)^2)	
	loc.genos<-sample(0:2,size=num.adults,prob=HWE,replace=TRUE)
	loc.genos
})

simindivgenoAdults<-cbind(simindivgenoAdults,adult.genos)

simindivgenoAll<-cbind(simindivgeno,
	simindivgenoAdults[match(simindivgeno$Indiv,simindivgenoAdults$Indiv),7:(nloci+6)])

rownames(simindivgenoAll)<-simindivgenoAll$Indiv

#simulate genotypes for nestlings via Mendelian transmission of alleles from parents
nest.years<-unique(simindivgenoNestlings$Year)

make.gametes<-function(g){
	gametes<-rep(NA,length(g))
	gametes[g==0]<-0
	gametes[g==2]<-1
	hets<- g==1 & !is.na(g)
	gametes[hets]<-sample(c(0,1),sum(hets),prob=c(0.5,0.5),replace=TRUE)
	gametes
}

for(year in nest.years){
	these.nestlings<-simindivgenoNestlings$Indiv[simindivgenoNestlings$Year==year]
	these.moms<-simindivgenoNestlings$Mom[simindivgenoNestlings$Year==year]
	these.dads<-simindivgenoNestlings$Dad[simindivgenoNestlings$Year==year]
	
	stopifnot(all(is.na(simindivgenoAll[these.nestlings,7:(nloci+6)])))
	
	mom.geno<-simindivgenoAll[these.moms,]
	dad.geno<-simindivgenoAll[these.dads,]

	stopifnot(all(!is.na(mom.geno[,7:(nloci+6)])))
	stopifnot(all(!is.na(dad.geno[,7:(nloci+6)])))

	Dad.gamete<-apply(dad.geno[,7:(nloci+6)],2,make.gametes)
	Mom.gamete<-apply(mom.geno[,7:(nloci+6)],2,make.gametes)
	
	nestling.geno<-Dad.gamete+Mom.gamete
	
	simindivgenoAll[these.nestlings,7:(nloci+6)] <- nestling.geno
}

simdataTrue<-merge(indivlist1,simindivgenoAll[,c(1,7:(nloci+6))],
	by.x='Indiv',by.y='Indiv',all.x=TRUE)	
	
#save
save(simdataTrue,file='simdataTrue.rdata')

#get number of all & genotyped indiv
counts<-ddply(indivlist,.(Year,Category),summarize,genotyped=2*sum(Genotyped=='Y'),
	total=2*length(Category))
	
countsAll<-ddply(indivlist,.(Year),summarize,genotyped=2*sum(Genotyped=='Y'),
	total=2*length(Category))

#calculate sample allele freq
#mimic sampling of genotyped indiv
simdataSample<-simdataTrue[simdataTrue$Genotyped=='Y',]

#get unique indiv
simdataTrueUnique<-simdataTrue[!duplicated(simdataTrue$Indiv),]
simdataSampleUnique<-simdataSample[!duplicated(simdataSample$Indiv),]

#calculate population (p) and sample allele freq (x)
#and the error in allele freq estimation due to sampling: err = x-p
#parallelize
simAlleleFreq<-data.frame(Year=integer(),category=character(),stringsAsFactors=F)

#parallelize snps
#cores=detectCores() #uncomment these two lines if you want to use more than 4 cores
#cl <- makeCluster(cores[1]-1) #not to overload your computer

cl <- makeCluster(4) #use 4 cores
registerDoParallel(cl)

#1998
sim<-foreach(i=names(simdataTrue)[7:(nloci+6)],.combine=cbind) %dopar% {
	tmp<-data.frame(Year=rep(year,each=4),Category=c('pt','xt','errT','Et'),
		stringsAsFactors=FALSE)
		
	frqYr1<-tmp$Year
	frqCat1<-tmp$category

	tmp[frqYr1==year & frqCat1=='pt',3]<-mean(simdataTrue[simdataTrue$Year==year,i])/2
	tmp[frqYr1==year & frqCat1=='xt',3]<-mean(simdataSample[simdataSample$Year==year,i])/2

	tmp[,3]
}
stopCluster(cl)

#save(sim,file=paste("SimAlleleFreqYr_",year,".rdata",sep=''))

simName<-data.frame(Year=rep(year,each=4),Category=c('pt','xt','errT','Et'),
	stringsAsFactors=FALSE)
	
sim1<-cbind(simName,sim)	
simAlleleFreq<-rbind(simAlleleFreq,sim1)

for(year in c(1999:2013))
{
	moms<-simdataTrue[simdataTrue$Year==year & simdataTrue$Category=='nestling','Mom']
	moms<-data.frame(Indiv=moms[!is.na(moms)],stringsAsFactors=FALSE)
	momgeno<-merge(moms,simdataTrueUnique[,c(1,7:(nloci+6))],by.x='Indiv',by.y='Indiv',
		all.x=TRUE)
		
	momgenoSample<-merge(moms,simdataSampleUnique[,c(1,7:(nloci+6))],by.x='Indiv',
		by.y='Indiv',all.x=TRUE)

	dads<-simdataTrue[simdataTrue$Year==year & simdataTrue$Category=='nestling','Dad']
	dads<-data.frame(Indiv=dads[!is.na(dads)],stringsAsFactors=FALSE)
	dadgeno<-merge(dads,simdataTrueUnique[,c(1,7:(nloci+6))],by.x='Indiv',by.y='Indiv',
		all.x=TRUE)
		
	dadgenoSample<-merge(dads,simdataSampleUnique[,c(1,7:(nloci+6))],by.x='Indiv',
		by.y='Indiv',all.x=TRUE)

	#parallelize snps
	cl <- makeCluster(4) #use 4 cores
	registerDoParallel(cl)

	sim<-foreach(i=names(simdataTrue)[7:(nloci+6)],.combine=cbind) %dopar% {
		tmp<-data.frame(Year=rep(year,each=41),Category=c('pt','xt','errT','Et','ps','xs',
			'errS','Es','pi','xi','errI','Ei','pb','xb','errB','Eb','pt1-pt','ps-pt',
			'xs-xt','errS-errT','pi-pt','xi-xt','errI-errT','pb-pt','xb-xt','errB-errT',
			'pm','xm','errM','pf','xf','errF','pfam','pmend','pfam-pt','xfam','xmend',
			'xfam-xt','errFAM','errMEND','errFAM-errT'),stringsAsFactors=FALSE)

		frqYr1<-tmp$Year
		frqCat1<-tmp$category

		tmp[frqYr1==year & frqCat1=='pt',3]<-mean(simdataTrue[simdataTrue$Year==year,i])/2
		
		tmp[frqYr1==year & frqCat1=='xt',3]<-
			mean(simdataSample[simdataSample$Year==year,i])/2

		tmp[frqYr1==year & frqCat1=='ps',3]<-mean(simdataTrue[simdataTrue$Year==year & 
			simdataTrue$Category=='survivor',i])/2
			
		tmp[frqYr1==year & frqCat1=='xs',3]<-mean(simdataSample[simdataSample$Year==year & 
			simdataSample$Category=='survivor',i])/2

		tmp[frqYr1==year & frqCat1=='pi',3]<-mean(simdataTrue[simdataTrue$Year==year & 
			simdataTrue$Category=='immigrant',i])/2
			
		tmp[frqYr1==year & frqCat1=='xi',3]<-
			ifelse(is.na(mean(simdataSample[simdataSample$Year==year & 
			simdataSample$Category=='immigrant',i])),0,
			mean(simdataSample[simdataSample$Year==year & 
			simdataSample$Category=='immigrant',i])/2)

		tmp[frqYr1==year & frqCat1=='pb',3]<-mean(simdataTrue[simdataTrue$Year==year & 
			simdataTrue$Category=='nestling',i])/2
			
		tmp[frqYr1==year & frqCat1=='xb',3]<-mean(simdataSample[simdataSample$Year==year &
			simdataSample$Category=='nestling',i])/2

		tmp[frqYr1==year & frqCat1=='pm',3]<-mean(momgeno[,i],na.rm=TRUE)/2
		tmp[frqYr1==year & frqCat1=='xm',3]<-mean(momgenoSample[,i],na.rm=TRUE)/2

		tmp[frqYr1==year & frqCat1=='pf',3]<-mean(dadgeno[,i],na.rm=TRUE)/2
		tmp[frqYr1==year & frqCat1=='xf',3]<-mean(dadgenoSample[,i],na.rm=TRUE)/2
		
		tmp[,3]
	}
	stopCluster(cl)
	
	#save(sim,file=paste("SimAlleleFreqYr_",year,".rdata",sep=''))
	
	simName<-data.frame(Year=rep(year,each=41),Category=c('pt','xt','errT','Et','ps','xs',
		'errS','Es','pi','xi','errI','Ei','pb','xb','errB','Eb','pt1-pt','ps-pt','xs-xt',
		'errS-errT','pi-pt','xi-xt','errI-errT','pb-pt','xb-xt','errB-errT','pm','xm',
		'errM','pf','xf','errF','pfam','pmend','pfam-pt','xfam','xmend','xfam-xt',
		'errFAM','errMEND','errFAM-errT'),stringsAsFactors=FALSE)
		
	sim1<-cbind(simName,sim)	
	simAlleleFreq<-rbind(simAlleleFreq,sim1)
}

#calculate error and allele freq differences between each category and the year before
#err = true error, E = hypergeometric error
frqYr<-simAlleleFreq$Year
frqCat<-simAlleleFreq$Category

#true error
simAlleleFreq[frqYr==1998 & frqCat=='errT',c(3:(nloci+2))]<-
	simAlleleFreq[frqYr==1998 & frqCat=='xt',c(3:(nloci+2))]-
	simAlleleFreq[frqYr==1998 & frqCat=='pt',c(3:(nloci+2))]

#hypergeometric
simAlleleFreq[frqYr==1998 & frqCat=='Et',c(3:(nloci+2))]<-
	(simAlleleFreq[frqYr==1998 & frqCat=='xt',c(3:(nloci+2))]*
	(1-simAlleleFreq[frqYr==1998 & frqCat=='xt',c(3:(nloci+2))])/
	(countsAll[countsAll$Year==1998,'genotyped']-1))*
	((countsAll[countsAll$Year==1998,'total']-
	countsAll[countsAll$Year==1998,'genotyped'])/
	(countsAll[countsAll$Year==1998,'total']-1))

for(year in c(1999:2013))
{
	#true error
	simAlleleFreq[frqYr==year & frqCat=='errT',c(3:(nloci+2))]<-
		simAlleleFreq[frqYr==year & frqCat=='xt',c(3:(nloci+2))]-
		simAlleleFreq[frqYr==year & frqCat=='pt',c(3:(nloci+2))]	
		
	simAlleleFreq[frqYr==year & frqCat=='errS',c(3:(nloci+2))]<-
		simAlleleFreq[frqYr==year & frqCat=='xs',c(3:(nloci+2))]-
		simAlleleFreq[frqYr==year & frqCat=='ps',c(3:(nloci+2))]
	
	simAlleleFreq[frqYr==year & frqCat=='errI',c(3:(nloci+2))]<-
		simAlleleFreq[frqYr==year & frqCat=='xi',c(3:(nloci+2))]-
		simAlleleFreq[frqYr==year & frqCat=='pi',c(3:(nloci+2))]	
		
	simAlleleFreq[frqYr==year & frqCat=='errB',c(3:(nloci+2))]<-
		simAlleleFreq[frqYr==year & frqCat=='xb',c(3:(nloci+2))]-
		simAlleleFreq[frqYr==year & frqCat=='pb',c(3:(nloci+2))]
		
	simAlleleFreq[frqYr==year & frqCat=='errM',c(3:(nloci+2))]<-
		simAlleleFreq[frqYr==year & frqCat=='xm',c(3:(nloci+2))]-
		simAlleleFreq[frqYr==year & frqCat=='pm',c(3:(nloci+2))]
		
	simAlleleFreq[frqYr==year & frqCat=='errF',c(3:(nloci+2))]<-
		simAlleleFreq[frqYr==year & frqCat=='xf',c(3:(nloci+2))]-
		simAlleleFreq[frqYr==year & frqCat=='pf',c(3:(nloci+2))]	

	#hypergeometric error
	simAlleleFreq[frqYr==year & frqCat=='Et',c(3:(nloci+2))]<-
		(simAlleleFreq[frqYr==year & frqCat=='xt',c(3:(nloci+2))]*
		(1-simAlleleFreq[frqYr==year & frqCat=='xt',c(3:(nloci+2))])/
		(countsAll[countsAll$Year==year,'genotyped']-1))*
		((countsAll[countsAll$Year==year,'total']-
		countsAll[countsAll$Year==year,'genotyped'])/
		(countsAll[countsAll$Year==year,'total']-1))
		
	simAlleleFreq[frqYr==year & frqCat=='Es',c(3:(nloci+2))]<-
		(simAlleleFreq[frqYr==year & frqCat=='xs',c(3:(nloci+2))]*
		(1-simAlleleFreq[frqYr==year & frqCat=='xs',c(3:(nloci+2))])/
		(counts[counts$Year==year & counts$Category=='survivor','genotyped']-1))*
		((counts[counts$Year==year & counts$Category=='survivor','total']-
		counts[counts$Year==year & counts$Category=='survivor','genotyped'])/
		(counts[counts$Year==year & counts$Category=='survivor','total']-1))
		
	simAlleleFreq[frqYr==year & frqCat=='Ei',c(3:(nloci+2))]<-
		(simAlleleFreq[frqYr==year & frqCat=='xi',c(3:(nloci+2))]*
		(1-simAlleleFreq[frqYr==year & frqCat=='xi',c(3:(nloci+2))])/
		(counts[counts$Year==year & counts$Category=='immigrant','genotyped']-1))*
		((counts[counts$Year==year & counts$Category=='immigrant','total']-
		counts[counts$Year==year & counts$Category=='immigrant','genotyped'])/
		(counts[counts$Year==year & counts$Category=='immigrant','total']-1))
		
	simAlleleFreq[frqYr==year & frqCat=='Eb',c(3:(nloci+2))]<-
		(simAlleleFreq[frqYr==year & frqCat=='xb',c(3:(nloci+2))]*
		(1-simAlleleFreq[frqYr==year & frqCat=='xb',c(3:(nloci+2))])/
		(counts[counts$Year==year & counts$Category=='nestling','genotyped']-1))*
		((counts[counts$Year==year & counts$Category=='nestling','total']-
		counts[counts$Year==year & counts$Category=='nestling','genotyped'])/
		(counts[counts$Year==year & counts$Category=='nestling','total']-1))

	#allele freq differences
	#true (population) freq
	simAlleleFreq[frqYr==year & frqCat=='pt1-pt',c(3:(nloci+2))]<-
		simAlleleFreq[frqYr==year & frqCat=='pt',c(3:(nloci+2))]-
		simAlleleFreq[frqYr==(year-1) & frqCat=='pt',c(3:(nloci+2))]

	simAlleleFreq[frqYr==year & frqCat=='ps-pt',c(3:(nloci+2))]<-
		simAlleleFreq[frqYr==year & frqCat=='ps',c(3:(nloci+2))]-
		simAlleleFreq[frqYr==(year-1) & frqCat=='pt',c(3:(nloci+2))]
		
	simAlleleFreq[frqYr==year & frqCat=='pi-pt',c(3:(nloci+2))]<-
		simAlleleFreq[frqYr==year & frqCat=='pi',c(3:(nloci+2))]-
		simAlleleFreq[frqYr==(year-1) & frqCat=='pt',c(3:(nloci+2))]
		
	simAlleleFreq[frqYr==year & frqCat=='pb-pt',c(3:(nloci+2))]<-
		simAlleleFreq[frqYr==year & frqCat=='pb',c(3:(nloci+2))]-
		simAlleleFreq[frqYr==(year-1) & frqCat=='pt',c(3:(nloci+2))]

	#sample freq
	simAlleleFreq[frqYr==year & frqCat=='xs-xt',c(3:(nloci+2))]<-
		simAlleleFreq[frqYr==year & frqCat=='xs',c(3:(nloci+2))]-
		simAlleleFreq[frqYr==(year-1) & frqCat=='xt',c(3:(nloci+2))]
		
	simAlleleFreq[frqYr==year & frqCat=='xi-xt',c(3:(nloci+2))]<-
		simAlleleFreq[frqYr==year & frqCat=='xi',c(3:(nloci+2))]-
		simAlleleFreq[frqYr==(year-1) & frqCat=='xt',c(3:(nloci+2))]
		
	simAlleleFreq[frqYr==year & frqCat=='xb-xt',c(3:(nloci+2))]<-
		simAlleleFreq[frqYr==year & frqCat=='xb',c(3:(nloci+2))]-
		simAlleleFreq[frqYr==(year-1) & frqCat=='xt',c(3:(nloci+2))]

	#true errors
	simAlleleFreq[frqYr==year & frqCat=='errS-errT',c(3:(nloci+2))]<-
		simAlleleFreq[frqYr==year & frqCat=='errS',c(3:(nloci+2))]-
		simAlleleFreq[frqYr==(year-1) & frqCat=='errT',c(3:(nloci+2))]
		
	simAlleleFreq[frqYr==year & frqCat=='errI-errT',c(3:(nloci+2))]<-
		simAlleleFreq[frqYr==year & frqCat=='errI',c(3:(nloci+2))]-
		simAlleleFreq[frqYr==(year-1) & frqCat=='errT',c(3:(nloci+2))]
		
	simAlleleFreq[frqYr==year & frqCat=='errB-errT',c(3:(nloci+2))]<-
		simAlleleFreq[frqYr==year & frqCat=='errB',c(3:(nloci+2))]-
		simAlleleFreq[frqYr==(year-1) & frqCat=='errT',c(3:(nloci+2))]
	
	#Mendelian noise
	simAlleleFreq[frqYr==year & frqCat=='pfam',c(3:(nloci+2))]<-0.5*
		(simAlleleFreq[frqYr==year & frqCat=='pm',c(3:(nloci+2))]+
		simAlleleFreq[frqYr==year & frqCat=='pf',c(3:(nloci+2))])
		
	simAlleleFreq[frqYr==year & frqCat=='pmend',c(3:(nloci+2))]<-
		simAlleleFreq[frqYr==year & frqCat=='pb',c(3:(nloci+2))]-
		simAlleleFreq[frqYr==year & frqCat=='pfam',c(3:(nloci+2))]
		
	simAlleleFreq[frqYr==year & frqCat=='pfam-pt',c(3:(nloci+2))]<-
		simAlleleFreq[frqYr==year & frqCat=='pfam',c(3:(nloci+2))]-
		simAlleleFreq[frqYr==(year-1) & frqCat=='pt',c(3:(nloci+2))]

	simAlleleFreq[frqYr==year & frqCat=='xfam',c(3:(nloci+2))]<-0.5*
		(simAlleleFreq[frqYr==year & frqCat=='xm',c(3:(nloci+2))]+
		simAlleleFreq[frqYr==year & frqCat=='xf',c(3:(nloci+2))])
		
	simAlleleFreq[frqYr==year & frqCat=='xmend',c(3:(nloci+2))]<-
		simAlleleFreq[frqYr==year & frqCat=='xb',c(3:(nloci+2))]-
		simAlleleFreq[frqYr==year & frqCat=='xfam',c(3:(nloci+2))]
		
	simAlleleFreq[frqYr==year & frqCat=='xfam-xt',c(3:(nloci+2))]<-
		simAlleleFreq[frqYr==year & frqCat=='xfam',c(3:(nloci+2))]-
		simAlleleFreq[frqYr==(year-1) & frqCat=='xt',c(3:(nloci+2))]

	simAlleleFreq[frqYr==year & frqCat=='errFAM',c(3:(nloci+2))]<-0.5*
		(simAlleleFreq[frqYr==year & frqCat=='errM',c(3:(nloci+2))]+
		simAlleleFreq[frqYr==year & frqCat=='errF',c(3:(nloci+2))])
		
	simAlleleFreq[frqYr==year & frqCat=='errMEND',c(3:(nloci+2))]<-
		simAlleleFreq[frqYr==year & frqCat=='errB',c(3:(nloci+2))]-
		simAlleleFreq[frqYr==year & frqCat=='errFAM',c(3:(nloci+2))]
		
	simAlleleFreq[frqYr==year & frqCat=='errFAM-errT',c(3:(nloci+2))]<-
		simAlleleFreq[frqYr==year & frqCat=='errFAM',c(3:(nloci+2))]-
		simAlleleFreq[frqYr==(year-1) & frqCat=='errT',c(3:(nloci+2))]

}

save(simAlleleFreq,file=paste("simAlleleFreq_",today,".rdata",sep=''))

#calculate variances and covariances
simVar<-data.frame(Year=rep(c(1999:2013),each=34),Category=rep(c('pt1-pt','ps-pt','xs-xt',
	'errS-errT','pspterrSerrT','pi-pt','xi-xt','errI-errT','pipterrIerrT','pb-pt','xb-xt',
	'errB-errT','pbpterrBerrT','pspi','xsxi','xserrI','errSxi','errSerrI','pspb','xsxb',
	'xserrB','errSxb','errSerrB','pipb','xixb','xierrB','errIxb','errIerrB','pmend',
	'xmend','errMEND','pfam-pt','xfam-xt','errFAM-errT'),15),stringsAsFactors=FALSE)

bsYr<-simVar$Year
bsCat<-simVar$Category
for(year in c(1999:2013))
{
	simVar[bsYr==year & bsCat=='pt1-pt',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pt1-pt',c(3:nloci+2)])^2)


	simVar[bsYr==year & bsCat=='ps-pt',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='ps-pt',c(3:nloci+2)])^2)
		
	simVar[bsYr==year & bsCat=='xs-xt',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xs-xt',c(3:nloci+2)])^2)
	
	simVar[bsYr==year & bsCat=='errS-errT',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errS-errT',c(3:nloci+2)])^2)
	
	simVar[bsYr==year & bsCat=='pspterrSerrT',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='ps-pt',c(3:nloci+2)])*
		as.numeric(simAlleleFreq[frqYr==year & frqCat=='errS-errT',c(3:nloci+2)]))


	simVar[bsYr==year & bsCat=='pi-pt',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pi-pt',c(3:nloci+2)])^2)
		
	simVar[bsYr==year & bsCat=='xi-xt',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xi-xt',c(3:nloci+2)])^2)
	
	simVar[bsYr==year & bsCat=='errI-errT',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errI-errT',c(3:nloci+2)])^2)
	
	simVar[bsYr==year & bsCat=='pipterrIerrT',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pi-pt',c(3:nloci+2)])*
		as.numeric(simAlleleFreq[frqYr==year & frqCat=='errI-errT',c(3:nloci+2)]))


	simVar[bsYr==year & bsCat=='pb-pt',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pb-pt',c(3:nloci+2)])^2)
		
	simVar[bsYr==year & bsCat=='xb-xt',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xb-xt',c(3:nloci+2)])^2)
		
	simVar[bsYr==year & bsCat=='errB-errT',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errB-errT',c(3:nloci+2)])^2)
		
	simVar[bsYr==year & bsCat=='pbpterrBerrT',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pb-pt',c(3:nloci+2)])*
		as.numeric(simAlleleFreq[frqYr==year & frqCat=='errB-errT',c(3:nloci+2)]))

	
	simVar[bsYr==year & bsCat=='pspi',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='ps-pt',c(3:(nloci+2))])*
		as.numeric(simAlleleFreq[frqYr==year & frqCat=='pi-pt',c(3:(nloci+2))]))

	simVar[bsYr==year & bsCat=='xsxi',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xs-xt',c(3:(nloci+2))])*
		as.numeric(simAlleleFreq[frqYr==year & frqCat=='xi-xt',c(3:(nloci+2))]))

	simVar[bsYr==year & bsCat=='xserrI',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xs-xt',c(3:(nloci+2))])*
		as.numeric(simAlleleFreq[frqYr==year & frqCat=='errI-errT',c(3:(nloci+2))]))

	simVar[bsYr==year & bsCat=='errSxi',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errS-errT',c(3:(nloci+2))])*
		as.numeric(simAlleleFreq[frqYr==year & frqCat=='xi-xt',c(3:(nloci+2))]))

	simVar[bsYr==year & bsCat=='errSerrI',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errS-errT',c(3:(nloci+2))])*
		as.numeric(simAlleleFreq[frqYr==year & frqCat=='errI-errT',c(3:(nloci+2))]))


	simVar[bsYr==year & bsCat=='pspb',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='ps-pt',c(3:(nloci+2))])*
		as.numeric(simAlleleFreq[frqYr==year & frqCat=='pb-pt',c(3:(nloci+2))]))

	simVar[bsYr==year & bsCat=='xsxb',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xs-xt',c(3:(nloci+2))])*
		as.numeric(simAlleleFreq[frqYr==year & frqCat=='xb-xt',c(3:(nloci+2))]))

	simVar[bsYr==year & bsCat=='xserrB',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xs-xt',c(3:(nloci+2))])*
		as.numeric(simAlleleFreq[frqYr==year & frqCat=='errB-errT',c(3:(nloci+2))]))

	simVar[bsYr==year & bsCat=='errSxb',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errS-errT',c(3:(nloci+2))])*
		as.numeric(simAlleleFreq[frqYr==year & frqCat=='xb-xt',c(3:(nloci+2))]))

	simVar[bsYr==year & bsCat=='errSerrB',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errS-errT',c(3:(nloci+2))])*
		as.numeric(simAlleleFreq[frqYr==year & frqCat=='errB-errT',c(3:(nloci+2))]))


	simVar[bsYr==year & bsCat=='pipb',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pi-pt',c(3:(nloci+2))])*
		as.numeric(simAlleleFreq[frqYr==year & frqCat=='pb-pt',c(3:(nloci+2))]))

	simVar[bsYr==year & bsCat=='xixb',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xi-xt',c(3:(nloci+2))])*
		as.numeric(simAlleleFreq[frqYr==year & frqCat=='xb-xt',c(3:(nloci+2))]))

	simVar[bsYr==year & bsCat=='xierrB',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xi-xt',c(3:(nloci+2))])*
		as.numeric(simAlleleFreq[frqYr==year & frqCat=='errB-errT',c(3:(nloci+2))]))

	simVar[bsYr==year & bsCat=='errIxb',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errI-errT',c(3:(nloci+2))])*
		as.numeric(simAlleleFreq[frqYr==year & frqCat=='xb-xt',c(3:(nloci+2))]))

	simVar[bsYr==year & bsCat=='errIerrB',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errI-errT',c(3:(nloci+2))])*
		as.numeric(simAlleleFreq[frqYr==year & frqCat=='errB-errT',c(3:(nloci+2))]))


	simVar[bsYr==year & bsCat=='pmend',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pmend',c(3:(nloci+2))])^2)

	simVar[bsYr==year & bsCat=='xmend',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xmend',c(3:(nloci+2))])^2)

	simVar[bsYr==year & bsCat=='errMEND',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errMEND',c(3:(nloci+2))])^2)


	simVar[bsYr==year & bsCat=='pfam-pt',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='pfam-pt',c(3:(nloci+2))])^2)

	simVar[bsYr==year & bsCat=='xfam-xt',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='xfam-xt',c(3:(nloci+2))])^2)

	simVar[bsYr==year & bsCat=='errFAM-errT',3]<-
		mean(as.numeric(simAlleleFreq[frqYr==year & frqCat=='errFAM-errT',c(3:(nloci+2))])^2)
}

#save output
save(simVar,file=paste("simVar_",today,".rdata",sep=''))
