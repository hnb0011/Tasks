#15 Feb 2018
#Nancy Chen
#script to analyze gene dropping results
#core demo nestlings
#multiple tests of selection

today<-format(Sys.Date(),format="%d%b%Y")
library(plyr)

snplist<-read.table('SNPlist.txt',header=TRUE,stringsAsFactors=FALSE)

#net change in allele freq between 1999-2013
pval_1999to2013<-data.frame(snp=snplist$SNP,stringsAsFactors=FALSE)

#look at obs vs sim change in allele freq between adjacent years over time
pval_change<-data.frame(snp=rep(snplist$SNP,each=14),year=rep(c(1999:2012),length(snplist$SNP)),stringsAsFactors=FALSE)

for (i in snplist$SNP)
{
	obs<-read.table(file=paste('batch_',i,'.',j,'.drop.data.txt',sep=''),header=TRUE)
	sim<-read.table(file=paste('batch_',i,'.',j,'.drop.sim.txt',sep=''),header=TRUE)
	simfreq<-mapply('/',sim[,11:25],obs[obs$allele==2 & obs$cohort_year>1998,'all_alleles_count'])

	#test for net selection between 1999-2013
	obsDelta1999<-obs[obs$allele==2 & obs$cohort_year==2013,'frequency_of_allele']-obs[obs$allele==2 & obs$cohort_year==1999,'frequency_of_allele']
	simDelta1999<-simfreq[,15]-simfreq[,1]
	pval_1999to2013[pval_1999to2013$snp==i,'obsChange']<-obsDelta1999
	pval_1999to2013[pval_1999to2013$snp==i,'simChange']<-median(simDelta1999)
	simdiff1999<-simDelta1999 - median(simDelta1999)
	obsdiff1999<-obsDelta1999 - median(simDelta1999)

	if (obsdiff1999 < 0)
	{
		pval_1999to2013[pval_1999to2013$snp==i,'dir']<-'-'
		pval_1999to2013[pval_1999to2013$snp==i,'pval']<-sum(simdiff1999<obsdiff1999)/1000000 + sum(simdiff1999==obsdiff1999)/2000000
	} else 
	{
		pval_1999to2013[pval_1999to2013$snp==i,'dir']<-'+'
		pval_1999to2013[pval_1999to2013$snp==i,'pval']<-sum(simdiff1999>obsdiff1999)/1000000 + sum(simdiff1999==obsdiff1999)/2000000
	}

	#test for selection in adjacent years
	for(x in c(1:14)) 
	{
		yr<-1998 + x
		obsDelta<-obs[obs$allele==2 & obs$cohort_year==(yr+1),'frequency_of_allele']-obs[obs$allele==2 & obs$cohort_year==yr,'frequency_of_allele']
		simDelta<-simfreq[,x+1]-simfreq[,x]
		pval_change[pval_change$snp==i & pval_change$year==yr,'obsChange']<-obsDelta
		pval_change[pval_change$snp==i & pval_change$year==yr,'simChange']<-median(simDelta)
		simdiff<-simDelta - median(simDelta)
		obsdiff<-obsDelta - median(simDelta)

		if (obsdiff < 0)
		{
			pval_change[pval_change$snp==i & pval_change$year==yr,'dir']<-'-'
			pval_change[pval_change$snp==i & pval_change$year==yr,'pval_change']<-sum(simdiff<obsdiff)/1000000 + sum(simdiff==obsdiff)/2000000
		} 
		else 
		{
			pval_change[pval_change$snp==i & pval_change$year==yr,'dir']<-'+'
			pval_change[pval_change$snp==i & pval_change$year==yr,'pval_change']<-sum(simdiff>obsdiff)/1000000 + sum(simdiff==obsdiff)/2000000
		}
	}

}

save(pval_change,pval_1999to2013,file=paste("genedropSelResults_",today,".rdata",sep=''))

