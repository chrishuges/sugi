# TODO: processing CIC mass spec data from IP samples processed in Proteome Discoverer
# 
# Author: cshughes
###############################################################################
library(RColorBrewer)
library(lattice)
library(gplots)
library(affy)
library(preprocessCore)
library(limma)
library(cluster)
library(vsn)
library(Biobase)
###############################################################################
#################################################
#data preprocessing
#################################################
#import the files
setwd(dir="/Users/cshughes/Documents/projects/sugi/")
pep<-as.data.frame(read.table("./ch_May2015_CIC_peptides.txt", header=TRUE, sep="\t", quote=""))
setwd(dir="/Users/cshughes/Documents/projects/sugi/Routput/")
dim(pep) #8787 129
#basic filtering of the data
diMProcessing <- function(x, ...){
	require(vsn)
	require(Biobase)
	#basic processing of the data to remove contaminants
	pep.s<-subset(x, !grepl('\\+',x$Potential.contaminant))
	pep.s<-subset(pep.s, !grepl('\\+',pep.s$Reverse))
	pep.s<-subset(pep.s, !grepl('CON',pep.s$Leading.razor.protein))
	print('Raw Dataset')
	print(dim(pep.s))
	#strip groups
	pro.names<-as.character(pep.s$Proteins)
	pro.names.parse<-sapply(strsplit(pro.names, ';'), '[', 1)
	pep.s$Proteins<-pro.names.parse
	#strip isoforms
	pro.names<-as.character(pep.s$Proteins)
	pro.names.parse<-sapply(strsplit(pro.names, '-'), '[', 1)
	pep.s$Proteins<-pro.names.parse
	#strip gene names
	pro.names<-as.character(pep.s$Gene.names)
	pro.names.parse<-sapply(strsplit(pro.names, ';'), '[', 1)
	pep.s$Gene.names<-pro.names.parse	
	#subset the matrix to remove non-informative columns
	pep.ss = pep.s[,c(1,34,35,38,40,103:104,106:107,109:110,112:113,115:116,118:119)]
	#filter any peptides without quant data looking only at medium and heavy
	#first convert 0 to NA
	pep.ss[pep.ss==0]<-NA
	colnames(pep.ss) = c('Sequence','Proteins','Accession','Gene','Unique','cH8l.1','cH8h.1','cHEKl.1','cHEKh.1','nH8l.1','nH8h.1','nHEKl.1','nHEKh.1','nH8l.2','nH8h.2','nHEKl.2','nHEKh.2')
	#output the data
	return(pep.ss)
}
#apply the function
pep.s = diMProcessing(pep)

#################################################
#analysis of nuclear IP only
#################################################
#make a new data frame with just nuc
pep.n = pep.s[,c(1:4,10:17)]
#filter out samples that are not present in both replicates, per cell line
pep.sn = subset(pep.n, rowSums(is.na(pep.n[,c(5:6,9:10)]))<3 & rowSums(is.na(pep.n[,c(7:8,11:12)]))<3)
dim(pep.sn) #2943 12
#normalize the values
exprsFile1<-as.matrix(pep.sn[,5:12])
xnorm1.1<-justvsn(exprsFile1)
xnorm1<-as.data.frame(justvsn(exprsFile1))
colnames(xnorm1)<- names(pep.sn[,5:12])
pep.q<-cbind(pep.sn[,c(1:4)],xnorm1)
#test the fit
pdf('ch_CIC_IP_peptides_vsnFit.pdf')
meanSdPlot(xnorm1.1, ranks=TRUE)
dev.off()
#get the heavy to light ratios
pep.q$rH8.1 = pep.q$nH8l.1 - pep.q$nH8h.1 
pep.q$rH8.2 = pep.q$nH8l.2 - pep.q$nH8h.2 
pep.q$rHEK.1 = pep.q$nHEKl.1 - pep.q$nHEKh.1 
pep.q$rHEK.2 = pep.q$nHEKl.2 - pep.q$nHEKh.2 
#aggregate into proteins
pro = aggregate(cbind(nH8l.1,nH8h.1,nH8l.2,nH8h.2,nHEKl.1,nHEKh.1,nHEKl.2,nHEKh.2,rH8.1,rH8.2,rHEK.1,rHEK.2)~Accession+Gene, data=pep.q, mean, trim=0.2, na.action=na.pass, na.rm=TRUE)
dim(pro) #698 14
#replace NaN with NA
pro[pro=='NaN'] = NA
#assign p-values
reps = pro[,c(11:14)]
fit <- lmFit(reps)
fit <- eBayes(fit)
p.value <- fit$p.value
p.adj = p.adjust(p.value, method="BH")
pro$pAdj <- p.adj
pro$MeanFC <- rowMeans(pro[,c(11:14)], na.rm=TRUE)
pro$Var <- apply(as.matrix(pro[,c(11:14)]),1,var,na.rm=TRUE)

########################################################
#data output
########################################################
#write the data out
write.table(pro,
		'ch_may2015_CIC-IP_proteins_diM.txt',
		quote=FALSE,
		sep='\t',
		col.names=TRUE,
		row.names=FALSE)













