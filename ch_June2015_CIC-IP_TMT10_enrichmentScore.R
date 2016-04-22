# TODO: Add comment
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
#############################################
##################################################
#data processing
##################################################
'HEKa1','HEKa2','HEKa3','E2HF12a1','E2HF12a2','E2HF12a3','E4HH8a1','E4HH8a2','E4HH8a3'
###data sets - CIC processed in Proteome Discoverer
setwd(dir="/Users/cshughes/Documents/projects/sugi/cicIP/")
cicpsm<-read.table("./ch_19June2015_CIC-IP_TMT10_Cyto_PSMs.txt", header=TRUE, sep='\t')
cicpro<-read.table("./ch_19June2015_CIC-IP_TMT10_Cyto_Proteins.txt", header=TRUE, sep='\t')
#cicpsm<-read.table("./ch_19June2015_CIC-IP_TMT10_Nuc_PSMs.txt", header=TRUE, sep='\t')
#cicpro<-read.table("./ch_19June2015_CIC-IP_TMT10_Nuc_Proteins.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/sugi/cicIP/Routput/")

##################################################
#identification processing
##################################################
#unfortunately they are different for some reason, so need two functions
processPSM <- function(psmFile, proteinFile, ... ){	
	require(limma)
	#initial processing to remove non-informative columns and rows from contaminant proteins or those with many missing values
	pep<-psmFile[,c(5,6,9,10,11,37:47)]
	pro<-proteinFile[,c(3:4,13)]
	print('Raw data frame dimensions')
	print(dim(pep))
	#remove non-unique peptides or those without quan values
	pep<-subset(pep, Number.of.Proteins==1)
	pep<-subset(pep, !grepl('NoQuanValues',pep$Quan.Info))
	print('Unique Peptides Only')
	print(dim(pep))
	#parse the protein accession
	pep$Accession = sapply(strsplit(as.character(pep$Master.Protein.Accessions), ';'),'[', 1)
	#merge the protein descriptions into the peptide file
	pep.m = merge(pep,pro,by='Accession')
	#get the gene name out
	pep.m$Gene<-sub(".*?GN=(.*?)( .*|$)", "\\1", pep.m$Description)
	#parse the peptide column for amino acids
	pep.m$Sequence = toupper(sub('.*?\\.(.*?)(\\..*|$)','\\1',pep.m$Annotated.Sequence))
	#filter information from Description
	pep.m$Descriptions = sub('(.*?)( OS=.*|$)','\\1',pep.m$Description)
	#remove specific proteins
	pep.m = subset(pep.m, !grepl('Keratin',pep.m$Descriptions))
	pep.m = subset(pep.m, !grepl('sp',pep.m$Accession, ignore.case=FALSE))
	pep.m = subset(pep.m, !grepl('ribosomal',pep.m$Descriptions))
	#filter the data columns
	print('removing contaminants')
	print(dim(pep.m))
	pep.r = pep.m[,c(1,20,22,21,19,8:17)]
	#aggregate the PSMs into peptides
	print('aggregating peptides')
	pep.a = aggregate(cbind(MW.in.kDa,X126,X127N,X127C,X128N,X128C,X129N,X129C,X130N,X130C,X131)~Accession+Gene+Descriptions+Sequence,data=pep.r,mean,na.action=na.pass,na.rm=TRUE)
	colnames(pep.a) = c('Accession','Gene','Descriptions','Sequence','MW','HEKa1','E2HF12a1','E4HH8a1','HEKa2','E2HF12a2','E4HH8a2','HEKa3','E2HF12a3','E4HH8a3','Neg')
	print(dim(pep.a))
	#replace NaN with NA
	pep.a[,6:15] = round(pep.a[,6:15],2)
	pep.a[,6:15][pep.a[,6:15]==0]<-NA
	pep.a[,6:15][is.na(pep.a[,6:15])]<-NA
	#output the data
	return(pep.a)	
}
#run function
cic.psm<-processPSM(cicpsm,cicpro)

##################################################
#expression processing
##################################################
#remove rows with too many missing values
cic.psm = subset(cic.psm, rowSums(is.na(cic.psm[,c(6:14)]))<6)
#data output
write.table(cic.psm,'ch_June2015_CIC-IP-MS_Cyto_peptideSet.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)
#set NA values to a baseline value
cic.psm$Neg = ifelse(is.na(cic.psm$Neg),2,cic.psm$Neg)
#make a new data fram with just annotation
xnorm = apply(cic.psm[,c(6:14)],2, function(x) x - cic.psm$Neg)
vnorm = scale(xnorm)
cic.q = cbind(cic.psm[,1:5],vnorm)

#check with a plot
pdf('ch_cic-IP_EnrichmentDensity.pdf')
plotDensities(cic.q[,6:14],
		legend = FALSE,
		main = 'Enrichment Score Density of Peptides per IP')
abline(v=0,col='red',lty=2,lwd=2)
dev.off()
#I think it looks ok here...the distributions look even


##################################################
#protein processing
##################################################
cic.q$pepNum = 1
cic.pnum = aggregate(MW~Accession+Gene+Descriptions,data=cic.q,mean,na.action=na.pass,na.rm=TRUE)
cic.pro1 = aggregate(cbind(pepNum,HEKa1,HEKa2,HEKa3,E2HF12a1,E2HF12a2,E2HF12a3,E4HH8a1,E4HH8a2,E4HH8a3)~Accession+Gene+Descriptions,data=cic.q,sum,na.action=na.pass,na.rm=TRUE)
cic.pro = merge(cic.pnum,cic.pro1,by=c('Accession','Gene','Descriptions'))
cic.pro$HEKfc = rowMeans(cic.pro[,c(6:8)],na.rm=TRUE)
cic.pro$E2HF12fc = rowMeans(cic.pro[,c(9:11)],na.rm=TRUE)
cic.pro$E4HH8fc = rowMeans(cic.pro[,c(12:14)],na.rm=TRUE)
#assign p-values
reps = cic.pro[,c(6,9,12)]
fit <- lmFit(reps)
fit <- eBayes(fit)
p.value <- fit$p.value
p.adj = p.adjust(p.value, method="BH")
cic.pro$HEK_pAdj <- p.adj
reps = cic.pro[,c(7,10,13)]
fit <- lmFit(reps)
fit <- eBayes(fit)
p.value <- fit$p.value
p.adj = p.adjust(p.value, method="BH")
cic.pro$E2HF12_pAdj <- p.adj
reps = cic.pro[,c(8,11,14)]
fit <- lmFit(reps)
fit <- eBayes(fit)
p.value <- fit$p.value
p.adj = p.adjust(p.value, method="BH")
cic.pro$E4HH8_pAdj <- p.adj

#data output
write.table(cic.pro,
		'ch_June2015_CIC-IP-MS_Cyto_proteinSet.txt',
		quote=FALSE,
		sep='\t',
		col.names=TRUE,
		row.names=FALSE)


##################################################
#data plotting
##################################################
ip = cic.pro[order(cic.pro$IP_FC),]
ipb = cic.pro[order(cic.pro$IPb_FC),]
ip.s = subset(ip, IP_FC>0)
ipb.s = subset(ipb, IPb_FC>0)
pro.s = subset(cic.pro, IP_FC>0 & IPb_FC>0)

pdf('ch_cic12-IP_293A_dataPlots.pdf')
plot(log2(ip.s$IP_FC),pch=19,cex = 1.75, ylab = 'log2(Enrichment Score)')
box(lwd=3)
plot(log2(ipb.s$IPb_FC),pch=19,cex = 1.75, ylab = 'log2(Enrichment Score)')
box(lwd=3)
heatscatter(log2(cic.pro$IP_FC),log2(cic.pro$IPb_FC),main=paste('spearman R = ',round(cor(log2(cic.pro$IP_FC),log2(cic.pro$IPb_FC),use='pairwise.complete.obs',method='spearman'),6),sep=''),pch=19,cex=1.5)
dev.off()



