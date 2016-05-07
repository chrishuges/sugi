# TODO: reprocessing data for Sugi
# 
# Author: cshughes
###############################################################################
library(lattice)
library(preprocessCore)
library(vsn)
library(Biobase)
###############################################################################
#load in files
setwd(dir="/Users/cshughes/Documents/projects/sugi/cicTotal/")
cpsm<-read.table("./ch_05Sept2015_Sugi_CICrep_HpH_All_PSMs.txt", header=TRUE, sep='\t')
cpro<-read.table("./ch_05Sept2015_Sugi_CICrep_HpH_All_Proteins.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/sugi/cicTotal/Routput/")

##################################################
#initial processing of the data to compile peptides and remove contaminants
##################################################
#dont split replicates yet or do any combination of them
processPSM <- function(psmFile, proteinFile, ... ){	
	require(limma)
	#initial processing to remove non-informative columns and rows from contaminant proteins or those with many missing values
	pep<-psmFile[,c(5,6,9:11,36:46)]
	pro<-proteinFile[,c(3:4)]
	print('Raw data frame dimensions')
	print(dim(pep))
	#remove non-unique peptides or those without quan values
	#pep<-subset(pep, Number.of.Proteins==1)
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
	#remove contaminant proteins
	print('removing contaminants')
	pep.m = subset(pep.m, !grepl('sp',pep.m$Accession, ignore.case=FALSE))
	#filter the data columns
	print(dim(pep.m))
	pep.r = pep.m[,c(1,19,21,20,8:17)]
	#aggregate the PSMs into peptides
	print('aggregating peptides')
	pep.a = aggregate(cbind(X126,X127N,X127C,X128N,X128C,X129N,X129C,X130N,X130C,X131)~Accession+Gene+Descriptions+Sequence,data=pep.r,median,na.action=na.pass,na.rm=TRUE)
	colnames(pep.a) = c('Accession','Gene','Descriptions','Sequence','HEKp11','F12p11','H8p11','A9p11','D10p11','D12p11','HEKp8','F12p8','A9p8','D10p8')
	print(dim(pep.a))
	#replace NaN with NA
	pep.a[,5:14] = round(pep.a[,5:14],2)
	pep.a[,5:14][is.na(pep.a[,5:14])]<-NA
	#filter based on S/N
	pep.f1 = subset(pep.a, rowMeans(pep.a[,5:14], na.rm=TRUE)>4)
	#filter based on NA
	pep.f2 = subset(pep.f1, rowSums(is.na(pep.f1[,5:14]))<9)
	#output the data
	return(pep.f2)	
}
#run the function
c.psm<-processPSM(cpsm, cpro)

#output the data objects
saveRDS(c.psm,'ch_CIC_TMT10_sept2015_processedPeptides.rds')



