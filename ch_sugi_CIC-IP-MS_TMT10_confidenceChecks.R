# TODO: Add comment
# 
# Author: cshughes
###############################################################################
#required libraries
library(limma)



###INPUT FILES
###this code will read in the input files specified by the user
###############################################################################
#read specified input files
setwd(dir="/Users/cshughes/Documents/projects/sugi/cicIP/")
#set PSMs.txt file
psm<-read.table("./ch_19June2015_CIC-IP_TMT10_Nuc_PSMs.txt", header=TRUE, sep='\t')
#set Proteins.txt file
pro<-read.table("./ch_19June2015_CIC-IP_TMT10_Nuc_Proteins.txt", header=TRUE, sep='\t')
#set working directory to an output directory (can be the same as above)
setwd(dir="/Users/cshughes/Documents/projects/sugi/cicIP/Routput/")


###PROCESS RAW DATA
###this code will process the raw PSMs to remove contaminants, redundant hits, and re-annotate the hits
###############################################################################
#assign the names of your 10plex sample set
plex = c('HEKa1','HEKa2','HEKa3','E2HF12a1','E2HF12a2','E2HF12a3','E4HH8a1','E4HH8a2','E4HH8a3','Neg')
#processing function that processes and annotates the peptide list
processPSM <- function(psmFile, proteinFile, sampleNames, ... ){	
	#reshape the PSM to remove non-informative columns
	pepCols = c('Annotated.Sequence','Modifications','Number.of.Protein.Groups','Master.Protein.Accessions','Quan.Info','X126','X127N','X127C','X128N','X128C','X129N','X129C','X130N','X130C','X131')
	pep<-psmFile[,pepCols]
	#grab more complete annotation information from the Proteins file
	proCols = c('Accession','Description','MW.in.kDa')
	pro<-proteinFile[,proCols]
	message('Raw number of PSMs')
	message(nrow(pep))
	#remove non-unique peptides or those without quan values
	pep<-subset(pep, Number.of.Protein.Groups==1)
	message('Unique Peptides Only')
	message(nrow(pep))
	#parse the protein accession to keep only the top hit if multiple groups assigned
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
	message('removing contaminants')
	pep.m = subset(pep.m, !grepl('Keratin',pep.m$Descriptions)) #this removes keratins...if you are interested in keratin, remove this statement!
	pep.m = subset(pep.m, !grepl('sp',pep.m$Accession, ignore.case=FALSE)) #this removes hits present in the common contaminants database
	pep.m = subset(pep.m, !grepl('Ig',pep.m$Descriptions, ignore.case=FALSE)) #this removes antibody hits
	pep.m = subset(pep.m, !grepl('ribosomal',pep.m$Descriptions, ignore.case=TRUE)) #this removes ribosomal hits
	message(nrow(pep.m))
	#reshape the data frame again to remove redundant columns we have already processed
	pCols = c('Accession','Gene','Descriptions','Sequence','MW.in.kDa','X126','X127N','X127C','X128N','X128C','X129N','X129C','X130N','X130C','X131')
	pep.r = pep.m[,pCols]
	#make a PSM counter
	pep.r$psms = 1
	#aggregate the PSMs into peptides...collapses instances where there are multiple PSMs for a single peptide sequence into a single entry
	message('collapsing redundant PSMs')
	pep.a1 = aggregate(cbind(MW.in.kDa,X126,X127N,X127C,X128N,X128C,X129N,X129C,X130N,X130C,X131)~Accession+Gene+Descriptions+Sequence,data=pep.r,median,na.rm=TRUE,na.action=na.pass)
	pep.a2 = aggregate(psms~Accession+Gene+Descriptions+Sequence,data=pep.r,sum,na.rm=TRUE,na.action=na.pass)
	colnames(pep.a1)[6:15] = sampleNames
	pep.a = merge(pep.a2,pep.a1,by=c('Accession','Gene','Descriptions','Sequence'),sort=FALSE)
	message(nrow(pep.a))
	#round all expression values to 2 places
	pep.a[,7:16] = round(pep.a[,7:16],2)
	#replaces zeroes with NA
	pep.a[,7:16][pep.a[,7:16]==0]<-NA
	#replaces NaN with NA
	pep.a[,7:16][is.na(pep.a[,7:16])]<-NA
	#remove any rows that are completely NA in the expression columns
	message('removing peptides with only NA values')
	pep.s1 = subset(pep.a, rowSums(is.na(pep.a[,7:16]))<9)
	message(nrow(pep.s1))
	#filter rows where the expression values are too low to get confident values
	message('filtering peptides by signal-to-noise')
	pep.s2 = subset(pep.s1, rowMeans(pep.s1[,7:16],na.rm=TRUE)>4)
	message(nrow(pep.s2))
	#output the data
	return(pep.s2)	
}
#run function
pepSet<-processPSM(psm,pro,plex)


###CALCULATE ENRICHMENT SCORE
###this code will calculate the enrichment score for your peptide hits
###you will need to change the column numbers around based on your experimental setup
###############################################################################
#what column are your IP samples in?
ipCols = 7:9
#what column are your negative samples in?
#negCols = 16 - should ignore the negative here...it is not good.
#prefilter the data
pepSet = subset(pepSet, rowSums(is.na(pepSet[,ipCols]))<2)
#need to change the column numbers to the numbers for your negative
#pepSet$Neg = pepSet[,negCols]
#set NA values in your negative to a base value of 1 for calculations
#pepSet$Neg = ifelse(is.na(pepSet$Neg),1,pepSet$Neg)
#calculate the signal in the IgG relative to the target protein signals
#pepSet$Neg = pepSet$Neg / rowMeans(pepSet[,ipCols],na.rm=TRUE)
#normalize the IP samples only (not the negative) using a median shift (limma package)
pepSet[,ipCols] = as.data.frame(normalizeMedianValues(as.matrix(pepSet[,ipCols])))
#calculate the signal score
eSignal = apply(pepSet[,ipCols],2,function(x) x * pepSet$psms)
#calculate the penalty score
pepSet$penalty = rowSums(!is.na(pepSet[,ipCols])) * sqrt(apply(pepSet[,ipCols],1,function(x) var(x,na.rm=TRUE))) #* pepSet$Neg
#calculate the score relative to the negative
eScore = apply(eSignal,2,function(x) x / pepSet$penalty)
#rebind the data
eSet = cbind(pepSet[,c(1:5,6)],eScore)


###CALCULATE ENRICHMENT SCORE
###roll into proteins
#PRMT5, RUVBL1, SMARCB1, ARID1A, RBBP4
###############################################################################
#add a peptide counter for the proteins
eSet$pepNum = 1
#aggregate with a sum for peptide number
pro1 = aggregate(cbind(psms,pepNum)~Accession+Gene+Descriptions, data=eSet, sum, na.rm=TRUE, na.action=na.pass)
#aggregate with a mean for the MW and enrichment score
pro2 = aggregate(cbind(MW.in.kDa,HEKa1,HEKa2,HEKa3)~Accession+Gene+Descriptions, data=eSet, mean, na.rm=TRUE, na.action=na.pass)
#merge the two aggregated sets
proSet = merge(pro1,pro2,by=c('Accession','Gene','Descriptions'))
proSet$mean_eScore = rowMeans(proSet[,7:9],na.rm=TRUE)
#sort by decreasing enrichment score
proSet = proSet[order(-proSet$mean_eScore),]
#output the result to a text file
write.table(proSet,'ch_test5.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)



