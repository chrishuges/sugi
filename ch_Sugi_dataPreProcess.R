# TODO: Add comment
# 
# Author: cshughes
###############################################################################
#load in files
setwd(dir="/Users/cshughes/Documents/projects/sugi/Rlessons/")
psm<-read.table("./SC_17Dec2015_Sugi_C__PSMs.txt", header=TRUE, sep='\t')
pro<-read.table("./SC_17Dec2015_Sugi_C__Proteins.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/sugi/Rlessons/Routput/")


###############################################################################
#process the raw data
processPSM <- function(psmFile, proteinFile, replicate, ... ){	
	require(limma)
	#initial processing to remove non-informative columns and rows from contaminant proteins or those with many missing values
	pep<-psmFile[,c(5,6,9:11,36:46)]
	pro<-proteinFile[,c(3:4,12)]
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
	pep.m = subset(pep.m, !grepl('sp',pep.m$Accession, ignore.case=FALSE))
	#filter the data columns
	print('removing contaminants')
	print(dim(pep.m))
	pep.r = pep.m[,c(1,20,22,21,19,8:17)]
	#aggregate the PSMs into peptides
	print('aggregating peptides')
	pep.a = aggregate(cbind(MW.in.kDa,X126,X127N,X127C,X128N,X128C,X129N,X129C,X130N)~Accession+Gene+Descriptions+Sequence,data=pep.r,median,na.action=na.pass,na.rm=TRUE)
	colnames(pep.a) = c('Accession','Gene','Descriptions','Sequence','MW','wt1','m1a','m2a','m3a','wt2','m1b','m2b','m3b')
	print(dim(pep.a))
	#replace NaN with NA
	pep.a[,6:13] = round(pep.a[,6:13],2)
	pep.a[,6:13][is.na(pep.a[,6:13])]<-NA
	#filter based on S/N
	pep.f1 = subset(pep.a, rowMeans(pep.a[,6:13], na.rm=TRUE)>5)
	#filter based on NA
	pep.f2 = subset(pep.f1, rowSums(is.na(pep.f1[,6:13]))<7)
	#add replicate counter
	pep.f2$Rep = replicate
	#output the data
	return(pep.f2)	
}
#run the function
psm.p<-processPSM(psm, pro, 'a1')


#output the data objects
saveRDS(psm.p,'ch_Sugi_processedPeptides.rds')
write.table(psm.p,'ch_Sugi_processedPeptides.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)

