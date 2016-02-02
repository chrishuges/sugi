###############################################################################
#starting with R
###############################################################################

###install a required package
#this tells R where to look for packages
source('http://bioconductor.org/biocLite.R')
#this tells R the package you want to install
biocLite('limma')
#this loads the package in R, you need to load a package every time you freshly start R
library(limma)


###load in your data files to be processed
#this changes what your home directory is
setwd(dir="/Users/cshughes/Documents/projects/sugi/Rlessons/Routput/")
#this tells R to read the file, saying the first row is the column names, and that it is tab-delimited 
pep<-read.table("./ch_Sugi_processedPeptides.txt", header=TRUE, sep='\t')
#or
pep<-readRDS("./ch_Sugi_processedPeptides.rds")
#it is always better to use an RDS file because a lot of times the gene names can get mixed up in a text file, or the import won't work properly
#if you have decided you don't like being cool and don't want to use RDS files anymore, you can just read it and export it to text
write.table(pep,'ch_Sugi_processedPeptides.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)


#look at the top part of the data, just to see that it is correct
head(pep)
#look how big the matrix is
dim(pep) #52071 8

#extract just the expression values
#first check which column numbers correspond to your data of interest
names(pep)
#looks like 6 through 13 are the columns we want
eset = as.matrix(pep[,6:13])
#or
eset = as.matrix(pep[,c(6:13)])
#these are technically the same command...but, in the first you can only use a set of columns that are in series (e.g. 6:13)
#in the second version, you can use any columns e.g. c(1,6:13)
#you can check it again to see if it is what you want
head(eset)
#you can also look at the characteristics of the total data set
summary(eset)

#now lets play with the data...
#lets log transform first
eset.log = log2(eset)
#median normalize using a premade function in Limma
eNorm1 = normalizeMedianValues(eset.log)
#you can also do this manually using the apply() functions

#we can also normalize the data using a premade function that works on signal intensity (e.g. VSN, Limma)
#luckily VSN is built into Limma already, so we don't need to install the package, but lets do it anyway
biocLite(vsn)
library(vsn)
#run the limma version
eNorm2 = normalizeVSN(eset)
#check the normalization
meanSdPlot(eNorm2)
#run the vsn version
eNorm3 = justvsn(eset)
meanSdPlot(eNorm3)
#the two plots should look identical...only need to use one version of this in the future
#lets save the data for use later
#as an RDS
saveRDS(eNorm3,'ch_Sugi_normalizedPeptides.rds')
#as a text file you can open in excel
write.table(eNorm3,'ch_Sugi_normalizedPeptides.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)



#lets do some stats of the conditions vs. controls using limma
#first we need to make a design matrix that tells limma what our samples are, 1 means yes, 0 means no based on columns
design <- cbind(wt=c(1,0,0,0,1,0,0,0),mut=c(0,1,1,1,0,1,1,1))
#calculate the linear models for each probe based on the sample groups we just made in design
fit <- lmFit(eNorm3,design)
#sometimes a warning message will appear: 'Partial NA coefficients for X probe(s)...this just means there wasn't enough data to build the model for 25 probes
#make a contrast matrix that tells limma how you want values calculated for your specific comparison
cont.matrix <- makeContrasts(wtvsmut=wt-mut, levels=design)
#this will recalculate the linear fits for the specific comparison of interest, for example if you had more than 2 groups in the original fit
fit2 <- contrasts.fit(fit, cont.matrix)
#now we calculate the statistical values, based on modeling of a population variance calculated by the eBayes function
eFit <- eBayes(fit2)
#the eFit object now has all of your data of interest
names(eFit)
#we can also correct the p-values for multiple testing
p.adj = p.adjust(eFit$p.value, method="BH")
#importantly, it is still the same size as the original data
dim(eFit) #52071 1
#the data right now is just expression and stats details, but we want to recombine it with the annotation details
ePep = cbind(pep[,c(1:5)], 'logFC' = eFit$coef, 'pVal' = eFit$p.value, 'adjPVal' = p.adj)
#cbind combines things that are the same length in terms of number of rows
#we need to change the column names for ePep
colnames(ePep) = c('Accession','Gene','Descriptions','Sequence','MW','logFC','pVal','adjPVal')
#this is our peptide set with all the annotation

#what about if we want a protein set?
#first lets make a peptide counter
ePep$pepNum = 1
#now aggregate the set into proteins calculating the number of peptides per protein
pro1 = aggregate(pepNum~Accession+Gene+Descriptions,data=ePep,sum,na.action=na.pass,na.rm=TRUE)
#but that only has peptide numbers, we want the expression data too
pro2 = aggregate(cbind(logFC,pVal,adjPVal)~Accession+Gene+Descriptions,data=ePep,median,na.action=na.pass,na.rm=TRUE)
#now we have two separate aggregated sets that both have data we want, lets combine them
pro = merge(pro1,pro2,by=c('Accession','Gene','Descriptions'))
#take a look at the data
head(pro)


#lets save the data for use later
#as an RDS
saveRDS(pro,'ch_Sugi_processedProteins.rds')
#as a text file you can open in excel
write.table(pro,'ch_Sugi_processedProteins.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)



