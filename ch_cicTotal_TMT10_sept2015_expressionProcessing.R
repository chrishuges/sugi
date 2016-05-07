# TODO: play around with the expression data
# 
# Author: cshughes
###############################################################################
library(limma)
library(PECA)
library(gplots)
library(vsn)
###############################################################################
#HEK, F12, D10 and A9 (getting rid of H8 and D12)
pep = readRDS('ch_CIC_TMT10_sept2015_processedPeptides.rds')

##get rid of the unwanted samples and rearrange the columns
pep = pep[,c(1:6,11:12,8:9,13:14)]

#normalize using a median transformation
normMed = normalizeMedianValues(as.matrix(pep[,5:12]))
normMed = as.data.frame(log2(normMed))

#or normalize using VSN
vsnMed = as.data.frame(justvsn(as.matrix(pep[,5:12])))
pdf('ch_test.pdf')
meanSdPlot(justvsn(as.matrix(pep[,5:12]))) #check the quality of the normalization
dev.off()
#################normalizeMedianValues Base Code
#x is your input matrix
#cmed <- log(apply(x, 2, median, na.rm = TRUE))
#cmed <- exp(cmed - mean(cmed))
#t(t(x)/cmed)
##################

######perform some statistical analysis of expression
#get the samples I want to compare
samplenames1 = colnames(pep[,5:8])
samplenames2 = colnames(pep[,9:12])
#make the design matrix
design = cbind(G1 = 1, G1vsG2 = c(rep(1, length(samplenames1)),rep(0, length(samplenames2))))
#calculate the gene expression value
#######just need to change the input matrix here to test different things (e.g. normMed)...values should be log scale here and normalized
fit = lmFit(as.matrix(vsnMed), design)
fit = eBayes(fit)
probeFC = fit$coefficients[, 2]
t <- fit$t[, 2]
df.total <- fit$df.residual[1] + fit$df.prior
#get some annotation details for aggregation
pGenes = paste(pep$Gene,'_',pep$Accession, sep='')
#perform the aggregation
gene.n <- tapply(t, pGenes, function(x) sum(!is.na(x)))
#calculate median fold change per protein
geneFC <- tapply(probeFC, pGenes, median, na.rm = TRUE)
#calculate median t-statistic per protein
t <- tapply(t, pGenes, median, na.rm = TRUE)
#calculate the p-values
gene.p <- 2 * pt(abs(t), df = df.total, lower.tail = FALSE)
#adjust the p-values
gene.p.fdr <- p.adjust(gene.p, method = "fdr")
#bind the data into an output set
result <- data.frame(cbind(logFC = geneFC, t = t, score = gene.p, p.fdr = gene.p.fdr))
result$Accession = row.names(result)
row.names(result) = NULL



#the VSN fit seems to give a bit better results, but in general the treatment of 'replicates' is hurting the analysis
#maybe combine the values across passages for the same cell line, or alternating cell lines and redo the stats? 
#the correlation matrix is quite good...hmm.



pdf('ch_cicTest_HeatMap.pdf')
#make the plot labels and boundaries
xLabels<- names(vsnMed)
mybreaks = seq(0,0.9,by=0.05)
ColSideColors = c(rep(brewer.pal(6,'Accent')[1],4),rep(brewer.pal(6,'Accent')[6],4))
#make the correlation heatmap
heatmap.2(
		cor(vsnMed, use='pairwise.complete.obs', method='pearson'),
		col= colorRampPalette(brewer.pal(9,"YlOrRd"))(length(mybreaks)-1),
		symkey=TRUE,
		Rowv=FALSE,
		Colv=FALSE,
		dendrogram="none",
		breaks=mybreaks,
		cellnote = round(cor(vsnMed, use='pairwise.complete.obs', method='pearson'),2),
		labRow = xLabels,
		labCol = xLabels,
		notecol = 'black',
		notecex = 1.2,
		colsep = 1:8,
		rowsep = 1:8,
		sepwidth = c(0.03, 0.03),
		sepcolor = 'white',
		ColSideColors=ColSideColors,
		## labels
		main='OvC Subtype Correlation',
		## color key
		key = TRUE,
		keysize = 1,
		density.info = "none",
		scale = "none",
		trace = "none",
		mar=c(8,8),
		cexRow=1.2,
		cexCol=1.2
)
dev.off()







