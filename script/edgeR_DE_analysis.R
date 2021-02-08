# Load packages, setup colours
library(edgeR)
library(gplots)
library(ggplot2)
library(dynamicTreeCut)
library(statmod)
install.packages("RColorBrewer")
library(RColorBrewer)

# Setup input (default is command line)
args <- commandArgs(trailingOnly = TRUE)
sub_analyse = paste(args[1])
FDR2use = as.numeric(paste(args[2]))

# example
sub_analyse <- 'Atab'
FDR2use  <- 0.05

datapath <- "~/Decaytrait_qtl/input/"
outpath <- paste("~/Decaytrait_qtl/output", sub_analyse, sep="")
dir.create(file.path(outpath))

annotation <- read.delim(file.path(datapath, "Atab_annotation.txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE) 

count <- read.table(file.path(datapath, paste(sub_analyse,'_count.txt', sep="")), header=T, row.names=1)
count <- round(count, digits=0)
design <- read.table(file.path(datapath, paste(sub_analyse,'_design.txt', sep="")), header=T)

model.formula <- as.formula("~0+group")
dmat <- model.matrix(model.formula,data=as.data.frame(design))
dgl <- DGEList(counts=count, group=design$group, genes=annotation)
dgl <- DGEList(counts=count, group=NULL, genes=annotation)
paste("all transcripts:", nrow(dgl))

dgl <- dgl[aveLogCPM(dgl) > 0,] # filter by average reads
#dgl <- dgl[rowSums(cpm(dgl)>1) >= 2,]

write(paste("dgl"), filter_file, append=T)
write(paste(nrow(dgl), "_", sep=""), filter_file, append=T)
write(summary(rowSums(cpm(dgl)/ncol(dgl))), filter_file, append=T, sep='\t', ncol=6)

summary(aveLogCPM(dgl))
summary(rowSums(dgl$count))

# estimate data normalisation factors and dispersion
xcpm <- mglmOneGroup(dgl$counts)     # computing a logCPM for making dispersion plot
dgl <- calcNormFactors(dgl)
dgl <- estimateGLMCommonDisp(dgl, dmat)
dgl <- estimateGLMTrendedDisp(dgl, dmat, min.n=1000)
dgl <- estimateGLMTagwiseDisp(dgl, dmat)

## dispersion plot
pdf(file.path(outpath, paste('Dispersion_', sub_analyse, '.pdf', sep="")), width=8, height=8)
par(mar=c(5,5,4,3))
plot(xcpm, dgl$tagwise.dispersion, pch=16, cex=0.5, xlab="log2CPM", ylab="Dispersion", main=paste(sub_analyse," dispersion", sep=""))
if(!is.null(dgl$trended.dispersion)) points(xcpm,dgl$trended.dispersion, pch=16, cex=0.5, col=cbgreen)
abline(h=dgl$common.dispersion ,col=cblblue, lwd=2)
legend("topright", c("Common","Trended","Tagwise"), pch=16, col=c(cblblue, cbgreen, "black"), title="Dispersion")
dev.off()

##  fit the data model ------------------------------------------
fitres <- glmFit(dgl, dmat, robust = TRUE, dispersion=0.7^2) ###command for samples without biological replicates

x <- read.delim(paste(datapath, sub_analyse, "_matrix.txt", sep=""), sep="\t", header=T)
sortedX <- data.frame(x[order(x$model_coefficients, decreasing=F),])
cmat <- as.matrix(sortedX[,-1])
colnames(cmat)[1] <- colnames(sortedX[2])
rownames(cmat) <- as.character(sortedX[,1])
# cmat

#Contrast fit and test results
lrtres <- list()
#for(k in 1:ncol(cmat)) lrtres[[k]] <- glmQLFTest(fitres, contrast=cmat[,k])
for(k in 1:ncol(cmat)) lrtres[[k]] <- glmLRT(fitres, contrast=cmat[,k])

logFC <- NULL
PV <- NULL
FDR <- NULL
for(k in 1:ncol(cmat)) {PV <- cbind(PV, lrtres[[k]]$table[,"PValue"])
FDR= cbind(FDR, p.adjust(PV[,k],method="BH"))
logFC= cbind(logFC, lrtres[[k]]$table[,"logFC"])
}

xcpm <- lrtres[[1]]$table[,"logCPM"]
allzeros <- which(rowSums(cpm(dgl)) < 3,)
allused <- which(rowSums(cpm(dgl)) >= 3,)

cname <-colnames(cmat)
colnames(logFC) <- paste("logFC", cname,sep=".")
colnames(PV) <- paste("PV", cname, sep=".")
colnames(FDR) <- paste("FDR", cname, sep=".")
rownames(logFC) <- rownames(PV) <- rownames(FDR) <- rownames(fitres$coefficients)

## Make results table
idxzeros <- allzeros
restab <- data.frame(rownames(dgl$counts),dgl$genes[,2], dgl$genes[,3], dgl$genes[,4], dgl$genes[,5],logCPM=xcpm,logFC,PV,FDR) # note start and end are the same in annotation file
colnames(restab)[1] <- 'gid'
colnames(restab)[2] <- 'gname'
colnames(restab)[3] <- 'chr'
colnames(restab)[4] <- 'start'
colnames(restab)[5] <- 'end'

write.table(restab$gid, file=file.path(outpath, paste(sub_analyse, "_used_gene_names.txt", sep="")), quote=F, row.names=F, sep="\t")

## pvalue histogram
pdf(file.path(outpath, paste('p_hist_',FDR2use, '_', sub_analyse,'.pdf', sep="")), width=8, height=8)
npanel <- ncol(logFC)
np <- ceiling(sqrt(npanel))
if(np*(np-1)>= npanel) mfcol <- c(np-1,np) else
{if((np-1)*(np+1)>= npanel) mfcol <- c(np+1,np-1) else mfcol <- c(np,np)}
par(mfcol=mfcol)
for(k in 1:npanel) {
  hist(PV[,k],n=100,xlab="P-value",main=colnames(logFC)[k])
}
dev.off()

## whinin group pairwise scatter plot
wx <- dgl$counts
wg <- as.character(dgl$samples$group)
wug <- unique(wg)
wn <- length(wug)
for (k in 1:wn) {
  ix <- wg %in% wug[k]
  xmat <- log2(wx[,ix])
  pdf(file.path(outpath, paste('pairwise_raw_count_', wug[k], '_', sub_analyse,'.pdf', sep="")), width=8, height=8)
  if (sum(ix) > 1 ) {
    pairs(xmat,pch=16,cex=0.4,main=wug[k])
  }
  #	dev.copy(pdf,file.path(outpath,'courtship_out',paste(sub_analyse,wug[k],'_pairwise_raw_count.pdf', sep="")), width=8, height=8)
  dev.off()
}


##Export number of DE genes table
restab_frame <- as.data.frame(restab)

for(logFC_use in c(3, 2, 1, 0) ) {
  de.yes.no <- FDR < FDR2use & abs(logFC) > logFC_use
  if (ncol(cmat) == 1) {
    de4 <- which((FDR[,1] < FDR2use) == 0 & abs(logFC[,1] > logFC_use) != 0) } else {
      de4 <- which(rowSums(FDR[,c(1:ncol(FDR))] < FDR2use) == 0 & rowSums(abs(logFC[,c(1:ncol(FDR))]) > logFC_use) != 0 )
    }

de.yes.no[de4,] <- FALSE
deidx <- ii <- rowSums(de.yes.no) > 0
delabel <- (sign(logFC)*de.yes.no)[ii,]
delabel[is.na(delabel)] <- 0
combinedp <- 1; for(k in 1:ncol(PV)) combinedp <- combinedp*PV[,k]
deidx <- ii

wtable_1 <- restab[deidx,] # this is the result table of DE genes.
demat <- as.matrix(logFC[deidx,])
write.table(wtable_1, file=file.path(outpath, paste('de_',FDR2use, '_',logFC_use, '_', sub_analyse,'.txt', sep="")), quote=F, row.names=F, sep='\t')

if (nrow(wtable_1) > 1) {
  if (ncol(cmat) == 1) {
    wtable_3 <- rbind(NUM_DE=sum(abs(delabel)), NUM_UP_DE =sum(delabel>0), NUM_DOWN_DE =sum(delabel<0)) } else {
      wtable_3 <- rbind(NUM_DE=colSums(abs(delabel)), NUM_UP_DE =colSums(delabel>0), NUM_DOWN_DE =colSums(delabel<0))
    }
wtable_3 <- cbind(DE_numbers=rownames(wtable_3), wtable_3)
write.table(wtable_3, file=file.path(outpath, paste('Number_de_',FDR2use, '_',logFC_use, '_', sub_analyse, '.txt', sep="")), quote=F, row.names=F, sep='\t')
  }
}

#logFC plot
fc <- logFC
fc[(fc)>10] <- 10
fc[ fc < -10] <- -10
par(mfcol=mfcol)

for(k in 1:ncol(cmat)) {
pdf(file.path(outpath, paste('FC-CPMplot_',FDR2use, '_', sub_analyse,'_', paste(colnames(cmat)[k]), '.pdf', sep="")), width=8, height=8)
par(mar=c(5,6,3,2)+0.1)
ylab <- colnames(logFC)[k]
deix <- which(de.yes.no[,k])
maPlot(x=NULL,y=NULL,ylim=c(-10,10),logAbundance= xcpm, logFC = fc[,k], xlab = bquote(paste(log^2, CPM)), ylab = paste(strsplit(colnames(cmat)[k], '\\.')[[1]][1], ' - ', strsplit(colnames(cmat)[k], '\\.')[[1]][2], sep=""), de.tags= deix, pch = 19, cex = 0.3, smearWidth = 0.5, cex.axis=1.8, cex.lab=2, panel.first = grid(), smooth.scatter = FALSE, lowess = FALSE, na.rm =TRUE, main = paste('LogFC plot ', strsplit(colnames(cmat)[k], '\\.')[[1]][1], ' vs ', strsplit(colnames(cmat)[k], '\\.')[[1]][2], sep=""))
dev.off()
}
