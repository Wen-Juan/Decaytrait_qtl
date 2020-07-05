#install and load qtl packages
install.packages("qtl")
library(qtl)

#also incoporate other R packages to make nice ggplot2 figures.
library(ggplot2)
library(ggrepel)
library(tidyr)
library(dplyr)

### download the Rqtl code from the tutorial.
url.show("http://www.rqtl.org/rqtltour.R")

############################################################
# Example 1: Hypertension
############################################################
#load data
qtl <- read.cross("csv", dir = "/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/input/", file = "genotype_linkagemap_mod3.csv", alleles=c("A","B"),na.string="na")
qtl <- jittermap(qtl)
summary(qtl) 

#some summary information of the dataset
nind(qtl)
nphe(qtl)
nchr(qtl)
totmar(qtl)
nmar(qtl)

# remove duplicate markers
print(dup <- findDupMarkers(qtl, exact.only=FALSE))
qtl <- drop.markers(qtl, unlist(dup))

# should report no problems
qtl <- est.rf(qtl)
checkAlleles(qtl, threshold=5)

#plot the summary
plot(qtl)

#plot individual perspect of information
plotMissing(qtl)
par(mfrow=c(1,1)) 
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/genentic_map_original.pdf", width=8, height=8)
plotMap(qtl,show.marker.names=FALSE,alternate.chrid=TRUE,horizontal=FALSE)
dev.off()

plotPheno(qtl, pheno.col=2)

#plot genetic map with marker names on certain chromosomes
plotMap(qtl, chr=c(7), show.marker.names=TRUE)

plotMissing(qtl, reorder=TRUE)

#remove markers having no genotype data
qtl <- drop.nullmarkers(qtl)
totmar(qtl)

#estimating recombination fractions between all pairs of markers and plot them.
qtl <- est.rf(qtl)
plotRF(qtl)
plotRF(qtl, chr=c(12))

plotRF(qtl, chr=6)
plotMissing(qtl, chr=6)

#Re-estimate the genetic map (keeping the order of markers fixed), and plot the original map against the newly estimated one.
newmap <- est.map(qtl, error.prob=0.001, m=0, p=0, map.function="kosambi")
plotMap(qtl, newmap)
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/genentic_map_originalvsnew.pdf", width=8, height=8)
plotMap(qtl, newmap,show.marker.names=FALSE,horizontal=FALSE,,xlab="Linkage group", ylab="Location (cM)", col = "dark green", lwd =2,cex.lab=5)
dev.off()

#If one wished to replace the genetic map with the estimated one, it could be done as follows:
qtl <- replace.map(qtl, newmap)
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/genentic_map_new.pdf", width=8, height=8)
plotMap(qtl,show.marker.names=FALSE,alternate.chrid=TRUE,horizontal=FALSE)
dev.off()

##error message here:Error in as.data.frame.default(x[[i]], optional = TRUE) : cannot coerce class ""A"" to a data.frame

#identification of genotyping errors
qtl <- calc.errorlod(qtl, error.prob=0.001)

#genotyping value <4 can be ignored.
top.errorlod(qtl)

#as top genotyping errors show bias in chr05, we plot it. ThefunctionplotGenomaybeusedtoinspecttheobservedgenotypesforachromosome,withlikelygenotypingerrors flagged.
plotGeno(qtl, chr=5, ind=c(24:34, 71:86))

#The function plotInfo plots a measure of the proportion of missing genotype information in the genotype data.
?plotInfo
plotInfo(qtl)
plotInfo(qtl, chr=c(1,2,3))
plotInfo(qtl, chr=c(1,4,15), method="entropy")
plotInfo(qtl, chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), method="variance")

########################################
########################################single-QTL genome scan
########################################
#The function calc.genoprob calculates QTL genotype probabilities, conditional on the available marker data
qtl <- calc.genoprob(qtl, step=1, error.prob = 0.001)
Cross <- data.frame(Cross=as.numeric(pull.pheno(qtl, "cross") == "F2"))
out.em <- scanone(qtl, pheno.col=1, addcovar=Cross, model=c("binary"),method="em")
out.hk <- scanone(qtl, pheno.col=1, addcovar=Cross, model=c("binary"),method="hk")


#The argument step indicates the step size (in cM) at which the probabilities are calculated, and determines the step size at which later LOD scores are calculated.
#We may also use the multiple imputation method of Sen and Churchill (2001). The n.draws indicates the number of imputations to perform. 
#step indicates the step size (in cM) in which the probability is calculated.
qtl <- sim.geno(qtl, step=1, n.draws=1000, error.prob=0.001) 
out.imp <- scanone(qtl, pheno.col=1, addcovar=Cross, model=c("binary"),method="imp") 
#as binary model is not comptatible with imp method, so decided not to use this method.

#thefunctionsummary.scanonedisplaysthemaximumLODscoreon each chromosome for which the LOD exceeds a specified threshold
summary(out.em,pvalues=TRUE) #for some i-donot-understand reason, this values changes each time i ran the code
summary(out.em, threshold=3)
#chr pos  lod
c12.loc9  12   9 39.6
#chr pos  lod
c12.loc8  12   8 41.3
##
chr pos  lod
c12.loc6  12   6 41.2
##
chr pos  lod
c12.loc7  12   7 41.3

summary(out.hk, threshold=3)
# chr pos lod
c12.loc10  12  10  39
#  chr pos  lod
c12.loc10  12  10 40.9
###chr pos  lod
c12.loc7  12   7 40.7
###
chr pos  lod
c12.loc8  12   8 40.8

summary(out.imp)
summary(out.imp, threshold=3)
# chr pos  lod
c12.loc9  12   9 39.6
#chr pos  lod
c12.loc2  12   2 50.8
###
chr pos  lod
c12.loc6  12   6 41.2
###
chr pos  lod
c12.loc7  12   7 41.3

max(out.em) 
#chr pos  lod
c12.loc8  12   8 41.3
#    chr pos  lod
c12.loc6  12   6 41.2
#  chr pos  lod
c12.loc7  12   7 41.3

max(out.hk) 
#chr pos  lod
c12.loc10  12  10 40.9
## chr pos  lod
c12.loc7  12   7 40.7
##chr pos  lod
c12.loc8  12   8 40.8

max(out.imp) 
# chr pos  lod
c12.loc2  12   2 50.8
#   chr pos  lod
c12.loc6  12   6 41.2
#chr pos  lod
c12.loc7  12   7 41.3

#plot.scanone can plot up to three genome scans at once, provided that they conform appropriately. Alternatively, one may use the argument add.
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/qtl_lod_genomewide_em.pdf", width=8, height=8)
plot(out.em, chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), col="blue", ylim=c(0,40),xlab="Linkage group", ylab="LOD",lwd = 2,alternate.chrid=TRUE)
abline(h = 3,col="blue", lwd=1, lty=2)
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/qtl_lod_genomewide_hk.pdf", width=8, height=8)
plot(out.hk, chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), col="blue",ylim=c(0,40),xlab="Linkage group", ylab="LOD",lwd = 2,alternate.chrid=TRUE)
abline(h = 3,col="blue", lwd=1, lty=2)
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/qtl_lod_genomewide_all2.pdf", width=12, height=8)
plot(out.hk, chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),col="red",ylim=c(0,40),xlab="Linkage group", ylab="LOD",lwd = 2,alternate.chrid=TRUE)
plot(out.em, chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), col="blue", add=TRUE,ylim=c(0,42),xlab="Linkage group", ylab="LOD",lwd = 2,alternate.chrid=TRUE)
abline(h = 3,col="blue", lwd=2, lty=2)
dev.off()

#Thefunctionscanonemayalsobeusedtoperformapermutationtesttogetagenome-wideLODsignificancethreshold.
operm.hk <- scanone(qtl, method="hk", n.perm=1000,addcovar=Cross,pheno.col=1, model ="binary")
operm.em <- scanone(qtl, method="em", n.perm=1000,addcovar=Cross,pheno.col=1, model ="binary")
#operm.imp <- scanone(qtl, method="imp", n.perm=100, addcovar=Cross,model ="binary") #Method imp not available for binary model; using em
summary(operm.hk, alpha=0.05)
#####
LOD thresholds (1000 permutations)
lod
5% 2.47
#####
LOD thresholds (1000 permutations)
lod
5% 2.4
###
LOD thresholds (1000 permutations)
lod
5% 2.36

lod_threshold <- summary(operm.hk, alpha = 0.05) 
labels_df <- as.data.frame(summary(out.hk, perms = operm.hk, alpha = 0.05, pvalues = TRUE))

#################
#################
#################
summary(out.hk, perms=operm.hk, lodcolumn = 1, alpha=0.05, pvalues=TRUE) ### alpha=0.05 5% is above this threshold. 95% are below. 
##Permutation = 1000
chr pos  lod pval
c12.loc10  12  10 40.9    0
##
chr pos  lod pval
c12.loc8  12   8 40.8    0
#################
#################
#################

summary(out.hk, perms=operm.hk, lodcolumn = 1,alpha=0.10, pvalues=TRUE)
##Permutation = 1000
chr pos  lod pval
c12.loc10  12  10 40.9    0
##
chr pos  lod pval
c12.loc8  12   8 40.8    0

summary(out.hk, perms=operm.hk, lodcolumn = 1,alpha=0.60, pvalues=TRUE)
##Permutation = 1000
chr pos   lod  pval
c5.loc47    5  47  1.56 0.364
c9.loc45    9  45  1.53 0.367
c12.loc10  12  10 40.89 0.000
##
chr pos   lod  pval
c5.loc40   5  40  1.55 0.346
c9.loc56   9  56  1.53 0.360
c12.loc8  12   8 40.80 0.000

summary(out.em, perms=operm.em, lodcolumn = 1,alpha=0.05, pvalues=TRUE)
##Permutation = 1000
chr pos  lod pval
c12.loc8  12   8 41.3    0
#     chr pos  lod pval
c12.loc7  12   7 41.3    0

summary(out.em, perms=operm.em, lodcolumn = 1,alpha=0.10, pvalues=TRUE)
##Permutation = 1000
chr pos  lod pval
c12.loc8  12   8 41.3    0
# chr pos  lod pval
c12.loc7  12   7 41.3    0

summary(out.em, perms=operm.em, lodcolumn = 1,alpha=0.60, pvalues=TRUE)
##Permutation = 1000
chr pos   lod  pval
c5.loc47   5  47  1.55 0.331
c9.loc45   9  45  1.53 0.344
c12.loc8  12   8 41.30 0.000
# chr pos   lod  pval
c5.loc40   5  40  1.55 0.336
c9.loc56   9  56  1.53 0.348
c12.loc7  12   7 41.28 0.000

mname1 <- find.marker(qtl, chr = 12, pos = 7)
mname1 
#[1] "Ajap26"
#[1] "Ajap26"
mname2 <- find.marker(qtl, chr = 12, pos = 8)
mname2 
#[1] "Ajap26"
#[1] "Ajap26"
mname3 <- find.marker(qtl, chr = 12, pos = 2)
mname3 
#[1] "Ajap26"
#[1] "Ajap26"
mname4 <- find.marker(qtl, chr = 12, pos = 11)
mname4 
#[1] "Ajap26"
#[1] "Ajap53"
effectplot(qtl, pheno.col = 1, mname1=mname1)
###heterozygosity is responsible for mating success.

effectplot(qtl, pheno.col = 1, mname1=mname4)
###heterozygosity is responsible for mating success.


pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/qtl_lod_2qtl_effectplotv2.pdf",  width=12, height=6)
par(mfrow=c(1,2)) 
par(mar=c(5,5,4,3))
effectplot(qtl, pheno.col = 1, mname1=mname1, draw=TRUE, xlab="Linkage group 12", ylab ="Proportion mating success",ylim = c(0,1))
effectplot(qtl, pheno.col = 1, mname1=mname4, draw=TRUE,xlab="Linkage group 12", ylab ="Proportion mating success",ylim = c(0,1))
dev.off()


CIchr12v1 <- bayesint(out.hk, chr = 12, prob = 0.95)
####
chr pos      lod
c12.loc7   12   7 40.02445
c12.loc10  12  10 40.88967
c12.loc13  12  13 40.03487
#####
chr pos      lod
c12.loc6   12   6 40.08095
c12.loc8   12   8 40.80095
c12.loc11  12  11 39.90136
####
####
plot(out.hk, chr=12, lodcolumn=1, main="Confidence interval for QTL at LG12", xlab="Linkage map", ylab="LOD", col = "blue", lwd =2)
lines(x=CIchr12v1[c(1,3),2], y=c(0,0), type="l", col="blue", lwd=4)

CIchr12v2 <- bayesint(out.em, chr = 12, prob = 0.95)
###
chr pos      lod
c12.loc5   12   5 40.59870
c12.loc8   12   8 41.29667
c12.loc12  12  12 40.70013
###

###
chr pos      lod
c12.loc4   12   4 40.49679
c12.loc7   12   7 41.28363
c12.loc10  12  10 40.69709

###
par(mfrow=c(1,1)) 
plot(out.em, chr=12, lodcolumn=1, main="Confidence interval for QTL at LG12", xlab="Linkage map", ylab="LOD", col = "dark green", lwd =2)
lines(x=CIchr12v2[c(1,3),2], y=c(0,0), type="l", col="dark green", lwd=4)


pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/qtl_lod_genomewide_all2.pdf", width=12, height=8)
plot(out.hk, chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),col="dark green",ylim=c(0,40),xlab="Linkage group", ylab="LOD",lwd = 2,alternate.chrid=TRUE)
plot(out.em, chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), col="blue", add=TRUE,ylim=c(0,42),xlab="Linkage group", ylab="LOD",lwd = 2,alternate.chrid=TRUE)
abline(h = 2.47,col="blue", lwd=2, lty=2)
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/qtl_lod_genomewide_lg12_confidenceinterval.pdf", width=10, height=8)
plot(out.em, chr=12, lodcolumn=1, main="Confidence interval for QTL at LG12", xlab="Linkage group 12", ylab="LOD", col = "dark green", lwd =2,lternate.chrid=TRUE,cex.lab=1.5)
plot(out.hk, chr=12, lodcolumn=1, xlab="Linkage map", ylab="LOD", col = "blue", lwd =2,add=TRUE,ylim=c(0,42),lternate.chrid=TRUE,cex.lab=1.5)
lines(x=CIchr12v2[c(1,3),2], y=c(0,0), type="l", col="dark green", lwd=4)
lines(x=CIchr12v1[c(1,3),2], y=c(0,0), type="l", col="blue", lwd=4)
abline(h = 2.4,col="black", lwd=1.5, lty=2)
dev.off()

###########
###########
###
#dataManual <- data.frame(out.em)
#qtl_plot(input=dataManual) This below qtl_plot function comes from URL:https://www.r-bloggers.com/conditional-ggplot2-geoms-in-functions-qtl-plots/
qtl_plot <- function(input,              # data frame input from scanone
                     mult.pheno = FALSE, # multiple phenotypes?
                     model = "binary",   # model used in scanone
                     chrs = NA,          # chromosomes to display
                     lod = 2.47,           # LOD threshold
                     rug = FALSE,        # plot marker positions as rug?
                     ncol = NA,          # number of columns for facetting
                     labels = NA         # optional dataframe to plot QTL labels
) {
  
  # if we have multiple phenotypes and/or a 2part model, gather input
  if (mult.pheno & model == "2part") {
    input <- gather(input, group, lod, grep("pheno", colnames(input)))
  } else if (mult.pheno) {
    input <- gather(input, group, lod, grep("pheno", colnames(input)))
  } else if (model == "2part") {
    input <- gather(input, method, lod, lod.p.mu:lod.mu)
  }
  
  # if not all chromosomes should be displayed, subset input
  if (!is.na(chrs)[1]) {
    input <- input[as.character(input$chr) %in% chrs, ]
  }
  
  # if there is more than one LOD column, gather input
  if (!any(colnames(input) == "lod")) {
    input$lod <- input[, grep("lod", colnames(input))]
  }
  
  # if no number of columns for facetting is defined, plot all in one row
  if (is.na(ncol)) {
    ncol <- length(unique(input$chr))
  }
  
  # if labels are set and there is no name column, set from rownames
  if (!is.na(labels)[1]) {
    if (is.null(labels$name)) {
      labels$name <- rownames(labels)
    }
  }
  
  # plot input data frame position and LOD score
  plot <- ggplot(input, aes(x = pos, y = lod)) + {
    
    # if LOD threshold is given, plot as horizontal line
    if (!is.na(lod)[1] & length(lod) == 1) geom_hline(yintercept = lod, linetype = "dashed")
  } + {
    
    if (!is.na(lod)[1] & length(lod) > 1) geom_hline(data = lod, aes(yintercept = lod, linetype = group))
  } + {
    
    # plot rug on bottom, if TRUE
    if (rug) geom_rug(size = 0.1, sides = "b")
  } + {
    
    # if input has column method but not group, plot line and color by method
    if (!is.null(input$method) & is.null(input$group)) geom_line(aes(color = method), size = 1, alpha = 0.8)
  } + {
    
    # if input has column group but not method, plot line and color by group
    if (!is.null(input$group) & is.null(input$method)) geom_line(aes(color = group), size = 1, alpha = 0.8)
  } + {
    
    # if input has columns method and group, plot line and color by method & linetype by group
    if (!is.null(input$group) & !is.null(input$method)) geom_line(aes(color = method, linetype = group), size = 1, alpha = 0.8)
  } + {
    
    # set linetype, if input has columns method and group
    if (!is.null(input$group) & !is.null(input$method)) scale_linetype_manual(values = c("solid", "twodash", "dotted"))
  } + {
    
    # if input has neither columns method nor group, plot black line
    if (is.null(input$group) & is.null(input$method)) geom_line(size = 1, alpha = 0.8)
  } + {
    
    # if QTL positions are given in labels df, plot as point...
    if (!is.na(labels)[1]) geom_point(data = labels, aes(x = pos, y = lod))
  } + {
    
    # ... and plot name as text with ggrepel to avoid overlapping
    if (!is.na(labels)[1]) geom_text_repel(data = labels, aes(x = pos, y = lod, label = "loc8"),nudge_y = 0.7, size=4.5) 
  } + 
    # facet by chromosome
   facet_wrap(~ chr, ncol = ncol, scales = "free_x") +
   theme(axis.title.x = element_text(size=15,colour = "black"),axis.title.y = element_text(size=15,colour = "black")) +
   theme(axis.text.x = element_text(colour="black",size=11,angle=90)) +
   theme(axis.text.y = element_text(colour="black",size=11)) +
    # minimal plotting theme
    #theme_minimal() +
    # increase strip title size
    theme(strip.text = element_text(face = "bold", size = 13)) +
    theme(legend.title = element_blank()) + 
    theme(legend.position = c(0.2, 0.8)) +
    scale_color_manual(values=c("#009900", "#0066FF")) +
    # use RcolorBrewer palette
 #   scale_color_brewer(palette = "Set1") +
    # Change plot labels
    labs(x = "Linkage map",
         y = "LOD",
         color = "black",
         linetype = "twodash")
  
  print(plot)
}

###to hide the x axis tick mark lables (but keep the marks)
#theme(axis.text.x = element_blank())

qtl_plot(input=out.em)

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/qtl_lod_genomewide_permutation1000_2modles.pdf", width=18, height=8)
qtl_plot(input = rbind(data.frame(out.em, method = "EM algorithm"), 
                       data.frame(out.hk, method = "Haley-Knott regression")), 
         lod = lod_threshold[1], 
         rug = TRUE, 
         labels = labels_df)
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/qtl_lod_lg12_permutation1000_2modles.pdf", width=10, height=8)
qtl_plot(input = rbind(data.frame(out.em, method = "EM algorithm"), 
                       data.frame(out.hk, method = "Haley-Knott regression")), 
         lod = lod_threshold[1], 
         chrs = c(12),
         rug = TRUE, 
         labels = labels_df) 
dev.off()

###fir the model for single QTL
#operm.em
multisim <- sim.geno(qtl, step=1, n.draws=10000, err=0.001)
multisim2 <- calc.genoprob(multisim, step=1, error.prob=0.001)

m1lg <- 12
m1cm <- 7
#m2lg <- 12
#m2cm <- 23

Cross <- data.frame(Cross=as.numeric(pull.pheno(qtl, "cross") == "BC"))
f1 <- "y ~ Cross + Q1 + Cross:Q1"

# Input number of loci for the qtl object
qtl.auto <- makeqtl(multisim2, chr=c(m1lg), pos=c(m1cm), what="draws")
fq.auto <- fitqtl(multisim2, cov=Cross,qtl=qtl.auto, formula=f1,model='binary')
summary(fq.auto, alpha=0.05, pvalues=TRUE, format="allpheno")

###
fitqtl summary

Method: multiple imputation 
Model:  binary phenotype
Number of observations : 253 

Full model result
----------------------------------  
  Model formula: y ~ Cross + Q1 + Cross:Q1 

df     LOD     %var Pvalue(Chi2)
Model  3 36.6254 48.65819            0


Drop one QTL at a time ANOVA table: 
  ----------------------------------  
  df   LOD  %var Pvalue(Chi2)    
Cross         2  0.00  0.00            1    
12@7.0        2 36.63 48.66       <2e-16 ***
  Cross:12@7.0  1  0.00  0.00            1    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
###
# check for interactions between qtl or with qtl and cross
qtl.auto.ai <- addint(multisim2, qtl=qtl.auto, formula=f1, pheno.col=1, covar=Cross)
summary(qtl.auto.ai, alpha=0.05, pvalues=TRUE, format="allpheno")

# refine positions
rqtl <- refineqtl(multisim2, qtl=qtl.auto, formula=f1, pheno.col=1, method=c("imp"), covar=Cross,model='binary')
final.fit <- fitqtl(multisim2, qtl=rqtl, formula=f1, pheno.col=1, covar=Cross, model='binary')
summary(final.fit, alpha=0.05, pvalues=TRUE, format="allpheno")
###
fitqtl summary

Method: multiple imputation 
Model:  binary phenotype
Number of observations : 253 

Full model result
----------------------------------  
  Model formula: y ~ Cross + Q1 + Cross:Q1 

df      LOD     %var Pvalue(Chi2)
Model  3 36.75003 48.77454            0


Drop one QTL at a time ANOVA table: 
  ----------------------------------  
  df   LOD  %var Pvalue(Chi2)    
Cross         2  0.00  0.00            1    
12@6.0        2 36.75 48.77       <2e-16 ***
  Cross:12@6.0  1  0.00  0.00            1    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
###

capture.output(summary(final.fit, alpha=0.05, pvalues=TRUE, format="allpheno"), file=file.path(outpath, paste("pheno",'2qtl_PVE_table.txt', sep="_")))


########################################
########################################two-QTL genome scan
########################################

# The function scantwo performs a two-dimensional genome scan with a two-QTL model. For every pair of positions, it calculates a LOD score for the full model (two QTL plus interaction) and a LOD score for the additive model (two QTL but no interaction). This be quite time consuming, and so you may wish to do the calculations on a coarser grid.
qtl <- calc.genoprob(qtl, step=1, error.prob=0.001)
Cross <- data.frame(Cross=as.numeric(pull.pheno(qtl, "cross") == "BC"))
###hk model
out2.hk <- scantwo(qtl, pheno.col=1, addcovar=Cross, model=c("binary"),method="hk")
summary(out2.hk)

plot(out2.hk) #here is the heatmap of LOD scores.
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/qtl_hk_scantwo_lod.pdf", width=12, height=12)
plot(out2.hk, chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16))
dev.off()
#The function max.scantwo returns the two-locus positions with the maximum LOD score for the full and additive models.
max(out2.hk)
###02062020
###
pos1f pos2f lod.full lod.fv1 lod.int     pos1a pos2a lod.add lod.av1
c2:c12    43    10     40.2    1.16   0.465        44    10    39.7   0.696
###

#One may also use scantwo to perform permutation tests in order to obtain genome-wide LOD significance thresholds.
operm2.hk <- scantwo(qtl, pheno.col=1, addcovar=Cross, method="hk", n.perm=5, model='binary')

summary(operm2.hk,alpha=0.05)
##
mating (1000 permutations)
full  fv1  int  add  av1  one
5%  4.90 3.83 3.61 3.87 1.89 2.34
10% 4.57 3.52 3.27 3.41 1.69 2.03
###permutation with 10
full  fv1  int  add  av1  one
5% 4.95 3.15 2.92 4.69 2.21 2.51
##
summary(out2.hk,thresholds=c(4.90, 3.83, 3.61, 3.87, 1.89),df=TRUE) ###5%
#There were no pairs of loci meeting the criteria.
summary(out2.hk,thresholds=c(4.57, 3.52, 3.27, 3.41, 1.69),df=TRUE) ###10%
###
pos1f pos2f lod.full lod.fv1 lod.int     pos1a pos2a lod.add lod.av1
c5:c9    49    41     4.56    2.77    1.02        48    42    3.54    1.75
###
summary(out2.hk, perms=operm2.hk,alphas=c(0.05, 0.05, 0, 0.05, 0.05), pvalues=TRUE)
##
pos1f pos2f lod.full pval lod.fv1 pval lod.int pval     pos1a pos2a lod.add pval lod.av1 pval
c5:c9    49    41     4.56  0.2    2.77    1    1.02    1        48    42    3.54    0    1.75  0.2
##
plot(operm2.hk) 
plot(out2.hk, chr=c(12),upper="cond-int")

###One may also restrict the summary to just the case of j = k, to look at evidence for linked QTL on each chromosome, by using allpairs=FALSE.
summary(out2.hk,allpairs=FALSE, df=TRUE)
######
pos1f pos2f lod.full lod.fv1 lod.int     pos1a pos2a lod.add lod.av1
c1 :c1     50    55    2.721   1.546  1.3432        27    87   1.378  0.2031
c2 :c2      0    51    1.598   0.451  0.3613         0    51   1.237  0.0894
c3 :c3      1     2    2.657   2.412  2.2577         6     7   0.400  0.1539
c4 :c4     51    53    0.522   0.100  0.0812        64    68   0.441  0.0191
c5 :c5     12    65    3.180   1.387  1.1686        19    58   2.011  0.2183
c6 :c6      0     4    1.071   0.717  0.5892         1     5   0.482  0.1276
c7 :c7      0   256    1.822   1.597  0.9610        26    27   0.861  0.6363
c8 :c8      0    82    1.229   0.311  0.1544         0    79   1.074  0.1569
c9 :c9      3    68    2.945   1.285  0.8772        47    50   2.068  0.4077
c10:c10    68    93    1.909   0.748  0.3087        17    56   1.600  0.4396
c11:c11    10    11    2.569   1.511  1.4377         7    12   1.132  0.0734
c12:c12     0    26   39.209   0.190  0.1521         6    28  39.057  0.0382
c13:c13     0     0    0.588   0.000  0.0000         0     0   0.588  0.0000
c14:c14     7    10    0.615   0.190  0.0807         1     9   0.535  0.1096
c15:c15     0     0    0.124   0.000  0.0000         0     0   0.124  0.0000
c16:c16     0     1    0.614   0.555  0.3490         0     1   0.265  0.2056
####

###########
###########
# end of rqtltour.R

