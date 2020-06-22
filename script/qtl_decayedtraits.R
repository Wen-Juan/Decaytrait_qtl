#install and load qtl packages
install.packages("qtl")
library(qtl)

### download the Rqtl code from the tutorial.
url.show("http://www.rqtl.org/rqtltour.R")

############################################################
# Example 1: Hypertension
############################################################
#load data
qtl <- read.cross("csv", dir = "/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/input/", file = "genotype_linkagemap_mod2.csv", alleles=c("A","B"),na.string="na")
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

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/genentic_map.pdf", width=8, height=8)
plotMap(qtl)
dev.off()

plotPheno(qtl, pheno.col=2)

#plot genetic map with marker names on certain chromosomes
plotMap(qtl, chr=c(1, 4), show.marker.names=TRUE)

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
newmap <- est.map(qtl, error.prob=0.001)
plotMap(qtl, newmap)

#If one wished to replace the genetic map with the estimated one, it could be done as follows:
qtl <- replace.map(qtl, newmap)
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
qtl <- sim.geno(qtl, step=1, n.draws=100, error.prob=0.001) 
out.imp <- scanone(qtl, pheno.col=1, addcovar=Cross, method="imp") #as binary model is not comptatible with imp method, so decided not to use this method.

#thefunctionsummary.scanonedisplaysthemaximumLODscoreon each chromosome for which the LOD exceeds a specified threshold
summary(out.em,pvalues=TRUE) #for some i-donot-understand reason, this values changes each time i ran the code
summary(out.em, threshold=3)
#chr pos  lod
c12.loc9  12   9 39.6
#chr pos  lod
c12.loc8  12   8 41.3

summary(out.hk, threshold=3)
# chr pos lod
c12.loc10  12  10  39
#  chr pos  lod
c12.loc10  12  10 40.9

summary(out.imp)
summary(out.imp, threshold=3)
# chr pos  lod
c12.loc9  12   9 39.6
#chr pos  lod
c12.loc2  12   2 50.8

max(out.em) 
#chr pos  lod
c12.loc8  12   8 41.3

max(out.hk) 
#chr pos  lod
c12.loc10  12  10 40.9

max(out.imp) 
# chr pos  lod
c12.loc2  12   2 50.8

#plot.scanone can plot up to three genome scans at once, provided that they conform appropriately. Alternatively, one may use the argument add.
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/qtl_lod_genomewide_em.pdf", width=8, height=8)
plot(out.em, chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), col="blue", ylim=c(0,40))
abline(h = 3,col="blue", lwd=1, lty=2)
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/qtl_lod_genomewide_hk.pdf", width=8, height=8)
plot(out.hk, chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), col="blue",ylim=c(0,40))
abline(h = 3,col="blue", lwd=1, lty=2)
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/qtl_lod_genomewide_all2.pdf", width=12, height=8)
plot(out.hk, chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),col="red",ylim=c(0,40))
plot(out.em, chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), col="blue", add=TRUE,ylim=c(0,42))
abline(h = 3,col="blue", lwd=1, lty=2)
dev.off()

#Thefunctionscanonemayalsobeusedtoperformapermutationtesttogetagenome-wideLODsignificancethreshold.
operm.hk <- scanone(qtl, method="hk", n.perm=100,addcovar=Cross,pheno.col=1, model ="binary")
operm.em <- scanone(qtl, method="em", n.perm=100,addcovar=Cross,pheno.col=1, model ="binary")
#operm.imp <- scanone(qtl, method="imp", n.perm=100, addcovar=Cross,model ="binary") #Method imp not available for binary model; using em

summary(operm.hk, alpha=0.05)
summary(out.hk, perms=operm.hk, alpha=0.05, pvalues=TRUE)
##
chr pos  lod pval
c12.loc10  12  10 40.9    0
##

summary(out.em, perms=operm.em, alpha=0.05, pvalues=TRUE)
#
chr pos  lod pval
c12.loc8  12   8 41.3    0

###the permutation results may also be used in teh summary.scanone function to calculate LOD thresholds and genome-scan-adjusted p values.
summary(out.hk, perms=operm.hk, alpha=0.05, pvalues=TRUE,format = "allpheno")
m1lg <- 12
m1cm <- 10
###
chr pos lod pval
c12.loc10  12  10  39    0
###


########################################
########################################two-QTL genome scan
########################################

# The function scantwo performs a two-dimensional genome scan with a two-QTL model. For every pair of positions, it calculates a LOD score for the full model (two QTL plus interaction) and a LOD score for the additive model (two QTL but no interaction). This be quite time consuming, and so you may wish to do the calculations on a coarser grid.
qtl <- calc.genoprob(qtl, step=1, error.prob=0.001)
Cross <- data.frame(Cross=as.numeric(pull.pheno(qtl, "cross") == "F2"))
###hk model
out2.hk <- scantwo(qtl, pheno.col=1, addcovar=Cross, model=c("binary"),method="hk")
summary(out2.hk)
max(out2.hk)

plot(out2.hk)
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/qtl_hk_scantwo_lod.pdf", width=12, height=12)
plot(out2.hk, chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16))
dev.off()
#The function max.scantwo returns the two-locus positions with the maximum LOD score for the full and additive models.
max(out2.hk)
###02062020
###
pos1f pos2f lod.full lod.fv1 lod.int     pos1a pos2a lod.add lod.av1
c2:c12    43     9     40.2    1.14   0.446        44     9    39.7   0.693
###
19062020
#       pos1f pos2f lod.full lod.fv1 lod.int     pos1a pos2a lod.add lod.av1
c3:c12    34    11     41.8   0.892   0.286        35    10    41.5   0.606

#One may also use scantwo to perform permutation tests in order to obtain genome-wide LOD significance thresholds.
operm2.hk <- scantwo(qtl, method="hk", n.perm=10, model='binary')

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
summary(out2.hk, perms=operm2.hk, pvalues=TRUE,
        alphas=0.05)
##
There were no pairs of loci meeting the criteria.
##
plot(operm2.hk)

###########
###########

# end of rqtltour.R

###pie chart
# Load ggplot2
library(ggplot2)

# Create Data
data <- data.frame(
  group=LETTERS[1:2],
  value=c(93.75,6.25)
)

# Basic piechart
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/piechartGO.pdf", width=12, height=12)
ggplot(data, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual(values=c("black", "white")) +
  theme(legend.position="none") 
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/piecharG1.pdf", width=12, height=12)
ggplot(data, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual(values=c("black", "white")) +
  theme(legend.position="none") 
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/piecharG2.pdf", width=12, height=12)
ggplot(data, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual(values=c("black", "white")) +
  theme(legend.position="none") 
dev.off()


pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/piecharG3.pdf", width=12, height=12)
ggplot(data, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual(values=c("black", "white")) +
  theme(legend.position="none") 
dev.off()


pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/piecharG4.pdf", width=12, height=12)
ggplot(data, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=1, color="black") +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual(values=c("black", "white")) +
  theme(legend.position="none") 
dev.off()
