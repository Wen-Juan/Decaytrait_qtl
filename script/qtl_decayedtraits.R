#install and load qtl packages
install.packages("qtl")
library(qtl)

### download the Rqtl code from the tutorial.
url.show("http://www.rqtl.org/rqtltour.R")

############################################################
# Example 1: Hypertension
############################################################
#load data
qtl <- read.cross("csv", dir = "/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/input/", file = "genotype_linkagemap_mod.csv", alleles=c("A","B"),na.string="na")
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

plotPheno(qtl, pheno.col=1)

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

#The function calc.genoprob calculates QTL genotype probabilities, conditional on the available marker data
#The argument step indicates the step size (in cM) at which the probabilities are calculated, and determines the step size at which later LOD scores are calculated.
qtl <- calc.genoprob(qtl, step=1, error.prob=0.001)

########################################
########################################single-QTL genome scan
########################################
# use the function scanone to perform a single-QTL genome scan with a binomial model.
out.em <- scanone(qtl,method="em",model='binary')
out.hk <- scanone(qtl, method="hk",model='binary')

#We may also use the multiple imputation method of Sen and Churchill (2001). The n.draws indicates the number of imputations to perform. 
#step indicates the step size (in cM) in which the probability is calculated.

qtl <- sim.geno(qtl, step=0.2, n.draws=100, error.prob=0.001)
out.imp <- scanone(qtl, method="imp",model='binary')

#thefunctionsummary.scanonedisplaysthemaximumLODscoreon each chromosome for which the LOD exceeds a specified threshold
summary(out.em,pvalues=TRUE) #for some i-donot-understand reason, this values changes each time i ran the code
summary(out.em, threshold=3)
#chr pos  lod
c12.loc9  12   9 39.6

summary(out.hk)
summary(out.hk, threshold=2)
# chr pos lod
c12.loc10  12  10  39

summary(out.imp)
summary(out.imp, threshold=3, alpha=0.05,pvalues = TRUE)
# chr pos  lod
c12.loc9  12   9 39.6

max(out.em) #c12.loc9  12   9 39.6
max(out.hk) #c12.loc10  12  10  39
max(out.imp) #c12.loc9  12   9 39.6

#plot.scanone can plot up to three genome scans at once, provided that they conform appropriately. Alternatively, one may use the argument add.
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/qtl_lod_genomewide_em.pdf", width=8, height=8)
plot(out.em, chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),ylim=c(0,42))
abline(h = 3,col="blue", lwd=1, lty=2)
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/qtl_lod_genomewide_hk.pdf", width=8, height=8)
plot(out.em, out.hk, out.imp, chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),ylim=c(0,42))
plot(out.hk, chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), col="blue", add=TRUE,ylim=c(0,42))
abline(h = 3,col="blue", lwd=1, lty=2)
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/qtl_lod_genomewide_imp.pdf", width=8, height=8)
plot(out.em, out.hk, out.imp, chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),ylim=c(0,42))
plot(out.imp, chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),ylim=c(0,50),col="red", add=TRUE)
abline(h = 3,col="blue", lwd=1, lty=2)
dev.off()

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/qtl_lod_genomewide_imp.pdf", width=8, height=8)
plot(out.em, out.hk, out.imp, chr=c(1,11,12),ylim=c(0,42))
plot(out.imp, chr=c(1,11,12),ylim=c(0,50),col="red", add=TRUE)
abline(h = 3,col="blue", lwd=1, lty=2)
dev.off()

#Thefunctionscanonemayalsobeusedtoperformapermutationtesttogetagenome-wideLODsignificancethreshold.
operm.hk <- scanone(qtl, method="hk", n.perm=10000,pheno.col=1, model ="binary")
operm.hk <- scantwo(qtl, method="hk", n.perm=10, pheno.col=1, model ="binary")
operm.imp <- scanone(qtl, method="imp", n.perm=3000)

summary(operm.hk, alpha=0.05)
summary(operm.imp, alpha=0.05)
##
LOD thresholds (10000 permutations)
lod
5%  2.42
10% 2.12
##

###the permutation results may also be used in teh summary.scanone function to calculate LOD thresholds and genome-scan-adjusted p values.
summary(out.hk, perms=operm.hk, alpha=0.05, pvalues=TRUE,format = "allpheno")
m1lg <- 12
m1cm <- 10
###
chr pos lod pval
c12.loc10  12  10  39    0
###

summary(out.imp, perms=operm.imp, alpha=0.05, pvalues=TRUE)
###
chr pos  lod pval
c12.loc4  12   4 54.5    0
###

###obtaining a bootstrap-based confidence interval for the location of a qtl with the function scanoneboot.
out.boot <- scanoneboot(qtl, chr = c(12), n.boot = 1000)

# The function scantwo performs a two-dimensional genome scan with a two-QTL model. For every pair of positions, it calculates a LOD score for the full model (two QTL plus interaction) and a LOD score for the additive model (two QTL but no interaction). This be quite time consuming, and so you may wish to do the calculations on a coarser grid.
qtl <- calc.genoprob(qtl, step=5, error.prob=0.01)

out2.hk <- scantwo(qtl, method="hk")

summary(out2.hk, thresholds=c(6.0, 4.7, 4.4, 4.7, 2.6)) #backcross

summary(out2.hk, thresholds=c(6.0, 4.7, Inf, 4.7, 2.6)) #intercross

plot(out2.hk)
plot(out2.hk, chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16))

#The function max.scantwo returns the two-locus positions with the maximum LOD score for the full and additive models.
max(out2.hk)
###02062020
pos1f pos2f lod.full lod.fv1 lod.int     pos1a pos2a lod.add lod.av1
c12:c14    10     0     50.9    0.97  0.0306        10     0    50.9   0.939
###
pos1f pos2f lod.full lod.fv1 lod.int     pos1a pos2a lod.add lod.av1
c2:c12    40     5     51.3       1    0.33        40     5      51   0.672
###

#One may also use scantwo to perform permutation tests in order to obtain genome-wide LOD significance thresholds.
operm2.hk <- scantwo(qtl, method="hk", n.perm=1000)

summary(operm2.hk)
##
mating (1000 permutations)
full  fv1  int  add  av1  one
5%  4.90 3.83 3.61 3.87 1.89 2.34
10% 4.57 3.52 3.27 3.41 1.69 2.03
##

summary(out2.hk, perms=operm2.hk, pvalues=TRUE,
        alphas=c(0.05, 0.05, 0, 0.05, 0.05))
##
There were no pairs of loci meeting the criteria.
##
##we consider the fit of multiple-QTL models.
chr <- c(1, 2, 4, 6, 15)
pos <- c(50, 76, 30, 70, 20)
qtl <- makeqtl(qtl, chr, pos)

my.formula <- y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q4:Q5
out.fitqtl <- fitqtl(qtl, qtl=qtl, formula=my.formula)
summary(out.fitqtl)

ls()
rm(list=ls())


########################################
########################################two-QTL genome scan
########################################
qtl <- calc.genoprob(qtl, step=1, error.prob=0.001)
qtl <-sim.geno(qtl,step=11,n.draws=100,err=0.001)

out2.em <- scantwo(qtl,method="em",model='binary')
out2.hk <- scanone(qtl, method="hk",model='binary')
out2.imp <- scanone(qtl, method="imp",model='binary')


############################################################
# Example 5: Multiple QTL mapping
###########################################################
out2.c4 <- scantwo(hyper, method="imp", addcovar=g, chr=c(1,6,7,15))

summary(out2.c4, thr=c(6.0, 4.7, Inf, 4.7, 2.6))
summary( subset(out2.c4, chr=1) )
summary( subset(out2.c4, chr=c(7,15)) )

plot(out2.c4)
plot(out2.c4, chr=1, lower="cond-int")
plot(out2.c4, chr=c(6,15), lower="cond-int")
plot(out2.c4, chr=c(7,15), lower="cond-int")

out2sub <- subset(out2, chr=c(1,6,7,15))
plot(out2.c4 - out2sub, allow.neg=TRUE, lower="cond-int")

qc <- c(1, 1, 4, 6, 15)
qp <- c(43.3, 78.3, 30.0, 62.5, 18.0)
qtl <- makeqtl(hyper, chr=qc, pos=qp)

myformula <- y ~ Q1+Q2+Q3+Q4+Q5 + Q4:Q5

out.fq <- fitqtl(hyper, qtl=qtl, formula = myformula)
summary(out.fq)

out.fq <- fitqtl(hyper, qtl=qtl, formula = myformula, drop=FALSE, get.ests=TRUE)
summary(out.fq)

revqtl <- refineqtl(hyper, qtl=qtl, formula = myformula)

revqtl

plot(revqtl)

out.fq2 <- fitqtl(hyper, qtl=revqtl, formula=myformula)
summary(out.fq2)

out1.c4r <- addqtl(hyper, qtl=revqtl, formula=y~Q3)

plot(out1.c4, out1.c4r, col=c("blue", "red"))

plot(out1.c4r - out1.c4, ylim=c(-1.7, 1.7))
abline(h=0, lty=2, col="gray")

out2.c4r <- addpair(hyper, qtl=revqtl, formula=y~Q3, chr=c(1,6,7,15))

plot(out2.c4r - out2.c4, lower="cond-int", allow.neg=TRUE)

out.1more <- addqtl(hyper, qtl=revqtl, formula=myformula)
plot(out.1more)

out.iw4 <- addqtl(hyper, qtl=revqtl, formula=y~Q1+Q2+Q3+Q4+Q5+Q4:Q5+Q6+Q5:Q6)
plot(out.iw4)

out.2more <- addpair(hyper, qtl=revqtl, formula=myformula, chr=c(2,5,7,15))

plot(out.2more, lower="cond-int")

out.ai <- addint(hyper, qtl=revqtl, formula=myformula)
out.ai

qtl2 <- addtoqtl(hyper, revqtl, 7, 53.6)
qtl2

qtl3 <- dropfromqtl(qtl2, index=2)
qtl3

qtl4 <- replaceqtl(hyper, qtl3, index=1, chr=1, pos=50)
qtl4

qtl5 <- reorderqtl(qtl4, c(1:3,5,4))
qtl5

stepout.a <- stepwiseqtl(hyper, additive.only=TRUE, max.qtl=6)
stepout.a

stepout.i <- stepwiseqtl(hyper, max.qtl=6)
stepout.i

############################################################
# Example 6: Internal data structure
############################################################
data(fake.bc)

class(fake.bc)

names(fake.bc)

fake.bc$pheno[1:10,]

names(fake.bc$geno)
sapply(fake.bc$geno, class)

names(fake.bc$geno[[3]])
fake.bc$geno[[3]]$data[1:5,]
fake.bc$geno[[3]]$map

names(fake.bc$geno[[3]])
fake.bc <- calc.genoprob(fake.bc, step=10, err=0.01)
names(fake.bc$geno[[3]])
fake.bc <- sim.geno(fake.bc, step=10, n.draws=8, err=0.01)
names(fake.bc$geno[[3]])
fake.bc <- argmax.geno(fake.bc, step=10, err=0.01)
names(fake.bc$geno[[3]])
fake.bc <- calc.errorlod(fake.bc, err=0.01)
names(fake.bc$geno[[3]])

names(fake.bc)
fake.bc <- est.rf(fake.bc)
names(fake.bc)

# end of rqtltour.R