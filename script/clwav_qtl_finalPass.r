library(qtl)

outpath <- '/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/'

operm.hk
#operm.em
multisim <- sim.geno(qtl, step=1, n.draws=1000, err=0.001)
# multisim <- calc.genoprob(multisim, step=1, error.prob=0.01)

m1lg <- 12
m1cm <- 3
m2lg <- 12
m2cm <- 21


f1 <- "y ~ Q1 + Q2"

# Input number of loci for the qtl object

qtl.auto <- makeqtl(multisim, chr=c(m1lg, m2lg), pos=c(m1cm, m2cm), what="draws")
fq.auto <- fitqtl(multisim, covar=NULL,qtl=qtl.auto, formula=f1,model='binary')
summary(fq.auto, alpha=0.05, pvalues=TRUE, format="allpheno")

# check for interactions between qtl
qtl.auto.ai <- addint(multisim, qtl=qtl.auto, formula=f1)
summary(qtl.auto.ai, alpha=0.05, pvalues=TRUE, format="allpheno") 
##p=0.21, so there is no significant interaction.

# refine positions
rqtl <- refineqtl(multisim, qtl=qtl.auto, formula=f1, verbose=FALSE) 
final.fit <- fitqtl(multisim, qtl=rqtl, formula=f1, model='binary')
summary(final.fit, alpha=0.05, pvalues=TRUE, format="allpheno")

# check interactions again, after the positions have beed refined. This allows to update both the position and the model with new interactions, in the next iteration
qtl.auto.ai3 <- addint(multisim, qtl=rqtl, formula=f1, pheno.col=1)
summary(qtl.auto.ai3, alpha=0.05, pvalues=TRUE, format="allpheno")
#p=0.198, so no interaction.

# try to add more qtl
add1q <- addqtl(multisim, qtl=rqtl, formula=f1, model='binary')

summary(add1q, perms=operm.hk, alpha=0.05, pvalues=TRUE, format="allpheno")
#summary(add1q, perms=operm.em, alpha=0.05, pvalues=TRUE, format="allpheno")
# These are the LOD thresholds used
summary(operm.hk)
#summary(operm.em)

# plots
pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/qtl_lod_hk90kperm_lodprofile.pdf",  width=12, height=6)
par(mfrow=c(1,2)) 
plot(rqtl)
plotLodProfile(rqtl)
dev.off()

# Identify the marker name closest to each QTL, to be used in plotting for quick reference
M1 <- find.marker(qtl, m1lg, m1cm)
M2 <- find.marker(qtl, m2lg, m2cm)
#M3 <- find.marker(qtl, 5, 0)

# print the final model, and manually change the interaction between QTL part of the plot
f1

pdf("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/qtl_lod_hk_perm90k_effectplot.pdf",  width=12, height=6)
par(mfrow=c(1,2)) 
par(mar=c(5,5,4,3))
effectplot(multisim, M1, main=paste("QTL 1 - LG", m1lg),  pheno.col=1)
effectplot(multisim, M2, main=paste("QTL 2 - LG", m2lg),  pheno.col=1)
# effectplot(multisim, mname1=M1, mname2=M3, main=paste("LG",m1lg," * LG",m3lg, sep=""))
dev.off()

capture.output(summary(final.fit, alpha=0.05, pvalues=TRUE, format="allpheno"), file=file.path(outpath, paste("pheno",'hk90kperm_PVE_table.txt', sep="_")))

M1
M2



