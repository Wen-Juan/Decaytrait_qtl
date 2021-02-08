library(qtl)

outpath <- '~/Decaytrait_qtl/output/'

#operm.em
multisim <- sim.geno(qtl, step=1, n.draws=10000, err=0.001)
multisim2 <- calc.genoprob(multisim, step=1, error.prob=0.001)

m1lg <- 12
m1cm <- 2
m2lg <- 12
m2cm <- 23

Cross <- data.frame(Cross=as.numeric(pull.pheno(qtl, "cross") == "F2"))
f1 <- "y ~ Cross + Q1 + Q2 + Cross:Q1 + Cross:Q2 + Q1:Q2"

# Input number of loci for the qtl object
qtl.auto <- makeqtl(multisim2, chr=c(m1lg, m2lg), pos=c(m1cm, m2cm), what="draws")
fq.auto <- fitqtl(multisim2, cov=Cross,qtl=qtl.auto, formula=f1,model='binary')
summary(fq.auto, alpha=0.05, pvalues=TRUE, format="allpheno")

# check for interactions between qtl or with qtl and cross
qtl.auto.ai <- addint(multisim2, qtl=qtl.auto, formula=f1, pheno.col=1, covar=Cross)
summary(qtl.auto.ai, alpha=0.05, pvalues=TRUE, format="allpheno")

# refine positions
rqtl <- refineqtl(multisim2, qtl=qtl.auto, formula=f1, pheno.col=1, method=c("imp"), covar=Cross,model='binary')
final.fit <- fitqtl(multisim2, qtl=rqtl, formula=f1, pheno.col=1, covar=Cross, model='binary')
summary(final.fit, alpha=0.05, pvalues=TRUE, format="allpheno")
# the positions are fine
###

# check interactions again, after the positions have beed refined. This allows to update both the position and the model with new interactions, in the next iteration
qtl.auto.ai3 <- addint(multisim, qtl=rqtl, formula=f1, pheno.col=1, covar=Cross)
summary(qtl.auto.ai3, alpha=0.05, pvalues=TRUE, format="allpheno")
#p=0.198, so no interaction.

# try to add more qtl
add1q <- addqtl(multisim2, qtl=rqtl, formula=f1, pheno.col=1, method=c("imp"),model='binary', covar=Cross)
summary(add1q, perms=operm.hk, alpha=0.05, pvalues=TRUE, format="allpheno")
#summary(add1q, perms=operm.em, alpha=0.05, pvalues=TRUE, format="allpheno")
# These are the LOD thresholds used
summary(operm.hk)
#summary(operm.em)

# plots
pdf("~/Decaytrait_qtl/output/qtl_lod_2qtl_lodprofile.pdf",  width=12, height=6)
par(mfrow=c(1,2)) 
plot(rqtl)
plotLodProfile(rqtl)
dev.off()

# Identify the marker name closest to each QTL, to be used in plotting for quick reference
M1 <- find.marker(qtl, m1lg, m1cm)
M2 <- find.marker(qtl, m2lg, m2cm)

# print the final model, and manually change the interaction between QTL part of the plot
f1

pdf("~/Decaytrait_qtl/output/qtl_lod_2qtl_effectplot.pdf",  width=12, height=6)
par(mfrow=c(1,2)) 
par(mar=c(5,5,4,3))
effectplot(multisim2, M1, main=paste("QTL 1 - LG", m1lg),  pheno.col=1)
effectplot(multisim2, M2, main=paste("QTL 2 - LG", m2lg),  pheno.col=1)
# effectplot(multisim, mname1=M1, mname2=M3, main=paste("LG",m1lg," * LG",m3lg, sep=""))
dev.off()

capture.output(summary(final.fit, alpha=0.05, pvalues=TRUE, format="allpheno"), file=file.path(outpath, paste("pheno",'2qtl_PVE_table.txt', sep="_")))

M1
M2



