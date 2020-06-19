
library(qtl)
# library(plotrix)

pheno <- "ul2p" # define the phenotype to use

outpath <- paste("~/git/Slat_SAqtl/output/", pheno, sep="")

# Load first script results
load(paste(outpath, '/',pheno,'.Rdata', sep=""))
outpath <- paste("~/git/Slat_SAqtl/output/", pheno, sep="")

# load function to automatically create auto qtl object based on number of loci
source("~/git/Slat_SAqtl/scripts/aqtl.r")

multisim <- sim.geno(mapthis2, step=1, n.draws=256, err=0.001)
multisim2 <- calc.genoprob(multisim, step=1, error.prob=0.001)

m1lg <- 1
m1cm <- 58.6
m2lg <- 9
m2cm <- 36
m3lg <- 1
m3cm <- 14.2
m4lg <- 4
m4cm <- 46.5
m5lg <- 8
m5cm <- 56.1
m6lg <- 7
m6cm <- 27
m7lg <- 4
m7cm <- 3

# Final model (edit iteratively)
f1 <- "y ~ Sex + Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q1:Q7 + Q3:Q4 + Q5:Q6 + Q2:Q4 + Q2:Q5 + Q3:Q7 + Sex:Q3 + Sex:Q6 "


Sex <- data.frame(Sex=as.numeric(pull.pheno(multisim2, "sex") == "m"))

# Input number of loci for the qtl object
qtl.auto <- autoqtl(7)
fq.auto <- fitqtl(multisim2, qtl=qtl.auto, formula=f1, pheno.col=pheno, cov=Sex, forceXcovar=T)
summary(fq.auto, alpha=0.05, pvalues=TRUE, format="allpheno")

# check for interactions between qtl or with qtl and sex 
qtl.auto.ai <- addint(multisim2, qtl=qtl.auto, formula=f1, pheno.col=pheno, covar=Sex)
summary(qtl.auto.ai, alpha=0.05, pvalues=TRUE, format="allpheno")

# refine positions
rqtl <- refineqtl(multisim2, qtl=qtl.auto, formula=f1, pheno.col=pheno, method=method, covar=Sex, forceXcovar=T )
final.fit <- fitqtl(multisim2, qtl=rqtl, formula=f1, pheno.col=pheno, covar=Sex, forceXcovar=T)
summary(final.fit, alpha=0.05, pvalues=TRUE, format="allpheno")
# the positions are fine

# check interactions again, after the positions have beed refined. This allows to update both the position and the model with new interactions, in the next iteration
qtl.auto.ai3 <- addint(multisim2, qtl=rqtl, formula=f1, pheno.col=pheno, covar=Sex)
summary(qtl.auto.ai3, alpha=0.05, pvalues=TRUE, format="allpheno")

# try to add more qtl
add1q <- addqtl(multisim2, qtl=rqtl, pheno.col=pheno, method=method, covar=Sex, formula=f1)

summary(add1q, perms=operm.a, alpha=0.05, pvalues=TRUE, format="allpheno")

# These are the LOD thresholds used
summary(operm.a)


# The folowing is used to make nice graphical summaries of the QTL and their associated loci.
# plots
pdf(file.path(outpath, paste(pheno,'lodProfile.pdf', sep="_")), width=12, height=6)
par(mfrow=c(1,2)) 
plot(rqtl)
plotLodProfile(rqtl)
dev.off()


# Identify the marker name closest to each QTL, to be used in plotting for quick reference
M1 <- find.marker(mapthis, m1lg, m1cm)
M2 <- find.marker(mapthis, m2lg, m2cm)
M3 <- find.marker(mapthis, m3lg, m3cm)
M4 <- find.marker(mapthis, m4lg, m4cm)
M5 <- find.marker(mapthis, m5lg, m5cm)
M6 <- find.marker(mapthis, m6lg, m6cm)
M7 <- find.marker(mapthis, m7lg, m7cm)


# print the final model, and manually change the interaction between QTL part of the plot
f1

pdf(file.path(outpath, paste(pheno,'effectPlot2.pdf', sep="_")), width=14, height=12)
par(mfrow=c(4,4)) 
par(mar=c(5,5,4,3))
effectplot(multisim2, M1, main=paste("QTL 1 - LG", m1lg),  pheno.col=pheno)
if(grepl("Sex:Q1", f1)) { effectplot(mapthis, mname1="sex", mname2=M1, main="",  pheno.col=pheno) }
effectplot(multisim2, M2, main=paste("QTL 2 - LG", m2lg),  pheno.col=pheno)
if(grepl("Sex:Q2", f1)) { effectplot(mapthis, mname1="sex", mname2=M2, main="",  pheno.col=pheno) }
effectplot(multisim2, M3, main=paste("QTL 3 - LG", m3lg),  pheno.col=pheno)
if(grepl("Sex:Q3", f1)) { effectplot(mapthis, mname1="sex", mname2=M3, main="",  pheno.col=pheno) }
effectplot(multisim2, M4, main=paste("QTL 4 - LG", m4lg),  pheno.col=pheno)
if(grepl("Sex:Q4", f1)) { effectplot(mapthis, mname1="sex", mname2=M4, main="",  pheno.col=pheno) }
effectplot(multisim2, M5, main=paste("QTL 5 - LG", m5lg),  pheno.col=pheno)
if(grepl("Sex:Q5", f1)) { effectplot(mapthis, mname1="sex", mname2=M5, main="",  pheno.col=pheno) }
effectplot(multisim2, M6, main=paste("QTL 6 - LG", m6lg),  pheno.col=pheno)
if(grepl("Sex:Q6", f1)) { effectplot(mapthis, mname1="sex", mname2=M6, main="",  pheno.col=pheno) }
effectplot(multisim2, M7 main=paste("QTL 7 - LG", m7lg),  pheno.col=pheno)
if(grepl("Sex:Q7", f1)) { effectplot(mapthis, mname1="sex", mname2=M7, main="",  pheno.col=pheno) }
effectplot(mapthis, mname1=M1, mname2=M7, main=paste("LG",m1lg," * LG",m7lg, sep=""),  pheno.col=pheno)
effectplot(mapthis, mname1=M3, mname2=M4, main=paste("LG",m3lg," * LG",m4lg, sep=""),  pheno.col=pheno)
effectplot(mapthis, mname1=M3, mname2=M7, main=paste("LG",m3lg," * LG",m7lg, sep=""),  pheno.col=pheno)
effectplot(mapthis, mname1=M5, mname2=M6, main=paste("LG",m5lg," * LG",m6lg, sep=""),  pheno.col=pheno)
effectplot(mapthis, mname1=M2, mname2=M4, main=paste("LG",m2lg," * LG",m4lg, sep=""),  pheno.col=pheno)
effectplot(mapthis, mname1=M2, mname2=M5, main=paste("LG",m2lg," * LG",m5lg, sep=""),  pheno.col=pheno)
dev.off()

sink(file=file.path(outpath, paste(pheno,'PVE_table.txt', sep="_")), split=T)
capture.output(summary(final.fit, alpha=0.05, pvalues=TRUE, format="allpheno"), file=file.path(outpath, paste(pheno,'PVE_table.txt', sep="_")))
sink()

M1
M2
M3
M4
M5
M6
M7

# par(mfrow=c(2,2))
# effectplot(mapthis, mname1="sex", mname2="100125767", main="Chromosome 3", add.legend=FALSE, pheno.col=pheno)
# effectplot(mapthis, mname1="sex", ylim=c(15.1, 17.2), mname2="X@58", main="X chromosome", add.legend=FALSE)



