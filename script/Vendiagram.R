### venns in R
install.packages("VennDiagram")

library("VennDiagram")
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

#for the shared female-biased genes across five stages.
setwd <-'~/Tvedora_dev_RNAseq/'
datapath <- '~/Tvedora_dev_RNAseq/input/sex_bias/'
results <- '~/Tvedora_dev_RNAseq/input/sex_bias/'

fivestages_fbias <- read.table("~/Tvedora_dev_RNAseq/input/sexbias_female_inclusr.txt", header = TRUE)
str(fivestages_fbias)

g43XY_fbias <- read.table("~/Tvedora_dev_RNAseq/input/XY_fbias.txt", header = TRUE)

fbis_g23 <- fivestages_fbias$trans[fivestages_fbias$stage=='G23']
fbis_g27 <- fivestages_fbias$trans[fivestages_fbias$stage=='G27']
fbis_g31 <- fivestages_fbias$trans[fivestages_fbias$stage=='G31']
fbis_g43 <- fivestages_fbias$trans[fivestages_fbias$stage=='G43']
fbis_g46 <- fivestages_fbias$trans[fivestages_fbias$stage=='G46']
fbis_g43xy <- fivestages_fbias$trans[fivestages_fbias$stage=='G43XY']
fbis_g46xx <- fivestages_fbias$trans[fivestages_fbias$stage=='G46XX']

#for shared female-biased genes among five stages
venn.plot <- venn.diagram(list(G23 = as.character(fbis_g23), G27 = as.character(fbis_g27), G31 = as.character(fbis_g31), G43=as.character(fbis_g43), G46=as.character(fbis_g46)), filename =NULL,
                          fill=c("orange","red","black","green","blue"),
                          ext.line.lwd = 3,
                          cex = 2,
                          cat.cex = 2.5,
                          rotation.degree = 60)

# to draw to the screen:
grid.arrange(gTree(children=venn.plot),ncol = 1 )

# to output to pdf
venn_fbias_out <- arrangeGrob(gTree(children=venn.plot),ncol = 1 )
ggsave(file="shared_fembias.pdf", venn_fbias_out, path = "~/Tvedora_dev_RNAseq/output/")

#overlap of G43 female-biased genes resulting from comparison between XX females, or XY0 females and XY0 males.
venn.plot1 <- venn.diagram(list(XX_f=as.character(fbis_g43), XY0_f=as.character(fbis_g43xy)), filename =NULL,
                           fill=c("orange","red"),
                           ext.line.lwd = 2,
                           cex = 1,
                           cat.cex = 1,rotation.degree = 60)

grid.arrange(gTree(children=venn.plot1),ncol = 1 )
venn_fbias_out <- arrangeGrob(gTree(children=venn.plot1),ncol = 1 )
ggsave(file="g43_XYXXfemale_fembias.pdf", venn_fbias_out, path = "~/Tvedora_dev_RNAseq/")

##overlap of G46 female-biased genes resulting from comparison between XX females and either XY0 males or XX male.
venn.plot2 <- venn.diagram(list(XY0_m=as.character(fbis_g46), XX_m=as.character(fbis_g46xx)), filename =NULL,
                           fill=c("orange","red"),
                           ext.line.lwd = 2,
                           cex = 1,
                           cat.cex = 1,rotation.degree = 60)

grid.arrange(gTree(children=venn.plot2),ncol = 1 )
venn_fbias_out <- arrangeGrob(gTree(children=venn.plot2),ncol = 1 )
ggsave(file="g46_XXmale_fembias.pdf", venn_fbias_out, path = "~/Tvedora_dev_RNAseq/output/")


#test whether the overlap of female-biased genes across five stages is greater than by chance.
f_bias <- list(as.character(fbis_g23), as.character(fbis_g27), as.character(fbis_g31),as.character(fbis_g43),as.character(fbis_g46))
list(f_bias)
str(f_bias)

total1 <- 9680

####### TO DO for all interections

res=supertest(f_bias, n=total1)
plot(res, sort.by="size")
plot(res, Layout="landscape", degree=2:4, sort.by="size")
summary(res)

write.csv(summary(res)$Table, file="~/Tvedora_dev_RNAseq/output/f_bias_testbychance.csv", row.names=FALSE)


#for shared male-biased genes
mbias <- read.table("~/Tvedora_dev_RNAseq/input/share_mbias_incluesr.txt", header = TRUE)
str(mbias)


mbis_g27 <- mbias$trans[mbias$stage=='G27']
mbis_g31 <- mbias$trans[mbias$stage=='G31']
mbis_g43 <- mbias$trans[mbias$stage=='G43']
mbis_g46 <- mbias$trans[mbias$stage=='G46']
mbis_g46xx <- mbias$trans[mbias$stage=='G46XX']


venn.plot.2 <- venn.diagram(list(G27 = as.character(mbis_g27), G31 = as.character(mbis_g31), G43=as.character(mbis_g43), G46=as.character(mbis_g46)), filename =NULL,
                            fill=c("green","blue"),
                            ext.line.lwd = 3,
                            cex = 2,
                            cat.cex = 2.5)

## to draw to the screen:
grid.arrange(gTree(children=venn.plot.2),ncol = 1 )

# to output to pdf
venn_mbias_out <- arrangeGrob(gTree(children=venn.plot.2),ncol = 1 )
ggsave(file="shared_mbias.pdf", venn_mbias_out, path = "~/Tvedora_dev_RNAseq/output/")

#overlap of male-biased genes resulting from the comparison between XX females and either XY0 males, or XX male.
venn.plot.3 <- venn.diagram(list(G46=as.character(mbis_g46), G46XX=as.character(mbis_g46xx)), filename =NULL,
                            fill=c("green","blue"),
                            ext.line.lwd = 2,
                            cex = 1,
                            cat.cex = 1,rotation.degree = 60)

grid.arrange(gTree(children=venn.plot.3),ncol = 1 )
venn_mbias_out <- arrangeGrob(gTree(children=venn.plot.3),ncol = 1 )
ggsave(file="g46_xxmale_mbias.pdf", venn_mbias_out, path = "~/Tvedora_dev_RNAseq/output/")


#test whether the overlap of male-biased genes across five stages is greater than by chance.
m_bias <- list(as.character(mbis_g27), as.character(mbis_g31),as.character(mbis_g43),as.character(mbis_g46))
list(m_bias)
str(m_bias)
