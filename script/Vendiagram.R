### venns in R
install.packages("VennDiagram")

library("VennDiagram")
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)


shared_genes <- read.table("/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/input/shared_genes.txt", header = TRUE)
str(shared_genes)

alloeum <- shared_genes$gene[shared_genes$species=='D._alloeum']
arisanus <- shared_genes$gene[shared_genes$species=='F._arisanus']
demolitor <- shared_genes$gene[shared_genes$species=='M._demolitor']

#for shared female-biased genes among five stages
venn.plot <- venn.diagram(list(D.alloeum = as.character(alloeum), F.arisanus = as.character(arisanus), M.demolitor = as.character(demolitor)), filename =NULL,
                          fill=c("orange","blue","green"),
                          ext.line.lwd = 3,
                          cex = 2,
                          cat.cex = 1.5,
                          rotation.degree = 0)

# to draw to the screen:
grid.arrange(gTree(children=venn.plot),ncol = 1 )

# to output to pdf
venn_fbias_out <- arrangeGrob(gTree(children=venn.plot),ncol = 1 )
ggsave(file="shared_genes.pdf", venn_fbias_out, path = "/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/")


