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

###

Diachasma <- read.table("/Users/Wen-Juan/Downloads/Diachasma_list1.txt", header = FALSE)
str(Diachasma)

Fopis <- read.table("/Users/Wen-Juan/Downloads/Fopis_list2.txt", header = FALSE)
str(Fopis)

Microplitis <- read.table("/Users/Wen-Juan/Downloads/Microplitis_list1.txt", header = FALSE)
str(Microplitis)

Ajap <-  read.table("/Users/Wen-Juan/Downloads/annotation_list_twocontigs_Ajap.txt", header = FALSE)
str(Ajap)


#for shared female-biased genes among five stages
venn.plot <- venn.diagram(list(D.alloeum = as.character(Diachasma$V1), F.arisanus = as.character(Fopis$V1), M.demolitor = as.character(Microplitis$V1), A.japonica = as.character(Ajap$V1)), filename =NULL,
                          fill=c("orange","blue","green","red"),
                          ext.line.lwd = 3,
                          cex = 2,
                          cat.cex = 1.5,
                          rotation.degree = 0)

# to draw to the screen:
grid.arrange(gTree(children=venn.plot),ncol = 1 )

# to output to pdf
venn_fbias_out <- arrangeGrob(gTree(children=venn.plot),ncol = 1 )
ggsave(file="total_annotation_4species.pdf", venn_fbias_out, path = "/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/")


##get overlap list
shared_anno <- list(D.alloeum = as.character(Diachasma$V1), F.arisanus = as.character(Fopis$V1), M.demolitor = as.character(Microplitis$V1), A.japonica = as.character(Ajap$V1))

Intersect <- function (x) {  
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], Intersect(x[-1]))
  }
}

Union <- function (x) {  
  # Multiple set version of union
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], Union(x[-1]))
  }
}

Setdiff <- function (x, y) {
  # Remove the union of the y's from the common x's. 
  # x and y are lists of characters.
  xx <- Intersect(x)
  yy <- Union(y)
  setdiff(xx, yy)
}

combs <- 
  unlist(lapply(1:length(shared_anno), 
                function(j) combn(names(shared_anno), j, simplify = FALSE)),
         recursive = FALSE)
names(combs) <- sapply(combs, function(i) paste0(i, collapse = ""))
str(combs)

elements <- 
  lapply(combs, function(i) Setdiff(shared_anno[i], shared_anno[setdiff(names(shared_anno), i)]))

n.elements <- sapply(elements, length)
print(n.elements)

ItemsList <- venn(shared_anno, show.plot = FALSE)

overlap <- calculate.overlap(list)

overlap
