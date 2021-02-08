#install R packages and load R libraries
library(ggplot2)
library(plyr)
install.packages("devtools")
library(devtools)
install_github("kassambara/easyGgplot2", force = TRUE)
library(easyGgplot2)
install.packages("gridExtra")
library(gridExtra)
require(cowplot)
library(Hmisc)
library(ggplot2)
library(gridExtra)
install.packages("tidyverse")
library(tidyverse)


#load the corresponding data files, mapping to 1064 a1 strain. 
am_35F <- read.table('~/AM_000035F_sorted_q20_variants_filtered_DP.vcf', sep="\t",header=F,blank.lines.skip=TRUE,
                                  comment.char = "#")
str(am_35F)

am_58F <- read.table('~/AM_000058F_sorted_q20_variants_filtered_DP.vcf', sep="\t",header=F,blank.lines.skip=TRUE,
                     comment.char = "#")
str(am_58F)

kg_35F <- read.table('~/KG_000035F_sorted_q20_variants_filtered_DP.vcf', sep="\t",header=F,blank.lines.skip=TRUE,
                     comment.char = "#")
str(kg_35F)

kg_58F <- read.table('~/KG_000058F_sorted_q20_variants_filtered_DP.vcf', sep="\t",header=F,blank.lines.skip=TRUE,
                     comment.char = "#")
str(kg_58F)

colnames(am_35F)<-c("chr","start","id","refallele","altallele","qual",
                                 "filter","info","format")
summary(am_35F)


colnames(am_58F)<-c("chr","start","id","refallele","altallele","qual",
                    "filter","info","format")
summary(am_58F)

colnames(kg_35F)<-c("chr","start","id","refallele","altallele","qual",
                    "filter","info","format")
summary(kg_35F)

colnames(kg_58F)<-c("chr","start","id","refallele","altallele","qual",
                    "filter","info","format")
summary(kg_58F)

###snps on promoter region and exon regions, 
ggplot(am_35F) + 
  geom_histogram(aes(x=start),binwidth=100, col= "dark grey") +  # pick a binwidth that is not too small 
  ggtitle("SNPs on homonal receptor 4 gene") +
  xlab("Position on contig") +
  ylim(0,10) +
  xlim(1506000,1535000) +
  ylab("SNP density") +
  theme_bw() + 
  geom_vline(xintercept =1506170,colour = "orange", size = 1,alpha=0.6) +
  geom_vline(xintercept =1511170,colour = "orange", size = 1,alpha=0.6) +
  geom_vline(xintercept =1511507,colour = "green", size =1,alpha=0.6) +
  geom_vline(xintercept =1519632,colour = "green", size =1,alpha=0.6) +
  geom_vline(xintercept =1522821,colour = "green", size =1,alpha=0.6) +
  geom_vline(xintercept =1523514,colour = "green", size =1,alpha=0.6) +
  geom_vline(xintercept =1523792,colour = "green", size =1,alpha=0.6) +
  geom_vline(xintercept =1524332,colour = "green", size = 1,alpha=0.6) +
  geom_vline(xintercept =1528355,colour = "green", size = 1,alpha=0.6) +
  geom_vline(xintercept =1529063,colour = "green", size = 1,alpha=0.6) +
  geom_vline(xintercept =1530801,colour = "green", size = 1,alpha=0.6) +
  geom_vline(xintercept =1531627,colour = "green", size = 1,alpha=0.6) +
  geom_vline(xintercept =1532586,colour = "green", size = 1,alpha=0.6) +
  geom_vline(xintercept =1533266,colour = "green", size = 1,alpha=0.6) +
  theme(axis.title.x = element_text(size=14,colour = "black"),axis.title.y = element_text(size=14,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12)) +
  theme(axis.text.y = element_text(colour="black",size=12)) +
  theme_bw()

pdf(file="~/Decaytrait_qtl/output/SNPdensity_AM35F_promoter_exons.pdf",width=10,height=10)
dev.off()

###SNP density on contig 000035F.
ggplot(am_35F) + 
  geom_histogram(aes(x=start),binwidth=1, col= "blue") +  # pick a binwidth that is not too small 
  ggtitle("SNPs on 000035F AM") +
  xlab("Position on contig") +
  ylim(0,5) +
  xlim(1510000,1541000) +
  ylab("SNP density") +
  theme_bw() + 
 # geom_vline(xintercept =1510000,colour = "yellow", size = 3,alpha=0.7) +
#  geom_vline(xintercept =1545000,colour = "yellow", size = 3,alpha=0.7) +
  theme(axis.title.x = element_text(size=14,colour = "black"),axis.title.y = element_text(size=14,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12)) +
  theme(axis.text.y = element_text(colour="black",size=12)) +
  theme_bw()
pdf(file="~/Decaytrait_qtl/output/SNPdensity_35F.pdf",width=10,height=16)
dev.off()


###SNP density on contig 000058F.
  ggplot(am_58F) + 
  geom_histogram(aes(x=start),binwidth=500, col= "blue") +  # pick a binwidth that is not too small 
  ggtitle("SNPs on 000058F AM") +
  xlab("Position on contig") +
  ylim(0,40) +
  xlim(0,1780000) +
  ylab("SNP density") +
  theme_bw() + 
  theme(axis.title.x = element_text(size=14,colour = "black"),axis.title.y = element_text(size=14,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12)) +
  theme(axis.text.y = element_text(colour="black",size=12)) +
  theme_bw()
pdf(file="~/Decaytrait_qtl/output/SNPdensity_58F.pdf",width=10,height=16)
dev.off()

