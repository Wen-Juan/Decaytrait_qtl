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
am_35F <- read.table('/Users/Wen-Juan/Downloads/AM_000035F_sorted_q20_variants_filtered_DP.vcf', sep="\t",header=F,blank.lines.skip=TRUE,
                                  comment.char = "#")
str(am_35F)

am_58F <- read.table('/Users/Wen-Juan/Downloads/AM_000058F_sorted_q20_variants_filtered_DP.vcf', sep="\t",header=F,blank.lines.skip=TRUE,
                     comment.char = "#")
str(am_58F)

kg_35F <- read.table('/Users/Wen-Juan/Downloads/KG_000035F_sorted_q20_variants_filtered_DP.vcf', sep="\t",header=F,blank.lines.skip=TRUE,
                     comment.char = "#")
str(kg_35F)

kg_58F <- read.table('/Users/Wen-Juan/Downloads/KG_000058F_sorted_q20_variants_filtered_DP.vcf', sep="\t",header=F,blank.lines.skip=TRUE,
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

y1 <-
  ggplot(am_35F) + 
  geom_histogram(aes(x=start),binwidth=500, col= "blue") +  # pick a binwidth that is not too small 
  ggtitle("SNPs on 000035F AM") +
  xlab("Position on contig") +
  ylim(0,40) +
  xlim(0,2500000) +
  ylab("SNP density") +
  theme_bw() + 
  geom_vline(xintercept =1510000,colour = "yellow", size = 3,alpha=0.7) +
  geom_vline(xintercept =1545000,colour = "yellow", size = 3,alpha=0.7) +
  theme(axis.title.x = element_text(size=14,colour = "black"),axis.title.y = element_text(size=14,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12)) +
  theme(axis.text.y = element_text(colour="black",size=12)) +
  theme_bw()

y2 <-
  ggplot(kg_35F) + 
  geom_histogram(aes(x=start),binwidth=500, col= "blue") +  # pick a binwidth that is not too small 
  ggtitle("SNPs on 000035F KG") +
  xlab("Position on contig") +
  ylim(0,40) +
  xlim(0,2500000) +
  ylab("SNP density") +
  theme_bw() + 
  geom_vline(xintercept =1510000,colour = "yellow", size = 3,alpha=0.7) +
  geom_vline(xintercept =1545000,colour = "yellow", size = 3,alpha=0.7) +
  theme(axis.title.x = element_text(size=14,colour = "black"),axis.title.y = element_text(size=14,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12)) +
  theme(axis.text.y = element_text(colour="black",size=12)) +
  theme_bw()

pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/SNPdensity_35F_withcandidategene.pdf",width=10,height=16)
grid.arrange(y1,y2, nrow = 2)
dev.off()

y3 <-
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


y4 <-
  ggplot(kg_58F) + 
  geom_histogram(aes(x=start),binwidth=500, col= "blue") +  # pick a binwidth that is not too small 
  ggtitle("SNPs on 000058F KG") +
  xlab("Position on contig") +
  ylim(0,40) +
  xlim(0,1780000) +
  ylab("SNP density") +
  theme_bw() + 
  theme(axis.title.x = element_text(size=14,colour = "black"),axis.title.y = element_text(size=14,colour = "black")) +
  theme(axis.text.x = element_text(colour="black",size=12)) +
  theme(axis.text.y = element_text(colour="black",size=12)) +
  theme_bw()

pdf(file="/Users/Wen-Juan/Dropbox (Amherst College)/Amherst_postdoc/github/Decaytrait_qtl/output/SNPdensity_58F.pdf",width=10,height=16)
grid.arrange(y3,y4, nrow = 2)
dev.off()
