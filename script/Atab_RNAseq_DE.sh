####Step1. trim the transcriptome a bit#####
################################################
################################################

###trim the trinity generated de novo transcriptome.
perl remove_short_transcript.sql 300 TCAGTT_Trinity.fasta > Atab_trinity_300bp.fa

###check the gene and transcript number using Trinity perl script
/usr/local/bio/trinityrnaseq/util/TrinityStats.pl Atab_trinity_300bp.fa

###################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':	21965
Total trinity transcripts:	26751
Percent GC: 44.58
###################################


####Step2. Quantify accounts of transcripts#####
################################################
################################################

####2.1 index the trascriptome
/usr/local/bio/kallisto/kallisto index -i Atab_trinity.idx Atab_trinity_300bp.fa

####2.2 run kallisto to quantify the expression
/usr/local/bio/kallisto/kallisto quant -i Atab_trinity.idx -o ./output/TCAGTT/ --single -l 350 -s 50 TCAGTT.fq -t 5
/usr/local/bio/kallisto/kallisto quant -i Atab_trinity.idx -o ./output/AGAGAT/ --single -l 350 -s 50 AGAGAT.fq -t 5

####2.3 build the matrix of transcript abundance
/usr/local/bio/kallisto/kallisto
/usr/local/bio/trinityrnaseq/util/R/

/usr/local/bio/trinityrnaseq/util/abundance_estimates_to_matrix.pl \
--est_method kallisto --out_prefix Atab_ \
--gene_trans_map none \
--name_sample_by_basedir \
output/AGAGAT/AGAGAT.tsv \
output/TCAGTT/TCAGTT.tsv

#####################
####################
