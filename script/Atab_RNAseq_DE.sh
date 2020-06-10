      ####Step1. trim the transcriptome a bit#####
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
####2.1 index the trascriptome
/usr/local/bio/kallisto/kallisto index -i Atab_trinity.idx Atab_trinity_300bp.fa

####2.2 run kallisto to quantify the expression
