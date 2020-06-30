####run kmer analysis from jellyfish to estimate genome size
/usr/local/bio/jellyfish/jellyfish count -m 31 -o KG_km31_fastq.counts -C -s 500000000 -t 5 KG_TCAGTT_R1_pairedR1.fastq KG_TCAGTT_R2_pairedR2.fastq

###generate histogram file
 /usr/local/bio/jellyfish/jellyfish histo KG_km31_fastq.counts > KG_km31_fastq.counts.histo
