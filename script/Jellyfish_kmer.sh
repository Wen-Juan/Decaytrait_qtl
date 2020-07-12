####run kmer analysis from jellyfish to estimate genome size
/usr/local/bio/jellyfish/jellyfish count -m 31 -o KG_km31_fastq.counts -C -s 500000000 -t 5 KG_TCAGTT_R1_pairedR1.fastq KG_TCAGTT_R2_pairedR2.fastq

###generate histogram file
 /usr/local/bio/jellyfish/jellyfish histo KG_km31_fastq.counts > KG_km31_fastq.counts.histo

###using genomescope v2.0 was recommended to use kmer 21 for most of genomes
###kmer 21
http://genomescope.org/genomescope2.0/analysis.php?code=pI28nXUlI0povpWlov7F

###kmer31
http://genomescope.org/genomescope2.0/analysis.php?code=wzOxPMMnlDutfLt3D7ys

###Kmer25
http://genomescope.org/genomescope2.0/analysis.php?code=93S8LzwWZNY9CJ99L05e

###kmer41
http://genomescope.org/genomescope2.0/analysis.php?code=T52ri5URhNsSLSYn7csY

###kmer21 AM population
http://genomescope.org/genomescope2.0/analysis.php?code=ULemyk0pi0OnvGbcEFuv
