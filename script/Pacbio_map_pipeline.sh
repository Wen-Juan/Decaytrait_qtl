#############################################################Step 1 i######################################
###########################################################################################################
###########################################################################################################
####################################################some rough stats about the sequencing coverage.
###########
###########
##paired-end 100bp sequencing

#Calculate the number of raw paired-end reads number.
zcat my.fastq.gz | echo $((`wc -l`/4)) 

#AM: 52250806*2=104501612, average coverage 42, using genome size 250Mb.
#KG:51756653*2=103513306, average coverage 42.


#############################################################Step 2 i######################################
###########################################################################################################
###########################################################################################################
trimming reads to remove adapter and lower sequening score, using Trimmomatic for paired end reads.
###########
###########
###using Amherst cluster, loading trimmomatic
/usr/local/bio/triommatic/triommatic 

###running the code for AM reads
/usr/local/bio/triommatic/triommatic
ADAPTERS="/usr/local/bio/triommatic/adapters/*.fa"

/usr/local/bio/triommatic/triommatic PE -phred33 -threads 5 AM_AGAGAT_R1.fq.gz AM_AGAGAT_R2.fq.gz \
AM_AGAGAT_R1_pairedR1.fastq AM_AGAGAT_R1_unpairedR1.fastq \
AM_AGAGAT_R2_pairedR2.fastq AM_AGAGAT_R2_unpairedR2.fastq \
ILLUMINACLIP:$ADAPTERS:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:36 &> AM_AGAGAT_R1_trim.log
###filtered AM paired read number: 52250806 Both Surviving: 51961649 (99.45%) 

### KG reads
/usr/local/bio/triommatic/triommatic
ADAPTERS="/usr/local/bio/triommatic/adapters/*.fa"

/usr/local/bio/triommatic/triommatic PE -phred33 -threads 5 KG_TCAGTT_R1.fq.gz KG_TCAGTT_R2.fq.gz \
KG_TCAGTT_R1_pairedR1.fastq KG_TCAGTT_R1_unpairedR1.fastq \
KG_TCAGTT_R2_pairedR2.fastq KG_TCAGTT_R2_unpairedR2.fastq \
ILLUMINACLIP:$ADAPTERS:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:36 &> KG_TCAGTT_R1_trim.log
###filtered KG paired read number: 51756653 Both Surviving: 51470421 (99.45%)

#######fastqc for trimmed paired reads
/usr/local/bio/fastqc/fastqc -o ./fastqc/ KG_TCAGTT_R1_pairedR1.fastq KG_TCAGTT_R2_pairedR2.fastq -t 5
/usr/local/bio/fastqc/fastqc -o ./fastqc/ AM_AGAGAT_R1_pairedR1.fastq AM_AGAGAT_R2_pairedR2.fastq -t 5

#############################################################Step 3 i######################################
###########################################################################################################
###########################################################################################################
####mapping the filtered paired-end reads against reference genome.

###creat genome index
/usr/local/bio/bwa.kit/bwa index Ajap_PacBio.fa

###mapping the paired end reads
##############BWA mem mapping
###mapping the paired end reads
/usr/local/bio/bwa.kit/bwa mem -M -t 16 ../Ajap_PacBio.fa  ../AM_AGAGAT_R1_pairedR1.fastq ../AM_AGAGAT_R2_pairedR2.fastq > AM_AGAGAT_paired.sam

/usr/local/bio/bwa.kit/bwa mem -M -t 16 ../Ajap_PacBio.fa  ../KG_TCAGTT_R1_pairedR1.fastq ../KG_TCAGTT_R2_pairedR2.fastq > KG_TCAGTT_paired.sam


####convert SAM file to BAM format
/usr/local/bio/bwa.kit/samtools view -bS AM_AGAGAT_paired.sam > AM_AGAGAT_paired.bam
/usr/local/bio/bwa.kit/samtools view -bS KG_TCAGTT_paired.sam > KG_TCAGTT_paired.bam

###sorting bam file
/usr/local/bio/bwa.kit/samtools sort -o AM_AGAGAT_paired_sorted.bam -T temporary1_sort.bam AM_AGAGAT_paired.bam
/usr/local/bio/bwa.kit/samtools sort -o KG_TCAGTT_paired_sorted.bam -T temporary2_sort.bam KG_TCAGTT_paired.bam

###index sorted bam file
/usr/local/bio/bwa.kit/samtools index AM_AGAGAT_paired_sorted.bam
/usr/local/bio/bwa.kit/samtools index KG_TCAGTT_paired_sorted.bam

##Count with flagstat for additional information:
/usr/local/bio/bwa.kit/samtools flagstat AM_AGAGAT_paired_sorted.bam > AM_AGAGAT_paired_sorted_stats.txt
#####
105124193 + 0 in total (QC-passed reads + QC-failed reads)
1200895 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
100539997 + 0 mapped (95.64% : N/A)
103923298 + 0 paired in sequencing
51961649 + 0 read1
51961649 + 0 read2
92220108 + 0 properly paired (88.74% : N/A)
98851700 + 0 with itself and mate mapped
487402 + 0 singletons (0.47% : N/A)
4209608 + 0 with mate mapped to a different chr
1849560 + 0 with mate mapped to a different chr (mapQ>=5)

######
/usr/local/bio/bwa.kit/samtools flagstat KG_TCAGTT_paired_sorted.bam >  KG_TCAGTT_paired_sorted_stats.txt
103481815 + 0 in total (QC-passed reads + QC-failed reads)
540973 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
100247694 + 0 mapped (96.87% : N/A)
102940842 + 0 paired in sequencing
51470421 + 0 read1
51470421 + 0 read2
96357952 + 0 properly paired (93.61% : N/A)
99369542 + 0 with itself and mate mapped
337179 + 0 singletons (0.33% : N/A)
2157956 + 0 with mate mapped to a different chr
870706 + 0 with mate mapped to a different chr (mapQ>=5)


###filtering the low quality of mapping reads
/usr/local/bio/bwa.kit/samtools view -q 20 -b KG_TCAGTT_paired_sorted.bam > KG_TCAGTT_paired_sorted_q20.bam
/usr/local/bio/bwa.kit/samtools view -q 20 -b  AM_AGAGAT_paired_sorted.bam > AM_AGAGAT_paired_sorted_q20.bam

#############################################################Step 4 i######################################
###########################################################################################################
###########################################################################################################
###substract a certain region 
/usr/local/bio/bwa.kit/samtools view -b -h AM_AGAGAT_paired_sorted_q20.bam "000035F" > AM_000035F_sorted_q20.bam
/usr/local/bio/bwa.kit/samtools view -b -h AM_AGAGAT_paired_sorted_q20.bam "000058F" > AM_000058F_sorted_q20.bam

/usr/local/bio/bwa.kit/samtools index AM_000035F_sorted_q20.bam
/usr/local/bio/bwa.kit/samtools index AM_000058F_sorted_q20.bam

####calling SNPs from BAM file
####Index the genome assembly (again!)
/usr/local/bio/bwa.kit/samtools faidx Ajap_PacBio.fa
/usr/local/bio/bwa.kit/samtools mpileup -g -f Ajap_PacBio.fa KG_TCAGTT_paired_sorted_q20.bam > KG_TCAGTT_paired_sorted_q20.bcf
/usr/local/bio/bwa.kit/samtools mpileup -g -f Ajap_PacBio.fa AM_AGAGAT_paired_sorted_q20.bam > AM_AGAGAT_paired_sorted_q20.bcf

/usr/local/bio/bwa.kit/samtools mpileup -g -f Ajap_PacBio.fa bwa_mapping/AM_000035F_sorted_q20.bam  > bwa_mapping/AM_000035F_sorted_q20.bcf
/usr/local/bio/bwa.kit/samtools mpileup -g -f Ajap_PacBio.fa bwa_mapping/AM_000058F_sorted_q20.bam  > bwa_mapping/AM_000058F_sorted_q20.bcf
/usr/local/bio/bwa.kit/samtools mpileup -g -f Ajap_PacBio.fa bwa_mapping/KG_000058F_sorted_q20.bam  > bwa_mapping/KG_000058F_sorted_q20.bcf
/usr/local/bio/bwa.kit/samtools mpileup -g -f Ajap_PacBio.fa bwa_mapping/KG_000035F_sorted_q20.bam  > bwa_mapping/KG_000035F_sorted_q20.bcf

###using bcftools
##second, Detect the single nucleotide polymorphisms (SNPs)
bcftools call -m -v -o AM_000035F_sorted_q20_variants.bcf AM_000035F_sorted_q20.bcf
bcftools call -m -v -o AM_000058F_sorted_q20_variants.bcf AM_000058F_sorted_q20.bcf
bcftools call -m -v -o  KG_000058F_sorted_q20_variants.bcf KG_000058F_sorted_q20.bcf
bcftools call -m -v -o  KG_000035F_sorted_q20_variants.bcf KG_000035F_sorted_q20.bcf

#Filter and report the SNP variants in variant calling format (VCF)
vcfutils.pl varFilter  AM_000035F_sorted_q20_variants.bcf >  AM_000035F_sorted_q20_variants_filtered.vcf
vcfutils.pl varFilter AM_000058F_sorted_q20_variants.bcf > AM_000058F_sorted_q20_variants_filtered.vcf
vcfutils.pl varFilter KG_000058F_sorted_q20_variants.bcf  > KG_000058F_sorted_q20_variants_filtered.vcf
vcfutils.pl varFilter KG_000035F_sorted_q20_variants.bcf > KG_000035F_sorted_q20_variants_filtered.vcf

###filter DP read depth 
bcftools view -i 'INFO/DP>3 & INFO/DP<80' AM_000035F_sorted_q20_variants_filtered.vcf  >AM_000035F_sorted_q20_variants_filtered_DP.vcf
bcftools view -i 'INFO/DP>3 & INFO/DP<80' AM_000058F_sorted_q20_variants_filtered.vcf > AM_000058F_sorted_q20_variants_filtered_DP.vcf
bcftools view -i 'INFO/DP>3 & INFO/DP<80' KG_000058F_sorted_q20_variants_filtered.vcf >KG_000058F_sorted_q20_variants_filtered_DP.vcf
bcftools view -i 'INFO/DP>3 & INFO/DP<80' KG_000035F_sorted_q20_variants_filtered.vcf > KG_000035F_sorted_q20_variants_filtered_DP.vcf
