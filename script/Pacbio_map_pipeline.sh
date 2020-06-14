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

###running the code AM reads
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
ls 
###filtering the low quality of mapping reads
/usr/local/bio/bwa.kit/samtools view -q 20 -b KG_TCAGTT_paired_sorted.bam > KG_TCAGTT_paired_sorted_q20.bam
/usr/local/bio/bwa.kit/samtools view -q 20 -b  AM_AGAGAT_paired_sorted.bam > AM_AGAGAT_paired_sorted_q20.bam

#############################################################Step 4 i######################################
###########################################################################################################
###########################################################################################################
####calling SNPs from BAM file
####Index the genome assembly (again!)
/usr/local/bio/bwa.kit/samtools faidx Ajap_PacBio.fa
/usr/local/bio/bwa.kit/samtools mpileup -g -f MvSl-1064-A2-R4.fa LAT3T_mappingto_MldSil1_1604a2_sorted_q20.bam > LAT3T_mappingto_MldSil1_1604a2_sorted_q20.bcf
/usr/local/bio/bwa.kit/samtools mpileup -g -f MvSl-1064-A2-R4.fa LAT3D_mappingto_MldSil1_1064a2_sorted_q20.bam > LAT3D_mappingto_MldSil1_1064a2_sorted_q20.bcf


###using bcftools
##(optional)first, calculate the read coverage of positions in the genome
bcftools mpileup -O b -o raw.bcf -f ref.fasta aligned.sorted.bam (or this step can be done using samtools as earlier step indicated to generate the raw .bcf file)

##second, Detect the single nucleotide polymorphisms (SNPs)
bcftools call --ploidy 1 -m -v -o LAT3T_mappingto_MldSil1_1604a2_sorted_q20_variants.bcf LAT3T_mappingto_MldSil1_1604a2_sorted_q20.bcf
bcftools call --ploidy 1 -m -v -o LAT3D_mappingto_MldSil1_1064a2_sorted_q20_variants.bcf LAT3D_mappingto_MldSil1_1064a2_sorted_q20.bcf
bcftools call --ploidy 1 -m -v -o LAT3T_mappingto_MsdSdi_1303D_sorted_q20_variants.bcf LAT3T_mappingto_MsdSdi_1303D_sorted_q20.bcf
bcftools call --ploidy 1 -m -v -o LAT3D_mappingto_MsdSdi_1303D_sorted_q20_variants.bcf LAT3D_mappingto_MsdSdi_1303D_sorted_q20.bcf

#Filter and report the SNP variants in variant calling format (VCF)
vcfutils.pl varFilter LAT3T_mappingto_MldSil1_1604a2_sorted_q20_variants.bcf > LAT3T_mappingto_MldSil1_1604a2_sorted_q20_variants_filtered.vcf
vcfutils.pl varFilter LAT3D_mappingto_MldSil1_1064a2_sorted_q20_variants.bcf > LAT3D_mappingto_MldSil1_1064a2_sorted_q20_variants_filtered.vcf
vcfutils.pl varFilter LAT3T_mappingto_MsdSdi_1303D_sorted_q20_variants.bcf > LAT3T_mappingto_MsdSdi_1303D_sorted_q20_variants_filtered.vcf
vcfutils.pl varFilter LAT3D_mappingto_MsdSdi_1303D_sorted_q20_variants.bcf > LAT3D_mappingto_MsdSdi_1303D_sorted_q20_variants_filtered.vcf

###filter DP read depth 
bcftools view -i 'INFO/DP>3 & INFO/DP<80' LAT3T_mappingto_MldSil1_1604a2_sorted_q20_variants_filtered.bcf > LAT3T_mappingto_MldSil1_1604a2_sorted_q20_variants_filtered_DP.bcf
bcftools view -i 'INFO/DP>3 & INFO/DP<80' LAT3D_mappingto_MldSil1_1064a2_sorted_q20_variants_filtered.bcf > LAT3D_mappingto_MldSil1_1064a2_sorted_q20_variants_filtered_DP.Vcf
bcftools view -i 'INFO/DP>3 & INFO/DP<80' LAT3T_mappingto_MsdSdi_1303D_sorted_q20_variants_filtered.vcf > LAT3T_mappingto_MsdSdi_1303D_sorted_q20_variants_filtered_DP.vcf
bcftools view -i 'INFO/DP>3 & INFO/DP<80' LAT3D_mappingto_MsdSdi_1303D_sorted_q20_variants_filtered.vcf > LAT3D_mappingto_MsdSdi_1303D_sorted_q20_variants_filtered_DP.vcf

bcftools view LAT3D_S12_L001_R1R2_mappingto_MldSil1_1064_a1_variantall_filterDP.bcf > LAT3D_S12_L001_R1R2_mappingto_MldSil1_1064_a1_variantall_filterDP.vcf

####plot R figures. 
#########################

#### extract individual sample file using bcftools
##bcftools view -s samplename,samplename file.vcf > samples_wanted.vcf

#### play with VCFtools for various parameters.
/usr/local/bio/vcftools/vcftools --SNPdensity 1000 --vcf LAT3D_S12_L001_R1R2_mappingto_MldSil1_1064_a1_final.vcf --out LAT3D_S12_L001_R1R2_mappingto_MldSil1_1064a1_snp.txt
/usr/local/bio/vcftools/vcftools --het --vcf LAT3D_S12_L001_R1R2_mappingto_MldSil1_1064_a1_final.vcf --out LAT3D_S12_L001_R1R2_mappingto_MldSil1_1064a1_het.txt ###this part #did not work, think that is due to heterozygosity only considered to occur in diploid, not haploid, individuals. 

###could use R directly plot SNP density, it works, check out the R code in /script folder.

###The SNP density of MAT looks very odd.
### selecting reads which were only properly paired /usr/local/bio/bwa.kit/samtools view -b -f 2 -F 524 LAT3T_S11_L001_R1R2_mappingto_MldSil1_1604_a1_sorted_q20.bam > LAT3T_S11_L001_R1R2_mappingto_MldSil1_1604_a1_sorted_q20_pairedonly.bam

###Results still does not look very promising..

### extract SNPs from only single copy gene between 1064a1 and a2 using vcftools.
vcftools --vcf LAT3D_S12_L001_R1R2_mappingto_MldSil1_1064_a1_final.vcf --bed 1064a1_singlecopy_cds_locations.bed --out 1064a1_singlecopy_cds_locations_LAT3D_snp.vcf --recode --keep-INFO-all

vcftools --vcf LAT3T_S11_L001_R1R2_mappingto_MsdSdi_1303T_variants_filtered.vcf --bed MsdSdi_1303T_singlecopy_coding_genomelocation_sorted1.bed --out 1303T_singlecopy_cds_locations_LAT3T.vcf --recode --keep-INFO-all

vcftools --vcf LAT3D_S12_L001_R1R2_mappingto_MsdsDI_1303T_variants_filtered.vcf --bed MsdSdi_1303T_singlecopy_coding_genomelocation_sorted1.bed --out 1303T_singlecopy_cds_locations_LAT3D.vcf --recode --keep-INFO-all

vcftools --vcf LAT3D_S12_L001_R1R2_mappingto_MldSil1_1064_a1_final.vcf --bed  MldSil1_hemi_80perc_genomeloc1.bed --out LAT3D_map1064a1_1064a1_hemilocation_snp.vcf --recode --keep-INFO-all

vcftools --vcf LAT3D_S12_L001_R1R2_mappingto_MsdsDI_1303T_variants_filtered.vcf --bed  MsdSdi2_hemi_80perc_genomeloc1.bed --out LAT3D_map1303T_MsdSdi2_hemi_80perc_snp.vcf --recode --keep-INFO-all

vcftools --vcf LAT3T_S11_L001_R1R2_mappingto_MldSil1_1604a1_variants_filtered.vcf --bed  MldSil1_hemi_80perc_genomeloc1.bed --out LAT3T_map1064a1_1064a1_hemi_80perc_snp.vcf --recode --keep-INFO-all

vcftools --vcf LAT3T_S11_L001_R1R2_mappingto_MsdSdi_1303T_variants_filtered.vcf --bed  MsdSdi2_hemi_80perc_genomeloc1.bed --out LAT3T_map1303T_1303T_hemi_80perc_snp.vcf --recode --keep-INFO-all

###convert bcf to vcf file
bcftools view my.bcf > my.vcf

###filter DP read depth 
bcftools view -i 'INFO/DP>3 & INFO/DP<80' LAT3D_S12_L001_R1R2_mappingto_MldSil1_1064_a1_variantall.bcf > LAT3D_S12_L001_R1R2_mappingto_MldSil1_1064_a1_variantall_filterDP.bcf
bcftools view LAT3D_S12_L001_R1R2_mappingto_MldSil1_1064_a1_variantall_filterDP.bcf > LAT3D_S12_L001_R1R2_mappingto_MldSil1_1064_a1_variantall_filterDP.vcf
