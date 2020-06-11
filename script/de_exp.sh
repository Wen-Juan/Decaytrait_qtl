/usr/local/bio/kallisto/kallisto
/usr/local/bio/trinityrnaseq/util/R/

/usr/local/bio/trinityrnaseq/Analysis/DifferentialExpression/run_DE_analysis.pl \
--matrix Atab_.isoform.counts.matrix \
--method edgeR \
--min_reps_min_cpm 3 \
--output DE_exp \
--dispersion 0.1
