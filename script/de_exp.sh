perl run_DE_analysis.pl \
--matrix ../input/Atab_.isoform.counts.matrix \
--method edgeR \
--min_reps_min_cpm 3 \
--output ../output/DE_exp \
--dispersion 0.1
