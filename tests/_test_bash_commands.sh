./build/admixld \
	--vcf tests/data/ancestry_blocks.vcf \
	--out test_excl_intra \
	--tile-size 1024 \
	--intra \
	--min-abs-r 0.5 \
	--chr scaffold-mi7 \
	--chr scaffold-mi8 \
	--chr scaffold-mi9 \
	--chr scaffold-mi10 \
	--hi-mode excl-focus

./build/admixld \
	--vcf tests/data/ancestry_blocks.vcf \
	--out test_wo_excl_intra \
	--tile-size 1024 \
	--intra \
	--min-abs-r 0.5 \
	--chr scaffold-mi7 \
	--chr scaffold-mi8 \
	--chr scaffold-mi9 \
	--chr scaffold-mi10

./build/admixld \
	--vcf tests/data/ancestry_blocks.vcf \
	--out test_excl_inter \
	--tile-size 1024 \
	--min-abs-r 0.5 \
	--chr scaffold-mi7 \
	--chr scaffold-mi8 \
	--chr scaffold-mi9 \
	--chr scaffold-mi10 \
	--hi-mode excl-focus

./build/admixld \
	--vcf tests/data/ancestry_blocks.vcf \
	--out test_excl_inter_distrib \
	--tile-size 1024 \
	--min-abs-r 0.5 \
	--chr scaffold-mi7 \
	--chr scaffold-mi8 \
	--chr scaffold-mi9 \
	--chr scaffold-mi10 \
	--distrib \
	--hi-mode excl-focus

./build/admixld \
	--vcf tests/data/ancestry_blocks.vcf \
	--out test_excl_target \
	--tile-size 1024 \
	--min-abs-r 0.5 \
	--target-chr scaffold-mi7 \
	--target-pos 12421 \
	--chr scaffold-mi7 \
	--chr scaffold-mi8 \
	--chr scaffold-mi9 \
	--chr scaffold-mi10 \
	--hi-mode excl-focus

./build/admixld \
	--vcf tests/data/ancestry_blocks.vcf \
	--out test_excl_sample_geno \
	--tile-size 1024 \
	--min-abs-r 0.3 \
	--sample-geno tests/data/concvir_mito.tsv \
	--chr scaffold-mi7 \
	--chr scaffold-mi8 \
	--chr scaffold-mi9 \
	--chr scaffold-mi10 \
	--hi-mode excl-focus

./build/admixld \
	--vcf tests/data/ancestry_blocks.vcf \
	--out test_excl_perm \
	--tile-size 1024 \
	--chr scaffold-mi7 \
	--chr scaffold-mi8 \
	--chr scaffold-mi9 \
	--chr scaffold-mi10 \
	--permute 4 \
	--threads 4 \
	--hi-mode excl-focus

./build/admixld \
	--vcf tests/data/ancestry_blocks.vcf \
	--out test_excl_perm_sample_geno \
	--tile-size 1024 \
	--sample-geno tests/data/concvir_mito.tsv \
	--chr scaffold-mi7 \
	--chr scaffold-mi8 \
	--chr scaffold-mi9 \
	--chr scaffold-mi10 \
	--permute 4 \
	--threads 4 \
	--hi-mode excl-focus

./build/admixld \
	--vcf tests/data/ancestry_blocks.vcf \
	--out test_compute_hi2 \
	--tile-size 1024 \
	--no-chr scaffold-ma4 \
	--no-chr scaffold-ma3 \
	--compute-hi

./build/admixld \
	--vcf tests/data/Orioles_filtered_final_093020_55individuals.recode.vcf.gz \
	--out test_missing \
	--max-windows 5000 \
	--intra \
	--min-callrate 0.9 \
	--min-abs-r 0.4 \
	--tile-size 1024