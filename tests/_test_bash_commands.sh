./build/adfinder \
	--vcf tests/data/ancestry_blocks.vcf \
	--out test_excl_intra \
	--block-size 1024 \
	--intra \
	--min-abs-r 0.5 \
	--chr scaffold-mi7 \
	--chr scaffold-mi8 \
	--chr scaffold-mi9 \
	--chr scaffold-mi10 \
	--hi-mode excl-focus

./build/adfinder \
	--vcf tests/data/ancestry_blocks.vcf \
	--out test_wo_excl_intra \
	--block-size 1024 \
	--intra \
	--min-abs-r 0.5 \
	--chr scaffold-mi7 \
	--chr scaffold-mi8 \
	--chr scaffold-mi9 \
	--chr scaffold-mi10

./build/adfinder \
	--vcf tests/data/ancestry_blocks.vcf \
	--out test_excl_inter \
	--block-size 1024 \
	--min-abs-r 0.5 \
	--chr scaffold-mi7 \
	--chr scaffold-mi8 \
	--chr scaffold-mi9 \
	--chr scaffold-mi10 \
	--hi-mode excl-focus

./build/adfinder \
	--vcf tests/data/ancestry_blocks.vcf \
	--out test_excl_target \
	--block-size 1024 \
	--min-abs-r 0.5 \
	--target-chr scaffold-mi7 \
	--target-pos 12421 \
	--chr scaffold-mi7 \
	--chr scaffold-mi8 \
	--chr scaffold-mi9 \
	--chr scaffold-mi10 \
	--hi-mode excl-focus

./build/adfinder \
	--vcf tests/data/ancestry_blocks.vcf \
	--out test_excl_sample_geno \
	--block-size 1024 \
	--min-abs-r 0.3 \
	--sample-geno tests/data/concvir_mito.tsv \
	--chr scaffold-mi7 \
	--chr scaffold-mi8 \
	--chr scaffold-mi9 \
	--chr scaffold-mi10 \
	--hi-mode excl-focus

./build/adfinder \
	--vcf tests/data/ancestry_blocks.vcf \
	--out test_excl_perm \
	--block-size 1024 \
	--chr scaffold-mi7 \
	--chr scaffold-mi8 \
	--chr scaffold-mi9 \
	--chr scaffold-mi10 \
	--permute 4 \
	--threads 4 \
	--hi-mode excl-focus

./build/adfinder \
	--vcf tests/data/ancestry_blocks.vcf \
	--out test_excl_perm_sample_geno \
	--block-size 1024 \
	--sample-geno tests/data/concvir_mito.tsv \
	--chr scaffold-mi7 \
	--chr scaffold-mi8 \
	--chr scaffold-mi9 \
	--chr scaffold-mi10 \
	--permute 4 \
	--threads 4 \
	--hi-mode excl-focus

./build/adfinder \
	--vcf tests/data/ancestry_blocks.vcf \
	--out test_compute_hi2 \
	--block-size 1024 \
	--no-chr scaffold-ma4 \
	--no-chr scaffold-ma3 \
	--compute-hi