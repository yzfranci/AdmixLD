# requirements
htslib
```bash
	sudo apt-get update
	sudo apt-get install -y libhts-dev pkg-config
```
eigen
```bash
	sudo apt-get install -y libeigen3-dev
```
# Build locally
```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

# Test
```bash
tests/test_smoke.sh
./build/adfinder --vcf tests/data/ancestry_blocks.vcf --out test --min-abs-r 0.5 --max-windows 200 --intra
```

# Run in parall
```bash
# interchrom permutation test
./build/adfinder --vcf tests/data/ancestry_blocks.vcf --out scan --max-windows 8000 --permute 10 --threads 4 --permute-sample 200000 --seed 123
```

# Other tests
```bash
# Empirical distribution test
./build/adfinder --vcf tests/data/ancestry_blocks.vcf --out distr --max-windows 8000 --distrib --seed 123 --min-abs-r 0.5
```


## As a reminder for the missing data handling
By default, adfinder requires complete genotype information across all samples within each ancestry block: blocks containing any missing genotypes are excluded from analysis, ensuring that all residualization and correlation calculations are performed on fully observed data. When the --min-callrate option is specified with a value less than 1.0, this strict requirement is relaxed. In this mode, ancestry blocks are first filtered to retain only those with a proportion of non-missing genotypes greater than or equal to the specified call-rate threshold. For the remaining blocks, any missing genotypes are then mean-imputed within each block prior to hybrid-index residualization. This imputation step is necessary to enable efficient matrix-based correlation scans while preserving the original hybrid-index correction framework. Mean imputation is applied only when --min-callrate < 1.0 is explicitly requested; the default behavior (--min-callrate 1.0) remains unchanged and excludes blocks with any missing data. This option is intended primarily for reduced-representation datasets (e.g., RADseq), where missing data are common and dense pairwise-complete correlation approaches would be computationally prohibitive.