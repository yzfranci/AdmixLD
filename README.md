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