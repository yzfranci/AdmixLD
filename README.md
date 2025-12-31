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
OMP_NUM_THREADS=4 ./build/adfinder --vcf tests/data/ancestry_blocks.vcf --out scan --max-windows 8000 --permute 20 --permute-sample 20000 --seed 123
```