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
./build/adfinder --vcf x.vcf.gz --out test --min-abs-r 0.5
```