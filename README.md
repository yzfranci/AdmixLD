# AdmixLD

A command-line tool for detecting **ancestry disequilibrium (AD)** in hybrid zones and admixed populations.

---

## Overview

AdmixLD scans for pairwise correlations between residualized ancestry markers — correlations that remain after conditioning on each sample's genome-wide hybrid index (ancestry proportion). Strong signals can indicate selection for conspecific loci (positive r values) or heterospecific loci (negative r values).

Each marker's dosage is residualized on the genome-wide hybrid index (HI), controlling for linkage disequilibrium due to admixture alone. Pearson correlations between residualized markers are then computed across the genome (interchromosomal by default, intrachromosomal with `--intra`).

---

## Requirements & Installation

### Linux (Ubuntu / Debian)

```bash
# HTSlib (VCF/BCF I/O)
sudo apt-get install -y libhts-dev pkg-config

# Eigen3 (linear algebra)
sudo apt-get install -y libeigen3-dev

# OpenMP — bundled with GCC; install GCC if not already present
sudo apt-get install -y gcc g++
```

#### Build

```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

#### Add to PATH

```bash
echo "export PATH=\"$(pwd)/build:\$PATH\"" >> ~/.bashrc
source ~/.bashrc
```

---

### macOS (Homebrew)

```bash
# HTSlib, pkg-config, Eigen3
brew install htslib pkg-config eigen

# GCC with OpenMP (replaces Apple Clang for this build)
brew install gcc
```

#### Build

macOS ships `gcc` as a symlink to Apple Clang, so CMake must be pointed at the Homebrew GCC.

```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=gcc-<GCC_VERSION> \
  -DCMAKE_CXX_COMPILER=g++-<GCC_VERSION>
cmake --build build -j
```

#### Add to PATH

```bash
echo "export PATH=\"$(pwd)/build:\$PATH\"" >> ~/.zshrc
source ~/.zshrc
```

---

## Input Files

AdmixLD accepts two input formats:

- **VCF/BCF** (`--vcf`): One record per ancestry marker. Each sample must have a `DS` (dosage, in [0, 2]) or `GT` (diploid genotype) FORMAT field. Compressed files (`.vcf.gz`, `.bcf`) are supported via HTSlib.
- **MSP** (`--msp`): `.msp.tsv` format produced by **RFMix2** or **Gnomix**. Each row is a local ancestry tract; per-sample dosage is the sum of the two haplotype ancestry calls (0, 1, or 2).

Optional inputs:

- **`--cov FILE`**: Pre-computed covariate matrix (TSV) — one or more columns (e.g., hybrid index or PCA scores) — used instead of the internally computed HI for residualization.
- **`--bed FILE`**: BED file of genomic intervals; only markers within these regions are retained.
- **`--sample-geno FILE`**: Per-sample numeric trait TSV; correlates the trait against all ancestry markers instead of running a pairwise scan.

---

## Usage examples

```bash
# Genome-wide interchromosomal scan writing hits with a minimum aboslute value of r
# Also outputs the genome-wide distribution of r values
./build/admixld --vcf example_data_files/data.vcf.gz --out results \
	--min-abs-r 0.5 \
	--distrib \
	--chr scaffold-ma1 \
	--chr scaffold-ma2

# Interchromosomal scan writing hits with a minimum aboslute value of r
# Filtered for chromosome 1 and chromosomme 2
./build/admixld --vcf example_data_files/data.vcf.gz --out results \
	--min-abs-r 0.5 \
	--chr scaffold-ma1 \
	--chr scaffold-ma2

# Intrachromosomal scan with a distance limit
# Filtered for chromosome 1
./build/admixld --vcf example_data_files/data.vcf.gz --out results \
	--intra \
	--max-dist 5000000 \
	--chr scaffold-ma1

# Genome-wide intterchromosomal scan with LOCO (leave-one-chromosome-out) residualization
./build/admixld --vcf example_data_files/data.vcf.gz --out results \
	--hi-mode excl-focus \
	--min-abs-r 0.5

# Permutation test
./build/admixld --vcf data.vcf.gz --out permutation \
	--permute 10 \
	--threads 8 --seed 123

# Compute hybrid index only and exit
./build/admixld --vcf data.vcf.gz --out results --compute-hi
```
---

## Outputs

| File | Written when | Description |
|------|-------------|-------------|
| `<prefix>.hits.tsv` | Always (unless `--compute-hi`) | Marker pairs exceeding the correlation threshold |
| `<prefix>.hi.tsv` | When HI is computed | Per-sample hybrid index |
| `<prefix>.perm.summary.tsv` | `--permute N` | Per-replicate permutation statistics |
| `<prefix>.scan.summary.tsv` | `--distrib` | Empirical scan distribution summary |

The hits file columns are: `wA`, `chrA`, `posA`, `wB`, `chrB`, `posB`, `r`, `n`.

See [documentation.md](documentation.md) for the full flag reference, file format specifications, and method details.
