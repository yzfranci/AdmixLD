# AdmixLD

A command-line tool for detecting **excess ancestry linkage disequilibrium** in hybrid zones and admixed populations.

---

## Overview

AdmixLD scans for pairwise correlations between residualized ancestry markers — correlations that remain after conditioning on each sample's genome-wide hybrid index (ancestry proportion). Strong signals can indicate selection for conspecific loci (positive r values) or heterospecific loci (negative r values).

Each marker's dosage is residualized on the genome-wide hybrid index (HI), controlling for linkage disequilibrium due to admixture alone. Pearson correlations between residualized markers are then computed across the genome (interchromosomal by default, intrachromosomal with `--intra`).

For more comprehensive documentation, see [documentation.md](documentation.md).

---

## Requirements & Installation

### Linux (Ubuntu / Debian)
In order to install on a Debian based system, you will need to install the following system level packages:

```bash
# HTSlib (VCF/BCF I/O)
sudo apt-get install -y libhts-dev pkg-config

# Eigen3 (linear algebra)
sudo apt-get install -y libeigen3-dev

# OpenMP — bundled with GCC; install GCC if not already present
sudo apt-get install -y gcc g++
```
#### Clone repository
After you have installed the necessary system level packages, clone the repository:

```bash
git clone https://github.com/yzfranci/AdmixLD.git
```

#### Build
Once you have cloned the repository, make sure to move into the new directory:
```bash
cd AdmixLD
```

Now that you are in that directory, run the following commands to build the program:
```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

#### Add to PATH
To make the executable callable from anywhere within the installed user, run the following to add the build directory to your PATH and restart your shell:
```bash
echo "export PATH=\"$(pwd)/build:\$PATH\"" >> ~/.bashrc
source ~/.bashrc
```

---

### macOS (Homebrew)
In order to install on a MacOS based system, you will need to install the following system level packages with brew. To install brew, follow the guide here: https://brew.sh/

```bash
# HTSlib, pkg-config, Eigen3
brew install htslib pkg-config eigen

# GCC with OpenMP (replaces Apple Clang for this build)
brew install gcc
```

#### Clone repository
After you have installed the necessary system level packages, clone the repository:

```bash
git clone https://github.com/yzfranci/AdmixLD.git
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
To make the executable callable from anywhere within the installed user, run the following to add the build directory to your PATH and restart your shell:

```bash
echo "export PATH=\"$(pwd)/build:\$PATH\"" >> ~/.zshrc
source ~/.zshrc
```

If you are using ``bash`` as your shell on MacOS, run this instead:
```bash
echo "export PATH=\"$(pwd)/build:\$PATH\"" >> ~/.bashrc
source ~/.bashrc
```

---

## Input Files

AdmixLD accepts two input formats:

- **VCF/BCF** (`--vcf`): One record per ancestry marker. Each sample must have a `DS` (dosage, in [0, 2]) or `GT` (diploid genotype) FORMAT field. Compressed files (`.vcf.gz`, `.bcf`) are supported via HTSlib.
- **MSP** (`--msp`): `.msp.tsv` format produced by **RFMix2** or **Gnomix**. Each row is a local ancestry tract; per-sample dosage is the sum of the two haplotype ancestry calls (0, 1, or 2).

Optional inputs:

- **`--cov FILE`**: External HI file (TSV: `sample<TAB>hi`, single column) — used instead of the internally computed HI for residualization.
- **`--ref-freq FILE`**: Parental allele frequency file (TSV: `chrom`, `pos`, `p1`, `p2`). When supplied with `--vcf`, enables frequency-based HI estimation and allele polarization (see [Frequency-Based HI](#frequency-based-hi-computation---ref-freq)).
- **`--bed FILE`**: BED file of genomic intervals; only markers within these regions are retained.
- **`--sample-haplo FILE`**: Per-sample mitochondrial haplotype TSV (values must be `0` or `1`); correlates the haplotype against all ancestry markers instead of running a pairwise scan.

---

## Usage examples

```bash
# Genome-wide interchromosomal scan writing hits with a minimum aboslute value of r
./build/admixld --vcf example_data_files/data.vcf.gz --out results \
	--min-abs-r 0.5 \
	--chr scaffold-ma1 \
	--chr scaffold-ma2

# Genome-wide interchromosomal scan writing hits with a minimum aboslute value of r
# Also outputs the genome-wide distribution of r values (with the sampled r values used to estimate the distribution)
./build/admixld --vcf example_data_files/data.vcf.gz --out results \
	--min-abs-r 0.5 \
	--distrib \
	--distrib-raw \
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

# Compute hybrid index only and exit
./build/admixld --vcf data.vcf.gz --out results --compute-hi
```
---

## Outputs

| File | Written when | Description |
|------|-------------|-------------|
| `<prefix>.hits.tsv` | Always (unless `--compute-hi`) | Marker pairs exceeding the correlation threshold |
| `<prefix>.hi.tsv` | When HI is computed | Per-sample hybrid index |
| `<prefix>.scan.summary.tsv` | `--distrib` | Empirical scan distribution summary (quantiles, mean, SD) |
| `<prefix>.scan.summary.reservoir.tsv` | `--distrib-raw` | Raw reservoir sample of r values (single column `r`) |

The hits file columns are: `wA`, `chrA`, `posA`, `wB`, `chrB`, `posB`, `r`, `n`.

See [documentation.md](documentation.md) for the full flag reference, file format specifications, and method details.

## Citation
TO BE COMPLETED