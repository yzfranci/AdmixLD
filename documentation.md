# AdmixLD — Full Documentation

## Overview

AdmixLD is a command-line tool for detecting **ancestry disequilibrium (AD)** in hybrid zones and admixed populations. It scans for pairwise correlations between residualized ancestry markers — correlations that remain after controlling for the genome-wide hybrid index (ancestry proportion). Strong signals can indicate selection for conspecific loci (positive r values) or heterospecific loci (negative r values).

**Core idea**: Each marker's dosage is residualized on the genome-wide hybrid index (HI), controlling for linkage disequilibrium due to admixture alone. Pearson correlations between residualized markers are then computed across the genome (interchromosomal by default, intrachromosomal with `--intra`).

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

OpenMP support is detected automatically by CMake. The build succeeds without it, but `--threads` will have no effect.

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
# Homebrew (if not already installed)
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# HTSlib, pkg-config, Eigen3
brew install htslib pkg-config eigen

# GCC with OpenMP support
brew install gcc
```

#### Build

macOS ships `gcc` as a symlink to Apple Clang, so CMake must be pointed at the real Homebrew GCC.

```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=gcc-<GCC_VER> \
  -DCMAKE_CXX_COMPILER=g++-<GCC_VER>
cmake --build build -j
```

#### Add to PATH

```bash
echo "export PATH=\"$(pwd)/build:\$PATH\"" >> ~/.zshrc
source ~/.zshrc
```

The binary is placed at `build/admixld`.

---

## Input File Formats

### VCF/BCF — Ancestry Markers (`--vcf`)

The primary input is a VCF or BCF file where each record is one **ancestry marker**. Each sample must have either:

- `DS` (dosage): a continuous value in [0, 2] representing ancestry dosage
- `GT` (genotype): a standard diploid genotype (e.g., `0|1`, `1/1`)

If both are present, `DS` takes precedence. Compressed input (`.vcf.gz`, `.bcf`) is supported via HTSlib.

**Position interpretation**: By default, the VCF `POS` field is interpreted as the **end** of the marker's inferred segment. Use `--pos-is-start` if your VCF records marker start positions instead. This affects inferred segment lengths used for weighted HI computation.

If the ancestry markers represent singular SNPs instead of ancestry tracts, make sure to use --unweighted-hi for hybird index computation.

**Missing data**: Missing genotypes are stored internally as NaN. By default, any marker with missing data in any sample is excluded. See [Missing Data](#missing-data) for relaxed handling.

### MSP File (`--msp`)

Alternative to `--vcf`. Accepts the `.msp.tsv` format produced by **RFMix2** and **Gnomix**. Each data row represents one local ancestry tract; the dosage for each sample is the sum of the two haplotype ancestry calls (0, 1, or 2).

The file must begin with exactly two comment header lines:

```
#Subpopulation order/codes: Pop0=0    Pop1=1
#chm    spos    epos    sgpos    egpos    n snps    sample1.0    sample1.1    sample2.0    sample2.1    ...
```

- **`chm`** — chromosome name
- **`spos` / `epos`** — integer start and end positions of the tract (both used directly; `--pos-is-start` has no effect for MSP input)
- **`sgpos` / `egpos` / `n snps`** — genetic positions and SNP count (read but not used)
- **`sample.0` / `sample.1`** — haplotype ancestry calls, must be `0` or `1`; any other value is treated as missing

### Covariate TSV (`--cov`)

Optional user-supplied covariate matrix used instead of the computed hybrid index for residualization. Tab-separated with an optional header row. Must contain all samples present in the input file.

Supports one or more covariate columns (e.g., a single hybrid index, or PC1–PC3 from a PCA, etc...). Cannot handle missing values. 

```
sample    cov1      cov2      cov3
CV0340    0.373028  0.021     -0.014
CV0359    0.213722  -0.003    0.008
CV0362    0.891145  0.041     0.002
```
### Sample Trait TSV (`--sample-geno`)

Per-sample numeric haplotype/trait for correlating against all ancestry markers. Tab-separated. Missing values can be encoded as `.`, `NA`, `NaN`, `nan`, or `NAN`.

This can be used to scan for correlation with the mitochondrial haplotype.

```
sample    trait
CV0340    1
CV0359    .
CV0362    0
```

### BED Region Filter (`--bed`)

A tab-separated BED file defining genomic intervals. Only markers whose position falls within a BED interval are retained.

```
scaffold-mi7    100000    200000
scaffold-mi8    50000     150000
```

### Sample Keep List (`--keep-indv`)

A plain-text file listing one sample ID per line. Only these samples are retained after loading. All other samples are dropped before any filtering or computation.

---

## Output File Formats

### Hits File — `<prefix>.hits.tsv`

Main output. Tab-separated. Each row is a marker pair exceeding the correlation threshold.

| Column | Description |
|--------|-------------|
| `wA`   | Index of marker A |
| `chrA` | Chromosome of marker A |
| `posA` | Position of marker A |
| `wB`   | Index of marker B |
| `chrB` | Chromosome of marker B |
| `posB` | Position of marker B |
| `r`    | Pearson correlation coefficient |
| `n`    | Number of samples used |

```
wA    chrA          posA      wB    chrB           posB      r           n
0     scaffold-mi9  421750    617   scaffold-mi10  103231    0.110808    104
0     scaffold-mi9  421750    636   scaffold-mi10  819774    0.125357    104
```

### Hybrid Index File — `<prefix>.hi.tsv`

Written when HI is computed (always if not pre-supplied). Tab-separated.

| Column   | Description |
|----------|-------------|
| `sample` | Sample name |
| `hi`     | Hybrid index (ancestry proportion, in [0, 1]) |

### Permutation Summary — `<prefix>.perm.summary.tsv`

Written when `--permute N` is used. One row per permutation replicate.

| Column      | Description |
|-------------|-------------|
| `rep`       | Replicate index |
| `max_r`     | Maximum correlation in replicate |
| `min_r`     | Minimum correlation in replicate |
| `p99`..`p01`| Percentiles of correlation distribution |
| `mean`      | Mean correlation |
| `sd`        | Standard deviation of correlations |
| `mean_r2`   | Mean of r² values |
| `sd_r2`     | SD of r² values |

### Scan Distribution Summary — `<prefix>.scan.summary.tsv`

Written when `--distrib` is used. Single row summarising the empirical distribution of all tested pairs.

| Column         | Description |
|----------------|-------------|
| `tested_pairs` | Total number of pairs tested |
| `max_r`, `min_r` | Global extremes |
| Percentiles    | `p99`, `p95`, `p75`, `median`, `p25`, `p05`, `p01` |
| `mean`, `sd`   | Mean and SD |
| `mean_r2`, `sd_r2` | R² statistics |

---

## Command-Line Reference

### Required

| Flag | Description |
|------|-------------|
| `--vcf FILE` | VCF/BCF input file (mutually exclusive with `--msp`) |
| `--msp FILE` | RFMix2/Gnomix `.msp.tsv` file (mutually exclusive with `--vcf`) |
| `--out PREFIX` | Output file prefix |

### Scanning Behaviour

| Flag | Default | Description |
|------|---------|-------------|
| `--min-abs-r FLOAT` | 0.0 | Minimum absolute Pearson r to report |
| `--min-neg-r FLOAT` | — | Keep only pairs with r ≤ −value (negative correlations) |
| `--min-pos-r FLOAT` | — | Keep only pairs with r ≥ value (positive correlations) |
| `--intra` | off | Scan intrachromosomal pairs (default: interchromosomal only) |
| `--phased` | off | Use phased haplotypes (2n rows); requires `--intra`; VCF GT or MSP only |
| `--max-dist INT` | unlimited | Maximum basepair distance for intrachromosomal pairs |
| `--max-windows N` | all | Load at most N markers (useful for quick testing) |

When either asymmetric filter (`--min-neg-r` or `--min-pos-r`) is set, `--min-abs-r` is ignored and a pair is kept if it satisfies either directional threshold.

### Hybrid Index / Covariates

| Flag | Default | Description |
|------|---------|-------------|
| `--cov FILE` | computed | Load covariate matrix from TSV (one or more columns); replaces computed HI |
| `--hi-mode STR` | `global` | HI correction mode: `global` or `excl-focus` (LOCO; incompatible with `--cov`) |
| `--compute-hi` | off | Compute and save HI, then exit without scanning (ignored with `--cov`) |
| `--unweighted-hi` | off | Use unweighted (simple mean) HI instead of length-weighted (ignored with `--cov`) |
| `--pos-is-start` | off | Treat VCF POS as marker start rather than marker end; no effect for `--msp` input |

LOCO (leave-one-chromosome-out) can be used to compute HI without the focus chromosome(s), to avoid confounding the focal chromosome's signal with its own contribution to the global HI estimate (Schumer and Brandvain, 2016).

### Filtering

| Flag | Default | Description |
|------|---------|-------------|
| `--chr STR` | all | Keep only this chromosome (repeatable) |
| `--no-chr STR` | none | Exclude this chromosome (repeatable) |
| `--bed FILE` | none | Keep only markers within BED intervals |
| `--min-callrate FLOAT` | 1.0 | Minimum call-rate to retain a marker; enables mean imputation if < 1.0 |
| `--keep-indv FILE` | all | Keep only samples listed in FILE (one ID per line) |

Be sure to use --min-callrate for loci that still have a decent amount of data (<10% missing), as it can introduce bias due to mean imputation.

### Target Marker

Used to scan correlations between a marker of interest against the rest of the markers.

| Flag | Description |
|------|-------------|
| `--target-chr STR` | Chromosome of single focus marker |
| `--target-pos INT` | Position of single focus marker (exact match); reports that marker vs all others |

### Sample Trait

Typically used to provide a mitochondrial haplotype that would not be contained in the VCF/MSP file.

| Flag | Description |
|------|-------------|
| `--sample-geno FILE` | Per-sample numeric trait TSV; correlates haplotype/trait against all markers |

### Permutation Testing

Estimates the genome-wide null distribution of r by repeatedly shuffling sample labels and re-scanning. Can be used to choose a r value threshold for outliers.

| Flag | Default | Description |
|------|---------|-------------|
| `--permute N` | 0 | Run N full-shuffle permutation replicates |
| `--permute-sample INT` | 200000 | Reservoir sample size for percentile estimation |
| `--seed INT` | 1 | RNG seed |

### Distribution Summary

| Flag | Default | Description |
|------|---------|-------------|
| `--distrib` | off | Write empirical scan distribution summary |
| `--distrib-sample INT` | 200000 | Reservoir sample size |

### Performance

| Flag | Default | Description |
|------|---------|-------------|
| `--threads INT` | 1 | OpenMP threads for parallel scanning |
| `--tile-size INT` | 1024 | Matrix tile size for cache efficiency |

---

## Methods

### Hybrid Index Computation

The hybrid index (HI) estimates the genome-wide ancestry proportion of each sample.

**Weighted HI** (default):

$$\text{HI}_i = \frac{\sum_j d_{ij} \cdot L_j}{2 \sum_j L_j}$$

where $d_{ij}$ is the dosage of sample $i$ at marker $j$, and $L_j$ is the inferred segment length in base pairs. For MSP input, $L_j$ is the exact tract length ($\text{epos}_j - \text{spos}_j$). For VCF input, lengths are inferred from consecutive marker positions on each chromosome.

**Unweighted HI** (`--unweighted-hi`):

$$\text{HI}_i = \frac{1}{2n} \sum_j d_{ij}$$

where $n$ is the number of markers. Equivalent to treating all markers equally regardless of physical length.

### Residualization

Each marker's dosage vector is residualized on the covariate matrix $\mathbf{H}$ (samples × $k$) via ordinary least squares. By default $\mathbf{H}$ is the single-column genome-wide hybrid index; with `--cov` it can be any number of columns (e.g., PC1–PC3 from a PCA).

1. Centre each covariate column: $H^c_{il} = H_{il} - \bar{H}_l$
2. Regress dosage of marker $j$ on covariates: $\hat{\boldsymbol{\beta}}_j = (\mathbf{H}^{cT} \mathbf{H}^c)^{-1} \mathbf{H}^{cT} \mathbf{d}_j$
3. Residuals: $\boldsymbol{\epsilon}_j = \mathbf{d}_j - \bar{d}_j \mathbf{1} - \mathbf{H}^c \hat{\boldsymbol{\beta}}_j$
4. Standardise: $Z_{ij} = (\epsilon_{ij} - \bar{\epsilon}_j) / \text{SD}(\boldsymbol{\epsilon}_j)$

The $k \times k$ system $(\mathbf{H}^{cT} \mathbf{H}^c)$ is factored once (LDLT) and reused across all markers. Markers with fewer than $k+2$ samples, a singular covariate matrix, or zero residual SD are excluded.

### Correlation Scanning

Pearson correlations are computed between all pairs of residualized, standardised markers:

$$r_{AB} = \frac{1}{n-1} \sum_i Z_{Ai} \cdot Z_{Bi}$$

**Interchromosomal** (default): all pairs from different chromosomes.  
**Intrachromosomal** (`--intra`): all pairs from the same chromosome, optionally within `--max-dist` bp.

Matrix multiplication is tiled (default tile size 1024) for cache efficiency with large marker sets.

### Permutation Testing

Genome-wide significance thresholds are estimated by full-shuffle permutation:

1. Residualize all markers on HI to produce the z-scored residual matrix Z (computed once)
2. Draw a single random permutation of sample indices per replicate; for every interchromosomal pair (A, B), apply it to the markers of A while leaving B in natural sample order — the same permutation is reused across all pairs within a replicate
3. Scan all interchromosomal pairs on the permuted Z and record the distribution of r values
4. Repeat N times (`--permute N`)

A reservoir sampler with configurable capacity (`--permute-sample`) tracks the distribution without storing all correlations.

The distribution as well as the min and max r values can be used to choose a r value threshold for outliers.

### LOCO Mode (`--hi-mode excl-focus`)

In leave-one-chromosome-out (LOCO) mode, the HI used to residualize markers on chromosome $c$ excludes that chromosome's contribution:

$$\text{HI}^{-c}_i = \frac{\sum_{j \notin c} d_{ij} L_j}{2 \sum_{j \notin c} L_j}$$

This prevents the focal chromosome's own ancestry signal from inflating its residuals, which would suppress true signals on that chromosome while producing artifactual correlations between it and others. LOCO mode is recommended when the data contain strong chromosomal-scale ancestry tracts.

---

## Missing Data

By default (`--min-callrate 1.0`), AdmixLD excludes any marker that has missing genotype data in any sample. All residualization and correlation calculations are therefore performed on fully observed data.

When `--min-callrate <value>` is set to a value less than 1.0:

1. Markers with a call-rate (proportion of non-missing genotypes) ≥ the threshold are **retained**.
2. Remaining missing genotypes within retained markers are **mean-imputed** using the per-marker mean dosage.
3. Residualization and scanning proceed on the imputed matrix.

Mean imputation is applied only when explicitly requested. This option is intended primarily for reduced-representation datasets (e.g., RADseq) where missing data are common and a strict complete-data requirement would exclude most markers. Note that mean imputation can attenuate correlation estimates.

---

## Notes and Caveats

- **Multithreading**: Parallel scanning writes per-thread temporary files that are merged after the scan. Thread count has no effect on results, only runtime.
- **Weighted vs unweighted HI**: Weighted HI is preferred when markers vary substantially in inferred segment length (typical for local ancestry tract data obtained via Loter, RFMIX etc...). Unweighted HI may be appropriate for SNP arrays.
- **LOCO and permutations**: When using `--hi-mode excl-focus` with `--permute`, HI components are re-used across permutations for efficiency; only the mapping of HI to shuffled samples changes.
- **`--msp` vs `--vcf`**: MSP input provides explicit tract start/end positions, so segment lengths for weighted HI are exact. VCF input infers lengths from consecutive marker positions (affected by `--pos-is-start`).
- **`--cov` and LOCO**: `--cov` is incompatible with `--hi-mode excl-focus`. LOCO requires per-chromosome HI components derived from the marker data itself and cannot be computed for externally supplied covariates.
- **`--cov` missing data**: No missing values are permitted in the covariate file. All samples in the input must be present.
