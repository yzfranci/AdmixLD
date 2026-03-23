# AdmixLD Documentation

## Overview

AdmixLD is a command-line tool for detecting **ancestry disequilibrium (AD)** in hybrid zones and admixed populations. It scans for statistically significant pairwise correlations between residualized ancestry blocks — correlations that remain after conditioning on global hybrid index (ancestry proportion). Such signals can indicate adaptive introgression, selection, or anomalous population structure at specific genomic loci.

**Core idea**: Each ancestry block's dosage is residualized on the genome-wide hybrid index (HI), removing background ancestry-linkage disequilibrium. Pearson correlations between residualized blocks across the genome then reveal statistically significant deviations.

---

## Requirements & Installation

### Dependencies

```bash
# HTSlib (VCF/BCF I/O)
sudo apt-get install -y libhts-dev pkg-config

# Eigen3 (linear algebra)
sudo apt-get install -y libeigen3-dev

# OpenMP (optional, for parallelization — typically included with GCC)
```

### Build

```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The binary is placed at `build/adfinder`.

---

## Quick Start

```bash
# Basic interchromosomal scan
./build/admixld --vcf data.vcf --out results --min-abs-r 0.5

# Run smoke tests
tests/test_smoke.sh
```

---

## Input File Formats

### VCF/BCF — Ancestry Blocks (`--vcf`)

The primary input is a VCF or BCF file where each record represents an **ancestry block** (e.g., from a local ancestry inference tool). Each sample must have either:

- `DS` (dosage): a continuous value in [0, 2] representing ancestry dosage
- `GT` (genotype): a standard diploid genotype (e.g., `0|1`, `1/1`)

If both are present, `DS` takes precedence.

**Position interpretation**: By default, the VCF `POS` field is interpreted as the **end** of the block. Use `--pos-is-start` if your VCF records block start positions instead. This affects inferred block lengths used for weighted HI computation.

**Missing data**: Missing genotypes are stored internally as NaN. By default, any block with missing data in any sample is excluded. See [Missing Data](#missing-data) for relaxed handling.

```
#CHROM          POS     ID  REF ALT ...  FORMAT      Sample1     Sample2
scaffold-mi9    421749  .   A   G   ...  GT:DS       0|0:0       1|0:0.5
scaffold-mi9    845000  .   A   G   ...  GT:DS       1|1:2       0|0:0
```

### Hybrid Index TSV (`--hi`)

Optional pre-computed hybrid index. Tab-separated with a header row. Must contain all samples present in the VCF.

```
sample    hi
CV0340    0.373028
CV0359    0.213722
CV0362    0.891145
```

### Sample Trait TSV (`--sample-geno`)

Per-sample numeric trait for correlating against all ancestry blocks. Tab-separated. Missing values can be encoded as `.`, `NA`, `NaN`, `nan`, or `NAN`. Requires at least 3 non-missing samples.

```
sample    trait
CV0340    1
CV0359    .
CV0362    0
```

### BED Region Filter (`--bed`)

A headerless, tab-separated BED file defining genomic intervals. Only blocks whose position falls within a BED interval are retained. Overlapping and adjacent intervals are merged automatically.

```
scaffold-mi7    100000    200000
scaffold-mi8    50000     150000
```

---

## Output File Formats

### Hits File — `<prefix>.hits.tsv`

Main output. Tab-separated. Each row is a block pair exceeding the correlation threshold.

| Column | Description |
|--------|-------------|
| `wA`   | Index of block A |
| `chrA` | Chromosome of block A |
| `posA` | Position of block A (1-based) |
| `wB`   | Index of block B |
| `chrB` | Chromosome of block B |
| `posB` | Position of block B (1-based) |
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
| `hi`     | Hybrid index (ancestry proportion, in [0,1]) |

### Permutation Summary — `<prefix>.perm.summary.tsv`

Written when `--permute N` is used. One row per permutation replicate. Contains:

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
| `--vcf FILE` | VCF/BCF input file |
| `--out PREFIX` | Output file prefix |

### Scanning Behaviour

| Flag | Default | Description |
|------|---------|-------------|
| `--min-abs-r FLOAT` | 0.0 | Minimum absolute Pearson r to report |
| `--intra` | off | Scan intrachromosomal pairs (default: interchromosomal only) |
| `--max-dist INT` | unlimited | Maximum basepair distance for intrachromosomal pairs |
| `--max-windows N` | all | Load at most N blocks (useful for testing) |

### Hybrid Index

| Flag | Default | Description |
|------|---------|-------------|
| `--hi FILE` | computed | Load pre-computed HI from TSV |
| `--hi-mode STR` | `global` | HI computation mode: `global` or `excl-focus` (LOCO) |
| `--compute-hi` | off | Compute and save HI, then exit without scanning |
| `--unweighted-hi` | off | Use unweighted (simple mean) HI instead of block-length weighted |
| `--pos-is-start` | off | Treat VCF POS as block start rather than block end |

### Filtering

| Flag | Default | Description |
|------|---------|-------------|
| `--chr STR` | all | Keep only this chromosome (repeatable) |
| `--no-chr STR` | none | Exclude this chromosome (repeatable) |
| `--bed FILE` | none | Keep only blocks within BED intervals |
| `--min-callrate FLOAT` | 1.0 | Minimum call-rate to retain a block; enables mean imputation if < 1.0 |

### Target Block

| Flag | Description |
|------|-------------|
| `--target-chr STR` | Chromosome of single focus block |
| `--target-pos INT` | Position of single focus block (exact match); reports that block vs all others |

### Sample Trait

| Flag | Description |
|------|-------------|
| `--sample-geno FILE` | Per-sample numeric trait TSV; correlates trait against all blocks |

### Permutation Testing

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

### Asymmetric Correlation Filters

| Flag | Description |
|------|-------------|
| `--min-neg-r FLOAT` | Keep only pairs with r ≤ −value (negative correlations) |
| `--min-pos-r FLOAT` | Keep only pairs with r ≥ value (positive correlations) |

When either asymmetric filter is set, `--min-abs-r` is ignored and a pair is kept if it satisfies either directional threshold.

---

## Methods

### Hybrid Index Computation

The hybrid index (HI) estimates the genome-wide ancestry proportion of each sample.

**Weighted HI** (default):

$$\text{HI}_i = \frac{\sum_j d_{ij} \cdot L_j}{2 \sum_j L_j}$$

where $d_{ij}$ is the dosage of sample $i$ at block $j$, and $L_j$ is the inferred block length in base pairs. Block lengths are inferred from consecutive VCF positions on each chromosome.

**Unweighted HI** (`--unweighted-hi`):

$$\text{HI}_i = \frac{1}{2n} \sum_j d_{ij}$$

where $n$ is the number of blocks. Equivalent to treating all blocks equally regardless of physical length.

### Residualization

Each block's dosage vector is residualized on HI via ordinary least squares:

1. Centre HI: $h_i^c = \text{HI}_i - \bar{\text{HI}}$
2. Regress: $\hat{\beta} = \text{cov}(\mathbf{d}_j, \mathbf{h}^c) / \text{var}(\mathbf{h}^c)$
3. Residuals: $\epsilon_{ij} = d_{ij} - \bar{d}_j - \hat{\beta} h_i^c$
4. Standardise: $Z_{ij} = (\epsilon_{ij} - \bar{\epsilon}_j) / \text{SD}(\boldsymbol{\epsilon}_j)$

Blocks with fewer than 3 samples, zero HI variance, or zero residual SD are excluded.

### Correlation Scanning

Pearson correlations are computed between all pairs of residualized, standardised blocks:

$$r_{AB} = \frac{1}{n-1} \sum_i Z_{Ai} \cdot Z_{Bi}$$

**Interchromosomal** (default): all pairs from different chromosomes.
**Intrachromosomal** (`--intra`): all pairs from the same chromosome, optionally within `--max-dist` bp.

Matrix multiplication is tiled (default tile size 1024) for cache efficiency with large block sets.

### Permutation Testing (Westfall-Young FWER)

Genome-wide significance thresholds are estimated by full-shuffle permutation:

1. Residualize all blocks on HI to produce the z-scored residual matrix Z (computed once)
2. Draw a single random permutation of sample indices per replicate; for every interchromosomal pair (A, B), apply it to the blocks of A while leaving B in natural sample order — the same permutation is reused across all pairs within a replicate
3. Scan all interchromosomal pairs on the permuted Z and record the distribution of r values
4. Repeat N times (`--permute N`)

**Interpretation**: The 95th percentile of `max_r` across replicates is a genome-wide threshold at α = 0.05 (Westfall-Young family-wise error rate). Equivalently, the 5th percentile of `min_r` gives a threshold for negative correlations.

A reservoir sampler with configurable capacity (`--permute-sample`) tracks the distribution without storing all correlations.

### LOCO Mode (`--hi-mode excl-focus`)

In leave-one-chromosome-out (LOCO) mode, the HI used to residualize blocks on chromosome $c$ excludes that chromosome's contribution:

$$\text{HI}^{-c}_i = \frac{\sum_{j \notin c} d_{ij} L_j}{2 \sum_{j \notin c} L_j}$$

This prevents the focal chromosome's own ancestry signal from inflating its residuals, which would suppress true signals on that chromosome while producing artifactual correlations between it and others. LOCO mode is recommended when the data contain strong chromosomal-scale ancestry blocks.

---

## Missing Data

By default (`--min-callrate 1.0`), AdmixLD excludes any block that has missing genotype data in any sample. All residualization and correlation calculations are therefore performed on fully observed data.

When `--min-callrate <value>` is set to a value less than 1.0:

1. Blocks with a call-rate (proportion of non-missing genotypes) ≥ the threshold are **retained**.
2. Remaining missing genotypes within retained blocks are **mean-imputed** using the block's per-block mean dosage.
3. Residualization and scanning proceed on the imputed matrix.

Mean imputation is applied only when explicitly requested. This option is intended primarily for reduced-representation datasets (e.g., RADseq) where missing data are common and a strict complete-data requirement would exclude most blocks. Note that mean imputation can attenuate correlation estimates.

---

## Usage Examples

### Compute and inspect the hybrid index

```bash
./build/admixld --vcf data.vcf --out results --compute-hi
# Writes results.hi.tsv and exits
```

### Basic interchromosomal scan

```bash
./build/admixld \
  --vcf data.vcf \
  --out results \
  --min-abs-r 0.5
```

### Intrachromosomal scan with distance limit

```bash
./build/admixld \
  --vcf data.vcf \
  --out results \
  --intra \
  --max-dist 5000000 \
  --min-abs-r 0.3
```

### Scan with pre-computed hybrid index

```bash
./build/admixld \
  --vcf data.vcf \
  --out results \
  --hi custom.hi.tsv \
  --min-abs-r 0.5
```

### LOCO (leave-one-chromosome-out) mode

```bash
./build/admixld \
  --vcf data.vcf \
  --out results \
  --hi-mode excl-focus \
  --min-abs-r 0.4
```

### Permutation test for genome-wide thresholds

```bash
./build/admixld \
  --vcf data.vcf \
  --out results \
  --permute 1000 \
  --seed 42 \
  --threads 8 \
  --permute-sample 200000
# Check results.perm.summary.tsv: p95 of max_r is the α=0.05 threshold
```

### Empirical distribution summary

```bash
./build/admixld \
  --vcf data.vcf \
  --out results \
  --distrib \
  --min-abs-r 0.0
```

### Focus block: correlate one block against the rest

```bash
./build/admixld \
  --vcf data.vcf \
  --out results \
  --target-chr scaffold-mi9 \
  --target-pos 421750
```

### Correlate a sample trait against all blocks

```bash
./build/admixld \
  --vcf data.vcf \
  --out results \
  --sample-geno phenotype.tsv
```

### Restrict to BED regions with relaxed call-rate

```bash
./build/admixld \
  --vcf data.vcf \
  --out results \
  --bed regions.bed \
  --min-callrate 0.9 \
  --min-abs-r 0.5
```

### Scan with directional filters (positive correlations only)

```bash
./build/admixld \
  --vcf data.vcf \
  --out results \
  --min-pos-r 0.5
```

---

## Notes and Caveats

- **Block count limit**: Hard-coded maximum of 300,000 blocks.
- **Sample size**: At least 3 samples with non-missing data are required for any block to be included.
- **Chromosome naming**: Chromosome names are matched exactly (case-sensitive) between VCF, BED, and any `--chr`/`--no-chr` arguments.
- **Multithreading**: Parallel scanning writes per-thread temporary files that are merged after the scan. Thread count has no effect on results, only runtime.
- **Weighted vs unweighted HI**: Weighted HI is preferred when blocks vary substantially in length (typical for phased local ancestry). Unweighted HI may be appropriate for evenly-spaced SNP data.
- **LOCO and permutations**: When using `--hi-mode excl-focus` with `--permute`, HI components are re-used across permutations for efficiency; only the HI mapping to shuffled samples changes.
