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

### External HI File (`--cov`)

Optional user-supplied hybrid index used instead of the computed HI for residualization. Tab-separated, two columns: `sample` and `hi`. A header row is auto-detected. Must contain all samples present in the input file. Exactly one value column is required.

```
sample    hi
CV0340    0.373028
CV0359    0.213722
CV0362    0.891145
```
### Sample Haplotype TSV (`--sample-haplo`)

Per-sample mitochondrial haplotype for correlating against all ancestry markers. Tab-separated, two columns: `sample` and haplotype value. Non-missing values must be **exactly `0` or `1`**; any other value is an error. Missing values can be encoded as `.`, `NA`, `NaN`, `nan`, or `NAN`.

```
sample    haplo
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

### Parental Allele Frequency File (`--ref-freq`)

A tab-separated file providing the allele frequency of the ancestry-1 allele (`p1`) and the ancestry-2 allele (`p2`) in the two parental populations, one row per marker. A header row is auto-detected.

```
chrom           pos       p1      p2
scaffold-ma1    10000     0.95    0.02
scaffold-ma1    25000     0.88    0.11
```

- `chrom` and `pos` must match a record in the VCF exactly (chromosome name and integer position).
- `p1` and `p2` are allele frequencies in [0, 1].
- Markers present in the VCF but absent from this file are dropped from the analysis (a count is reported).
- Only used with `--vcf`; ignored (with a warning) when `--msp` is provided.

See [Frequency-Based HI Computation](#frequency-based-hi-computation---ref-freq) for details on how these frequencies are used.

### Sample Keep List (`--keep-indv`)

A plain-text file listing one sample ID per line. Only these samples are included in the **correlation scan**. HI estimation, LOCO component building, and marker residualization are performed on the full sample set; the filter is applied afterwards so that it does not bias covariate estimation.

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
| `--ref-freq FILE` | none | Parental allele frequency TSV (`chrom`, `pos`, `p1`, `p2`); enables frequency-based HI and allele polarization (VCF only) |
| `--min-delta-afd FLOAT` | 0.0 | Minimum allele frequency difference \|p1 − p2\| to retain a marker when `--ref-freq` is used |
| `--hi-mode STR` | `global` | HI correction mode: `global` or `excl-focus` (LOCO; incompatible with `--cov`) |
| `--compute-hi` | off | Compute and save HI, then exit without scanning (ignored with `--cov`) |
| `--unweighted-hi` | off | Use unweighted (simple mean) HI instead of length-weighted (ignored with `--cov` and `--ref-freq`) |
| `--pos-is-start` | off | Treat VCF POS as marker start rather than marker end; no effect for `--msp` input |

LOCO (leave-one-chromosome-out) can be used to compute HI without the focus chromosome(s), to avoid confounding the focal chromosome's signal with its own contribution to the global HI estimate (Schumer and Brandvain, 2016).

### Filtering

| Flag | Default | Description |
|------|---------|-------------|
| `--chr STR` | all | Keep only this chromosome (repeatable) |
| `--no-chr STR` | none | Exclude this chromosome (repeatable) |
| `--bed FILE` | none | Keep only markers within BED intervals |
| `--min-callrate FLOAT` | 1.0 | Minimum call-rate to retain a marker; enables mean imputation if < 1.0 |
| `--keep-indv FILE` | all | Keep only samples listed in FILE (one ID per line) for the correlation scan; residualization uses the full set |

Be sure to use --min-callrate for loci that still have a decent amount of data (<10% missing), as it can introduce bias due to mean imputation.

### Target Marker

Used to scan correlations between a marker of interest against the rest of the markers.

| Flag | Description |
|------|-------------|
| `--target-chr STR` | Chromosome of single focus marker |
| `--target-pos INT` | Position of single focus marker (exact match); reports that marker vs all others |

### Sample Haplotype

Typically used to provide a mitochondrial haplotype that would not be contained in the VCF/MSP file.

| Flag | Description |
|------|-------------|
| `--sample-haplo FILE` | Per-sample mitochondrial haplotype TSV (0 or 1 only); correlates against all markers |

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

### Frequency-Based HI Computation (`--ref-freq`)

When `--ref-freq` is supplied together with `--vcf`, AdmixLD uses parental allele frequencies to estimate HI via a δ²-weighted estimator instead of the default dosage-average approach.

**Allele polarization**

Before computing HI, each marker is polarized so that p1 ≥ p2 (i.e., so that the ancestry-1 allele is always the higher-frequency allele in parental population 1). For any marker where p2 > p1 in the input file, the dosage column is flipped:

$$d'_{ij} = d_{\max} - d_{ij}$$

and p1 and p2 are swapped, where $d_{\max}$ is 2 for diploid and 1 for phased haplotype data. This step ensures that δ = p1 − p2 > 0 for every retained marker.

Markers absent from the `--ref-freq` file or with |δ| < `--min-delta-afd` are dropped before HI estimation (a count is printed to stdout).

**Formula**

$$\text{HI}_i = \frac{\displaystyle\sum_j \left(\frac{d_{ij}}{2} - p_{2j}\right) \delta_j}{\displaystyle\sum_j \delta_j^2}$$

where $\delta_j = p_{1j} - p_{2j}$ is the allele frequency difference. The estimator projects the mean-centred dosage onto the δ axis and weights each SNP by its informativeness (δ²), yielding a value in [0, 1] for pure-ancestry individuals. Each SNP contributes equally regardless of its genomic position; no segment-length weighting is applied.

**LOCO with `--ref-freq`**

When `--hi-mode excl-focus` is used together with `--ref-freq`, per-chromosome numerator and denominator components are accumulated separately, so the LOCO HI for chromosome $c$ can be computed without re-scanning all markers:

$$\text{HI}^{-c}_i = \frac{\displaystyle\sum_{j \notin c} \left(\frac{d_{ij}}{2} - p_{2j}\right) \delta_j}{\displaystyle\sum_{j \notin c} \delta_j^2}$$

### Residualization

Each marker's dosage vector is residualized on the hybrid index $h$ using a **plug-in** formula derived from the expected dosage under Hardy–Weinberg equilibrium in an admixed population. The slope relating dosage to ancestry is fixed by the biological model rather than estimated by regression, making the method more robust at small sample sizes.

**Without `--ref-freq`** (markers assumed ancestry-diagnostic, $p_1 \approx 1$, $p_2 \approx 0$):

$$\epsilon_{iw} = d_{iw} - 2 h_i \quad \text{(diploid)} \qquad \epsilon_{iw} = d_{iw} - h_i \quad \text{(phased)}$$

**With `--ref-freq`** (per-marker parental frequencies known):

$$\epsilon_{iw} = d_{iw} - 2\,(p_{2w} + \delta_w \cdot h_i) \quad \text{(diploid)} \qquad \epsilon_{iw} = d_{iw} - (p_{2w} + \delta_w \cdot h_i) \quad \text{(phased)}$$

where $\delta_w = p_{1w} - p_{2w}$ after polarization ($\delta_w \geq 0$ always). In both cases the residuals are then mean-centred and standardised:

$$Z_{iw} = \frac{\epsilon_{iw} - \bar{\epsilon}_w}{\text{SD}(\boldsymbol{\epsilon}_w)}$$

Markers with zero residual SD are excluded as uninformative.

**For `--sample-haplo`** (mitochondrial haplotype, haploid 0/1): the plug-in residual is $r_i = y_i - h_i$ (scale=1), then z-scored over non-missing entries.

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

When used with `--ref-freq`, the same exclusion principle applies using the frequency-based components; see [Frequency-Based HI Computation](#frequency-based-hi-computation---ref-freq).

---

## Missing Data

By default (`--min-callrate 1.0`), AdmixLD excludes any marker that has missing genotype data in any sample. All residualization and correlation calculations are therefore performed on fully observed data.

When `--min-callrate <value>` is set to a value less than 1.0:

1. Markers with a call-rate (proportion of non-missing genotypes) ≥ the threshold are **retained**.
2. Remaining missing genotypes within retained markers are **mean-imputed** using the per-marker mean dosage.
3. Residualization and scanning proceed on the imputed matrix.

Mean imputation is applied only when explicitly requested. This option is intended primarily for reduced-representation datasets (e.g., RADseq) where a strict complete-data requirement would exclude most markers. Note that mean imputation can attenuate correlation estimates.

---

## Notes and Caveats

- **Multithreading**: Parallel scanning writes per-thread temporary files that are merged after the scan. Thread count has no effect on results, only runtime.
- **Weighted vs unweighted HI**: Weighted HI is preferred when markers vary substantially in inferred segment length (typical for local ancestry tract data obtained via Loter, RFMIX etc...). Unweighted HI may be appropriate for SNP arrays.
- **Frequency-based HI (`--ref-freq`)**: This mode is recommended when markers are individual SNPs rather than ancestry tracts, since the δ²-weighting emphasises informative diagnostic markers. It requires a parental allele frequency file and is only compatible with `--vcf`.
- **Allele polarization**: Polarization is applied internally before HI computation and residualization. The dosage values stored in memory and used throughout the analysis reflect the polarized allele. The input VCF is not modified.
- **`--ref-freq` and `--msp`**: `--ref-freq` is silently ignored when `--msp` is provided. MSP dosages are already coded as ancestry counts and do not require frequency-based polarization.
- **LOCO and permutations**: When using `--hi-mode excl-focus` with `--permute`, HI components are re-used across permutations for efficiency; only the mapping of HI to shuffled samples changes.
- **`--msp` vs `--vcf`**: MSP input provides explicit tract start/end positions, so segment lengths for weighted HI are exact. VCF input infers lengths from consecutive marker positions (affected by `--pos-is-start`).
- **`--cov` and LOCO**: `--cov` is incompatible with `--hi-mode excl-focus`. LOCO requires per-chromosome HI components derived from the marker data itself and cannot be computed for externally supplied covariates.
- **`--cov` format**: `--cov` now accepts exactly one value column (`sample<TAB>hi`). Files with multiple numeric columns are rejected with an error.
- **Plug-in residualization**: AdmixLD uses a model-based plug-in formula (expected dosage = $2 h_i$ or $2(p_{2w} + \delta_w h_i)$) rather than OLS regression. This eliminates per-marker slope estimation and is more robust when sample sizes are small or markers are not perfectly ancestry-informative.
- **`--sample-haplo` values**: Non-missing values must be exactly 0 or 1. Any other value (including 0.5) is rejected at load time. Use missing-value tokens (`.`, `NA`) for samples without haplotype data.
- **`--keep-indv` and residualization**: The sample filter is applied **after** residualization. HI estimation, LOCO component building, and OLS residualization all use the full sample set to avoid biasing covariate estimation with a non-representative subset. Only the correlation scan (and permutations) use the filtered subset.
