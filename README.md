# AdmixLD Documentation

## Overview

AdmixLD is a command-line tool for detecting **ancestry disequilibrium (AD)** in hybrid zones and admixed populations. It scans for statistically significant pairwise correlations between residualized ancestry markers — correlations that remain after conditioning on global hybrid index (ancestry proportion). Such signals can indicate adaptive introgression, selection, or anomalous population structure at specific genomic loci.

**Core idea**: Each ancestry marker's dosage is residualized on the genome-wide hybrid index (HI), removing background ancestry-linkage disequilibrium. Pearson correlations between residualized markers across the genome then reveal statistically significant deviations.

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

The binary is placed at `build/admixld`.

---

### Input files

### Usage
admixld (load blocks into matrix)
Usage:
  admixld --vcf input.vcf[.gz] --out output_prefix
  admixld --msp input.msp.tsv   --out output_prefix
  --min-abs-r FLOAT      Keep pairs with |r| >= value (default)
  --min-neg-r FLOAT      Keep pairs with r <= -value (asymmetric)
  --min-pos-r FLOAT      Keep pairs with r >= value (asymmetric)
  --intra                Scan intrachromosomal pairs (default: interchrom only)
  --phased               Use phased haplotypes (2n rows); requires --intra; MSP/VCF GT only
  --max-dist INT         Intra only: max end-distance (bp) between block pairs
  --max-windows N        Load at most N blocks from the VCF (used solely for test runs and debugging)
  --min-callrate FLOAT   Minimum block call rate; values <1.0 enable within-block mean imputation (beta; might introduce bias).
  --cov FILE             Covariate file (TSV: for HI: sample<TAB>hi; For more covariates (PCA etc...) sample<TAB>cov1 [cov2 ...]; header optional); overrides computed HI
  --hi-mode STR          HI correction: global (default) | excl-focus (LOCO/LOCO2; incompatible with --cov)
  --compute-hi           Compute HI using FILTERED blocks only, write out.hi.tsv, and exit (ignored with --cov)
  --unweighted-hi        Use unweighted HI (mean(dosage)/2; legacy behavior; ignored with --cov)
  --pos-is-start         For weighted HI: interpret VCF block pos as START (ignored with --cov)
  --tile-size INT                tile size for processing (default: 1024)
  --threads INT          Number of OpenMP threads for scan steps (default: 1)
  --distrib              Write empirical scan distribution summary
  --distrib-sample INT   Reservoir sample size for distribution summary (default: 200000)
  --permute N            Run N interchrom full-shuffle permutations (summary stats)
  --permute-sample INT   Reservoir sample size for percentile estimates (default: 200000)
  --seed INT             RNG seed for permutations (default: 1)
  --chr STR              Keep only this chromosome (repeatable)
  --no-chr STR           Exclude this chromosome (repeatable; opposite of --chr)
  --bed FILE             Keep blocks whose position is within BED intervals (chr start end; no header)
  --target-chr STR       Scan one target block/pos vs all others (target chromosome)
  --target-pos INT       Scan one target block/pos vs all others (target position; matches single position)
  --sample-geno FILE     Per-sample numeric mito haplotype/genotype/trait TSV: sample<TAB>value
  --keep-indv FILE       Keep only samples listed in FILE (one ID per line)

### Outputs