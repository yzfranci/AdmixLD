#!/usr/bin/env bash
set -euo pipefail

BIN=./build/adfinder
VCF=tests/data/ancestry_blocks.vcf

OUT=$(mktemp -d)

$BIN \
	--vcf "$VCF" \
	--out "$OUT/test" \
	--max-windows 200 \
	--block-size 64 \
	--intra \
	--min-abs-r 0.2 \
	> "$OUT/log.txt"

# Extract numbers
TESTED=$(awk '/tested_pairs/ {for (i=1;i<=NF;i++) if ($i ~ /^[0-9]+$/) {print $i; exit}}' "$OUT/log.txt")
KEPT=$(awk '/kept_pairs/ {for (i=1;i<=NF;i++) if ($i ~ /^[0-9]+$/) {print $i; exit}}' "$OUT/log.txt")

# Expected values (you fill these once)
EXPECTED_TESTED=19900
EXPECTED_KEPT=9288

if [[ "$TESTED" != "$EXPECTED_TESTED" ]]; then
	echo "FAILED: tested_pairs = $TESTED (expected $EXPECTED_TESTED)"
	exit 1
fi

if [[ "$KEPT" != "$EXPECTED_KEPT" ]]; then
	echo "FAILED: kept_pairs = $KEPT (expected $EXPECTED_KEPT)"
	exit 1
fi

echo "OK: regression test passed"
