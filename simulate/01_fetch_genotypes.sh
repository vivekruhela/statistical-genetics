#!/usr/bin/env bash
# Step 1: pull phased GRCh38 genotypes for the cohort from the public HGDP+1kG panel
# (gnomAD resource, SHAPEIT5-phased, built on the 1000 Genomes 30x calls).
# Keeps biallelic SNVs with minor allele frequency >= 1% in the cohort.
# Usage: CHROMS="22" bash simulate/01_fetch_genotypes.sh data/eur_yri_samples.txt data/genotypes
set -euo pipefail
SAMPLES=${1:?sample list, one ID per line}
OUT=${2:?output directory}
CHROMS=${CHROMS:-$(seq 1 22)}
THREADS=${THREADS:-2}
PANEL=https://storage.googleapis.com/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2
mkdir -p "$OUT"
for c in $CHROMS; do
  f="$OUT/chr$c.bcf"
  [[ -s "$f.csi" ]] && { echo "chr$c already done"; continue; }
  bcftools view --threads "$THREADS" -S "$SAMPLES" --force-samples -v snps -m2 -M2 -Ou \
      "$PANEL/hgdp1kgp_chr$c.filtered.SNV_INDEL.phased.shapeit5.bcf" 2>/dev/null \
    | bcftools annotate -x INFO,^FORMAT/GT -Ou \
    | bcftools view -q 0.01:minor -Ob -o "$f.tmp"
  mv "$f.tmp" "$f"
  bcftools index "$f"
  echo "chr$c $(bcftools index -n "$f") variants"
done
