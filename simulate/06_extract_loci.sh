#!/usr/bin/env bash
# Step 6: the released genotypes. Genome-wide PLINK 2 files and genotype PCs for the cohort,
# and the cohort cut to each teaching locus (all variants with MAF >= 1%).
#   OUT/genome/cohort.{pgen,pvar,psam}, OUT/genome/cohort_pca.eigenvec
#   OUT/loci/Lxx.{pgen,pvar,psam}  PLINK 2
#   OUT/loci/Lxx.vcf.gz            bgzipped VCF, for tools that want VCF
#   OUT/loci/Lxx.traw              ALT dosages, read by steps 7 and 9 and by the posts
# Usage: bash simulate/06_extract_loci.sh data/work/loci/loci.tsv data/work/plink data/cohort/genotypes
set -euo pipefail
LOCI=${1:?loci.tsv from step 5}
PLINK=${2:?PLINK directory from step 4}
OUT=${3:?output directory}
mkdir -p "$OUT/genome" "$OUT/loci"
cp "$PLINK"/cohort.{pgen,pvar,psam} "$PLINK"/cohort_pca.{eigenvec,eigenval} "$OUT/genome/"
P="plink2 --threads ${THREADS:-2} --silent"
tail -n +2 "$LOCI" | while IFS=$'\t' read -r locus scenario chr start end rest; do
  c=${chr#chr}
  $P --pfile "$PLINK/chr$c" --chr "$chr" --from-bp "$start" --to-bp "$end" --make-pgen --out "$OUT/loci/$locus"
  $P --pfile "$OUT/loci/$locus" --export vcf bgz id-paste=iid --out "$OUT/loci/$locus"
  $P --pfile "$OUT/loci/$locus" --export A-transpose --out "$OUT/loci/$locus"
  echo "$locus ($scenario) $chr:$start-$end: $(grep -vc '^#' "$OUT/loci/$locus.pvar") variants"
done
