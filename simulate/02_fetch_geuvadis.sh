#!/usr/bin/env bash
# Step 2: GEUVADIS LCL RNA-seq (ArrayExpress E-GEUV-1): gene counts, RPKM and sample sheet.
# NOT YET TESTED: written from the E-GEUV-1 file names; check them on
# https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-GEUV-1 if a download fails.
set -euo pipefail
OUT=${1:?output directory}
BASE=https://www.ebi.ac.uk/biostudies/files/E-GEUV-1
mkdir -p "$OUT"
for f in E-GEUV-1.sdrf.txt \
         GD660.GeneQuantCount.txt.gz \
         GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz; do
  [[ -s "$OUT/$f" ]] || curl -fL --retry 3 -o "$OUT/$f" "$BASE/$f"
  echo "got $f"
done
