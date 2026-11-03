#!/usr/bin/env bash
# Step 2: GEUVADIS LCL RNA-seq (ArrayExpress / BioStudies E-GEUV-1): sample sheet and gene counts.
# The sdrf sits at the top of the study; the processed files are under E-GEUV-1/analysis_results/.
# Checked on 3 Nov 2026: GD660.GeneQuantCount.txt.gz is 51,372,443 bytes, 660 libraries.
# Usage: bash simulate/02_fetch_geuvadis.sh data/geuvadis
set -euo pipefail
OUT=${1:?output directory}
BASE=https://ftp.ebi.ac.uk/biostudies/fire/E-GEUV-/001/E-GEUV-1/Files
mkdir -p "$OUT"
# file, path under BASE, md5 (from the study's processed-data_filelist.tsv; "-" = not checked)
while read -r f path md5; do
  if [[ ! -s "$OUT/$f" ]]; then
    curl -fL --retry 3 -o "$OUT/$f.tmp" "$BASE/$path"
    mv "$OUT/$f.tmp" "$OUT/$f"
  fi
  if [[ "$md5" != "-" ]]; then
    got=$(md5sum "$OUT/$f" | cut -d' ' -f1)
    [[ "${got,,}" == "${md5,,}" ]] || { echo "md5 mismatch for $f: delete it and rerun" >&2; exit 1; }
  fi
  echo "got $f"
done <<'LIST'
E-GEUV-1.sdrf.txt E-GEUV-1.sdrf.txt -
GD660.GeneQuantCount.txt.gz E-GEUV-1/analysis_results/GD660.GeneQuantCount.txt.gz 9F91E5548BEF1CD703E53B773FA9D628
LIST
