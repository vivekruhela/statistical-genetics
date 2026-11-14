#!/usr/bin/env bash
# Post R: every analysis and figure, in order. Steps skip work already done where they can.
#   TRACK=python   run the Python twins (default: r)
#   STATGEN_CONFIG the config the data were built with (default: config/cohort.json)
#   THREADS        threads for Beagle (default: 2)
# Needs the built data (simulate/run_all.sh), and on the PATH: Rscript or python, bcftools and
# beagle. The in-sample LD (step 4) has to be made by the language that built the GWAS cohort,
# so step 4 picks R or Python from data/gwas/EUR/summary.tsv, whatever TRACK says.
# Usage: bash posts/R-reference-panels/run_all.sh data results/post-R
set -euo pipefail
D=${1:-data}
OUT=${2:-results/post-R}
TRACK=${TRACK:-r}
P=posts/R-reference-panels
run() {  # run the R or Python version of a step
  local s=$1; shift
  echo "== $s ($TRACK)"
  if [[ $TRACK == r ]]; then Rscript "$P/$s.R" "$@"; else python "$P/$s.py" "$@"; fi
}
mkdir -p "$OUT"
run 01_ld_maps "$D" "$OUT"
echo "== 02_imputation (Beagle)"
bash "$P/02_imputation.sh" "$D" "$OUT/imputation"
run 03_imputation_accuracy "$D" "$OUT/imputation" "$OUT"
built=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == "built_with") c = i} NR == 2 {print $c}' "$D/gwas/EUR/summary.tsv")
echo "== 04_gwas_ld (the GWAS cohort was built with $built)"
if [[ $built == R ]]; then Rscript "$P/04_gwas_ld.R" "$D" "$OUT"; else python "$P/04_gwas_ld.py" "$D" "$OUT"; fi
run 05_finemap_panels "$D" "$OUT" "${SHOW:-L05_trait_q}"
echo "== done: tables and figures in $OUT"
