#!/usr/bin/env bash
# Post R: run the R track, then the Python twins on the same inputs, and check that every table
# matches (to 1e-6). The Python track reuses the R track's Beagle runs and in-sample LD, which
# are the same inputs for both. Needs Rscript and python (and bcftools, beagle) on the PATH.
# Usage: bash posts/R-reference-panels/check.sh data results/post-R-check
set -euo pipefail
D=${1:-data}
OUT=${2:-results/post-R-check}
P=posts/R-reference-panels
TRACK=r bash "$P/run_all.sh" "$D" "$OUT/r"
mkdir -p "$OUT/py/ld"
cp "$OUT/r/windows.tsv" "$OUT/py/"
ln -sfn "$(cd "$OUT/r/ld/insample" && pwd)" "$OUT/py/ld/insample"
python "$P/01_ld_maps.py" "$D" "$OUT/py"
python "$P/03_imputation_accuracy.py" "$D" "$OUT/r/imputation" "$OUT/py"
python "$P/05_finemap_panels.py" "$D" "$OUT/py" "${SHOW:-L05_trait_q}"
echo "== R and Python give the same numbers"
python tests/compare_tables.py "$OUT/r" "$OUT/py" ld_decay.tsv ld_map.tsv imputation_summary.tsv \
  imputation_r2.tsv.gz finemap_summary.tsv finemap_sets.tsv finemap_variants.tsv.gz
