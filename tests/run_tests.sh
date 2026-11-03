#!/usr/bin/env bash
# All checks. Usage: bash tests/run_tests.sh [data]
#  1. unit tests of lib/ in R and in Python
#  2. if the data are built: the planted truth is in the data (tests/test_data.py)
#  3. if the data are built: the Python steps 3 and 5 reproduce the R outputs exactly
# Needs both environments on PATH, or run the R and Python parts from each env.
set -euo pipefail
D=${1:-data}
echo "== unit tests (R)";      Rscript -e 'res <- testthat::test_file("tests/test_lib.R"); if (any(as.data.frame(res)$failed > 0)) quit(status = 1)'
echo "== unit tests (Python)"; python -m pytest -q tests/test_lib.py
if [[ -s "$D/truth/truth.tsv" ]]; then
  echo "== planted truth";     STATGEN_DATA="$D" python -m pytest -q tests/test_data.py
  echo "== R and Python give the same numbers (steps 3 and 5)"
  C="$D/work/py_check"
  python simulate/03_prepare_expression.py "$D/geuvadis" "$D/meta" "$C/expression_real" > /dev/null
  python tests/compare_tables.py "$D/work/expression_real" "$C/expression_real"
  python simulate/05_pick_loci.py "$D/work/expression_real" "$D/work/plink" "$C/loci" > /dev/null
  python tests/compare_tables.py "$D/work/loci" "$C/loci"
else
  echo "(no data in $D yet: skipped the data checks)"
fi
echo "== all checks passed"
