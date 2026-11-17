#!/usr/bin/env bash
# All checks. Usage: bash tests/run_tests.sh [data]
#  1. unit tests of lib/ in R and in Python
#  2. if the data are built: the planted truth is in the data (tests/test_data.py)
#  3. if the data are built: the Python steps 3 and 5 reproduce the R outputs exactly
# R and Python come from different conda envs (statgen-r, statgen-py). Run it from either one:
# whichever language is missing from the active env is run through `conda run -n <env>`.
# Override with RSCRIPT="..." and PYTHON="..." if your envs have other names.
# The data checks rerun the genome-wide scan of step 5 in Python: run on a compute node.
set -euo pipefail
D=${1:-data}

if [[ -z ${PYTHON:-} ]]; then
  if python -c "import pytest, numpy, scipy, pandas, statsmodels" 2>/dev/null; then PYTHON=python
  elif command -v conda > /dev/null && conda run -n statgen-py python -c "import pytest" 2>/dev/null; then
    PYTHON="conda run --no-capture-output -n statgen-py python"
  else
    echo "error: no Python with pytest, numpy, scipy, pandas and statsmodels." >&2
    echo "Create the env (mamba env create -f env-python.yml) or set PYTHON=/path/to/python." >&2
    exit 1
  fi
fi
if [[ -z ${RSCRIPT:-} ]]; then
  if command -v Rscript > /dev/null && Rscript -e 'library(testthat); library(data.table)' > /dev/null 2>&1; then RSCRIPT=Rscript
  elif command -v conda > /dev/null && conda run -n statgen-r Rscript -e 'library(testthat)' > /dev/null 2>&1; then
    RSCRIPT="conda run --no-capture-output -n statgen-r Rscript"
  else
    echo "error: no Rscript with testthat and data.table." >&2
    echo "Create the env (mamba env create -f env-r.yml) or set RSCRIPT=/path/to/Rscript." >&2
    exit 1
  fi
fi
echo "(R: $RSCRIPT; Python: $PYTHON)"

echo "== unit tests (R)";      $RSCRIPT -e 'res <- testthat::test_file("tests/test_lib.R"); if (any(as.data.frame(res)$failed > 0)) quit(status = 1)'
echo "== unit tests (Python)"; $PYTHON -m pytest -q tests/test_lib.py
if [[ -s "$D/truth/truth.tsv" ]]; then
  echo "== planted truth";     STATGEN_DATA="$D" $PYTHON -m pytest -q tests/test_data.py
  echo "== R and Python give the same numbers (steps 3 and 5)"
  C="$D/work/py_check"
  $PYTHON simulate/03_prepare_expression.py "$D/geuvadis" "$D/meta" "$C/expression_real" > /dev/null
  $PYTHON tests/compare_tables.py "$D/work/expression_real" "$C/expression_real"
  $PYTHON simulate/05_pick_loci.py "$D/work/expression_real" "$D/work/plink" "$C/loci" > /dev/null
  $PYTHON tests/compare_tables.py "$D/work/loci" "$C/loci"
else
  echo "(no data in $D yet: skipped the data checks)"
fi
echo "== all checks passed"
