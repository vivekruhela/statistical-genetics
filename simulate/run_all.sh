#!/usr/bin/env bash
# Build the whole dataset, steps 0 to 9, into one directory (default: data).
#   TRACK=python   use the Python twins for steps 3, 5, 7 and 9 (default: r)
#   CHROMS="22"    limit genotypes to some chromosomes (for a quick trial run)
#   STATGEN_CONFIG another config file (default: config/cohort.json)
# Every step skips work whose output already exists, so a failed run can be restarted.
# Usage: bash simulate/run_all.sh [data]
set -euo pipefail
D=${1:-data}
TRACK=${TRACK:-r}
SAMPLES=data/eur_yri_samples.txt
step() {  # run the R or Python version of a step
  local s=$1; shift
  echo "== step ${s%%_*} ($TRACK): $s"
  if [[ $TRACK == r ]]; then Rscript "simulate/$s.R" "$@"; else python "simulate/$s.py" "$@"; fi
}
echo "== step 00: metadata";   bash simulate/00_fetch_metadata.sh "$SAMPLES" "$D/meta"
echo "== step 01: genotypes";  bash simulate/01_fetch_genotypes.sh "$SAMPLES" "$D/genotypes"
if [[ ! -s "$D/geuvadis/GD660.GeneQuantCount.txt.gz" ]]; then
  echo "== step 02: GEUVADIS"; bash simulate/02_fetch_geuvadis.sh "$D/geuvadis"
fi
step 03_prepare_expression "$D/geuvadis" "$D/meta" "$D/work/expression_real"
echo "== step 04: PLINK";      bash simulate/04_genotypes_to_plink.sh "$D/genotypes" "$D/work/expression_real/samples.txt" "$D/work/plink"
step 05_pick_loci "$D/work/expression_real" "$D/work/plink" "$D/work/loci"
echo "== step 06: loci";       bash simulate/06_extract_loci.sh "$D/work/loci/loci.tsv" "$D/work/plink" "$D/cohort/genotypes"
step 07_simulate_layers "$D/work" "$D/cohort" "$D/truth" "$D/meta"
echo "== step 08: templates";  bash simulate/08_fetch_gwas_templates.sh "$D/work/loci/loci.tsv" "$D/meta" "$D/cohort/genotypes/loci" "$D/work/gwas_templates"
step 09_gwas_cohort "$D/cohort" "$D/truth" "$D/work/gwas_templates" "$D/gwas"
echo "== done: cohort in $D/cohort, answer key in $D/truth, GWAS summary statistics in $D/gwas"
