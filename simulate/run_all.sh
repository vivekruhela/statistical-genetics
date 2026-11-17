#!/usr/bin/env bash
# Build the whole dataset, steps 0 to 9, into one directory (default: data).
#   TRACK=python   use the Python twins for steps 3, 5, 7 and 9 (default: r)
#   CHROMS="22"    limit genotypes to some chromosomes (for a quick trial run)
#   STATGEN_CONFIG another config file (default: config/cohort.json)
#   FROM=05        start at this step, keeping what steps before it already built
#                  (e.g. after step 5 was killed: FROM=05 bash simulate/run_all.sh data)
# Steps 1, 2 and 4 also skip downloads and files that already exist, so a plain rerun is safe too.
# The full build is heavy (step 5 scans ~6.9 million variants; it peaks around 3 GB, PLINK in
# step 4 takes more): run it on a compute node, not a login/head node, which may kill it with a
# bare "Terminated". simulate/run_all.sbatch is a ready-made Slurm job.
# Usage: bash simulate/run_all.sh [data]
set -euo pipefail
D=${1:-data}
TRACK=${TRACK:-r}
FROM=${FROM:-00}
SAMPLES=data/eur_yri_samples.txt
todo() { (( 10#$1 >= 10#$FROM )); }  # is step $1 at or after FROM?
step() {  # run the R or Python version of a step
  local s=$1; shift
  todo "${s%%_*}" || return 0
  echo "== step ${s%%_*} ($TRACK): $s"
  if [[ $TRACK == r ]]; then Rscript "simulate/$s.R" "$@"; else python "simulate/$s.py" "$@"; fi
}
todo 00 && { echo "== step 00: metadata";   bash simulate/00_fetch_metadata.sh "$SAMPLES" "$D/meta"; }
todo 01 && { echo "== step 01: genotypes";  bash simulate/01_fetch_genotypes.sh "$SAMPLES" "$D/genotypes"; }
if todo 02 && [[ ! -s "$D/geuvadis/GD660.GeneQuantCount.txt.gz" ]]; then
  echo "== step 02: GEUVADIS"; bash simulate/02_fetch_geuvadis.sh "$D/geuvadis"
fi
step 03_prepare_expression "$D/geuvadis" "$D/meta" "$D/work/expression_real"
todo 04 && { echo "== step 04: PLINK";      bash simulate/04_genotypes_to_plink.sh "$D/genotypes" "$D/work/expression_real/samples.txt" "$D/work/plink"; }
step 05_pick_loci "$D/work/expression_real" "$D/work/plink" "$D/work/loci"
todo 06 && { echo "== step 06: loci";       bash simulate/06_extract_loci.sh "$D/work/loci/loci.tsv" "$D/work/plink" "$D/cohort/genotypes"; }
step 07_simulate_layers "$D/work" "$D/cohort" "$D/truth" "$D/meta"
todo 08 && { echo "== step 08: templates";  bash simulate/08_fetch_gwas_templates.sh "$D/work/loci/loci.tsv" "$D/meta" "$D/cohort/genotypes/loci" "$D/work/gwas_templates"; }
step 09_gwas_cohort "$D/cohort" "$D/truth" "$D/work/gwas_templates" "$D/gwas"
echo "== done: cohort in $D/cohort, answer key in $D/truth, GWAS summary statistics in $D/gwas"
