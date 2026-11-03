#!/usr/bin/env bash
# Step 4: genotypes for the final cohort in PLINK 2 format, genotype PCs, and per-chromosome
# dosage tables for the cis-eQTL scan in step 5.
#   OUT/chrN.{pgen,pvar,psam}   per chromosome, variant IDs chr:pos:ref:alt
#   OUT/cohort.{pgen,pvar,psam} all chromosomes merged
#   OUT/cohort_pca.eigenvec     genotype PCs from LD-pruned common variants
#   OUT/scan/chrN.traw          ALT dosages, MAF >= scan_maf (config), for step 5
# Usage: bash simulate/04_genotypes_to_plink.sh data/genotypes data/work/expression_real/samples.txt data/work/plink
set -euo pipefail
GENO=${1:?directory with chrN.bcf from step 1}
SAMPLES=${2:?final sample list from step 3}
OUT=${3:?output directory}
CHROMS=${CHROMS:-$(seq 1 22)}
THREADS=${THREADS:-2}
# settings from the config (honours STATGEN_CONFIG and "_base"), via the R helpers
cfg() { Rscript -e "suppressWarnings(source('lib/statgen.R')); cat(read_config()\$genotype\$$1)" 2>/dev/null; }
MAF=$(cfg scan_maf)
NPC=$(cfg n_pcs)
[[ -n "$MAF" && -n "$NPC" ]] || { echo "could not read the config" >&2; exit 1; }
mkdir -p "$OUT/scan"
P="plink2 --threads $THREADS --silent"
: > "$OUT/merge_list.txt"
for c in $CHROMS; do
  if [[ ! -s "$OUT/chr$c.pgen" ]]; then
    $P --bcf "$GENO/chr$c.bcf" --keep "$SAMPLES" --output-chr chrM \
       --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 200 truncate \
       --make-pgen --out "$OUT/chr$c"
  fi
  echo "$OUT/chr$c" >> "$OUT/merge_list.txt"
  [[ -s "$OUT/scan/chr$c.traw" ]] || \
    $P --pfile "$OUT/chr$c" --maf "$MAF" --export A-transpose --out "$OUT/scan/chr$c"
  echo "chr$c: $(grep -vc '^#' "$OUT/chr$c.pvar") variants, $(($(wc -l < "$OUT/scan/chr$c.traw")-1)) for the scan"
done

if [[ $(wc -l < "$OUT/merge_list.txt") -gt 1 ]]; then
  $P --pmerge-list "$OUT/merge_list.txt" pfile --make-pgen --out "$OUT/cohort"
else
  $P --pfile "$(cat "$OUT/merge_list.txt")" --make-pgen --out "$OUT/cohort"
fi
$P --pfile "$OUT/cohort" --maf 0.05 --indep-pairwise 1000kb 0.1 --out "$OUT/prune"
$P --pfile "$OUT/cohort" --extract "$OUT/prune.prune.in" --pca "$NPC" --out "$OUT/cohort_pca"
echo "PCA on $(wc -l < "$OUT/prune.prune.in") pruned variants -> $OUT/cohort_pca.eigenvec"
