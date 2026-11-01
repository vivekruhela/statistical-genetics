#!/usr/bin/env bash
# Step 8: phased haplotypes of the GWAS template people at each teaching locus.
# The large GWAS cohort (step 9) is built by copying segments of these haplotypes, so it has
# realistic LD at exactly the cohort's variants but contains none of the cohort's people.
#   OUT/POP/Lxx.hap.gz   one row per cohort variant (id), two 0/1 columns per template person
# Usage: bash simulate/08_fetch_gwas_templates.sh data/work/loci/loci.tsv data/meta data/cohort/genotypes/loci data/work/gwas_templates
set -euo pipefail
LOCI=${1:?loci.tsv from step 5}
META=${2:?metadata directory from step 0}
LOCUS_GENO=${3:?locus genotypes from step 6 (for the variant list)}
OUT=${4:?output directory}
THREADS=${THREADS:-2}
PANEL=https://storage.googleapis.com/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2
for pop in EUR AFR; do
  mkdir -p "$OUT/$pop"
  tail -n +2 "$LOCI" | while IFS=$'\t' read -r locus scenario chr start end rest; do
    f="$OUT/$pop/$locus.hap.gz"
    [[ -s "$f" ]] && { echo "$pop $locus already done"; continue; }
    c=${chr#chr}
    grep -v '^#' "$LOCUS_GENO/$locus.pvar" | cut -f3 > "$OUT/$pop/$locus.ids"
    bcftools view --threads "$THREADS" -S "$META/gwas_templates_$pop.txt" --force-samples \
        -r "$chr:$start-$end" -v snps -m2 -M2 -Ob -o "$OUT/$pop/$locus.tmp.bcf" \
        "$PANEL/hgdp1kgp_chr$c.filtered.SNV_INDEL.phased.shapeit5.bcf" 2>/dev/null
    {
      # header from the subset file itself, so columns follow its sample order
      bcftools query -l "$OUT/$pop/$locus.tmp.bcf" 2>/dev/null \
        | awk 'BEGIN {printf "id"} {printf "\t%s_1\t%s_2", $1, $1} END {printf "\n"}'
      bcftools query -f '%CHROM:%POS:%REF:%ALT[\t%GT]\n' "$OUT/$pop/$locus.tmp.bcf" 2>/dev/null \
        | awk -F'\t' 'NR==FNR {keep[$1]=1; next} ($1 in keep)' "$OUT/$pop/$locus.ids" - \
        | tr '|' '\t'
    } | gzip > "$f.tmp"
    rm "$OUT/$pop/$locus.tmp.bcf"
    mv "$f.tmp" "$f"
    rm "$OUT/$pop/$locus.ids"
    echo "$pop $locus: $(($(zcat "$f" | wc -l)-1)) of $(grep -vc '^#' "$LOCUS_GENO/$locus.pvar") cohort variants, $(($(zcat "$f" | head -1 | awk -F'\t' '{print NF}')-1)) haplotypes"
  done
done
