#!/usr/bin/env bash
# Post R, analysis 2: imputation with a matched, a mismatched and a combined reference panel.
#
# Pretend our cohort was typed on a genotyping chip: keep only every 10th common variant at each
# teaching locus, then let Beagle 5 fill the rest back in from a reference panel. Three panels,
# each made of people who are NOT in the cohort (the GWAS template people from step 8):
#   EUR   every European-ancestry template person (about 258)
#   AFR   the same number of African-ancestry template people (the first ones in the list), so
#         the two panels differ in ancestry, not in size
#   BOTH  the two together
# 03_imputation_accuracy.{R,py} then compares the filled-in dosages with the real genotypes.
#   OUT/chip/Lxx.sites.txt     the variants kept on the pretend chip
#   OUT/imputed/Lxx.PANEL.ds.tsv.gz   imputed ALT dosage (DS) for every variant, every person
# Needs bcftools and beagle (bioconda: beagle) on the PATH.
# Usage: bash posts/R-reference-panels/02_imputation.sh data results/post-R/imputation
set -euo pipefail
D=${1:?data directory}
OUT=${2:?output directory}
THREADS=${THREADS:-2}
CHIP_EVERY=${CHIP_EVERY:-10}     # keep 1 in this many common variants on the "chip"
TMPL=$D/work/gwas_templates
mkdir -p "$OUT"/{chip,ref,target,imputed,log}

# a template haplotype table (id, then two 0/1 columns per person) as a phased VCF;
# $2 = first person, $3 = last person (1-based) to keep, $4 = chromosome name to write
hap_to_vcf() {
  zcat "$1" | awk -F'\t' -v OFS='\t' -v first="$2" -v last="$3" -v chrom="$4" '
    NR == 1 {
      print "##fileformat=VCFv4.2"
      print "##contig=<ID=" chrom ">"
      print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotype\">"
      printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
      for (p = first; p <= last; p++) { s = $(2 * p); sub(/_1$/, "", s); printf "\t%s", s }
      printf "\n"; next }
    { split($1, v, ":")
      printf "%s\t%s\t%s\t%s\t%s\t.\tPASS\t.\tGT", chrom, v[2], $1, v[3], v[4]
      for (p = first; p <= last; p++) printf "\t%s|%s", $(2 * p), $(2 * p + 1)
      printf "\n" }'
}

tail -n +2 "$D/cohort/loci.tsv" | while IFS=$'\t' read -r locus scenario chr rest; do
  chrom=${chr#chr}   # the cohort VCFs, written by plink2, call it "22", not "chr22"; the panel must match
  if [[ -s "$OUT/imputed/$locus.BOTH.ds.tsv.gz" ]]; then echo "$locus already done"; continue; fi
  # the pretend chip: every CHIP_EVERY-th variant with MAF >= 5% in the cohort
  # (2>/dev/null: plink2 VCFs have no contig lines, which bcftools warns about)
  bcftools view -q 0.05:minor -Ou "$D/cohort/genotypes/loci/$locus.vcf.gz" 2>/dev/null \
    | bcftools query -f '%CHROM\t%POS\t%ID\n' | awk -v k="$CHIP_EVERY" '(NR - 1) % k == 0' > "$OUT/chip/$locus.sites.txt"
  bcftools view -T "$OUT/chip/$locus.sites.txt" -Oz -o "$OUT/target/$locus.vcf.gz" "$D/cohort/genotypes/loci/$locus.vcf.gz" 2>/dev/null

  n_people() { { gzip -cd "$1" || true; } | awk -F'\t' 'NR == 1 {print (NF - 1) / 2; exit}'; }
  N_REF=$(( $(n_people "$TMPL/EUR/$locus.hap.gz") < $(n_people "$TMPL/AFR/$locus.hap.gz") ?
            $(n_people "$TMPL/EUR/$locus.hap.gz") : $(n_people "$TMPL/AFR/$locus.hap.gz") ))
  hap_to_vcf "$TMPL/EUR/$locus.hap.gz" 1 "$N_REF" "$chrom" | bgzip > "$OUT/ref/$locus.EUR.vcf.gz"
  hap_to_vcf "$TMPL/AFR/$locus.hap.gz" 1 "$N_REF" "$chrom" | bgzip > "$OUT/ref/$locus.AFR.vcf.gz"
  for p in EUR AFR; do bcftools index -f "$OUT/ref/$locus.$p.vcf.gz"; done
  bcftools merge -Oz -o "$OUT/ref/$locus.BOTH.vcf.gz" "$OUT/ref/$locus.EUR.vcf.gz" "$OUT/ref/$locus.AFR.vcf.gz"

  for panel in EUR AFR BOTH; do
    beagle ref="$OUT/ref/$locus.$panel.vcf.gz" gt="$OUT/target/$locus.vcf.gz" out="$OUT/imputed/$locus.$panel" \
      impute=true gp=false ap=false seed=20261121 nthreads="$THREADS" > "$OUT/log/$locus.$panel.log" 2>&1
    {
      bcftools query -l "$OUT/imputed/$locus.$panel.vcf.gz" 2>/dev/null | awk 'BEGIN {printf "id"} {printf "\t%s", $1} END {printf "\n"}'
      bcftools query -f '%ID[\t%DS]\n' "$OUT/imputed/$locus.$panel.vcf.gz" 2>/dev/null
    } | gzip > "$OUT/imputed/$locus.$panel.ds.tsv.gz"
    rm "$OUT/imputed/$locus.$panel.vcf.gz"
  done
  echo "$locus: $(wc -l < "$OUT/chip/$locus.sites.txt") chip variants; imputed from EUR and AFR panels of $N_REF people each, and both"
done
