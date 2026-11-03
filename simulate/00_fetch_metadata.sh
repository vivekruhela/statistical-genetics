#!/usr/bin/env bash
# Step 0: sample metadata for the HGDP+1kG panel and a GRCh38 gene annotation.
#   data/meta/panel_samples.tsv     every panel sample: population, region, sex, relatedness flag
#   data/meta/cohort_samples.tsv    the 654 EUR+YRI candidates with population and sex
#   data/meta/gwas_templates_*.txt  unrelated people NOT in the cohort, used in step 6 to build
#                                   the large GWAS cohort (so GWAS and eQTL come from different people)
#   data/meta/gencode.v39.annotation.gtf.gz   GENCODE v39 (GRCh38), mirrored by gnomAD on GCS
# Usage: bash simulate/00_fetch_metadata.sh data/eur_yri_samples.txt data/meta
set -euo pipefail
SAMPLES=${1:?sample list, one ID per line}
OUT=${2:?output directory}
GCS=https://storage.googleapis.com/gcp-public-data--gnomad
META_URL=$GCS/release/3.1/secondary_analyses/hgdp_1kg_v2/metadata_and_qc/gnomad_meta_updated.tsv
GTF_URL=$GCS/resources/grch38/gencode/gencode.v39.annotation.gtf.gz
mkdir -p "$OUT"
[[ -s "$OUT/gnomad_meta_updated.tsv" ]] || curl -fsSL --retry 3 -o "$OUT/gnomad_meta_updated.tsv" "$META_URL"
[[ -s "$OUT/gencode.v39.annotation.gtf.gz" ]] || curl -fsSL --retry 3 -o "$OUT/gencode.v39.annotation.gtf.gz" "$GTF_URL"

# Pick columns by name so a reordered file still works.
awk -F'\t' -v OFS='\t' '
  NR==1 { for (i=1;i<=NF;i++) col[$i]=i
          print "sample","project","population","region","sex","related"; next }
  { sex=$col["project_meta.sex"]
    if (sex!="male" && sex!="female") {
      k=$col["sex_imputation.sex_karyotype"]; sex=(k=="XY")?"male":(k=="XX")?"female":"NA" }
    rel=$col["sample_filters.release_related"]; if (rel=="") rel="NA"
    print $col["s"],$col["hgdp_tgp_meta.Project"],$col["hgdp_tgp_meta.Population"],
          $col["hgdp_tgp_meta.Genetic.region"],sex,rel }' \
  "$OUT/gnomad_meta_updated.tsv" > "$OUT/panel_samples.tsv"

awk -F'\t' -v OFS='\t' 'NR==FNR {want[$1]=1; next}
  FNR==1 {print; next} ($1 in want) {print}' "$SAMPLES" "$OUT/panel_samples.tsv" \
  > "$OUT/cohort_samples.tsv"

# GWAS templates: unrelated, not in the cohort, not admixed (ASW, ACB) and not the deeply
# diverged hunter-gatherer groups (San, Biaka, Mbuti), whose LD differs a lot from YRI.
awk -F'\t' 'NR==FNR {want[$1]=1; next} FNR>1 && !($1 in want) && $6=="false" && $4=="EUR" {print $1}' \
  "$SAMPLES" "$OUT/panel_samples.tsv" > "$OUT/gwas_templates_EUR.txt"
awk -F'\t' 'NR==FNR {want[$1]=1; next} FNR>1 && !($1 in want) && $6=="false" && $4=="AFR" &&
  $3!~/^(ASW|ACB|San|BiakaPygmy|MbutiPygmy)$/ {print $1}' \
  "$SAMPLES" "$OUT/panel_samples.tsv" > "$OUT/gwas_templates_AFR.txt"

echo "cohort candidates found in panel metadata: $(($(wc -l < "$OUT/cohort_samples.tsv")-1)) of $(wc -l < "$SAMPLES")"
echo "GWAS templates: EUR $(wc -l < "$OUT/gwas_templates_EUR.txt"), AFR $(wc -l < "$OUT/gwas_templates_AFR.txt")"
