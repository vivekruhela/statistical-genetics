# Data

Everything here except `eur_yri_samples.txt` is built by `simulate/run_all.sh` and is not in git.
The finished bundle will be published on Zenodo with a DOI.

## What a build contains

```
cohort/
  samples.tsv                   sample, population, superpopulation, sex
  clinical/clinical.tsv         age, bmi, trait_q, disease (0/1)
  expression/counts.tsv.gz      gene counts, real GEUVADIS genes plus SIMG01-SIMG10
  expression/expression_int.tsv.gz      inverse-normal expression (filtered genes)
  expression/expression_log2norm.tsv.gz log2(size-factor-normalised count + 1)
  expression/genes.tsv, size_factors.tsv, expression_pcs.tsv
  methylation/methylation_m.tsv.gz      CpG M-values (beta = 2^M / (2^M + 1))
  methylation/cpgs.tsv          CpG positions; promoter_of names the teaching gene
  genotypes/genome/             genome-wide PLINK 2 files and genotype PCs
  genotypes/loci/Lxx.*          each teaching locus: PLINK 2, VCF, ALT dosages (.traw)
  loci.tsv                      the 10 teaching loci and their scenarios
gwas/EUR, gwas/AFR/
  trait_q.sumstats.tsv.gz, disease.sumstats.tsv.gz   per-variant beta, se, z, -log10 p
  summary.tsv
truth/                          the answer key; do not peek before the exercises
  truth.tsv                     one row per causal variant effect, with its path to the trait
  edges.tsv, nodes.tsv, model.json      the full causal graph (used to build the GWAS cohort)
  causal_variants.tsv, de_genes.tsv, teaching_genes.tsv, cpg_parameters.tsv
  latent_small_cohort.tsv.gz    simulated latent expression and trait values
work/                           intermediate files (can be deleted after the build)
```

Coordinates are GRCh38. Variant IDs are `chr:pos:ref:alt`; dosages and effects count the ALT allele.

## Sources and terms

| Data | Source | Terms |
|---|---|---|
| Genotypes | HGDP+1kG phased haplotypes v2, gnomAD (Koenig et al. 2024, Genome Research) | open access; cite the 1000 Genomes Project, HGDP (Bergström et al. 2020) and Koenig et al. 2024 |
| Sample metadata | gnomAD HGDP+1kG metadata (`gnomad_meta_updated.tsv`) | open access |
| Expression | GEUVADIS, ArrayExpress E-GEUV-1 (Lappalainen et al. 2013, Nature) | open access; cite Lappalainen et al. 2013 |
| Genes | GENCODE v39 (mirrored by gnomAD) | open access |
| Simulated layers | this repository | CC BY 4.0 |

Licence details are to be checked against each source when the bundle is published.
