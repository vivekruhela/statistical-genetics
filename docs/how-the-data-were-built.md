# How the data were built

This is the reference for the advanced track of post 0. Every number below is a setting in
`config/cohort.json`.

## Real layers

**People.** The GEUVADIS people (CEU, FIN, GBR, TSI from Europe; YRI from Nigeria) who are also in
the HGDP+1kG phased panel. Sex and population come from the panel metadata. Where GEUVADIS has
more than one RNA-seq library for a person, the deepest library is kept.

**Genotypes.** Biallelic SNVs with minor allele frequency of at least 1% in the cohort, from the
gnomAD HGDP+1kG phased panel (GRCh38, phased with SHAPEIT5).

**Expression.** GEUVADIS gene counts, matched to GENCODE v39 by Ensembl gene ID. Genes are kept if
they are protein coding or lncRNA on an autosome and have at least 6 reads and 0.5 counts per
million in at least 20% of people. Counts are divided by DESeq2-style median-of-ratios size
factors, transformed with log2(x + 1), then rank-transformed to a normal distribution per gene
(inverse normal transform). The top 10 principal components of that matrix are the expression PCs.

## Choosing the loci

A cis-eQTL scan runs over every gene: each variant within 750 kb of the gene's transcription
start site with MAF of at least 5% is tested against expression, after removing sex, 5 genotype
PCs and 10 expression PCs from both. The scan keeps the lead variant per gene and the best
variant after conditioning on the lead.

The 10 loci are the strongest eQTL genes (lead -log10 p of at least 15) whose 1.5 Mb window has
2,000 to 15,000 scanned variants, avoiding the MHC and the 8p23 and 17q21.31 inversions, at most
two per chromosome and at least 1 Mb apart. Each locus then takes the next scenario from the list.

## Simulated layers: one causal graph

Everything simulated is a node in one causal graph (`data/truth/edges.tsv` and `nodes.tsv`):

- observed inputs: genotype dosages of the causal variants, and age, sex and BMI on a fixed
  standard scale;
- latent nodes, each with variance 1: a regional factor per locus, CpGs, teaching-gene
  expression, the trait and the disease liability.

Each edge has a target share of variance `r2`. Its weight is `sign * sqrt(r2 / var(parent))`,
where `var(parent)` is `2p(1-p)` for a variant (p is the allele frequency in the whole cohort) and
1 for a latent node. A node is the sum of its weighted, centred parents plus normal noise with
variance `1 - sum(r2)`. The weights are fixed once, from the small cohort, and reused for the
GWAS cohorts, so both cohorts share exactly one model. Variant effects therefore carry over across
cohorts while their variance explained changes with allele frequency, as in real data.

**Causal variants.** At each locus the teaching gene sits 15 kb from the anchor gene. Causal
variants are SNVs with MAF of at least 10% in EUR and 5% in YRI, eQTL variants within 200 kb of
the teaching gene. A second independent variant has r² below 0.1 with the first (in EUR); the
linkage scenario's trait variant has r² between 0.4 and 0.75 with the eQTL variant; the ancestry
scenario's variant has MAF of 3 to 12% in EUR and at least 30% in YRI.

**Effect sizes.** An eQTL variant explains 20% of its gene's expression (10% for a second eQTL);
an mQTL variant explains 35% of its CpG, and that CpG explains 45% of the gene's expression,
always lowering it. For the trait, the target is the share explained by one causal variant: 1.2%,
or 8% for the strong direct locus that even 460 people can detect. When a variant acts through
expression, the expression-to-trait edge is scaled up to hit that target. The disease variant
explains 2% of the liability. Age, sex and BMI explain 5%, 3% and 4% of both the trait and the
liability.

**Clinical.** Age is uniform from 20 to 75; BMI is log-normal around 26 and rises slightly with
age. `trait_q = 50 + 10 * node`. The disease is liability above the 70th percentile of a standard
normal, so about 30% of people are cases.

**Methylation.** Each teaching gene has 8 promoter CpGs within 2 kb of its start, sharing a
regional factor (25% of variance). Another 400 null CpGs are spread over the loci. A random 30% of
CpGs change with age (10% of variance). M-values are `mean + scale * node`, with promoter means
from -4 to -1 (mostly unmethylated) and null means from -3 to 3.

**Teaching-gene expression.** Counts are negative binomial with mean
`2^(7 + 0.6 * node) * size factor` and dispersion 0.05, using each person's real size factor.

**Disease effects on expression.** 100 well-expressed real genes go up and 100 go down in cases
(log2 fold change about 0.7, genes at least 1 Mb from any locus): in cases, an up gene gains
Poisson extra reads and a down gene's reads are binomially thinned. Expression is then
renormalised with the teaching genes included.

## The GWAS cohorts

New people are mosaics of phased haplotypes from panel people who are not in the cohort:
261 unrelated Europeans (1000 Genomes IBS and HGDP Europeans) and 466 unrelated Africans
(1000 Genomes ESN, GWD, LWK, MSL and HGDP African groups, leaving out the admixed ASW and ACB and
the deeply diverged San, Biaka and Mbuti). Each new haplotype copies one template and switches to
another at random breakpoints, on average every 50 kb, miscopying 0.1% of alleles, a simple
copying model in the spirit of HAPGEN. This keeps real local LD at exactly the cohort's variants.
Mosaics from a few hundred templates share more long segments than real people do, which is fine
at single loci but is why the GWAS cohort only covers the teaching loci.

The people go through the same causal graph, with new age, sex and BMI. Every locus variant is
tested against `trait_q` by least squares and against the disease by a logistic score test (beta
is the one-step log odds ratio U / I), both adjusted for age and sex.

## Checks

`tests/test_data.py` confirms that the planted eQTLs, mQTLs and CpG-to-expression effects are
visible in the small cohort with the right sign, that the gwas-only gene has no eQTL, that the
disease prevalence and the DE genes' directions are right, and that the GWAS finds every trait
variant (except the ancestry variant in EUR, which is rare there by design) while the eQTL-only
locus stays quiet.
