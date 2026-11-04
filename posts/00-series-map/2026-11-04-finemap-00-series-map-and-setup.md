---
title: "Statistical fine-mapping, post 0: series map and setup"
date: 2026-11-04
type: method
summary: "What fine-mapping is for, a map of the 20-post series, and how to build the one cohort every post uses: real genotypes and expression with a simulated answer key."
tags: ["fine-mapping", "GWAS", "eQTL", "statistical genetics", "tutorial"]
draft: false
---

[![Code](https://img.shields.io/badge/code-GitHub-181717?logo=github)](https://github.com/vivekruhela/statistical-genetics)
[![Languages](https://img.shields.io/badge/code-R%20%2B%20Python-2f6f9f)](https://github.com/vivekruhela/statistical-genetics)
[![License](https://img.shields.io/badge/license-MIT-blue)](https://github.com/vivekruhela/statistical-genetics/blob/main/LICENSE)

## Contents

- [Start here: what fine-mapping is for](#start-here-what-fine-mapping-is-for)
- [The series map](#the-series-map)
- [Why one cohort with an answer key](#why-one-cohort-with-an-answer-key)
- [Set up](#set-up)
- [Build the data](#build-the-data)
- [Test the build](#test-the-build)
- [How the data were built](#how-the-data-were-built)
- [The ten teaching loci](#the-ten-teaching-loci)
- [Concepts in this post](#concepts-in-this-post)
- [References](#references)

## Start here: what fine-mapping is for

A genome-wide association study (GWAS) tests millions of variants, one at a time, for a link with
a trait. When it works, it does not give you one variant. It gives you a tower of significant
variants, sometimes hundreds, all inside a stretch of DNA that is inherited as a block. Those
variants are correlated with each other through **linkage disequilibrium (LD)**, so each one
carries some of the signal of the variant that actually does something.

The variant with the smallest p-value is the natural first guess, but it is often not the causal
one. Its p-value depends on how well it is measured, on its allele frequency, and on chance. The
causal variant may sit a few places down the list.

**Statistical fine-mapping** asks a sharper question: given the association statistics and the LD
between the variants, how probable is it that each variant is causal? The answer comes as a
**posterior inclusion probability (PIP)** for each variant, and as **credible sets**: the smallest
groups of variants that contain a causal variant with, say, 95% probability. A locus with one
clean signal may give a credible set of two or three variants. A locus with messy LD may give
thirty.

Fine-mapping is the step between "this region matters" and "this gene, through this variant,
in this tissue". The later posts join it to expression data, which is where the gene comes in:
colocalization asks whether a GWAS signal and an expression signal share one causal variant, and
SMR and TWAS ask whether the gene's expression is the route to the trait.

The series builds up to that, one method per post, on one dataset.

## The series map

Each circle is a post. The genotype lane (GWAS, then runs of homozygosity) runs above the shared
trunk and the expression lane runs below it, and the two meet in fine-mapping. Filled circles are
published and open when clicked. Hover over a circle to see the question that post answers.

<iframe src="/finemap/00/series-map.html" title="Series map: 20 posts in six parts" width="100%" height="480" style="border:0;" loading="lazy"></iframe>

The same map as a list, in reading order:

| Post | Part | Question | Planned |
|---|---|---|---|
| 0 | Setup | What is fine-mapping for, and how do I set up the code and data? | 4 Nov 2026 |
| R | Setup | Which reference panel should I use for PCA, imputation, LD or eQTL lookups? | 21 Nov 2026 |
| 1.1 | Data | What files and formats describe one person across genotypes, expression, clinical data and methylation? | 28 Nov 2026 |
| 1.2 | Data | Which variants and people should I remove, and why? | 5 Dec 2026 |
| 1.3 | Data | How does ancestry confound association, and how do I correct for it? | 12 Dec 2026 |
| 2.1 | GWAS | How do I test every variant for a trait and read a Manhattan plot? | 19 Dec 2026 |
| 2.2 | GWAS | Why does one causal variant light up a whole tower of hits? | 26 Dec 2026 |
| 2.3 | GWAS | How do I turn significant variants into independent loci? | 9 Jan 2027 |
| 3.1 | Expression | How do raw counts become a matrix I can test? | 16 Jan 2027 |
| 3.2 | Expression | Which genes differ between cases and controls, and does that make them causal? | 23 Jan 2027 |
| 3.3 | Expression | Which variants change a gene's expression? | 6 Feb 2027 |
| 4.1 | Fine-mapping | With one causal variant, how sure can I be which one it is? | 13 Feb 2027 |
| 4.2 | Fine-mapping | How do I fine-map a locus with several causal variants? | 27 Feb 2027 |
| 4.3 | Fine-mapping | Can methylation and annotations sharpen the credible sets? | 6 Mar 2027 |
| 5.1 | Integration | Do the GWAS and eQTL signals share one causal variant? | 13 Mar 2027 |
| 5.2 | Integration | Does expression of the gene cause the trait, or is it linkage? | 20 Mar 2027 |
| 5.3 | Integration | Can predicted expression and methylation trace the chain from variant to trait? | 27 Mar 2027 |
| 6 | Capstone | Can the whole pipeline find the causal gene at a new locus? | 10 Apr 2027 |
| 2.4 | ROH | What is a run of homozygosity, and how do callers compare? | 17 Apr 2027 |
| 2.5 | ROH | Why does ROH length depend on ancestry, and can ROH burden explain a trait? | 24 Apr 2027 |

The two posts on runs of homozygosity branch off the GWAS lane and come last.

## Why one cohort with an answer key

Every post uses the same few hundred people, so each new method can be read against what the
earlier posts found. And every method can be **scored**, because part of the data is simulated
from a causal model that is written down.

This matters for fine-mapping in particular. On real data nobody knows which variant is causal,
so a credible set can only be judged by whether it looks reasonable. Here the causal variants are
known, so you can see when a credible set misses, when a colocalization is fooled by LD, and when
a method is right for the wrong reason.

The cohort is a hybrid:

| Layer | Where it comes from |
|---|---|
| Genotypes | **Real.** The HGDP + 1000 Genomes phased panel, for the GEUVADIS people of European (CEU, FIN, GBR, TSI) and Yoruba (YRI) ancestry |
| Expression | **Real** GEUVADIS lymphoblastoid RNA-seq counts, plus 10 simulated teaching genes `SIMG01` to `SIMG10` |
| Clinical | Sex is real. Age, BMI, a quantitative trait and a disease are simulated |
| Methylation | Simulated CpGs around each teaching gene, plus null CpGs |
| GWAS summary statistics | Simulated cohorts of 10,000 European and 4,000 African-ancestry people, built from *other* panel people |
| Answer key | Every causal variant, its effect and the path it acts through |

The full build has 443 people with both genotypes and expression. The large GWAS cohorts exist
because 443 people are too few to detect small genetic effects on a trait; real studies have the
same problem, which is why fine-mapping usually starts from summary statistics of a large GWAS
and uses LD from a smaller reference.

## Set up

The code is in R, with a full Python version of every step next to it. You need conda or mamba
(Miniforge is the simplest install). There are two environments, because the R and Python stacks
pin different NumPy versions:

```bash
git clone https://github.com/vivekruhela/statistical-genetics
cd statistical-genetics
mamba env create -f env-r.yml        # R, PLINK 2 and BCFtools
mamba env create -f env-python.yml   # Python
```

The repository layout:

```
config/cohort.json      every setting and effect size, in one place
lib/statgen.{R,py}      shared helpers: normalisation, cis scan, causal graph, association tests
simulate/               the data build, steps 0 to 9, and run_all.sh
tests/                  unit tests, data tests, and an R versus Python comparison
docs/                   how the data were built, in full
posts/                  code for each post, starting with this one
```

## Build the data

One command runs every step. Each step skips work that is already done, so a failed run can be
started again. Try one chromosome first, with a test configuration that loosens the locus rules
so chromosome 22 alone still gives ten loci:

```bash
conda activate statgen-r
STATGEN_CONFIG=tests/config_test.json CHROMS="22" bash simulate/run_all.sh data_chr22
```

Then the real build, on all autosomes:

```bash
bash simulate/run_all.sh data
```

On a server with a good connection this takes about an hour: around 20 minutes to download the
genotypes, 20 to 30 minutes for the genome-wide eQTL scan, and a few minutes for everything else.
To run the Python versions instead, use `TRACK=python bash simulate/run_all.sh data` from the
`statgen-py` environment.

## Test the build

```bash
bash tests/run_tests.sh data
```

This needs both `Rscript` and `python` on the PATH. It runs the unit tests of the shared helpers
in both languages, checks that the planted causal effects can be found again in the built data,
and reruns the expression and eQTL steps in Python to confirm that they give the same numbers as
R. The steps that draw random numbers (the simulation and the GWAS cohorts) give the same design
but different draws in the two languages; the R build is the one the posts use.

## How the data were built

<details open>
<summary><strong>Beginner track</strong>: the idea in five steps</summary>

1. **Take real people.** About 450 people have both whole-genome genotypes and RNA-seq from blood
   cell lines. Their genotypes and gene expression are kept as they are.
2. **Find genes with strong real eQTLs.** A scan over every gene finds the ones whose expression
   depends strongly on a nearby variant. Ten of them, on different chromosomes, become the
   **teaching loci**.
3. **Plant a story at each locus.** Next to each real gene, a simulated teaching gene is added.
   Chosen real variants are made to change it, and through it, or directly, a trait and a disease.
   Each locus has a different story (the table below).
4. **Simulate the rest of the person.** Age, BMI, the trait, the disease and methylation are drawn
   from one causal graph, so every effect has a known size and direction.
5. **Build a big GWAS.** New genomes are stitched together from short stretches of real
   haplotypes from other panel people, so their LD looks like real LD. The same causal graph gives
   them traits, and their association statistics become the GWAS summary statistics.

</details>

<details>
<summary><strong>Advanced track</strong>: the details</summary>

**Genotypes.** Biallelic SNVs with minor allele frequency (MAF) of at least 1% in the cohort,
from the gnomAD HGDP + 1000 Genomes panel, GRCh38, phased with SHAPEIT5. People listed in the
panel metadata without a relatedness call are not in the genotype files and are dropped.

**Expression.** GEUVADIS gene counts matched to GENCODE v39. A gene is kept if it is protein
coding or lncRNA, on an autosome, with at least 6 reads and 0.5 counts per million in at least
20% of people. Counts are divided by median-of-ratios size factors, log2 transformed, and rank
transformed to a normal distribution gene by gene. The first 10 principal components of that
matrix are the expression PCs.

**The eQTL scan.** Every variant within 750 kb of a gene's start, with MAF of at least 5%, is
tested against that gene's expression after removing sex, 5 genotype PCs and 10 expression PCs.
The scan keeps the lead variant per gene and the best variant after conditioning on the lead.
Loci are the strongest genes (lead $-\log_{10} p \ge 15$) with 2,000 to 15,000 variants in their
1.5 Mb window, avoiding the MHC and two large inversions, at most two per chromosome.

**The causal graph.** Each edge from a parent to a child has a target share of variance $r^2$.
Its weight is

$$
w = \pm\sqrt{\frac{r^2}{\operatorname{var}(\text{parent})}},
$$

where the variance of a variant is $2p(1-p)$ for allele frequency $p$, and 1 for a simulated node.
A node is the sum of its weighted, centred parents plus normal noise with variance
$1 - \sum r^2$, so it has variance 1. The weights are fixed once, in the small cohort, and reused
in the GWAS cohorts. A variant's effect is then the same in both, while its variance explained
changes with its allele frequency, as it does in real data.

**Effect sizes.** An eQTL variant explains 20% of its gene's expression (10% for a second one).
A trait variant explains 1.2% of the trait, or 8% at the one strong locus that 443 people can
detect. A disease variant explains 2% of the disease liability, and the disease is liability
above the 70th percentile, so about 30% of people are cases. Age, sex and BMI explain a few
percent each.

**Methylation.** Each teaching gene has 8 promoter CpGs sharing a regional factor, and 400 null
CpGs are spread over the loci. A third of all CpGs drift with age. At the mediated loci, a
variant changes a CpG and the CpG lowers expression.

**Disease effects on expression.** The disease also shifts 200 real genes (100 up, 100 down),
all far from the loci. These are consequences of disease, not causes, which is the trap the
differential expression post is about.

**The GWAS cohorts.** Each new genome is a mosaic: for each 50 kb stretch, both haplotypes are
copied from randomly chosen panel people of the same ancestry who are not in the small cohort,
with a small miscopy rate. LD inside a stretch is real LD. The full build has 10,000 European
people (2,923 cases) and 4,000 African-ancestry people (1,340 cases).

Every number here is a setting in `config/cohort.json`, and the full description is in
`docs/how-the-data-were-built.md` in the repository.

</details>

## The ten teaching loci

Each locus is a 1.5 Mb window around a real gene with a strong eQTL, and each carries one
scenario. Which variants are causal is not given here; finding them is the work of the later
posts. The capstone locus is held back for post 6.

| Locus | Window around | Scenario | What is planted | Main posts |
|---|---|---|---|---|
| L01 | *ERAP2* (chr5) | mediated | a variant changes a CpG, the CpG lowers expression, expression changes the trait | 4.3, 5.1, 5.3 |
| L02 | *SLFN5* (chr17) | two signals | two independent eQTL variants; the trait follows expression | 4.2 |
| L03 | *SNHG5* (chr6) | linkage | an eQTL at one variant, a trait effect at a different variant in moderate LD with it | 5.1, 5.2 |
| L04 | *FLVCR1-DT* (chr1) | pleiotropy | one variant acts on expression and on the trait through separate paths | 5.2 |
| L05 | *ZNF266* (chr19) | GWAS only | a strong direct trait variant, visible even in the small cohort | 2.1, 5.1 |
| L06 | *POMZP3* (chr7) | eQTL only | an eQTL with no effect on any trait | 5.1 |
| L07 | *CENPU* (chr4) | two GWAS signals | two trait variants; only one is also an eQTL | 4.2, 5.1 |
| L08 | *ARL16* (chr17) | binary | an eQTL whose gene raises disease risk | 2.1, 5.1 |
| L09 | a lncRNA (chr5) | ancestry | an eQTL variant common in YRI and rare in Europeans | 4.3 |
| L10 | held back | capstone | several signals at once | 6 |

The scenarios are the ways real loci mislead. A GWAS hit next to an eQTL looks like a causal gene
in L01, L03 and L04, but only L01 is the simple story. L06 is an eQTL that a careless analysis
would link to the trait. L09 is the case where a credible set from European data is wide and one
from African-ancestry data is narrow, because LD blocks are shorter in African populations.

## Concepts in this post

- **Linkage disequilibrium (LD):** the correlation between alleles at nearby variants, because
  they are inherited together. It is why one causal variant makes many variants significant.
- **Fine-mapping:** estimating, for each variant in a locus, the probability that it is causal,
  given the association statistics and the LD.
- **Posterior inclusion probability (PIP):** that probability, for one variant.
- **Credible set:** the smallest set of variants that contains a causal variant with a chosen
  probability, often 95%.
- **eQTL:** a variant associated with the expression level of a gene.
- **Summary statistics:** per-variant effect estimates and standard errors from a GWAS, shared in
  place of individual data. Most fine-mapping starts from them.
- **Reference panel:** a set of sequenced people used to estimate LD, ancestry or missing
  genotypes for another study. Post R is about choosing one.
- **Haplotype mosaic:** a simulated genome copied in stretches from real haplotypes, so that its
  local LD is realistic.

## References

1. Schaid DJ, Chen W, Larson NB. From genome-wide associations to candidate causal variants by
   statistical fine-mapping. *Nature Reviews Genetics* 19, 491–504 (2018).
   [doi:10.1038/s41576-018-0016-z](https://doi.org/10.1038/s41576-018-0016-z)
2. Wang G, Sarkar A, Carbonetto P, Stephens M. A simple new approach to variable selection in
   regression, with application to genetic fine mapping. *Journal of the Royal Statistical
   Society B* 82, 1273–1300 (2020).
   [doi:10.1111/rssb.12388](https://doi.org/10.1111/rssb.12388)
3. Lappalainen T et al. Transcriptome and genome sequencing uncovers functional variation in
   humans. *Nature* 501, 506–511 (2013).
   [doi:10.1038/nature12531](https://doi.org/10.1038/nature12531)
4. Koenig Z et al. A harmonized public resource of deeply sequenced diverse human genomes.
   *Genome Research* 34, 796–809 (2024).
   [doi:10.1101/gr.278378.123](https://doi.org/10.1101/gr.278378.123)
5. The 1000 Genomes Project Consortium. A global reference for human genetic variation.
   *Nature* 526, 68–74 (2015).
   [doi:10.1038/nature15393](https://doi.org/10.1038/nature15393)
6. Giambartolomei C et al. Bayesian test for colocalisation between pairs of genetic association
   studies using summary statistics. *PLoS Genetics* 10, e1004383 (2014).
   [doi:10.1371/journal.pgen.1004383](https://doi.org/10.1371/journal.pgen.1004383)
7. Zhu Z et al. Integration of summary data from GWAS and eQTL studies predicts complex trait
   gene targets. *Nature Genetics* 48, 481–487 (2016).
   [doi:10.1038/ng.3538](https://doi.org/10.1038/ng.3538)
8. Gusev A et al. Integrative approaches for large-scale transcriptome-wide association studies.
   *Nature Genetics* 48, 245–252 (2016).
   [doi:10.1038/ng.3506](https://doi.org/10.1038/ng.3506)

*Figures in this post were made for it with the code in `posts/00-series-map/`.*
