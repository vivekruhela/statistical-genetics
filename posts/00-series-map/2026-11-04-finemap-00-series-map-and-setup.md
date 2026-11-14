---
title: "Statistical fine-mapping, post 0: the series map and our playground"
date: 2026-11-04
type: method
summary: "Why the top GWAS hit is often not the cause, what fine-mapping does about it, and the one cohort with a hidden answer key that we'll use for the whole series."
tags: ["fine-mapping", "GWAS", "eQTL", "statistical genetics", "tutorial"]
draft: false
---

[![Code](https://img.shields.io/badge/code-GitHub-181717?logo=github)](https://github.com/vivekruhela/statistical-genetics)
[![Languages](https://img.shields.io/badge/code-R%20%2B%20Python-2f6f9f)](https://github.com/vivekruhela/statistical-genetics)
[![License](https://img.shields.io/badge/license-MIT-blue)](https://github.com/vivekruhela/statistical-genetics/blob/main/LICENSE)

## Contents

- [A window breaks on a crowded street](#a-window-breaks-on-a-crowded-street)
- [The same thing happens in our DNA](#the-same-thing-happens-in-our-dna)
- [What fine-mapping does about it](#what-fine-mapping-does-about-it)
- [Where we are going](#where-we-are-going)
- [Our playground: real people, hidden answers](#our-playground-real-people-hidden-answers)
- [Try it yourself](#try-it-yourself)
- [Concepts in this post](#concepts-in-this-post)
- [References](#references)

## A window breaks on a crowded street

Picture a busy street. A window breaks, and a camera on the corner shows you everyone who was
standing near it at that moment. Most of them did nothing. But there is a catch: these people
always walk around in the same little groups. Wherever the culprit goes, the same five friends are
right next to them. So every time a window breaks, the camera shows the same faces.

If you only count how often each face shows up near broken windows, the culprit and their friends
look equally guilty. One of the friends might even score a little higher, just by luck of where
they were standing. Catching the right person takes a different question. Not "who was near?",
but "given who always walks with whom, who is most likely to have thrown the stone?"

That second question is what this whole series is about.

## The same thing happens in our DNA

In our genomes, the people on the street are genetic variants, and the groups that walk together
are stretches of DNA that we inherit as blocks. Variants inside a block travel together from
parent to child, so knowing one of them tells you a lot about its neighbours. Geneticists call
this **linkage disequilibrium**, or **LD**.

A genome-wide association study (GWAS) is the camera. It tests each variant on its own and asks
whether people who carry it differ in some trait. When one variant really changes the trait, all
its travelling companions light up too. You don't get one hit. You get a tower of them.

Here is a small example I simulated, so we know the truth. There are 120 variants, 1,500 people,
and exactly one variant that changes the trait. I circled it in orange.

<iframe src="/finemap/00/ld-tower.html" title="One causal variant lights up a tower of GWAS hits; fine-mapping narrows it to three" width="100%" height="570" style="border:0;" loading="lazy"></iframe>

Look at the top panel first. Nineteen variants pass the usual genome-wide significance line, and
the darker the dot, the more tightly it travels with the real cause. The variant with the
smallest p-value, the "top hit", is not the one I planted. It sits right next door and simply
got a little luckier in this sample. If you stopped here and wrote "variant 63 causes this trait"
in a paper, you would be wrong.

## What fine-mapping does about it

Fine-mapping asks the second question from the street: knowing the association signal *and* who
travels with whom, how likely is each variant to be the real cause?

The answer is one number per variant, the **posterior inclusion probability (PIP)**, shown in the
bottom panel. Nineteen suspects shrink to a short list. The blue bars are the **95% credible
set**: the smallest group of variants that, together, contain the real cause with 95%
probability. Here it has three variants, and the one I planted is among them.

Notice what fine-mapping does *not* do. It doesn't point at a single winner, because the data
can't honestly tell these three apart. It tells you how uncertain you should be. A clean locus
gives a credible set of one or two variants; a locus with a big LD block can give thirty. Both
are useful answers.

For a single causal variant, the PIP of variant $j$ is just its share of the evidence:

$$
\text{PIP}_j = \frac{\text{BF}_j}{\sum_k \text{BF}_k},
$$

where $\text{BF}_j$, the Bayes factor, says how much better the data are explained if variant
$j$ is the cause than if nothing is. Post 4.1 builds this up slowly, and post 4.2 handles loci
with more than one cause.

Fine-mapping is also where genes come in. If the short list of variants sits inside a stretch of
DNA that turns a gene up or down, we can start asking whether that gene is the route from DNA to
disease. That is what the expression half of the series is for.

## Where we are going

Each circle below is one post. The top lane is about genotypes (GWAS, and later runs of
homozygosity), the bottom lane is about gene expression, and the two meet in fine-mapping.
Filled circles are published and open when clicked. Hover over or tap a circle to see the
question that post answers. On a phone, swipe the map sideways to see all of it.


<iframe src="/finemap/00/series-map.html" title="Series map: 20 posts in six parts" width="100%" height="490" style="border:0;" loading="lazy"></iframe>

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


## Our playground: real people, hidden answers

On real data, nobody knows which variant is the cause. So you can never be sure whether a
credible set is right, or just looks reasonable. I wanted a dataset where we *can* check.

So I built one, like a mystery novel where the author hid the answers at the back of the book.
The people are real: 443 people from Europe and from the Yoruba in Nigeria, with their real
genomes (from the HGDP and 1000 Genomes projects) and the real activity of their genes in blood
cells (from the GEUVADIS project). On top of that, I planted clues. At ten places in the genome I
chose real variants and made them change a trait, a disease, or the activity of a nearby gene,
each in a different way. Everything I planted is written down in an answer key, so every method
in this series can be graded.

Ten places means ten small mysteries, each modelled on a way real loci fool people:

| Mystery | What I planted | Where we solve it |
|---|---|---|
| L01 | A variant changes the chemical marks on DNA, which turns a gene down, which changes the trait | posts 4.3, 5.1, 5.3 |
| L02 | Two separate variants both change a gene; the trait follows the gene | post 4.2 |
| L03 | One variant changes a gene, a *different* nearby variant changes the trait | posts 5.1, 5.2 |
| L04 | One variant changes both, through two separate routes | post 5.2 |
| L05 | A strong effect on the trait, no gene involved | posts 2.1, 5.1 |
| L06 | A variant changes a gene that has nothing to do with the trait | post 5.1 |
| L07 | Two variants change the trait, only one of them through a gene | posts 4.2, 5.1 |
| L08 | A variant changes a gene that raises disease risk | posts 2.1, 5.1 |
| L09 | A variant that is common in Nigeria and rare in Europe | post 4.3 |
| L10 | Held back for the capstone, post 6 | post 6 |

Each mystery sits around a real gene (L01 at *ERAP2*, L02 at *SLFN5*, L03 at *SNHG5*, L04 at
*FLVCR1-DT*, L05 at *ZNF266*, L06 at *POMZP3*, L07 at *CENPU*, L08 at *ARL16*, L09 at a lncRNA on
chromosome 5), but the stories are made up, so don't read biology into them.

Two more pieces make the playground feel like real research. First, 443 people are far too few
to find small genetic effects, which is also true in real life. So, as in real studies, the big
GWAS comes from *other* people: 10,000 European and 4,000 African-ancestry genomes that I stitched
together from short real stretches of DNA, so their LD looks real. Second, the disease also
changes the activity of 200 genes that have nothing to do with its cause. Those are traps for the
differential expression post, and I'll let you fall into them first.

## Try it yourself

Everything is in R, with a Python version of every step next to it. You can reproduce the figure
above in a few seconds, with no data download:

```bash
git clone https://github.com/vivekruhela/statistical-genetics
cd statistical-genetics
python posts/00-series-map/ld_tower.py ld-tower.html
```

Change `SEED` or `RHO` at the top of the script and watch the tower and the credible set change.
With weaker LD (try `RHO = 0.9`) the tower gets thinner and the credible set shrinks. That one
knob is most of the story of this series.

To build the whole playground, you need conda or mamba, a machine with a decent connection, and
about an hour:

```bash
mamba env create -f env-r.yml          # R, PLINK 2 and BCFtools
mamba env create -f env-python.yml     # Python
conda activate statgen-r
STATGEN_CONFIG=tests/config_test.json CHROMS="22" bash simulate/run_all.sh data_chr22   # a quick trial
bash simulate/run_all.sh data          # the real build, all chromosomes
bash tests/run_tests.sh data           # checks the answer key, and that R and Python agree
```

Each step skips work it has already done, so a failed run can simply be started again. The
`data/truth/` folder is the answer key. Try not to peek until the posts ask you to.

<details>
<summary><strong>For the curious</strong>: exactly how the data were built</summary>

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

**The GWAS cohorts.** Each new genome is a mosaic: for each 500 kb stretch (on average), both haplotypes are
copied from randomly chosen panel people of the same ancestry who are not in the small cohort,
with a small miscopy rate. LD inside a stretch is real LD. The full build has 10,000 European
people (2,923 cases) and 4,000 African-ancestry people (1,340 cases).

Every number here is a setting in `config/cohort.json`, and the full description is in
`docs/how-the-data-were-built.md` in the repository.

</details>

## Concepts in this post

- **Linkage disequilibrium (LD):** nearby variants are inherited together, so they are
  correlated. It is why one cause makes a whole tower of hits.
- **GWAS:** testing every variant, one at a time, for a link with a trait.
- **Fine-mapping:** using the association signal and the LD together to estimate how likely each
  variant is to be the cause.
- **Posterior inclusion probability (PIP):** that likelihood, for one variant.
- **Credible set:** the smallest group of variants that contains the cause with a chosen
  probability, often 95%.
- **eQTL:** a variant that changes how active a gene is.
- **Summary statistics:** per-variant results of a GWAS, shared instead of people's data. Most
  fine-mapping starts from them.

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

*Both figures were made for this post with the code in `posts/00-series-map/`.*
