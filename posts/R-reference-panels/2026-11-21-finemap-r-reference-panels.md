---
title: "Statistical fine-mapping, post R: reference panels, or borrowing other people's genomes"
date: 2026-11-21
type: method
summary: "What a reference panel is, the five jobs it does, and what happens when you borrow from the wrong one: a matched panel fills in genomes and fine-maps a locus well, a mismatched one quietly invents answers."
tags: ["fine-mapping", "reference panels", "LD", "imputation", "statistical genetics", "tutorial"]
draft: true
---

<!-- Numbers in this draft come from the chromosome 22 test build. Refresh them from
     results/post-R/*.tsv after the full run (look for the figures quoted in figs 2, 3 and 5). -->

[![Code](https://img.shields.io/badge/code-GitHub-181717?logo=github)](https://github.com/vivekruhela/statistical-genetics/tree/main/posts/R-reference-panels)
[![Languages](https://img.shields.io/badge/code-R%20%2B%20Python-2f6f9f)](https://github.com/vivekruhela/statistical-genetics)
[![License](https://img.shields.io/badge/license-MIT-blue)](https://github.com/vivekruhela/statistical-genetics/blob/main/LICENSE)

## Contents

- [Borrowed notes](#borrowed-notes)
- [What a reference panel is](#what-a-reference-panel-is)
- [Five jobs a panel does](#five-jobs-a-panel-does)
- [Why ancestry matters: every population has its own LD map](#why-ancestry-matters-every-population-has-its-own-ld-map)
- [Job one in action: filling in the blanks](#job-one-in-action-filling-in-the-blanks)
- [Job two in action: LD for fine-mapping](#job-two-in-action-ld-for-fine-mapping)
- [How to tell when a panel doesn't fit](#how-to-tell-when-a-panel-doesnt-fit)
- [Which panel should I use?](#which-panel-should-i-use)
- [Try it yourself](#try-it-yourself)
- [Concepts in this post](#concepts-in-this-post)
- [References](#references)

## Borrowed notes

Say you missed three weeks of lectures and the exam is on Monday. You borrow a classmate's notes.
Whether that saves you depends on three things. Did they take the *same course*? Notes from a
different section with a different teacher will look useful and still mislead you. Did they write
*enough*? Two pages of scribbles leave you guessing. And did they write in a *notation you can
read*? If they numbered the chapters differently, you'll put the right facts in the wrong places.

Almost every genetic study borrows notes like this. The notes are called a **reference panel**,
and the same three questions decide whether borrowing helps you or quietly hurts you: is the panel
from people like yours, is it big enough, and does it use the same coordinates and alleles as your
data. This post is about those three questions. It is a short detour before the series really
starts, because every later post leans on a reference panel somewhere.

## What a reference panel is

A reference panel is a group of people whose genomes were read completely and carefully, and then
shared so that anyone can use them. The best known are the **1000 Genomes Project** (about 2,500
people from 26 populations, later sequenced again at high depth), the **Human Genome Diversity
Project** (HGDP, about 900 people from 54 populations), and the combined, cleaned-up release of the
two that gnomAD made, which is what our own playground is built from.

Most studies don't read whole genomes. They use a cheap genotyping chip that measures a few hundred
thousand positions, or they only have *summary statistics*: one line per variant from somebody
else's GWAS, with no people attached. Either way, something is missing, and the panel is where you
borrow it from. You are trusting that the panel's people carry the same stretches of DNA, in the
same combinations, as your own people do.

## Five jobs a panel does

| Job | What you borrow | Where it comes up in this series |
|---|---|---|
| **Ancestry map (PCA)** | a fixed map of world ancestry to place your people on | post 1.3 |
| **Imputation** | the variants your chip didn't measure, filled in from haplotypes in the panel | below, and post 1.2 |
| **LD for summary statistics** | how variants travel together, when you have z-scores but no people | below, and posts 2.2, 4.1 to 5.2 |
| **Local ancestry** | which ancestry each stretch of an admixed genome came from | post 2.5 |
| **QTL lookups** | what a variant does to genes or methylation, from studies like GTEx or eQTLGen | posts 3.3, 5.1 to 5.3 |

The last one is a reference panel in a looser sense: you borrow *results* rather than genomes. But
the same three questions apply, as we'll see.

The two jobs that matter most for fine-mapping are imputation and LD, so those are the two I'll
show you on our own data. First, though, the reason ancestry matters so much.

## Why ancestry matters: every population has its own LD map

In post 0 we met **linkage disequilibrium (LD)**: variants that sit close together are inherited in
blocks, so they are correlated. How long those blocks are depends on a population's history. Every
generation, recombination cuts the blocks a little shorter. African populations are the oldest and
have stayed large, so their genomes have been cut up for longer and their blocks are short. People
who left Africa went through a small founding population, which reset a lot of that shuffling, so
European and Asian genomes carry longer blocks.

Here is that difference in our playground. On the left is one 150 kb stretch near locus L02: each
square is the r² between two common variants, European people below the diagonal and Yoruba people
above it. On the right is the average r² by distance, over all ten teaching loci. I used the same
number of people from each group, because LD estimated from few people is inflated by chance.

<iframe src="/finemap/R/ld-maps.html" title="LD map of one region in European and Yoruba people, and LD decay with distance" width="100%" height="450" style="border:0;" loading="lazy"></iframe>

The European half has big dark blocks, the Yoruba half has smaller, paler ones, and on the right
the Yoruba curve sits below the European one at every distance. Neither is "better". They are two
different maps of the same territory. A panel from one population hands you that population's map,
and the trouble starts when you use it to navigate somewhere else.

## Job one in action: filling in the blanks

Let's pretend our 443 people were typed on a chip. At each of the ten loci I kept only one in every
ten common variants, threw the other nine away, and asked a standard imputation program (Beagle) to
fill them back in. It does this by matching each person's chip variants to stretches of haplotypes
in the panel and copying in what those stretches carry. I tried three panels, all made of people who
are *not* in our cohort: 258 Europeans, 258 Africans, and the two together. Then I compared the
filled-in genotypes with the real ones we threw away.

<iframe src="/finemap/R/imputation.html" title="Imputation accuracy by panel, for European and Yoruba people" width="100%" height="440" style="border:0;" loading="lazy"></iframe>

The y-axis is the imputation r²: 1 means the panel brought the variant back perfectly, 0 means it
had no idea. With a matched panel, common variants come back almost perfectly (r² around 0.9) in
both groups. With the wrong panel, things fall apart, and they fall apart worst for the Yoruba
genomes imputed from Europeans: the average r² drops to between 0.1 and 0.5. Those short African
LD blocks are not in the European notes at all.

Two more things are worth noticing. Rarer variants are always harder, whatever the panel, because
fewer panel haplotypes carry them. And the combined panel does as well as the matched one: a good
imputation program picks the haplotypes that fit each person, so adding more people rarely hurts.
That is why the big modern panels, like TOPMed, mix many ancestries on purpose. Hold on to that
thought, because the next job behaves differently.

## Job two in action: LD for fine-mapping

Remember the formula for fine-mapping from post 0: the evidence for each variant, *together with*
how variants travel together. When you have the people, you compute that LD from them. But most of
the time you only have a GWAS's summary statistics, and you borrow the LD from a panel. The method
then assumes the panel's LD is the LD of the GWAS people. If it isn't, it tries to explain the
mismatch, and the only tool it has is more causal variants.

So let's take our European GWAS (10,000 people) and fine-map the same z-scores with SuSiE (post
4.2 explains how SuSiE works; here it is a black box) and LD from five sources:

- **the GWAS people themselves** (in-sample LD): the ideal, and rarely shared;
- **a European panel**: the 373 European people of our cohort, who are real Europeans but not the
  same Europeans as the GWAS (much like using 1000 Genomes for a UK study);
- **a small European panel**: 50 of them;
- **an African panel**: 466 African-ancestry people;
- **a mixed panel**: the Europeans and the Africans together, like using all of 1000 Genomes.

This is locus L05, where one planted variant has a strong effect on the trait. The top row is the
GWAS. Each row below is one panel: dots are variants, their height is the PIP, the coloured ones
are in a 95% credible set, and the orange ring marks the real cause from the answer key.

<iframe src="/finemap/R/finemap-panels.html" title="Fine-mapping the same GWAS with LD from five panels" width="100%" height="880" style="border:0;" loading="lazy"></iframe>

With in-sample LD, fine-mapping does exactly what we want: one credible set, one variant, and it
is the real cause. The European panel also finds it, but adds a second credible set out of nothing.
The small panel finds the cause and *eight* invented signals. The African panel loses the cause
completely and puts ten confident credible sets in the wrong places. And the mixed panel, which
did so well at imputation, misses the cause and invents nine.

That last one surprised me the first time I saw it. For imputation, more people is more choice.
For LD, a mixed panel is simply a different population, one that matches nobody: its LD is a blend
that no real GWAS has.

One locus is one story, so here is the scorecard across every locus that had a GWAS hit:

<iframe src="/finemap/R/scorecard.html" title="Across all loci: real causes caught and false credible sets, by LD panel" width="100%" height="340" style="border:0;" loading="lazy"></iframe>

The pattern holds. In-sample LD catches 7 of the 8 real causes with no false signals. A large
matched panel comes close. Everything else loses real causes and, worse, hands you confident false
ones. And look at the last row: perfect in-sample LD, with just one variant whose alleles were
swapped in the summary statistics. That single bookkeeping error is enough to break most loci.
This is why fine-mapping papers ask for in-sample LD, and why the stronger your GWAS, the more an
imperfect panel hurts: tiny LD errors don't matter for z = 3, but they matter a lot for z = 20.

## How to tell when a panel doesn't fit

The good news is that a mismatch leaves fingerprints, and you can look for them before you trust
a credible set. The idea is simple. If the LD matrix is right, each variant's z-score should be
predictable from its neighbours' z-scores: a variant in perfect LD with a z = 10 neighbour should
itself have z close to 10. So, for every variant, predict its z-score from all the others using the
panel's LD, and compare it with what the GWAS actually saw.

<iframe src="/finemap/R/z-check.html" title="Observed versus predicted z-scores, for three LD sources" width="100%" height="410" style="border:0;" loading="lazy"></iframe>

With in-sample LD, every point sits on the diagonal. With the African panel, the predictions are
all over the place, and many variants (orange) are more than four standard deviations from their
prediction. With the flipped allele, everything fits except one point that sits far off on the
other side of the diagonal: the flipped variant, with its sign the wrong way round.

Each panel also gets one number, *s*, which measures how much of the z-scores the panel's LD can't
explain: 0 means a perfect fit, and anything above about 0.1 is a warning. The susieR package
calls these tools `estimate_s_rss` and `kriging_rss`, and DENTIST is a stand-alone tool built on
the same idea. Flagged variants are worth a look: often they are allele flips, strand problems or
low-quality variants, and removing them helps. A high *s* across the board means the panel itself
is the wrong one, and no amount of cleaning fixes that.

## Which panel should I use?

Pick what you need the panel for and who is in your study, and the chart gives you a place to
start, with the reason and the usual trap.

<iframe src="/finemap/R/decision-chart.html" title="Which reference panel: pick a task and an ancestry" width="100%" height="580" style="border:0;" loading="lazy"></iframe>

If you remember one thing from this post, make it this: **for imputation, bigger and more diverse
is better; for LD, matched is everything.** In-sample LD beats any panel. A large panel from the
same population is the next best. A small, mixed or mismatched panel will give you confident
answers that are wrong, and the only protection is to check the fit before you believe them.

## Try it yourself

All of this runs on the playground from post 0. The whole post, every table and figure, in R:

```bash
conda activate statgen-r
bash posts/R-reference-panels/run_all.sh data results/post-R
```

and the same in Python with `TRACK=python` in front. It takes 15 to 30 minutes on a laptop-sized
machine, most of it fine-mapping. `bash posts/R-reference-panels/check.sh data` runs both and
checks that R and Python give the same numbers.

The heart of the fine-mapping experiment is only a few lines. In R, for one locus:

```r
library(susieR); library(data.table)
source("lib/statgen.R")

ss <- fread("data/gwas/EUR/trait_q.sumstats.tsv.gz")[locus == "L05" & af >= 0.01 & af <= 0.99]
ss <- ss[abs(pos - pos[which.max(abs(z))]) <= 250000]           # 250 kb either side of the top hit
samples <- fread("data/cohort/samples.tsv")
geno <- read_traw("data/cohort/genotypes/loci/L05.traw")
panel <- geno$g[ss$variant_id, samples[superpopulation == "EUR", sample]]

R <- ld_matrix(panel)                       # LD borrowed from the European panel
keep <- !is.na(diag(R))                     # drop variants that don't vary in the panel
fit <- susie_rss(ss$z[keep], R[keep, keep], n = 10000)
fit$sets$cs                                 # the credible sets
check <- ld_check(ss$z[keep], R[keep, keep], n = 10000)
check$s[1]                                  # how badly the panel fits: 0 is perfect
```

And in Python, with the SuSiE port from the repository:

```python
import sys
import numpy as np, pandas as pd
sys.path.insert(0, "lib")
import statgen as sg
from susie_rss import susie_rss

ss = pd.read_csv("data/gwas/EUR/trait_q.sumstats.tsv.gz", sep="\t")
ss = ss[(ss.locus == "L05") & ss.af.between(0.01, 0.99)]
ss = ss[(ss.pos - ss.pos[ss.z.abs().idxmax()]).abs() <= 250000]
samples = pd.read_csv("data/cohort/samples.tsv", sep="\t")
info, g, people = sg.read_traw("data/cohort/genotypes/loci/L05.traw")
rows = info.reset_index().set_index("id").loc[ss.variant_id, "index"]
cols = [people.index(s) for s in samples["sample"][samples.superpopulation == "EUR"]]

R = sg.ld_matrix(g[np.ix_(rows, cols)])
keep = ~np.isnan(np.diag(R))
fit = susie_rss(ss.z.to_numpy()[keep], R[np.ix_(keep, keep)], n=10000)
print([s["variants"] for s in fit["sets"]])
print(sg.ld_check(ss.z.to_numpy()[keep], R[np.ix_(keep, keep)], n=10000).s[0])
```

A few things to try:

1. Swap the European panel for the Yoruba people (`superpopulation == "AFR"`) and watch the number
   of credible sets and *s* go up.
2. Use `L = 1` in `susie_rss`. With only one causal variant allowed, the LD barely matters (post 4.1
   shows why), so even the wrong panel gives a sensible answer at a single-signal locus. LD becomes
   critical exactly when a locus can have more than one cause.
3. In the imputation step, change `CHIP_EVERY=10` to 30 (a sparser chip) and see which panel copes.

<details>
<summary><strong>For the curious</strong>: builds, alleles, access, and one thing I fixed in our data</summary>

**The notation problem: builds, chromosome names and alleles.** The same variant has different
coordinates in genome build GRCh37 (hg19) and GRCh38 (hg38), so data and panel must be on the same
build; tools such as `bcftools +liftover`, CrossMap or UCSC liftOver convert between them, and lose
a few variants on the way. Even within a build, files disagree about names: our cohort VCFs, written
by PLINK 2, call the chromosome `22`, while the panel files call it `chr22`, and Beagle refused to
run until I made them match. Then alleles: a GWAS reports effects for one allele, the panel counts
another, and if you don't align them the sign of the LD flips, which is exactly the flipped-allele
row above. The worst cases are A/T and C/G variants, where the two strands look the same and only
the allele frequency can tell you which way round they are; many pipelines simply drop them.

**The panels, side by side.** Sizes are for the releases most people use at the time of writing.

| Panel | People | Who | Build | How you get it | Mostly for |
|---|---|---|---|---|---|
| 1000 Genomes (phase 3) | 2,504 | 26 populations, 5 continents | GRCh37 | open download | PCA, small LD panels |
| 1000 Genomes 30x (NYGC) | 3,202, including 602 trios | same 26 populations | GRCh38 | open download | imputation, LD, phasing |
| HGDP | 929 | 54 populations | GRCh38 | open download | diversity, PCA, local ancestry |
| HGDP + 1000 Genomes (gnomAD) | about 4,100 | about 80 populations | GRCh38 | open download | PCA, local ancestry; this series |
| gnomAD v4 | about 807,000 (76,000 genomes) | many ancestries | GRCh38 | allele frequencies only | frequency lookups |
| TOPMed | about 97,000 (r2) to 133,000 (r3) | many ancestries | GRCh38 | through the imputation server only | imputation |
| HRC | 32,470 | mostly European | GRCh37 | imputation servers; EGA with approval | imputation |
| UK Biobank / Pan-UKB LD | LD matrices only | European; six groups in Pan-UKB | GRCh37 | open download | LD for summary statistics |
| GTEx v8 | 838 donors, 49 tissues | about 85% European | GRCh38 | results open; genotypes by application | eQTL lookups |
| eQTLGen (phase 1) | 31,684 | mostly European, blood | GRCh37 | results open | blood eQTL lookups |
| eQTL Catalogue | dozens of studies | mixed, reprocessed uniformly | GRCh38 | results open | eQTL lookups, colocalization |
| GoDMC | about 30,000 | European, blood | GRCh37 | results open | methylation QTL lookups |

**How our panels were built.** The European and African "panels" in this post are real people from
the HGDP + 1000 Genomes release who are not in our cohort: the same people step 9 uses as templates
for the big GWAS. The European GWAS templates are mostly Iberian and HGDP Europeans, while our
cohort's Europeans are British, Finnish, Tuscan and CEPH-Utah, so the "European panel" is matched
the way 1000 Genomes is matched to a typical European study: same continent, not the same people.

**One thing this post made me fix.** The GWAS genomes in our playground are stitched together from
real haplotypes. The first build switched haplotypes every 50 kb on average, and that quietly cut
LD between variants more than about 50 kb apart: at 50 to 200 kb, the GWAS people had an average
|r| of 0.02 against 0.11 in real Europeans. No real panel could match that, and even the European
panel made six to eight false credible sets per locus. Switching every 500 kb instead keeps real
LD across the whole fine-mapping window, so step 9 was rebuilt that way. It is a nice example of
this post's point, from the other side: the LD in your summary statistics and the LD in your panel
have to describe the same people.

**What the numbers here are.** SuSiE-RSS was run with its defaults (up to 10 effects, 95% credible
sets, purity of at least 0.5) and the GWAS sample size. A credible set is "false" if it contains
none of the planted causal variants for that trait. The fit check uses the null-likelihood estimate
of *s* and flags a variant when its z-score is more than 4 standard deviations from the prediction.

</details>

## Concepts in this post

- **Reference panel:** a set of fully sequenced, openly shared genomes that you borrow information
  from: missing variants, LD, or an ancestry map.
- **Imputation:** filling in unmeasured variants by matching a person's measured variants to
  haplotypes in a panel.
- **Imputation r²:** the squared correlation between imputed and true genotypes; the standard
  measure of imputation accuracy.
- **In-sample LD:** LD computed from the same people as the summary statistics.
- **LD mismatch:** a difference between the LD of the panel and the LD of the people behind the
  summary statistics; it makes fine-mapping invent or lose signals.
- **Allele alignment:** making sure the summary statistics and the panel count the same allele of
  every variant.
- **Genome build:** the version of the reference genome that coordinates refer to (GRCh37 or
  GRCh38).

## References

1. The 1000 Genomes Project Consortium. A global reference for human genetic variation.
   *Nature* 526, 68–74 (2015). [doi:10.1038/nature15393](https://doi.org/10.1038/nature15393)
2. Byrska-Bishop M et al. High-coverage whole-genome sequencing of the expanded 1000 Genomes
   Project cohort including 602 trios. *Cell* 185, 3426–3440 (2022).
   [doi:10.1016/j.cell.2022.08.004](https://doi.org/10.1016/j.cell.2022.08.004)
3. Bergström A et al. Insights into human genetic variation and population history from 929
   diverse genomes. *Science* 367, eaay5012 (2020).
   [doi:10.1126/science.aay5012](https://doi.org/10.1126/science.aay5012)
4. Koenig Z et al. A harmonized public resource of deeply sequenced diverse human genomes.
   *Genome Research* 34, 796–809 (2024).
   [doi:10.1101/gr.278378.123](https://doi.org/10.1101/gr.278378.123)
5. Taliun D et al. Sequencing of 53,831 diverse genomes from the NHLBI TOPMed Program.
   *Nature* 590, 290–299 (2021). [doi:10.1038/s41586-021-03205-y](https://doi.org/10.1038/s41586-021-03205-y)
6. McCarthy S et al. A reference panel of 64,976 haplotypes for genotype imputation.
   *Nature Genetics* 48, 1279–1283 (2016). [doi:10.1038/ng.3643](https://doi.org/10.1038/ng.3643)
7. Browning BL, Zhou Y, Browning SR. A one-penny imputed genome from next-generation reference
   panels. *American Journal of Human Genetics* 103, 338–348 (2018).
   [doi:10.1016/j.ajhg.2018.07.015](https://doi.org/10.1016/j.ajhg.2018.07.015)
8. Zou Y, Carbonetto P, Wang G, Stephens M. Fine-mapping from summary data with the "Sum of
   Single Effects" model. *PLoS Genetics* 18, e1010299 (2022).
   [doi:10.1371/journal.pgen.1010299](https://doi.org/10.1371/journal.pgen.1010299)
9. Benner C et al. Prospects of fine-mapping trait-associated genomic regions by using summary
   statistics from genome-wide association studies. *American Journal of Human Genetics* 101,
   539–551 (2017). [doi:10.1016/j.ajhg.2017.08.012](https://doi.org/10.1016/j.ajhg.2017.08.012)
10. Kanai M et al. Meta-analysis fine-mapping is often miscalibrated at single-variant resolution.
    *Cell Genomics* 2, 100210 (2022). [doi:10.1016/j.xgen.2022.100210](https://doi.org/10.1016/j.xgen.2022.100210)
11. Chen W et al. Improved analyses of GWAS summary statistics by reducing data heterogeneity and
    errors. *Nature Communications* 12, 7117 (2021).
    [doi:10.1038/s41467-021-27438-7](https://doi.org/10.1038/s41467-021-27438-7)
12. Weissbrod O et al. Functionally informed fine-mapping and polygenic localization of complex
    trait heritability. *Nature Genetics* 52, 1355–1363 (2020).
    [doi:10.1038/s41588-020-00735-5](https://doi.org/10.1038/s41588-020-00735-5)
13. Karczewski KJ et al. The mutational constraint spectrum quantified from variation in 141,456
    humans. *Nature* 581, 434–443 (2020). [doi:10.1038/s41586-020-2308-7](https://doi.org/10.1038/s41586-020-2308-7)
14. The GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human
    tissues. *Science* 369, 1318–1330 (2020). [doi:10.1126/science.aaz1776](https://doi.org/10.1126/science.aaz1776)
15. Võsa U et al. Large-scale cis- and trans-eQTL analyses identify thousands of genetic loci and
    polygenic scores that regulate blood gene expression. *Nature Genetics* 53, 1300–1310 (2021).
    [doi:10.1038/s41588-021-00913-z](https://doi.org/10.1038/s41588-021-00913-z)
16. Kerimov N et al. A compendium of uniformly processed human gene expression and splicing
    quantitative trait loci. *Nature Genetics* 53, 1290–1299 (2021).
    [doi:10.1038/s41588-021-00924-w](https://doi.org/10.1038/s41588-021-00924-w)
17. Min JL et al. Genomic and phenotypic insights from an atlas of genetic effects on DNA
    methylation. *Nature Genetics* 53, 1311–1321 (2021).
    [doi:10.1038/s41588-021-00923-x](https://doi.org/10.1038/s41588-021-00923-x)
18. Maples BK, Gravel S, Kenny EE, Bustamante CD. RFMix: a discriminative modeling approach for
    rapid and robust local-ancestry inference. *American Journal of Human Genetics* 93, 278–288
    (2013). [doi:10.1016/j.ajhg.2013.06.020](https://doi.org/10.1016/j.ajhg.2013.06.020)
19. Pan-UKB team. Pan-UK Biobank: summary statistics and LD matrices for six ancestry groups.
    [pan.ukbb.broadinstitute.org](https://pan.ukbb.broadinstitute.org)

*Every figure was made for this post with the code in `posts/R-reference-panels/`.*
