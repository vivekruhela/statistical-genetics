# statistical-genetics

Code and data recipes for the statistical genetics blog series on vivekruhela.github.io:
GWAS, differential expression, eQTL, fine-mapping (SuSiE), colocalization, SMR, TWAS and ROH,
all run on one cohort with genotypes, expression, clinical traits and methylation for the same people.

Every step has an R version (the main track) and a Python version next to it. Where a step is
deterministic, both produce the same numbers, and the tests check that they do.

Status: the data build (steps 0 to 9) is built genome-wide and passes all tests, in R and in Python. The blog posts'
own code will be added post by post.

## What the data are

A hybrid cohort: real genotypes and real expression, with simulated layers on top so every
method can be scored against a known answer.

| Layer | Source | People |
|---|---|---|
| Genotypes | HGDP+1kG phased panel (gnomAD, GRCh38, SHAPEIT5) | about 462 GEUVADIS people (EUR and YRI) |
| Expression | GEUVADIS LCL RNA-seq (E-GEUV-1), plus 10 simulated teaching genes `SIMG01`-`SIMG10` | same |
| Methylation | simulated CpG M-values around each teaching gene, plus null CpGs | same |
| Clinical | sex is real; age, BMI, a quantitative trait `trait_q` and a disease are simulated | same |
| GWAS summary statistics | simulated cohorts of 10,000 EUR and 4,000 AFR people built from *other* panel people | different people, as in real studies |
| Answer key | `data/truth/`: every causal variant, its effect and the path it acts through | |

The build picks 10 teaching loci (1.5 Mb windows around genes with strong real eQTLs) and gives
each one scenario, so that every method in the series has a case to find:

| Locus | Scenario | What is planted | Main posts |
|---|---|---|---|
| L01 | mediated | variant A changes a CpG, the CpG lowers expression, expression changes the trait | 4.3, 5.1, 5.3 |
| L02 | two_signals | two independent eQTL variants; the trait follows expression | 4.2 |
| L03 | linkage | eQTL at A, trait at a different variant B in moderate LD | 5.1, 5.2 |
| L04 | pleiotropy | A acts on expression and on the trait through separate paths | 5.2 |
| L05 | gwas_only | a strong direct trait variant, visible even in the small cohort | 2.1, 5.1 |
| L06 | eqtl_only | eQTL only, no trait effect | 5.1 |
| L07 | gwas_two | two trait variants; only one is also an eQTL | 4.2, 5.1 |
| L08 | binary | eQTL variant whose expression raises disease risk | 2.1, 5.1 |
| L09 | ancestry | eQTL variant common in YRI, rare in EUR | 4.3 |
| L10 | capstone | CpG-mediated eQTL, a second eQTL and a direct trait variant, held out for post 6 | 6 |

The disease also shifts the expression of 200 real genes (100 up, 100 down; listed in
`data/truth/de_genes.tsv`), so the differential expression post has true positives that are
consequences of disease, not causes. How the simulation works is in
[docs/how-the-data-were-built.md](docs/how-the-data-were-built.md).

## 1. Set up (once)

Needs conda or mamba (Miniforge recommended). Two environments, because the R and Python stacks
pin different NumPy versions:

```bash
mamba env create -f env-r.yml        # R track + PLINK 2 + BCFtools
mamba env create -f env-python.yml   # Python track
```

## 2. Build the data

One command runs every step; each step skips work already done, so a failed run can be restarted.
Try one chromosome first:

```bash
conda activate statgen-r
STATGEN_CONFIG=tests/config_test.json CHROMS="22" bash simulate/run_all.sh data_chr22
```

`tests/config_test.json` loosens the locus rules so chromosome 22 alone still yields 10 loci.
Then the real build, all autosomes, with `config/cohort.json`:

```bash
bash simulate/run_all.sh data
```

Expected time on a server with a good connection: about 20 minutes to download genotypes,
20 to 30 minutes for the genome-wide cis-eQTL scan in step 5, a few minutes for the rest.
To run the Python versions of steps 3, 5, 7 and 9 instead: `TRACK=python bash simulate/run_all.sh data`
(from the `statgen-py` environment, with `plink2` and `bcftools` on the PATH).

The steps, if you want to run them one by one:

| Step | Script | Output |
|---|---|---|
| 0 | `00_fetch_metadata.sh` | sample sex and population, GWAS template people, GENCODE v39 genes |
| 1 | `01_fetch_genotypes.sh` | `data/genotypes/chrN.bcf`, biallelic SNVs with MAF >= 1% |
| 2 | `02_fetch_geuvadis.sh` | GEUVADIS sample sheet and gene counts (from EBI BioStudies, md5 checked) |
| 3 | `03_prepare_expression.{R,py}` | one library per person, gene filter, size factors, log2, inverse normal, expression PCs |
| 4 | `04_genotypes_to_plink.sh` | PLINK 2 files, genotype PCs, dosage tables for the scan |
| 5 | `05_pick_loci.{R,py}` | genome-wide cis-eQTL scan (lead and conditional second signal per gene), `loci.tsv` |
| 6 | `06_extract_loci.sh` | released genotypes: genome-wide PLINK 2 and per-locus PLINK 2, VCF and dosages |
| 7 | `07_simulate_layers.{R,py}` | causal variants, clinical, methylation, teaching genes, DE spike, renormalised expression, truth |
| 8 | `08_fetch_gwas_templates.sh` | phased haplotypes of the GWAS template people at each locus |
| 9 | `09_gwas_cohort.{R,py}` | GWAS cohorts and summary statistics for `trait_q` and `disease` |

Steps 3 and 5 give identical numbers in R and Python. Steps 7 and 9 draw random numbers, so the
two languages build the same design with different draws; the released dataset is the R build.

## 3. Test

```bash
bash tests/run_tests.sh data     # needs Rscript and python on the PATH
```

This runs unit tests of `lib/` in both languages, checks that the planted truth can be found in
the built data (`tests/test_data.py`), and reruns steps 3 and 5 in Python to confirm they match
the R outputs. Without real GEUVADIS access, `tests/make_mock_geuvadis.py` writes a mock count
table in the same layout, which is how this code was tested.

## Layout

```
config/cohort.json      every setting and effect size, in one place
lib/statgen.{R,py}      shared helpers (normalisation, cis scan, causal graph, association tests)
simulate/               the data build, steps 0 to 9, and run_all.sh
tests/                  unit tests, data tests, R vs Python comparison, mock GEUVADIS
docs/                   how the data were built
data/                   built data (not in git; see data/README.md)
```

## Licences

Code: MIT (see LICENSE). Data keep their own terms, listed in [data/README.md](data/README.md).
