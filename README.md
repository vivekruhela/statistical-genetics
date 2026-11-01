# statistical-genetics

Code and data recipes for the statistical genetics blog series on vivekruhela.github.io:
GWAS, differential expression, eQTL, fine-mapping (SuSiE), colocalization, SMR, TWAS and ROH,
all run on one cohort with genotypes, expression, clinical traits and methylation for the same people.

Status: work in progress. Only steps 1 and 2 of the data build exist so far.

## 1. Set up (once)

Needs conda or mamba (Miniforge recommended). Two environments, because the R and Python stacks
pin different NumPy versions:

```bash
mamba env create -f env-r.yml        # R track + PLINK 2 + BCFtools
mamba env create -f env-python.yml   # Python track + stdpopsim/SLiM simulation
```

## 2. Build the data

```bash
conda activate statgen-r

# Step 1: genotypes (HGDP+1kG phased panel, GRCh38). Test on one chromosome first:
CHROMS="22" bash simulate/01_fetch_genotypes.sh data/eur_yri_samples.txt data/genotypes
# then all autosomes (about 15-20 min on a good connection, a few GB downloaded, ~0.5 GB kept):
bash simulate/01_fetch_genotypes.sh data/eur_yri_samples.txt data/genotypes

# Step 2: GEUVADIS expression (not yet tested)
bash simulate/02_fetch_geuvadis.sh data/geuvadis
```

`data/eur_yri_samples.txt` lists the 654 1000 Genomes people from CEU, FIN, GBR, TSI and YRI,
the populations GEUVADIS sequenced. It is a superset of the GEUVADIS cohort and will be cut to the
people with expression data in step 3.

Tested so far (29 Oct 2026): step 1 on chr22 gives 654 samples and 159,612 SNVs at MAF >= 1%
in about 30 seconds.

## Next steps (to be written)

3. Match genotype and expression samples; normalise expression.
4. Pick teaching loci from the data (strong cis-eQTL genes, some with two independent signals).
5. Simulate clinical traits, methylation and planted causal variants; write `data/truth/`.
6. Simulate the large GWAS cohort (stdpopsim) and its summary statistics.
7. Tests: every step checked against the truth files, and R vs Python results compared.
