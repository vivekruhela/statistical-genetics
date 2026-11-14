#!/usr/bin/env python
"""Post R, analysis 3a: the fine-mapping windows, and the in-sample LD of the GWAS people.

Python twin of 04_gwas_ld.R. A window is 250 kb either side of the strongest GWAS variant at a
teaching locus, for each trait (trait_q, disease) whose strongest variant passes p < 5e-8 in the
European GWAS; the capstone locus (L10) is left out on purpose. For each window this rebuilds the
European GWAS genomes exactly as step 9 made them and writes their LD matrix (the "in-sample"
LD). The windows are the same in R and Python; the in-sample LD is only the same for the same
build, because the two languages draw different random genomes. This script reads a Python
build (TRACK=python); the R build (the released data) needs 04_gwas_ld.R.
  OUT/windows.tsv
  OUT/ld/insample/<locus>_<trait>.ld.{bin,ids}
Usage: python posts/R-reference-panels/04_gwas_ld.py data results/post-R
"""
import io
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "lib")
import statgen as sg  # noqa: E402

D, OUT = sys.argv[1:3]
HALF_WINDOW, MIN_MLOG10P, MIN_AF = 250000, -np.log10(5e-8), 0.01
cfg = sg.read_config()
gcc = cfg["gwas_cohort"]
loci = pd.read_csv(f"{D}/cohort/loci.tsv", sep="\t")

windows = []
for trait in ["trait_q", "disease"]:
    ss = pd.read_csv(f"{D}/gwas/EUR/{trait}.sumstats.tsv.gz", sep="\t")
    for lid in loci.locus[loci.scenario != "capstone"]:
        x = ss[(ss.locus == lid) & (ss.af >= MIN_AF) & (ss.af <= 1 - MIN_AF)]
        lead = x.loc[x.z.abs().idxmax()]
        if lead.mlog10p < MIN_MLOG10P:
            continue
        w = x[(x.pos - lead.pos).abs() <= HALF_WINDOW]
        windows.append(dict(locus=lid, trait=trait, lead_id=lead.variant_id, lead_pos=lead.pos,
                            lead_mlog10p=lead.mlog10p, start=w.pos.min(), end=w.pos.max(),
                            n_variants=len(w), n_gwas=lead.n))
windows = pd.DataFrame(windows).sort_values(["locus", "trait"]).reset_index(drop=True)
sg.write_tsv(windows, f"{OUT}/windows.tsv")
print(windows.to_string(index=False))

built = pd.read_csv(f"{D}/gwas/EUR/summary.tsv", sep="\t").iloc[0]
if built.built_with != "Python":
    sys.exit(f"these GWAS genomes were built with {built.built_with}; for the in-sample LD run 04_gwas_ld.R instead")
if built.n != gcc["n"]["EUR"]:
    sys.exit(f"the GWAS has {built.n} people but the config says {gcc['n']['EUR']}: "
             "set STATGEN_CONFIG to the config the data were built with")


def build_locus(pop_i, i, n):
    """The GWAS genomes at one locus, exactly as in simulate/09_gwas_cohort.py."""
    rng = np.random.default_rng([cfg["seed"], pop_i, i])
    with open(f"{D}/cohort/genotypes/loci/{loci.locus[i]}.pvar") as fh:
        body = "".join(line for line in fh if not line.startswith("##"))
    pv = pd.read_csv(io.StringIO(body), sep="\t").iloc[:, :5]
    pv.columns = ["chr", "pos", "id", "ref", "alt"]
    H = sg.read_haplotypes(f"{D}/work/gwas_templates/{list(gcc['n'])[pop_i]}/{loci.locus[i]}.hap.gz", pv.id)
    return pv, sg.mosaic_dosage(H, pv.pos.to_numpy(), n, gcc["segment_kb"] * 1000, gcc["miscopy_rate"], rng)


pop_i = list(gcc["n"]).index("EUR")
for w in windows.itertuples():
    prefix = f"{OUT}/ld/insample/{w.locus}_{w.trait}"
    if os.path.exists(prefix + ".ld.bin"):
        print(w.locus, w.trait, "already done")
        continue
    pv, G = build_locus(pop_i, int(np.flatnonzero(loci.locus == w.locus)[0]), gcc["n"]["EUR"])
    ss = pd.read_csv(f"{D}/gwas/EUR/{w.trait}.sumstats.tsv.gz", sep="\t")
    ss = ss[(ss.locus == w.locus) & (ss.af >= MIN_AF) & (ss.af <= 1 - MIN_AF) & (ss.pos >= w.start) & (ss.pos <= w.end)]
    row = {v: k for k, v in enumerate(pv.id)}
    g = G[[row[v] for v in ss.variant_id]].astype(float)
    if np.max(np.abs(g.mean(1) / 2 - ss.af.to_numpy())) > 1e-9:
        sys.exit("rebuilt genomes do not match the summary statistics: was the data built with another config or seed?")
    sg.write_ld(sg.ld_matrix(g), list(ss.variant_id), prefix)
    print(f"{w.locus} {w.trait}: in-sample LD for {len(ss)} variants in {g.shape[1]} people")
