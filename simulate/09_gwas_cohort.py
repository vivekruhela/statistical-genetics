#!/usr/bin/env python
"""Step 9 (Python): a large GWAS cohort at the teaching loci, and its summary statistics.

Python twin of simulate/09_gwas_cohort.R: same model and outputs, NumPy random numbers.
Usage: python simulate/09_gwas_cohort.py data/cohort data/truth data/work/gwas_templates data/gwas
"""
import io
import json
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "lib")
import statgen as sg  # noqa: E402

cohort_dir, truth_dir, tmpl_dir, out = sys.argv[1:5]
cfg = sg.read_config()
gcc = cfg["gwas_cohort"]
loci = pd.read_csv(f"{cohort_dir}/loci.tsv", sep="\t")
edges = pd.read_csv(f"{truth_dir}/edges.tsv", sep="\t")
nodes = pd.read_csv(f"{truth_dir}/nodes.tsv", sep="\t")
with open(f"{truth_dir}/model.json") as fh:
    model = json.load(fh)
causal = set(edges["from"][edges["from"].str.startswith("var:")].str.replace("^var:", "", regex=True))


def locus_variants(lid):
    with open(f"{cohort_dir}/genotypes/loci/{lid}.pvar") as fh:
        body = "".join(line for line in fh if not line.startswith("##"))
    pv = pd.read_csv(io.StringIO(body), sep="\t").iloc[:, :5]
    pv.columns = ["chr", "pos", "id", "ref", "alt"]
    return pv.reset_index(drop=True)


def build_locus(pop_i, i, n):
    rng = np.random.default_rng([cfg["seed"], pop_i, i])
    pv = locus_variants(loci.locus[i])
    H = sg.read_haplotypes(f"{tmpl_dir}/{list(gcc['n'])[pop_i]}/{loci.locus[i]}.hap.gz", pv.id)
    return pv, sg.mosaic_dosage(H, pv.pos.to_numpy(), n, gcc["segment_kb"] * 1000, gcc["miscopy_rate"], rng)


for pop_i, (pop, n) in enumerate(gcc["n"].items()):
    print(f"{pop}: building {n} people")
    inputs = {}
    for i in range(len(loci)):
        pv, G = build_locus(pop_i, i, n)
        for j in np.flatnonzero(pv.id.isin(causal)):
            inputs["var:" + pv.id[j]] = G[j].astype(float)
    rng = np.random.default_rng([cfg["seed"], pop_i])
    covs = sg.simulate_covariates(n, rng)
    inputs.update(sg.standard_covariates(covs.age.to_numpy(), covs.male.to_numpy(), covs.bmi.to_numpy()))
    vals = sg.simulate_dag(nodes, edges, inputs, n, rng)
    trait_q = model["trait_q"]["offset"] + model["trait_q"]["scale"] * vals["trait_q"]
    disease = (vals["disease_liability"] > model["disease_threshold"]).astype(float)
    X = np.column_stack([np.ones(n), covs.age, covs.male])

    res = []
    for i in range(len(loci)):
        pv, Gall = build_locus(pop_i, i, n)
        af = Gall.mean(1) / 2
        for s0 in range(0, len(pv), 1000):
            s = slice(s0, min(len(pv), s0 + 1000))
            for tr, y in (("trait_q", trait_q), ("disease", disease)):
                st = sg.assoc_stats(Gall[s], y, X, binary=tr == "disease")
                block = pv.iloc[s].rename(columns={"id": "variant_id", "alt": "effect_allele"})
                block = block[["variant_id", "chr", "pos", "ref", "effect_allele"]].reset_index(drop=True)
                block.insert(0, "locus", loci.locus[i])
                block.insert(0, "trait", tr)
                block["af"], block["n"] = af[s], n
                res.append(pd.concat([block, st], axis=1))
        print(f"  {pop} {loci.locus[i]} done")
    res = pd.concat(res, ignore_index=True)
    res = res[(res.af > 0) & (res.af < 1)]
    for tr in ("trait_q", "disease"):
        x = res[res.trait == tr].drop(columns="trait")
        if tr == "disease":
            x = x.assign(n_cases=int(disease.sum()), n_controls=int(n - disease.sum()))
        sg.write_tsv(x, f"{out}/{pop}/{tr}.sumstats.tsv.gz")
    summ = pd.DataFrame([{"population": pop, "n": n, "n_cases": int(disease.sum()), "trait_mean": trait_q.mean(),
                          "trait_sd": trait_q.std(ddof=1), "mean_age": covs.age.mean(),
                          "male_fraction": covs.male.mean(), "built_with": "Python"}])
    sg.write_tsv(summ, f"{out}/{pop}/summary.tsv")
    print(summ.to_string(index=False))
