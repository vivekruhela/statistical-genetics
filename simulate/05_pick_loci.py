#!/usr/bin/env python
"""Step 5 (Python): genome-wide cis-eQTL scan on the real expression, then pick the teaching loci.

Usage: python simulate/05_pick_loci.py data/work/expression_real data/work/plink data/work/loci
R twin: simulate/05_pick_loci.R (same outputs, same numbers).
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "lib")
import statgen as sg  # noqa: E402

expr_dir, plink_dir, out = sys.argv[1:4]
cfg = sg.read_config()

samples = pd.read_csv(f"{expr_dir}/samples.tsv", sep="\t")
genes = pd.read_csv(f"{expr_dir}/genes.tsv", sep="\t")
expr = pd.read_csv(f"{expr_dir}/expression_int.tsv.gz", sep="\t", index_col=0)[samples["sample"]]
epcs = pd.read_csv(f"{expr_dir}/expression_pcs.tsv", sep="\t", index_col=0)
gpcs = pd.read_csv(f"{plink_dir}/cohort_pca.eigenvec", sep="\t", index_col=0)
no_geno = sorted(set(samples["sample"]) - set(gpcs.index))
if no_geno:
    sys.exit("no genotypes (so no genotype PCs) for: " + ", ".join(no_geno) + "\nRerun from step 0 so step 3 leaves them out.")
covar = np.column_stack([
    (samples["sex"] == "male").astype(float).to_numpy(),
    gpcs.loc[samples["sample"]].iloc[:, :cfg["genotype"]["n_pcs"]].to_numpy(),
    epcs.loc[samples["sample"]].iloc[:, :cfg["expression"]["n_expression_pcs"]].to_numpy(),
])
eres = pd.DataFrame(sg.residualise(expr.to_numpy(), covar), index=expr.index)

hw = cfg["loci"]["half_window_bp"]
scan = []
for f in sorted(glob.glob(f"{plink_dir}/scan/*.traw")):
    info, g, ids = sg.read_traw(f)
    g = g[:, [ids.index(s) for s in samples["sample"]]]
    chr_name = info["chr"].iloc[0]
    gi = genes[genes.chr == chr_name].sort_values("gene_id")
    if gi.empty:
        continue
    gres = sg.residualise(g, covar)
    res = sg.cis_scan(eres.loc[gi.gene_id].to_numpy(), gi.gene_id.to_list(), gi.tss.to_numpy(),
                      gres, info["id"].to_list(), info["pos"].to_numpy(), hw, covar.shape[1])
    print(f"{chr_name}: {len(res)} genes scanned")
    scan.append(res)
scan = pd.concat(scan)
scan = genes[["gene_id", "gene_name", "chr", "tss"]].merge(scan, on="gene_id")
scan = scan.sort_values(["lead_mlog10p", "gene_id"], ascending=[False, True]).reset_index(drop=True)
sg.write_tsv(scan, f"{out}/cis_scan_lead.tsv")

lc = cfg["loci"]
cand = scan[(scan.lead_mlog10p >= lc["min_lead_mlog10p"]) & (scan.n_cis >= lc["min_variants"]) &
            (scan.n_cis <= lc["max_variants"])].copy()
cand["start"] = np.maximum(1, cand.tss - hw)
cand["end"] = cand.tss + hw
excl = pd.DataFrame(lc["exclude_regions"])
bad = [((excl.chr == r.chr) & (excl.start <= r.end) & (excl.end >= r.start)).any() for r in cand.itertuples()]
cand = cand[~np.array(bad, dtype=bool)]
chosen = []
for r in cand.itertuples():
    if len(chosen) == lc["n_loci"]:
        break
    same = [c for c in chosen if c.chr == r.chr]
    if len(same) >= lc["max_loci_per_chr"]:
        continue
    if any(c.start - lc["min_gap_bp"] <= r.end and c.end + lc["min_gap_bp"] >= r.start for c in same):
        continue
    chosen.append(r)
if len(chosen) < lc["n_loci"]:
    sys.exit(f"only {len(chosen)} loci passed the filters")
scen = [s["name"] for s in cfg["scenarios"]]
loci = pd.DataFrame({
    "locus": [f"L{i + 1:02d}" for i in range(len(chosen))],
    "scenario": scen[:len(chosen)],
    "chr": [c.chr for c in chosen], "start": [c.start for c in chosen], "end": [c.end for c in chosen],
    "anchor_gene_id": [c.gene_id for c in chosen], "anchor_gene": [c.gene_name for c in chosen],
    "anchor_tss": [c.tss for c in chosen], "anchor_lead_id": [c.lead_id for c in chosen],
    "anchor_lead_mlog10p": [c.lead_mlog10p for c in chosen],
    "anchor_second_mlog10p": [c.second_mlog10p for c in chosen],
    "n_scan_variants": [c.n_cis for c in chosen],
})
sg.write_tsv(loci, f"{out}/loci.tsv")
print(loci[["locus", "scenario", "chr", "start", "anchor_gene", "anchor_lead_mlog10p"]])
