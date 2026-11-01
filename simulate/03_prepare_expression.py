#!/usr/bin/env python
"""Step 3 (Python): match GEUVADIS libraries to genotyped people, filter and normalise expression.

Usage: python simulate/03_prepare_expression.py data/geuvadis data/meta data/work/expression_real
R twin: simulate/03_prepare_expression.R (same outputs, same numbers).
"""
import os
import sys

import pandas as pd

sys.path.insert(0, "lib")
import statgen as sg  # noqa: E402

geuv_dir, meta_dir, out = sys.argv[1:4]
cfg = sg.read_config()

counts = sg.read_geuvadis_counts(os.path.join(geuv_dir, "GD660.GeneQuantCount.txt.gz"))
cohort = pd.read_csv(os.path.join(meta_dir, "cohort_samples.tsv"), sep="\t")
keep = sorted(set(counts.columns) & set(cohort["sample"]))
print(f"GEUVADIS individuals: {counts.shape[1]}; with genotypes and metadata: {len(keep)}")
counts = counts[keep]
samples = cohort.set_index("sample").loc[keep].reset_index()
print(samples.groupby(["population", "sex"]).size())

genes = sg.read_genes(os.path.join(meta_dir, "gencode.v39.annotation.gtf.gz"))
norm = sg.normalise_expression(counts, genes, cfg)
pcs = sg.sample_pcs(norm["int"], cfg["expression"]["n_expression_pcs"])
print(f"genes kept after filtering: {norm['int'].shape[0]} of {counts.shape[0]}")

os.makedirs(out, exist_ok=True)
sg.write_tsv(samples, f"{out}/samples.tsv")
with open(f"{out}/samples.txt", "w") as fh:
    fh.write("\n".join(samples["sample"]) + "\n")
sg.write_tsv(sg.matrix_table(counts[counts.index.isin(genes.gene_id)], "gene_id"), f"{out}/counts_all.tsv.gz")
sg.write_tsv(genes.set_index("gene_id").loc[norm["int"].index].reset_index(), f"{out}/genes.tsv")
sg.write_tsv(pd.DataFrame({"sample": norm["size_factors"].index, "size_factor": norm["size_factors"].to_numpy()}),
             f"{out}/size_factors.tsv")
sg.write_tsv(sg.matrix_table(norm["int"], "gene_id"), f"{out}/expression_int.tsv.gz")
sg.write_tsv(sg.matrix_table(pcs, "sample"), f"{out}/expression_pcs.tsv")
