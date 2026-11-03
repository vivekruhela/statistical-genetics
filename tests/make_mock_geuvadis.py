#!/usr/bin/env python
"""Write a small GEUVADIS-shaped count table for testing steps 3 to 7 without the real download.

Real GEUVADIS needs EBI access. This mock has the same layout (TargetID, Gene_Symbol, Chr, Coord,
then one column per library), about 460 people, duplicate libraries for a few of them, lab batch
effects, and planted cis-eQTLs on chr22 read from the real step 1 genotypes.

Usage: python tests/make_mock_geuvadis.py data/meta data/genotypes/chr22.bcf data/geuvadis
"""
import gzip
import os
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "lib")
import statgen as sg  # noqa: E402

meta_dir, bcf, out = sys.argv[1:4]
rng = np.random.default_rng(7)
cohort = pd.read_csv(f"{meta_dir}/cohort_samples.tsv", sep="\t")
in_bcf = set(subprocess.run(["bcftools", "query", "-l", bcf], capture_output=True, text=True).stdout.split())
cohort = cohort[cohort["sample"].isin(in_bcf)]
eur = cohort[cohort.region == "EUR"]["sample"].to_numpy()
yri = cohort[cohort.region == "AFR"]["sample"].to_numpy()
people = np.concatenate([rng.choice(eur, 373, replace=False), rng.choice(yri, 89, replace=False)])

genes = sg.read_genes(f"{meta_dir}/gencode.v39.annotation.gtf.gz")
genes = genes[genes.gene_type.isin(["protein_coding", "lncRNA"])]
chr22 = genes[genes.chr == "chr22"]
other = genes[(genes.chr != "chr22") & genes.chr.str.match(r"^chr\d+$")].sample(1500, random_state=1)
g = pd.concat([chr22, other]).reset_index(drop=True)

n = len(people)
log_mu = rng.normal(5, 2, len(g))
low = rng.random(len(g)) < 0.25
log_mu[low] = rng.normal(0, 1, low.sum())
lab = rng.integers(0, 7, n)
lab_eff = rng.normal(0, 0.15, (len(g), 7))
depth = rng.lognormal(0, 0.25, n)
eta = log_mu[:, None] + lab_eff[:, lab] + rng.normal(0, 0.3, (len(g), n)) + np.log(depth)[None, :]

# planted cis-eQTLs: 60 chr22 genes, one variant within 200 kb of the TSS
samples_arg = ",".join(people)
eg = rng.choice(np.flatnonzero(g.chr.to_numpy() == "chr22"), 60, replace=False)
for i in eg:
    tss = int(g.loc[i, "tss"])
    q = subprocess.run(["bcftools", "query", "-s", samples_arg, "--force-samples", "-r",
                        f"chr22:{max(1, tss - 200000)}-{tss + 200000}", "-i", "MAF>0.1",
                        "-f", "[%GT\t]\n", bcf], capture_output=True, text=True).stdout.strip().split("\n")
    if not q or q == [""]:
        continue
    row = q[rng.integers(len(q))].rstrip("\t").split("\t")
    dos = np.array([sum(int(a) for a in gt.replace("|", "/").split("/")) for gt in row], float)
    eta[i] += rng.choice([-1, 1]) * rng.uniform(0.2, 0.5) * (dos - dos.mean())

counts = rng.negative_binomial(10, 10 / (10 + np.exp(eta)))
libs = [f"{p}.1.M_1111{lab[j]:02d}_{j % 8}" for j, p in enumerate(people)]
df = pd.DataFrame(counts, columns=libs)
for j in range(5):  # a shallower second library for five people
    df[f"{people[j]}.2.M_120000_1"] = (counts[:, j] * 0.3).round().astype(int)
info = pd.DataFrame({"TargetID": g.gene_id + ".1", "Gene_Symbol": g.gene_id + ".1",
                     "Chr": g.chr.str.replace("chr", ""), "Coord": g.tss})
os.makedirs(out, exist_ok=True)
with gzip.open(f"{out}/GD660.GeneQuantCount.txt.gz", "wt") as fh:
    pd.concat([info, df], axis=1).to_csv(fh, sep="\t", index=False)
print(f"mock GEUVADIS: {len(g)} genes x {df.shape[1]} libraries ({n} people)")
