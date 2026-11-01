#!/usr/bin/env python
"""Step 7 (Python): plant the simulated layers on the real cohort and write the truth.

Python twin of simulate/07_simulate_layers.R: the same model, rules and outputs, built with
NumPy's random numbers. The released dataset comes from the R script; this version shows the
same build in Python and can make an independent replicate cohort (same design, new draws).

Usage: python simulate/07_simulate_layers.py data/work data/cohort data/truth data/meta
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, "lib")
import statgen as sg  # noqa: E402

work, cohort_dir, truth_dir, meta_dir = sys.argv[1:5]
cfg = sg.read_config()
ef, cv, me = cfg["effects"], cfg["causal_variants"], cfg["methylation"]
rng = np.random.default_rng(cfg["seed"])

samples = pd.read_csv(f"{work}/expression_real/samples.tsv", sep="\t")
n = len(samples)
is_eur = (samples["region"] == "EUR").to_numpy()
loci = pd.read_csv(f"{work}/loci/loci.tsv", sep="\t")

# ---- covariates (sex is real) -------------------------------------------------------
covs = sg.simulate_covariates(n, rng, male=(samples["sex"] == "male").astype(int).to_numpy())
inputs = sg.standard_covariates(covs.age, covs.male, covs.bmi)


def maf(p):
    return np.minimum(p, 1 - p)


def pick1(idx, what, locus):
    idx = np.flatnonzero(idx)
    if not len(idx):
        sys.exit(f"{locus}: no variant fits the rule for {what}")
    return int(idx[rng.integers(len(idx))])


edges, nodes, genes_sim, cpgs, variants = [], [], [], [], []


def add_edge(locus, scenario, frm, to, r2, sign, center=0.0, var_ref=1.0, role=""):
    edges.append(dict(locus=locus, scenario=scenario, **{"from": frm}, to=to, role=role, r2=r2, sign=sign,
                      center=center, var_ref=var_ref, weight=sign * np.sqrt(r2 / var_ref)))


def rsign():
    return float(rng.choice([-1, 1]))


for i, L in enumerate(loci.itertuples(), start=1):
    lid, sc = L.locus, L.scenario
    info, g, ids = sg.read_traw(f"{cohort_dir}/genotypes/loci/{lid}.traw")
    g = g[:, [ids.index(s) for s in samples["sample"]]]
    info["af"] = g.mean(1) / 2
    info["af_eur"] = g[:, is_eur].mean(1) / 2
    info["af_yri"] = g[:, ~is_eur].mean(1) / 2
    tss = L.anchor_tss + ef["teaching_gene_offset_bp"]
    if tss > L.end:
        tss = L.anchor_tss - ef["teaching_gene_offset_bp"]
    snv = (info.ref.str.len() == 1).to_numpy() & (info.alt.str.len() == 1).to_numpy()
    ok = snv & (maf(info.af_eur) >= cv["min_maf_eur"]).to_numpy() & (maf(info.af_yri) >= cv["min_maf_yri"]).to_numpy()
    near = (np.abs(info.pos - tss) <= ef["cis_window_bp"]).to_numpy()
    ge = sg.unit_rows(g[:, is_eur])

    def r2_eur(a):  # variants monomorphic in EUR get 0, but they already fail the MAF rules
        return (ge @ ge[a]) ** 2

    if sc == "ancestry":
        A = pick1(snv & near & (maf(info.af_eur) >= cv["ancestry_maf_eur"][0]).to_numpy()
                  & (maf(info.af_eur) <= cv["ancestry_maf_eur"][1]).to_numpy()
                  & (maf(info.af_yri) >= cv["ancestry_maf_yri_min"]).to_numpy(), "A (ancestry)", lid)
    elif sc == "gwas_only":
        A = pick1(ok, "A", lid)
    else:
        A = pick1(ok & near, "A", lid)
    rA = np.nan_to_num(r2_eur(A), nan=1.0)
    B = C = None
    if sc in ("two_signals", "capstone"):
        B = pick1(ok & near & (rA < cv["low_ld_r2_max"]), "B", lid)
    if sc == "gwas_two":
        B = pick1(ok & (rA < cv["low_ld_r2_max"]), "B", lid)
    if sc == "linkage":
        B = pick1(ok & (np.abs(info.pos - info.pos[A]) > 1000).to_numpy()
                  & (rA >= cv["linkage_r2"][0]) & (rA <= cv["linkage_r2"][1]), "B (linkage)", lid)
    if sc == "capstone":
        rB = np.nan_to_num(r2_eur(B), nan=1.0)
        C = pick1(ok & (rA < cv["low_ld_r2_max"]) & (rB < cv["low_ld_r2_max"]), "C", lid)
    roles = {k: v for k, v in {"A": A, "B": B, "C": C}.items() if v is not None}

    def vnode(k):
        return "var:" + info.id[k]

    def vcenter(k):
        return 2 * info.af[k]

    def vvar(k):
        return 2 * info.af[k] * (1 - info.af[k])

    for r, k in roles.items():
        inputs[vnode(k)] = g[k]
        variants.append({"locus": lid, "scenario": sc, "role": r,
                         **info.loc[k, ["id", "chr", "pos", "ref", "alt", "af", "af_eur", "af_yri"]].to_dict()})

    gene = f"SIMG{i:02d}"
    expr, fac = f"expr:{gene}", f"factor:{lid}"
    genes_sim.append(dict(gene_id=gene, gene_name=gene, gene_type="teaching", chr=L.chr,
                          start=tss, end=tss + 20000, strand="+", tss=tss))
    cpg_ids = [f"simcg_{lid}_{j + 1:02d}" for j in range(me["cpgs_per_locus"])]
    pos = np.sort(tss + rng.choice(np.arange(-2000, 2001), me["cpgs_per_locus"], replace=False))
    cpgs.append(pd.DataFrame({"cpg_id": cpg_ids, "chr": L.chr, "pos": pos, "locus": lid, "cluster": gene}))
    nodes.append(dict(node=fac, type="factor", order=1))
    nodes.append(dict(node=expr, type="expression", order=3))
    for cg in cpg_ids:
        add_edge(lid, sc, fac, f"cpg:{cg}", me["regional_factor_r2"], 1.0, role="regional")
    cpg1 = f"cpg:{cpg_ids[0]}"

    def eq(k, r2, role):
        add_edge(lid, sc, vnode(k), expr, r2, rsign(), vcenter(k), vvar(k), role)

    def direct(k, r2, to, role):
        add_edge(lid, sc, vnode(k), to, r2, rsign(), vcenter(k), vvar(k), role)

    def via_expr(target_r2, chain_r2, to):
        add_edge(lid, sc, expr, to, target_r2 / chain_r2, rsign(), role="expr_to_trait")

    def mediated_chain(k):
        add_edge(lid, sc, vnode(k), cpg1, ef["mqtl_r2"], rsign(), vcenter(k), vvar(k), "mqtl")
        add_edge(lid, sc, cpg1, expr, ef["cpg_to_expr_r2"], -1.0, role="cpg_to_expr")
        return ef["mqtl_r2"] * ef["cpg_to_expr_r2"]

    if sc == "mediated":
        via_expr(ef["trait_r2"], mediated_chain(A), "trait_q")
    elif sc == "two_signals":
        eq(A, ef["eqtl_r2"], "eqtl"); eq(B, ef["eqtl_second_r2"], "eqtl")  # noqa: E702
        via_expr(ef["trait_r2"], ef["eqtl_r2"], "trait_q")
    elif sc == "linkage":
        eq(A, ef["eqtl_r2"], "eqtl"); direct(B, ef["trait_r2"], "trait_q", "direct")  # noqa: E702
    elif sc == "pleiotropy":
        eq(A, ef["eqtl_r2"], "eqtl"); direct(A, ef["trait_r2"], "trait_q", "direct")  # noqa: E702
    elif sc == "gwas_only":
        direct(A, ef["trait_strong_r2"], "trait_q", "direct")
    elif sc == "eqtl_only":
        eq(A, ef["eqtl_r2"], "eqtl")
    elif sc == "gwas_two":
        eq(A, ef["eqtl_r2"], "eqtl"); via_expr(ef["trait_r2"], ef["eqtl_r2"], "trait_q")  # noqa: E702
        direct(B, ef["trait_r2"], "trait_q", "direct")
    elif sc == "binary":
        eq(A, ef["eqtl_r2"], "eqtl"); via_expr(ef["disease_r2"], ef["eqtl_r2"], "disease_liability")  # noqa: E702
    elif sc == "ancestry":
        eq(A, ef["eqtl_r2"], "eqtl"); via_expr(ef["trait_r2"], ef["eqtl_r2"], "trait_q")  # noqa: E702
    elif sc == "capstone":
        ch = mediated_chain(A)
        eq(B, ef["eqtl_second_r2"], "eqtl")
        via_expr(ef["trait_r2"], ch, "trait_q")
        direct(C, ef["trait_r2"], "trait_q", "direct")
    else:
        sys.exit(f"unknown scenario {sc}")

nnull = me["null_cpgs"]
li = rng.integers(0, len(loci), nnull)
cpgs.append(pd.DataFrame({
    "cpg_id": [f"simcg_null_{j + 1:04d}" for j in range(nnull)], "chr": loci.chr.to_numpy()[li],
    "pos": np.round(loci.start.to_numpy()[li] + rng.random(nnull) * (loci.end - loci.start).to_numpy()[li]).astype(int),
    "locus": loci.locus.to_numpy()[li], "cluster": "none"}))
cpgs = pd.concat(cpgs, ignore_index=True)
null = (cpgs.cluster == "none").to_numpy()
cpgs["m_mean"] = np.where(null, rng.uniform(-3, 3, len(cpgs)), rng.uniform(-4, -1, len(cpgs)))
cpgs["m_scale"] = rng.uniform(0.3, 0.9, len(cpgs))
cpgs = cpgs.sort_values(["chr", "pos"]).reset_index(drop=True)
age_cpg = cpgs.cpg_id[rng.random(len(cpgs)) < me["age_cpg_fraction"]]
for cg in age_cpg:
    add_edge(cpgs.locus[cpgs.cpg_id == cg].iloc[0], "methylation_age", "cov:age", f"cpg:{cg}", me["age_r2"], rsign(), role="age")
nodes.append(pd.DataFrame({"node": "cpg:" + cpgs.cpg_id, "type": "cpg", "order": 2}))

for to in ("trait_q", "disease_liability"):
    for cn, r2 in ef["covariate_r2"].items():
        add_edge("global", "covariates", f"cov:{cn}", to, r2, 1.0, role="covariate")
nodes.append(pd.DataFrame({"node": ["trait_q", "disease_liability"], "type": ["trait", "liability"], "order": 4}))

edges = pd.DataFrame(edges)
nodes = pd.concat([pd.DataFrame([x]) if isinstance(x, dict) else x for x in nodes], ignore_index=True)
r2_in = edges.groupby("to")["r2"].sum()
nodes["r2_in"] = nodes.node.map(r2_in).fillna(0.0)
if (nodes.r2_in >= 1).any():
    sys.exit("variance explained >= 1 for: " + ", ".join(nodes.node[nodes.r2_in >= 1]))
nodes["resid_sd"] = np.sqrt(1 - nodes.r2_in)
nodes = nodes.sort_values(["order", "node"]).reset_index(drop=True)

# ---- simulate the graph for the cohort ------------------------------------------------
vals = sg.simulate_dag(nodes, edges, inputs, n, rng)
threshold = stats.norm.ppf(1 - ef["disease_prevalence"])
clinical = pd.DataFrame({
    "sample": samples["sample"], "population": samples.population, "superpopulation": samples.region,
    "sex": samples.sex, "age": covs.age.to_numpy(), "bmi": covs.bmi.to_numpy(),
    "trait_q": np.round(50 + 10 * vals["trait_q"], 2),
    "disease": (vals["disease_liability"] > threshold).astype(int)})
print(f"disease cases: {clinical.disease.sum()} of {n}")
meth = pd.DataFrame(np.vstack([cpgs.m_mean[k] + cpgs.m_scale[k] * vals["cpg:" + c] for k, c in enumerate(cpgs.cpg_id)]),
                    index=cpgs.cpg_id, columns=samples["sample"])

# ---- expression: teaching genes and the disease spike -----------------------------------
cmat = pd.read_csv(f"{work}/expression_real/counts_all.tsv.gz", sep="\t", index_col=0)[samples["sample"]]
sf = pd.read_csv(f"{work}/expression_real/size_factors.tsv", sep="\t", index_col=0)["size_factor"].loc[samples["sample"]].to_numpy()
genes_sim = pd.DataFrame(genes_sim)
size = 1 / ef["teaching_gene_dispersion"]
sim_counts = []
for gn in genes_sim.gene_id:
    mu = 2 ** (ef["teaching_gene_log2_mean"] + ef["teaching_gene_log2_sd"] * vals[f"expr:{gn}"]) * sf
    sim_counts.append(rng.negative_binomial(size, size / (size + mu)))
sim_counts = pd.DataFrame(np.vstack(sim_counts), index=genes_sim.gene_id, columns=samples["sample"])

genes = sg.read_genes(f"{meta_dir}/gencode.v39.annotation.gtf.gz")
genes_real = pd.read_csv(f"{work}/expression_real/genes.tsv", sep="\t")
in_loci = np.array([((loci.chr == r.chr) & (loci.start - 1e6 <= r.tss) & (loci.end + 1e6 >= r.tss)).any()
                    for r in genes_real.itertuples()])
well = genes_real.gene_id.isin(cmat.index[cmat.median(axis=1) >= 30]).to_numpy()
pool = np.sort(genes_real.gene_id[~in_loci & well].to_numpy())
de = cfg["de"]
de_ids = rng.choice(pool, de["n_up"] + de["n_down"], replace=False)
de_tab = pd.DataFrame({"gene_id": de_ids, "direction": ["up"] * de["n_up"] + ["down"] * de["n_down"]})
de_tab["log2fc"] = np.where(de_tab.direction == "up", 1, -1) * np.maximum(0.2, rng.normal(de["log2fc_mean"], de["log2fc_sd"], len(de_tab)))
cases = np.flatnonzero(clinical.disease.to_numpy() == 1)
cvals = cmat.to_numpy().copy()
row_of = {gid: k for k, gid in enumerate(cmat.index)}
for gid, fc2 in zip(de_tab.gene_id, de_tab.log2fc):
    x = cvals[row_of[gid], cases]
    fc = 2 ** fc2
    cvals[row_of[gid], cases] = x + rng.poisson(x * (fc - 1)) if fc > 1 else rng.binomial(x.astype(int), fc)
cmat = pd.DataFrame(cvals, index=cmat.index, columns=cmat.columns)
de_tab = de_tab.merge(genes_real[["gene_id", "gene_name", "chr", "tss"]], on="gene_id").sort_values(["direction", "gene_id"])

all_counts = pd.concat([cmat, sim_counts])
all_genes = pd.concat([genes, genes_sim], ignore_index=True)
cfg["expression"]["gene_types"] = cfg["expression"]["gene_types"] + ["teaching"]
norm = sg.normalise_expression(all_counts, all_genes, cfg)
pcs = sg.sample_pcs(norm["int"], cfg["expression"]["n_expression_pcs"])

# ---- write -------------------------------------------------------------------------------
print("writing outputs")
sg.write_tsv(clinical, f"{cohort_dir}/clinical/clinical.tsv")
sg.write_tsv(samples[["sample", "population", "region", "sex"]].rename(columns={"region": "superpopulation"}),
             f"{cohort_dir}/samples.tsv")
sg.write_tsv(sg.matrix_table(meth.round(4), "cpg_id"), f"{cohort_dir}/methylation/methylation_m.tsv.gz")
sg.write_tsv(cpgs.assign(promoter_of=np.where(cpgs.cluster == "none", "", cpgs.cluster))[["cpg_id", "chr", "pos", "promoter_of"]],
             f"{cohort_dir}/methylation/cpgs.tsv")
ex = f"{cohort_dir}/expression"
sg.write_tsv(sg.matrix_table(all_counts[all_counts.index.isin(all_genes.gene_id)], "gene_id"), f"{ex}/counts.tsv.gz")
sg.write_tsv(all_genes.set_index("gene_id").loc[norm["int"].index].reset_index(), f"{ex}/genes.tsv")
sg.write_tsv(sg.matrix_table(norm["int"], "gene_id"), f"{ex}/expression_int.tsv.gz")
sg.write_tsv(sg.matrix_table(norm["log2norm"].round(6), "gene_id"), f"{ex}/expression_log2norm.tsv.gz")
sg.write_tsv(pd.DataFrame({"sample": norm["size_factors"].index, "size_factor": norm["size_factors"].to_numpy()}), f"{ex}/size_factors.tsv")
sg.write_tsv(sg.matrix_table(pcs, "sample"), f"{ex}/expression_pcs.tsv")
sg.write_tsv(loci, f"{cohort_dir}/loci.tsv")

variants = pd.DataFrame(variants)
sg.write_tsv(edges, f"{truth_dir}/edges.tsv")
sg.write_tsv(nodes, f"{truth_dir}/nodes.tsv")
sg.write_tsv(variants, f"{truth_dir}/causal_variants.tsv")
sg.write_tsv(de_tab, f"{truth_dir}/de_genes.tsv")
sg.write_tsv(cpgs, f"{truth_dir}/cpg_parameters.tsv")
sg.write_tsv(genes_sim, f"{truth_dir}/teaching_genes.tsv")
latent = pd.DataFrame({k: v for k, v in vals.items() if k.startswith(("expr", "trait", "disease"))})
latent.insert(0, "sample", samples["sample"])
sg.write_tsv(latent, f"{truth_dir}/latent_small_cohort.tsv.gz")
with open(f"{truth_dir}/model.json", "w") as fh:
    json.dump({"seed": cfg["seed"], "disease_threshold": float(threshold),
               "disease_prevalence": ef["disease_prevalence"], "trait_q": {"offset": 50, "scale": 10},
               "teaching_gene": {"log2_mean": ef["teaching_gene_log2_mean"], "log2_sd": ef["teaching_gene_log2_sd"],
                                 "dispersion": ef["teaching_gene_dispersion"]},
               "built_with": "Python"}, fh, indent=2)

v_edges = edges[edges["from"].str.startswith("var:")].copy()
v_edges["id"] = v_edges["from"].str.replace("^var:", "", regex=True)
truth = v_edges.merge(variants.drop_duplicates("id")[["id", "role", "chr", "pos", "ref", "alt", "af", "af_eur", "af_yri"]],
                      on="id", suffixes=("_edge", "_variant"))
down = edges[edges.role.isin(["cpg_to_expr", "expr_to_trait"])]


def follow(to):
    """Path down the cpg -> expression -> trait chain, its endpoint and the product of weights."""
    out, cur, gain = to, to, 1.0
    while True:
        nx = down[down["from"] == cur]
        if nx.empty:
            return out, cur, gain
        gain *= nx.weight.iloc[0]
        cur = nx.to.iloc[0]
        out += " -> " + cur


fol = truth.to.map(follow)
truth["path"] = fol.str[0]
truth["endpoint"] = fol.str[1]
truth["endpoint_effect"] = truth.weight * fol.str[2]
truth = pd.DataFrame({
    "locus": truth.locus, "scenario": truth.scenario, "variant_id": truth.id, "variant_role": truth.role_variant,
    "chr": truth.chr, "pos": truth.pos, "ref": truth.ref, "effect_allele": truth.alt, "af": truth.af,
    "af_eur": truth.af_eur, "af_yri": truth.af_yri, "acts_on": truth.to, "effect_type": truth.role_edge,
    "path": truth.path, "r2_target": truth.r2, "per_allele_effect": truth.weight,
    "endpoint": truth.endpoint, "endpoint_effect": truth.endpoint_effect,
}).sort_values(["locus", "variant_role", "acts_on"])
sg.write_tsv(truth, f"{truth_dir}/truth.tsv")
print(truth[["locus", "scenario", "variant_role", "variant_id", "effect_type", "path"]].to_string(index=False))
