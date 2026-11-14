"""Shared helpers for the data build (Python track).

lib/statgen.R is the R twin: the same functions, the same numbers. Every function here is
deterministic; all randomness lives in the simulate/ scripts.
"""
import gzip
import json
import math
import os
import re

import numpy as np
import pandas as pd
from scipy import stats


def _merge(base, over):
    out = dict(base)
    for k, v in over.items():
        out[k] = _merge(base[k], v) if isinstance(v, dict) and isinstance(base.get(k), dict) else v
    return out


def read_config(path=None):
    """Read the build config; a config may name a "_base" config and override some values."""
    path = path or os.environ.get("STATGEN_CONFIG", "config/cohort.json")
    with open(path) as fh:
        cfg = json.load(fh)
    if "_base" in cfg:
        cfg = _merge(read_config(cfg.pop("_base")), cfg)
    return cfg


# ---- genes -------------------------------------------------------------------

def read_genes(gtf):
    """One row per gene from a GENCODE GTF: id without version, name, type, position, TSS."""
    rows = []
    pat = {k: re.compile(k + r' "([^"]+)"') for k in ("gene_id", "gene_name", "gene_type")}
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "gene":
                continue
            full_id = pat["gene_id"].search(f[8]).group(1)
            if "_PAR_Y" in full_id:
                continue
            rows.append((full_id.split(".")[0], pat["gene_name"].search(f[8]).group(1),
                         pat["gene_type"].search(f[8]).group(1), f[0], int(f[3]), int(f[4]), f[6]))
    g = pd.DataFrame(rows, columns=["gene_id", "gene_name", "gene_type", "chr", "start", "end", "strand"])
    g["tss"] = np.where(g["strand"] == "+", g["start"], g["end"])
    return g.drop_duplicates("gene_id").reset_index(drop=True)


# ---- expression --------------------------------------------------------------

def read_geuvadis_counts(path):
    """GEUVADIS count table -> DataFrame genes x individuals (rounded counts).

    Library columns look like HG00096.1.M_111124_6; the individual is the part before the
    first dot. Where one person has several libraries we keep the deepest one.
    """
    d = pd.read_csv(path, sep="\t")
    info_cols = [c for c in ("TargetID", "Gene_Symbol", "Chr", "Coord") if c in d.columns]
    libs = [c for c in d.columns if c not in info_cols]
    m = d[libs].astype(float)
    m.index = d.iloc[:, 0].astype(str).str.replace(r"\..*", "", regex=True)
    m.index.name = "gene_id"
    m = m[~m.index.duplicated()]
    ind = pd.Series([c.split(".")[0] for c in libs])
    depth = m.sum(axis=0).to_numpy()
    order = sorted(range(len(libs)), key=lambda j: (ind[j], -depth[j]))
    keep, seen = [], set()
    for j in order:
        if ind[j] not in seen:
            seen.add(ind[j])
            keep.append(j)
    m = m.iloc[:, keep]
    m.columns = ind[keep].to_list()
    m = m[sorted(m.columns)]
    return m.round()


def inverse_normal(x):
    """Rank-based inverse normal transform (Blom offsets), ties get average ranks."""
    r = stats.rankdata(x, method="average")
    return stats.norm.ppf((r - 3 / 8) / (len(x) + 1 / 4))


def size_factors(counts):
    """DESeq2-style median-of-ratios size factors, using genes with no zero count."""
    c = counts.to_numpy()
    pos = (c == 0).sum(axis=1) == 0
    logc = np.log(c[pos])
    return pd.Series(np.exp(np.median(logc - logc.mean(axis=1, keepdims=True), axis=0)),
                     index=counts.columns)


def normalise_expression(counts, genes, cfg):
    """Filter genes, normalise and transform. Returns a dict of DataFrames (genes x samples)."""
    ecfg = cfg["expression"]
    autosomes = [f"chr{i}" for i in range(1, 23)]
    ann = genes[genes.gene_type.isin(ecfg["gene_types"]) & genes.chr.isin(autosomes)]
    counts = counts[counts.index.isin(ann.gene_id)]
    cpm = counts / counts.sum(axis=0) * 1e6
    ok = (counts >= ecfg["min_count"]) & (cpm >= ecfg["min_cpm"])
    counts = counts[ok.mean(axis=1) >= ecfg["min_sample_fraction"]].sort_index()
    sf = size_factors(counts)
    log2norm = np.log2(counts / sf + 1)
    intm = log2norm.apply(lambda row: pd.Series(inverse_normal(row.to_numpy()), index=row.index), axis=1)
    return {"counts": counts, "size_factors": sf, "log2norm": log2norm, "int": intm}


# ---- linear algebra -----------------------------------------------------------

def sample_pcs(mat, k):
    """Top k PC scores of the samples (columns of mat); largest-magnitude score is positive."""
    x = mat.to_numpy().T
    x = x - x.mean(axis=0)
    u, d, _ = np.linalg.svd(x, full_matrices=False)
    scores = u[:, :k] * d[:k]
    for j in range(k):
        if scores[np.argmax(np.abs(scores[:, j])), j] < 0:
            scores[:, j] = -scores[:, j]
    return pd.DataFrame(scores, index=mat.columns, columns=[f"PC{j + 1}" for j in range(k)])


def residualise(y, covar):
    """Residuals of each row of y (features x samples) after regressing on covariates."""
    c = np.column_stack([np.ones(y.shape[1]), np.asarray(covar, dtype=float)])
    q, _ = np.linalg.qr(c)
    yt = np.asarray(y, dtype=float).T
    return (yt - q @ (q.T @ yt)).T


def unit_rows(x):
    """Centre each row and scale it to unit length; constant rows become all zero."""
    x = x - x.mean(axis=1, keepdims=True)
    nrm = np.sqrt((x ** 2).sum(axis=1, keepdims=True))
    nrm[nrm == 0] = np.inf
    return x / nrm


def cor_mlog10p(r, df):
    """-log10 two-sided p-value of a correlation r with df residual degrees of freedom."""
    r = np.minimum(np.abs(r), 1 - 1e-15)
    t = r * np.sqrt(df / (1 - r ** 2))
    return -(stats.t.logsf(t, df) + math.log(2)) / math.log(10)


# ---- genotypes ----------------------------------------------------------------

def read_traw(path):
    """PLINK 2 .traw -> (variant table, ALT dosage array variants x samples, sample IDs)."""
    d = pd.read_csv(path, sep="\t")
    d.columns = ["chr", "id", "cm", "pos", "counted", "other"] + list(d.columns[6:])
    samples = [re.sub(r"^[^_]*_", "", c) for c in d.columns[6:]]
    g = np.array(d.iloc[:, 6:].to_numpy(dtype=float), copy=True)
    alt = d["id"].str.split(":").str[-1]
    flip = (d["counted"] != alt).to_numpy()
    g[flip] = 2 - g[flip]
    parts = d["id"].str.split(":", expand=True)
    info = pd.DataFrame({"id": d["id"], "chr": parts[0], "pos": parts[1].astype(int),
                         "ref": parts[2], "alt": parts[3]})
    return info, g, samples


# ---- cis-eQTL scan ---------------------------------------------------------------

def cis_scan(expr_res, gene_ids, gene_tss, geno_res, snp_ids, snp_pos, half_window, n_covar):
    """Lead cis variant per gene, and the best variant after conditioning on the lead."""
    n = expr_res.shape[1]
    df = n - 2 - n_covar
    e = unit_rows(np.asarray(expr_res, dtype=float))
    g = unit_rows(np.asarray(geno_res, dtype=float))
    snp_pos = np.asarray(snp_pos)
    assert np.all(np.diff(snp_pos) >= 0)
    lo = np.searchsorted(snp_pos, np.asarray(gene_tss) - half_window, side="left")
    hi = np.searchsorted(snp_pos, np.asarray(gene_tss) + half_window, side="right")
    out = []
    for i in range(e.shape[0]):
        if hi[i] <= lo[i]:
            continue
        gi = g[lo[i]:hi[i]]
        r = gi @ e[i]
        best = int(np.argmax(np.abs(r)))
        gl = gi[best]
        e2 = e[i] - (e[i] @ gl) * gl
        g2 = gi - np.outer(gi @ gl, gl)
        r2 = unit_rows(g2) @ (e2 / np.sqrt(e2 @ e2))
        r2[best] = 0
        best2 = int(np.argmax(np.abs(r2)))
        out.append({
            "gene_id": gene_ids[i], "n_cis": hi[i] - lo[i],
            "lead_id": snp_ids[lo[i] + best], "lead_r": r[best],
            "lead_mlog10p": float(cor_mlog10p(r[best], df)),
            "second_id": snp_ids[lo[i] + best2],
            "second_mlog10p": float(cor_mlog10p(r2[best2], df - 1)),
        })
    return pd.DataFrame(out)


# ---- simulation helpers ------------------------------------------------------------

def ld_r2(g):
    u = unit_rows(np.asarray(g, dtype=float))
    return (u @ u.T) ** 2


def write_tsv(df, path, index=False):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    df.to_csv(path, sep="\t", index=index)


def matrix_table(m, id_name):
    out = m.copy()
    out.index.name = id_name
    return out.reset_index()


# ---- the causal model (shared by steps 7 and 9) ----------------------------------

def simulate_dag(nodes, edges, inputs, n, rng):
    """Simulate every latent node of the causal graph for n people.

    nodes:  DataFrame node, order, resid_sd (one row per latent node)
    edges:  DataFrame from, to, weight, center; value(to) += weight * (value(from) - center)
    inputs: dict of observed arrays (variant dosages "var:<id>", covariates "cov:<name>")
    Nodes are drawn in `order`, one standard normal vector each.
    """
    vals = dict(inputs)
    nodes = nodes.sort_values(["order", "node"])
    by_to = {k: v for k, v in edges.groupby("to")}
    for nd, sd in zip(nodes["node"], nodes["resid_sd"]):
        v = sd * rng.standard_normal(n)
        if nd in by_to:
            for f, w, c in zip(by_to[nd]["from"], by_to[nd]["weight"], by_to[nd]["center"]):
                v = v + w * (vals[f] - c)
        vals[nd] = v
    return {k: vals[k] for k in nodes["node"]}


def standard_covariates(age, male, bmi):
    """Covariates on the fixed standard scale shared by both cohorts."""
    return {"cov:age": (np.asarray(age) - 47.5) / (55 / math.sqrt(12)),
            "cov:sex": (np.asarray(male) - 0.5) / 0.5,
            "cov:bmi": (np.log(bmi) - math.log(26)) / 0.15}


def simulate_covariates(n, rng, male=None):
    age = np.round(rng.uniform(20, 75, n))
    if male is None:
        male = rng.binomial(1, 0.5, n)
    bmi = np.round(np.exp(rng.normal(math.log(26) + 0.004 * (age - 47.5), 0.14)), 1)
    return pd.DataFrame({"age": age, "male": male, "bmi": bmi})


# ---- the GWAS cohort (step 9) -------------------------------------------------------

def read_haplotypes(path, ids):
    """Template haplotype table from step 8, rows aligned to the variant IDs `ids`."""
    h = pd.read_csv(path, sep="\t", index_col=0)
    return h.loc[list(ids)].to_numpy(dtype=np.int8)


def mosaic_dosage(H, pos, n_people, segment_bp, miscopy, rng, batch=500):
    """New people as mosaics of template haplotypes (copying model, as in HAPGEN).

    Each new haplotype copies one template, switches to a random template at breakpoints that
    fall on average every segment_bp, and miscopies an allele with probability miscopy.
    Returns ALT dosages, variants x n_people (int8).
    """
    m, nt = H.shape
    pos = np.asarray(pos)
    span = pos.max() - pos.min() + 1
    out = np.zeros((m, n_people), dtype=np.int8)
    rows = np.arange(m)
    for b0 in range(0, n_people, batch):
        people = np.arange(b0, min(n_people, b0 + batch))
        nh = 2 * len(people)
        hap = np.empty((m, nh), dtype=np.int8)
        for h in range(nh):
            k = rng.poisson(span / segment_bp)
            brk = np.sort(rng.uniform(pos.min(), pos.max(), k))
            tmpl = rng.integers(0, nt, k + 1)
            hap[:, h] = H[rows, tmpl[np.searchsorted(brk, pos, side="right")]]
        nflip = rng.binomial(hap.size, miscopy)
        if nflip:
            f = rng.choice(hap.size, nflip, replace=False)
            flat = hap.reshape(-1)
            flat[f] = 1 - flat[f]
        out[:, people] = hap[:, 0::2] + hap[:, 1::2]
    return out


def logistic_null(y, X, iters=50):
    beta = np.zeros(X.shape[1])
    for _ in range(iters):
        mu = 1 / (1 + np.exp(-(X @ beta)))
        w = mu * (1 - mu)
        step = np.linalg.solve(X.T @ (X * w[:, None]), X.T @ (y - mu))
        beta += step
        if np.max(np.abs(step)) < 1e-10:
            break
    return 1 / (1 + np.exp(-(X @ beta)))


def assoc_stats(G, y, X, binary=False):
    """Per-variant association statistics adjusting for X (n x q, with intercept).

    Quantitative y: least squares (Frisch-Waugh). Binary y: score test from the null logistic
    fit, with beta = U / I, the one-step approximation to the log odds ratio.
    """
    G = np.asarray(G, dtype=float)
    if not binary:
        q, _ = np.linalg.qr(X)
        yr = y - q @ (q.T @ y)
        gy = G @ yr
        CG = q.T @ G.T
        gg = (G ** 2).sum(1) - (CG ** 2).sum(0)
        gg[gg <= 1e-8] = np.nan  # monomorphic in this cohort
        beta = gy / gg
        df = len(y) - X.shape[1] - 1
        s2 = ((yr ** 2).sum() - beta * gy) / df
        se = np.sqrt(s2 / gg)
    else:
        mu = logistic_null(y, X)
        w = mu * (1 - mu)
        U = G @ (y - mu)
        inv = np.linalg.inv(X.T @ (X * w[:, None]))
        CWG = (X * w[:, None]).T @ G.T
        I = (G ** 2) @ w - (CWG * (inv @ CWG)).sum(0)
        I[I <= 1e-8] = np.nan
        beta = U / I
        se = 1 / np.sqrt(I)
    z = beta / se
    return pd.DataFrame({"beta": beta, "se": se, "z": z,
                         "mlog10p": -(stats.norm.logsf(np.abs(z)) + math.log(2)) / math.log(10)})


# ---- LD and reference panels (post R onward) ---------------------------------------

def ld_matrix(g):
    """Correlation (r, not r^2) between rows of a dosage array, variants x people.

    A variant that does not vary in these people gets NaN: the panel cannot say anything about it.
    """
    u = unit_rows(np.asarray(g, dtype=float))
    r = u @ u.T
    bad = (u ** 2).sum(1) == 0
    r[bad, :] = np.nan
    r[:, bad] = np.nan
    r[np.diag_indices_from(r)] = np.where(bad, np.nan, 1.0)
    return r


def write_ld(r, ids, prefix):
    """An LD matrix on disk: PREFIX.ld.bin (float32, little endian, p x p) and PREFIX.ld.ids."""
    os.makedirs(os.path.dirname(prefix) or ".", exist_ok=True)
    with open(prefix + ".ld.ids", "w") as fh:
        fh.write("\n".join(ids) + "\n")
    np.asarray(r, dtype="<f4").T.tofile(prefix + ".ld.bin")   # column-major, as R writes it


def read_ld(prefix):
    with open(prefix + ".ld.ids") as fh:
        ids = fh.read().split()
    r = np.fromfile(prefix + ".ld.bin", dtype="<f4").astype(float).reshape(len(ids), len(ids)).T
    return r, ids


def ld_check(z, R, n, r_tol=1e-8):
    """Do these z-scores fit this LD matrix?

    For each variant, the z-score its neighbours predict (the conditional mean under
    z ~ N(0, R_s), with R_s = (1 - s) R + s I) and how far the observed z is from it, in standard
    deviations. s is estimated by maximum likelihood, as in susieR::estimate_s_rss(method =
    "null-mle"); a value near 0 means R fits the z-scores well. This is the idea behind
    susieR::kriging_rss and DENTIST.
    """
    from scipy.optimize import minimize_scalar
    z = np.nan_to_num(np.asarray(z, dtype=float))
    z = np.sqrt((n - 1) / (z ** 2 + n - 2)) * z
    d, v = np.linalg.eigh(R)
    d = np.where(d < r_tol, 0.0, d)
    ztv = v.T @ z
    nll = lambda s: 0.5 * np.sum(np.log((1 - s) * d + s)) + 0.5 * np.sum(ztv ** 2 / ((1 - s) * d + s))
    with np.errstate(divide="ignore", invalid="ignore"):
        s = minimize_scalar(nll, bounds=(0, 1), method="bounded", options={"xatol": 1.5e-8}).x
        dinv = 1 / ((1 - s) * d + s)
    dinv[~np.isfinite(dinv)] = 0
    P = (v * dinv) @ v.T
    pd_ = np.diag(P)
    Pz = P @ z
    return pd.DataFrame({"z": z, "expected": z - Pz / pd_, "sd": np.sqrt(1 / pd_),
                         "z_std_diff": Pz / np.sqrt(pd_), "s": s})
