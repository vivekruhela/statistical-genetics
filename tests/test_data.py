"""Checks that the built cohort carries the planted truth (pytest).

Run after the data build:  STATGEN_DATA=data pytest tests/test_data.py
Each check uses only simple marginal statistics, so a failure points at the build, not at a
fancy method. Thresholds are loose enough for the 462-person cohort.
"""
import os
import sys

import numpy as np
import pandas as pd
import pytest
from scipy import stats

sys.path.insert(0, "lib")
import statgen as sg  # noqa: E402

DATA = os.environ.get("STATGEN_DATA", "data")
pytestmark = pytest.mark.skipif(not os.path.exists(f"{DATA}/truth/truth.tsv"), reason="data not built")


def table(path, **kw):
    return pd.read_csv(f"{DATA}/{path}", sep="\t", **kw)


@pytest.fixture(scope="module")
def cohort():
    samples = table("cohort/samples.tsv")["sample"].to_list()
    truth = table("truth/truth.tsv")
    geno = {}
    for locus in truth.locus.unique():
        info, g, ids = sg.read_traw(f"{DATA}/cohort/genotypes/loci/{locus}.traw")
        g = g[:, [ids.index(s) for s in samples]]
        geno.update({v: g[i] for i, v in enumerate(info["id"]) if v in set(truth.variant_id)})
    return {
        "samples": samples, "truth": truth, "geno": geno,
        "expr": table("cohort/expression/expression_int.tsv.gz", index_col=0)[samples],
        "meth": table("cohort/methylation/methylation_m.tsv.gz", index_col=0)[samples],
        "clin": table("cohort/clinical/clinical.tsv").set_index("sample").loc[samples],
    }


def assoc(x, y):
    r, p = stats.pearsonr(x, y)
    return r, p


def test_eqtls_are_planted(cohort):
    t = cohort["truth"]
    for row in t[t.effect_type == "eqtl"].itertuples():
        gene = row.acts_on.split(":")[1]
        r, p = assoc(cohort["geno"][row.variant_id], cohort["expr"].loc[gene])
        assert p < 1e-4, f"{row.locus} {gene}: eQTL at {row.variant_id} p={p:.2g}"
        assert np.sign(r) == np.sign(row.per_allele_effect)


def test_mqtls_are_planted(cohort):
    t = cohort["truth"]
    for row in t[t.effect_type == "mqtl"].itertuples():
        cpg = row.acts_on.split(":")[1]
        r, p = assoc(cohort["geno"][row.variant_id], cohort["meth"].loc[cpg])
        assert p < 1e-10 and np.sign(r) == np.sign(row.per_allele_effect), f"{row.locus} mQTL p={p:.2g}"


def test_methylation_lowers_expression(cohort):
    t = cohort["truth"]
    for row in t[t.effect_type == "mqtl"].itertuples():
        cpg = row.acts_on.split(":")[1]
        gene = row.path.split(" -> ")[1].split(":")[1]
        r, p = assoc(cohort["meth"].loc[cpg], cohort["expr"].loc[gene])
        assert r < 0 and p < 1e-6, f"{row.locus}: CpG-expression r={r:.2f}"


def test_strong_trait_locus_is_visible(cohort):
    t = cohort["truth"]
    strong = t[(t.scenario == "gwas_only") & (t.effect_type == "direct")].iloc[0]
    r, p = assoc(cohort["geno"][strong.variant_id], cohort["clin"]["trait_q"])
    assert p < 1e-5 and np.sign(r) == np.sign(strong.per_allele_effect)


def test_gwas_only_gene_has_no_eqtl(cohort):
    t = cohort["truth"]
    row = t[t.scenario == "gwas_only"].iloc[0]
    gene = "SIMG%02d" % int(row.locus[1:])
    _, p = assoc(cohort["geno"][row.variant_id], cohort["expr"].loc[gene])
    assert p > 1e-4


def test_disease_prevalence(cohort):
    prev = cohort["clin"]["disease"].mean()
    assert 0.2 < prev < 0.4


def test_de_genes_move_the_right_way(cohort):
    de = table("truth/de_genes.tsv")
    logn = table("cohort/expression/expression_log2norm.tsv.gz", index_col=0)[cohort["samples"]]
    assert de.gene_id.isin(logn.index).mean() > 0.95, "DE genes lost in gene filtering"
    de = de[de.gene_id.isin(logn.index)]
    case = cohort["clin"]["disease"].to_numpy() == 1
    diff = logn.loc[de.gene_id].to_numpy()[:, case].mean(1) - logn.loc[de.gene_id].to_numpy()[:, ~case].mean(1)
    agree = np.mean(np.sign(diff) == np.sign(de.log2fc.to_numpy()))
    assert agree > 0.9, f"only {agree:.0%} of DE genes move in the planted direction"


def test_ancestry_variant_frequencies(cohort):
    t = cohort["truth"]
    row = t[t.scenario == "ancestry"].iloc[0]
    maf = lambda p: min(p, 1 - p)  # noqa: E731
    assert maf(row.af_eur) < 0.15 and maf(row.af_yri) > 0.25


def test_no_missing_values(cohort):
    for k in ("expr", "meth"):
        assert not cohort[k].isna().any().any()
    assert not cohort["clin"][["age", "bmi", "trait_q", "disease"]].isna().any().any()


@pytest.mark.parametrize("pop", ["EUR", "AFR"])
def test_gwas_finds_trait_variants(cohort, pop):
    path = f"{DATA}/gwas/{pop}/trait_q.sumstats.tsv.gz"
    if not os.path.exists(path):
        pytest.skip("GWAS cohort not built")
    ss = pd.read_csv(path, sep="\t").set_index("variant_id")
    t = cohort["truth"]
    hits = t[t.endpoint == "trait_q"]
    z = ss.loc[hits.variant_id, "z"].to_numpy()
    agree = np.sign(z) == np.sign(hits.endpoint_effect.to_numpy())
    if pop == "EUR":  # the main GWAS: every trait variant clearly found, right direction
        # (except the ancestry locus, whose variant is rare in Europeans by design)
        bad = hits[~(agree & (np.abs(z) > 4)) & (hits.scenario != "ancestry").to_numpy()]
        assert bad.empty, f"missed: {list(zip(bad.locus, bad.variant_id))}"
    else:  # smaller, different LD and frequencies: most directions still right
        assert agree.mean() >= 0.75


def test_gwas_disease_locus(cohort):
    path = f"{DATA}/gwas/EUR/disease.sumstats.tsv.gz"
    if not os.path.exists(path):
        pytest.skip("GWAS cohort not built")
    ss = pd.read_csv(path, sep="\t").set_index("variant_id")
    t = cohort["truth"]
    row = t[t.scenario == "binary"].iloc[0]
    z = ss.loc[row.variant_id, "z"]
    assert abs(z) > 2.5 and np.sign(z) == np.sign(row.endpoint_effect)


def test_eqtl_only_locus_is_quiet_in_gwas(cohort):
    path = f"{DATA}/gwas/EUR/trait_q.sumstats.tsv.gz"
    if not os.path.exists(path):
        pytest.skip("GWAS cohort not built")
    ss = pd.read_csv(path, sep="\t")
    lid = cohort["truth"].query("scenario == 'eqtl_only'").locus.iloc[0]
    # null locus: with thousands of correlated tests the best p-value should stay modest
    assert ss[ss.locus == lid].mlog10p.max() < 7
