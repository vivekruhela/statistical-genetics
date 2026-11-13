"""Unit tests for lib/statgen.py (pytest), mirroring tests/test_lib.R."""
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

sys.path.insert(0, "lib")
import statgen as sg  # noqa: E402


def test_inverse_normal():
    x = np.array([5, 1, 3, 3, 10.0])
    z = sg.inverse_normal(x)
    assert z[2] == z[3] and np.argmax(z) == 4
    assert abs(sg.inverse_normal(np.arange(7.0)).sum()) < 1e-12


def test_size_factors():
    m = pd.DataFrame(np.array([[10, 20, 30], [20, 40, 60], [5, 10, 15]]).T)
    expect = np.array([1, 2, 0.5]) / np.exp(np.mean(np.log([1, 2, 0.5])))
    assert np.allclose(sg.size_factors(m), expect)


def test_cis_scan_matches_ols():
    rng = np.random.default_rng(1)
    n = 80
    g = rng.binomial(2, 0.4, (5, n)).astype(float)
    cv = rng.normal(size=(n, 2))
    e = (0.5 * g[2] + cv[:, 0] + rng.normal(size=n))[None, :]
    res = sg.cis_scan(sg.residualise(e, cv), ["gene"], [300], sg.residualise(g, cv),
                      [f"v{i}" for i in range(1, 6)], [100, 200, 300, 400, 500], 1000, 2)
    fit = sm.OLS(e[0], sm.add_constant(np.column_stack([g[2], cv]))).fit()
    assert res.lead_id[0] == "v3"
    assert np.isclose(res.lead_mlog10p[0], -np.log10(fit.pvalues[1]), rtol=1e-8)


def test_linear_assoc_matches_ols():
    rng = np.random.default_rng(2)
    n = 200
    G = rng.binomial(2, 0.3, (3, n)).astype(float)
    X = np.column_stack([np.ones(n), rng.normal(size=n), rng.binomial(1, 0.5, n)])
    y = 0.3 * G[1] + X[:, 1] + rng.normal(size=n)
    st = sg.assoc_stats(G, y, X)
    fit = sm.OLS(y, np.column_stack([X, G[1]])).fit()
    assert np.isclose(st.beta[1], fit.params[-1]) and np.isclose(st.se[1], fit.bse[-1])


def test_logistic_score_matches_statsmodels():
    rng = np.random.default_rng(3)
    n = 500
    G = rng.binomial(2, 0.3, (1, n)).astype(float)
    X = np.column_stack([np.ones(n), rng.normal(size=n)])
    y = rng.binomial(1, 1 / (1 + np.exp(-(-0.5 + 0.4 * G[0] + 0.3 * X[:, 1])))).astype(float)
    st = sg.assoc_stats(G, y, X, binary=True)
    full = sm.Logit(y, np.column_stack([X, G[0]]))
    score = full.score(np.r_[sm.Logit(y, X).fit(disp=0).params, 0.0])
    info = full.hessian(np.r_[sm.Logit(y, X).fit(disp=0).params, 0.0])
    rao = score @ np.linalg.solve(-info, score)
    assert np.isclose(st.z[0] ** 2, rao, rtol=1e-6)


def test_mosaic_shape():
    rng = np.random.default_rng(4)
    H = rng.binomial(1, 0.5, (50, 20)).astype(np.int8)
    d = sg.mosaic_dosage(H, np.arange(1000, 50001, 1000), 30, 10000, 0, rng)
    assert d.shape == (50, 30) and set(np.unique(d)) <= {0, 1, 2}


def test_dag_unit_variance():
    rng = np.random.default_rng(5)
    nodes = pd.DataFrame({"node": ["a", "b"], "order": [1, 2], "resid_sd": [1, np.sqrt(0.7)]})
    edges = pd.DataFrame({"from": ["a"], "to": ["b"], "weight": [np.sqrt(0.3)], "center": [0.0]})
    v = sg.simulate_dag(nodes, edges, {}, 20000, rng)
    assert abs(v["b"].var() - 1) < 0.03 and abs(stats.pearsonr(v["a"], v["b"])[0] ** 2 - 0.3) < 0.03


def test_read_geuvadis_counts(tmp_path):
    # GEUVADIS layout: versioned TargetID, three info columns, one column per library;
    # HG1 has two libraries and the deeper one is kept
    p = tmp_path / "counts.txt"
    pd.DataFrame({"TargetID": ["ENSG1.5", "ENSG2.1"], "Gene_Symbol": ["ENSG1.5", "ENSG2.1"],
                  "Chr": ["1", "2"], "Coord": [10, 20],
                  "HG1.1.M_1": [1, 2], "HG1.2.M_2": [10, 20], "HG2.1.M_3": [3, 4]}).to_csv(p, sep="\t", index=False)
    m = sg.read_geuvadis_counts(p)
    assert m.index.name == "gene_id" and list(m.index) == ["ENSG1", "ENSG2"]
    assert list(m.columns) == ["HG1", "HG2"] and m.loc["ENSG1", "HG1"] == 10


# ---- LD and reference panels (post R) ----
def toy_ld():
    p, n = 30, 1000
    i = np.arange(p)
    R = 0.9 ** np.abs(np.subtract.outer(i, i))
    b = np.zeros(p)
    b[7], b[21] = 0.15, -0.12
    return R, n, np.sqrt(n) * R @ b + 0.8 * np.sin(np.arange(1, p + 1))


def test_ld_matrix_and_round_trip(tmp_path):
    g = np.array([[0, 1, 2, 1, 0], [0, 1, 2, 2, 0], [1, 1, 1, 1, 1]], dtype=float)
    r = sg.ld_matrix(g)
    assert abs(r[0, 1] - np.corrcoef(g[0], g[1])[0, 1]) < 1e-12
    assert np.isnan(r[2]).all()
    R, _, _ = toy_ld()
    sg.write_ld(R, [f"v{k}" for k in range(30)], str(tmp_path / "toy"))
    R2, ids = sg.read_ld(str(tmp_path / "toy"))
    assert np.allclose(R2, R, atol=1e-7) and ids[3] == "v3"


def test_susie_rss_matches_susieR():
    from susie_rss import susie_rss
    R, n, z = toy_ld()
    fit = susie_rss(z, R, n)
    # susieR 0.14.2: susie_rss(z, R, n = 1000)$pip on the same toy data
    expect = [0.000200, 0.000257, 0.000234, 0.000456, 0.000674, 0.000657, 0.029513, 0.958001, 0.012094, 0.000954,
              0.001967, 0.001140, 0.000451, 0.000361, 0.000480, 0.002573, 0.018182, 0.021978, 0.006018, 0.002763,
              0.007286, 0.193364, 0.547232, 0.176512, 0.009185, 0.000940, 0.000483, 0.000797, 0.002323, 0.002740]
    assert np.allclose(fit["pip"], expect, atol=2e-6)
    assert [s["variants"] for s in fit["sets"]] == [[7]]


def test_ld_check_on_matching_ld():
    R, n, z = toy_ld()
    chk = sg.ld_check(z, R, n)
    assert chk.s[0] < 1e-6        # the z-scores were made from this R, so it fits
    assert chk.z_std_diff.abs().max() < 4
