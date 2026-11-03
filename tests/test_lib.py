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
