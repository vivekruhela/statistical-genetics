"""SuSiE from summary statistics (SuSiE-RSS), a compact Python port of susieR::susie_rss.

The R track calls susieR directly. This file is the Python twin: it follows susieR 0.14's
susie_rss(z, R, n) with its defaults (L = 10, prior variance estimated by Brent's method on the
log scale, residual variance fixed at 1, 95% credible sets kept only if their purity, the
smallest |r| between members, is at least 0.5), so the two give the same PIPs and credible
sets on the same input. tests/test_lib.R and tests/test_lib.py check that they do.

    fit = susie_rss(z, R, n)
    fit["pip"], fit["sets"]  # PIP per variant; list of credible sets (0-based indices)

Wang et al. 2020 (JRSS-B) describe SuSiE; Zou et al. 2022 (PLoS Genetics) describe SuSiE-RSS.
"""
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import norm

_EPS = np.sqrt(np.finfo(float).eps)


class SusieError(RuntimeError):
    """Raised where susieR stops with an error (for example, z and R that disagree badly)."""


def _loglik(V, betahat, shat2, prior_weights):
    lbf = norm.logpdf(betahat, 0, np.sqrt(V + shat2)) - norm.logpdf(betahat, 0, np.sqrt(shat2))
    lpo = lbf + np.log(prior_weights + _EPS)
    m = lpo.max()
    return np.log(np.exp(lpo - m).sum()) + m


def _optimize_prior_variance(betahat, shat2, prior_weights, V_init):
    f = lambda lv: -_loglik(np.exp(lv), betahat, shat2, prior_weights)
    lv = minimize_scalar(f, bounds=(-30, 15), method="bounded", options={"xatol": 1.5e-8}).x
    with np.errstate(divide="ignore"):
        lv_init = np.log(V_init)
    if f(lv) > f(lv_init):
        lv = lv_init
    V = np.exp(lv)
    if _loglik(0, betahat, shat2, prior_weights) >= _loglik(V, betahat, shat2, prior_weights):
        V = 0.0
    return V


def _single_effect(XtR, d, V, sigma2, prior_weights):
    betahat = XtR / d
    shat2 = sigma2 / d
    V = _optimize_prior_variance(betahat, shat2, prior_weights, V)
    lbf = norm.logpdf(betahat, 0, np.sqrt(V + shat2)) - norm.logpdf(betahat, 0, np.sqrt(shat2))
    lpo = lbf + np.log(prior_weights + _EPS)
    m = lpo.max()
    w = np.exp(lpo - m)
    alpha = w / w.sum()
    post_var = 1 / (1 / V + d / sigma2) if V > 0 else np.zeros_like(d)
    post_mean = post_var * XtR / sigma2
    return alpha, post_mean, post_var + post_mean ** 2, V, m + np.log(w.sum())


def _credible_sets(alpha, V, R, coverage=0.95, min_abs_corr=0.5):
    sets, seen = [], []
    for l in range(alpha.shape[0]):
        order = np.argsort(-alpha[l], kind="stable")
        k = int(np.searchsorted(np.cumsum(alpha[l, order]), coverage)) + 1
        cs = sorted(order[:k].tolist())
        duplicate = cs in seen
        seen.append(cs)
        if V[l] <= 1e-9 or duplicate:
            continue
        sub = np.abs(R[np.ix_(cs, cs)])
        purity = 1.0 if len(cs) == 1 else sub[np.triu_indices(len(cs), 1)].min()
        if purity >= min_abs_corr:
            sets.append({"effect": l, "variants": cs, "purity": float(purity),
                         "coverage": float(alpha[l, cs].sum())})
    return sorted(sets, key=lambda s: -s["purity"])


def susie_rss(z, R, n, L=10, max_iter=100, tol=1e-3, coverage=0.95, min_abs_corr=0.5,
              scaled_prior_variance=0.2, check_prior=True):
    """Fit SuSiE to z-scores z with LD (correlation) matrix R from n people.

    Returns a dict with pip, alpha (L x p), V (prior variance per effect), sets (credible sets,
    0-based indices, purest first), elbo, niter and converged. Raises SusieError where susieR
    would stop.
    """
    z = np.nan_to_num(np.asarray(z, dtype=float))
    R = np.asarray(R, dtype=float)
    p = len(z)
    adj = (n - 1) / (z ** 2 + n - 2)
    z = np.sqrt(adj) * z
    XtX = (n - 1) * R
    XtX = (XtX + XtX.T) / 2
    Xty = np.sqrt(n - 1) * z
    yty = n - 1
    d = np.diag(XtX).copy()
    sigma2 = 1.0
    L = min(L, p)
    pi = np.full(p, 1 / p)
    alpha = np.full((L, p), 1 / p)
    mu = np.zeros((L, p))
    mu2 = np.zeros((L, p))
    V = np.full(L, scaled_prior_variance)
    KL = np.zeros(L)
    XtXr = np.zeros(p)
    zm = np.max(np.abs(Xty / d / np.sqrt(sigma2 / d)))
    elbo = [-np.inf]
    converged = False
    for it in range(max_iter):
        for l in range(L):
            XtXr -= XtX @ (alpha[l] * mu[l])
            XtR = Xty - XtXr
            alpha[l], mu[l], mu2[l], V[l], lbf_model = _single_effect(XtR, d, V[l], sigma2, pi)
            Eb, Eb2 = alpha[l] * mu[l], alpha[l] * mu2[l]
            KL[l] = -lbf_model - 0.5 / sigma2 * (-2 * Eb @ XtR + d @ Eb2)
            XtXr += XtX @ (alpha[l] * mu[l])
        if check_prior and np.any(V > 100 * zm ** 2):
            raise SusieError("the estimated prior variance is unreasonably large: "
                             "the z-scores and the LD matrix disagree")
        B = alpha * mu
        betabar = B.sum(0)
        ER2 = (yty - 2 * betabar @ Xty + betabar @ XtX @ betabar
               - np.sum((B @ XtX) * B) + np.sum(d * (alpha * mu2).sum(0)))
        e = -n / 2 * np.log(2 * np.pi * sigma2) - ER2 / (2 * sigma2) - KL.sum()
        if not np.isfinite(e):
            raise SusieError("the objective became infinite")
        elbo.append(e)
        if elbo[-1] - elbo[-2] < tol:
            converged = True
            break
    keep = V > 1e-9
    pip = 1 - np.prod(1 - alpha[keep], axis=0) if keep.any() else np.zeros(p)
    return {"pip": pip, "alpha": alpha, "V": V, "sets": _credible_sets(alpha, V, R, coverage, min_abs_corr),
            "elbo": np.array(elbo[1:]), "niter": it + 1, "converged": converged}
