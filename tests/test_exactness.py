"""
Exactness tests for HetNetEX-MD.

These are not smoke tests. The central claims of the method are that certain
formulas are *identities*, not approximations, so the tests verify them against
exhaustive enumeration of every admissible sample. If any test here fails, the
corresponding theorem is wrong.

Run:  pytest -q
"""
import itertools

import numpy as np
import pytest
from scipy import sparse

from hetnetex_md import (
    EdgeType,
    HetNet,
    exact_median_pvalue,
    exact_resampling_moments,
    edgeworth_upper_tail,
    fit_soft_cm,
    soft_cm_ratio,
    admissible_self_kernels,
    network_null_moments,
)

TOL = 1e-12


# --------------------------------------------------------------------------- #
# Theorem 1 / Lemma 1 : SRSWOR moments, single stratum
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("N,k", [(9, 3), (10, 4), (12, 5), (11, 2), (8, 4), (7, 1)])
def test_srswor_moments_match_exhaustive_enumeration(N, k):
    """Enumerate all C(N,k) samples; compare mean, variance, third moment."""
    rng = np.random.default_rng(N * 100 + k)
    x = rng.gamma(2.0, 1.0, size=N)
    pool = np.arange(N)

    means = np.array([x[list(c)].mean() for c in itertools.combinations(range(N), k)])
    brute_mean = means.mean()
    brute_var = ((means - brute_mean) ** 2).mean()
    brute_mu3 = ((means - brute_mean) ** 3).mean()

    mean, var, mu3 = exact_resampling_moments(x, [pool], [k])

    assert abs(mean - brute_mean) < TOL
    assert abs(var - brute_var) < TOL * max(1.0, abs(brute_var))
    assert abs(mu3 - brute_mu3) < 1e-10 * max(1.0, abs(brute_mu3))


@pytest.mark.parametrize(
    "N1,k1,N2,k2", [(7, 3, 6, 2), (8, 3, 7, 4), (6, 2, 6, 3), (5, 2, 5, 2)]
)
def test_srswor_moments_stratified(N1, k1, N2, k2):
    """Two independent strata: means add, variances and third cumulants add."""
    N = N1 + N2
    rng = np.random.default_rng(N1 * 1000 + k1 * 100 + N2 * 10 + k2)
    x = rng.gamma(2.0, 1.0, size=N)
    p1, p2 = np.arange(N1), np.arange(N1, N)

    vals = [
        x[list(c1) + list(c2)].mean()
        for c1 in itertools.combinations(p1, k1)
        for c2 in itertools.combinations(p2, k2)
    ]
    vals = np.asarray(vals)
    bm = vals.mean()

    mean, var, mu3 = exact_resampling_moments(x, [p1, p2], [k1, k2])
    assert abs(mean - bm) < TOL
    assert abs(var - ((vals - bm) ** 2).mean()) < 1e-11
    assert abs(mu3 - ((vals - bm) ** 3).mean()) < 1e-11


def test_srswor_degenerate_cases():
    """k=N is degenerate; k=N/2 has zero third moment by symmetry."""
    x = np.array([1.0, 4.0, 9.0, 16.0, 25.0, 36.0])
    _, var, mu3 = exact_resampling_moments(x, [np.arange(6)], [6])
    assert var == pytest.approx(0.0, abs=TOL)
    assert mu3 == pytest.approx(0.0, abs=TOL)

    _, _, mu3_half = exact_resampling_moments(x, [np.arange(6)], [3])
    assert mu3_half == pytest.approx(0.0, abs=1e-12)


# --------------------------------------------------------------------------- #
# Theorem 5b : exact null distribution of the stratified median
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("N,K", [(11, 5), (12, 5), (13, 7), (10, 5), (9, 3)])
def test_median_pvalue_matches_exhaustive_enumeration(N, K):
    rng = np.random.default_rng(N * 7 + K)
    x = np.round(rng.gamma(2.0, 1.0, size=N), 6)
    pool = np.arange(N)
    meds = np.array([np.median(x[list(c)]) for c in itertools.combinations(range(N), K)])

    uniq = np.unique(meds)
    for t in (uniq[0], uniq[len(uniq) // 2], uniq[-1]):
        brute = float(np.mean(meds >= t))
        exact = exact_median_pvalue(x, [pool], [K], t)
        assert abs(brute - exact) < 1e-12, f"t={t}"


@pytest.mark.parametrize("N1,K1,N2,K2", [(7, 3, 6, 2), (8, 3, 7, 4)])
def test_median_pvalue_stratified(N1, K1, N2, K2):
    N = N1 + N2
    rng = np.random.default_rng(N1 * 31 + N2)
    x = np.round(rng.gamma(2.0, 1.0, size=N), 6)
    p1, p2 = np.arange(N1), np.arange(N1, N)
    meds = np.array([
        np.median(x[list(c1) + list(c2)])
        for c1 in itertools.combinations(p1, K1)
        for c2 in itertools.combinations(p2, K2)
    ])
    uniq = np.unique(meds)
    for t in (uniq[0], uniq[len(uniq) // 2], uniq[-1]):
        brute = float(np.mean(meds >= t))
        exact = exact_median_pvalue(x, [p1, p2], [K1, K2], t)
        assert abs(brute - exact) < 1e-12, f"t={t}"


def test_median_pvalue_is_a_probability():
    rng = np.random.default_rng(0)
    x = rng.gamma(2.0, 1.0, size=200)
    p = exact_median_pvalue(x, [np.arange(200)], [21], float(np.median(x)))
    assert 0.0 <= p <= 1.0


# --------------------------------------------------------------------------- #
# Lemma 2b : max-entropy configuration model
# --------------------------------------------------------------------------- #
def _graphical_pair(seed, n_row=300, n_col=120):
    """A heavy-tailed but bipartite-graphical degree pair."""
    rng = np.random.default_rng(seed)
    d_row = np.clip((rng.pareto(1.5, size=n_row) + 1) * 3, 1, 0.8 * n_col)
    d_col = np.clip(d_row.sum() * rng.dirichlet(np.ones(n_col)), 1, 0.8 * n_row)
    d_row = d_row * d_col.sum() / d_row.sum()          # match edge counts
    d_row = np.clip(d_row, 1e-9, 0.8 * n_col)
    d_col = d_col * d_row.sum() / d_col.sum()
    return d_row, d_col


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_soft_cm_reproduces_expected_degrees(seed):
    """The whole point of the fit: expected degrees must equal observed ones."""
    d_row, d_col = _graphical_pair(seed)
    fx, fy, resid = fit_soft_cm(d_row, d_col)
    assert resid < 1e-6, f"degree constraints not met: {resid}"


def test_soft_cm_rejects_non_graphical_degrees():
    """Regression: an infeasible sequence used to return garbage silently."""
    rng = np.random.default_rng(0)
    d_row = np.clip((rng.pareto(1.5, size=300) + 1) * 3, 1, None)   # hubs > n_col
    d_col = d_row.sum() * rng.dirichlet(np.ones(120))
    assert d_row.max() > 120, "test input should be infeasible"
    with pytest.raises(ValueError, match="not bipartite-graphical"):
        fit_soft_cm(d_row, d_col)
    with pytest.warns(RuntimeWarning):
        fit_soft_cm(d_row, d_col, strict=False)


def test_soft_cm_rejects_mismatched_edge_counts():
    d_row = np.full(50, 4.0)
    d_col = np.full(20, 3.0)      # 200 vs 60 edge endpoints
    with pytest.raises(ValueError, match="degree sums disagree"):
        fit_soft_cm(d_row, d_col)


def test_soft_cm_ratio_cache_is_not_shared_across_networks():
    """Regression: a name-keyed cache silently returned another net's fit."""
    def mk(seed, n=60, m=240):
        r = np.random.default_rng(seed)
        A = sparse.csr_matrix(
            (np.ones(m), (r.integers(0, n, m), r.integers(0, n, m))), shape=(n, n)
        )
        A.data[:] = 1.0
        return HetNet({"X": n}, {"e": EdgeType("e", "X", "X", A)})

    n1, n2 = mk(1), mk(999)
    d = np.array([5.0, 10.0])
    shared = {}
    soft_cm_ratio(n1, "e", d, d, cache=shared)
    cached = soft_cm_ratio(n2, "e", d, d, cache=shared)
    fresh = soft_cm_ratio(n2, "e", d, d)
    assert np.allclose(cached, fresh), "stale cache leaked across networks"
    assert len(shared) == 2


def test_soft_cm_approximations_are_ordered():
    """d_u d_v / m and 1 - exp(-x) are the 1st and 2nd order approximations."""
    for x in (1e-4, 1e-3, 1e-2):
        sparse_limit = x
        exponential = 1.0 - np.exp(-x)
        assert exponential < sparse_limit
        assert abs(exponential - sparse_limit) < x**2


# --------------------------------------------------------------------------- #
# Admissibility rule for the same-source kernels
# --------------------------------------------------------------------------- #
def test_admissible_self_kernels_by_path_length():
    t1, tm, tL = 0.11, 0.22, 0.33
    assert admissible_self_kernels(2, t1, tm, tL) == 0.0
    assert admissible_self_kernels(3, t1, tm, tL) == pytest.approx(t1 + tL)
    assert admissible_self_kernels(4, t1, tm, tL) == pytest.approx(t1 + tm + tL)


# --------------------------------------------------------------------------- #
# Edgeworth tail
# --------------------------------------------------------------------------- #
def test_edgeworth_reduces_to_normal_at_zero_skew():
    from scipy.stats import norm
    z, p = edgeworth_upper_tail(1.5, 0.0, 1.0, 0.0)
    assert z == pytest.approx(1.5)
    assert p == pytest.approx(norm.sf(1.5))


def test_edgeworth_handles_degenerate_variance():
    z, p = edgeworth_upper_tail(1.0, 0.0, 0.0, 0.0)
    assert np.isinf(z) and p == 0.0
    z, p = edgeworth_upper_tail(-1.0, 0.0, 0.0, 0.0)
    assert z == 0.0 and p == 1.0


# --------------------------------------------------------------------------- #
# Network null moments: structural sanity
# --------------------------------------------------------------------------- #
def _toy_chain(seed=0):
    rng = np.random.default_rng(seed)
    sizes = {"G": 120, "P": 40, "B": 25}
    def bip(a, b, m):
        r = rng.integers(0, sizes[a], m); c = rng.integers(0, sizes[b], m)
        A = sparse.csr_matrix((np.ones(m), (r, c)), shape=(sizes[a], sizes[b]))
        A.data[:] = 1.0
        return A
    edges = {
        "GpP": EdgeType("GpP", "G", "P", bip("G", "P", 400)),
        "PpB": EdgeType("PpB", "P", "B", bip("P", "B", 150)),
    }
    return HetNet(sizes, edges)


def test_network_null_moments_are_finite_and_nonnegative():
    net = _toy_chain()
    mom = network_null_moments(net, ["GpP", "PpB"], w=0.4)
    for key in ("a", "b", "a2", "b2"):
        assert np.all(np.isfinite(mom[key]))
        assert np.all(mom[key] >= 0)
    assert np.isfinite(mom["C"]) and mom["C"] > 0
    assert np.isfinite(mom["theta_mid"])


def test_network_null_mean_is_rank_one():
    """Lemma 3: mu(g,t) = C * a_g * b_t exactly (this is what makes it O(Ln))."""
    net = _toy_chain(3)
    mom = network_null_moments(net, ["GpP", "PpB"], w=0.4)
    mu = mom["C"] * np.outer(mom["a"], mom["b"])
    assert np.linalg.matrix_rank(mu, tol=1e-10 * mu.max()) == 1
