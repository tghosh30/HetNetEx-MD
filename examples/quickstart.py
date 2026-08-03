"""
Quickstart: exact inference for the Multi-DWPC gene-set resampling null.

Runs in a couple of seconds on a toy problem and prints the exact null moments,
the effect size, the analytic p-value, and -- for comparison -- what B=1000
resampling gives for the same input.

    python examples/quickstart.py
"""
import time

import numpy as np

from hetnetex_md import exact_median_pvalue, hetnetex_md_resampling

rng = np.random.default_rng(0)

# ---------------------------------------------------------------------------
# A per-gene score vector, already computed from the observed graph.
# In production these are arcsinh-transformed DWPC values to a fixed target;
# here we just simulate something with the right shape.
# ---------------------------------------------------------------------------
N = 5_000                                   # genes in the universe
scores = np.abs(rng.gamma(0.6, 1.0, size=N))
scores[rng.random(N) < 0.55] = 0.0          # most genes have no path

# Degree strata: 10 quantile bins of total degree.
degree = rng.gamma(2.0, 5.0, size=N)
order = np.argsort(np.argsort(degree))
bin_ids = np.minimum((order / N * 10).astype(int), 9)

# The query gene set, with implanted signal.
K = 150
gene_idx = rng.choice(N, size=K, replace=False)
scores[gene_idx[:60]] += rng.gamma(1.5, 0.8, size=60)

# ---------------------------------------------------------------------------
# Exact inference -- one pass, no replicates.
# ---------------------------------------------------------------------------
t0 = time.perf_counter()
res = hetnetex_md_resampling(scores, gene_idx, bin_ids)
t_exact = time.perf_counter() - t0

print("HetNetEX-MD, exact (Theorem 2 + Edgeworth)")
print(f"  observed T      {res['t_obs']:.6f}")
print(f"  null mean       {res['null_mean']:.6f}")
print(f"  null SD         {res['null_sd']:.6f}")
print(f"  null skewness   {res['null_skew']:+.4f}   (exact, not estimated)")
print(f"  z               {res['z']:.3f}")
print(f"  p (Edgeworth)   {res['p_edgeworth']:.3e}")
print(f"  p (normal)      {res['p_normal']:.3e}")
print(f"  time            {t_exact*1000:.2f} ms")

# ---------------------------------------------------------------------------
# The same quantity by resampling, for comparison.
# ---------------------------------------------------------------------------
B = 1000
pools, counts = [], []
for b in np.unique(bin_ids[gene_idx]):
    k = int(np.sum(bin_ids[gene_idx] == b))
    universe = np.ones(N, bool)
    universe[gene_idx] = False
    pools.append(np.flatnonzero((bin_ids == b) & universe))
    counts.append(k)

t0 = time.perf_counter()
draws = np.empty(B)
for i in range(B):
    total = sum(scores[rng.choice(p, size=k, replace=False)].sum()
                for p, k in zip(pools, counts))
    draws[i] = total / K
p_emp = (1 + int(np.sum(draws >= res["t_obs"]))) / (1 + B)
t_mc = time.perf_counter() - t0

print(f"\nResampling, B={B}")
print(f"  null mean       {draws.mean():.6f}   "
      f"(exact: {res['null_mean']:.6f})")
print(f"  null SD         {draws.std(ddof=1):.6f}   "
      f"(exact: {res['null_sd']:.6f})")
print(f"  p (empirical)   {p_emp:.3e}   floor = {1/(1+B):.3e}")
print(f"  time            {t_mc*1000:.2f} ms   -> {t_mc/t_exact:.0f}x slower")

print("\nNote: the p-values differ by construction -- one counts exceedances")
print("among finitely many draws, the other evaluates a limiting tail.")
print("They agree in ranking, which is what a discovery procedure needs.")

# ---------------------------------------------------------------------------
# The median aggregate also has an exact null -- and needs no limit theorem.
# ---------------------------------------------------------------------------
m_obs = float(np.median(scores[gene_idx]))
p_med = exact_median_pvalue(scores, pools, counts, m_obs)
print(f"\nMedian aggregate (Theorem 5b, distribution-free)")
print(f"  observed median {m_obs:.6f}")
print(f"  p (exact)       {p_med:.3e}")
if m_obs == 0.0:
    print("  (median is exactly zero -- the test is vacuous here, which is")
    print("   common on sparse metapaths; prefer the mean aggregate)")
