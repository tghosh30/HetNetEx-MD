"""Vignette 1: the exact resampling null, verified against brute force."""
import itertools
import time

import numpy as np

from hetnetex_md import exact_resampling_moments, hetnetex_md_resampling

print("## 1. The formulas are identities\n")
rng = np.random.default_rng(11)
N, k = 12, 5
x = rng.gamma(2.0, 1.0, size=N)
means = np.array([x[list(c)].mean() for c in itertools.combinations(range(N), k)])
bm = means.mean()
brute = (bm, ((means - bm) ** 2).mean(), ((means - bm) ** 3).mean())
exact = exact_resampling_moments(x, [np.arange(N)], [k])
print(f"Enumerated all C({N},{k}) = {len(means)} samples.\n")
print(f"{'':14s}{'brute force':>22s}{'closed form':>22s}{'difference':>14s}")
for nm, b, e in zip(("mean", "variance", "third moment"), brute, exact):
    print(f"{nm:14s}{b:22.15f}{e:22.15f}{abs(b - e):14.2e}")

print("\n\n## 2. Stratification\n")
N1, k1, N2, k2 = 7, 3, 6, 2
y = rng.gamma(2.0, 1.0, size=N1 + N2)
p1, p2 = np.arange(N1), np.arange(N1, N1 + N2)
vals = np.array([
    y[list(c1) + list(c2)].mean()
    for c1 in itertools.combinations(p1, k1)
    for c2 in itertools.combinations(p2, k2)
])
vm = vals.mean()
e = exact_resampling_moments(y, [p1, p2], [k1, k2])
print(f"Two strata: C({N1},{k1}) x C({N2},{k2}) = {len(vals)} samples.\n")
print(f"{'':14s}{'brute force':>22s}{'closed form':>22s}{'difference':>14s}")
for nm, b, ee in zip(
    ("mean", "variance", "third moment"),
    (vm, ((vals - vm) ** 2).mean(), ((vals - vm) ** 3).mean()),
    e,
):
    print(f"{nm:14s}{b:22.15f}{ee:22.15f}{abs(b - ee):14.2e}")

print("\n\n## 3. At realistic scale\n")
rng = np.random.default_rng(0)
N = 21_000
scores = np.abs(rng.gamma(0.6, 1.0, size=N))
scores[rng.random(N) < 0.55] = 0.0
deg = rng.gamma(2.0, 5.0, size=N)
bins = np.minimum((np.argsort(np.argsort(deg)) / N * 10).astype(int), 9)
K = 200
gene_idx = rng.choice(N, size=K, replace=False)
scores[gene_idx[:80]] += rng.gamma(1.5, 0.8, size=80)

t0 = time.perf_counter()
res = hetnetex_md_resampling(scores, gene_idx, bins)
t_ex = time.perf_counter() - t0

B = 1000
pools, counts = [], []
for b in np.unique(bins[gene_idx]):
    kk = int(np.sum(bins[gene_idx] == b))
    uni = np.ones(N, bool)
    uni[gene_idx] = False
    pools.append(np.flatnonzero((bins == b) & uni))
    counts.append(kk)

t0 = time.perf_counter()
draws = np.array([
    sum(scores[rng.choice(p, size=kk, replace=False)].sum()
        for p, kk in zip(pools, counts)) / K
    for _ in range(B)
])
p_emp = (1 + int((draws >= res["t_obs"]).sum())) / (1 + B)
t_mc = time.perf_counter() - t0

print(f"N = {N:,} genes, K = {K}, B = {B}\n")
print(f"{'':16s}{'exact':>16s}{'resampling':>16s}{'rel. diff':>12s}")
print(f"{'null mean':16s}{res['null_mean']:16.6f}{draws.mean():16.6f}"
      f"{abs(draws.mean() - res['null_mean']) / res['null_mean']:11.2%}")
print(f"{'null SD':16s}{res['null_sd']:16.6f}{draws.std(ddof=1):16.6f}"
      f"{abs(draws.std(ddof=1) - res['null_sd']) / res['null_sd']:11.2%}")
print(f"{'skewness':16s}{res['null_skew']:16.6f}{'(unavailable)':>16s}")
print(f"\n{'z':16s}{res['z']:16.3f}")
print(f"{'p-value':16s}{res['p_edgeworth']:16.3e}{p_emp:16.3e}"
      f"   <- floor {1 / (1 + B):.3e}")
print(f"{'time':16s}{t_ex * 1000:14.2f}ms{t_mc * 1000:14.2f}ms{t_mc / t_ex:11.0f}x")
