"""Vignette 2: the median aggregate has an exact null too -- and no limit theorem."""
import itertools
import time

import numpy as np

from hetnetex_md import exact_median_pvalue

print("## 1. Verified against every possible sample\n")
rng = np.random.default_rng(4)
N, K = 12, 5
x = np.round(rng.gamma(2.0, 1.0, size=N), 6)
pool = np.arange(N)
meds = np.array([np.median(x[list(c)]) for c in itertools.combinations(range(N), K)])
uniq = np.unique(meds)
print(f"Enumerated all C({N},{K}) = {len(meds)} samples.\n")
print(f"{'threshold t':>14s}{'brute force':>20s}{'convolution':>20s}{'difference':>14s}")
for t in (uniq[0], uniq[len(uniq) // 2], uniq[-1]):
    b = float(np.mean(meds >= t))
    e = exact_median_pvalue(x, [pool], [K], t)
    print(f"{t:14.6f}{b:20.15f}{e:20.15f}{abs(b - e):14.2e}")

print("\n\n## 2. Stratified, also exact\n")
N1, K1, N2, K2 = 7, 3, 6, 2
y = np.round(rng.gamma(2.0, 1.0, size=N1 + N2), 6)
p1, p2 = np.arange(N1), np.arange(N1, N1 + N2)
meds2 = np.array([
    np.median(y[list(c1) + list(c2)])
    for c1 in itertools.combinations(p1, K1)
    for c2 in itertools.combinations(p2, K2)
])
u2 = np.unique(meds2)
print(f"{'threshold t':>14s}{'brute force':>20s}{'convolution':>20s}{'difference':>14s}")
for t in (u2[0], u2[len(u2) // 2], u2[-1]):
    b = float(np.mean(meds2 >= t))
    e = exact_median_pvalue(y, [p1, p2], [K1, K2], t)
    print(f"{t:14.6f}{b:20.15f}{e:20.15f}{abs(b - e):14.2e}")

print("\n\n## 3. No floor, and no CLT\n")
rng = np.random.default_rng(7)
N, K = 20_000, 151
scores = np.abs(rng.gamma(0.8, 1.0, size=N))
deg = rng.gamma(2.0, 5.0, size=N)
bins = np.minimum((np.argsort(np.argsort(deg)) / N * 10).astype(int), 9)
gene_idx = rng.choice(N, size=K, replace=False)
scores[gene_idx] += rng.gamma(2.2, 0.9, size=K)      # strong implanted signal

pools, counts = [], []
for b in np.unique(bins[gene_idx]):
    kk = int(np.sum(bins[gene_idx] == b))
    uni = np.ones(N, bool)
    uni[gene_idx] = False
    pools.append(np.flatnonzero((bins == b) & uni))
    counts.append(kk)

m_obs = float(np.median(scores[gene_idx]))
t0 = time.perf_counter()
p_exact = exact_median_pvalue(scores, pools, counts, m_obs)
t_ex = time.perf_counter() - t0

B = 1000
t0 = time.perf_counter()
draws = np.array([
    np.median(scores[np.concatenate([rng.choice(p, size=kk, replace=False)
                                     for p, kk in zip(pools, counts)])])
    for _ in range(B)
])
p_mc = (1 + int((draws >= m_obs).sum())) / (1 + B)
t_mc = time.perf_counter() - t0

print(f"N = {N:,}, K = {K}, B = {B}\n")
print(f"observed median      {m_obs:.6f}")
print(f"exact p              {p_exact:.3e}      (no distributional assumption)")
print(f"resampling p         {p_mc:.3e}      floor = {1 / (1 + B):.3e}")
print(f"time                 {t_ex * 1000:.2f} ms vs {t_mc * 1000:.2f} ms"
      f"  -> {t_mc / t_ex:.0f}x")

print("\n\n## 4. The caution: the median is often degenerate\n")
rng = np.random.default_rng(3)
sparse_scores = np.abs(rng.gamma(0.6, 1.0, size=N))
sparse_scores[rng.random(N) < 0.55] = 0.0
gi = rng.choice(N, size=K, replace=False)
frac_zero = float(np.mean(sparse_scores[gi] == 0))
print(f"fraction of the gene set with no path : {frac_zero:.0%}")
print(f"observed median                       : {np.median(sparse_scores[gi]):.6f}")
print("\nOn sparse metapaths more than half a gene set typically has no path to")
print("the target, so the median is exactly zero and the test is vacuous. This")
print("argues for the mean aggregate on scientific grounds, not computational ones.")
