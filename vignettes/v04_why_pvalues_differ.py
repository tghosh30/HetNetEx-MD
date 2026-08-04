"""Vignette 4: why permutation and analytical p-values cannot agree."""
import numpy as np

from hetnetex_md import hetnetex_md_resampling

rng = np.random.default_rng(42)

print("## 1. A permutation p-value is a random variable\n")
p_star = 0.002                      # the true tail probability
B = 1000
reps = 20
counts = rng.binomial(B, p_star, size=reps)
phat = (1 + counts) / (1 + B)
print(f"True tail probability p* = {p_star}\n")
print(f"Twenty independent runs of the SAME analysis, B = {B}:\n")
print("  " + "  ".join(f"{v:.4f}" for v in phat[:10]))
print("  " + "  ".join(f"{v:.4f}" for v in phat[10:]))
print(f"\n  range      {phat.min():.4f} to {phat.max():.4f}")
print(f"  mean       {phat.mean():.5f}   (theory: {p_star + (1-p_star)/(1+B):.5f})")
print(f"  SD         {phat.std(ddof=1):.5f}   "
      f"(theory: {np.sqrt(p_star*(1-p_star)/B):.5f})")
print("\nSame data, same code, different seed -> different p-value.")
print("The analytical p-value is deterministic: one number, every time.")

print("\n\n## 2. Relative error blows up exactly where it matters\n")
print(f"{'p*':>12s}{'rel. std. error':>20s}{'reportable at B=1000?':>26s}")
for ps in (0.5, 0.1, 0.01, 1e-3, 1e-4, 1e-6):
    rse = np.sqrt((1 - ps) / (B * ps))
    ok = "yes" if ps >= 1 / (1 + B) else "NO - below the floor"
    print(f"{ps:12.0e}{rse:19.1%}{ok:>26s}")
print("\nThe smaller the true p-value -- i.e. the stronger the finding -- the")
print("worse resampling estimates it, until it cannot represent it at all.")

print("\n\n## 3. The floor meets multiple-testing correction\n")
alpha = 0.05
print(f"BH needs p_(k) <= alpha*k/m, but the floor gives p >= 1/(B+1).")
print(f"So NO discovery is possible unless  B >= m/alpha - 1.\n")
print(f"{'metapaths m':>14s}{'B required':>14s}{'B=200?':>12s}{'B=1000?':>12s}")
for m in (5, 11, 25, 50, 100, 200, 500):
    need = m / alpha - 1
    print(f"{m:14d}{need:14.0f}{'ok' if 200 >= need else 'BLOCKED':>12s}"
          f"{'ok' if 1000 >= need else 'BLOCKED':>12s}")

print("\n\n## 4. What this looks like on real features\n")
N = 20_000
scores = np.abs(rng.gamma(0.6, 1.0, size=N))
scores[rng.random(N) < 0.5] = 0.0
deg = rng.gamma(2.0, 5.0, size=N)
bins = np.minimum((np.argsort(np.argsort(deg)) / N * 10).astype(int), 9)

m_tests = 11
B = 200
print(f"{m_tests} metapaths, B = {B}, FDR alpha = {alpha}")
print(f"BH threshold for the top feature = alpha/m = {alpha/m_tests:.5f}")
print(f"Permutation floor                = 1/(B+1) = {1/(1+B):.5f}")
print(f"-> the floor is {'ABOVE' if 1/(1+B) > alpha/m_tests else 'below'} "
      f"the threshold, so the best possible permutation p-value "
      f"{'CANNOT' if 1/(1+B) > alpha/m_tests else 'can'} clear it.\n")

rows = []
for j in range(m_tests):
    gi = rng.choice(N, size=150, replace=False)
    boost = [0, 0, 0, 0, 0, 0, 0, 0, 1.1, 1.4, 1.8][j]
    sc = scores.copy()
    if boost:
        sc[gi[:70]] += rng.gamma(boost, 0.9, size=70)
    r = hetnetex_md_resampling(sc, gi, bins)
    rows.append((j, r["z"], r["p_edgeworth"]))

print(f"{'feature':>9s}{'z':>10s}{'analytic p':>14s}"
      f"{'best possible perm p':>24s}")
for j, z, p in sorted(rows, key=lambda r: -r[1])[:6]:
    print(f"{j:9d}{z:10.2f}{p:14.2e}{1/(1+B):24.5f}")

sig = sum(1 for _, z, p in rows if p < alpha / m_tests)
print(f"\nanalytic, clearing the strictest BH bar : {sig} of {m_tests}")
print(f"permutation at B={B}, best case          : 0 of {m_tests} "
      f"(floor exceeds the bar)")

print("\n\n## 5. What DOES agree\n")
print("Both p-values are monotone in the same z. The ranking transfers even")
print("though the numbers never match:")
zs = np.array([r[1] for r in rows])
ps = np.array([r[2] for r in rows])
order_z = np.argsort(-zs)
order_p = np.argsort(ps)
print(f"\n  ranking by z      : {[int(i) for i in order_z]}")
print(f"  ranking by p      : {[int(i) for i in order_p]}")
print(f"  identical         : {np.array_equal(order_z, order_p)}")
print("\nThe scientific decision is 'which metapaths are most promising',")
print("not 'is this p-value 0.0031 or 0.0028'. That decision is preserved.")
