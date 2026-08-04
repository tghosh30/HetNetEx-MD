"""Vignette 3: the network-permutation null, and why the edge probability matters."""
import numpy as np
from scipy import sparse

from hetnetex_md import (
    EdgeType,
    HetNet,
    admissible_self_kernels,
    fit_soft_cm,
    network_null_moments,
)


def powerlaw(n, gamma, rng):
    x = (rng.random(n) + 1e-3) ** (-1.0 / (gamma - 1.0))
    return x / x.sum()


def bipartite(n_src, n_dst, m, rng, gamma=2.1):
    ps, pd = powerlaw(n_src, gamma, rng), powerlaw(n_dst, gamma, rng)
    seen, rows, cols = set(), [], []
    while len(seen) < m:
        r = rng.choice(n_src, size=2 * (m - len(seen)), p=ps)
        c = rng.choice(n_dst, size=2 * (m - len(seen)), p=pd)
        for i, j in zip(r, c):
            key = i * n_dst + j
            if key not in seen:
                seen.add(key)
                rows.append(i)
                cols.append(j)
                if len(seen) >= m:
                    break
    A = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n_src, n_dst))
    A.data[:] = 1.0
    return A


rng = np.random.default_rng(11)
sizes = {"G": 1500, "P": 220, "B": 110}
net = HetNet(sizes, {
    "GpP": EdgeType("GpP", "G", "P", bipartite(1500, 220, 4500, rng)),
    "PpB": EdgeType("PpB", "P", "B", bipartite(220, 110, 1200, rng)),
})
metapath = ["GpP", "PpB"]
W = 0.4

print("## 1. The null mean is rank one (Lemma 3)\n")
mom = network_null_moments(net, metapath, w=W)
mu = mom["C"] * np.outer(mom["a"], mom["b"])
print(f"mu has shape {mu.shape} = {mu.size:,} entries")
print(f"numerical rank: {np.linalg.matrix_rank(mu, tol=1e-10 * mu.max())}")
print("\nSo one O(Ln) pass gives every gene against every target. The whole")
print("matrix is reconstructed from a scalar C and two vectors a, b:")
print(f"  C            = {mom['C']:.6e}")
print(f"  a (sources)  = {mom['a'].shape[0]:,} values")
print(f"  b (targets)  = {mom['b'].shape[0]:,} values")
print(f"  stored       = {1 + mom['a'].size + mom['b'].size:,} numbers "
      f"instead of {mu.size:,}")

print("\n\n## 2. The edge-probability ladder (Lemma 4 -> 4b)\n")
e = net.edges["GpP"]
fx, fy, resid = fit_soft_cm(e.d_src, e.d_dst)
print(f"max-entropy fit: median relative degree residual = {resid:.2e}")
print("\nCompare the three models on representative degree pairs:\n")
print(f"{'d_u':>6s}{'d_v':>6s}{'sparse d_u d_v/m':>20s}"
      f"{'1 - exp(-x)':>16s}{'max-entropy':>16s}")
m = e.m
for du, dv in [(2, 3), (10, 20), (40, 60), (100, 90), (200, 150)]:
    base = du * dv / m
    expo = 1.0 - np.exp(-base)
    xu = float(np.atleast_1d(fx(np.array([float(du)])))[0])
    yv = float(np.atleast_1d(fy(np.array([float(dv)])))[0])
    soft = xu * yv / (1 + xu * yv)
    print(f"{du:6d}{dv:6d}{base:20.6f}{expo:16.6f}{soft:16.6f}")
print("\nThe three agree when d_u d_v / m is small and separate as it grows.")
print("Hub-to-hub pairs are exactly where they differ -- and where the")
print("sparse limit can exceed 1, which is not a probability.")

print("\n\n## 3. Inter-source dependence is not a correction (Theorem 6)\n")
theta = mom["theta_mid"] + float(np.atleast_1d(mom["theta_L"])[0])
print(f"theta (one scalar per metapath+target) = {theta:.6e}\n")
print(f"{'gene set size K':>16s}{'independent SD':>18s}{'with theta':>16s}{'ratio':>10s}")
target = 0
for K in (10, 50, 100, 200, 400):
    g = np.arange(K)
    mu_g = mom["C"] * mom["a"][g] * mom["b"][target]
    mu2 = mom["C2"] * mom["a2"][g] * mom["b2"][target]
    th1 = np.atleast_1d(mom["theta_1"])[g]
    self_th = admissible_self_kernels(len(metapath), th1, mom["theta_mid"],
                                      float(np.atleast_1d(mom["theta_L"])[target]))
    vg = np.clip(mu2 + mu_g ** 2 * self_th, 0, None)
    sd_i = np.sqrt(vg.sum()) / K
    sd_f = np.sqrt((vg.sum() + theta * (mu_g.sum() ** 2 - (mu_g ** 2).sum()))) / K
    print(f"{K:16d}{sd_i:18.6e}{sd_f:16.6e}{sd_f / sd_i:10.2f}")
print("\nThe ratio grows with K: assuming genes are independent understates the")
print("null spread more and more as the gene set gets larger.")

print("\n\n## 4. The admissibility rule\n")
print(f"{'L':>4s}{'admissible same-source kernels':>36s}")
for L in (2, 3, 4, 5):
    v = admissible_self_kernels(L, 0.11, 0.22, 0.33)
    which = {2: "none", 3: "theta_1 + theta_L",
             4: "theta_1 + theta_mid + theta_L"}.get(L, "all")
    print(f"{L:4d}{which:>36s}   value = {v:.2f}")
print("\nForcing layer l shared pins interior positions l-1 and l. If no interior")
print("position is left free the two paths coincide, so the term must not be")
print("counted. At L=2 nothing is admissible and Var(Y_g) = mu^(2) exactly.")
