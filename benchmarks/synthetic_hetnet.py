"""
run_experiments.py
==================
Benchmarks HetNetEX-MD against the resampling-based Multi-DWPC null and against
a degree-preserving network-permutation (XSwap) null on a synthetic hetnet whose
layer structure mirrors Hetionet.

Outputs: results CSVs + figures under ./figs, and a JSON of headline numbers.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr, pearsonr
from statsmodels.stats.multitest import multipletests

from hetnetex_md.core import (
    EdgeType,
    HetNet,
    network_null_moments,
    aggregate_network_null,
    exact_resampling_moments,
    edgeworth_upper_tail,
)

RNG = np.random.default_rng(20260727)


def stable_seed(*parts) -> int:
    """Deterministic 64-bit seed from arbitrary parts.

    Python's built-in hash() is randomised per process (PYTHONHASHSEED), so
    seeding an RNG with it makes results irreproducible across runs. BLAKE2b is
    stable across processes, machines, and Python versions.
    """
    import hashlib
    h = hashlib.blake2b(repr(parts).encode("utf-8"), digest_size=8)
    return int.from_bytes(h.digest(), "big")


OUT = Path("figs")
OUT.mkdir(exist_ok=True)
W = 0.4

NODE_COUNTS = dict(G=2000, C=400, PW=250, A=150, D=200, BP=120)

EDGE_SPEC = [
    ("GpPW", "G", "PW", 6000),
    ("GaA", "G", "A", 4500),
    ("GbC", "G", "C", 5000),
    ("GdD", "G", "D", 4000),
    ("CpPW", "C", "PW", 1600),
    ("CaA", "C", "A", 1100),
    ("AaD", "A", "D", 900),
    ("DaPW", "D", "PW", 1300),
    ("DbC", "D", "C", 900),
    ("PWpBP", "PW", "BP", 1400),
    ("AeBP", "A", "BP", 800),
    ("DdBP", "D", "BP", 1000),
    ("CbBP", "C", "BP", 900),
]

METAPATHS = {
    "GpPW-PWpBP": ["GpPW", "PWpBP"],
    "GaA-AeBP": ["GaA", "AeBP"],
    "GdD-DdBP": ["GdD", "DdBP"],
    "GbC-CbBP": ["GbC", "CbBP"],
    "GbC-CpPW-PWpBP": ["GbC", "CpPW", "PWpBP"],
    "GaA-AaD-DdBP": ["GaA", "AaD", "DdBP"],
    "GdD-DaPW-PWpBP": ["GdD", "DaPW", "PWpBP"],
    "GbC-CaA-AeBP": ["GbC", "CaA", "AeBP"],
    "GbC-CaA-AaD-DdBP": ["GbC", "CaA", "AaD", "DdBP"],
    "GaA-AaD-DaPW-PWpBP": ["GaA", "AaD", "DaPW", "PWpBP"],
    "GdD-DbC-CpPW-PWpBP": ["GdD", "DbC", "CpPW", "PWpBP"],
}


# --------------------------------------------------------------------------- #
# Network construction
# --------------------------------------------------------------------------- #
def powerlaw_weights(n: int, gamma: float, rng) -> np.ndarray:
    x = (rng.random(n) + 1e-3) ** (-1.0 / (gamma - 1.0))
    return x / x.sum()


def make_bipartite(n_src: int, n_dst: int, m: int, rng, gamma=2.1) -> sparse.csr_matrix:
    ps = powerlaw_weights(n_src, gamma, rng)
    pd_ = powerlaw_weights(n_dst, gamma, rng)
    seen = set()
    rows, cols = [], []
    while len(seen) < m:
        need = m - len(seen)
        r = rng.choice(n_src, size=need * 2, p=ps)
        c = rng.choice(n_dst, size=need * 2, p=pd_)
        for i, j in zip(r, c):
            key = i * n_dst + j
            if key not in seen:
                seen.add(key)
                rows.append(i)
                cols.append(j)
                if len(seen) >= m:
                    break
    A = sparse.csr_matrix(
        (np.ones(len(rows)), (rows, cols)), shape=(n_src, n_dst), dtype=np.float64
    )
    return A


def build_net(rng) -> HetNet:
    edges = {}
    for name, s, d, m in EDGE_SPEC:
        A = make_bipartite(NODE_COUNTS[s], NODE_COUNTS[d], m, rng)
        edges[name] = EdgeType(name, s, d, A)
    return HetNet(NODE_COUNTS, edges)


# --------------------------------------------------------------------------- #
# DWPC column for a single target, all genes  (backward sparse matvecs)
# --------------------------------------------------------------------------- #
def _weighted_ops(net: HetNet, metapath, w=W):
    ops = []
    for e in net.edges_for(metapath) if hasattr(net, "edges_for") else [net.edges[x] for x in metapath]:
        ds = np.where(e.d_src > 0, e.d_src, 1.0) ** (-w)
        dt = np.where(e.d_dst > 0, e.d_dst, 1.0) ** (-w)
        ops.append((sparse.diags(ds) @ e.A @ sparse.diags(dt)).tocsr())
    return ops


def dwpc_column(net: HetNet, metapath, target_pos: int, w=W) -> np.ndarray:
    ops = _weighted_ops(net, metapath, w)
    v = np.zeros(ops[-1].shape[1])
    v[target_pos] = 1.0
    for M in reversed(ops):
        v = M @ v
    return v


def dwpc_raw_mean(net: HetNet, metapath, w=W) -> float:
    ops = _weighted_ops(net, metapath, w)
    out = ops[0]
    for M in ops[1:]:
        out = out @ M
    out = np.asarray(out.todense())
    return float(out.mean())


# --------------------------------------------------------------------------- #
# Degree-preserving XSwap on a bipartite edge list
# --------------------------------------------------------------------------- #
def xswap(A: sparse.csr_matrix, rng, attempts_per_edge=10) -> sparse.csr_matrix:
    coo = A.tocoo()
    src = coo.row.copy()
    dst = coo.col.copy()
    n_dst = A.shape[1]
    edge_set = set((src * n_dst + dst).tolist())
    m = len(src)
    n_att = attempts_per_edge * m
    ia = rng.integers(0, m, size=n_att)
    ib = rng.integers(0, m, size=n_att)
    for k in range(n_att):
        i, j = ia[k], ib[k]
        if i == j:
            continue
        s1, d1 = src[i], dst[i]
        s2, d2 = src[j], dst[j]
        if s1 == s2 or d1 == d2:
            continue
        k1 = s1 * n_dst + d2
        k2 = s2 * n_dst + d1
        if k1 in edge_set or k2 in edge_set:
            continue
        edge_set.discard(s1 * n_dst + d1)
        edge_set.discard(s2 * n_dst + d2)
        edge_set.add(k1)
        edge_set.add(k2)
        dst[i], dst[j] = d2, d1
    return sparse.csr_matrix(
        (np.ones(m), (src, dst)), shape=A.shape, dtype=np.float64
    )


def permute_net(net: HetNet, rng) -> HetNet:
    edges = {}
    for name, e in net.edges.items():
        edges[name] = EdgeType(name, e.src_type, e.dst_type, xswap(e.A, rng))
    return HetNet(net.node_counts, edges)


# --------------------------------------------------------------------------- #
# Gene sets
# --------------------------------------------------------------------------- #
def build_gene_sets(net: HetNet, rng, n_sets=6):
    nG = NODE_COUNTS["G"]
    deg = np.zeros(nG)
    for e in net.edges.values():
        if e.src_type == "G":
            deg += e.d_src
        if e.dst_type == "G":
            deg += e.d_dst
    order = np.argsort(np.argsort(deg + rng.random(nG) * 1e-6))
    bin_ids = np.minimum((order / nG * 10).astype(int), 9)

    sizes = [40, 60, 90, 120, 160, 200]
    sets = []
    for i in range(n_sets):
        K = sizes[i]
        target = int(rng.integers(0, NODE_COUNTS["BP"]))
        genes = rng.choice(nG, size=K, replace=False)
        sets.append(dict(set_id=f"LV{i+1:02d}", genes=np.sort(genes), target=target, K=K))
    return sets, bin_ids, deg


def implant_signal(net: HetNet, sets, rng, frac=0.5, n_hubs=6):
    """Add edges so that some gene sets have genuine connectivity to their target."""
    for si, s in enumerate(sets):
        if si % 2 == 1:
            continue  # leave half as pure nulls
        t = s["target"]
        pw = np.flatnonzero(np.asarray(net.edges["PWpBP"].A[:, t].todense()).ravel() > 0)
        if pw.size == 0:
            pw = rng.choice(NODE_COUNTS["PW"], size=3, replace=False)
            A = net.edges["PWpBP"].A.tolil()
            for p in pw:
                A[p, t] = 1.0
            net.edges["PWpBP"].A = A.tocsr()
        hubs = pw[: min(n_hubs, pw.size)]
        A = net.edges["GpPW"].A.tolil()
        chosen = rng.choice(s["genes"], size=max(2, int(frac * s["K"])), replace=False)
        for g in chosen:
            A[g, rng.choice(hubs)] = 1.0
        net.edges["GpPW"].A = A.tocsr()
    return net


# --------------------------------------------------------------------------- #
# Experiment A : gene-set resampling null
# --------------------------------------------------------------------------- #
def experiment_resampling(net, sets, bin_ids, B_list=(200, 1000)):
    nG = NODE_COUNTS["G"]
    pooled = np.unique(np.concatenate([s["genes"] for s in sets]))
    pooled_mask = np.zeros(nG, dtype=bool)
    pooled_mask[pooled] = True

    rows = []
    timing = {f"emp_B{B}": 0.0 for B in B_list}
    timing["exact"] = 0.0

    score_cache = {}
    for mp_name, mp in METAPATHS.items():
        c = dwpc_raw_mean(net, mp)
        for s in sets:
            key = (mp_name, s["target"])
            if key not in score_cache:
                score_cache[key] = np.arcsinh(dwpc_column(net, mp, s["target"]) / c)
    for mp_name, mp in METAPATHS.items():
        for s in sets:
            x = score_cache[(mp_name, s["target"])]
            gidx = s["genes"]
            K = s["K"]
            t_obs = float(x[gidx].mean())

            for null_type in ("random", "permuted"):
                if null_type == "random":
                    universe = np.ones(nG, dtype=bool)
                    universe[gidx] = False
                else:
                    universe = pooled_mask.copy()

                bins = np.unique(bin_ids[gidx])
                pools, counts = [], []
                for b in bins:
                    k = int(np.sum(bin_ids[gidx] == b))
                    pool = np.flatnonzero((bin_ids == b) & universe)
                    if pool.size <= k:
                        pool = np.flatnonzero(bin_ids == b)
                    pools.append(pool)
                    counts.append(k)

                # ---- exact (HetNetEX-MD) ---- #
                t0 = time.perf_counter()
                mean_e, var_e, mu3_e = exact_resampling_moments(x, pools, counts)
                z_e, p_edge = edgeworth_upper_tail(t_obs, mean_e, var_e, mu3_e)
                timing["exact"] += time.perf_counter() - t0
                sd_e = float(np.sqrt(var_e))
                z_plain = (t_obs - mean_e) / sd_e if sd_e > 0 else np.inf

                rec = dict(
                    metapath=mp_name, set_id=s["set_id"], K=K, L=len(mp),
                    null_type=null_type, t_obs=t_obs,
                    ex_mean=mean_e, ex_sd=sd_e,
                    ex_skew=float(mu3_e / var_e ** 1.5) if var_e > 0 else 0.0,
                    ex_z=z_plain, ex_p=p_edge,
                )

                # ---- empirical (Multi-DWPC) ---- #
                for B in B_list:
                    rng = np.random.default_rng(stable_seed(mp_name, s["set_id"], null_type, B))
                    t0 = time.perf_counter()
                    draws = np.empty(B)
                    for bi in range(B):
                        tot = 0.0
                        for pool, k in zip(pools, counts):
                            tot += x[rng.choice(pool, size=k, replace=False)].sum()
                        draws[bi] = tot / K
                    exceed = int(np.sum(draws >= t_obs))
                    p_emp = (1 + exceed) / (1 + B)
                    timing[f"emp_B{B}"] += time.perf_counter() - t0
                    rec[f"emp_mean_B{B}"] = float(draws.mean())
                    rec[f"emp_sd_B{B}"] = float(draws.std(ddof=1))
                    rec[f"emp_z_B{B}"] = float(
                        (t_obs - draws.mean()) / draws.std(ddof=1)
                    ) if draws.std(ddof=1) > 0 else np.inf
                    rec[f"emp_p_B{B}"] = p_emp
                rows.append(rec)
    return pd.DataFrame(rows), timing


# --------------------------------------------------------------------------- #
# Experiment B : network-permutation null
# --------------------------------------------------------------------------- #
def experiment_network(net, sets, B_net=100):
    c_scale = {mp: dwpc_raw_mean(net, path) for mp, path in METAPATHS.items()}

    obs = {}
    for mp_name, mp in METAPATHS.items():
        for s in sets:
            x = np.arcsinh(dwpc_column(net, mp, s["target"]) / c_scale[mp_name])
            obs[(mp_name, s["set_id"])] = float(x[s["genes"]].mean())

    # ---- analytical ---- #
    t0 = time.perf_counter()
    ana = {}
    for mp_name, mp in METAPATHS.items():
        mom = network_null_moments(net, mp, w=W)
        for s in sets:
            m, v = aggregate_network_null(mom, s["genes"], s["target"], c_scale[mp_name])
            ana[(mp_name, s["set_id"])] = (m, np.sqrt(v))
    t_ana = time.perf_counter() - t0

    # ---- empirical ---- #
    rng = np.random.default_rng(7)
    draws = {k: [] for k in obs}
    t0 = time.perf_counter()
    for b in range(B_net):
        pnet = permute_net(net, rng)
        for mp_name, mp in METAPATHS.items():
            for s in sets:
                x = np.arcsinh(dwpc_column(pnet, mp, s["target"]) / c_scale[mp_name])
                draws[(mp_name, s["set_id"])].append(float(x[s["genes"]].mean()))
    t_emp = time.perf_counter() - t0

    rows = []
    for (mp_name, sid), d in draws.items():
        d = np.asarray(d)
        m_a, sd_a = ana[(mp_name, sid)]
        t_obs = obs[(mp_name, sid)]
        sd_e = d.std(ddof=1)
        rows.append(dict(
            metapath=mp_name, set_id=sid, L=len(METAPATHS[mp_name]),
            t_obs=t_obs,
            emp_mean=float(d.mean()), emp_sd=float(sd_e),
            emp_z=float((t_obs - d.mean()) / sd_e) if sd_e > 0 else np.inf,
            emp_p=float((1 + np.sum(d >= t_obs)) / (1 + B_net)),
            ana_mean=float(m_a), ana_sd=float(sd_a),
            ana_z=float((t_obs - m_a) / sd_a) if sd_a > 0 else np.inf,
        ))
    from scipy.stats import norm
    df = pd.DataFrame(rows)
    df["ana_p"] = norm.sf(df["ana_z"])
    return df, t_emp, t_ana, B_net


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #
def concordance(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    return pearsonr(a[ok], b[ok])[0], spearmanr(a[ok], b[ok])[0]


def main():
    rng = np.random.default_rng(11)
    net = build_net(rng)
    sets, bin_ids, deg = build_gene_sets(net, rng)
    net = implant_signal(net, sets, rng)

    print("== Experiment A: gene-set resampling null ==")
    dfA, timing = experiment_resampling(net, sets, bin_ids)
    dfA.to_csv("results_resampling.csv", index=False)

    print("== Experiment B: network permutation null ==")
    dfB, t_emp, t_ana, B_net = experiment_network(net, sets, B_net=60)
    dfB.to_csv("results_network.csv", index=False)

    summary = {}

    # ---- Experiment A summaries ---- #
    for B in (200, 1000):
        r_mean = concordance(dfA["ex_mean"].values, dfA[f"emp_mean_B{B}"].values)
        r_sd = concordance(dfA["ex_sd"].values, dfA[f"emp_sd_B{B}"].values)
        r_z = concordance(dfA["ex_z"].values, dfA[f"emp_z_B{B}"].values)
        lp_e = -np.log10(np.clip(dfA["ex_p"].values, 1e-300, 1))
        lp_m = -np.log10(np.clip(dfA[f"emp_p_B{B}"].values, 1e-300, 1))
        r_p = concordance(lp_e, lp_m)
        rel_mean = np.abs(dfA[f"emp_mean_B{B}"] - dfA["ex_mean"]) / np.abs(dfA["ex_mean"])
        rel_sd = np.abs(dfA[f"emp_sd_B{B}"] - dfA["ex_sd"]) / dfA["ex_sd"]
        floor_hits = int(np.sum(np.isclose(dfA[f"emp_p_B{B}"], 1.0 / (1 + B))))
        summary[f"A_B{B}"] = dict(
            r_mean=r_mean[0], rho_mean=r_mean[1],
            r_sd=r_sd[0], rho_sd=r_sd[1],
            r_z=r_z[0], rho_z=r_z[1],
            r_logp=r_p[0], rho_logp=r_p[1],
            med_rel_err_mean=float(np.median(rel_mean)),
            med_rel_err_sd=float(np.median(rel_sd)),
            p_floor_hits=floor_hits, n_features=len(dfA),
            floor=1.0 / (1 + B),
        )
        # selection agreement (BH within set_id x null_type)
        def bh(col):
            out = np.full(len(dfA), np.nan)
            for _, idx in dfA.groupby(["set_id", "null_type"]).groups.items():
                ii = dfA.index.get_indexer(idx)
                out[ii] = multipletests(dfA.loc[idx, col].values, method="fdr_bh")[1]
            return out
        sel_e = bh("ex_p") < 0.05
        sel_m = bh(f"emp_p_B{B}") < 0.05
        inter = int(np.sum(sel_e & sel_m))
        union = int(np.sum(sel_e | sel_m))
        summary[f"A_B{B}"]["n_sel_exact"] = int(sel_e.sum())
        summary[f"A_B{B}"]["n_sel_emp"] = int(sel_m.sum())
        summary[f"A_B{B}"]["jaccard"] = inter / union if union else 1.0
        summary[f"A_B{B}"]["n_disagree"] = union - inter

    summary["A_timing"] = {k: float(v) for k, v in timing.items()}
    summary["A_speedup_B200"] = timing["emp_B200"] / timing["exact"]
    summary["A_speedup_B1000"] = timing["emp_B1000"] / timing["exact"]

    # ---- Experiment B summaries ---- #
    r_mean = concordance(dfB["ana_mean"].values, dfB["emp_mean"].values)
    r_sd = concordance(dfB["ana_sd"].values, dfB["emp_sd"].values)
    r_z = concordance(dfB["ana_z"].values, dfB["emp_z"].values)
    summary["B"] = dict(
        r_mean=r_mean[0], rho_mean=r_mean[1],
        r_sd=r_sd[0], rho_sd=r_sd[1],
        r_z=r_z[0], rho_z=r_z[1],
        med_rel_err_mean=float(np.median(np.abs(dfB.emp_mean - dfB.ana_mean) / np.abs(dfB.ana_mean))),
        med_rel_err_sd=float(np.median(np.abs(dfB.emp_sd - dfB.ana_sd) / dfB.ana_sd)),
        t_emp=t_emp, t_ana=t_ana, speedup=t_emp / t_ana, B=B_net,
        n_features=len(dfB),
    )
    for L in sorted(dfB.L.unique()):
        sub = dfB[dfB.L == L]
        summary["B"][f"rho_z_L{L}"] = float(spearmanr(sub.ana_z, sub.emp_z)[0])

    with open("summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=float)
    print(json.dumps(summary, indent=2, default=float))

    make_figures(dfA, dfB, summary)


def make_figures(dfA, dfB, summary):
    plt.rcParams.update({"font.size": 8, "figure.dpi": 200})

    # Fig 1: exact vs empirical null mean / sd / z / p (resampling null)
    fig, ax = plt.subplots(1, 4, figsize=(11, 2.7))
    B = 1000
    sub = dfA
    colors = {"random": "#2b7bba", "permuted": "#d95f02"}
    for nt, g in sub.groupby("null_type"):
        ax[0].scatter(g["ex_mean"], g[f"emp_mean_B{B}"], s=8, alpha=.7, c=colors[nt], label=nt)
        ax[1].scatter(g["ex_sd"], g[f"emp_sd_B{B}"], s=8, alpha=.7, c=colors[nt])
        ax[2].scatter(g["ex_z"], g[f"emp_z_B{B}"], s=8, alpha=.7, c=colors[nt])
        ax[3].scatter(-np.log10(np.clip(g["ex_p"], 1e-300, 1)),
                      -np.log10(np.clip(g[f"emp_p_B{B}"], 1e-300, 1)),
                      s=8, alpha=.7, c=colors[nt])
    for a, t in zip(ax, ["null mean", "null SD", "z", r"$-\log_{10}p$"]):
        lo = min(a.get_xlim()[0], a.get_ylim()[0]); hi = max(a.get_xlim()[1], a.get_ylim()[1])
        a.plot([lo, hi], [lo, hi], "k--", lw=.7)
        a.set_xlabel(f"HetNetEX-MD (exact) {t}")
        a.set_ylabel(f"Multi-DWPC B={B} {t}")
    ax[3].axhline(np.log10(B + 1), color="r", ls=":", lw=.8)
    ax[0].legend(frameon=False, fontsize=6)
    fig.tight_layout()
    fig.savefig(OUT / "fig1_resampling_agreement.png")
    plt.close(fig)

    # Fig 2: Monte-Carlo error of the null SD vs B, and p-value floor
    fig, ax = plt.subplots(1, 2, figsize=(6.4, 2.6))
    for B, col in zip((200, 1000), ("#7570b3", "#1b9e77")):
        rel = np.abs(dfA[f"emp_sd_B{B}"] - dfA["ex_sd"]) / dfA["ex_sd"]
        ax[0].hist(rel, bins=30, alpha=.6, color=col, label=f"B={B}")
    ax[0].set_xlabel("relative error of Monte-Carlo null SD")
    ax[0].set_ylabel("features")
    ax[0].legend(frameon=False, fontsize=6)

    ax[1].scatter(-np.log10(np.clip(dfA["ex_p"], 1e-300, 1)),
                  -np.log10(np.clip(dfA["emp_p_B200"], 1e-300, 1)),
                  s=8, alpha=.6, c="#7570b3", label="B=200")
    ax[1].scatter(-np.log10(np.clip(dfA["ex_p"], 1e-300, 1)),
                  -np.log10(np.clip(dfA["emp_p_B1000"], 1e-300, 1)),
                  s=8, alpha=.6, c="#1b9e77", label="B=1000")
    ax[1].axhline(np.log10(201), color="#7570b3", ls=":", lw=.9)
    ax[1].axhline(np.log10(1001), color="#1b9e77", ls=":", lw=.9)
    ax[1].set_xlabel(r"HetNetEX-MD $-\log_{10}p$")
    ax[1].set_ylabel(r"Multi-DWPC $-\log_{10}p$")
    ax[1].legend(frameon=False, fontsize=6)
    fig.tight_layout()
    fig.savefig(OUT / "fig2_resolution.png")
    plt.close(fig)

    # Fig 3: network null agreement by path length
    fig, ax = plt.subplots(1, 3, figsize=(8.4, 2.7))
    ms = {2: "o", 3: "s", 4: "^"}
    for L, g in dfB.groupby("L"):
        ax[0].scatter(g.ana_mean, g.emp_mean, s=12, marker=ms[L], alpha=.75, label=f"L={L}")
        ax[1].scatter(g.ana_sd, g.emp_sd, s=12, marker=ms[L], alpha=.75)
        ax[2].scatter(g.ana_z, g.emp_z, s=12, marker=ms[L], alpha=.75)
    for a, t in zip(ax, ["null mean", "null SD", "z"]):
        lo = min(a.get_xlim()[0], a.get_ylim()[0]); hi = max(a.get_xlim()[1], a.get_ylim()[1])
        a.plot([lo, hi], [lo, hi], "k--", lw=.7)
        a.set_xlabel(f"HetNetEX-MD {t}")
        a.set_ylabel(f"XSwap {t}")
    ax[0].legend(frameon=False, fontsize=6)
    fig.tight_layout()
    fig.savefig(OUT / "fig3_network_agreement.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
