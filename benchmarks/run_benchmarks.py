"""
run_final.py
============
Final integrated benchmark for HetNetEX-MD.

  Experiment A : gene-set resampling null (Multi-DWPC LV track)      -> exact
  Experiment B : network permutation null (Multi-DWPC year track)    -> asymptotic
  Experiment C : finite-size scaling of the configuration-model error
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
from scipy.stats import spearmanr, pearsonr, norm
from statsmodels.stats.multitest import multipletests

from hetnetex_md.core import (
    network_null_moments, correction_factors,
    exact_resampling_moments, edgeworth_upper_tail,
)
from synthetic_hetnet import (
    stable_seed,
    METAPATHS, W, NODE_COUNTS, EDGE_SPEC, build_net, build_gene_sets,
    implant_signal, dwpc_column, dwpc_raw_mean, permute_net, HetNet, EdgeType,
    make_bipartite,
)

OUT = Path("figs")
OUT.mkdir(exist_ok=True)
B_NET = 100


def conc(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    return float(pearsonr(a[ok], b[ok])[0]), float(spearmanr(a[ok], b[ok])[0])


# --------------------------------------------------------------------------- #
def analytic_raw(mom, corr, genes, t):
    gm, g1, gL = corr
    mu = mom["C"] * mom["a"][genes] * mom["b"][t] * (gm * g1[genes] * gL[t])
    mu2 = mom["C2"] * mom["a2"][genes] * mom["b2"][t] * (gm * g1[genes] * gL[t])
    th = mom["theta_mid"] + float(np.atleast_1d(mom["theta_L"])[t])
    th1 = np.atleast_1d(mom["theta_1"])[genes]
    var_g = np.clip(mu2 + mu ** 2 * (th1 + th), 0, None)
    K = len(genes)
    var = (var_g.sum() + th * (mu.sum() ** 2 - (mu ** 2).sum())) / K ** 2
    var_ind = var_g.sum() / K ** 2
    return float(mu.mean()), float(np.sqrt(max(var, 0))), float(np.sqrt(max(var_ind, 0))), float(th)


# --------------------------------------------------------------------------- #
def experiment_A(net, sets, bin_ids, B_list=(200, 1000)):
    nG = NODE_COUNTS["G"]
    pooled = np.unique(np.concatenate([s["genes"] for s in sets]))
    pooled_mask = np.zeros(nG, bool)
    pooled_mask[pooled] = True

    cache = {}
    for mp, path in METAPATHS.items():
        c = dwpc_raw_mean(net, path)
        for s in sets:
            cache[(mp, s["target"])] = np.arcsinh(dwpc_column(net, path, s["target"]) / c)

    rows = []
    timing = {"exact": 0.0, **{f"emp_B{B}": 0.0 for B in B_list}}
    for mp, path in METAPATHS.items():
        for s in sets:
            x = cache[(mp, s["target"])]
            gidx, K = s["genes"], s["K"]
            t_obs = float(x[gidx].mean())
            for null_type in ("random", "permuted"):
                universe = np.ones(nG, bool)
                if null_type == "random":
                    universe[gidx] = False
                else:
                    universe = pooled_mask.copy()
                pools, counts = [], []
                for b in np.unique(bin_ids[gidx]):
                    k = int(np.sum(bin_ids[gidx] == b))
                    pool = np.flatnonzero((bin_ids == b) & universe)
                    if pool.size <= k:
                        pool = np.flatnonzero(bin_ids == b)
                    pools.append(pool)
                    counts.append(k)

                t0 = time.perf_counter()
                m_e, v_e, m3_e = exact_resampling_moments(x, pools, counts)
                z_e, p_e = edgeworth_upper_tail(t_obs, m_e, v_e, m3_e)
                timing["exact"] += time.perf_counter() - t0
                sd_e = float(np.sqrt(v_e))
                rec = dict(metapath=mp, set_id=s["set_id"], K=K, L=len(path),
                           null_type=null_type, t_obs=t_obs, ex_mean=m_e, ex_sd=sd_e,
                           ex_skew=float(m3_e / v_e ** 1.5) if v_e > 0 else 0.0,
                           ex_z=(t_obs - m_e) / sd_e if sd_e > 0 else np.inf,
                           ex_p=p_e,
                           ex_p_normal=float(norm.sf((t_obs - m_e) / sd_e)) if sd_e > 0 else 0.0)
                for B in B_list:
                    rng = np.random.default_rng(
                        stable_seed(mp, s["set_id"], null_type, B))
                    t0 = time.perf_counter()
                    draws = np.empty(B)
                    for bi in range(B):
                        tot = 0.0
                        for pool, k in zip(pools, counts):
                            tot += x[rng.choice(pool, size=k, replace=False)].sum()
                        draws[bi] = tot / K
                    exceed = int(np.sum(draws >= t_obs))
                    timing[f"emp_B{B}"] += time.perf_counter() - t0
                    sd = draws.std(ddof=1)
                    rec[f"emp_mean_B{B}"] = float(draws.mean())
                    rec[f"emp_sd_B{B}"] = float(sd)
                    rec[f"emp_z_B{B}"] = float((t_obs - draws.mean()) / sd) if sd > 0 else np.inf
                    rec[f"emp_p_B{B}"] = (1 + exceed) / (1 + B)
                rows.append(rec)
    return pd.DataFrame(rows), timing


# --------------------------------------------------------------------------- #
def experiment_B(net, sets):
    obs, ana = {}, {}
    t0 = time.perf_counter()
    for mp, path in METAPATHS.items():
        mom = network_null_moments(net, path, w=W)
        corr = correction_factors(net, path, w=W)
        for s in sets:
            ana[(mp, s["set_id"])] = analytic_raw(mom, corr, s["genes"], s["target"])
    t_ana = time.perf_counter() - t0

    for mp, path in METAPATHS.items():
        for s in sets:
            y = dwpc_column(net, path, s["target"])
            obs[(mp, s["set_id"])] = float(y[s["genes"]].mean())

    draws = {k: [] for k in obs}
    rng = np.random.default_rng(707)
    t0 = time.perf_counter()
    for b in range(B_NET):
        pnet = permute_net(net, rng)
        for mp, path in METAPATHS.items():
            for s in sets:
                y = dwpc_column(pnet, path, s["target"])
                draws[(mp, s["set_id"])].append(float(y[s["genes"]].mean()))
    t_emp = time.perf_counter() - t0

    rows = []
    for (mp, sid), d in draws.items():
        d = np.asarray(d)
        m_a, sd_a, sd_i, th = ana[(mp, sid)]
        t_obs = obs[(mp, sid)]
        sd_e = float(d.std(ddof=1))
        rows.append(dict(metapath=mp, set_id=sid, L=len(METAPATHS[mp]), theta=th,
                         t_obs=t_obs, emp_mean=float(d.mean()), emp_sd=sd_e,
                         emp_z=(t_obs - d.mean()) / sd_e if sd_e > 0 else np.inf,
                         emp_p=float((1 + np.sum(d >= t_obs)) / (1 + B_NET)),
                         ana_mean=m_a, ana_sd=sd_a, ana_sd_indep=sd_i,
                         ana_z=(t_obs - m_a) / sd_a if sd_a > 0 else np.inf))
    df = pd.DataFrame(rows)
    df["ana_p"] = norm.sf(df["ana_z"])
    return df, t_emp, t_ana


# --------------------------------------------------------------------------- #
def experiment_C(scales=(0.5, 1.0, 2.0, 4.0), mp="GbC-CpPW-PWpBP", B=50):
    """Finite-size scaling of |E_emp - mu_ana| / mu_ana."""
    rows = []
    path = METAPATHS[mp]
    for sc in scales:
        rng = np.random.default_rng(int(1000 * sc) + 5)
        nc = {k: max(30, int(v * sc)) for k, v in NODE_COUNTS.items()}
        edges = {}
        for name, s, d, m in EDGE_SPEC:
            edges[name] = EdgeType(name, s, d,
                                   make_bipartite(nc[s], nc[d], int(m * sc), rng))
        net = HetNet(nc, edges)
        genes = rng.choice(nc["G"], size=min(120, nc["G"] // 4), replace=False)
        t = int(rng.integers(0, nc["BP"]))
        mom = network_null_moments(net, path, w=W)
        corr = correction_factors(net, path, w=W)
        m_a, sd_a, _, _ = analytic_raw(mom, corr, genes, t)
        r2 = np.random.default_rng(99)
        d = np.array([float(dwpc_column(permute_net(net, r2), path, t)[genes].mean())
                      for _ in range(B)])
        rows.append(dict(scale=sc, n_gene=nc["G"], ana=m_a, emp=float(d.mean()),
                         rel_err=abs(float(d.mean()) - m_a) / m_a))
        print(f"   scale {sc}: nG={nc['G']} rel_err={rows[-1]['rel_err']:.4f}")
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
def main():
    rng = np.random.default_rng(11)
    net = build_net(rng)
    sets, bin_ids, deg = build_gene_sets(net, rng)
    net = implant_signal(net, sets, rng)

    print("== A: resampling null ==")
    dfA, timing = experiment_A(net, sets, bin_ids)
    dfA.to_csv("resA.csv", index=False)

    print("== B: network null ==")
    dfB, t_emp, t_ana = experiment_B(net, sets)
    dfB.to_csv("resB.csv", index=False)

    print("== C: finite-size scaling ==")
    dfC = experiment_C()
    dfC.to_csv("resC.csv", index=False)

    S = {}
    for B in (200, 1000):
        s = {}
        s["r_mean"], s["rho_mean"] = conc(dfA.ex_mean, dfA[f"emp_mean_B{B}"])
        s["r_sd"], s["rho_sd"] = conc(dfA.ex_sd, dfA[f"emp_sd_B{B}"])
        s["r_z"], s["rho_z"] = conc(dfA.ex_z, dfA[f"emp_z_B{B}"])
        lp_e = -np.log10(np.clip(dfA.ex_p, 1e-300, 1))
        lp_m = -np.log10(np.clip(dfA[f"emp_p_B{B}"], 1e-300, 1))
        s["r_logp"], s["rho_logp"] = conc(lp_e, lp_m)
        s["med_rel_mean"] = float(np.median(np.abs(dfA[f"emp_mean_B{B}"] - dfA.ex_mean) / np.abs(dfA.ex_mean)))
        s["med_rel_sd"] = float(np.median(np.abs(dfA[f"emp_sd_B{B}"] - dfA.ex_sd) / dfA.ex_sd))
        s["max_rel_sd"] = float(np.max(np.abs(dfA[f"emp_sd_B{B}"] - dfA.ex_sd) / dfA.ex_sd))
        s["floor"] = 1.0 / (1 + B)
        s["n_at_floor"] = int(np.sum(np.isclose(dfA[f"emp_p_B{B}"], 1.0 / (1 + B))))
        s["n_p_identical"] = int(np.sum(np.isclose(dfA.ex_p, dfA[f"emp_p_B{B}"], rtol=1e-6)))

        def bh(col):
            out = np.full(len(dfA), np.nan)
            for _, idx in dfA.groupby(["set_id", "null_type"]).groups.items():
                ii = dfA.index.get_indexer(idx)
                out[ii] = multipletests(dfA.loc[idx, col].values, method="fdr_bh")[1]
            return out
        se, sm = bh("ex_p") < 0.05, bh(f"emp_p_B{B}") < 0.05
        s["n_sel_exact"], s["n_sel_emp"] = int(se.sum()), int(sm.sum())
        s["jaccard"] = float((se & sm).sum() / max((se | sm).sum(), 1))
        s["n_discordant"] = int((se | sm).sum() - (se & sm).sum())
        S[f"A_B{B}"] = s

    S["A_timing"] = {k: float(v) for k, v in timing.items()}
    S["A_speedup_B200"] = timing["emp_B200"] / timing["exact"]
    S["A_speedup_B1000"] = timing["emp_B1000"] / timing["exact"]
    S["A_n_features"] = int(len(dfA))
    S["A_med_skew"] = float(np.median(np.abs(dfA.ex_skew)))
    S["A_edgeworth_vs_normal"] = float(np.median(np.abs(
        -np.log10(np.clip(dfA.ex_p, 1e-300, 1)) + np.log10(np.clip(dfA.ex_p_normal, 1e-300, 1)))))

    b = dict(B=B_NET, n_features=int(len(dfB)), t_emp=t_emp, t_ana=t_ana,
             speedup=t_emp / t_ana)
    b["r_mean"], b["rho_mean"] = conc(dfB.ana_mean, dfB.emp_mean)
    b["r_sd"], b["rho_sd"] = conc(dfB.ana_sd, dfB.emp_sd)
    b["r_z"], b["rho_z"] = conc(dfB.ana_z, dfB.emp_z)
    b["med_rel_mean"] = float(np.median(np.abs(dfB.emp_mean - dfB.ana_mean) / dfB.ana_mean))
    b["med_rel_sd"] = float(np.median(np.abs(dfB.emp_sd - dfB.ana_sd) / dfB.ana_sd))
    b["ratio_sd_indep"] = float(np.median(dfB.emp_sd / dfB.ana_sd_indep))
    b["ratio_sd_full"] = float(np.median(dfB.emp_sd / dfB.ana_sd))
    b["theta_inflation"] = float(np.median(dfB.ana_sd / dfB.ana_sd_indep))
    for L in sorted(dfB.L.unique()):
        sub = dfB[dfB.L == L]
        b[f"L{L}"] = dict(n=int(len(sub)),
                          rho_mean=float(spearmanr(sub.ana_mean, sub.emp_mean)[0]),
                          rho_z=float(spearmanr(sub.ana_z, sub.emp_z)[0]),
                          med_rel_mean=float(np.median(np.abs(sub.emp_mean - sub.ana_mean) / sub.ana_mean)),
                          theta_inflation=float(np.median(sub.ana_sd / sub.ana_sd_indep)))
    S["B"] = b
    S["C"] = dfC.to_dict(orient="records")

    with open("summary_final.json", "w") as f:
        json.dump(S, f, indent=2, default=float)
    print(json.dumps(S, indent=2, default=float))
    figures(dfA, dfB, dfC)


def figures(dfA, dfB, dfC):
    plt.rcParams.update({"font.size": 7.5, "figure.dpi": 220,
                         "axes.spines.top": False, "axes.spines.right": False})
    col = {"random": "#1f78b4", "permuted": "#e6550d"}

    # Figure 1
    fig, ax = plt.subplots(1, 4, figsize=(10.6, 2.5))
    B = 1000
    for nt, g in dfA.groupby("null_type"):
        ax[0].scatter(g.ex_mean, g[f"emp_mean_B{B}"], s=7, alpha=.75, c=col[nt], label=nt, lw=0)
        ax[1].scatter(g.ex_sd, g[f"emp_sd_B{B}"], s=7, alpha=.75, c=col[nt], lw=0)
        ax[2].scatter(g.ex_z, g[f"emp_z_B{B}"], s=7, alpha=.75, c=col[nt], lw=0)
        ax[3].scatter(-np.log10(np.clip(g.ex_p, 1e-300, 1)),
                      -np.log10(np.clip(g[f"emp_p_B{B}"], 1e-300, 1)),
                      s=7, alpha=.75, c=col[nt], lw=0)
    labs = ["null mean", "null SD", "z", r"$-\log_{10}p$"]
    for a, t in zip(ax, labs):
        lo = min(a.get_xlim()[0], a.get_ylim()[0]); hi = max(a.get_xlim()[1], a.get_ylim()[1])
        a.plot([lo, hi], [lo, hi], "k--", lw=.6)
        a.set_xlabel("HetNetEX-MD " + t); a.set_ylabel(f"Multi-DWPC (B={B}) " + t)
    ax[3].axhline(np.log10(B + 1), color="k", ls=":", lw=.8)
    ax[0].legend(frameon=False, fontsize=6, loc="upper left")
    fig.tight_layout(); fig.savefig(OUT / "fig1_resampling.png"); plt.close(fig)

    # Figure 2
    fig, ax = plt.subplots(1, 3, figsize=(8.4, 2.5))
    for B, c in zip((200, 1000), ("#7570b3", "#1b9e77")):
        rel = np.abs(dfA[f"emp_sd_B{B}"] - dfA.ex_sd) / dfA.ex_sd
        ax[0].hist(rel, bins=28, alpha=.62, color=c, label=f"B={B}")
    ax[0].set_xlabel("relative error of Monte-Carlo null SD"); ax[0].set_ylabel("features")
    ax[0].legend(frameon=False, fontsize=6)
    for B, c in zip((200, 1000), ("#7570b3", "#1b9e77")):
        ax[1].scatter(-np.log10(np.clip(dfA.ex_p, 1e-300, 1)),
                      -np.log10(np.clip(dfA[f"emp_p_B{B}"], 1e-300, 1)),
                      s=7, alpha=.6, c=c, lw=0, label=f"B={B}")
        ax[1].axhline(np.log10(B + 1), color=c, ls=":", lw=.9)
    ax[1].set_xlabel(r"HetNetEX-MD $-\log_{10}p$")
    ax[1].set_ylabel(r"Multi-DWPC $-\log_{10}p$")
    ax[1].legend(frameon=False, fontsize=6)
    ax[2].scatter(dfA.ex_skew, -np.log10(np.clip(dfA.ex_p, 1e-300, 1))
                  + np.log10(np.clip(dfA.ex_p_normal, 1e-300, 1)), s=7, alpha=.6, c="k", lw=0)
    ax[2].axhline(0, color="r", ls="--", lw=.7)
    ax[2].set_xlabel("null skewness  $\\gamma_1$")
    ax[2].set_ylabel(r"$\log_{10}(p_{\rm normal}/p_{\rm Edgeworth})$")
    fig.tight_layout(); fig.savefig(OUT / "fig2_resolution.png"); plt.close(fig)

    # Figure 3
    fig, ax = plt.subplots(1, 3, figsize=(8.4, 2.5))
    ms = {2: "o", 3: "s", 4: "^"}
    for L, g in dfB.groupby("L"):
        ax[0].scatter(g.ana_mean, g.emp_mean, s=13, marker=ms[L], alpha=.8, lw=0, label=f"L={L}")
        ax[1].scatter(g.ana_sd, g.emp_sd, s=13, marker=ms[L], alpha=.8, lw=0)
        ax[2].scatter(g.ana_z, g.emp_z, s=13, marker=ms[L], alpha=.8, lw=0)
    for a, t in zip(ax, ["null mean", "null SD", "z"]):
        lo = min(a.get_xlim()[0], a.get_ylim()[0]); hi = max(a.get_xlim()[1], a.get_ylim()[1])
        a.plot([lo, hi], [lo, hi], "k--", lw=.6)
        a.set_xlabel("HetNetEX-MD " + t); a.set_ylabel(f"XSwap (B={B_NET}) " + t)
    ax[0].set_xscale("log"); ax[0].set_yscale("log")
    ax[0].legend(frameon=False, fontsize=6)
    fig.tight_layout(); fig.savefig(OUT / "fig3_network.png"); plt.close(fig)

    # Figure 4
    fig, ax = plt.subplots(1, 2, figsize=(6.2, 2.4))
    ax[0].loglog(dfC.n_gene, dfC.rel_err, "o-", color="#d95f02")
    ref = dfC.rel_err.iloc[0] * dfC.n_gene.iloc[0] / dfC.n_gene
    ax[0].loglog(dfC.n_gene, ref, "k--", lw=.8, label=r"$O(1/n)$")
    ax[0].set_xlabel("number of gene nodes $n$")
    ax[0].set_ylabel("relative error of null mean")
    ax[0].legend(frameon=False, fontsize=6)
    infl = dfB.ana_sd / dfB.ana_sd_indep
    ax[1].boxplot([infl[dfB.L == L] for L in sorted(dfB.L.unique())],
                  labels=[f"L={L}" for L in sorted(dfB.L.unique())], widths=.55)
    ax[1].axhline(1, color="r", ls="--", lw=.7)
    ax[1].set_ylabel(r"SD inflation from $\vartheta$")
    fig.tight_layout(); fig.savefig(OUT / "fig4_scaling.png"); plt.close(fig)


if __name__ == "__main__":
    main()
