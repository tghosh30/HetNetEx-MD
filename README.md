# HetNetEX-MD

**Exact and asymptotic inference for Multi-DWPC — removing the resampling null from gene-set connectivity search.**

[Multi-DWPC](https://github.com/lagillenwater/multi-dwpc) scores how strongly a *set* of source genes converges on a shared phenotype along a metapath of a heterogeneous knowledge graph, and tests that score against two degree-aware null models built by resampling.

This package replaces the resampling stage with closed-form calculations.

| | Multi-DWPC | HetNetEX-MD |
|---|---|---|
| Gene-set resampling null (N1) | `B` replicates | **exact**, one `O(N)` pass |
| Network permutation null (N2) | `B` XSwap rewirings | asymptotic, `O(Ln)` |
| Median aggregate | `B` replicates | **exact** null CDF, distribution-free |
| Smallest reportable *p* | `1/(B+1)` | unbounded |
| Reproducibility | seed-dependent | deterministic |
| Measured speed-up | — | 195× (N1), 967× (N2), 126× (median) |

It computes the *same* null hypothesis by a different route. It is not a different test.

---

## Install

```bash
git clone https://github.com/tghosh30/hetnetex-md.git
cd hetnetex-md
pip install -e ".[dev]"
pytest -q          # 31 tests, ~1 s
```

Core dependencies are only `numpy` and `scipy`. The benchmark scripts additionally need `pandas`, `matplotlib`, and `statsmodels` (`pip install -e ".[benchmarks]"`).

---

## Quickstart

```python
import numpy as np
from hetnetex_md import hetnetex_md_resampling

# scores[i] = arcsinh(DWPC(gene_i -> target) / c_metapath), computed once
# bin_ids[i] = degree-stratum label for gene i
res = hetnetex_md_resampling(scores, gene_idx, bin_ids)

res["null_mean"]     # exact, not estimated
res["null_sd"]       # exact
res["null_skew"]     # exact third moment -- free Edgeworth correction
res["z"]
res["p_edgeworth"]   # no floor
```

A complete runnable version, including a side-by-side comparison against `B=1000` resampling, is in [`examples/quickstart.py`](examples/quickstart.py).

---

## What the functions correspond to

| Function | Implements | Status |
|---|---|---|
| `exact_resampling_moments` | Lemma 1 / Theorem 2 — stratified SRSWOR moments | **exact** (identity) |
| `edgeworth_upper_tail` | Theorem 4 — Edgeworth-corrected tail | asymptotic, `O(K⁻¹)` |
| `exact_median_pvalue` | Theorem 5b — hypergeometric convolution | **exact**, distribution-free |
| `network_null_moments` | Lemma 3, Theorems 6–7 — configuration-model moments | asymptotic, `O(1/n)` |
| `fit_soft_cm` / `soft_cm_ratio` | Lemma 4b — max-entropy edge probability | asymptotic |
| `correction_factors` | Lemma 4 — mean-field simple-graph correction | superseded by `fit_soft_cm` |
| `admissible_self_kernels` | admissibility rule for same-source kernels | exact rule |

The distinction in the last column is load-bearing. **The N1 results are identities**: `exact_resampling_moments` and `exact_median_pvalue` return what infinite resampling would converge to, and the test suite verifies this against exhaustive enumeration of every subset. **The N2 results are approximations** whose error is measured and reported, not assumed.

---

## Reproducibility

Two things were fixed specifically to make results reproducible, and both are worth knowing about:

1. **Seeding.** Benchmark replicate seeds are derived via BLAKE2b (`stable_seed`), not Python's built-in `hash()`. Python randomises string hashing per process (`PYTHONHASHSEED`), so `hash()`-derived seeds silently differ between runs. Verified stable across processes.
2. **Fit caching.** `soft_cm_ratio` keys its cache on the edge's degree sequence, not its name, so two networks sharing an edge-type name cannot collide. The earlier name-keyed version returned another network's fitted multipliers, a silent ~6% error.

The solver also **refuses to fail quietly**: `fit_soft_cm` raises on a degree sequence that admits no simple bipartite graph, and on failure to converge. Both conditions previously returned plausible-looking numbers. Pass `strict=False` to downgrade to a warning.

To reproduce the published benchmark numbers exactly, install from `requirements-lock.txt`. Ordinary use does not need pinned versions.

```bash
pip install -r requirements-lock.txt
python benchmarks/run_benchmarks.py     # writes CSVs + figures
```

**Timings are ratios, not absolutes.** All reported speed-ups were measured single-threaded in a Linux container and describe the *null stage only* — computing the observed DWPC costs the same under both methods.

---

## Known limitations

Stated plainly, because they affect how results should be read.

- **Synthetic networks.** All measurements are on synthetic hetnets with Hetionet-like layer structure. The N1 results do not depend on this (they are identities for any fixed score vector), but **every N2 number does**. Real Hetionet benchmarks via `hetmatpy` are the main outstanding task.
- **Variance overshoot at `L=2`.** The analytical null SD is exact at `L=4` (ratio 1.044) but roughly 20% too large at `L=2`. Two explanations were tested and refuted — the canonical/microcanonical ensemble gap, and double-counting in `Var(Y_g)`. The cause is open. The direction is conservative: it inflates *p*-values, costing power rather than creating false positives.
- **Canonical vs microcanonical.** `fit_soft_cm` is the expected-degree ensemble; XSwap samples the exact-degree one. The gap measured below the noise floor here, but Hetionet's hub degrees are far more extreme and this must be rechecked.
- **Transform under N2.** No accurate arcsinh-scale network-null moments at `L≥3`; use the raw-scale aggregate (`mean_dwpc`).
- **Aggregate choice.** Mean and median are covered exactly. Trimmed means and weighted quantiles are not — use resampling for those.
- **Chain metapaths.** Metapaths revisiting a node type need a self-avoidance correction, which is `O(1/n)` but not zero.
- **Extreme tails.** A `z` of 49 yields `p ≈ 1e-300`, an extrapolation of a limit theorem far past where it was validated. Report `z`, or truncate `p`, beyond roughly `z > 6`.

---

## Integrating with the Multi-DWPC pipeline

HetNetEX-MD emits the same schema the existing statistics assembly consumes (`null_mean`, `null_std`, `p_empirical`, `n_eff`, `d_median`), so downstream code needs no changes. Two columns lose their meaning and should **not** be used as quality filters in analytical mode:

- `n_eff` → `inf` (there are no replicates)
- `d_iqr` → `0` (deterministic; no spread across replicates)

If a spread diagnostic is wanted, report the exact null skewness `null_skew` instead.

---

## Repository layout

```
src/hetnetex_md/     core.py          the method; every public function
tests/               test_exactness.py  exhaustive-enumeration verification
benchmarks/          synthetic_hetnet.py  network builder, XSwap, DWPC
                     run_benchmarks.py    reproduces the paper's tables/figures
examples/            quickstart.py
```

---

## Citing

See [`CITATION.cff`](CITATION.cff). The accompanying manuscript describes the theory in full, including proofs of every lemma and theorem referenced above.

## License

BSD 3-Clause. See [`LICENSE`](LICENSE).

## Acknowledgments

Multi-DWPC and its reference implementation are the work of Lucas A. Gillenwater, with James C. Costello and Casey S. Greene. This package builds directly on that framework and on the XSwap permutation method of Himmelstein et al. Supported by NIH R01 HD109765.
