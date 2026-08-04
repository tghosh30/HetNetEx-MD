# Vignettes

Each vignette is a runnable script plus a write-up. Every number shown in the
markdown is real output from the script beside it — nothing is hand-written.

```bash
pip install -e .
python vignettes/v01_exact_resampling.py
```

| | Vignette | Script | Covers |
|---|---|---|---|
| 1 | [The exact resampling null](01-exact-resampling-null.md) | `v01_exact_resampling.py` | Lemma 1, Theorem 2 — verified against exhaustive enumeration |
| 2 | [The exact median null](02-exact-median-null.md) | `v02_exact_median.py` | Theorem 5b — hypergeometric convolution, distribution-free |
| 3 | [The network null](03-network-null.md) | `v03_network_null.py` | Lemma 3, Lemma 4b, Theorems 6–7 — and what is still open |
| 4 | [Why the $p$-values differ](04-why-p-values-differ.md) | `v04_why_pvalues_differ.py` | The floor, the FDR bound, and what does agree |

For the underlying statements and proofs, see the
[theory notes](../docs/theory.md).
