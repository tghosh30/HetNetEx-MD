# Vignette 4 — Why the two $p$-values cannot agree

**Run it:** `python vignettes/v04_why_pvalues_differ.py`

This is the first question anyone asks, and the answer is not "they nearly do."

**They cannot agree, and expecting them to misunderstands what each one is.**
A permutation $p$-value is *measured by counting*; an analytical $p$-value is
*computed from a formula*. They target the same quantity by different routes.

---

## 1. A permutation $p$-value is a random variable

```
True tail probability p* = 0.002

Twenty independent runs of the SAME analysis, B = 1000:

  0.0040  0.0030  0.0050  0.0040  0.0010  0.0060  0.0040  0.0040  0.0010  0.0030
  0.0020  0.0050  0.0030  0.0040  0.0030  0.0020  0.0030  0.0010  0.0040  0.0030

  range      0.0010 to 0.0060
  mean       0.00325   (theory: 0.00300)
  SD         0.00137   (theory: 0.00141)

Same data, same code, different seed -> different p-value.
The analytical p-value is deterministic: one number, every time.
```

Twenty runs of identical code on identical data, differing only in seed, spread
over a sixfold range. Both the conservative bias and the standard deviation match
theory. The analytical $p$-value has neither — it is the same number every time.

## 2. Relative error blows up exactly where it matters

```
          p*     rel. std. error     reportable at B=1000?
       5e-01               3.2%                       yes
       1e-01               9.5%                       yes
       1e-02              31.5%                       yes
       1e-03              99.9%                       yes
       1e-04             316.2%      NO - below the floor
       1e-06            3162.3%      NO - below the floor

The smaller the true p-value -- i.e. the stronger the finding -- the
worse resampling estimates it, until it cannot represent it at all.
```

This is the uncomfortable part. The relative standard error of a permutation
$p$-value is $\sqrt{(1-p^*)/(Bp^*)}$, which **diverges as $p^*\to0$**. The
stronger your finding, the worse resampling estimates it — until, below
$1/(B+1)$, it cannot represent it at all.

## 3. The floor meets multiple-testing correction

```
BH needs p_(k) <= alpha*k/m, but the floor gives p >= 1/(B+1).
So NO discovery is possible unless  B >= m/alpha - 1.

   metapaths m    B required      B=200?     B=1000?
             5            99          ok          ok
            11           219     BLOCKED          ok
            25           499     BLOCKED          ok
            50           999     BLOCKED          ok
           100          1999     BLOCKED     BLOCKED
           200          3999     BLOCKED     BLOCKED
           500          9999     BLOCKED     BLOCKED
```

Benjamini–Hochberg demands a *smaller* $p$-value as $m$ grows; the floor sets a
*minimum*. Push them together and there is a hard requirement on $B$ before any
discovery is arithmetically possible.

Note the practical consequence: with 100 metapaths — routine for a
connectivity-search query — even $B=1000$ is blocked.

## 4. What this looks like on real features

```
11 metapaths, B = 200, FDR alpha = 0.05
BH threshold for the top feature = alpha/m = 0.00455
Permutation floor                = 1/(B+1) = 0.00498
-> the floor is ABOVE the threshold, so the best possible permutation p-value CANNOT clear it.

  feature         z    analytic p    best possible perm p
       10     14.94      1.34e-48                 0.00498
        9     10.40      7.12e-24                 0.00498
        8      8.86      1.39e-17                 0.00498
        3      0.82      2.03e-01                 0.00498
        7      0.77      2.16e-01                 0.00498
        6      0.10      4.41e-01                 0.00498

analytic, clearing the strictest BH bar : 3 of 11
permutation at B=200, best case          : 0 of 11 (floor exceeds the bar)
```

Three features carry overwhelming evidence ($z$ of 14.9, 10.4, 8.9). The
analytical method reports them. The permutation method reports **nothing** — not
because it disagrees, but because its best possible output, 0.00498, sits above
the bar of 0.00455.

This is a *resolution* failure, not a statistical one.

## 5. What does agree: the ranking

```
Both p-values are monotone in the same z. The ranking transfers even
though the numbers never match:

  ranking by z      : [10, 9, 8, 3, 7, 6, 4, 5, 0, 2, 1]
  ranking by p      : [10, 9, 8, 3, 7, 6, 4, 5, 0, 2, 1]
  identical         : True

The scientific decision is 'which metapaths are most promising',
not 'is this p-value 0.0031 or 0.0028'. That decision is preserved.
```

Both $p$-values are monotone in the same standardised statistic $z$, so the
ordering is preserved even though the numbers never coincide. In the
132-feature benchmark, Spearman $\rho$ was 0.997–0.999 on $z$, and **zero
features** had matching $p$-values.

---

## The bottom line

| | Permutation | Analytical |
|---|---|---|
| Nature | empirical — counts exceedances | distributional — evaluates a tail |
| Support | $\{1/(B+1),\dots,1\}$ | $(0,1)$ |
| Bias | $+(1-p^*)/(1+B)$ | $O(K^{-1})$ (mean); **none** (median) |
| Variance | $\approx p^*(1-p^*)/B$ | none |
| Reproducible | no — seed-dependent | yes — bit-for-bit |
| Assumes a limit law | no | yes for the mean, **no** for the median |

The last row is the honest cost, and it is why
[Theorem 5b](../docs/theory.md#theorem-5b) matters: for the median aggregate,
even that cost disappears.

---

**Back to:** [Vignette 1](01-exact-resampling-null.md) ·
[Theory notes](../docs/theory.md)
