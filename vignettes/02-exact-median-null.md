# Vignette 2 — The exact median null

**Run it:** `python vignettes/v02_exact_median.py`

The moment route of Theorem 2 needs a *linear* aggregate. The median is not
linear, and it is tempting to conclude that resampling is unavoidable for it.

It is not. The median has its own exact null, and in one respect a **stronger**
one — it delivers the whole distribution function and needs no limit theorem at
all ([Theorem 5b](../docs/theory.md#theorem-5b)).

The identity: the sample median is $\ge t$ **iff** at least
$\lceil(K+1)/2\rceil$ of the sampled values are $\ge t$. That count is a sum of
*independent hypergeometrics*, one per stratum, so its distribution is an exact
convolution.

---

## 1. Verified against every possible sample

```
Enumerated all C(12,5) = 792 samples.

   threshold t         brute force         convolution    difference
      0.386696   1.000000000000000   1.000000000000000      0.00e+00
      1.276011   0.500000000000000   0.500000000000000      1.11e-16
      2.207246   0.045454545454545   0.045454545454545      0.00e+00
```

## 2. Stratified, also exact

```
   threshold t         brute force         convolution    difference
      0.386696   1.000000000000000   1.000000000000000      0.00e+00
      1.276011   0.500000000000000   0.500000000000000      1.11e-16
      2.207246   0.045454545454545   0.045454545454545      0.00e+00


## 2. Stratified, also exact

   threshold t         brute force         convolution    difference
      0.834309   1.000000000000000   1.000000000000000      0.00e+00
      1.662692   0.594285714285714   0.594285714285714      0.00e+00
      3.222365   0.028571428571429   0.028571428571429      3.47e-18
```

## 3. No floor, and no CLT

```
N = 20,000, K = 151, B = 1000

observed median      2.175119
exact p              8.528e-44      (no distributional assumption)
resampling p         9.990e-04      floor = 9.990e-04
time                 6.78 ms vs 203.81 ms  -> 30x
```

The exact calculation reports $p\approx8.5\times10^{-44}$ where resampling is
pinned at its floor of $9.99\times10^{-4}$. And unlike the mean route, this
number involves **no distributional assumption** — it is a finite combinatorial
sum, so the extreme-tail extrapolation problem that affects a CLT-based $p$-value
simply does not arise.

## 4. The caution: the median is often degenerate

```
fraction of the gene set with no path : 56%
observed median                       : 0.000000

On sparse metapaths more than half a gene set typically has no path to
the target, so the median is exactly zero and the test is vacuous. This
argues for the mean aggregate on scientific grounds, not computational ones.
```

This is worth knowing before choosing the median. In our 132-feature benchmark,
122 features had an observed median of exactly zero. The degeneracy is a property
of sparse metapaths, not of the method — but it means the mean aggregate is
usually the better scientific choice, and Theorem 5b matters when the median is
genuinely what you want.

---

**Next:** [Vignette 3 — the network null](03-network-null.md)
