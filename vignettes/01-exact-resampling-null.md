# Vignette 1 — The exact resampling null

**Run it:** `python vignettes/v01_exact_resampling.py`

Multi-DWPC's gene-set null holds the graph fixed and redraws the *set*. Because
the graph never moves, every gene's score is a fixed number, and the only
randomness is which $K$ of $N$ fixed numbers get averaged. That is a
finite-population sampling problem with exact closed-form moments
([Theorem 2](../docs/theory.md#lemma-1--theorem-2)).

---

## 1. The formulas are identities, not approximations

The cleanest way to check a closed form is to enumerate *every* possible sample
and compare. For small $(N,k)$ that is feasible.

```
Enumerated all C(12,5) = 792 samples.

                         brute force           closed form    difference
mean               2.347850952073050     2.347850952073050      0.00e+00
variance           0.148776724875777     0.148776724875777      2.78e-17
third moment       0.003321046122407     0.003321046122407      3.25e-17
```

Agreement to fifteen decimal places. The closed form is not an approximation of
resampling — it is what resampling converges to.

## 2. Stratification

Multi-DWPC bins genes by degree and samples within bins, so the real statistic is
a *stratified* SRSWOR mean. Strata are independent, so means add, variances add
with weights $(K_r/K)^2$, and third cumulants add with weights $(K_r/K)^3$.

```
Two strata: C(7,3) x C(6,2) = 525 samples.

                         brute force           closed form    difference
mean               1.658853713597177     1.658853713597177      0.00e+00
variance           0.169299349105733     0.169299349105733      8.33e-17
third moment       0.002217686761734     0.002217686761734      8.24e-17
```

## 3. At realistic scale

Hetionet has roughly 21,000 genes. Here is a query of 200 genes against a
degree-stratified null, compared with $B=1000$ resampling.

```
N = 21,000 genes, K = 200, B = 1000

                           exact      resampling   rel. diff
null mean               0.268382        0.264332      1.51%
null SD                 0.041733        0.042511      1.86%
skewness                0.264046   (unavailable)

z                          8.109
p-value                6.239e-15       9.990e-04   <- floor 9.990e-04
time                      7.68ms        156.14ms         20x
```

Three things to notice.

**The Monte-Carlo moments are converging to the exact ones**, not to something
else — the 1.5% and 1.9% gaps shrink as $B^-0.5$.

**The skewness is simply unavailable to resampling.** Third moments converge
slowly, so a Monte-Carlo estimate would be too noisy to use. Here it is exact and
free, which is what makes the Edgeworth correction worth applying.

**The $p$-values differ by orders of magnitude** — and this is not an error.
The resampling $p$ has hit its floor of $1/(B+1)$: it cannot report anything
smaller, however strong the evidence. See
[vignette 4](04-why-p-values-differ.md).

---

**Next:** [Vignette 2 — the exact median null](02-exact-median-null.md)
