# Vignette 3 — The network-permutation null

**Run it:** `python vignettes/v03_network_null.py`

Under the network null the *graph* moves and the gene set stays fixed. This is
the harder case: the results are asymptotic rather than exact, and the residual
error is measured rather than assumed.

---

## 1. The null mean is rank one

[Lemma 3](../docs/theory.md#lemma-3) says the configuration-model null mean
factorises into a metapath constant, a source factor, and a target factor.

```
mu has shape (1500, 110) = 165,000 entries
numerical rank: 1

So one O(Ln) pass gives every gene against every target. The whole
matrix is reconstructed from a scalar C and two vectors a, b:
  C            = 4.839282e-04
  a (sources)  = 1,500 values
  b (targets)  = 110 values
  stored       = 1,611 numbers instead of 165,000
```

This is the structural reason the method is fast: **one $O(Ln)$ pass gives every
gene against every target**, which is exactly what a many-to-one query needs. A
permutation pipeline cannot amortise this way, because every replicate is a
different graph and the matrix chain must be walked again.

<a name="2-the-edge-probability-ladder"></a>
## 2. The edge-probability ladder

Three models for $p_{uv}$, in increasing order of correctness
([Lemma 4b](../docs/theory.md#lemma-4b)):

```
max-entropy fit: median relative degree residual = 1.56e-12

Compare the three models on representative degree pairs:

   d_u   d_v    sparse d_u d_v/m     1 - exp(-x)     max-entropy
     2     3            0.001333        0.001332        0.000918
    10    20            0.044444        0.043471        0.048232
    40    60            0.533333        0.413354        0.555946
   100    90            2.000000        0.864665        0.683948
   200   150            6.666667        0.998727        0.823117

The three agree when d_u d_v / m is small and separate as it grows.
Hub-to-hub pairs are exactly where they differ -- and where the
sparse limit can exceed 1, which is not a probability.
```

Look at the last two rows. **The sparse limit returns 2.0 and 6.67** — those are
not probabilities. The exponential form saturates near 1, which is better but
still wrong. Only the max-entropy form respects the degree constraints, and it is
fitted here to a relative residual of $1.6\times10^{-12}$.

This is not a cosmetic difference. Switching from the sparse limit to
max-entropy cut the median error of the null mean from 22.9% to 4.7% and moved
the bias ratio from 0.771 to 1.004 — at which point the residual equals the
Monte-Carlo noise of the reference itself.

> `fit_soft_cm` **raises** on a degree sequence that admits no simple bipartite
> graph, and on failure to converge. Both conditions previously returned
> plausible-looking garbage. Real hetnets with extreme hubs are exactly where
> this matters.

## 3. Inter-source dependence is not a correction

Genes in a set share paths and all end at the same target, so their scores rise
and fall together ([Theorem 6](../docs/theory.md#theorems-67--inter-source-dependence-and-the-aggregate-variance)).

```
theta (one scalar per metapath+target) = 1.110607e+00

 gene set size K    independent SD      with theta     ratio
              10      2.931001e-03    3.115371e-03      1.06
              50      1.381125e-03    1.842150e-03      1.33
             100      1.008244e-03    1.680667e-03      1.67
             200      6.940542e-04    1.422911e-03      2.05
             400      4.873191e-04    1.358340e-03      2.79

The ratio grows with K: assuming genes are independent understates the
null spread more and more as the gene set gets larger.
```

**The ratio grows with the gene-set size.** At $K=10$ ignoring dependence costs
6%; at $K=400$ it understates the null spread by nearly threefold. Multi-DWPC
runs on sets of tens to hundreds of genes, so this is a leading effect — and a
correction to standard practice that holds whether or not you adopt the
analytical null.

## 4. The admissibility rule

```
   L      admissible same-source kernels
   2                                none   value = 0.00
   3                   theta_1 + theta_L   value = 0.44
   4       theta_1 + theta_mid + theta_L   value = 0.66
   5                                 all   value = 0.66

Forcing layer l shared pins interior positions l-1 and l. If no interior
position is left free the two paths coincide, so the term must not be
counted. At L=2 nothing is admissible and Var(Y_g) = mu^(2) exactly.
```

An exact combinatorial rule, though numerically small: it moved the $L=2$
standard-deviation ratio from 0.786 to 0.802.

---

## What is still open

The analytical standard deviation is exact at $L=4$ (ratio 1.044) but roughly
20% too large at $L=2$. Two explanations were tested and **refuted** — the
canonical/microcanonical ensemble gap, and double-counting in
$\operatorname{Var}(Y_g)$. The cause is not yet identified. The direction is
conservative: it inflates $p$-values, costing power rather than creating false
positives.

---

**Next:** [Vignette 4 — why the $p$-values differ](04-why-p-values-differ.md)
