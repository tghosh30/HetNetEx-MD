# Theory notes

Complete statements and proofs for HetNetEX-MD. Every public function implements
a specific result below, and each is labelled **exact** or **asymptotic** — that
distinction is the honest core of the method and is never blurred.

**Contents**

- [Setting and notation](#setting-and-notation) · [Assumptions](#assumptions)
- **Part I — exact theory for N1**: [Lemma 1](#lemma-1) · [Theorem 2](#theorem-2) · [Corollary 3](#corollary-3) · [Theorem 4](#theorem-4) · [Theorem 5](#theorem-5) · [Theorem 5b](#theorem-5b)
- **Part II — asymptotic theory for N2**: [Lemma 3](#lemma-3) · [Lemma 4](#lemma-4) · [Lemma 4b](#lemma-4b) · [Theorem 6](#theorem-6) · [Theorem 7](#theorem-7) · [Corollary 7b](#corollary-7b) · [Admissibility](#admissibility) · [Theorem 8](#theorem-8) · [Theorem 9](#theorem-9) · [Theorem 10](#theorem-10)
- **Part III — the two $p$-values**: [Prop 11](#proposition-11) · [Prop 12](#proposition-12) · [Theorem 13](#theorem-13) · [Prop 14](#proposition-14) · [Prop 15](#proposition-15)
- [Results and figures](results.md) · [Limitations](#limitations-and-open-problems)

| Function | Implements | Status |
|---|---|---|
| `exact_resampling_moments` | [Lemma 1](#lemma-1) / [Theorem 2](#theorem-2) | **exact (identity)** |
| `edgeworth_upper_tail` | [Theorem 4](#theorem-4) | asymptotic, $O(K^{-1})$ |
| `exact_median_pvalue` | [Theorem 5b](#theorem-5b) | **exact, distribution-free** |
| `network_null_moments` | [Lemma 3](#lemma-3), [Theorems 6–7](#theorem-6) | asymptotic, $O(1/n)$ |
| `fit_soft_cm` / `soft_cm_ratio` | [Lemma 4b](#lemma-4b) | asymptotic |
| `correction_factors` | [Lemma 4](#lemma-4) | superseded by `fit_soft_cm` |
| `admissible_self_kernels` | [admissibility rule](#admissibility) | **exact rule** |

---

## Setting and notation

A metapath $\mathcal{P}=(e_1,\dots,e_L)$ visits node types $V_0,\dots,V_L$.
Degrees are **metaedge-specific**: $d^{(e)}_v$ counts only edges of type $e$, and
$m_e$ is the edge count of type $e$. For a path $\pi=(v_0,\dots,v_L)$ the DWPC
weight damps every edge endpoint,

$$w_\pi=\prod_{\ell=1}^{L}\bigl(d^{(e_\ell)}_{v_{\ell-1}}\bigr)^{-w}\bigl(d^{(e_\ell)}_{v_\ell}\bigr)^{-w},$$

so an **interior** node contributes $d^{-2w}$ in total. Then

$$Y_g=\mathrm{DWPC}(g,t;\mathcal{P})=\sum_{\pi:g\rightsquigarrow t}w_\pi I_\pi,
\qquad x_g=\operatorname{arcsinh}(Y_g/c_\mathcal{P}),
\qquad T(S)=\frac{1}{K}\sum_{g\in S}x_g.$$

Write $N$ for genes in the universe, $K=|S|$, $R$ for degree strata, $L$ for path
length, $n$ for nodes of the relevant type, $B$ for replicates. Set

$$a^{(r)}_v=\bigl(d^{(e_1)}_v\bigr)^{1-rw},\quad
b^{(r)}_v=\bigl(d^{(e_L)}_v\bigr)^{1-rw},\quad
\psi^{(r)}_j(v)=\bigl(d^{(e_j)}_v\bigr)^{1-rw}\bigl(d^{(e_{j+1})}_v\bigr)^{1-rw},$$

with $S^{(r)}_j=\sum_{v\in V_j}\psi^{(r)}_j(v)$. Abbreviate $a_v=a^{(1)}_v$,
$\psi_j=\psi^{(1)}_j$, $S_j=S^{(1)}_j$.

Multi-DWPC tests $T$ against two nulls:

- **N1 — gene-set resampling.** Graph fixed; the *set* is redrawn by stratified
  sampling without replacement from degree-binned pools.
- **N2 — network permutation.** Set fixed; the *graph* is rewired by XSwap,
  preserving every degree.

These are mathematically very different situations, which is why Part I is exact
and Part II is not.

### Assumptions

**A1 (Configuration model).** For each edge type, the network null draws uniformly
from simple bipartite graphs with the observed degree sequence. This is the null
XSwap targets.

**A2 (Edge-type independence).** Wirings of distinct edge types are independent
given node identities. Hetionet's edge types come from distinct resources; the
assumption concerns wiring, not node overlap.

**A3 (Sparsity).** $d^{(e)}_{\max}\ll\sqrt{m_e}$, so
$\Pr(u\sim v)=1-e^{-d_ud_v/m}+O(1/n)=d_ud_v/m+O((d_ud_v/m)^2)$.

**A4 (Fixed transform).** $c_\mathcal{P}$ is computed once from the observed graph
and held fixed across replicates, so $\psi$ is a fixed strictly increasing map.

**A5 (Exchangeability within strata).** Under N1 the null set is a stratified
simple random sample without replacement. **This is not an approximation** — it
restates what the reference implementation's sampler does.

---

# Part I — exact theory for N1

Everything here rests on one observation:

> **Under N1, the graph does not move.**

Each $x_g$ was computed once from the observed graph and does not change while the
set is reshuffled. The only randomness is *which $K$ of $N$ fixed numbers get
averaged* — a finite-population sampling problem with a complete moment theory.

<a name="lemma-1"></a>
## Lemma 1 — SRSWOR moments of the sample mean

Let $x_1,\dots,x_N$ be fixed reals with mean $\bar X$ and central moments
$m_q=N^{-1}\sum_i(x_i-\bar X)^q$. Let $\bar x_k$ be the mean of a simple random
sample of size $1\le k\le N$ drawn **without replacement**. Then

$$\mathbb{E}[\bar x_k]=\bar X\tag{1}$$

$$\operatorname{Var}(\bar x_k)=\frac{m_2}{k}\cdot\frac{N}{N-1}\Bigl(1-\frac{k}{N}\Bigr)\tag{2}$$

$$\mathbb{E}\bigl[(\bar x_k-\bar X)^3\bigr]=\frac{(N-k)(N-2k)}{k^2(N-1)(N-2)}\,m_3,\qquad N\ge3.\tag{3}$$

<details>
<summary><b>Proof</b></summary>

Put $\delta_i=x_i-\bar X$, so $\sum_i\delta_i=0$, and let $\mathcal{S}$ be a
uniform random $k$-subset with $D=\sum_{i\in\mathcal{S}}\delta_i$, so that
$\bar x_k-\bar X=D/k$. The SRSWOR inclusion probabilities are

$$\pi_i=\frac{k}{N},\qquad \pi_{ij}=\frac{k(k-1)}{N(N-1)},\qquad \pi_{ijl}=\frac{k(k-1)(k-2)}{N(N-1)(N-2)}$$

for distinct $i,j,l$.

**First moment.** $\mathbb{E}[D]=\sum_i\pi_i\delta_i=(k/N)\sum_i\delta_i=0$, giving (1).

**Second moment.**

$$\mathbb{E}[D^2]=\sum_i\pi_i\delta_i^2+\sum_{i\ne j}\pi_{ij}\delta_i\delta_j=\frac{k}{N}\sum_i\delta_i^2+\frac{k(k-1)}{N(N-1)}\Bigl[\Bigl(\sum_i\delta_i\Bigr)^2-\sum_i\delta_i^2\Bigr]$$

$$=k\,m_2-\frac{k(k-1)}{N-1}m_2=k\,m_2\,\frac{N-k}{N-1}.$$

Dividing by $k^2$ gives (2).

**Third moment.** Grouping by the number of distinct indices,

$$\mathbb{E}[D^3]=\sum_i\pi_i\delta_i^3+3\sum_{i\ne j}\pi_{ij}\delta_i^2\delta_j+\sum_{i\ne j\ne l}\pi_{ijl}\delta_i\delta_j\delta_l.$$

Using $\sum_i\delta_i=0$ we have $\sum_{i\ne j}\delta_i^2\delta_j=-\sum_i\delta_i^3=-Nm_3$
and $\sum_{i\ne j\ne l}\delta_i\delta_j\delta_l=(\sum\delta)^3-3(\sum\delta)(\sum\delta^2)+2\sum\delta^3=2Nm_3$.
Therefore

$$\mathbb{E}[D^3]=k\,m_3\Bigl[1-\frac{3(k-1)}{N-1}+\frac{2(k-1)(k-2)}{(N-1)(N-2)}\Bigr]=k\,m_3\,\frac{(N-k)(N-2k)}{(N-1)(N-2)},$$

the last step an algebraic identity verified by expanding both sides. Divide by
$k^3$. $\blacksquare$
</details>

Three consistency checks pin the constants: $k=1\Rightarrow m_3$;
$k=N\Rightarrow0$ (degenerate); $k=N/2\Rightarrow0$ by the symmetry
$\bar x_k-\bar X=-(\bar x_{N-k}-\bar X)$.

<a name="theorem-2"></a>
## Theorem 2 — exact moments of the Multi-DWPC statistic

Fix a target $t$ and metapath $\mathcal{P}$, and let $x_g$ be **any** per-gene
score computed from the observed graph. Let the null set be drawn by independent
SRSWOR of $K_r$ genes from pool $\mathcal{U}_r$, $r=1,\dots,R$, with
$K=\sum_rK_r$ and $N_r=|\mathcal{U}_r|$. Then $\tilde T=K^{-1}\sum_{g\in\tilde S}x_g$
satisfies, **exactly and with no asymptotic approximation**,

$$\mathbb{E}[\tilde T]=\sum_{r=1}^{R}\frac{K_r}{K}\bar X_r\tag{4}$$

$$\operatorname{Var}(\tilde T)=\sum_{r=1}^{R}\Bigl(\frac{K_r}{K}\Bigr)^{2}\frac{m_{2,r}}{K_r}\cdot\frac{N_r}{N_r-1}\Bigl(1-\frac{K_r}{N_r}\Bigr)\tag{5}$$

$$\mathbb{E}\bigl[(\tilde T-\mathbb{E}\tilde T)^3\bigr]=\sum_{r=1}^{R}\Bigl(\frac{K_r}{K}\Bigr)^{3}\frac{(N_r-K_r)(N_r-2K_r)}{K_r^{2}(N_r-1)(N_r-2)}m_{3,r}\tag{6}$$

All three are computable in $O(\sum_rN_r)\le O(N)$ time and $O(1)$ extra memory per
feature.

<details>
<summary><b>Proof</b></summary>

Write $\tilde T=\sum_r(K_r/K)\,\bar x^{(r)}_{K_r}$ where $\bar x^{(r)}_{K_r}$ is
the SRSWOR sample mean within stratum $r$. Strata are sampled independently, so
means add; variances add with weights $(K_r/K)^2$; and third **cumulants** — which
coincide with third central moments — add with weights $(K_r/K)^3$. Apply Lemma 1
within each stratum. $\blacksquare$
</details>

```python
mean, var, mu3 = exact_resampling_moments(scores, pools, counts)
```

<a name="corollary-3"></a>
## Corollary 3 — transform invariance

Theorem 2 holds verbatim for $x_g=\psi(Y_g)$ with **any** fixed $\psi$ (A4): the
arcsinh transform, raw DWPC, ranks, loading-weighted scores, anything.

This is precisely why N1 admits an exact treatment while N2 does not. **Under N1
the transform is applied to constants, so Jensen's inequality never bites; under
N2 it is applied to a random variable and it does.** Weighted gene sets therefore
follow immediately — put the weights into $x_g$.

<a name="theorem-4"></a>
## Theorem 4 — Edgeworth-corrected tail

Let $\gamma_1=\mathbb{E}[(\tilde T-\mathbb{E}\tilde T)^3]/\operatorname{Var}(\tilde T)^{3/2}$.
Under the conditions of [Theorem 5](#theorem-5) with finite third moment,
uniformly on compact $z$-sets,

$$\Pr(Z\le z)=\Phi(z)-\phi(z)\frac{\gamma_1}{6}(z^2-1)+O(K^{-1}),$$

so the upper-tail $p$-value is

$$p_{\mathrm{HX}}=1-\Phi(z)+\phi(z)\frac{\gamma_1}{6}(z^2-1)+O(K^{-1}).\tag{7}$$

<details>
<summary><b>Proof sketch</b></summary>

The cumulants of $Z$ satisfy $\kappa_1=0$, $\kappa_2=1$,
$\kappa_3=\gamma_1=O(K^{-1/2})$, $\kappa_4=O(K^{-1})$ under the Lindeberg
condition. The formal Edgeworth series for a lattice-free standardised sum is
$F_Z(z)=\Phi(z)-\phi(z)\frac{\gamma_1}{6}He_2(z)+O(K^{-1})$ with $He_2(z)=z^2-1$.
Validity for sampling without replacement, including the required Cramér-type
condition, is established by Robinson (1978) and Babu & Singh (1985); the
stratified case follows by applying the expansion to the independent stratum sum.
Since $\gamma_1$ here is computed **exactly** from (6) rather than estimated, no
additional estimation error enters the correction term. $\blacksquare$
</details>

This is the one asymptotic step in Part I, and the only reason a Part I $p$-value
is not exact even though its moments are. Measured: median $|\gamma_1|=0.28$,
median $|\log_{10}(p_{\mathrm{normal}}/p_{\mathrm{Edgeworth}})|=0.0082$.

Resampling cannot match this — third moments converge slowly, so a Monte-Carlo
$\gamma_1$ would be too noisy to use.

```python
z, p = edgeworth_upper_tail(t_obs, mean, var, mu3)
```

<a name="theorem-5"></a>
## Theorem 5 — finite-population CLT (Hájek)

Suppose $\{x_g\}$ satisfies the Hájek–Lindeberg condition

$$\frac{1}{N_rm_{2,r}}\sum_{g\in\mathcal{U}_r}(x_g-\bar X_r)^2\,\mathbf{1}\bigl\{|x_g-\bar X_r|>\varepsilon\sqrt{N_rm_{2,r}}\bigr\}\longrightarrow0\quad\forall\varepsilon>0$$

for each stratum with $K_r\to\infty$ and $K_r/N_r\to f_r\in[0,1)$. Then

$$Z=\frac{\tilde T-\mathbb{E}[\tilde T]}{\sqrt{\operatorname{Var}(\tilde T)}}\xrightarrow{d}\mathcal{N}(0,1).$$

<details>
<summary><b>Proof</b></summary>

Hájek's CLT for SRSWOR applies within each stratum; the strata are independent,
and a finite sum of independent asymptotically normal terms with non-degenerate
variance ratios is asymptotically normal. $\blacksquare$
</details>

<a name="theorem-5b"></a>
## Theorem 5b — exact null distribution of the stratified median

The moment route requires a linear aggregate. The median is not linear — but it
has its own exact null, and in one respect a **stronger** one: it delivers the
entire distribution function and needs no limit theorem at all.

Let $\tilde M=\operatorname{median}\{x_g:g\in\tilde S\}$. For any threshold $t$,
put $c_r(t)=|\{g\in\mathcal{U}_r:x_g\ge t\}|$. Then

$$\Pr\bigl(\tilde M\ge t\bigr)=\Pr\Bigl(\sum_{r=1}^{R}A_r\ge\bigl\lceil(K{+}1)/2\bigr\rceil\Bigr),\qquad A_r\sim\mathrm{Hypergeom}\bigl(N_r,c_r(t),K_r\bigr)\ \text{independent},\tag{8}$$

and the right-hand side is the exact convolution of $R$ hypergeometric pmfs,
computable in $O(K\sum_rK_r)$ directly or $O(K\log K)$ by FFT. Setting
$t=M_{\mathrm{obs}}$ gives the exact upper-tail $p$-value.

<details>
<summary><b>Proof</b></summary>

The sample median is $\ge t$ if and only if at least $\lceil(K+1)/2\rceil$ of the
$K$ sampled values are $\ge t$ (for even $K$, at least $K/2+1$ under the
lower-median convention). Let $A_r$ count sampled genes in stratum $r$ with score
$\ge t$. Within stratum $r$ we draw $K_r$ of $N_r$ items without replacement, of
which $c_r(t)$ are "successes", so $A_r\sim\mathrm{Hypergeom}(N_r,c_r(t),K_r)$
exactly. Strata are independent, so the law of $\sum_rA_r$ is the convolution.
$\blacksquare$
</details>

**No CLT, no Edgeworth, no distributional assumption whatever.** Exact in finite
samples, exact in the far tail, no $p$-value floor. This sidesteps the extreme-tail
extrapolation problem entirely.

> **Caution.** On sparse metapaths more than half a typical gene set has no path to
> the target, so the observed median is exactly zero and the test is vacuous. In
> our benchmark this happened for 122 of 132 features. Prefer the mean aggregate on
> *scientific* grounds; Theorem 5b matters when the median is genuinely the
> statistic of interest.

```python
p = exact_median_pvalue(scores, pools, counts, t_obs)
```

---

# Part II — asymptotic theory for N2

Now the graph moves. Both the marginal law of each $Y_g$ **and the dependence
between $Y_g$ and $Y_h$** enter. Everything is leading order in $1/n$, and the
residual error is measured rather than assumed.

<a name="lemma-3"></a>
## Lemma 3 — chain factorisation of the null mean

Under A1–A3, for a chain metapath with distinct interior node types and any
$r\ge0$,

$$A_r(g,t):=\sum_{\pi:g\rightsquigarrow t}w_\pi^{\,r}\,\mathbb{E}[I_\pi]=\underbrace{C^{(r)}_\mathcal{P}}_{\text{metapath}}\underbrace{a^{(r)}_g}_{\text{source}}\underbrace{b^{(r)}_t}_{\text{target}}+O(1/n),\qquad C^{(r)}_\mathcal{P}=\Bigl(\prod_{\ell=1}^{L}\frac{1}{m_\ell}\Bigr)\prod_{j=1}^{L-1}S^{(r)}_j.\tag{9}$$

In particular the null mean is $\mu_{g,t}=A_1(g,t)$, the expected path *count* is
$\Lambda_{g,t}=A_0(g,t)$, and the weighted second moment is
$\mu^{(2)}_{g,t}=A_2(g,t)$. Each is computable for all sources and all targets
simultaneously in $O(Ln)$.

<details>
<summary><b>Proof</b></summary>

$\mathbb{E}[I_\pi]=\prod_\ell d^{(e_\ell)}_{v_{\ell-1}}d^{(e_\ell)}_{v_\ell}/m_\ell$
by A1–A3 and asymptotic independence of distinct node pairs. Multiplying by
$w_\pi^r$ converts the exponent of every endpoint from $1$ to $1-rw$. Every
interior node $v_j$ appears once as the destination of $e_j$ and once as the source
of $e_{j+1}$, so its total contribution is $\psi^{(r)}_j(v_j)$, and the sum over
interior nodes separates into $\prod_jS^{(r)}_j$. The endpoints $g$ and $t$ appear
once each, contributing $a^{(r)}_g$ and $b^{(r)}_t$. $\blacksquare$
</details>

**This is what makes the method fast.** The null mean is a *rank-one* object:
source effect × target effect × metapath constant. One $O(Ln)$ pass yields the null
mean of every gene against every target — exactly what a many-to-one query needs. A
permutation pipeline cannot amortise this way, because every replicate is a
different graph and the matrix chain must be walked again.

<a name="lemma-4"></a>
## Lemma 4 — mean-field simple-graph correction *(superseded)*

Because XSwap generates **simple** graphs, the edge probability is
$\tilde p_{uv}=1-e^{-d_ud_v/m}$ rather than $p_{uv}=d_ud_v/m$. Writing
$g(x)=(1-e^{-x})/x\le1$,

$$\tilde\mu_{g,t}=\mu_{g,t}\cdot\bar g_1(g)\Bigl(\prod_{\ell=2}^{L-1}\bar g_\ell\Bigr)\bar g_L(t)\bigl(1+O(\operatorname{Var}g)\bigr),\tag{10}$$

where $\bar g_\ell$ is the path-measure-weighted mean of $g(d_ud_v/m_\ell)$ over the
endpoints of layer $\ell$, computable in $O(n+b^2)$ with $b$ degree bins.

<details>
<summary><b>Proof</b></summary>

Writing $\mu=\sum_\pi\lambda_\pi$ and $\tilde\mu=\sum_\pi\lambda_\pi\prod_\ell g_\ell(\pi)$,
we have $\tilde\mu/\mu=\mathbb{E}_\lambda[\prod_\ell g_\ell]$, an expectation under
the normalised path measure $\lambda_\pi/\mu$. Under that measure the layers are
independent to leading order because $\lambda$ factorises (Lemma 3), so
$\mathbb{E}_\lambda[\prod_\ell g_\ell]=\prod_\ell\mathbb{E}_\lambda[g_\ell]+O(\mathrm{Cov})$.
Each $\mathbb{E}_\lambda[g_\ell]=\bar g_\ell$ is a double sum over the layer's
endpoint degrees weighted by $\psi_{\ell-1}$ and $\psi_\ell$, evaluated on $b$
quantile bins. For $\ell=1$ the source endpoint is fixed at $g$ and for $\ell=L$ the
target is fixed at $t$, giving the vector-valued factors. $\blacksquare$
</details>

`correction_factors` implements this. It is **retained for reproducibility of the
published comparison and superseded in practice** by Lemma 4b.

The mean-field step is not exact, because adjacent layers share an interior node. An
exact backward transfer-matrix pass in $O(Lb^2)$ moved the bias ratio from 0.850 to
0.889 but left the median error essentially unchanged at 15.4% — a real but minor
error next to the choice of $p_{uv}$ itself.

<a name="lemma-4b"></a>
## Lemma 4b — maximum-entropy simple-graph probability

For a bipartite edge type with degree sequences $\mathbf d^{\mathrm{row}}$,
$\mathbf d^{\mathrm{col}}$, the maximum-entropy distribution over simple bipartite
graphs subject to expected-degree constraints has independent edges with

$$p_{uv}=\frac{x_uy_v}{1+x_uy_v},\qquad\sum_vp_{uv}=d^{\mathrm{row}}_u,\quad\sum_up_{uv}=d^{\mathrm{col}}_v.\tag{11}$$

The multipliers are the unique positive solution, obtained by the fixed-point
iteration $x_u\leftarrow d^{\mathrm{row}}_u\big/\sum_vy_v/(1+x_uy_v)$ and
symmetrically for $y$, at $O(b^2)$ per iteration on $b$ degree bins. Moreover

$$p_{uv}=\frac{d_ud_v}{m}\Bigl(1+O\bigl(\tfrac{d_ud_v}{m}\bigr)\Bigr)=1-e^{-d_ud_v/m}+O\Bigl(\bigl(\tfrac{d_ud_v}{m}\bigr)^{3}\Bigr),$$

so the sparse-limit and exponential forms are its first- and second-order
approximations. Substituting (11) throughout Lemma 3 gives the corrected null mean,
again multiplicatively.

The three models agree when $d_ud_v/m$ is small and separate badly for hub pairs —
where the sparse limit can **exceed 1** and stop being a probability at all. See
[vignette 3](../vignettes/03-network-null.md#2-the-edge-probability-ladder) for a
worked table.

**Measured effect** across 66 network-null features:

| $p_{uv}$ model | median rel. error | bias ratio |
|---|---|---|
| sparse limit $d_ud_v/m$ | 22.9% | 0.771 |
| exponential $1-e^{-d_ud_v/m}$ | 15.5% | 0.850 |
| exact chain (transfer matrix) | 15.4% | 0.889 |
| **max-entropy** | **4.7%** | **1.004** |
| *Monte-Carlo noise floor of the reference* | *4.7%* | — |

The residual now **equals the sampling error of the $B=100$ reference itself**.

> Up to the canonical/microcanonical distinction, (11) is the object the XSwap edge
> prior (Zietz et al. 2024) estimates by sampling. Adopting it makes this a *reuse*
> of that result rather than a competing approximation to it.

```python
f_row, f_col, residual = fit_soft_cm(d_row, d_col)   # raises if infeasible
g = soft_cm_ratio(net, "GpP", du, dv, cache={})
```

`fit_soft_cm` **raises** on a degree sequence admitting no simple bipartite graph
and on failure to converge — both previously returned plausible-looking garbage.
Pass `strict=False` to downgrade to a warning.

<a name="theorem-6"></a>
## Theorem 6 — inter-source covariance

Two genes of the same set are not independent: their paths can share edges, and
both must terminate at the same target. For $g\ne h$, to leading order in the
number of shared edges,

$$\operatorname{Cov}(Y_g,Y_h)=\vartheta_t\,\mu_{g,t}\mu_{h,t}\bigl(1+O(1/n)\bigr),\qquad\vartheta_t=\sum_{\ell=2}^{L}\Theta_\ell,\tag{12}$$

where, with $\rho_\ell(u,x)=m_\ell/(d^{(e_\ell)}_ud^{(e_\ell)}_x)-1$,

$$\Theta_\ell=\frac{\sum_{u\in V_{\ell-1}}\sum_{x\in V_\ell}\psi_{\ell-1}(u)^2\psi_\ell(x)^2\rho_\ell(u,x)}{S_{\ell-1}^2S_\ell^2},\quad2\le\ell\le L-1;\qquad\Theta_L=\frac{\sum_{u\in V_{L-1}}\psi_{L-1}(u)^2\rho_L(u,t)}{S_{L-1}^2}.\tag{13}$$

**Crucially $\vartheta_t$ does not depend on $g$ or $h$**: to leading order the null
correlation structure across a gene set is exchangeable, a single scalar per
(metapath, target) rather than a $K\times K$ matrix. Both are $O(n)$ per layer via

$$\sum_{u,x}f(u)h(x)\rho_\ell(u,x)=m_\ell\Bigl(\sum_u\tfrac{f(u)}{d_u}\Bigr)\Bigl(\sum_x\tfrac{h(x)}{d_x}\Bigr)-\Bigl(\sum_uf(u)\Bigr)\Bigl(\sum_xh(x)\Bigr).\tag{14}$$

<details>
<summary><b>Proof</b></summary>

Let $\pi$ be a path from $g$ to $t$ and $\pi'$ from $h$ to $t$, and let
$F=F(\pi,\pi')$ be the set of edge slots where the two paths use the identical
edge. Under A1, distinct node pairs are asymptotically independent, so

$$\mathbb{E}[I_\pi I_{\pi'}]=\frac{\mathbb{E}[I_\pi]\mathbb{E}[I_{\pi'}]}{\prod_{e\in F}p_e},\qquad\operatorname{Cov}(I_\pi,I_{\pi'})=\mathbb{E}[I_\pi]\mathbb{E}[I_{\pi'}]\Bigl(\prod_{e\in F}p_e^{-1}-1\Bigr).$$

Since $g\ne h$, $1\notin F$. Retain $|F|=1$; a set of size $q$ carries $q$
additional degree constraints and contributes $O(n^{-(q-1)})$ relative to $|F|=1$.

Fix $F=\{\ell\}$ with $2\le\ell\le L-1$. Summing over all $(\pi,\pi')$ with
$v_{\ell-1}=v'_{\ell-1}=u$ and $v_\ell=v'_\ell=x$ and all other interior nodes
free, the separability of Lemma 3 applies at every position except $\ell-1$ and
$\ell$: at a free position $j$ the two paths each contribute $S_j$, giving $S_j^2$;
at the two forced positions they contribute $\psi_{\ell-1}(u)^2$ and
$\psi_\ell(x)^2$ instead of $S_{\ell-1}^2$ and $S_\ell^2$. Multiplying by
$\rho_\ell(u,x)=p_{ux}^{-1}-1$ and normalising by $\mu_{g,t}\mu_{h,t}$ — whose
expression contains $S_{\ell-1}^2S_\ell^2$ at these positions — yields the first
part of (13). The source factors $a_g,a_h$ and target factor $b_t$ cancel in the
normalisation, which is why $\Theta_\ell$ is free of $g$ and $h$.

For $F=\{L\}$ the destination is $t$ for both paths automatically, so only
$u=v_{L-1}$ is forced; the sum over $x$ does not appear and the normalisation
carries only $S_{L-1}^2$, giving the second part. Finally $\rho_\ell+1$ is itself
separable, which gives (14). $\blacksquare$
</details>

Note the sign: higher-order terms enter the expansion of
$\prod_{e\in F}p_e^{-1}-1$ with **positive** sign, so truncation cannot produce an
*over*estimate. This is relevant to the open problem below.

<a name="theorem-7"></a>
## Theorem 7 — aggregate null variance

With $\Theta_1(g)=S_1^{-2}\sum_{x\in V_1}\psi_1(x)^2\rho_1(g,x)$ the same-source
layer-1 kernel,

$$\operatorname{Var}(Y_g)=\mu^{(2)}_g+\mu_g^2\bigl(\Theta_1(g)+\vartheta_t\bigr)=:(1+\kappa_g)\mu_g\tag{15}$$

$$\operatorname{Var}(T_{\mathrm{raw}})=\frac{1}{K^2}\Bigl[\sum_{g\in S}\operatorname{Var}(Y_g)+\vartheta_t\Bigl\{\Bigl(\sum_{g\in S}\mu_g\Bigr)^2-\sum_{g\in S}\mu_g^2\Bigr\}\Bigr].\tag{16}$$

<details>
<summary><b>Proof</b></summary>

For (15),
$\operatorname{Var}(Y_g)=\sum_\pi w_\pi^2\operatorname{Var}(I_\pi)+\sum_{\pi\ne\pi'}w_\pi w_{\pi'}\operatorname{Cov}(I_\pi,I_{\pi'})$.
The diagonal term is
$\sum_\pi w_\pi^2\mathbb{E}[I_\pi](1-\mathbb{E}[I_\pi])=A_2(g,t)+O(1/n)$. The
off-diagonal term repeats the argument of Theorem 6, but now layer 1 **can** be
shared because both paths start at $g$; this adds $\Theta_1(g)$. Equation (16) is
the bilinear expansion of $\operatorname{Var}(\sum_gY_g)$ with off-diagonal entries
from (12). $\blacksquare$
</details>

<a name="corollary-7b"></a>
## Corollary 7b — reduction to single-source HetNetEX

Setting $K=1$ in (16) recovers (15), i.e. $\operatorname{Var}(Y)=(1+\kappa)\mu$ with
$\kappa=\mu^{(2)}/\mu+\mu(\Theta_1+\vartheta)-1$. The multi-source theory contains
the single-source theorem as a special case, and makes $\kappa$ **explicit** rather
than parameterised by degree-heterogeneity ratios.

> **A leading effect, not a correction.** The ratio of the covariance contribution
> to the independent contribution is $\approx\vartheta_tK\bar\mu/(1+\bar\kappa)$ —
> **linear in the gene-set size**. Multi-DWPC runs on sets of tens to hundreds of
> genes. Measured: ignoring it understates the null SD by a median factor of 1.474
> overall and 2.39× at $L=4$. This is a correction to standard practice that holds
> whether or not the analytical null is adopted.

<a name="admissibility"></a>
## Admissibility of the same-source kernels

The shared-layer terms in $\operatorname{Var}(Y_g)$ sum over path **pairs**
$\pi\ne\pi'$ from the same source. Forcing layer $\ell$ shared pins interior
positions $\ell-1$ and $\ell$; if that leaves no interior position free, the two
paths are forced to **coincide** and the term must not be counted. The number of
free interior positions is $(L-1)$ minus the number pinned, so

| $L$ | admissible |
|---|---|
| 2 | none — $\operatorname{Var}(Y_g)=\mu^{(2)}_g$ exactly |
| 3 | only $\Theta_1$ and $\Theta_L$; drop the middle kernel $\Theta_2$ |
| $\ge4$ | all |

Cross-source covariance is unaffected: paths from $g\ne h$ already differ at
position 0. This is an **exact combinatorial rule**, though numerically small — it
moved the $L=2$ SD ratio from 0.786 to 0.802.

```python
self_theta = admissible_self_kernels(L, theta_1, theta_mid, theta_L)
```

<a name="theorem-8"></a>
## Theorem 8 — CLT for the aggregate under N2

Suppose $K\to\infty$, no single source dominates
($\max_g\operatorname{Var}(Y_g)/\sum_h\operatorname{Var}(Y_h)\to0$), and
$\vartheta_tK\bar\mu$ is bounded. Then

$$\frac{T_{\mathrm{raw}}-\mathbb{E}T_{\mathrm{raw}}}{\sqrt{\operatorname{Var}(T_{\mathrm{raw}})}}\xrightarrow{d}\mathcal{N}(0,1),$$

and the upper-tail $p$-value is $1-\Phi(z)$.

<details>
<summary><b>Proof sketch</b></summary>

The $Y_g$ form a triangular array with exchangeable leading-order correlation
$\vartheta_t\mu_g\mu_h/(\sigma_g\sigma_h)$. Writing $T_{\mathrm{raw}}$ as a common
factor plus an idiosyncratic remainder and applying Lindeberg–Feller to the
remainder gives the result; the common factor is itself asymptotically normal by
Theorem 9-type arguments for individual DWPCs. $\blacksquare$
</details>

<a name="theorem-9"></a>
## Theorem 9 — Poisson regime and the transform

For short metapaths with $\Lambda_{g,t}=O(1)$ and bounded dependency
neighbourhoods, the path count $N_g$ satisfies
$d_{\mathrm{TV}}(N_g,\mathrm{Poisson}(\Lambda_{g,t}))=O(n^{-1})$ by Stein–Chen, and
$Y_g=\sum_{i=1}^{N_g}W_i$ is compound Poisson with
$\mathbb{E}[W^r]=A_r(g,t)/A_0(g,t)$. Under hub-induced clustering, $N_g$ is better
matched by a negative binomial with mean $\Lambda$ and variance
$\Lambda+\Lambda^2\vartheta^{(0)}_t$, where $\vartheta^{(0)}$ is $\vartheta$
evaluated at $w=0$. Hence $\mathbb{E}[\psi(Y_g)]$ and
$\operatorname{Var}(\psi(Y_g))$ can be evaluated by one-dimensional quadrature
against this compound law rather than by a delta expansion.

> **When the delta method fails, and what to do.** For sparse metapaths most $Y_g$
> are exactly zero and $\psi$ is strongly concave, so the naive second-order delta
> approximation is badly biased — **+126% at $L=2$** in a pilot run. The
> compound-Poisson quadrature reduced this to +4% at $L=2$ but retained 50–60% bias
> at $L=3$–4.
>
> **We therefore recommend the raw-scale aggregate under N2**, where Theorems 6–7
> apply directly. That is exactly the `mean_dwpc` statistic the snapshot track
> already reports. Under N1 the issue does not arise at all (Corollary 3).
>
> We do **not** claim accurate arcsinh-scale network-null moments at $L\ge3$.

<a name="theorem-10"></a>
## Theorem 10 — equivalence with the XSwap stationary law

Let $q_{ij}$ and $r_{ij}$ be the forward and reverse XSwap rates for a candidate
edge $(i,j)$. Then $\lim_{m\to\infty}q_{ij}/(q_{ij}+r_{ij})=d_id_j/m$, the
configuration-model edge probability of A1.

<details>
<summary><b>Proof</b></summary>

$q_{ij}\propto d_id_j$ (the number of half-edge pairs that create $(i,j)$) and
$r_{ij}\propto m-d_i-d_j+1$; the ratio converges as stated. $\blacksquare$
</details>

Hence Theorems 6–7 **target the same null distribution** the Multi-DWPC network
permutation samples from. The discrepancy at finite $B$ is sampling error plus the
$O(1/n)$ terms, not a different hypothesis.

---

# Part III — why the two $p$-values cannot agree

<a name="proposition-11"></a>
## Proposition 11 — the permutation $p$-value is a random variable

Let $p^\ast=\Pr(\tilde T\ge T_{\mathrm{obs}})$ be the exact resampling $p$-value and
$\hat p=(1+\#\{T_b\ge T_{\mathrm{obs}}\})/(1+B)$. Then

1. $\hat p$ is supported on $\{1/(B{+}1),2/(B{+}1),\dots,1\}$; in particular
   $\hat p\ge1/(B{+}1)$ **always**, whatever the evidence;
2. $\mathbb{E}[\hat p]=p^\ast+\dfrac{1-p^\ast}{1+B}$ — conservatively biased upward;
3. $\operatorname{Var}(\hat p)=\dfrac{Bp^\ast(1-p^\ast)}{(1+B)^2}\approx\dfrac{p^\ast(1-p^\ast)}{B}$,
   so the relative standard error $\sqrt{(1-p^\ast)/(Bp^\ast)}$ **diverges** as
   $p^\ast\to0$;
4. $\hat p$ depends on the random seed.

<details>
<summary><b>Proof</b></summary>

$\#\{T_b\ge T_{\mathrm{obs}}\}\sim\mathrm{Binomial}(B,p^\ast)$; (1)–(3) follow
immediately, and (4) restates that the $T_b$ come from a pseudo-random sampler.
$\blacksquare$
</details>

<a name="proposition-12"></a>
## Proposition 12 — the analytical $p$-value is deterministic

$p_{\mathrm{HX}}$ is a deterministic continuous function of the observed statistic,
the degree sequence, and the stratum definitions. It has full support on $(0,1)$, is
reproducible bit-for-bit, and carries no Monte-Carlo variance. Its error relative to
$p^\ast$ is $O(K^{-1})$ under N1 and $O(K^{-1/2})+O(n^{-1})$ under N2 — and **zero**
for the median (Theorem 5b).

<a name="theorem-13"></a>
## Theorem 13 — decomposition of the discrepancy

$$\underbrace{\hat p-p_{\mathrm{HX}}}_{\text{observed}}=\underbrace{(\hat p-p^\ast)}_{\text{(a) Monte-Carlo}}+\underbrace{\Bigl(p^\ast-\tfrac{\lceil(1{+}B)p^\ast\rceil}{1{+}B}\Bigr)}_{\text{(b) discretisation / floor}}+\underbrace{(p^\ast_{\text{grid}}-p_{\mathrm{HX}})}_{\text{(c) model / expansion}}\tag{17}$$

Term (a) has mean $(1-p^\ast)/(1+B)$ and SD $\approx\sqrt{p^\ast(1-p^\ast)/B}$; term
(b) is bounded by $1/(B{+}1)$; term (c) is $O(K^{-1})$ under N1. Terms (a) and (b)
vanish only as $B\to\infty$; (c) only as $K\to\infty$. **There is no regime in which
all three vanish at finite $B$ and finite $K$**, so $\hat p=p_{\mathrm{HX}}$ is a
measure-zero event.

Measured: **0 of 132** features returned numerically identical $p$-values at either
$B=200$ or $B=1000$.

<a name="proposition-14"></a>
## Proposition 14 — rank equivalence

Both $\hat p$ and $p_{\mathrm{HX}}$ are non-increasing in the standardised statistic
$z$ (for $\hat p$ up to grid ties; for $p_{\mathrm{HX}}$ strictly, provided
$|\gamma_1z|<3$). Hence

$$z_{f_1}>z_{f_2}\implies p_{\mathrm{HX},f_1}<p_{\mathrm{HX},f_2}\ \text{ and }\ \hat p_{f_1}\le\hat p_{f_2}\ \text{w.p. }1-O(B^{-1/2}).$$

The two procedures induce the **same ordering** up to Monte-Carlo reversals and
ties, even though they never agree on the numbers. The correct comparison metric is
therefore rank concordance and agreement of the *selected sets*, not $p$-value
equality.

<a name="proposition-15"></a>
## Proposition 15 — permutation resolution bound for BH-FDR

Let $m$ hypotheses be tested with permutation $p$-values at $B$ replicates and
corrected by Benjamini–Hochberg at level $\alpha$. Since $p_{(k)}\ge1/(B{+}1)$ for
every $k$, **no** discovery is possible unless

$$\frac{1}{B+1}\le\frac{\alpha}{m}\qquad\Longleftrightarrow\qquad B\ \ge\ \frac{m}{\alpha}-1.\tag{18}$$

HetNetEX-MD is unaffected, since $p_{\mathrm{HX}}$ has no floor.

<details>
<summary><b>Proof</b></summary>

Immediate from Proposition 11(1) and the BH step-up rule. $\blacksquare$
</details>

**The bound is active in practice.** With $m=11$ metapaths and $\alpha=0.05$, (18)
requires $B\ge219$. In our benchmark at $B=200$ the resampling null made **zero**
FDR-significant discoveries even though six features were strongly signalled; at
$B=1000$ it recovered exactly the same six HetNetEX-MD reports. Realistic queries
enumerate far more than eleven metapaths — with $m=200$ the bound demands
$B\ge3999$.

---

## Limitations and open problems

- **Synthetic networks.** Part I does not depend on this — Lemma 1 and Theorem 5b
  are identities for any fixed score vector — but **every N2 number does**. Real
  Hetionet benchmarks via `hetmatpy` are the main outstanding task.
- **Variance overshoot at $L=2$.** The analytical SD is exact at $L=4$ (ratio 1.044)
  but ~20% too large at $L=2$. Two explanations were tested and **refuted**: the
  canonical/microcanonical ensemble gap (the two references agree with each other to
  5% while the formula overshoots both), and double-counting in
  $\operatorname{Var}(Y_g)$ (the admissibility rule is correct but numerically
  small). Since higher-order shared-edge terms enter with positive sign, truncation
  cannot explain an overshoot either. The cause is open; the direction is
  conservative.
- **Canonical vs microcanonical.** Lemma 4b is the expected-degree ensemble; XSwap
  samples the exact-degree one. Measured below the noise floor here, but Hetionet's
  hubs are far more extreme.
- **Transform under N2.** No accurate arcsinh-scale moments at $L\ge3$; use the
  raw-scale aggregate (Theorem 9).
- **Aggregate choice.** Mean and median are covered exactly; trimmed means and
  weighted quantiles are not.
- **Chain metapaths.** Metapaths revisiting a node type need a self-avoidance
  correction, $O(1/n)$ but not zero.
- **Extreme tails.** $z=49\Rightarrow p\approx10^{-300}$ extrapolates a limit
  theorem far past validation. Report $z$, or truncate $p$, beyond $z\approx6$.
  Theorem 5b avoids this entirely.

---

## References

1. Himmelstein & Baranzini (2015), *PLoS Comput Biol* **11**(7):e1004259 — DWPC and XSwap.
2. Himmelstein et al. (2017), *eLife* **6**:e26726 — Rephetio.
3. Himmelstein et al. (2023), *GigaScience* **12**:giad047 — Hetionet connectivity search.
4. Zietz et al. (2024), *GigaScience* **13**:giae001 — the XSwap edge prior.
5. Hanhijärvi, Garriga & Puolamäki (2009), *SDM* — randomization techniques for graphs.
6. Hájek (1960), *Publ. Math. Inst. Hungar. Acad. Sci.* **5**:361–374 — CLT for SRSWOR.
7. Robinson (1978), *Ann. Statist.* **6**(5):1005–1011 — Edgeworth expansion, finite population.
8. Babu & Singh (1985), *J. Multivariate Anal.* **17**:261–278 — expansions without replacement.
9. Chen (1975), *Ann. Probab.* **3**(3):534–545 — Poisson approximation for dependent trials.
10. Barbour, Holst & Janson (1992), *Poisson Approximation*, Oxford.
11. Molloy & Reed (1995), *Random Struct. Algorithms* **6**:161–180 — configuration model.
12. Benjamini & Hochberg (1995), *JRSS-B* **57**(1):289–300 — FDR control.
