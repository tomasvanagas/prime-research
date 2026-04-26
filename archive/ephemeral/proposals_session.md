# Proposals Session — Fresh Attack Vectors for Polylog p(n)

Date: 2026-04-26
Author: proposal-session
Constraint: written WITHOUT consulting CLOSED_PATHS.md to avoid bias toward
prior negative results. Some proposals may overlap with closed paths in
spirit; the *test designs* are still useful because they sharpen what
specifically about prior attempts failed.

Target: compute p(n) — the nth prime — exactly in O(polylog n).
Reduction: p(n) reduces to a single inversion of pi(x), so all proposals
attack pi(x) and the residual delta(n) = p(n) - R^{-1}(n).

---

## Proposal 1 — Compressed Sensing on the Zeta-Zero Contribution Matrix

### Idea
From Riemann's explicit formula
    pi(x) = R(x) - sum_{rho} R(x^rho) - small corrections.
Build, for a target window of x-values, the "contribution matrix"
    M[i, j] = R(x_i^{rho_j}) + R(x_i^{conj(rho_j)})
                     for x_i in [N, 2N], j = 1..K zeros.

The known fact is that the *full* sum is large (each row sums to a
correction comparable in magnitude to log x); the *novel question* is
whether M itself has very fast singular value decay (numerical
low-rank-ness) or sparse representation in some basis (Fourier in t,
wavelets in x, etc.). If sigma_k(M) decays like 2^{-k} or faster, then
O(polylog) zero contributions on a clever basis suffice to reconstruct
pi(x) up to 1/2.

### Pseudocode
```
def build_M(x_lo, x_hi, K, rhos):
    xs = arange(x_lo, x_hi)
    M = zeros((len(xs), K))
    for j, t in enumerate(rhos[:K]):
        rho = 0.5 + 1j*t
        col = R(xs**rho) + R(xs**conj(rho))
        M[:, j] = col.real
    return M

def test():
    M = build_M(10_000, 20_000, K=2000, rhos=load_zeros(8000))
    U, S, Vt = svd(M, full_matrices=False)
    # do singular values decay?
    # also check sparsity in DCT-2D, FFT-2D, Haar wavelet bases
```

### Time complexity
If the *m*-truncated SVD captures pi(x) to within 1/2 over the window
and m grows polylog in N, the cost to evaluate pi at any single x in
window is O(m * polylog) — *after* an offline precomputation of the m
basis vectors (which is shared across all x in window). For asymptotic
polylog this needs the basis to be data-independent (e.g., a fixed
Fourier basis) AND the coefficient of pi(x) in that basis to be
computable from x alone in polylog time.

### Key assumption / conjecture
**Compressibility conjecture.** The matrix M (or some normalized
variant) has effective rank or sparsity polylog in the window size,
in some basis admitting closed-form coefficient extraction.

This is a *concrete* version of "the explicit formula has secret
structure beyond GUE statistics that lets us avoid the O(sqrt(x))-zero
truncation cost."

### Test for n < 10000
1. Build M for x in [1, 10000], K = 8000 zeros.
2. Compute SVD; report sigma_k/sigma_0 for k = 1, 2, 4, 8, …, 256.
3. If sigma_k <= 2^{-k}, **the conjecture survives a first test**.
4. Try Fourier-2D and Haar wavelet representations, count significant
   coefficients (|c| > 0.01 * ||M||_F).

**Falsifier:** flat singular spectrum (sigma_k ~ 1/sqrt(k), GUE-like)
or no basis with sub-linear sparsity. This is the most likely outcome
and is directly informative.

---

## Proposal 2 — Targeted PSLQ Hunt for delta(n) on Sparse Subsequences

### Idea
Brute PSLQ on the full sequence delta(n) = p(n) - R^{-1}(n) is doomed
because the sequence is incompressible globally. But maybe:

- delta(2^k) for k = 1..20 has a closed form.
- delta(p_k) restricted to k a Fibonacci index.
- delta(p_{k!} mod something).

These are "sparse probes" — enough samples to find low-degree algebraic
relations against a dictionary of fundamental constants. We're not
trying to predict ALL primes; we're looking for *any* infinite
subsequence on which delta has structure. Even one such subsequence
would be a major finding.

### Pseudocode
```
import mpmath
mpmath.mp.dps = 60

dictionary = [
    mpmath.mpf(1),
    mpmath.zeta(2), mpmath.zeta(3), mpmath.zeta(5),
    mpmath.euler,
    mpmath.log(2), mpmath.log(3), mpmath.log(2*mpmath.pi),
    mpmath.catalan,
]

for k in range(1, 21):
    n = 2**k
    d = delta_at(n)
    for scale in [1, log(n), sqrt(n)]:
        target = d * scale
        rel = mpmath.pslq([target] + dictionary, tol=1e-30)
        if rel and max(abs(r) for r in rel) < 1e6:
            print(f"k={k}, scale={scale}: {rel}")
```

### Time complexity
PSLQ trivial for a single search. The *resulting* algorithm, if a
relation is found, would be O(polylog n) on the subsequence.

### Key assumption / conjecture
**Subsequence integrability.** There exists an infinite, polylog-
indexable subsequence S = {n_k} and a closed-form formula F(k, fundamental
constants) such that delta(n_k) = F(k).

Much weaker than "delta is integrable globally" — even local
integrability on n = 2^k would be a structural breakthrough.

### Test for n < 10000
PSLQ at k = 1..13 (n_k = 2^k). 60 digits of precision. Try several
scalings and dictionary sizes. Success criterion: any nontrivial
relation with small integer coefficients (< 100) holding to within 1e-25.

---

## Proposal 3 — Tensor Train (TT) Decomposition of the Prime Indicator

### Idea
Treat chi_P : {0, 1}^L -> {0, 1} (the indicator of "is x prime" with
x in binary) as an order-L tensor of shape (2, 2, …, 2). Compute its
tensor train (TT) decomposition. The TT-rank r tells us how many
parameters are needed to encode primes up to 2^L *deterministically and
exactly*.

Crucially: prefix sums of a TT can be computed in O(L * r^3) time. If
r = polylog(L), pi(2^L) is computable in polylog(2^L) time.

This is a measurement of an information-theoretic quantity not
previously run *as TT-decomposition* on the prime indicator. (Related
work has measured MPS bond dim asymptotics as a theorem; this is a
direct numerical experiment that tests the same question with a
deterministic construction.)

### Pseudocode
```python
import numpy as np
from sympy import isprime

L = 18
chi = np.array([1.0 if isprime(i) else 0.0 for i in range(2**L)])
T = chi.reshape((2,)*L)
ranks = [1]
A = T.reshape(2, -1)
for ell in range(L-1):
    U, S, Vt = np.linalg.svd(A, full_matrices=False)
    eps = 1e-12
    r = int(np.sum(S > eps * (S[0] if S[0] > 0 else 1)))
    ranks.append(r)
    A = (np.diag(S[:r]) @ Vt[:r]).reshape(r*2, -1)
print("TT-ranks:", ranks)
```

### Time complexity
- TT decomposition itself: O(2^L * max_rank^2) per SVD step.
- Once TT obtained, prefix-sum (= pi(x)): O(L * max_rank^3).
- IF max_rank = polylog L, polylog pi(x) achieved.

### Key assumption / conjecture
**Low TT-rank conjecture.** chi_P on [0, 2^L) has TT-rank polylog L
when sites are ordered by binary digit (most significant first).

A *random* binary string has TT-rank min(2^ell, 2^{L-ell}) — exponential
in L/2. So "primes have low TT-rank" would directly contradict the
GUE-randomness picture and would be a major positive finding.

### Test for n < 10000
L = 13 (n < 8192) or L = 14 (n < 16384) is easy. Compute TT-ranks
across all bipartitions. Compare with:
- random Bernoulli(p = pi(2^L)/2^L),
- a known-low-TT-rank function (n divisible by k for small k),
- a known-high-TT-rank function (random binary).

**Strong success:** max TT-rank < L^2.
**Mild success:** max TT-rank significantly below random-Bernoulli baseline.
**Likely outcome:** TT-rank scales as ~2^{L/2} (random-like).

This experiment is *cheap* and gives a clean numerical verdict
on whether tensor-network descriptions of primes admit any compression
at all.

---

## Proposal 4 — GUE-Aware Adaptive Importance Sampling of Zeros

### Idea
The contribution of zero rho_j to pi(x) is roughly cos(t_j * log x) /
(t_j * log t_j). Standard truncation uses the FIRST K zeros — but for
*specific* x the most important zeros may be elsewhere (those whose
phase t_j * log x is near 0 mod 2 pi, and whose 1/t-decay places them
near a "resonance").

Treat zero positions as a determinantal point process (the GUE
sine-kernel limit), and compute pi(x) as an *expectation* under this
DPP, using:
1. The first K_low << sqrt(x) zeros (deterministically), giving O(polylog)
   "smooth" contribution; AND
2. A statistical correction based on GUE moments — the analytic
   *variance* of the truncation tail.

If the variance is provably below 1/4, rounding to nearest integer
gives the *exact* pi(x) with high probability. With repeated trials
this becomes a Las Vegas algorithm with exponential confidence.

### Pseudocode
```
def pi_gue(x, K_low):
    # Deterministic prefix
    main = R(x) - sum(2*Re(R(x**rho_j)) for j in range(K_low))
    # Estimate variance from j > K_low using GUE pair-correlation
    sigma2 = analytic_tail_variance(K_low, x)
    if sigma2 < 1/8:
        return round(main)
    else:
        return None  # Need more zeros
```

### Time complexity
O(K_low * polylog) = O(polylog x) if the GUE tail-variance estimate
is below 1/4 with K_low = polylog. Whether this is true depends on
cancellation rates tied to moments of the zero-zero correlation.

### Key assumption / conjecture
**Variance-truncation conjecture.** Under GUE statistics, the
truncation tail of the explicit formula has variance O(1/T^c) for some
c > 0, with K_low = polylog chosen such that the tail variance is < 1/4.

The standard "K = sqrt(x)" comes from worst-case bounds; this attempts
to exploit GUE *cancellation* via random matrix moment results
(Keating-Snaith, etc.).

### Test for n < 10000
For x in {100, 200, 500, 1000, 2000, 5000}:
1. Compute pi_gue(x) with K_low = 10, 20, 50, 100.
2. Compare to ground truth.
3. Measure empirical variance of truncation tail vs GUE prediction.

---

## Proposal 5 — Modular / Theta Function Bridge (theoretical)

### Idea
Theta functions admit O(polylog precision) evaluation via the modular
transformation theta(-1/tau) = sqrt(-i tau) theta(tau). Riemann's xi is
essentially a Mellin transform of theta, which is why zeta has a
functional equation. The chain
    theta -> Mellin -> zeta -> Perron -> pi(x)
goes through O(sqrt(x)) zeros at the last step.

**Question:** can we write pi(x) as the Mellin transform of a
modular-or-near-modular function f(tau), without going through the
zeta zeros at all? If so, evaluate f via modular acceleration, then
take its Mellin transform via Gauss-Legendre quadrature in O(polylog).

### Key assumption / conjecture
**Off-line evaluability.** log zeta(s)/s is computable at any s in
the critical strip in O(polylog |s|). Currently: only known on the
critical line via Riemann-Siegel. Off-line evaluation requires
Euler-Maclaurin, costing ~O(|s|^c).

### Test for n < 10000 (preliminary)
Numerical contour integration of the Perron formula with truncation
T = 50, 100, 200; use Riemann-Siegel for zeta(1/2 + it). Measure error
vs T scaling. If error decreases like 1/T^c with c > 1/2, that's
encouraging (cheaper than naive sqrt(x) truncation).

This proposal is the most speculative; included for completeness.

---

## Cross-cutting test plan

For each tested proposal:
1. Save runnable Python to `experiments/proposals/<name>.py`.
2. Save findings to `experiments/proposals/<name>_results.md`.
3. Report verdict: SURVIVES / CLOSED / INCONCLUSIVE.
4. If CLOSED, classify failure mode: C (circularity), E (equivalence),
   I (information loss), R (resource), N (no structure).

## Priority order for testing in this session

1. **Proposal 3 (TT-rank)** — cheapest, structurally informative.
2. **Proposal 1 (compressed sensing)** — moderate cost, clean spectrum.
3. **Proposal 2 (PSLQ subsequence)** — cheap, low expected payoff.
4. **Proposal 4 (GUE truncation)** — needs variance analysis.
5. **Proposal 5 (theta bridge)** — too theoretical for this session.

---

# === SECOND PROPOSAL BATCH (same date, fresh angles) ===

The prior batch above (compressed sensing on zero contribution matrix,
PSLQ on subsequences, TT-rank, GUE truncation, theta bridge) overlaps
heavily with what's already runnable under `experiments/proposals/`.
The five proposals below intentionally pivot to angles that are
**structurally distinct** from those frames, with different mathematical
objects of attack.

| # | Frame                                            | Object of attack                               | Win condition                                    |
| - | ------------------------------------------------ | ---------------------------------------------- | ------------------------------------------------ |
| B1 | Hierarchical-matrix kernel compression          | The explicit-formula sum as a low-rank operator | HSS rank << #zeros for the relevant kernel slab  |
| B2 | Reservoir-computing dynamics on delta(n) bits   | The bit-residual sequence as a chaotic signal  | ESN beats a random predictor on bit 0 of delta   |
| B3 | Determinant-representation hunt for pi(x)       | A small auxiliary matrix M(x) with det = pi(x) | Empirical agreement on first 10^4 values         |
| B4 | Modular-form / Hecke-trace coupling             | pi(x) as a difference of two trace formulas    | Algebraic cancellation reduces to polylog work   |
| B5 | Self-correcting recursive halving p(n)->p(n/2)  | A recursion with provable polylog corrector    | A working recursion with measurable shrink rate  |

---

## Proposal B1. Hierarchical-matrix (H-matrix) compression of the explicit-formula kernel

### Idea

The Riemann explicit formula gives, schematically,
    psi_0(x) = x - sum_rho x^rho/rho - log(2 pi) - (1/2) log(1 - x^{-2}).
For x evaluated at a vector of test points {x_1, ..., x_M} this is the
matrix-vector product [psi_0(x_i)]_i  =  K v with K_{i,k} = x_i^{rho_k}/rho_k.
Naive cost of M evaluations using N zeros is `O(M N)`.

**Compressed sensing attacks the *vector* v** (sparsity in zero basis,
prior batch's Proposal 1). **This proposal attacks the *kernel* K itself.**
If K admits an H-matrix / HODLR / HSS decomposition with hierarchical
ranks r << N, then evaluating it on M new points costs `O((M+N) r log)`.

The novel question: does K(x, gamma) = x^{1/2 + i gamma}/(1/2 + i gamma)
have the **off-diagonal low-rank structure** that H-matrices exploit? The
GUE statistics of zeros are local; H-matrices care about pair (i,k)
"distance" |x_i - x_k| or |gamma_i - gamma_k|. The kernel oscillates in
the second slot at rate log x, so for nearby gammas the kernel rows are
nearly proportional --- which is exactly the H-matrix favorable case.

### Pseudocode

```python
import numpy as np
from mpmath import zetazero
def kernel_block(x_grid, gamma_grid):
    s = 0.5 + 1j*gamma_grid
    return np.exp(np.log(x_grid)[:,None] * s[None,:]) / s[None,:]
def hssranks(K, eps=1e-10):
    n = K.shape[1]
    out = []
    for level in range(int(np.log2(n))):
        m = 2**level
        block = K[:, :m]
        s = np.linalg.svd(block, compute_uv=False)
        out.append(int(np.sum(s > eps*s[0])))
    return out
gammas = np.array([float(zetazero(k).imag) for k in range(1, 65)])
xs = np.linspace(1e3, 1e4, 64)
print(hssranks(kernel_block(xs, gammas)))
```

### Complexity

If empirical HSS ranks are O(polylog), one-time `O(N r^2)` precomputation
gives O(N r log N) per evaluation. **Crucial caveat:** the win is contingent
on shrinking N from O(sqrt(x)) to polylog(x) using H-matrix-style
hierarchical cancellation. That extra step is the conjecture.

### Key conjecture

The off-diagonal blocks of K(x, gamma), restricted to any fixed dyadic
window of x and zeros up to height T, have epsilon-ranks
O(log T  log(1/eps)).

### Test for n < 10000

64 test points x in [10^3, 10^4], 64 zeros, build K, compute rank-vs-block
size curve. Plateau at O(log) -> genuine effect; linear growth -> closed.

Code: `experiments/proposals/hmatrix_kernel_compression.py`.

---

## Proposal B2. Reservoir-computing dynamics on delta(n) bit residuals

### Idea

Reservoir computing (echo state networks, ESN) is the **dequantized cousin
of recurrent quantum systems**: a fixed random recurrent network with only
the readout layer trained. Recent work shows ESNs learn chaotic attractors
with state size polylog in prediction horizon. Riemann zero phases are
GUE-random but **deterministic** --- they *are* a chaotic dynamical system
in the Hilbert-Polya picture.

Train an ESN on bit-0 of delta(n) := p(n) - round(R^{-1}(n)) for n in a
training window. Input at step n: log-features (log n, R^{-1}(n) mod 1,
n mod q for small q). Output: bit-0 of delta(n). If the ESN beats 50% on
held-out window, **the entropy barrier is partially broken**.

Distinct from gradient-boosting attempts because the ESN has **internal
recurrent dynamics** that hold positional memory across n -> n+1
transitions.

### Pseudocode

```python
import numpy as np
class ESN:
    def __init__(self, in_dim, res_dim=512, sr=0.95, sparsity=0.1):
        self.W_in = np.random.randn(res_dim, in_dim) * 0.1
        W = np.random.randn(res_dim, res_dim)
        W *= (np.random.random(W.shape) < sparsity)
        W *= sr / max(abs(np.linalg.eigvals(W)))
        self.W = W
        self.h = np.zeros(res_dim)
    def step(self, u):
        self.h = np.tanh(self.W @ self.h + self.W_in @ u)
        return self.h
    def fit(self, U, y):
        H = np.array([self.step(u) for u in U])
        self.W_out = np.linalg.lstsq(H, y, rcond=None)[0]
```

### Complexity

Training `O(T res_dim^2)`, inference `O(res_dim^2)`. Polylog inference if
res_dim is polylog.

### Key conjecture

There exists a finite-dimensional dynamical system whose orbit, projected
to its first coordinate, agrees with bit-0 of delta(n).

### Test for n < 10000

Train on n in [2000, 7000], evaluate on n in [7001, 10000]. Win = >2%
above majority-class baseline.

Code: `experiments/proposals/esn_delta_bits.py`.

---

## Proposal B3. Determinant-representation hunt for pi(x)

### Idea

Many counting problems admit determinant representations: spanning trees
(matrix-tree), tilings (Kasteleyn), non-intersecting paths (LGV),
constructible-sheaf characteristics (Frobenius). **Hypothesis to test**:
there exists a polynomial-size matrix M(x) with polylog-time computable
entries such that det M(x) = pi(x) (or a simple algebraic variant).

Empirical attack: for k=2, 3, ..., parameterise M(x) by k^2 functions
chosen from a feature library {1, log x, R^{-1}(x), li(x), x^{1/2},
x mod q, ...}. Score by L2 error of det M(x) - pi(x) over x = 10..10^4.

### Pseudocode

```python
import itertools, numpy as np
def hunt(k=2, F_lib=...):
    best = (np.inf, None)
    for choice in itertools.product(F_lib, repeat=k*k):
        M = build(choice)  # shape (k, k, T)
        d = np.array([np.linalg.det(M[..., t]) for t in range(T)])
        err = np.linalg.norm(d - pi_target)
        if err < best[0]:
            best = (err, choice)
    return best
```

### Complexity

If a positive hit at small k exists, det M(x) is `O(k^omega)` = polylog,
giving directly an O(polylog) algorithm.

### Key conjecture

pi(x) is in the algebraic-circuit class VP rather than VNP --- i.e.,
expressible by polynomial-size determinants over polylog-computable entries.

### Test for n < 10000

k=2 with 12-feature library = 12^4 = 20736 candidates. Top-10 by L2
error vs pi on 1000 evaluation points.

Code: `experiments/proposals/det_representation_hunt.py`.

---

## Proposal B4. Modular-form / Hecke-trace coupling (sketch only)

### Idea

For weight-2 newform f on Gamma_0(N), Eichler-Selberg gives
sum_{p<=x} g(a_p(f)) in time *independent of x* for polynomial g.
If we find f and g such that g(a_p) = 1 for all p, then pi(x) inherits
fast evaluation. Hard part: a_p has Sato-Tate distribution; no fixed
polynomial g equals 1 on it. A *family* of forms {f_N} with growing
levels gives more spectrum to play with --- a CRT-style combination
might work.

This is an organising sketch, not yet a runnable experiment. The "build a
polynomial that equals 1 on Hecke spectra" is itself research-grade.

### Test for n < 10000

Pick N = 11 (smallest level with cusp form), enumerate a_p for p < 10^4
(can use mpmath / sage), fit polynomials g of deg <= 4 to minimize
sum |g(a_p) - 1|. If best fit error stays bounded as p grows -> open;
diverges -> closed.

Code: `experiments/proposals/hecke_trace_sketch.py` (data-only this session).

---

## Proposal B5. Self-correcting recursive halving p(n) -> p(n/2)

### Idea

Define theta(n) := p(n) - 2 p(n/2). By PNT, the leading O(log n) "smooth"
terms in p(n) and 2 p(n/2) cancel (modulo a `n log 2` shift). The remaining
*chaotic* part of theta might be **smaller than the chaotic part of p(n)
alone** because of GUE phase coherence between near-multiplicative scales.
If theta has a polylog-predictable structure, recurse: depth log n levels,
polylog work each, total polylog^2.

### Pseudocode

```python
def p_via_recursion(n, base_table):
    if n < len(base_table):
        return base_table[n]
    a = p_via_recursion(n//2, base_table)
    theta_est = theta_predictor(n, a)
    cand = 2*a + theta_est
    return correct_locally(cand, n)  # local prime walk
```

### Complexity

If theta_predictor is `O(log n)` and correct_locally shrinks search to
`O(log^c n)`, total is `O(log^{c+1} n)`.

### Key conjecture

theta(n) := p(n) - 2 p(n/2) is *more predictable* than p(n) directly,
because the leading 50% of digits cancel and the remaining 50% partially
re-cancel via near-multiplicative GUE coherence.

### Test for n < 10000

Compute theta(n) for n = 100, 200, ..., 10000. Fit theta(n) - n log 2
to polynomial in log n. If residual / sqrt(n) -> 0: open, if not: closed
(reduces to direct estimation problem).

Code: `experiments/proposals/halving_recursion.py`.

---

## Honest priors

| Proposal      | Prior of new info | Risk of duplicate-shape          |
| ------------- | ----------------- | -------------------------------- |
| B1 H-matrix   | medium-high       | low (kernel-side, not signal)    |
| B2 Reservoir  | low-medium        | medium (ML-on-residuals family)  |
| B3 Det-hunt   | medium            | low (algebraic-circuit family)   |
| B4 Modular    | low (sketch only) | medium (Hecke family)            |
| B5 Halving    | medium            | low (recursion is fresh)         |

I will run B1, B2, B3, B5 as live experiments and leave B4 as a sketch.
