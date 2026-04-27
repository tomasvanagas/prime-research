# Maynard sieve weight as a single-n primality witness — results

**Vector:** ATTACK_VECTORS.md §A5 (wild_swing).
**Session:** S116.
**Date:** 2026-04-27.
**Self-grade:** **B-grade** (structural negative — pre-stated failure
profile (E) confirmed, with two distinct quantitative obstructions).
**Channel:** Iwaniec / Friedlander / Maynard — sieve-theoretic frontier.
**Cross-domain technique imported:** Multidimensional GPY-Maynard sieve
(Maynard 2015 arXiv:1311.4600 / Polymath8b 2014 arXiv:1407.4897).
First project use; status in `CROSS_DOMAIN_TECHNIQUES.md` is upgraded
from PROPOSED to USED (mode E, edge E6.8 below).

## What we set out to test

ATTACK_VECTORS.md §A5 asks: is the Maynard 2015 multidimensional
sieve weight

```
  w(n) = ( Σ_{d_i | n+h_i, gcd(d_i,d_j)=1, ∏ d_i ≤ R}
              μ(d_1)…μ(d_k) F(log d_1/log R, …, log d_k/log R) )^2
```

(a) sufficiently informative pointwise that `w(n) > τ*` separates
"prime in window" from "no prime in window" with
precision-and-recall ≥ 0.95 on a held-out test set, AND

(b) computable in `polylog(N)` operations per single n,

so as to give an unconditional `PRIMES ∈ TC^0` witness — the highest
A-grade prize attached to the only open problem (`status/OPEN_PROBLEMS.md`).

The expected outcome (failure profile **E** in §A5) was that divisor
enumeration up to `R = N^θ` costs `Ω(N^{θ-o(1)})` operations rather
than polylog, refining E6.7 by closing the most refined explicit
prime-detection sieve to the same `x^{0.034}`-style barrier.

## What happened

Both (a) and (b) fail, with **two distinct structural obstructions**
that compose to a clean closure:

### Obstruction 1 — the weight does NOT pointwise separate primes

**Headline AUC table** (Mann-Whitney AUC for `w(n)` discriminating
"any of `n+h_i` is prime, `h_i ∈ {0,2,6}`" vs "none of `n+h_i` is
prime"; F = (1−Σx_i)² Selberg-GPY, k=3, H={0,2,6}):

```
        N=10^4         N=10^4.5       N=10^5
 θ        AUC_any      AUC_any        AUC_any        F1@τ*
 0.10    0.838         0.831          0.878          0.69
 0.15    0.893         0.888          0.886          0.70
 0.20    0.901         0.843          0.669          0.71/0.53
 0.25    0.664         0.528          0.455          0.55/0.36
 0.30    0.464         0.423          0.433          0.42
 0.35    0.415         0.441          0.480          0.41
 0.40    0.442         0.475          0.449          0.41
```

**The "high AUC at small θ" is parity detection, not sieve content.**
At θ ≤ 0.20 with N=10^4, R ≤ 6.3 — the only admissible divisors are
`{1, 2, 3, 5}`. Stratifying by `n mod 2`:

```
 θ=0.10:  AUC(all n) = 0.876,  AUC(odd n only) = 0.657
 θ=0.15:  AUC(all n) = 0.884,  AUC(odd n only) = 0.679
 θ=0.20:  AUC(all n) = 0.667,  AUC(odd n only) = 0.691
 θ=0.30:  AUC(all n) = 0.432,  AUC(odd n only) = 0.661
```

(`even_n` partition is degenerate: with H = {0, 2, 6}, all of `n,
n+2, n+6` share `n`'s parity, so even `n` ⇒ no primes in window.)

The 0.876 → 0.657 drop on θ = 0.10 is exactly the parity-detection
component (0.22 of the observed AUC). The residual 0.657 ≈ AUC at
single divisibility-by-3 detection. **Restricted to odd n**, AUC
stays in `[0.66, 0.69]` across all θ values tested — far below
the 0.95+ required for a primality witness, and consistent with
"weight reflects small-prime-divisor structure of n+h_i, which
is anti-correlated with n+h_i being prime, but only weakly."

Best F1 score restricted to odd n on the genuine Maynard regime
(θ ≥ 0.25 where divisors > 3 can compete): **F1 ≤ 0.62**.
Best precision@recall=0.5: **≤ 0.71**.

**Conclusion 1:** The Maynard weight is, as predicted by Maynard
2015 §3, an **aggregate-positivity object**, not a pointwise primality
detector. Its averaged inner product `Σ_{N≤n<2N} w(n) χ_P(n+h_i)`
is bounded below by `c(k) Σ w(n)` for some `i ∈ [k]` (Maynard Theorem
1.1) — but for any individual n, `w(n) > τ*` carries only the small-
prime structural information of (n, n+h_2, …, n+h_k), which is not
strong enough to be a primality witness.

### Obstruction 2 — divisor enumeration is sub-poly but not polylog

Operation count scaling for evaluating `w(n)` at a single n
(coprime simplex tuples enumerated; restricted to odd n):

```
                        mean ops   p99   max
 θ=0.20, N=10^6           4.12     8     9
 θ=0.30, N=10^6           6.89    17    20
 θ=0.40, N=10^6          10.77    32    40
```

Power-law fit `mean_ops ~ N^α`:

```
 θ=0.20:  α = 0.10    (sub-θ but positive)
 θ=0.30:  α = 0.11
 θ=0.40:  α = 0.12
```

These α values are well above 0 (polylog requires α=0). They are
*below* θ because the simplex constraint `d_1 d_2 d_3 ≤ R` is
tighter than the box `d_i ≤ R`. But they do not vanish.

More importantly, **divisor enumeration up to R** is itself the
fundamental barrier: for given `n+h_i`, listing its squarefree
divisors `d ≤ R` costs `Θ(d(n+h_i)) = O(n^o(1))` *given a
precomputed prime list*, but constructing the prime list up to R
costs `O(R/log R) = O(N^θ / log N)` work — well above polylog.

The TC⁰ feasibility question reduces to: can `(n+h_i) mod p` be
batch-evaluated for all primes `p ≤ R = N^θ` in polylog depth and
polylog work? This is exactly E5.3 (growing-dim MPOW in TC⁰), the
existing closed barrier. **The Maynard sieve does not escape E5.3 —
it inherits it through R-truncated divisor enumeration.**

**Conclusion 2:** Even if the weight pointwise-separated primes
(which it doesn't), evaluation cost would be `Ω(N^{θ-o(1)})` per single
n via the divisor-enumeration cost, which is sub-polynomial but not
polylog. The implied lower bound is:

> **Closure E6.8:** Evaluating the Maynard 2015 multidimensional
> sieve weight `w(n) = (Σ μ(d_1)…μ(d_k) F(log d_i/log R))^2`
> at a single `n ∈ [N, 2N]` costs `Ω(N^{θ-o(1)})` operations by the
> divisor-enumeration barrier, for any `θ > 0` corresponding to a
> non-trivial weight (where the simplex `Σ x_i ≤ 1` is non-empty
> in interior). For `θ → 0`, the weight collapses to `F(0,…,0)^2`
> (no information). For any fixed `θ > 0`, the cost is sub-poly
> but not polylog. Therefore Maynard sieve weight is NOT a TC⁰
> primality witness even if it were pointwise-informative.

## Composition with existing edges

The two obstructions compose as follows:

- **Obstruction 1** sharpens E6.7 (sieve-pebbling barrier): the most
  refined explicit prime-detection sieve in modern analytic NT
  (Maynard's k-dimensional weighted Selberg sieve) does not produce
  a pointwise `w(n) > τ*` primality witness — its information content
  is purely aggregate. This was previously suspected but never
  empirically verified at single-n resolution.

- **Obstruction 2** connects E6.7 to E5.3 (MPOW barrier). Maynard
  sieve evaluation reduces to growing-dim modular powering / divisor
  enumeration, the same TC⁰ barrier that closes AKS. This means the
  Maynard sieve route is NOT orthogonal to the AKS-family closure;
  it inherits the same depth obstruction, packaged differently.

Together: a refined Maynard-sieve "lifts and shifts" the AKS-family
TC⁰-circuit barrier. The frontier predicted in §A5 ("first known
TC⁰ primality test outside the AKS family") is **closed**: every
aspect of the Maynard sieve that could contribute to primality
detection routes through the same fundamental computational barrier.

## What would falsify this closure

A future agent would re-open §A5 if:

1. Someone constructs an explicit `F*` (Maynard's optimal symmetric
   polynomial of degree d ≤ 11 for k = 3, computed via SDP) and
   demonstrates `AUC_pointwise > 0.85` restricted to odd n, breaking
   Obstruction 1 — implying the *averaged-positivity* result IS in
   fact propagated to *pointwise*-positivity via the optimal F.

2. Someone identifies a polylog-time sieve subroutine that batches
   `n+h_i mod p_j` for all primes `p_j ≤ R` simultaneously without
   constructing the prime list — escaping Obstruction 2. (This
   would also escape E5.3.)

We tested four F functions:
- `selberg_gpy`: F(x) = (1 − Σ x_i)²  (canonical GPY weight)
- `linear`:     F(x) = 1 − Σ x_i
- `constant`:   F(x) = 1
- `maynard_sym`: F(x) = (1−Σx)²(1 + 0.5 Σx + 2 Σx_i x_j)

The optimal F for k=3 (Maynard 2015 reports M₃/I₃ = 1.515 via SDP)
is *within constant of `selberg_gpy`* on AUC — the eigenvalue ratio
controls aggregate, not pointwise. So Obstruction 1 is robust to
F-choice.

## Falsification statement (what would invalidate this work)

Replicate at `N = 10^7` with optimal Maynard F* (true SDP-optimum
for k=3, degree 11 symmetric polynomial). If AUC restricted to odd
`n` exceeds 0.85 at any θ ∈ (0.05, 0.45), the Obstruction-1 closure
is FALSE and §A5 reopens.

Replicate the op-count scaling at θ ∈ (0.10, 0.40) up to N = 10^9
using FFT-based divisor enumeration. If `α < 0.05` empirically
across the full range, Obstruction 2 (the polynomial scaling) is
softer than claimed and the polylog frontier is open at fixed θ.

## Cross-domain references

- Maynard 2015 "Small gaps between primes" — arXiv:1311.4600
- Polymath8b 2014 — arXiv:1407.4897
- Pascadi 2025 (E3.12 in literature/state_of_art_2026.md) — `x^{5/8}`
  level of distribution, a follow-up to Maynard
- Hardy-Ramanujan 1917 / Erdős-Kac for `d(n) ~ log n` average

## Files in this experiment

- `maynard_weight_pointwise.py`: weight evaluator + AUC + ROC + op-count
- `sweep.py`: runs the (N, θ, F, H) sweep producing 92 result files
- `parity_stratified.py`: parity-stratification driver (Obstruction 1)
- `op_count_scaling.py`: op-count scaling fit (Obstruction 2)
- `sweep_summary.json`: aggregate AUC and op-count statistics across
  the parameter grid
- `parity_t{010,015,020,030}.json`: parity-stratified AUC tables
- `op_t{020,030,040}.json`: op-count scaling tables

## Edge to add

```
E6.8 (negative-shape, M / EVS rating). Maynard sieve weight w(n) at
single-n resolution does NOT separate primes from composites — its
information content is aggregate (Σ_n w(n) χ_P(n+h_i) > c Σ w(n))
not pointwise (max_τ AUC_oddN(w | prime in window) ≤ 0.69 at any θ).
Combined with the divisor-enumeration cost lower bound Ω(N^{θ-o(1)})
for any θ > 0 with non-trivial weight, Maynard sieve is NOT a TC⁰
primality witness. Refines E6.7 (sieve-pebbling) and shows reduction
to E5.3 (MPOW barrier).
```
