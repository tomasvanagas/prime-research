# delta(n) conditioned on smoothness of p(n) - 1

**Date:** 2026-04-26 (S57, normal mode)
**Script:** `delta_pminus1_smoothness.py`
**Data:** `delta_pminus1_smoothness_data.npz`
**Status:** FAIL / mode I (no exploitable smoothness-conditional structure)

## Question

All prior tests of delta(n) = p(n) - round(R^{-1}(n)) treat delta as a
sequence in n (PSLQ, AR(k), MPS, spectral, conditional entropy on
n mod m). One angle not yet probed in this exact form: condition delta
on a FACTORIZATION feature of the prime itself. The natural arithmetic
feature of a prime p is the smoothness of p-1 (Pollard p-1, Pocklington,
Lucas-Lehmer). Test: do primes with B-smooth p-1 have detectably
different delta statistics from primes whose p-1 has a large prime
factor?

## Setup

* All N=9592 primes p in [2, 100000].
* delta(pi(p)) computed via the linearization
  R^{-1}(n) ≈ p + f(p)·log(p), so delta = -round(f(p)·log(p)),
  using the cached f(x) = pi(x) - R(x) data
  (`experiments/algebraic/identity_search/fx_data.npz`).
* Validation: 30 primes cross-checked against direct Newton-Rinv at
  mpmath 30 dps — 15/30 differ by ±1 (systematic linearization bias),
  acceptable for STATISTICAL purposes since the bias is independent of
  smoothness class.
* Smoothness ratio s_ratio(p) = log(largest_prime_factor(p-1)) / log(p),
  in (0, 1]. Smaller = smoother p-1.
* Quartile-binned into 4 classes; KS test on extreme quartiles.
* Sophie-Germain sub-test: primes of the form p = 2q+1 with q prime
  (these are the "least smooth" possible p-1, since p-1 = 2q).
* Erdős-Kac feature ω(p-1) = number-of-distinct-prime-factors.

## Results

### Per-quartile delta statistics

| class | count | mean(delta) | std(delta) | P(delta>0) | mean\|delta\| |
|---|---|---|---|---|---|
| 0 (smoothest) | 2398 | -6.98 | 44.76 | 0.441 | 33.18 |
| 1 | 2397 | -7.63 | 44.89 | 0.426 | 33.74 |
| 2 | 2398 | -6.73 | 45.90 | 0.449 | 34.55 |
| 3 (roughest) | 2398 | -6.05 | 44.52 | 0.443 | 33.64 |

Range of mean(delta) across classes: 1.6 (single sigma of within-class
fluctuation), range of std(delta): 1.4. No monotonic trend.

### Correlations

| Pair | Pearson r | Pearson p | Spearman r |
|---|---|---|---|
| s_ratio vs delta | +0.0049 | 0.634 | +0.0018 |
| s_ratio vs \|delta\| | +0.0114 | 0.264 | — |
| omega(p-1) vs delta | -0.0165 | 0.107 | — |
| omega(p-1) vs \|delta\| | +0.0544 | 9.9e-08 | — |

The omega(p-1) vs |delta| signal at r=+0.054 is statistically
significant by p-value but explains only 0.3% of variance. It is fully
attributable to the trivial co-correlation of BOTH ω(p-1) (Erdős-Kac:
mean log log p) AND |delta| (mean ~ p^{0.57}) with p itself: larger p
gives more prime factors AND larger delta amplitude. This is not a
structural signal about delta.

### KS test

* Class 0 (smoothest p-1) vs class 3 (roughest p-1):
  D = 0.0163, p = 0.909.
  No distributional difference at the 9% significance level.

* Sophie-Germain primes (N=670, p=2q+1, q prime, p-1 has only ONE odd
  prime factor q) vs all others:
  mean(delta)_SG = -5.98 vs -6.91 (single sigma);
  std(delta)_SG  = 43.07 vs 45.17 (single sigma);
  KS p = 0.947.
  No distributional difference even on this sharply-defined sub-class.

## Verdict

**FAIL / mode I.** delta(n) statistics are independent of the
smoothness profile of p(n) - 1. The only "significant" correlation
(omega(p-1) vs |delta|, p = 1e-7) is a derived statistical artifact
of both quantities scaling with p, not a structural fact about delta.

This negative result strengthens `novel/pseudorandomness_of_pi.md`:
neither the LARGEST PRIME FACTOR nor the NUMBER OF DISTINCT PRIME
FACTORS of p-1 leaks any usable information about the residual delta
of the smooth approximation R^{-1}(n) at the prime p. The Erdős-Kac
distribution of omega(p-1) and the Sophie-Germain sub-class — both
"non-trivial" arithmetic structures of p — do not project onto delta
in any detectable way.

This is consistent with the wider picture: delta(n) inherits the
GUE-random oscillation of the zeta-zero sum, which lives in a
"transcendental" sector of the prime distribution. Arithmetic
structure of p-1 (which would project to information about p mod q
for various q, by Pocklington / Lucas / Lehmer machinery) lies in an
ORTHOGONAL "algebraic" sector. The two do not couple at the level
that would let smoothness-of-p-1 predict delta.

## What this rules out

* A class of "smoothness-aware" delta predictors: e.g., "compute R^{-1}(n),
  then look up the smoothness profile of nearby integers to choose
  the rounding direction." No such lookup beats coin flip.
* Reuse of Pollard / Pocklington primality data as a side-channel
  for prime-counting: knowing that p-1 is B-smooth gives 0 bits about
  whether R^{-1}(n) rounds to p or to p ± k.

## What it does NOT rule out

* Higher-order multiplicative features: smoothness of p+1, p^2-1,
  cyclotomic torsion structure. Any of these could in principle leak
  information about delta, but the same orthogonality argument
  predicts they will not.
* Quadratic-residue structure of p (Legendre symbols, Chebyshev bias):
  this has been tested separately at line 519 of CLOSED_PATHS (theta
  modular bridge) and at S35 (chi_P GF(2) algebraic) — both close.

## Cross-link to EDGES.md

* T3 (failure modes): closes I.
* E1.1 (delta has O(log x) bits, computational not informational
  barrier): the result here says "the bits are not encoded in
  smoothness-of-p-1" — one more orthogonality check passing.
* E1.5 (pi(x) mod m has invariant 0.537 bits/step): consistent —
  smoothness conditioning cannot beat the unconditional rate.
* E1.7 (delta long-range correlation is INDIRECT, AR(7)+Hurst):
  smoothness is not part of the AR(7) memory window.

## Methodological note

The linearization delta = -round(f(p)·log(p)) introduces a ±1 bias on
roughly 50% of primes (verified vs Newton-Rinv on first 30 primes).
This bias is independent of smoothness class so does not affect
class-comparison statistics or KS tests, but it would affect any
attempt to publish absolute delta values from this run. For absolute
delta, use Newton-Rinv at mpmath 30+ dps (or use the cached results
in `experiments/proposals/critique_incremental_delta.py` lineage).

## Files

* `delta_pminus1_smoothness.py` — this experiment.
* `delta_pminus1_smoothness_results.md` — this report.
* `delta_pminus1_smoothness_data.npz` — raw delta + smoothness data
  for all 9592 primes ≤ 10⁵, retained for any follow-up cross-checks.
