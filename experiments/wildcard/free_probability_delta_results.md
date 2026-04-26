# Free-Probability Moment Test on δ(x) — Results

**Script:** `free_probability_delta.py`
**Date:** 2026-04-26 (Session 53, fresh-perspective)
**Verdict:** CLOSED — classical CLT (Gaussian) wins, no asymptotic freeness.
**Failure mode:** **E** (equivalence — empirical distribution equates to a
generic limit law from which no algorithmic structure is derivable).

## Goal

The 2026-04-26 morning brainstorm
(`fresh_perspective_session_2026_04_26.md`, item 3) proposed: treat
the contribution of each non-trivial zeta zero ρ to δ(x) at log-scale
as a free random variable. The R-transform composes additively under
free convolution, so if `δ(x)/√x` over a dyadic window matches a known
free distribution (free Wigner / semicircle, free Poisson, Marchenko–
Pastur), one could in principle invert via R-transform and bypass
explicit zero summation.

## What was computed

For `n_max = 2¹⁸ = 262,144`, sieve `π(n)` and evaluate `Li(n)` on the
integer grid (mpmath, ~9.5 s). Form the Selberg-scaled fluctuation

  z_n = (π(n) − Li(n)) · log(n) / √n     — variance ~ O(1).

Standardize z within each dyadic window `[2^j, 2^{j+1})`, compute
sample moments `m_k = E[zₛ^k]` for `k = 1..6`, and compare against
the standardized moments of four reference distributions.

## Reference moment table (standardized: mean 0, variance 1)

| distribution         | m1 | m2 | m3 | m4 | m5 | m6 |
|----------------------|---:|---:|---:|---:|---:|---:|
| Gaussian (CLT)       | 0  | 1  | 0  | 3  | 0  | 15 |
| Wigner / semicircle  | 0  | 1  | 0  | 2  | 0  | 5  |
| Free Poisson(λ=1)    | 0  | 1  | 1  | 3  | 6  | 15 |
| Marchenko–Pastur(c=1)| 0  | 1  | 1  | 3  | 6  | 15 |

The discriminator is `m4` (kurtosis) and `m6`: Gaussian = 3, 15;
semicircle = 2, 5.

## Empirical moments

| window           | n      | m1     | m2 | m3      | m4    | m5      | m6     |
|------------------|-------:|-------:|---:|--------:|------:|--------:|-------:|
| [256,    512)    |    256 | 0.0000 | 1  | +0.2423 | 2.734 | +1.683  | 11.677 |
| [512,   1024)    |    512 | 0.0000 | 1  | −0.0951 | 2.739 | −0.847  | 10.837 |
| [1024,  2048)    |   1024 | 0.0000 | 1  | −0.2501 | 3.154 | −2.443  | 15.851 |
| [2048,  4096)    |   2048 | 0.0000 | 1  | −0.2631 | 2.171 | −1.270  |  6.602 |
| [4096,  8192)    |   4096 | 0.0000 | 1  | −0.1225 | 3.211 | −1.042  | 15.309 |
| [8192, 16384)    |   8192 | 0.0000 | 1  | −0.1334 | 2.380 | −0.740  |  7.811 |
| [16384,32768)    |  16384 | 0.0000 | 1  | +0.1838 | 3.272 | +2.443  | 18.748 |
| [32768,65536)    |  32768 | 0.0000 | 1  | +0.5227 | 2.746 | +3.328  | 11.528 |
| [65536,131072)   |  65536 | 0.0000 | 1  | +0.0400 | 3.085 | +0.098  | 14.965 |
| [131072,262144)  | 131072 | 0.0000 | 1  | −0.1818 | 3.081 | −1.608  | 14.595 |

m1, m2 are forced by standardization. The interesting columns are
**m4 ≈ 3** and **m6 ≈ 15**, both squarely on the Gaussian curve and
far from semicircle (2 / 5). Odd moments m3, m5 hover near zero with
window-to-window jitter consistent with finite-sample noise on a
symmetric-mean distribution; they show no systematic offset toward
+1 / +6 (free-Poisson / MP).

## Closeness verdict (L_∞ over m3..m6 at largest window)

| reference        | L_∞    |
|------------------|-------:|
| **Gaussian**     | **1.61** (dominated by m5 sample noise) |
| Free Poisson(1)  | 7.61   |
| Marchenko–Pastur(c=1) | 7.61 |
| Wigner           | 9.59   |

m4 = 3.08 and m6 = 14.59 confirm the largest-window distribution is
Gaussian to ~3% on the leading even moments. The Wigner / free hypotheses
are off by factors of ~1.5 and ~3 on m4 and m6 respectively — clear
falsification.

## Interpretation

The fluctuation `δ(x)/√x` follows a classical Gaussian limit law on
dyadic windows up to `n = 2¹⁸`. This is consistent with Cramér's
heuristic and with Selberg's CLT for Riemann-zero contributions: the
sum `Σ_ρ x^ρ/ρ` distributes as a sum of (nearly) independent random
phases, and standard CLT applies — not the *free* CLT.

**Algorithmic implication:** there is no free-convolution / R-transform
shortcut for δ(x). Free probability gives non-trivial computational
leverage exactly when joint moments fail Wick's theorem (i.e., are
*non-Gaussian* in the free sense, governed by non-crossing instead of
all pair partitions). Since the empirical second-moment-determines-
fourth pattern (m4 = 3·m2² as in Gaussian, not m4 = 2·m2² as in
semicircle) holds, all higher moments are forced by the variance, and
δ(x) is effectively a *commutative* Gaussian random function.

## Why this is informative (not just yet-another-pseudorandomness check)

The project's `pseudorandomness_of_pi.md` already collects 22+
independent measures showing `π(x) mod 2` and related sequences look
random. Those measures test for *deviation from i.i.d. binary noise*.
This experiment tests a strictly stronger object: the **joint moment
structure** of `δ(x)/√x` as a real-valued process. A semicircle
match would have been compatible with the standing pseudorandomness
findings AND would have implied a non-classical (free) generative
mechanism amenable to R-transform inversion. The Gaussian match
*excludes* that family — closing not just one shortcut but the entire
class of free-probability-style algorithms.

## File / cleanup

- Script: `experiments/wildcard/free_probability_delta.py`
- This results file: `experiments/wildcard/free_probability_delta_results.md`
- Brainstorm doc (this session): `experiments/wildcard/fresh_perspective_session53.md`
- Total runtime: ~9.5 s (mpmath-bound on Li evaluation).
