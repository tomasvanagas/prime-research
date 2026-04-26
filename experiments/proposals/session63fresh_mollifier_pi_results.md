# P2 — Mollifier-corrected explicit formula

**Script:** `session63fresh_mollifier_pi.py`

## What was tested

The mollifier idea: replace ζ(s) by ζ(s)·M(s) where M(s) =
Σ_{n≤Y} a_n n^{-s} is a Dirichlet polynomial whose values M(½+iγ_j)
vanish for j = 1..K. If the explicit formula for π(x) acquired a
factor M(ρ)/M(1) at each ζ-zero ρ, then suppressing the first K zeros
should let us truncate at much smaller T.

## Method

1. Load 1000 ζ-zeros from `data/zeta_zeros_1000.txt`.
2. For each (K, Y) ∈ {(5,30), (10,50), (20,100), (40,200)}: solve a
   complex least-squares problem for a_1=1, a_2..a_Y minimising
   Σ_{j=1..K} |M(½+iγ_j)|² + λ‖a‖² with λ = 10⁻⁸.
3. Verify: indeed max|M(ρ_j)| < 10⁻⁶ for j ≤ K — the mollifier zeros
   the targeted zeros effectively.
4. Compute *baseline* sharp partial sum
   S_T(x) = R(x) − Σ_{|γ_j|≤T} 2 Re Ei(ρ_j log x)
   and *mollified* partial sum where each term gets a weight
   (M(ρ_j)/M(1)).

## Result

Baseline (sharp) at T = 50, x = 100, 1000, 10000:
errors are **{−0.02, +0.21, +0.44}** — already extremely good.

Mollified with K = 20, Y = 100 at the same T:
errors are **{+0.66, +0.36, −2.07}** — strictly *worse*.

| K, Y | T=50 max |error| | T=1000 max |error| |
|------|------------------|---------------------|
| sharp baseline | 0.44 | 0.29 |
| K=5, Y=30 | 1.81 | 3.21 |
| K=10, Y=50 | 2.07 | 1.73 |
| K=20, Y=100 | 2.07 | 1.93 |
| K=40, Y=200 | 2.07 | 1.72 |

The mollified version is *systematically worse* at every T.

## Why it fails

Multiplying ζ by M does suppress the K targeted zeros, but the explicit
formula for π(x) under ζ·M is **not** "weight each ζ-zero contribution
by M(ρ)/M(1)". The actual formula has TWO contributions:

1. ζ-zeros weighted by M(ρ), AND
2. zeros and poles of M itself, with their own weights.

The naive heuristic "weight the zero terms by M(ρ)/M(1) and stop"
double-counts/distorts the answer because it ignores M's own
contribution — which, since M is a length-Y Dirichlet polynomial,
introduces O(Y log Y) extra terms. To balance the books, you would
need to sum those Y terms anyway: the trade between K (zeros killed)
and Y (mollifier length) is roughly even, with no net polylog gain.

## Verdict

**CLOSED.** Failure mode: **Equivalence (E)** — suppressing K ζ-zeros
via a length-Y mollifier introduces a length-Y compensating sum that
costs essentially as much as those K zeros would have. There is no
free lunch.

## What would change the verdict

A *non-Dirichlet* mollifier (e.g., one based on continuous integral
operators, finite-element bases on the critical line, or compactly
supported wavelets) might break the Y-vs-K duality. But every
example tried in the literature seems to have the same equality of
costs (Bombieri-Friedlander, Conrey-Ghosh).

## One-line summary

The mollifier successfully kills M(ρ_j) for the first K ζ-zeros, but
"weighting each zero contribution by M(ρ)/M(1)" is not the explicit
formula for π under ζ·M; the resulting estimator is strictly worse
than the sharp truncation at every T tested, because it ignores the
length-Y compensation term.
