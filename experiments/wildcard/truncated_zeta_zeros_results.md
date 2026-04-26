# Truncated zeta polynomial zeros vs Riemann zeros — Results

## Question
Riemann's explicit formula uses the nontrivial zeros `rho` of `zeta(s)`. Computing each true zero requires Riemann-Siegel-style evaluation, costing `O(sqrt(t))` per zero. The partial sum `zeta_N(s) = sum_{n=1..N} n^{-s}` is trivial to evaluate (`O(N)` per point). **If the zeros of `zeta_N` on the critical line approximated the true zeros**, we'd have a cheap source of zero approximations and could plug them into the explicit formula.

## Setup
For each `N ∈ {5, 10, 20, 50, 100, 200, 500}`:
1. Sweep `t ∈ (0, 60)` with 1200 samples on `s = 1/2 + it`.
2. Find local minima of `|zeta_N(1/2 + it)|` below threshold 0.2 → candidate zero locations.
3. For each of the first 10 true Riemann zeros, find the closest candidate. Report avg and min error.

## Numbers

| N   | # minima | avg err | min err |
|----:|---------:|--------:|--------:|
| 5   | 3        | 3.78    | 0.13    |
| 10  | 12       | 0.28    | 0.03    |
| 20  | 13       | 0.75    | 0.03    |
| 50  | 10       | 1.42    | 0.06    |
| 100 | 15       | 0.83    | 0.08    |
| 200 | 8        | 4.29    | 0.17    |
| 500 | 11       | 1.34    | 0.23    |

True Riemann zeros (heights): 14.13, 21.02, 25.01, 30.42, 32.94, 37.59, 40.92, 43.33, 48.01, 49.77.

## Interpretation
**Error does NOT converge to zero as N grows.** It bounces in the range 0.3–4.3 with no monotone trend. The min-error column shows occasional close hits (≈ 0.03), but these are coincidental — adjacent N values give very different "nearby" minima.

This is the expected mathematical behavior:
- `zeta(s)` has its nontrivial zeros only via *analytic continuation* from `Re(s) > 1` to the critical strip.
- The partial sum `zeta_N(s)` is a finite Dirichlet polynomial; on `Re(s) = 1/2` it has no relationship to the analytic-continuation zeros.
- In particular, by the Lindelöf-class theory of Dirichlet polynomials, `zeta_N(1/2 + it)` oscillates with `O(N^{1/2})` magnitude and `O(log N)` zero density per unit t-interval — completely different from `zeta`'s zero density `(1/2π) log(t/2π)`.

A truncated Dirichlet polynomial cannot recover analytic-continuation behavior. The functional equation `zeta(s) = chi(s) zeta(1-s)` is what produces the critical-line zeros, and truncation breaks it.

## What about Riemann-Siegel-style truncation?
The Riemann-Siegel formula `Z(t) = 2 sum_{n=1..M} cos(theta(t) - t log n) / sqrt(n) + R(t)` with `M = floor(sqrt(t/2π))` *does* approximate `zeta(1/2 + it)` accurately, including its zeros. But:
- `M = sqrt(t/2π)` means we need `M` up to `sqrt(t)` zeros, where `t` ranges to `x log^3 x` for integer-precision pi(x).
- So Riemann-Siegel costs `O(sqrt(t))` per zero, same as direct evaluation, no speedup.

A "polynomial-zero shortcut" requires a polynomial whose zeros approximate Riemann zeros — and no such polynomial of polylog size is known.

## Verdict
**CLOSED — failure mode (I) Information loss.** Truncating the Dirichlet series destroys the analytic continuation, and the partial-sum zeros bear no useful resemblance to Riemann zeros. The truncated polynomial cannot serve as a polylog-time substitute.

## What would change the verdict
A *non-truncation* polynomial approximation: e.g., a polynomial in `s` of polylog degree whose roots approximate the first `T` Riemann zeros within `1/T`. Such a "compressed zero polynomial" would be a major analytic discovery. The Hardy `Z`-function has the right zeros but is transcendental, not polynomial; rational-function approximations (Padé) for `Z` show `O(N)` poles needed near each zero cluster — no compression.
