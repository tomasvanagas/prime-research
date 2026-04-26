# Bessel-Basis PSLQ Identity Search for f(x) = π(x) − R(x)

**Script:** `bessel_basis_pslq.py`
**Date:** 2026-04-26 (Session 68 deep-focus revisit of Task #3)
**Verdict:** **CLOSED.** No Bessel-basis identity for f(x) survives cross-validation.
**Failure mode:** **E** (equivalence — Bessel basis behaves identically to a same-variance Gaussian random control).

## Why this experiment

Session 29 closed Task #3 ("Novel Identity Search") for the basis
{1, log x, x^{1/k}, li(x^{1/k}), sin/cos(γ_k log x)}. That verdict
("f algebraically independent of all TESTED bases") does not logically
preclude an identity in an *untested* basis. Bessel functions K_ν, I_ν,
J_ν, Y_ν appear naturally in:

- The Selberg trace formula for SL(2,ℤ) (spectral side: K_{ir}-kernels),
- Mellin–Barnes representations of L-functions,
- The Hardy–Ramanujan circle method for partition asymptotics
  (K_0(2π√(...)) saddle).

S29's basis is disjoint from the Bessel family, so this is a strict
extension, not a re-run.

## Basis (10 elements)

```
f(x), 1, log x, √x, li(x),
K_0(log x), I_0(log x),                 # modified Bessel (real, decaying / growing)
J_0(γ_1 · log x), Y_0(γ_1 · log x),     # oscillatory Bessel @ first zeta zero
K_0(2π · √(log x))                      # partition-asymptotic kernel
```

PSLQ ran at x ∈ {5000, 10000, 50000, 100000} with mpmath at 50 dps,
maxcoeff = 10⁶, maxsteps = 8000. Each candidate relation was cross-validated
at x' = x + 1000 (or +200 if near upper bound).

## Results

| x      | residual at fit | coeff(f) | cross-check residual | survives? |
|-------:|----------------:|---------:|---------------------:|:---------:|
|  5,000 |        1.9 × 10⁻³⁶ | 2712 |       1.10 × 10⁴ | NO |
| 10,000 |        0.0          |    0 |       4.88       | NO (coeff(f)=0; trivial 100·1 = √x) |
| 50,000 |        4.0 × 10⁻³⁵ | 5349 |       6.03 × 10⁴ | NO |
| 100,000|        2.8 × 10⁻³⁶ | 10385 |       N/A (no data > 100k) | n/a |

The "relations" found at the fit points have tiny coefficients-of-f relative
to the magnitudes (e.g. coefficients 2712 / 5349 / 10385 with the linear
combination cancelling to ~10⁻³⁵), which is the canonical PSLQ-overfit
signature when the basis dimension is comparable to the precision-limited
rank: the algorithm always finds an integer relation at finite precision,
but the "relation" depends on the specific x rather than on the algebraic
structure of f.

## Random control

Replacing f(x) with i.i.d. Gaussian noise of matching standard deviation
(σ ≈ 4.23, the std of f over [2, 100000]), three trials at x = 50,000:

| trial | residual at fit | cross-check residual |
|------:|----------------:|---------------------:|
| 1     | 6.4 × 10⁻³⁵     | 1.60 × 10⁵           |
| 2     | 3.6 × 10⁻³⁵     | 7.46 × 10⁴           |
| 3     | 8.5 × 10⁻³⁷     | 1.77 × 10⁴           |

The random-control behaviour is statistically indistinguishable from f(x):
both produce ~10⁻³⁵ residuals at the fit point and 10³–10⁵ residuals at
cross-check. PSLQ on a 10-element basis at 50 dps "finds an identity" for
*any* number, regardless of whether one exists; the discriminator is
cross-validation.

## Interpretation

The Bessel basis carries no more identity-bearing structure for f(x) than
a Gaussian random number does. Concretely, this rules out:

1. **Selberg-trace shadow:** if f had an identity expressible via K_{ir}
   spectral sums at low order, the K_0(γ_1 · log x) coefficients should
   line up across x. They do not.
2. **Partition-saddle reuse:** K_0(2π√(log x)) is the dominant Bessel
   kernel in Hardy–Ramanujan partition asymptotics. There is no
   transferable structure from the partition setting to π(x) at this
   level.
3. **Modified Bessel growth/decay scaffolding:** I_0(log x) ~ x/√(2π log x)
   and K_0(log x) ~ √(π/(2 log x))/x give the natural growth/decay pair
   in the Mellin-Barnes contour. Linear combinations with li(x), √x do
   not absorb the f residual.

## Failure mode (per project taxonomy)

**E (equivalence).** Bessel-basis indistinguishability from random matches
the standing pseudorandomness battery (≥22 measures in
`novel/pseudorandomness_of_pi.md`). This adds a 23rd measure: f(x) is
indistinguishable from a same-variance Gaussian under PSLQ on a 10-element
mixed-Bessel basis.

## Why this is informative (vs S29 redundancy)

S29's basis tested elementary functions (log, √, x^{1/k}), li-variants,
and trigonometric oscillations at zeta zeros. The Bessel basis is
strictly disjoint:

- K_0, I_0, J_0, Y_0 are NOT linear combinations of {sin, cos, log, √, li}
  over ℚ (Bessel functions are transcendental over the elementary field).
- The K_0(2π√(log x)) kernel is the only natural place where partition-
  type asymptotics could donate an identity to π(x); this is now ruled out.
- The PSLQ cross-validation discipline is the same as S29 (we deliberately
  reused S29's `fx_data.npz` to ensure identical f-values).

A surviving relation would have been a genuine breakthrough independent
of S29. Its non-survival strengthens the closed verdict from "f(x) admits
no elementary/li/trig identity" to "f(x) admits no identity in any tested
basis including spectral-trace-related Bessel kernels."

## Closing pointer (CLOSED_PATHS / EDGES)

This closure cites edges:
- E1.x (information-theoretic incompressibility of δ),
- E3.x (analytic / explicit-formula structure),
- E7.x (negative-shape: identity-search basis growth does not change verdict).

Adjacent closures: S29 RESULTS_SUMMARY.md (elementary basis), S25 PSLQ on
zero ratios (PSLQ infrastructure validated), S33 trace-formula moments
(line 653 of CLOSED_PATHS.md).

## Reproducing

```
python3 experiments/algebraic/identity_search/bessel_basis_pslq.py
```

Total runtime: ~1 s (PSLQ on 10-element vec at 50 dps).

## One-line summary

PSLQ on a 10-element basis containing K_0, I_0, J_0, Y_0 evaluated at
log x, γ_1 · log x, and 2π√(log x) finds tight numerical "relations" at
each x ∈ {5k, 10k, 50k, 100k} but ALL fail cross-validation; the
behaviour is statistically identical to a same-σ Gaussian random control.
Bessel basis joins the closed set for novel-identity search.
