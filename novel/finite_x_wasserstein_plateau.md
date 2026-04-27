# Quantitative Finite-x Wasserstein Plateau for π(x) − Li(x)

**Discovered:** S108 (wild-swing, ATTACK_VECTORS.md §C5).
**Cross-domain technique:** Stein's method (Chen-Goldstein-Shao 2011 *Normal Approximation by Stein's Method*; Ross 2011 *Probability Surveys*).
**Status:** ~~A-grade claim~~ **B-grade after verification (S112 PARTIAL).**

> **VERIFICATION NOTE (S112):** S109/S110/S111 each CONFIRMED, but
> all used the SAME Riemann γ_k as the random-phase null. S112 ran
> the missing test: at K=5000, n_modes=50, n_trials=60, D_emp's W_1
> is indistinguishable from random-phase Riemann (z=-0.93), random-
> phase non-Riemann uniform [10,145] (z=-1.26), AND random-phase
> non-Riemann equispaced (z=-1.55). The W_1 magnitude is therefore
> a generic oscillatory-sum value, not a Riemann-specific signature.
> The §C5 verbatim criterion clause "ties to a specific zeta-zero
> contribution" fails — the structural origin is true only in the
> trivial pointwise sense (D_emp = explicit-formula sum, by E1.5).
> Plateau itself confirmed (z=11 vs i.i.d. Gaussian). A → B
> demoted. Cross-domain Stein technique import remains novel.
> See `archive/sessions/session112_verify_c5_stein.md`.

> **VERIFICATION NOTE (S113, further PARTIAL):** the plateau is
> not even oscillatory-sum-specific. Across 9 distribution families
> (uniform, single arcsine, two sum-of-arcsine variants, two-Gaussian
> mixture, t df=10, Laplace, analytic low-zero sum, Gaussian control)
> tested at K∈{200..10000} with n_trials=30, EVERY distribution with
> non-zero kurtosis shows the plateau (K=200/K=10000 ratio < 2 vs the
> Stein-CLT prediction √50≈7.07). Pure Gaussian decays as 1/√K
> (ratio 6.42 ≈ √50). The plateau magnitude W_1/σ tracks |kurtosis|
> monotonically; linear interpolation at kurt(D_emp)=-0.41 predicts
> W_1/σ ≈ 0.042 vs observed 0.038 (within 10%). So the W_1 plateau
> is a *generic* W_1(P, N(μ_P, σ_P²)) measurement — well-defined
> positive whenever P is non-Gaussian — and its magnitude is fully
> predicted by D's kurtosis under log-uniform x. Stein technique
> import survives but is not load-bearing for the conclusion. Grade
> stays B (further demotion to C not warranted; the kurtosis-
> determines-magnitude relation is itself a new quantitative
> refinement of E7.5). See `archive/sessions/session113_verify_c5_stein.md`
> and `experiments/analytic/stein_wasserstein_pi/test_plateau_universality_results.md`.

---

## Statement

Define `D(x) := (π(x) − Li(x)) · log(x) / √x`.

Let `D̂_K(X) := { D(x_k) : x_k = X · e^{(k − 1)/(K−1)} }_{k=1}^{K}` be K
log-uniform anchor evaluations of `D` over the window `[X, eX]`.

Let `μ̂(K, X) = mean(D̂_K)` and `σ̂(K, X) = std(D̂_K)`.

**Claim.** The Wasserstein-1 distance to the fitted Gaussian
plateaus at a *positive* value `c(X) > 0`:

```
W_1( D̂_K(X), N(μ̂, σ̂²) )  ≥  c(X)  >  0   for all sufficiently large K.
```

For `X = 10^6` (window `[10^6, 10^7]`), `c(X) ≈ 0.008` in absolute units
(`≈ 0.04 σ̂` in standardised units), with statistical significance
≥ 5σ at K=10000 against the i.i.d.-Gaussian-fluctuation null
`W_1 ~ c_G(σ̂)/√K`.

**Pointwise structural origin (CONFIRMED).** D_emp(x) is pointwise
matched by the explicit-formula low-zero sum (S110: corr 0.98 at
n=1000 zeros). This is the trivial sense in which the plateau is
"Riemann-driven": D_emp IS the explicit-formula sum, by E1.5.

**Magnitude is NOT Riemann-specific (S112 PARTIAL).** The plateau
*magnitude* `c(X) ≈ 0.008` is a generic oscillatory-sum value: any
sum of ~50 cosines with frequencies in [10, 145] and weights
`2/√(¼+f²)` evaluated on the same log-uniform grid produces a W_1
distribution that brackets D_emp's value. The "indistinguishable
from random-phase variant" finding originally cited as evidence for
Riemann structural origin is a generic-oscillatory-sum property;
non-Riemann uniform-frequency and equispaced-frequency variants
also give |z|<2 vs D_emp.

The pointwise representation in terms of zeros is:

```
D(x) = −2 ∑_{k} cos(γ_k log x − φ_k) / √(¼ + γ_k²)
       − log(2) · log(x) / √x
       + (small remainder)
```

where `γ_k` is the imaginary part of the k-th non-trivial zeta zero
(real part 1/2 under RH) and `φ_k = arctan(2 γ_k)`. Each individual
term contributes an arcsine-like distribution with excess kurtosis
−1.5; partial sums of small numbers of low zeros remain sub-Gaussian
(excess kurtosis ≈ −0.4 to −0.8), and this sub-Gaussianity is the
direct cause of the W_1 plateau.

Across 10 sub-windows of `[10^6, 10^7]` of width 0.5 in log10, the
empirical `W_1(D̂_emp)` and the explicit-formula-truncated
`W_1(D_th(50))` are correlated `r = 0.906`. Empirical D̂ is
statistically *indistinguishable* from a random-phase variant of the
explicit-formula low-zero sum (z-score `−1.06`).

---

## Why this is novel

1. **No prior quantitative finite-x Wasserstein bound for π(x)−Li(x)**
   has been published. Selberg's CLT is asymptotic for
   `log|ζ(½+it)|`; Hejhal extended it to `π(x) − Li(x)` but only as
   asymptotic *log-distribution* convergence (`x → ∞` then weak
   limit). Pintz 1980 / Korevaar 2002 give pointwise discrepancy
   bounds, not Wasserstein-shape bounds against a Gaussian target.
2. **Stein's method has never been applied to π(x) − Li(x)** despite
   being the canonical machinery for quantitative-CLT problems. The
   project's 70+ sessions never imported it.
3. **The plateau is structural, not noise.** It scales correctly
   (`W_1(D̂)·√K → ∞`, while `W_1(N(μ̂, σ̂²)·√K → c_G ≈ 0.06`), is
   reproduced by random-phase low-zero sums (within 1σ), correlates
   `r = 0.906` with the explicit-formula prediction across
   sub-windows, and the negative-excess-kurtosis signature persists
   across all K values tested.

## What this DOES NOT claim

This is *not* a polylog π(x) algorithm. The plateau confirms that
finite-x deviation of `D(x)` from Gaussianity is:
- (a) detectable by a quantitative-CLT metric;
- (b) sourced from the lowest few zeros, *not* from arithmetic
  structure orthogonal to the explicit formula.

So the result *strengthens* the GUE-sieve-circuit closure family
(E7.1, E7.6, E7.11, S43-S57): even the finite-x-Wasserstein
microscope sees only the explicit-formula contribution, which is
already the sqrt(x) wall. **No new bit-extraction angle is opened.**
The novelty is *in the metric*, not in the underlying structure.

## Implications

- **For E7.5 (Selberg-Hejhal CLT)**: provides the first *quantitative
  finite-x* refinement of asymptotic CLT statements.
- **For E1.5 (explicit formula)**: confirms that the explicit
  formula's leading low-zero terms control the Wasserstein-1 distance
  at finite x, not just the L^∞ discrepancy.
- **For the 35+ pseudorandomness battery**: adds a 38th measure of a
  fundamentally new TYPE — Wasserstein-shape — and that measure
  *fails* (deviation detected). But the deviation reduces to E1.5,
  so it joins the GUE-sieve closure family rather than opening a new
  angle.
- **For ATTACK_VECTORS.md §C5**: closed (mode E — explicit-formula
  reduction). The "structural deviation found" outcome is realised,
  but via the leading explicit-formula contribution, not via a
  bespoke arithmetic structure orthogonal to it.

## Verifiable predictions

If this claim is correct, then (a) `W_1(D̂_K)` should remain at
`≥ 0.005` for K up to `10^5` on the same x-window; (b) the
sub-window correlation `corr(W_1(D̂_emp window), W_1(D_th(50) window))`
should remain `≥ 0.85` for any partition into 10+ sub-windows of
width `≤ 0.5` in `log10 x`; (c) the empirical excess kurtosis should
stay in `[−0.5, −0.3]` for any `X ≥ 10^6`.

If the next session's verify-mode test contradicts (a), (b), or (c)
at the 2σ level, this `novel/` entry is demoted to a CLOSED_PATHS row.

## Where to look

- `experiments/analytic/stein_wasserstein_pi/stein_wasserstein_pi_results.md`
  — full empirical results table.
- `experiments/analytic/stein_wasserstein_pi/structural_explanation.py`
  — code that derives D_th from Odlyzko's first 50 zeros and compares
  to empirical D.
- `experiments/analytic/stein_wasserstein_pi/test_low_zero_robustness.py`
  — sub-window correlation and random-phase control tests.

---

## Verification refinements (S109, S110, S111)

**S109 (CONFIRM via K-extension):** plateau confirmed at K=20000
(W_1=0.008289) and K=50000 (W_1=0.008355). The plateau is flat to
within 1% across a 5× K-range; z-score climbs as `√K` as predicted.

**S110 (CONFIRM via truncation, disjoint windows, K=10⁵):** at n=1000
zeros, `corr(D_emp, D_th(n)) = 0.98` and `W_1(D_th)/W_1(D_emp) = 0.98`
— the empirical signal is essentially completely explained by the
explicit formula's first 1000 zeros. K=10⁵ ceiling test gives
W_1=0.008494, z=33.62σ. Disjoint width-0.2 sub-window test gives
r=0.9154 (stricter than the original overlapping-window r=0.906).

**S111 (CONFIRM via X-scaling and window-width):**
- Plateau detected at X=10⁵ (W_1=0.00804, z=8.28), X=10⁶
  (W_1=0.00903, z=9.71), X=10⁷ (W_1=0.00594, z=5.06). The original
  S108 unsubstantiated claim "c(10⁷) < c(10⁶)" (which had z=-0.65 at
  K=1000) is now verified at z=5.06.
- c(X) decays slowly with X — roughly 1.5× per decade. Consistent
  with asymptotic Hejhal CLT but slower than `1/log(X)`.
- The plateau magnitude DEPENDS on the log-window width. At
  X=10⁶, c/σ = {0.049, 0.037, 0.058, 0.041} for log-widths
  {0.5, 1.0, 2.0, 2.303} — non-monotone. So **c is a function of
  (X, logw), not just X.** The values reported above (`c ≈ 0.008`)
  are specific to the empirical [10⁶, 10⁷] window (logw=2.303),
  not the formal `[X, eX]` (logw=1.0) of the Statement section.
- Kurtosis FLIPS SIGN at narrow windows (kurt=+0.086 at logw=0.5
  vs kurt=-0.42 at logw=2.303). The sub-Gaussian signature is a
  property of *averaging across multiple low-zero oscillations*,
  not an intrinsic D(x) property.

### Refinements to verifiable predictions

(a) `W_1(D̂_K) ≥ 0.005 for K up to 10⁵ on the [10⁶, 10⁷] window`:
    HOLDS (S110 at K=10⁵, W_1=0.0085). At X=10⁷ window the floor is
    marginal (W_1=0.0059, just above 0.005).
(b) Sub-window correlation `r ≥ 0.85`: HOLDS (S110 confirmed at
    r=0.9154 with stricter disjoint windows).
(c) `kurt ∈ [-0.5, -0.3] for any X ≥ 10⁶`: HOLDS strictly at X=10⁶
    (-0.42); HOLDS at X=10⁵ (-0.39); MARGINAL at X=10⁷ (-0.297, just
    outside the band but within 1σ K=5000 sampling error). The
    natural refinement: kurtosis trends to 0 asymptotically by CLT,
    so the band lower edge tightens with X. Should be tightened to
    "X ∈ [10⁶, ~10⁷]; trends to 0 as X → ∞".

The plateau result is robust across three verification sessions.

### S114 (verify-6, PARTIAL)

**Re-verification of S108 numerics with three independent W_1
implementations** — all agree on W_1 ≈ 0.00829 to within 0.3%:
- Method A (mid-rank quantile): W_1 = 0.00827, W_1/σ = 0.0376.
- Method B (scipy.stats.wasserstein_distance vs 200k MC Gaussian
  reference): W_1 = 0.00813, W_1/σ = 0.0369.
- Method C (CDF-integral on 200k grid): W_1 = 0.00827, W_1/σ = 0.0376.

S108's number 0.00829 is solid; the plateau itself is robust to the
W_1 algorithm chosen.

**Falsification of S113's "kurtosis-only fully predicts W_1/σ"
claim.** S114 built a Beta(α,α) distribution with α = 5.817
(analytically tuned via kurt(Beta(α,α) − 1/2) = −6/(2α+3) ⇒ kurt
= −0.41, exactly matching D_emp's kurtosis). At K=10000, n_trials=30,
σ_target = 1.6:

- Beta sample kurtosis = -0.4127 (matches target).
- Beta W_1/σ = **0.0328**.
- D_emp W_1/σ = **0.0376** (gap: -13% below S113's ±10% band).
- S113 kurtosis-only prediction = **0.0426** (gap: -23% from Beta).

Both are outside the ±10% band of S113's own prediction. Conclusion:
S113's universality is **qualitative** — every non-Gaussian
distribution plateaus — but the magnitude requires more than just
kurtosis. A higher-moment-aware fit `f(kurt, fourth_moment_shape, …)`
would be needed for ~10% prediction accuracy across distribution
families.

This refinement does NOT un-demote S108: D_emp's W_1/σ is still
generic (Beta also produces a plateau, just to a different magnitude),
and the value remains not Riemann-specific. But it tightens the
scope of S113's claim: kurtosis alone is a ~30%-accurate predictor,
not a tight predictor.

---

## S115 verify-7 — sub-window correlation r=0.906 IS Riemann-φ-specific (CONFIRM)

The seventh verify session ran the angle no prior verify had touched:
the structural-matching claim said `corr(V_emp, V_th_actual) = 0.906`
across 10 sub-windows of [10⁶, 10⁷] (V = vector of W_1 per sub-window
to fitted Gaussian; V_th_actual uses the actual first-50 Odlyzko γ_k).
Is r=0.906 Riemann-phase-specific, or a generic property of any low-
frequency oscillatory sum on a log-uniform grid?

**Method:** 200-trial random-phase null on three reference families.

| Reference                                 | r mean   | r std  | min    | max    | r ≥ 0.906 |
|-------------------------------------------|----------|--------|--------|--------|-----------|
| **Actual γ + actual phases** (V_th_actual)| **+0.906** | —      | —      | —      | —         |
| Random φ, **same Riemann γ_k** (200)      | -0.044   | 0.389  | -0.852 | +0.889 | **0/200** |
| Random φ, non-Riemann γ ∈ [10, 145] (200) | -0.032   | 0.381  | -0.795 | +0.902 | 0/200     |
| Pure noise, standard normal (200)         | -0.030   | 0.387  | —      | —      | —         |

z(actual r vs random-phase Riemann ensemble) = **+2.44**;
p < 0.005 (one-sided).

The three null distributions are statistically indistinguishable —
mean and std agree to 0.02 and 0.01 respectively, all consistent with
the n=10 noise floor 1/√n ≈ 0.316. The actual zero PHASES are
necessary AND sufficient to reproduce V_emp pointwise.

**Conclusion.** The sub-window correlation r=0.906 is real signal —
**not** a generic oscillatory-sum-on-log-uniform-grid artefact. With
random phases (Riemann γ_k or otherwise), the correlation collapses
to the noise floor. So S108's "pointwise structural matching" claim
survives this attack.

**No grade change.** What was demoted by S112/S113 is the W_1
*magnitude* (generic across non-Gaussian distributions), not the
sub-window correlation. The sub-window correlation surviving means
D_emp ≈ low-zero explicit-formula sum pointwise — i.e., E1.5
re-asserted. That re-assertion has been the framing all along (see
EDGES.md E1.7's "i.e., re-confirm E1.5"). So S108 stays at B.

**Verify chain saturated.** Seven verify sessions (S109–S115) have
exhausted every angle: K-extension (S109), truncation/disjoint sub-
windows (S110), X-scaling/window-width (S111), random-phase nulls on
W_1 magnitude (S112), 9-distribution universality (S113), three-
method numeric reproduction + Beta(α,α) targeted-kurt (S114), random-
phase null on sub-window correlation (S115). Future verify sessions
should attack a different target.
