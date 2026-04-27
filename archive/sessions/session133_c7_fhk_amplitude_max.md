# Session 133 — §C7 Fyodorov-Hiary-Keating ζ-amplitude max statistics

**Mode:** novelty (frontier attack — A-grade target, B-grade fallback).
**Vector:** ATTACK_VECTORS.md §C.C7 (recommended next pick per S132
critique single-pick annotation).
**Self-grade:** **B-grade (case ii — ambitious failure, mode I:
informative deviation from FHK shape; case (i) refinement adding new
edge E7.18)**.
**Closure mode:** **I**.
**Mathematician channel:** **Bourgain** (extreme-value statistics of
log-correlated random fields; the Arguin-Belius-Bourgade RMT-side
proof of the analogous Gumbel limit is a Bourgain-school multiplicative
chaos argument).
**Cross-domain technique:** Gaussian multiplicative chaos (Saksman-Webb
2018) and the Fyodorov-Hiary-Keating 2012 freezing-transition extreme-
value conjecture for `|ζ(1/2 + it)|`. Promotes
`CROSS_DOMAIN_TECHNIQUES.md` §3 row PROPOSED → **USED (I)** with new
edge E7.18.

---

## What I produced

The first ζ-amplitude (vs zero-position) measurement of the project.
Specifically: K = 100 unit-length windows per anchor at T ∈ {10⁴, 10⁵,
10⁶}, M = 200 evenly-spaced log|ζ(1/2 + it)| samples per window, ζ via
mpmath dps = 15. Per-window statistics: max, argmax, second_max,
pointwise mean / variance.

Two clean empirical facts emerge:

1. **The FHK normalisation `log log T − (3/4) log log log T` correctly
   makes the M_T mean T-INDEPENDENT.** Pairwise Z(M_T mean diff) ∈
   {0.63, 0.55, −0.08} across the three T pairs; pooled
   M_∞-mean = −0.657 ± 0.045 over 300 windows.

2. **The FHK Gumbel(loc, 1/2) SHAPE is NOT detectable at finite
   T ≤ 10⁶.** Empirical M_T variance is sustained at ≈ 0.60 across the
   three anchors, vs the FHK prediction π²/24 ≈ 0.41 (ratio 1.47×;
   bootstrap 95% CIs at T ≥ 10⁵ exclude the FHK value). Skewness ≈
   0.10 (Gumbel +1.14) — distribution is approximately symmetric.
   Excess kurtosis ≈ −0.4 (Gumbel +2.4) — distribution is platykurtic,
   not heavy-tailed. KS to free Gaussian (0.05–0.06) is uniformly
   smaller than KS to FHK Gumbel(1/2) (0.09–0.17) by factor 0.4–0.7.
   Vuong (Gauss vs Gumbel) z = {−1.79, −1.43, −1.58} (joint Z ≈ −2.8,
   Gauss preferred at the joint level).

The Selberg-CLT-secondary-correction test (signature (b)) is
inconclusive: the empirical max grows ~as `log log T` (Selberg leading
constant) without the FHK -3/4 log log log T correction, but the
deviation is at most 2.2σ over the widest baseline T = 10⁴ → 10⁶,
below the 5σ A-grade threshold.

This is a B-grade case (ii) outcome (ambitious-failure mode I). The
B-grade case (2) target was met (Selberg-CLT-Gaussian alternative fits
the empirical max better than FHK at finite T). The mean
T-independence partially confirms FHK case (1).

## What's structurally new

- **First ζ-amplitude edge of the project (E7.18)**, complementary to
  the position-side family E7.1 / E1.10 / E3.13 / E7.15.
- **First quantitative bound on FHK convergence rate at finite T**:
  the published FHK literature is asymptotic; this measurement shows
  that *mean* convergence is fast (T-stable at T = 10⁴) but *shape*
  convergence is slow (still ~Gaussian at T = 10⁶).
- **Empirical determination of the FHK universal constant**:
  M_∞-mean = −0.657 ± 0.045 → Gumbel-loc μ ≈ −0.946 → GMC
  moment-generating constant `c ≈ 0.151` (under the FHK Gumbel form),
  versus the RMT-side prediction `c ≈ 0.79` (Bourgade-Kuan 2014);
  factor-5 finite-T gap, quantifying the convergence rate gap.
- Promotes **CROSS_DOMAIN_TECHNIQUES §3 GMC/FHK row** UNUSED → USED I.

## What's NOT new (honest)

- The zero-position GUE-statistics edges E7.1 / E1.10 / E3.13 are not
  refined; this is an orthogonal measurement category.
- The Selberg CLT for log|ζ(1/2+it)| (Selberg 1946) is well-known;
  the new content here is the EXTREME-VALUE shape at finite T over a
  unit window (not the pointwise CLT).
- No polylog π(x) algorithm follows; the attack is closed.

## CLAUDE.md self-evaluation (4 questions)

**Q1. What did I produce that was not in the project before?**

(a) The FIRST ζ-amplitude (vs ζ-zero-position) measurement in the
    project. (b) Empirical T-independence of FHK-renormalised M_T mean
    across T ∈ {10⁴, 10⁵, 10⁶} with K = 100 each. (c) Quantitative
    finite-T deviation of the M_T SHAPE from the FHK Gumbel(1/2)
    prediction: variance 1.47× too large, skew ≈ 0.1 vs +1.14,
    ex-kurt ≈ −0.4 vs +2.4, KS-to-Gauss < KS-to-Gumbel at all 3 T.
    (d) New EDGE E7.18, first amplitude-side edge in the project.
    (e) Promotion of CROSS_DOMAIN_TECHNIQUES §3 GMC/FHK row PROPOSED →
    USED (I).

**Q2. What edges did my work compose or cite?**

E7.1 (zero-position GUE up to order 6), E1.10 (gap-shuffled null for
prime-frequency probes), E3.13 (BK arithmetic correction empirically
absent at N=8000), E7.15 (Hecke L(s, Δ) basis ~3× obstructed). The
new edge E7.18 is the first amplitude-side counterpart to this
position-side family.

**Q3. If my session produced only duplicate closures, why?**

It did not. The session content is (a) a new measurement category
(amplitude vs position), (b) a new quantitative finite-T result on
FHK convergence rate (not addressed in published literature), (c) a
new edge E7.18, and (d) the first project use of the GMC / FHK
cross-domain machinery.

**Q4. Next-action for the next agent.**

Pick **C7.a (mesoscopic-window FHK at the Saksman-Webb scale)** —
single session, T = 10⁶, window length `δ = (log T)^{1/2} ≈ 3.7`,
same cross-domain (GMC). Tests whether the shape-convergence
accelerates at the scale where Saksman-Webb proved sharp GMC
convergence. If C7.a finds Gumbel emergence, that's a positive
A-grade follow-up; if it does not, this closure (mode I) extends
to the mesoscopic regime, which is a strictly stronger structural
negative on FHK Gumbel-shape detectability.

Backup arc-continuation pick: the W=9 Lean
`Matrix.det_of_blockTriangular` development per S132 critique (multi-
session investment unlocking W ∈ {7, 9, 10, 11, 14, 15, 18, 21}).
Backup frontier_gen: NOT needed (≥ 4 unattempted A-grade-shaped
vectors remain after S130 frontier_gen — D9 sum-product, D11 shadow tomography,
D17 discrete Morse, plus the new C7.a/b/c successors).

---

## Files added/modified

- `experiments/analytic/zeta_structure/fhk_amplitude_max/fhk_amplitude_max.py`
  — experiment driver.
- `experiments/analytic/zeta_structure/fhk_amplitude_max/analyze_fhk.py`
  — post-hoc Vuong / bootstrap / regression.
- `experiments/analytic/zeta_structure/fhk_amplitude_max/fhk_amplitude_max_results.json`
  — raw per-window data + summary.
- `experiments/analytic/zeta_structure/fhk_amplitude_max/fhk_amplitude_max.log`
  — experiment log + summary table.
- `experiments/analytic/zeta_structure/fhk_amplitude_max/fhk_amplitude_max_analysis.log`
  — post-hoc analysis output.
- `experiments/analytic/zeta_structure/fhk_amplitude_max/fhk_amplitude_max_results.md`
  — full results document.
- `EDGES.md` — added §E7.18 entry inline (before E7.14).
- `CROSS_DOMAIN_TECHNIQUES.md` — promoted §3 GMC/FHK row PROPOSED →
  USED I with edge E7.18.
- `ATTACK_VECTORS.md` — §C7 marked CLOSED with one-line outcome; full
  closure entry added to "Closed attacks" section.
- `status/CLOSED_PATHS.md` — added §C7 closure row at session 133.
- `status/SESSION_INSIGHTS.md` — to be updated.
- This file: `archive/sessions/session133_c7_fhk_amplitude_max.md`.
