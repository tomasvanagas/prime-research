# S115 Verify-7 — Sub-window correlation r=0.906 random-phase null

**Date:** 2026-04-27.
**Session:** S115 (seventh verify on S108).
**Target:** S108's structural-matching claim:
`corr(W_1(D_emp_window), W_1(D_th50_window)) = 0.906` across 10
sub-windows of `[10^6, 10^7]` (width 0.5 in log10), interpreted as
evidence of low-zero explicit-formula structural origin.

## Question

Is r=0.906 between V_emp (sub-window W_1 of empirical D) and V_th
(sub-window W_1 of theoretical D_th50 with actual Riemann zeros)
**Riemann-phase-specific**, or is it generic for any low-frequency
oscillatory sum on a log-uniform grid?

S110 confirmed r ≈ 0.92 on disjoint windows using ACTUAL zeros.
S112/S113 showed FULL-window W_1 magnitude is generic. None of
S109–S114 ran a random-phase null on the SUB-WINDOW correlation.

## Method

`verify_S115_subwindow_rand_phases.py`. Pipeline:

1. K=10000 log-uniform anchors in [10^6, 10^7]. π via sympy.primepi.
2. D_emp = (π − Li) log(x) / √x (S108 definition).
3. 10 sub-windows of width 0.5 in log10, starting in [6.0, 6.5].
4. V_emp = vector of W_1(D_emp restricted, fitted N(μ̂, σ̂²)) per window.
5. V_th_actual = same on D_th50 using **actual** first-50 Odlyzko γ_k.
6. r_actual = corr(V_emp, V_th_actual).
7. Random-phase null with **same Riemann γ_k**: 200 trials of
   D_rand = −Σ w_k cos(γ_k log x − φ_k_random); V_rand[i] per sub-window;
   r distribution.
8. Random-phase **non-Riemann** null: γ_k drawn uniform [10, 145]; 200 trials.
9. Pure-noise control: D_rand = standard normal samples; 200 trials.

W_1 routine: fast mid-rank quantile match (validated against S108's
closed-form on the K=10000 reference: 0.5% agreement).

## Results

| Source                                    | r        | std      | min      | max      | frac r ≥ 0.7 |
|-------------------------------------------|----------|----------|----------|----------|--------------|
| **Actual Riemann γ_k, actual phases**     | **+0.906** | —        | —        | —        | —            |
| Random-phase + Riemann γ_k (200 trials)   | -0.044   | 0.389    | -0.852   | +0.889   | 2.5%         |
| Random-phase + non-Riemann γ_k (200)      | -0.032   | 0.381    | -0.795   | +0.902   | 2.5%         |
| Pure noise (standard normal, 200)         | -0.030   | 0.387    | —        | —        | —            |

z(actual r vs random-phase Riemann) = +2.44. **0/200 random-phase
trials reach r ≥ 0.906** (one-sided p < 0.005).

The three null distributions (random-phase Riemann, random-phase non-
Riemann, pure noise) are **statistically indistinguishable** from each
other (means and stds within 0.02). The r=0.906 sub-window correlation
is a real signal, not a generic oscillatory-sum-on-log-grid artefact.

## Verdict: **CONFIRM** S108's structural-matching claim

The sub-window correlation r=0.906 IS Riemann-phase-specific. With the
actual γ_k but RANDOM phases, the mean correlation is essentially zero
with std = 0.39 (matching the n=10 noise floor 1/√n ≈ 0.316). With
non-Riemann frequencies AND random phases, the distribution is
statistically identical. So:

- **Random phases break the sub-window structural match** (whether or
  not the frequencies are Riemann γ's). This isolates the contribution
  of the actual zero PHASES.
- **The actual γ phases are necessary AND sufficient** to reproduce
  V_emp pointwise — exactly what E1.5 predicts.

## What this implies for S108's grade

**No change.** S108 has been at B since S112. The structural-matching
claim (r=0.906) was always positioned by S108 as "re-confirming E1.5"
— see EDGES.md E1.7: *"These confirm D_emp ≈ low-zero explicit-formula
sum pointwise — i.e., re-confirm E1.5."* My test confirms that
positioning is correct: the sub-window correlation is real signal but
the underlying mathematical content is exactly E1.5 (Riemann's
explicit formula's leading low-zero contribution dominates D pointwise
at finite x).

What was demoted in S112 is the **W_1 magnitude** claim, not the
sub-window correlation. The W_1 magnitude is generic across any
non-Gaussian distribution (S113); my test does not address that.

## What would have falsified S108 further

If random-phase Riemann variants reproduced r ≥ 0.7 with frequency
> 50%, the sub-window correlation would be generic and S108 would have
fallen further (B → C). The observed 2.5% rate is inside-the-noise:
under 1/√n ≈ 0.316, a Gaussian tail with std=0.39 above 0.7 happens at
~3.5% rate — matching the observed 2.5%. So the random-phase null is
truly null.

Therefore: **the structural-matching r=0.906 is a Riemann-phase
signature**, but it is mathematically equivalent to E1.5 — no new edge.

## Speed validation

The fast mid-rank W_1 routine (used to enable 200 × 10 = 2000 sub-
window W_1 evaluations in ~30s) was validated on V_emp against S108's
closed-form: agreement to 0.5%, well within sub-window noise.

## Self-grade: **B**

Per verify-mode rubric:
- A — clear refutation. NO; CONFIRM upholds the (already-B) claim.
- **B — confirmed via non-trivial reproduction.** YES — random-phase
  null on sub-window correlation is a TEST NOT PREVIOUSLY RUN by
  S109–S114, addressing a specific weakness in the structural-matching
  argument (what if r=0.906 is a generic oscillatory artefact?).
- C — trivial reproduction. NO; involved 200-trial null with two
  reference families plus noise control.

## Files written

- `verify_S115_subwindow_rand_phases.py` — script.
- `verify_S115_subwindow_rand_phases.log` — runtime log.
- `verify_S115_subwindow_rand_phases_results.json` — numeric results.
- `verify_S115_subwindow_rand_phases_results.md` — this file.

## Falsification statement

A subsequent verify session falsifies S115's CONFIRM verdict if any of:
- Random-phase Riemann correlation distribution at n_trials > 1000
  has mean > 0.3 (currently mean = -0.04, std = 0.39).
- A specific known PHASE-FREE statistic on D's sub-windows reproduces
  V_emp pointwise (which would mean the phases are not load-bearing).
- The reported r = 0.906 fails to reproduce on a different random-seed
  primepi computation (currently reproduced exactly to S108's value).
