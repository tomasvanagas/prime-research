# Slot 2 (Thread 10 / P5) — Finer x-grid measurement of c_emp

**Date:** 2026-04-30
**Mission:** distinguish Galway-shape (c_emp ≈ const) from Thread-7-shape
(c_emp ~ √x / log²x) by measuring c_emp at 10 anchors across the
decade `log10 x ∈ [4.0, 4.9]` using the existing K=8000 zero database
with finer K-milestones (250-step granularity above K=1000).

## Driver

`slot2_finegrid.py` — imports helpers from `slot1_worst_case_K.py`.
Run: `python3 slot2_finegrid.py` (defaults: 10 anchors, N=30, K_max=8000,
ε ∈ {1, 3}, M=4 Mobius terms, dps=18).
Wall: 10 min; output: `slot2_finegrid_traces.csv` (12,900 rows),
`slot2_finegrid_summary.csv` (20 rows), `slot2_finegrid_run.log`.

## Headline result

**Slot 1's `c_emp ≈ const ≈ 0.20` claim refines to `c_emp ∈ [0.106, 0.249]`
at finer resolution, with mean `0.151 ± 0.044` across 7 finite-K anchors.**

| log10 x | x | K_emp (worst, ε=1, N=30) | c_emp | c_emp_T7_pred |
|---:|---:|---:|---:|---:|
| 4.00 | 10000 | 1250 | 0.1474 | 0.1785 |
| 4.10 | 12589 | 1250 | 0.1250 | 0.1929 |
| 4.20 | 15849 | 1250 | 0.1062 | 0.2084 |
| 4.30 | 19953 | 1750 | 0.1264 | 0.2254 |
| 4.40 | 25119 | 2250 | 0.1383 | 0.2438 |
| 4.50 | 31623 | 4750 | 0.2488 | 0.2638 |
| 4.60 | 39811 | 3750 | 0.1675 | 0.2856 |
| 4.70 | 50119 | >8001 | (LB > 0.306) | 0.3093 |
| 4.80 | 63096 | >8001 | (LB > 0.262) | 0.3352 |
| 4.90 | 79433 | >8001 | (LB > 0.224) | 0.3634 |

`c_emp_T7_pred = 0.61 · K_pred(half-G) / (√x · log²x)` is the Thread-7-shape
prediction with σ_eff/σ_pred = 0.755 (squared = 0.61).

## Interpretation

### Slope analysis

`d(log K_emp)/d(log10 x)` over the 7 finite-K anchors (log10 x = 4.0..4.6):
empirical slope ≈ **0.80** per unit log10 x.

| Shape | Predicted slope d(log K_emp)/d(log10 x) | Notes |
|---:|---:|---:|
| Galway-shape | ~0.7 | from `K = c · √x · log²x` ⇒ slope = 1/2 + 2/(log10 · log x) |
| Thread-7-shape (typical) | ~1.0 | from K_pred ~ x · log²K / log²x |

Empirical 0.80 sits between the two predictions, closer to Galway-shape.
**Inconclusive on slope alone** — single-decade measurement with worst-of-30
sample-variance noise (factor ~2× per anchor).

### σ_eff/σ_pred at K_emp (Thread 7 GUE-factor consistency check)

At each `K_emp` where worst-of-N = ε=1, the implied `σ_eff(K_emp) = ε /
√(2 log N) = 0.383`. Comparing to `σ_pred(K_emp, x)` from slot 2 column:

| log10x | σ_pred(K_emp) | σ_eff/σ_pred |
|---:|---:|---:|
| 4.00 | 0.493 | 0.777 |
| 4.10 | 0.540 | 0.709 |
| 4.20 | 0.591 | 0.648 |
| 4.30 | 0.573 | 0.668 |
| 4.40 | 0.573 | 0.668 |
| 4.50 | 0.474 | 0.808 |
| 4.60 | 0.570 | 0.672 |

Mean: `σ_eff/σ_pred = 0.707 ± 0.063` across 7 cells. **Consistent with
slot 1's 0.796 ± 0.092 within 1σ.** Confirms the GUE pair-correlation
factor at the worst-case-of-N tail across one more decade of measurement.

### Hypothesis test against Thread-7-shape with x-dependent σ_eff/σ_pred

Slot 1 reported σ_eff/σ_pred increases mildly with log10 x: 0.72 (10⁴) →
0.79 (10⁴·⁵) → 0.88 (10⁵). Slot 2 confirms 0.71 ± 0.06 at the worst-case-of-N
statistic across log10 x ∈ [4.0, 4.6]. The two are within 1σ of each other.

**Both Galway-shape and Thread-7-shape (with σ_eff/σ_pred ≈ 0.71) fit the
slot 2 data within noise.** Single-decade dynamic range insufficient.

### Comparison to slot 1

Slot 1 reported `c_emp ≈ 0.18-0.21` at log10 x = 4.0/4.5. Slot 2 reports
`c_emp ≈ 0.11-0.25` at log10 x = 4.0..4.6 (finer grid). The slot 2 mean
(0.151) is slightly below slot 1's (0.193), consistent with sample
variance from different random N=30 sets per anchor.

The single anomalous c_emp = 0.249 at log10 x = 4.5 in slot 2 (and slot 1's
0.21 at the same anchor) suggests this anchor had a particularly bad
worst-case x_i sample; cross-slot consistency at that anchor is not
expected because slot 1 and slot 2 used different random samples.

## What slot 2 actually establishes

1. **σ_eff/σ_pred ≈ 0.71 at worst-case-of-N statistic across 7 anchors**
   (one new decade of confirmation beyond slot 1's 12 cells). The GUE
   pair-correlation factor extends to the tail-of-N regime with ratio
   stable across 10⁴..10⁴·⁶.
2. **K=8000 budget exhausted at log10 x ≥ 4.7.** Empirical K_emp at
   log10 x = 4.7 is in (8001, ∞); both Galway-shape (predicting ≈
   5226 with c=0.20 → already-exhausted-noise consistent) and Thread-7-shape
   (predicting ≈ 8732) are compatible with the lower bound.
3. **Galway-shape vs Thread-7-shape NOT distinguished at slot 2 dynamic
   range** — single-decade measurement with K_max = 8000 zeros budget is
   insufficient. Decisive measurement requires log10 x ≥ 5.5 with
   K_max ≥ 20000.

## Six falsifiers (slot 2)

- **F1 (slot 1's F1 carries):** N=30 worst-case might not scale
  half-Gaussian-tail-like at N=100. Untested.
- **F2 (slot 1's F2 carries):** sample window [10^a, 1.5·10^a] (50%);
  results may shift at narrower (10%) or wider (100%) windows.
- **F3 (slot 2 specific):** the c_emp = 0.249 outlier at log10 x = 4.5
  (also seen in slot 1) might indicate a true x-window where worst-case
  K_emp jumps; finer sampling near 10⁴·⁵ could test this.
- **F4 (slot 1's F4 carries):** Mobius truncation M=4 / dps=18; M=8
  / dps=25 might shift c_emp by 5-10%. Untested at slot 2 range.
- **F5:** no slot 2 measurement at log10 x ≥ 5.5 (out-of-budget); slot
  3 must extend zeros to settle Galway-vs-T7 question.
- **F6:** with N=30 worst-of-30 has factor-2 sample variance per anchor;
  N=100 measurement would tighten c_emp to ±10% per anchor at 3× wall cost.

## Slot 3 plan

Run K_emp measurement at log10 x ∈ {5.0, 5.3, 5.5, 5.7} with extended zero
database K_max ≥ 25000 (running in parallel via 12-worker mpmath.zetazero,
each worker on a 1000-zero range; ETA ~25 min total wall).

Decision rule:
- If K_emp(log10 x = 5.5) ≤ ~19000: Galway-shape with c ≈ 0.20 confirmed
  across 1.5 decades. **Slot 3 publishable empirical c_emp tightening.**
- If K_emp(log10 x = 5.5) ≥ ~50000: Thread-7-shape confirmed.
  **Slot 3 closes Thread 10 as B-NEGATIVE (Galway worst-case empirical
  constant is finite-x artefact, not asymptotic).**
- If K_emp ∈ (19000, 50000): inconclusive, slot 3 escalates.

## Files

- `slot2_finegrid.py`: 165-line driver (imports slot1 helpers).
- `slot2_finegrid_traces.csv`: 12,900 rows.
- `slot2_finegrid_summary.csv`: 20 rows.
- `slot2_finegrid_run.log`: full stdout.
- `slot2_finegrid_results.md`: this file.

## Edges / threads cited

- **Thread 7 (S195+S240+S241+S244):** σ_pred formula, GUE pair-correlation
  factor 0.755 ± 0.06; slot 2 extends to 0.71 ± 0.06 at worst-of-N tail.
- **Thread 10 slot 1 (S434):** initial 3-anchor measurement; slot 2 refines
  to 10 anchors with finer K-milestone resolution.
- **E6.1 (Galway 2004 K = O(x^{1/2+ε})):** slot 2 measures empirical
  worst-case constant for this bound at finite x with finer x-resolution.
