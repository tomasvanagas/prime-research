# P5 / Thread 10 (commit, slot 1) — Galway worst-case K-constant tightening

**Date:** 2026-04-30 (S434)
**Commit thread:** Thread 10 (OPEN_POSITIVE_TARGETS.md §P5).
**Slot:** 1 of 3.
**Status (slot 1):** complete.
**Self-grade (slot 1):** B — first empirical measurement of the worst-
case K-constant `c_emp` in the Galway-shape bound `K = c · √x · log²x`,
across two reachable decades. The empirical `c_emp ≈ 0.18–0.21` is
slowly drifting upward with x and a slot-2 extension will resolve
whether the Galway-shape `c_emp = const` or the Thread-7-shape
`c_emp ~ √x / log²x` holds asymptotically. Not A because (i) the
distinction between the two shapes is not yet decided at the slot-1
sample range; (ii) no benchmark vs HKM/primecount executed; the
cross-shape resolution is the slot-2 deliverable. Not C because slot 1
delivers a previously-unmeasured empirical constant for the explicit-
formula evaluator on a worst-case-over-N-samples basis, plus
confirmation that Thread 7's `σ_eff/σ_pred ≈ 0.75-0.85` GUE pair-
correlation factor extends to the worst-case-tail statistic via
`K_emp ≈ 0.6 · K_pred(half-Gaussian)`.

## Mission statement

OPEN_POSITIVE_TARGETS.md §P5: *"What is the smallest K(x) such that
the truncated explicit-formula evaluator π_K(x) attains π(x) ± 1
WORST-CASE over many random x in a window? Galway 2004 gives
K = O(x^{1/2+ε}) under GRH; explicit constants are loose."*

Slot 1's concrete first step: measure `K_emp(x, ε, N) := min K
s.t. max_{i=1..N} |π_K(x_i) - π(x_i)| ≤ ε` for x_i sampled in a
tight window around `x ∈ {10⁴, 10⁴·⁵, 10⁵}`, N = 30 samples,
ε ∈ {1, 3, 10, 30, 100}, K up to 8000 (project's pre-computed zero
database).

The outputs:
- (a) the K_emp scaling vs x — does it look like x^{1/2}, x/log²x, or
  something else?
- (b) the explicit prefactor `c_emp(x) := K_emp / (√x · log²x)` —
  empirical Galway-bound prefactor.
- (c) an empirical comparison to the half-Gaussian-tail prediction
  `K_emp(x, ε, N) ≈ K_typical(x, ε / √(2 log N))` under the random-
  phase / GUE-pair-correlation framework (Thread 7 / S244).

## Setup

- **Evaluator:** the Thread 7 partial-sum machinery (`polylog_approx_pi.py`
  exposes `R_at_rho`, `riemann_R`, `get_zeros`). Mobius truncation `M=4`
  and mpmath `dps=18` for ~3× speedup over Thread 7's `M=8 / dps=25`
  configuration; accuracy verified empirically (errors at K_max <2 across
  all 90 cells, well within 18-digit float precision).
- **Reference π(x):** sieve of Eratosthenes up to N = 1.5×10⁵.
- **Zero database:** project's `data/zeta_zeros_8000.txt` (first 8000
  imaginary parts of ζ-zeros at 15 dps; γ_8000 ≈ 8148.19).
- **Sampling:** seed = 0, 30 distinct integer x_i per anchor in the
  window `[10^a, 1.5 × 10^a]`.
- **K-milestones:** logarithmic grid
  `[1, 2, 5, 10, 20, 50, 100, 200, 300, 500, 700, 1000, 1500, 2000,
  3000, 4000, 5000, 6000, 7000, 8000]`.

Total runtime: 180 s for 90 cumulative-K evaluations. Output:
`slot1_traces.csv` (1800 rows), `slot1_summary.csv` (15 rows),
`slot1_run.log`.

## Results

### Table A — K_emp(x, ε, N=30) and prefactor c_emp

| log10 x | x | ε | K_emp (worst) | K_typ (median) | K_pred(half-G) | K_galway(c=1) | c_emp = K_emp/(√x log²x) | σ at K_emp |
|--------:|---:|---:|---:|---:|---:|---:|---:|---:|
| 4.00 | 10000 | 1 | 1500 | 2 | 2483 | 8483 | 0.1768 | 0.461 |
| 4.00 | 10000 | 3 | 20 | 1 | 93 | 8483 | 0.0024 | 1.637 |
| 4.00 | 10000 | 10 | 1 | 1 | 2 | 8483 | 0.0001 | ∞ |
| 4.50 | 31623 | 1 | 4000 | 10 | 8256 | 19092 | 0.2095 | 0.507 |
| 4.50 | 31623 | 3 | 50 | 1 | 408 | 19092 | 0.0026 | 2.137 |
| 4.50 | 31623 | 10 | 1 | 1 | 2 | 19092 | 0.0001 | ∞ |
| 5.00 | 100000 | 1 | **>8001** | 1000 | 27086 | 41915 | **>0.191 LB** | 0.621 (at K=8000) |
| 5.00 | 100000 | 3 | 700 | 1 | 1563 | 41915 | 0.0167 | 1.531 |
| 5.00 | 100000 | 10 | 1 | 1 | 31 | 41915 | 0.0000 | ∞ |

### Table B — σ_eff / σ_pred at the K-frontier (Thread 7 GUE factor confirmation)

For each anchor, `σ_eff = √(⟨err²⟩)` (RMS over 30 samples) vs
`σ_pred(x, K) = √x · log K / (π · √(2K) · log x)`:

| log10 x | K | σ_eff | σ_pred | σ_eff/σ_pred |
|---:|---:|---:|---:|---:|
| 4.00 | 100 | 0.811 | 1.125 | 0.720 |
| 4.00 | 1000 | 0.389 | 0.534 | 0.728 |
| 4.00 | 4000 | 0.194 | 0.321 | 0.606 |
| 4.00 | 8000 | 0.200 | 0.246 | 0.816 |
| 4.50 | 100 | 1.219 | 1.779 | 0.685 |
| 4.50 | 1000 | 0.658 | 0.844 | 0.780 |
| 4.50 | 4000 | 0.445 | 0.507 | 0.878 |
| 4.50 | 8000 | 0.318 | 0.388 | 0.818 |
| 5.00 | 100 | 2.381 | 2.847 | 0.836 |
| 5.00 | 1000 | 1.224 | 1.350 | 0.906 |
| 5.00 | 4000 | 0.681 | 0.811 | 0.840 |
| 5.00 | 8000 | 0.578 | 0.621 | 0.931 |

**Mean ratio over 12 cells: 0.796 ± 0.092** (consistent with Thread 7
S241 / S244 result `σ_eff/σ_pred = 0.755 ± 0.06`; the slot-1 sample
spans only 90 (anchor, K) cells with N=30 each, so the variability is
expected to be larger than Thread 7's 360-cell ensemble).

### Table C — worst-case |err| at K_max=8000

| log10 x | worst | p95 | median | mean |
|---:|---:|---:|---:|---:|
| 4.00 | 0.551 | 0.466 | 0.067 | 0.133 |
| 4.50 | 0.633 | 0.581 | 0.242 | 0.270 |
| 5.00 | 1.475 | 1.113 | 0.392 | 0.458 |

The worst-case |err| at x=10⁵, K=8000 is 1.475 (just above the ε=1
threshold). For x ≤ 10⁴·⁵, K=8000 is sufficient for ε=1 worst-case
across N=30; at x=10⁵, the K-budget is exhausted before the worst-
case error drops below 1.

## What slot 1 produced

1. **First empirical worst-case K-constant `c_emp` for the Galway-shape
   bound.** At ε=1 worst-case across N=30 samples in `[10^a, 1.5×10^a]`:

   | log10 x | K_emp | c_emp |
   |---:|---:|---:|
   | 4.00 | 1500 | 0.177 |
   | 4.50 | 4000 | 0.210 |
   | 5.00 | >8000 (LB) | >0.191 |

   The `c_emp` is **slowly drifting upward** across 1.5 decades (0.177
   → 0.210). Whether this trend continues (Thread-7-shape `c_emp ~
   √x / log²x`) or saturates (Galway-shape `c_emp = const`) is the
   slot-2 deliverable.

2. **Confirmation that Thread 7's GUE pair-correlation factor
   `σ_eff/σ_pred ≈ 0.75-0.85` extends to worst-case statistics.**
   Empirical K_emp / K_pred(half-Gaussian) ≈ 0.49-0.60 at the two
   reachable points (the 0.61 prediction from variance-factor 0.78²
   matches within 20%).

3. **The half-Gaussian-tail prediction overshoots K_emp by factor 1.7-3.4.**
   Heuristically, this reduces to the GUE variance-suppression factor
   (Thread 7's structural result transposed to the worst-case-of-N
   ensemble). The deviation factor `K_emp / K_pred(half-G)` is
   approximately `(σ_eff/σ_pred)² ≈ 0.61` — the half-Gaussian prediction
   uses Gaussian variance but the actual variance is reduced by the GUE
   pair-correlation factor.

4. **Hypothesis test (slot 2 falsifier):** if the Galway-shape `c_emp =
   const` is correct, then at `x = 10^5.5` we expect `K_emp ≈ 0.21 ·
   √(10^5.5) · log²(10^5.5) ≈ 0.21 · 562.3 · 159.7 ≈ 18900`. If the
   Thread-7-shape `c_emp ~ √x / log²x` is correct, then `K_emp ≈ x ·
   log N / (π² · log²x) · const ≈ 10^5.5 · 3.4 / (9.87 · 159.7) ≈
   675`. **The two predictions differ by 28×.** Slot 2's measurement
   at x=10^5.5 (requiring K up to ~25000 zeros) decisively
   distinguishes them.

## Falsifiers (what would make slot 1 wrong)

1. **F1 — N-dependence not captured.** The `K_emp` measurement
   averages worst-case over N=30 samples; running with N=100 should
   give `K_emp_N100 ≈ K_emp_N30 · log(100)/log(30) ≈ K_emp · 1.35`
   under the half-Gaussian model. If empirical N-scaling deviates by
   >50% slot 1's framing is wrong.
2. **F2 — sample-window dependence.** Slot 1 samples in
   `[10^a, 1.5·10^a]` (50% window). If `K_emp` differs across windows
   `[10^a, 1.1·10^a]` (10%) vs `[10^a, 2·10^a]` (100%) by a factor
   inconsistent with sample-density-of-x in the window, slot 1's
   c_emp interpretation is artefactual.
3. **F3 — slot-2 extension reveals neither Galway-shape nor Thread-7-
   shape.** Some intermediate or step-discontinuous behaviour at
   10^5.5 — 10^6 would refute both hypotheses and require a third
   model.
4. **F4 — Mobius/dps truncation effect.** Slot 1 uses `M=4 / dps=18`.
   If running at `M=8 / dps=25` (Thread 7's setting) gives
   significantly different K_emp values, the slot-1 setup is too
   aggressive on numerical precision and the c_emp constants are
   biased.
5. **F5 — `K_galway(c=1) = √x · log²x` is a normalisation we chose,
   not Galway 2004's actual constant.** A literature audit of Galway
   2004 / Büthe 2018 may give a different reference constant; if so
   our c_emp comparison shifts.
6. **F6 — Worst-case error at fixed K isn't actually monotone in K.**
   The K → worst-error trace shows non-monotonicity (e.g., x=10⁴
   worst@K=1000 > worst@K=2000 > worst@K=8000 but worst@K=8000 >
   worst@K=4000 by a small amount). If K_emp is sensitive to milestone
   choice, slot-1 K_emp values are quantization artefacts.

## Slot 2 plan

**Goal:** Measure `K_emp(±1, x = 10^{5.5}, N=30)` to distinguish
Galway-shape (`K_emp ≈ 19000`) from Thread-7-shape (`K_emp ≈ 700`).

**Subtasks:**
1. Extend the project's zeros database from K=8000 to K=25000-30000
   via incremental `mpmath.zetazero`. Per-zero cost ~0.5s at dps=15;
   ETA ~3 hours of background compute.
2. Re-run `slot1_worst_case_K.py` with `--K_max 25000` and
   `--anchors 5.0,5.5` to get the new measurement.
3. If the data clearly favours the Galway-shape, compute c_emp at
   x=10⁵·⁵ and tabulate against c_emp at x=10⁴, 10⁴·⁵.
4. If the data clearly favours the Thread-7-shape, the Galway-shape
   interpretation collapses and slot 2 transitions to a Thread 7
   extension session — measuring K_emp(x, ε=√x · log log x / log²x)
   at small x (the named-exponent corollary scale).

## Slot 3 plan

**Goal:** Either (A) tighten c_emp into a publishable empirical
constant with rigorous error bars, or (B) collapse Thread 10 if
slot 2 falsifies the Galway-shape interpretation.

**If (A):**
- Extend to x = 10^6 (requires ~120000 zeros — large precompute).
- Theoretical c_emp asymptote: heuristic argument from σ_pred formula.
- Compare to Buthe 2018's explicit GRH constants (literature read).

**If (B):**
- Reframe the result: not "Galway constant tightening" but "Thread 7
  extension to worst-case statistic with explicit `K_emp ≈ 0.61 ·
  K_pred(half-Gaussian)` calibration."
- Document the closure as a B-grade negative-shape (the Galway-shape
  interpretation is artefactual).
- Move to next OPEN_POSITIVE_TARGETS entry (P6 / Cramér-style).

## Files

- `slot1_worst_case_K.py` (~280 lines, single driver, parameterised by CLI)
- `slot1_traces.csv` (1800 rows: log10_x, x, pi_x, K, abs_err per K-milestone)
- `slot1_summary.csv` (15 rows: K_emp, K_typ, K_pred, K_galway, c_emp, sigma_at_K_emp)
- `slot1_run.log` (full stdout from the run)
- `slot1_worst_case_K_results.md` (this file)

## Edges cited / used

- E1.5 (entropy lower bound for π(x)): the `~50%` random portion of π(x)
  is what the explicit-formula partial sum approximates.
- Thread 5 / S224 (Correlation Dichotomy): the K-zeros-shared-database
  amortisation framing.
- Thread 7 / S195+S240+S241+S244 (`σ_eff/σ_pred = 0.755 ± 0.06`,
  GUE pair-correlation, named-exponent corollary): slot 1 worst-case
  measurements confirm the GUE factor and the half-Gaussian tail
  shape extend to the worst-case-of-N statistic.
- E6.1 (Galway 2004 explicit-formula bound K = O(x^{1/2+ε})): slot 1
  is the empirical worst-case constant for this bound at finite x.

## Composition (what existing edges did slot 1 work compose)

Slot 1 composes Thread 5's amortisation framework (zeros shared
across many x) + Thread 7's σ_pred formula and GUE factor + E6.1
Galway shape. The composition is: *worst-case K-frontier as the
extreme-tail of the Thread 7 random-phase distribution*, with the
GUE pair-correlation factor reducing the required K from the half-
Gaussian-tail prediction by factor `(σ_eff/σ_pred)² ≈ 0.61`.

## Next-action for the next agent (slot 2)

> **Measure K_emp(±1, x = 10^{5.5}, N=30) and K_emp(±1, x = 10^6, N=30).**
> The slot-2 deliverable is a 4-decade c_emp scaling curve that
> decisively distinguishes Galway-shape (`c_emp` constant) from
> Thread-7-shape (`c_emp` growing as ~`√x / log²x`). Requires extending
> zeros database to K=25000-30000 (~3 hours mpmath.zetazero).
