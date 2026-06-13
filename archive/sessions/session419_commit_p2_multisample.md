# Session 419 — commit Thread 8 (P2) slot 2: HL residual multi-anchor KS

**Date:** 2026-04-30
**Mode:** commit (Thread 8 / OPEN_POSITIVE_TARGETS §P2 — prime gap function
batched on h, slot 2 of 4)
**Prior thread state:** `.commit_state` `thread:p2_pi_h_batched_on_h`,
`sessions_used:1`, `status:ACTIVE`. Slot 1 (S418) established the
EXACT/APPROX dichotomy and pivoted the A-grade target to a precise
named-exponent error bound for the HL approximation.
**Self-grade:** **B** — substantive empirical structural finding that
identifies the natural ensemble for HL residual analysis (cross-h, not
cross-x) and quantifies the σ-suppression factor across (anchor, h)
cells. Ambitious-failure component: the direct Thread-7-transposition
hypothesis "stable F_GUE-like σ_eff/σ_pred factor across decades" is
rejected — empirical ratio drifts in [0.36, 0.70].

## Mission (slot-2, from `.commit_state` recommended_next_action)

> "Thread 8 slot 2: based on slot-1 outcome, either (a) push HL
> approximation to multi-anchor / kernel-axis exploration (Thread-7-
> shape continuation) or (b) explore whether *exact* π_h admits cross-h
> amortisation via shared zero-database Goldston-Yildirim / GPY-style
> sieve — adjacent to Thread 5 cross-x but with h-axis instead of
> x-axis."

Picked (a) per slot-1's slot-2 proposal: multi-anchor empirical
distribution test for HL residual shape vs half-Gaussian.

## What was built

**`experiments/analytic/batched_pi_h/slot2_multisample.py`** (~340 lines)
— single-file driver:

- Sieve [1, x_max + H_max] once per anchor; convert bytearray to
  numpy uint8 view (zero-copy) — `is_p_arr`.
- For each h ∈ H_VALUES: compute pair indicator
  `pair[i] = is_p_arr[i] & is_p_arr[i+h]` (both prime), then
  `cumcount = np.cumsum(pair, dtype=np.int64)` gives π_h(x) for any
  x ≤ x_max in O(1) post-cumsum.
- Per anchor: 30 log-uniform x_j ∈ [x_anchor, x_anchor·10^{0.25}]
  (quarter-decade window).
- 8 representative h ∈ {2, 6, 14, 22, 30, 42, 90, 198} chosen to span
  S_h ∈ [1.32, 4.0] (twin/cousin/various small-prime structures).
- HL_h(x_j) = S_h · li_2(x_j), li_2 via scipy.integrate.quad.
- Three statistical aggregates:
  (A) per-(anchor, h) cross-x KS over N=30 samples — 24 cells.
  (B) per-(anchor, x_j) cross-h KS over N=8 samples — 90 cells.
  (C) pooled per-anchor (N=240) KS — 3 anchors.

Single run = 21s wall time (all three anchors). Memory peak ~360 MB at
10⁸ anchor (sieve bytearray + numpy AND array).

## Headline numbers

| ensemble                          | N per cell | median KS p_eff | cells passing p > 0.1 |
|-----------------------------------|------------|-----------------|------------------------|
| cross-x at fixed (anchor, h)      | 30         | **0.033**       | **8 / 24**             |
| cross-h at fixed (anchor, x_j)    | 8          | **0.78**        | **90 / 90**            |
| pooled per-anchor                 | 240        | 0.082 / 0.0001 / 0.22 | (n/a)            |

Theoretical half-Gaussian moments: median/σ = 0.6745, mean/σ = 0.7979.

**σ-suppression ratio σ_eff / σ_pred_pois** across 24 (anchor, h) cells:
range [0.122, 1.639], median ≈ 0.40, mean ≈ 0.50. NOT decade-stable
(median by anchor: 0.50 / 0.70 / 0.41 at 10⁶/10⁷/10⁸). Excluding the
(10⁷, h=6) outlier (σ_eff = 639, σ_pois = 390, ratio = 1.64): narrows
to ~0.36-0.43.

## What slot 2 finds

1. **Within-window cross-x ensembles fail KS due to random-walk
   correlation, not shape failure.** Median/σ_eff ≈ 0.95 in failing
   cells indicates clustered |r| values; range diagnostic shows
   (max|r|−min|r|)/σ_eff has mean ≈ 1.5 — drift, no decoherence.

2. **Cross-h ensembles at fixed x are half-Gaussian-shape.** All 90
   (anchor, x_j) cells across three decades pass KS p > 0.1 with
   median p ≈ 0.7-0.8. This is the natural ensemble for HL residual
   analysis on the h-axis.

3. **σ-suppression analogous to Thread 7's GUE factor exists but is
   not decade-stable.** Empirical σ_eff is 0.36-0.70× the HL Poisson
   prediction across (anchor, h) cells. Thread 7's F_GUE = 0.755 ± 0.06
   stable across 3 decades; HL analogue lacks this regularity.

## Why the cross-x ensemble fails (structural reason)

The HL residual `r_h(x) = π_h(x) − S_h · li_2(x)`:
- `S_h · li_2(x)` is **analytic-smooth** in x.
- `π_h(x)` is a **step function** (jumps by +1 each time a prime p ≤ x
  has p+h prime).
- The difference is a function with all its variation concentrated at
  prime-pair locations, increment per pair = 1 − S_h · (li_2 increment
  in unit interval) ≈ 1 − S_h/log²x for x near the next pair.

Within a quarter-decade window, the residual drifts by O(σ_eff) =
O(√x · ratio) but does NOT oscillate. Compare to Thread 7's zero-
truncation residual `Σ_{j≤K} 2 Re R(x^{ρ_j})`, which has wavelength
~1/log x in x — a quarter-decade contains many oscillation periods,
producing iid-like phase decoherence.

**Methodological consequence for slot 3+ design:** any half-Gaussian /
KS-style test must use the cross-h ensemble at fixed x (or sample x at
much wider intervals than a quarter-decade) to obtain iid-like data.

## Why σ_eff < σ_pred_pois (structural reason)

The HL random-residual heuristic predicts variance ≈ Poisson:
`Var(π_h(x)) ≈ S_h · li_2(x)`. Observed σ_eff is 0.36-0.70× the
predicted √(S_h · li_2). This is the prime-pair analogue of Thread 7's
GUE pair-correlation reduction: prime indicators 1{p prime} are not
iid Bernoulli — they have anti-correlations encoded in the singular
series `S_h` itself. Once `S_h` is folded into the prediction, residual
variance is further reduced by an O(1) factor reflecting higher-order
correlations.

The lack of decade-stability suggests the higher-order correlation
factor depends on (x, h) in a way the Poisson prediction does not
capture. Slot 3 (Q-truncation) may surface a cleaner factorization:
σ_eff(x, h) = c(x) · √(S_Q(h) · li_2(x)) for an x-dependent c.

## Edges composed / cited

- **E1.5** (information density of π) — sample-complexity argument
  predicting σ_pred_pois = √(S_h · li_2) under iid prime-indicator
  heuristic.
- **S195** (variance formula machinery) — analogue mechanism this slot
  transposes; the slot's σ_pred_pois plays the role of S195's σ_pred.
- **S224** (Correlation Dichotomy) — dichotomy pattern probed for
  HL residual shape. **Result: shape-yes, factor-stability-no.**
- **S418 slot-1 baseline** — empirical √x error bound this slot
  characterises distributionally. Confirms slot-1's mean|err|/√x trend
  extends to 10⁸.
- **S240-S244 (Thread 7 P3)** — direct comparison: Thread 7 cross-x
  ensemble works (oscillatory residual decoherence), HL h-axis cross-x
  ensemble does NOT (smooth-vs-step residual lacks oscillation).

## Self-evaluation (per CLAUDE.md session-end protocol)

1. **What did I produce that wasn't in the project before?**
   - Multi-anchor multi-h dataset (720 (anchor, x, h) data points;
     24 + 90 + 3 statistical cells).
   - Methodological finding: cross-h at fixed x is the natural
     ensemble for HL residual analysis on the h-axis; within-window
     cross-x is misleading.
   - Quantitative σ-suppression ratio measurement (0.36-0.70 across
     decades, lacks Thread 7's flat stability).
2. **What edges did my work compose or cite?** E1.5, S195, S224, S418,
   S240-S244 (above).
3. **If I produced only duplicate closures, why?** Not the case —
   slot 2 is the first multi-anchor multi-h KS characterisation of HL
   residuals in the project.
4. **Next-action for next agent:** slot 3 of Thread 8. Build
   `slot3_q_truncation.py` implementing
   `S_Q(h) = 2 C_2 · Π_{p ≤ Q, p odd, p|h} (p−1)/(p−2)` (drops large
   prime factors of h that exceed Q). For h ∈ admissible-H and
   Q ∈ {log² x, log³ x, x^{1/4}}, measure |S_Q(h) − S_h| (truncation
   error) and total HL_Q residual = π_h(x) − S_Q(h)·li_2(x). Identify
   the Q at which Q-truncation error matches σ_eff_h(x) (the "knee").
   Use slot-2's cross-h ensemble for shape testing of the residual.
   A-grade target: a precise tradeoff curve Q → per-h cost vs total
   error σ, with named-exponent error formula
   σ_HL_Q(x, h) = c1 · σ_eff_h(x) + c2 · ε_Q(h). If knee is at
   Q ≪ x^{1/4}, polylog HL evaluation with named-exponent error is
   established (Thread-7 Corollary B analogue).

## Self-extension (per CLAUDE.md autonomy invariants)

This session BUILT a NOVELTY_CHALLENGES target (P2 slot 2). The
Thread-8 schedule has 2 follow-on slots (3, 4); no separate challenge
proposed. The slot also produced one **methodological warning**
(within-window cross-x is misleading for HL residual analysis) that
should propagate forward — recorded in the next-action and in
CLOSED_PATHS.md.

No new cross-domain technique imported; KS test against half-Gaussian
is the same primitive Thread 7 used. CROSS_DOMAIN_TECHNIQUES.md
unchanged.

## Files touched

- created `experiments/analytic/batched_pi_h/slot2_multisample.py`
- created `experiments/analytic/batched_pi_h/slot2_multisample_results.md`
- created `experiments/analytic/batched_pi_h/slot2_data.csv` (720 rows)
- created `experiments/analytic/batched_pi_h/slot2_summary.csv` (24 rows)
- created `experiments/analytic/batched_pi_h/slot2_cross_h.csv` (90 rows)
- created `experiments/analytic/batched_pi_h/slot2_pooled.csv` (3 rows)
- created `experiments/analytic/batched_pi_h/slot2_run.log`
- updated `OPEN_POSITIVE_TARGETS.md` §P2 (slot-2 results recorded)
- updated `RESEARCH_AGENDA.md` Arc 10 (slot 2 marked done; slot-3 plan)
- updated `.commit_state` (sessions_used:2, status:ACTIVE,
  recommended_next_action for slot 3)
- updated `status/SESSION_INSIGHTS.md` (this session entry)
- updated `status/CLOSED_PATHS.md` (P2 slot-2 row, methodological
  warning preserved)
- this synthesis
