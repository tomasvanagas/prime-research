# Session 242 — commit Thread 7 slot 3: smoothed kernels do not beat hard cutoff for π_K(x)

**Date:** 2026-04-30
**Mode:** commit (Thread 7 / OPEN_POSITIVE_TARGETS.md §P3 polylog
approximate π(x) with named ε)
**Slot:** 3 of 5 (continuing fresh thread from S240, S241)
**Self-grade:** **B** — substantive empirical refinement closing the
S202-wrap §"Non-Gaussian smoothing kernels" legitimate-falsifier
listing. Negative-shape result. See full detail in
`experiments/analytic/polylog_approx_pi/smoothed_kernels_results.md`.

## Mission

S202 wrap explicitly listed compactly-supported smoothing kernels
(sinc, sech², Hann, raised-cosine, etc.) as a *legitimate falsifier*
of the Thread 3 closure that S196 had only addressed for the
log-Gaussian (multiplicative damping) case. The structural prediction
in S202: "Var_TAIL is determined by truncation, not kernel" — i.e.,
no compactly-supported kernel can beat hard cutoff. But this had
*not* been empirically tested.

Slot 3 task per `.commit_state`: *smoothing kernel selection — do
Gaussian / raised-cosine / Hann / Hamming / Triangle / Riesz / Tukey
kernels reduce ε beyond the unsmoothed sum at fixed K_compute?*

## What slot 3 produced

1. **Smoothed-kernel evaluator (`smoothed_kernels.py`, ~340 lines)**:
   per-sample cumulative loop computing c_k = 2 Re R(x^{ρ_k}) for
   k = 1..K_max = 4000 (one R_at_rho call per (sample, k)), then
   numpy-dot post-processing for each (kernel, K_compute) combination
   sharing the same R_at_rho work. 9 kernels, 4 K-values, 20 paired
   samples, 3 anchors → 108 (anchor, K, kernel) summary cells from
   2160 (anchor, x, K, kernel) per-sample data triples.

2. **Headline closure**: ZERO of 96 (anchor, K, kernel != hard) cells
   show σ_eff(kernel) < σ_eff(hard) at paired sign-test p < 0.05
   (one-sided). Mean σ_eff/σ_eff(hard) by kernel across 12 (anchor, K)
   cells:
   - tukey25 = **1.04**, tukey50 = 1.11, cosine = 1.14, riesz = 1.12,
   - triangle = **1.23**, riesz4 = 1.21, hamming = 1.20, hann = **1.23**.
   All > 1. Decisive hard wins (p < 0.01) at 10⁸/K = 2000+4000 and
   10⁹/K = 2000+4000.

3. **Asymmetric-pattern verification**: 16 of 24 ratios ≥ 1.05 are
   1-σ above 1.0; 0 of 9 sub-1 ratios are 1-σ below 1.0. The smallest
   "kernel beats hard" p-value is 0.13 at 10⁷/K=1000/tukey25 — not
   significant at 5%.

4. **Paired-sign-test post-processor (`smoothed_kernels_paired_analysis.py`)**:
   for each (anchor, K, kernel), wins/N = #{i : |err_kernel(x_i)| <
   |err_hard(x_i)|}, with one-sided p-values for "kernel beats hard"
   and "hard beats kernel" under Binom(N, 0.5) null. Surfaces the
   asymmetric pattern explicitly.

5. **Structural-reason confirmation**: under random-phase
   Var(error) = Σ_j (1−w_j)² Var_j is minimised at hard cutoff
   (w_j = 1 ∀ j ≤ K, 0 else) — the L2-optimality argument is exact.
   The 0.74 GUE factor scales σ_eff *uniformly* across kernels:
   σ_eff(hard)/σ_pred clusters at 0.62-0.85 across 12 cells; GUE
   corrections do NOT discriminate between kernels at the second-
   moment level. **The S202 prediction holds.**

## What slot 3 did NOT produce

- An A-grade. The closure is *negative-shape* — no kernel beats hard,
  consistent with the theoretical L2-optimality. CLAUDE.md grading:
  "negative-shape closure of a legitimate falsifier" is B-grade.
- A non-symmetric / position-correlated kernel test. Slot 3 covered
  symmetric kernels w_k = w(u) where u = (k-1)/(K-1). Paired-weight
  kernels (w_{2j-1} = w_{2j}) exploiting Wigner repulsion remain
  open and are the slot-4 candidate.
- An empirical anchor at x = 10¹². Deferred from slot 2 with same
  reasoning — K_max=8000 caps the policy range and another decade
  would not change the headline.
- A rigorous variance bound. Slot 5 territory.

## Edges composed / cited

- **E1.5** (information-theoretic per-query barrier) — the
  L2-optimality argument shows that within the linear partial-sum
  framework, no compactly-supported kernel can extract more
  information per zero than hard cutoff.
- **E3.1** (Connes-Consani-Moscovici spectral triple) — Thread 3
  closure-of-closure transitivity; this slot adds the kernel-axis
  closure complementing S196's bandwidth-axis closure.
- **S195 σ-formula** — the prediction validated empirically; its
  L2-optimality is the structural reason behind slot 3's outcome.
- **S196** (log-Gaussian smoothing closure) — extended to the
  compactly-supported case here.
- **S202 unified Connes-Galway theorem** — §"Non-Gaussian smoothing
  kernels" legitimate-falsifier closed.
- **S240** (Thread 7 slot 1 named-exponent corollary) — confirmed
  kernel-optimal up to GUE-uniform 0.74 scale factor.
- **S241** (Thread 7 slot 2 distribution shape test) — half-Gaussian
  shape preserved under all tested kernels (KS p-values uniformly
  > 0.1, median ~0.7 across the kernel family).

## Cross-domain ingredient

Window-function theory (signal processing): Bartlett, Hann, Hamming,
Tukey, Riesz, raised-cosine windows are standard in spectral
estimation for trade-offs between mainlobe width and sidelobe
suppression. Their use here is a direct import: each window is
applied as a multiplicative weight on the explicit-formula partial
sum, and the L2-optimality of hard cutoff in this context is a
consequence of the orthogonal-decomposition variance formula under
random-phase.

`CROSS_DOMAIN_TECHNIQUES.md` should mark "compactly-supported window
functions for explicit-formula partial-sum smoothing — USED-E (S242)
— confirmed no improvement over hard cutoff for the partial-sum
approximation π_K(x)".

## Falsifiability

Slot 3's claim is falsified by:

1. A compactly-supported kernel w on [1, K_compute] yielding
   σ_eff(kernel) / σ_eff(hard) < 0.95 with paired sign-test p < 0.01
   at any tested (anchor, K) cell. *Directly tested across 96 cells;
   ruled out.*
2. A non-symmetric / position-correlated kernel structure (paired
   or position-dependent weights) yielding σ_eff < σ_eff(hard) at
   meaningful significance. *Slot 4 candidate; not ruled out.*
3. A second-moment GUE correction calculation showing that the
   pair-correlation term σ²_{GUE-cov} is *kernel-dependent* in such
   a way that some kernel has it negative. *Theoretical, outside
   current machinery.*

## What was built

1. `experiments/analytic/polylog_approx_pi/smoothed_kernels.py`
   (~340 lines): smoothed-kernel evaluator with 9 kernels, paired
   evaluation across 20 samples × 4 K_targets × 9 kernels, summary
   CSV output, headline ratio table, predicted-ratio (uniform-variance)
   sanity column.
2. `experiments/analytic/polylog_approx_pi/smoothed_kernels_paired_analysis.py`
   (~80 lines): paired-sign-test post-processor on the per-sample CSV.
3. `experiments/analytic/polylog_approx_pi/smoothed_kernels_data.csv`:
   2160 (anchor, x, K, kernel) per-sample rows.
4. `experiments/analytic/polylog_approx_pi/smoothed_kernels_summary.csv`:
   108 (anchor, K, kernel) summary rows with σ_pred, σ_eff,
   ratio_obs/hard, KS p-value.
5. `experiments/analytic/polylog_approx_pi/smoothed_kernels_main.log`:
   full run log.
6. `experiments/analytic/polylog_approx_pi/smoothed_kernels_paired.log`:
   full paired-sign-test grid.
7. `experiments/analytic/polylog_approx_pi/smoothed_kernels_results.md`:
   slot-3 write-up with full headline grid, per-kernel summary,
   structural reason, falsifiability, self-grade.
8. `status/CLOSED_PATHS.md`: §P.P3 slot-3 row appended.
9. `OPEN_POSITIVE_TARGETS.md` §P3: slot-3 negative-shape closure
   and slot-4 next-action options added.
10. `status/SESSION_INSIGHTS.md`: Session 242 entry appended.
11. `.commit_state`: sessions_used 2 → 3, session_history S240,S241
    → S240,S241,S242, last_synthesis updated, slot-4 next-action
    detailed.
12. `archive/sessions/session242_commit_p3_smoothed_kernels.md`:
    this file.

## Self-evaluation (per CLAUDE.md 4-question)

1. **What did I produce that was not in the project before this session?**
   - Empirical kernel-axis closure of the Thread 7 P3 partial-positive
     bound. S202 explicitly listed compactly-supported smoothing
     kernels as untested; slot 3 tests them with 2160 paired data
     triples and confirms hard cutoff is L2-optimal across the
     compactly-supported family.
   - Paired-sign-test methodology for kernel comparison: a binary
     pairs-test that controls for sample-level variance and
     extracts kernel signal cleanly even at N=20.
   - The asymmetric pattern (16 of 24 ratios ≥ 1.05 are 1-σ above
     1.0; 0 of 9 sub-1 ratios are 1-σ below 1.0) confirms the GUE
     0.74 factor scales σ_eff uniformly across kernels.
   - 9 kernel weight functions implemented as numpy primitives
     (Hann, Hamming, Tukey-25, Tukey-50, Triangle, Riesz, Riesz⁴,
     Cosine) for use in any future slot 4-5 experiments.

2. **What edges did my work compose or cite?**
   - E1.5, E3.1, S195, S196, S202, S240, S241.

3. **If my session produced only duplicate closures, why?**
   - Did not. The compactly-supported smoothing test was an explicitly-
     listed open falsifier in S202, and slot 3 resolves it with
     primary empirical data (96 paired comparisons). The empirical
     test is genuinely new content; the structural prediction was
     known but unverified.

4. **What is the next-action for the next agent?**
   - **Slot 4 of Thread 7**: option (a) non-symmetric / position-
     correlated kernels (paired weights w_{2j-1} = w_{2j} exploiting
     Wigner repulsion in adjacent zero pairs). This is the ONLY
     remaining kernel-family direction not yet tested. If σ_eff
     paired < σ_eff hard with paired sign-test p < 0.05, that's
     an A-grade-shaped result (GUE pair correlation IS computationally
     exploitable for the partial-sum approximation). If it fails,
     pivot to option (c) slot-5 theoretical wrap. See `.commit_state`
     `recommended_next_action`.

## Honest summary

Slot 3 closes the kernel-axis falsifier of the Thread 3 closure that
S196 had only addressed for log-Gaussian smoothing (multiplicative
damping). The compactly-supported family — Hann, Hamming, Triangle,
Riesz, Tukey, Cosine — was the open piece of the kernel-axis story
per S202 wrap. Slot 3 tests 9 such kernels at matched K_compute and
finds none beat hard cutoff in any (anchor, K) cell at meaningful
significance. The asymmetric pattern (no significant sub-1 ratio,
many significant > 1 ratios) directly confirms the L2-optimality
argument predicted in S202.

The closure is negative-shape — no positive algorithmic improvement
over the slot-1 named-exponent bound. The slot's contribution is to
document that the only kernel-shape lever within the symmetric
compactly-supported family does not move σ_eff. To reach A-grade for
Thread 7 there are now three remaining levers: non-symmetric kernels
(slot 4 option a), non-linear post-processing (outside the kernel
framework), or a rigorous variance proof (slot 5 territory).

The session produces no new theorem, no algorithm faster than slot 1's
heuristic O(polylog(x)) named-ε partial-sum evaluator. But the kernel-
axis is now closed structurally + empirically, which is the right
methodology to keep Thread 7 honest. **B-grade negative-shape
closure**, not A.
