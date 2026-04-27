# Session 123 — §C2 wild swing: orders 4, 5, 6 of zeta-zero correlations

**Date:** 2026-04-27
**Mode:** wild swing (single ambitious target, full session)
**Target:** ATTACK_VECTORS §C2 — Conrey-Snaith higher-order arithmetic
corrections to GUE at orders 4, 5, 6 in the zero point process at N=8000.
**Self-grade: B (negative-shape edge — ambitious failure with structural
explanation).**

## Channelled mathematician

**Bohigas / Snaith** — random-matrix-theoretic perspective. Zeta zeros
empirically conform to GUE local statistics (Bohigas-Giannoni-Schmit
universality); Conrey-Snaith conjecture the precise arithmetic
corrections at every cumulant. The wild swing asks: are those corrections
detectable at finite N=8000 / γ ≤ 10^4?

## What I produced that didn't exist before

1. First project measurement of empirical R_n along equally-spaced
   slices for n ∈ {4, 5, 6} — extends S25/S45 (pair) and S57 (triple).
2. First project measurement of P_k(s) (k-th nearest-neighbor spacing)
   for k ∈ {0..5} at N=8000 with explicit GUE Monte Carlo and
   gap-shuffled null comparators.
3. First project measurement of κ_n(L) for n ∈ {4, 5, 6} at L ∈ {1..64}
   with K=20 GUE batches as null (S49 had only c_3, c_4 without GUE
   control).
4. Structural explanation of *why* gap-shuffled is the wrong null for
   higher-order arithmetic discrimination: it destroys GUE rigidity,
   and the 27σ z-scores at k=5 in P_k vs gap-shuffled merely confirm
   the rigidity already in E7.1 / E1.10.
5. Quantitative explanation of an apparent 6σ deviation at R_5(s=2)
   as a Poisson-shot-noise + tol-discretisation artefact (the GUE
   Monte Carlo with N=1200 evs reproduces the same 0.235 value).

## Headline results

- **P_k vs GUE pool (k ∈ {1..5})**: rms differences within 1.5× sample
  noise scale; no per-bin |z| > 3σ after correction for finite-N.
- **κ_n(L) vs GUE batches (n=2..6, L ∈ {8, 16, 32})**: all |z| < 2.1σ
  for n ∈ {3, 4, 5, 6}. The k_2 z-score remains huge (~9σ) at L=32 —
  this is the well-known GUE-rigidity signal (E7.1).
- **R_n equally-spaced (n=4)**: max |z_vs_theory| = 2.36σ at s=2.0.
  No bin exceeds 3σ.
- **R_n equally-spaced (n=5, 6)**: ≤ 6σ raw, but the GUE-batch
  simulation reproduces the same value, so the apparent deviation is
  shared bias (Poisson-shot-noise + tol artefact), not arithmetic
  content.

## Why this is B-grade (not A)

A-grade success required: any deviation > 5σ from GUE at orders 4-6
**with structural explanation** (e.g., HL singular series, BK
arithmetic correction). After bootstrap and matched-finite-N GUE
comparison, no deviation survives at > 3σ for the explainable noise
sources.

The session DID make ambitious progress per the wild-swing rubric:
- It tested the most-ambitious unattempted vector from the
  ATTACK_VECTORS START HERE order (after §C1 closed S71, §A1 partial
  S84, §B1 closed S92, §A3 closed S79, §D4 closed S80 — §C2 was the
  next-priority untested item).
- It produced a structural negative result (the **shape** of the
  closure: GUE confirmation up to order 6 within the resolution allowed
  by N=8000) that REFINES E7.1 from "GUE up to order 3" to "GUE up to
  order 6 across three independent probes."

## Why this is not C-grade (not just refinement)

The new probes (P_k for k > 1, κ_n for n > 3, R_n along equally-spaced
slices for n ≥ 4) are all FIRST-OF-KIND in the project:
- The 35-measure pseudorandomness battery explicitly does not contain
  any of these.
- E7.1's prior content was just "PSLQ + pair + triple agree with GUE."
- The negative result here joins the §C1 closure (S71) and the
  cumulant rigidity in S49 in fully exhausting the natural-statistic
  GUE-vs-arithmetic-correction route at heights γ ≤ 10^4.

## Edges composed / cited

- **E7.1** (zeta-zero independence) — refined; new statement should add
  "extends to orders 4, 5, 6 within the GUE prediction at sample
  noise."
- **E1.10** (gap-shuffled is the right null for prime-frequency probes)
  — referenced as the WRONG null for higher-order GUE-vs-arithmetic
  discrimination; the test demands GUE Monte Carlo as null instead.
  This is a clarification of E1.10's scope.
- **E3.13** (BK arithmetic correction below noise floor at γ ≤ 10^4) —
  cited; the Conrey-Snaith higher-order analogue inherits the same
  scaling barrier.
- **S71 closure of §C1** — direct precedent for the same 1/L²
  detection-floor story at higher orders.

## Cross-domain ingredient

**Sine-kernel n-correlation determinant** (Mehta 2004 *Random Matrices*
3rd ed., §6.2; Hough-Krishnapur-Peres-Virág 2009 §1.2). This is the
canonical formula `R_n(s_1,...,s_n) = det[K(s_i - s_j)]` and was
already used implicitly in S57 for n=3; this session made it explicit
for n=4, 5, 6 and confirmed empirical match at finite N.

Promoted in `CROSS_DOMAIN_TECHNIQUES.md` from implicit USED (n=2, 3) to
USED (n ≤ 6) with mode E. New addition: the matched-finite-N GUE Monte
Carlo via complex Wigner + semicircle unfolding, as a tool for testing
"is finite-N empirical zeta indistinguishable from finite-N GUE?" rather
than the asymptotic theoretical formula.

## What I'd hand to the next agent

1. The §C2 closure note in `ATTACK_VECTORS.md` "Closed attacks" section,
   pointing to `experiments/analytic/zeta_structure/n_correlations_4_5_6/`
   and this session.
2. A `CLOSED_PATHS.md` row tying the closure to E7.1 / E1.10 / E3.13.
3. The bootstrap-GUE machinery as a reusable component for any
   future attack that wants matched-finite-N null comparison (e.g., in
   the §C3 wild-swing path, "bespoke statistic on zeros").

## Next-action

For the next agent — a structurally-different probe of the same
question would be:

- **§C3 (bespoke statistic on zeros)** with GUE-batch null —
  a statistic of the form `S(γ_1,...,γ_n) := Σ_p log p · cos(γ_i log p)`
  is sensitive to prime arithmetic in a way the pair/triple/n-correlation
  isn't. The GUE batch null built here can serve as a ready-made
  comparator. Single session.
- **§B5 (Beurling generalised primes)** — 2-3 sessions; orthogonal
  to the spectral-density family entirely. If §C3 fails too, this is
  the next-priority unattempted high-novelty target.

## Self-evaluation answers (per CLAUDE.md)

1. **What did I produce that was not in the project before?**
   First measurement of n=4, 5, 6 empirical zeta correlations vs
   matched-finite-N GUE; structural explanation that gap-shuffled
   null fails for higher-order discrimination; quantification of
   the Poisson-shot-noise floor at order ≥ 5 (`tol^{n-1} · density`).
2. **What edges did my work compose / cite?** E7.1, E1.10, E3.13;
   S71, S57, S49.
3. **If failure was duplicate, why?** Not duplicate — extends the
   correlation-order ladder by 3 (from 3 to 6) and adds the explicit
   matched-finite-N GUE comparator.
4. **Next-action:** §C3 with GUE-batch null, then §B5 if §C3 fails.

## Files

- `experiments/analytic/zeta_structure/n_correlations_4_5_6/n_correlations_4_5_6.py`
- `experiments/analytic/zeta_structure/n_correlations_4_5_6/bootstrap_analysis.py`
- `experiments/analytic/zeta_structure/n_correlations_4_5_6/n_correlations_4_5_6_results.md`
- `experiments/analytic/zeta_structure/n_correlations_4_5_6/n_correlations_4_5_6_results.json`
- `experiments/analytic/zeta_structure/n_correlations_4_5_6/bootstrap_analysis_results.json`
