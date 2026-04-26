# Session 74 — C2: Free cumulants of the chi_P MPS unfolding operator

**Date:** 2026-04-26
**Mode:** construction
**Target:** NOVELTY_CHALLENGES.md §1 C2 — MPS bond-dim × free-probability moments
**Composes:** E2.1 (TT/MPS bond-dim identity) × free probability framework
**Cross-domain:** Voiculescu free probability (Mingo & Speicher 2017,
Ch. 1–2, Ch. 4)

## Self-grade: B

Substantive refinement of E2.1 with a cross-domain technique that the
project had not used at the operator level. The construction produced
three new artifacts (Mertens-product = MP-rate identification; spike-count
regularity `k* ∝ R^{0.85}`; verified MP(c) standardized cumulant identity
`κ_r = c^{r-1}`), all falsifiable, all empirically verified across
W ∈ {2, 3, 5, 6, 30}. No polylog opening — the spike-count regularity
recovers a polynomial-in-N spectral compression barrier parallel to
Lagarias-Odlyzko, which is why it lands at B and not A.

## What I produced that was not in the project before this session

1. **Free-cumulant evaluator** for the chi_P MPS unfolding operator
   (`free_cumulants_chi_p.py`). Computes moments + free cumulants under
   several normalizations (full spectrum / nonzero-only / drop-k bulk),
   plus matched iid Bernoulli baselines (full-shape and active-block).

2. **MP-bulk identification.** Empirically, the bulk free cumulants of
   the chi_P MPS spectrum match Marchenko-Pastur with rate
   `c = φ(W)/W = ∏_{p ≤ W}(1 − 1/p)`. The Mertens product is exactly
   the free-Poisson rate of the bulk. This is a **non-trivial
   identification** of a free-probabilistic invariant with a sieve
   density — and as far as I can find, it does not appear in the
   project's prior work or in `pseudorandomness_of_pi.md`'s
   35-measure battery (which is about pi(x) mod 2 specifically, not
   about the operator spectrum of chi_P).

3. **Spike-count regularity.** Define `k*(W, d)` = smallest `k` such
   that the drop-k bulk free cumulants of chi_P match MP(c=φ(W)/W)
   within 10% relative on `(κ_2, κ_3, κ_4)`. Empirically,
   `k*(W, d) ≈ R^{0.85}` (fitted on W=2, d=14..22), where `R` is the
   rank from E2.1. The matched-active-Bernoulli baseline has `k* ≤ 1`
   for every (W, d) tested. So the spike-band excess is a deterministic
   feature of chi_P that grows polynomially in N — a structural fact
   beyond E2.1's rank statement.

4. **Verification of the MP(c) standardized cumulant identity
   `κ_r = c^{r-1}`** on a Gaussian iid (4000 × 1000) ensemble. (The
   κ_r convention I had initially written was wrong by an exponent
   sign; the empirical sanity check caught it.)

5. **Annotation of E2.1** (EDGES.md) recording the moment-level
   refinement.

## What edges did my work compose / cite?

- **E2.1** (composed) — TT/MPS bond-dim identity. The construction takes
  the same unfolding `M^(j)` and probes its singular value distribution
  beyond the rank.
- **Free probability framework** (cross-domain). Reference: Mingo &
  Speicher (2017). The S53 closure on `δ(x)/√x` was a *scalar* free-
  moment test; this is the *operator* free-cumulant test, a different
  object.
- Cross-references / citations: E2.1 (refined), Mertens product (used
  by name), Marchenko-Pastur (used by name), BBP transition / deformed
  Wishart (referenced as related random-matrix-theory regimes).

## Honest failure / partial-result reporting

The construction is **NOT** A-grade because:

- The spike-count regularity `k* ∝ R^{0.85}` reveals a polynomial-in-N
  spectral compression barrier, not a polylog one. The barrier is
  parallel to Lagarias-Odlyzko `O(N^{2/3})` and the explicit-formula
  `O(N^{1/2+ε})` bound — not an improvement.
- The Mertens-product identification, while clean and (I believe) new
  to this project, reduces algorithmically to E2.1 + a known fact
  about MP(c). It refines but does not transcend E2.1.
- The cross-domain import succeeded (free probability gave a clean
  bulk identification) but did NOT produce a polylog handle that the
  project previously lacked. Free-Poisson distributions admit
  R-transform inversion for sampling, but inversion of `M^(j)` doesn't
  buy us anything because we already KNOW chi_P (it's a sieve).

The construction is **NOT** F-grade because:

- The MP-bulk = Mertens-product identification is a new, falsifiable,
  empirically-verified identity tying free probability to the
  wheel-W sieve.
- The spike-count regularity is a new empirical fact with a clean
  fitted exponent and reproduces across 5 values of W.
- The verified κ_r = c^{r-1} identity is itself a small piece of
  rigor that future agents working with MP-style baselines will need.

So **B-grade** is the honest call.

## Next actions for the next agent

Three concrete next-actions, in declining order of value:

1. **Characterize the spike eigenvectors.** The spike band of
   `O(N^{0.42})` outliers must correspond to specific structured
   directions in `(C^W)^{⊗d}`. Candidates: small-prime indicator
   eigenfunctions, Selberg eigenfunctions, cross-block prime-density-
   variation modes, or twin-prime / prime-constellation correlations.
   The first few spike eigenvectors are computable by SVD truncation
   (just `M^(j) v_i = σ_i u_i`); compare them to known structured
   vectors. If the spike eigenvectors are exactly small-prime
   indicator functions, this would be a clean structural theorem.
   Effort: 1 session.

2. **Push the W=2 sweep to d=24, 26.** Confirm or refute the
   `k* ∝ R^{0.85}` exponent. The exponent could converge to 0.5
   (clean √rank), 1.0 (rank-deficit growth), or settle at 0.85.
   Effort: 0.5 sessions (SVD on 4096 × 4096 takes ~30s).

3. **Try other tensor reshapes.** The base-W reshape gives MP-bulk +
   spikes. What about CP / Tucker / hierarchical Tucker / MERA
   reshapes? N1 (negative-shape conjecture in NOVELTY_CHALLENGES.md
   §4) hypothesises a unified "no polylog tensor compression"
   barrier. The free-cumulant evaluator built here is reusable for
   testing whether other tensor structures of chi_P also produce an
   MP-bulk-plus-spike spectrum.
   Effort: 1-2 sessions for code reuse + 4-5 reshapes.

## RESEARCH_AGENDA.md update

Arc 4 advanced: C2 marked BUILT. Three new sub-action candidates added
(see "Next actions" above). C3, C4, C6 remain open and substantive.

## Session-end self-evaluation

1. **Q: What did I produce that was not in the project before this
   session?** A: Free-cumulant evaluator + Mertens-product=MP-rate
   identification + spike-count regularity + verified κ_r = c^{r-1}.
2. **Q: What edges did my work compose or cite?** A: E2.1 (composed
   and refined inline); cross-domain Mingo-Speicher 2017.
3. **Q: If only duplicate closures, why?** A: Not applicable — this
   session produced new content.
4. **Q: Next-action for the next agent?** A: characterize the spike
   eigenvectors of `M^(j)` (top 1 sub-action above), or push the W=2
   sweep to d=24, 26 (top 2 sub-action above).

## Files

- `experiments/constructions/free_cumulants_chi_p/free_cumulants_chi_p.py` — evaluator
- `experiments/constructions/free_cumulants_chi_p/free_cumulants_chi_p_results.md` — results
- `experiments/constructions/free_cumulants_chi_p/free_cumulants_chi_p_results.json` — JSON output
- `experiments/constructions/free_cumulants_chi_p/definition.md` — formal construction spec
- `experiments/constructions/free_cumulants_chi_p/run_full.log` — console transcript
- `EDGES.md` — E2.1 annotated with S74 free-probability refinement
- `NOVELTY_CHALLENGES.md` — C2 marked BUILT
- `RESEARCH_AGENDA.md` — Arc 4 advanced (C2 milestone done)
- `status/CLOSED_PATHS.md` — row 748 added
