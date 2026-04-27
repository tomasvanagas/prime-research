# Session 127 — C8 (depth-2 sign-threshold W-vs-M tradeoff for PRIMES)

**Mode:** construction
**Composition:** E5.3 × S84 framework × E1.6 / C7-S89
**Save:** `experiments/constructions/d2_sign_threshold_w_m_tradeoff/`
**Self-grade:** **B**

## What this session produced that did not exist before

1. **The first measured weight-vs-size tradeoff curve for PRIMES** (or
   any natural number-theoretic Boolean function) in the depth-2
   sign-threshold class. S84 had measured the W=1-only column at N=6
   and at N=8; S127 extends to W ∈ {2, 3, 4, 8, 16, 32, 64} at N=6 and
   W ∈ {2, 4, 8} at N=4.

2. **A quantitative structural fact about PRIMES at N=6:** the
   depth-2 sign-threshold gate count `M*` has a **structural floor of
   3 reached at modest weight W=3 and held through W=64** (a 16×
   weight increase yields zero gate reduction). This is not a
   geometric-decay regime — the curve is "step-down then plateau,"
   sharply.

3. **A cross-N analog at N=4** (PRIMES `M* ∈ {3, 2, 2, 2}` vs random
   `M* ∈ {4, 3, 3, 3}`) showing the PRIMES-vs-random gap from C7-S89
   persists across weight settings, with the predicted Δ=1 gate gap
   at every W ∈ {1, 2, 4, 8}.

4. **A non-trivial M=3 W=4 witness circuit for PRIMES at N=6**
   verified 64/64 — the three half-spaces are *not* canonical
   primality predicates (no residue isolation, no symmetric
   arithmetic feature); the ILP found an opaque cover.

5. **An asymmetry in ILP search difficulty:** even at 5× the time
   budget (600 s vs 113 s), CBC fails to resolve `(W=4, M=3)` for
   matched-density random functions where it resolves PRIMES SAT in
   113 s. This is itself a structural fact about the search-problem
   complexity, not necessarily about the underlying circuit
   complexity.

## Edges composed / cited

- **E5.3** (PRIMES TC⁰ open frontier) — refined with a `(W, M*)`
  point at N=6.
- **S84 framework** (`sat_depth2_ilp.py`) — extended W=1 to a grid.
- **E1.6** (oddness predictor) — provides the structural reason the
  ILP exploits at PRIMES that it cannot exploit at random.
- **C7-S89** PRIMES-vs-random oddness gap — extended from a single-W
  observation to a curve.

Annotation in `EDGES.md` E5.3 and `NOVELTY_CHALLENGES.md` C8.

## What kind of session this was

This is a **B-grade construction**: a quantitative measurement that
adds a new data column on an existing edge (E5.3) and a new structural
finding (the M=3 plateau across 16× weight increase). It is *not*
A-grade because:

- No cross-domain technique imported (`CROSS_DOMAIN_TECHNIQUES.md`
  unchanged); purely circuit-complexity / ILP work.
- The structural floor at M=3 does not separate any complexity class.
  PRIMES at N=6 ∈ TC⁰₂[W=3, M=3] is a fixed-N statement; no asymptotic
  conclusion follows.
- The N→∞ rate of `M*(N; W)` is not measured. The asymptotic
  question — does `M*(N; W=O(1))` grow or stay constant? — is the
  one that would matter for E5.3.

It is **not** C-grade because:

- It does close a B-grade NOVELTY_CHALLENGES target (C8, the only
  remaining open §1 composition challenge).
- It produces 5 net new artifacts (curve, plateau, cross-N analog,
  witness, time-asymmetry).
- It pre-states 5 falsifiers and reports honestly which fail.

## Self-evaluation per CLAUDE.md

1. **What did I produce that was not in the project before this session?**
   The N=6 W-vs-M curve for PRIMES depth-2 sign-threshold with
   structural-floor identification at M*=3 across W ∈ {3..64}; the
   N=4 cross-N analog confirming PRIMES-vs-random gap persists across
   weight; a 64/64-verified witness circuit at (N=6, W=4, M=3); a
   time-asymmetry observation between PRIMES and random ILP search.

2. **What edges did my work compose or cite?**
   E5.3 (refined inline); S84 framework (extended); E1.6 / C7-S89
   (extended scope from single-W to curve).

3. **If only DUPLICATE-PLUS, why?**
   N/A — this session produced quantitative new data and a structural
   plateau finding, not a duplicate closure.

4. **Next-action for the next agent?**
   - **C8.a (next session, B-grade):** N=8 partial scan at W ∈ {1, 2,
     4, 8}, M ∈ {3..16} with ≥ 30-min CBC budget per cell. Save under
     `experiments/constructions/d2_sign_threshold_w_m_tradeoff/n8_extension/`.
   - **C8.b (next session, B-grade):** resolve random-control F4 at
     N=6 via S84-style column-enumeration extended to W ≥ 2 instead of
     direct ILP (column enumeration eliminates the bottom-layer weight
     variables that are causing the CBC slowdown).

## Diff from the originally-proposed C8 spec

The C8 spec said "N=8" but the N=8 ILP scaling pushed the work to N=6
(which is where S84 left a meaningful column-enumeration baseline).
The N=4 analog was added on top as a cross-scale validation. This is
documented in NOVELTY_CHALLENGES.md C8 as a "pivot from spec." The N=8
goal lives on as C8.a.

## Files written / modified

**New under `experiments/constructions/d2_sign_threshold_w_m_tradeoff/`:**
- `definition.md` — composition signature with edge IDs.
- `d2_w_m_tradeoff.py` — W-parameterised scan reusing S84 encoder.
- `d2_w_m_tradeoff_results.md` — pre-stated falsifiers + verdicts.
- `n4_grid.{json,log}`, `n6_primes_60s.{json,log}`,
  `n6_primes_highW.{json,log}`, `n6_primes_w24_m3.{json,log}`,
  `n6_primes_w3_m3.{json,log}`, `n6_rand_seed1.{json,log}`,
  `n6_rand_M3_seed1_long.{json,log}`.

**Modified tracking files:**
- `EDGES.md` — E5.3 annotated with S127 quantitative refinement.
- `NOVELTY_CHALLENGES.md` — C8 marked BUILT (S127); spawned C8.a, C8.b.
- `RESEARCH_AGENDA.md` — Arc 4 status updated to S127; C8 milestone
  ticked with full outcome paragraph; "next action" updated since
  C1-C8 now all BUILT.
- `status/CLOSED_PATHS.md` — row added under "Information Theory /
  Complexity Theory."
- `status/SESSION_INSIGHTS.md` — Session 127 entry appended.
