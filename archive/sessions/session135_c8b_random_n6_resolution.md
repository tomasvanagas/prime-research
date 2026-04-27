# Session 135 — C8.b Construction: Random-Control F4 Resolution at N=6

**Date:** 2026-04-27
**Mode:** construction
**Target:** §1 NOVELTY_CHALLENGES C8.b
**Self-grade:** **B**

## Pick rationale

CLAUDE.md "Read ATTACK_VECTORS.md FIRST" preamble plus session-prompt
"Pick C_x (C1, C2, C3, C4, C5, or C6)". All §1 composition targets
C1-C8 are BUILT (S70, S71, S74, S81, S89, S105, S120, S127). The
remaining §1 work-units are *successors* of BUILT challenges: C3.a,
C4.a, C7.a, C8.a, C8.b. Among these, **C8.b** is the cleanest pickup
because (i) S127 left an *unresolved* F4 verdict at N=6 W ≥ 4 — the
random ILP cells were UNKNOWN at 600 s; (ii) the spec named a concrete
alternative encoding (column-enumeration extending S84's
`enum_d2_smart`) which the joint-ILP couldn't try; (iii) cost
estimate was ≤ 1 session.

## Composition

E5.3 × S84 column-enum framework × E1.6 / C7-S89 × C8/S127. The
question: does the S127 PRIMES-vs-random gate gap survive at W ≥ 2
once the joint-ILP UNKNOWN cells are resolved by an alternative
encoding?

## Object built

**Extended column enumeration** `Θ(N=6, W) = {distinct sign-threshold
truth tables on 6 bits with weights in {-W..W}}` extending S84's
W=1-only `enumerate_w1_thresholds_pruned`. Catalog sizes verified:
W=1: 1,458; W=2: 30,898; W=3: 218,066. Then
S84's `depth2_search` ILP at K=30,898 on PRIMES and density-matched
random (seeds 1, 5, 42) at M ∈ {2, 3, 4}.

## Pre-stated falsifiers (in `definition.md`)

F0/F0' sanity (PRIMES reproduces S84/S127); F1 (random easier at W=2
→ rejects F4); F2 (random ties PRIMES at W=2 → Δ = 0); F3 (random
strictly harder at W=2 → Δ ≥ 1, requires M=4 UNSAT proof); F4
(cross-seed robustness).

## Outcome

- **F0 HOLDS.** PRIMES W=1 K=1458: M=5 UNSAT (7 s), M=6 SAT (5 s) —
  matches S84.
- **F0' HOLDS.** PRIMES W=2 K=30898: M=3 UNSAT (157 s), M=4 SAT
  (181 s, gates=4, verified 64/64) — matches S127's joint-ILP `M*(W=2;
  6) = 4`. Column-enum 1.8× faster on M=3 UNSAT (157 vs 277 s) than
  S127's joint-ILP.
- **F1 REJECTED on every seed.** Random N=6 W=2 has M=2 UNSAT
  (130–196 s) AND M=3 UNSAT (147–230 s) at all three seeds {1, 5, 42}.
- **`M*(rand_s; W=2) ≥ 4 = M*(PRIMES; W=2)`** for s ∈ {1, 5, 42}.
- **F2 partially established** (lower bound matches PRIMES). **F3
  UNRESOLVED**: random M=4 W=2 seed=42 returned UNKNOWN at 618 s
  (CBC time-limit at 600 s, marginally exceeded).
- **F4 (cross-seed robustness) HOLDS** at all three tested seeds.

## Net new content

1. **First quantitative resolution of S127's open random N=6 cell at
   W ≥ 2.** Column enumeration reduces the random `(W=2, M=3)` UNSAT
   proof from "UNKNOWN at 600 s" (S127 joint ILP) to **UNSAT proven
   in 147–230 s**. The cross-encoding shift (column-enum vs joint)
   eliminates alpha-bilinear constraints `v[k] = sel[k] AND beta[k]`
   that dominate joint-ILP search.

2. **F4 holds at W=2 with Δ ≥ 0 robustly across 3 seeds.** Density-
   matched random is *never easier* than PRIMES at W=2; the strict-
   inequality magnitude (Δ = 0 vs Δ ≥ 1) remains open.

3. **Catalog size empirical fact:** `#Θ(6, W=2) = 30,898`,
   `#Θ(6, W=3) = 218,066`. Column-enum at W=2 is tractable at M ≤ 3
   in 130–230 s; intractable at M = 4.

4. **Methodological side-finding for the S84 framework.** Column-enum
   dominates joint-ILP for low-M UNSAT proofs (1.8× faster on PRIMES
   W=2 M=3; resolves random M=3 cells joint-ILP couldn't touch in
   600 s). Joint-ILP remains the right tool for M ≥ 4 SAT search
   (PRIMES W=2 M=4 SAT in 17 s joint vs 181 s column-enum).

## Refines

- **E5.3** with the new statement: PRIMES `M*(W=2; N=6) = 4` is *not*
  breakable by replacing PRIMES with density-matched random — the
  random class is at least as hard at W=2.
- **C7-S89 (E1.6)** with: even *outside* the calibrated-oddness regime
  (where density+bit_0 control absorbs the gap at W=1), the W=2 gap
  holds in direction across multiple seeds.
- **C8/S127** F4 verdict at N=6: from "unresolved" to "direction-
  confirmed across 3 seeds; magnitude between 0 and 1 gate."

## Closure verdict

Mode E (extended measurement of S127 with new ILP encoding). BUILT, no
polylog opening. F4 lower-bound established at W=2 across three random
seeds; strict-magnitude question unresolved.

## Cross-domain ingredient

NONE. The technique innovation is *cross-encoding* (column-enum vs
joint-ILP), not cross-domain. Per CLAUDE.md "Cross-Domain Imports"
this is a B-grade target by construction.

## Successors proposed

- **C8.b.i.** Greedy/LNS SAT search on random M=4 W=2. CBC's UNSAT proof
  at this scale is intractable; a local search may find a witness if
  one exists. Empirical evidence for `M*(rand; W=2) ≥ 5` (Δ ≥ 1) if
  no SAT in ~1 h CPU. Cost: 1 session.
- **C8.b.ii.** W=3 column-enum on random K=218,066. PRIMES W=3 M=3 SAT
  (S127, 65 s); does random M=3 stay UNSAT at W=3? Strengthens F4 by
  one weight-doubling step. Cost: 1-2 sessions.
- **C8.b.iii.** Seed-distribution histogram at N=6 W=2 M=3 across 100+
  seeds in 32-core parallel ≈ 10 min wall-clock. Tests "always ≥ 4" vs
  "sometimes 3" empirically. Cost: 1 session.

## Self-evaluation

1. **What did I produce that was not in the project before?**
   - The catalog `Θ(N=6, W=2)` of 30,898 distinct sign-threshold truth
     tables (and the W=3 catalog at 218,066, computed but not run).
   - The ILP-proven lower bound `M*(rand_s; W=2) ≥ 4` for three random
     seeds, resolving an S127 UNKNOWN cell.
   - The cross-encoding methodological observation (column-enum
     dominates joint-ILP for low-M UNSAT bounds).
2. **What edges did my work compose or cite?** E5.3 (refined inline),
   E1.6 / C7-S89 (refined direction at W=2), C8/S127 (F4 verdict
   sharpened). Edge IDs cited in CLOSED_PATHS row and EDGES.md
   annotation.
3. **Honest grade:** B. The session resolved a previously-open ILP
   cell with a real technique change (cross-encoding), produced a new
   empirical lower bound across three seeds, and identified concrete
   methodological dominance (column-enum > joint-ILP for low-M UNSAT).
   The strict-magnitude F3 question (Δ = 0 vs Δ ≥ 1) is left open;
   that's a real limitation, not a packaging issue. Not A-grade —
   this is a refinement of an existing edge with a precise new
   statement, not a frontier-attack partial positive. Not C-grade —
   the cross-encoding shift IS new structural content, and the
   M*(rand) ≥ 4 result robustly across 3 seeds is non-obvious *ex
   ante* (and S127 explicitly couldn't establish it).
4. **Next action for next agent:** Pick C8.b.iii (seed-distribution
   histogram, 32-core parallelisable, ≈ 10 min wall-clock cost). If
   the histogram shows even one random seed at M=3, F4 itself is
   fragile and we need to reconsider whether the gap is structural or
   seed-dependent. Or, if a stronger statement is preferred, C8.b.ii
   at W=3 — substantially more compute but a tighter F4 statement.
