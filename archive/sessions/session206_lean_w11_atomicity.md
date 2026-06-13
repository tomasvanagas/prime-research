# Session 206 — Lean Arc 2: W=11 corner pre-search reveals atomic 5×5 odd block

**Mode:** ARC continuation (Arc 2 — Lean Formalisation Track).
**Date:** 2026-04-29.
**Targeted next-action (post-S159):** "W=11 closure via `1 + 5 + 5`
BlockTriangular split with local `Fin.prod_univ_five` lemma. Single-
session if pre-search yields ≤ 6 new prime helpers."

## Bottom line

Arc-anticipated single-session closure is **structurally invalidated**.
The W=11 11×11 matrix admits a parity-block decomposition into a 6×6
even block (which decomposes further as 1+5 with the inner 5×5 leading-
row triangulable) and a 5×5 odd block, and the 5×5 odd block is
**atomically irreducible**: zero block-triangular decompositions exist
across all 15 ordered partitions of 5 with parts ≤ 4. The arc next-
action has been replaced with three multi-session paths.

## What the session produced

1. **Parity-block decomposition** of the W=11 11×11 matrix. Under the
   permutation `(0, 2, 4, 6, 8, 10, 1, 3, 5, 7, 9)`, the matrix is
   BlockTriangular with `6×6 even-row × even-col` upper-left and `5×5
   odd-row × odd-col` lower-right, lower-left 5×6 zero by parity, single
   off-diagonal nonzero at `(row 0, col 1) = chiP(2) = 1`.

2. **Atomicity verification** of the W=11 odd 5×5 block. Exhaustive
   enumeration over all ordered partitions of 5 with parts ≤ 4 yields
   ZERO BlockTriangular decompositions. The block has det = 1, rank = 5,
   no row with fewer than 2 nonzero entries → no leading-row
   triangulation, no further sub-block decomposition.

3. **Single-session candidate for j ≥ 2 closure.** Inner 10×10 leading-
   row triangulation
   ```
   ρ ↦ (3, 2, 1, 6, 5, 7, 16, 10, 19, 18)
   σ ↦ (9, 8, 7, 6, 3, 5, 4, 2, 1, 0)
   ```
   over rows [1, 22), max row 19 (forces j ≥ 2 since 11^j ≥ 19 ⇔ j ≥ 2).
   Diagonal primes `{19, 31, 43, 59, 73, 83, 113, 181, 199, 211}`, 6 new
   prime helpers needed: `{67, 71, 73, 79, 83, 113, 181}`. Single-session
   viable for the j ≥ 2 part (matches arc's 6-helper budget) but yields
   a strictly weaker theorem (j ≥ 1 separate sub-problem).

4. **Three multi-session paths registered** in `RESEARCH_AGENDA.md` Arc 2
   next-action queue:
   * **Path A: j ≥ 2 corner + separate j = 1 cofactor expansion.**
     Single-session for j ≥ 2; ~2-3 additional sessions for j = 1 via
     parity decomposition + 5×5 cofactor expansion.
   * **Path B: develop reusable `det_fin_five` lemma**, analogous to
     mathlib's `det_fin_three`. Estimated 2 sessions; reusable for all
     future W needing 5×5 dets.
   * **Path C: pivot to W=14** (R=7, composite, lacks parity-atomicity
     issue). **Recommended for S207.**

## Files

* `experiments/formalisations/E2_1_mps_bond_dim/w11_blocktriangular_search.py`
* `experiments/formalisations/E2_1_mps_bond_dim/w11_inner_triangulation.py`
* `experiments/formalisations/E2_1_mps_bond_dim/w11_general_search.py`
* `experiments/formalisations/E2_1_mps_bond_dim/w11_nested_search.py`
* `experiments/formalisations/E2_1_mps_bond_dim/w11_odd_block_atomicity.py`
* `experiments/formalisations/E2_1_mps_bond_dim/w11_blocktriangular_search_results.md`

## Edges referenced

* **E2.1** (MPS bond-dim): the arc's L1 statement; the orthogonal corner
  family `d = j + 1` is its Lean instance set.

(No new edges; the atomicity finding is sub-edge content.)

## Self-grade: B-grade

Ambitious failure at the arc's stated single-session target, with a
clear structural diagnostic explaining the failure mode. Per CLAUDE.md
"Ambitious Failure is Encouraged" — a failed arc-target attempt with
structural rather than time-out failure is B-grade.

**Why not C.** The atomicity of the 5×5 odd block is a NEW structural
fact about W=11 (and likely all primes ≥ 11) that none of the prior
twelve W-closures encountered. It refines the post-S144 "block-
triangular-required" set into a strict sub-class.

**Why not A.** No Lean file produced, no theorem proven, no actual
closure. The contribution is to arc planning + future-agent guidance,
not the project's body of mathematical content.

## Self-evaluation (CLAUDE.md required)

1. **What did I produce that was not in the project before this
   session?** The parity-block decomposition of the W=11 11×11 matrix
   and the atomicity proof for the 5×5 odd block. Both are genuinely
   new structural facts, exhaustively verified by computer search. The
   new arc next-action queue (three concrete multi-session paths) is
   actionable for the next agent.

2. **What edges did my work compose or cite?** E2.1 (the arc's L1
   target). No edge added — the contribution is sub-E2.1.

3. **If my session produced only duplicate closures, why?** The session
   did not produce duplicate closures; it produced a structural
   refinement of the arc's anticipated approach.

4. **What is the next-action for the next agent?** Path C: pre-search
   for W=14 BlockTriangular partitions (use `w7_blocktriangular_search.py`
   adapted to W=14, R=7) and execute the corresponding Lean closure if a
   single-session candidate exists with ≤ 6 new prime helpers. If W=14
   also hits an obstruction, pivot to Path B (`det_fin_five` lemma).
