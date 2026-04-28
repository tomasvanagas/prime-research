# Session 144 — Arc 2 Lean Formalisation: W=10 corner + comprehensive leading-row enumeration

**Mode:** arc continuation (Arc 2 — Lean Formalisation Track).
**Date:** 2026-04-27.
**Self-grade:** **B-grade** (substantive arc step: closes a previously-
believed-obstructed corner via DP-based pre-search and a clean Lean
proof; provides a comprehensive leading-row enumeration that maps the
exhaustion of the existing technique).

## What was produced that wasn't in the project before

1. **Three new theorems in `experiments/formalisations/E2_1_mps_bond_dim/`:**
   - `chiP_ninety_seven_eq_one : chiP 97 = 1`
   - `exists_invertible_submatrix_W_eq_10_d_eq_j_plus_1 : ∀ j ≥ 1, ∃ ρ σ, IsUnit (...)`
   - `mps_bond_dim_W_eq_10_d_eq_j_plus_1 : ∀ j ≥ 1, (unfolding 10 (j+1) j).rank = 5`
   All sorry-free; `#print axioms` confirms only
   `[propext, Classical.choice, Quot.sound]`.

2. **Eleventh unconditional `mps_bond_dim` instance.** Closed W's are now
   `{2, 3, 4, 5, 6, 8, 10, 12, 18, 20}`. **First instance refuting an
   entry on the S128/S129 "structurally obstructed" list** — the prior
   sessions claimed W=10 was obstructed by a "multiplicity-2 residue
   pattern", but the S144 DP-based search found a clean triangulation
   that the earlier search apparently restricted out (by limiting to
   row prefixes `{0, 1, 2, 3, 4}` rather than `[0, W)`).

3. **DP-based leading-row enumeration script** at
   `experiments/formalisations/E2_1_mps_bond_dim/leading_row_search.py`
   (plus `leading_row_search_results.txt`). The script enumerates
   triangulations via bitmask DP over column subsets: `O(2^R · R · W)`
   per W. Comprehensive scan over `W ∈ [2, 72]` with `R ≤ 22` shows the
   leading-row family closes **exactly** the 10-element set above.
   Structurally obstructed at this parameter range:
   `{7, 9, 11, 13, 14, 15, 16, 17, 19, 21, 22, 24, 25, 26, 27, 28, 30,
   32, 33, 34, 36, 38, 40, 42, 44, 48, 50, 54, 60, 66}`.

4. **Negative shape edge documented**: leading-row + dead-col upper-
   triangulation closes `mps_bond_dim` orthogonal corner exactly for
   `W ∈ {2, 3, 4, 5, 6, 8, 10, 12, 18, 20}` at `R ≤ 22`. Filed in
   CLOSED_PATHS.md as `§L1 EXHAUSTED` row with full DP justification.

5. **Linter cleanup**: added the missing comment for the W=20
   `set_option maxHeartbeats 2000000 in` line per the linter's request.

## Edges composed / cited

- **E2.1** (the MPS bond dimension identity, the subject of L1).
- The closure depends on the existing `upper_bound` lemma (general)
  for the upper bound and on the corner-case prime exhibit (with rows
  `{0, 1, 3, 4, 9}` and dead column `1`) for the lower bound.

## Cross-domain technique

- Bitmask DP over column subsets for triangulation existence (algorithmic
  enumeration). The DP is a complete decision procedure for
  "does an `R × R` 0/1 matrix admit a row+col permutation to upper
  triangular with 1's on the diagonal?" — `reach[full] ⇔ ∃ permutation`.
- Lean side: same `det_of_upperTriangular` + `BlockTriangular id`
  skeleton as W=8/W=12.

## Triangulation data (W=10)

- **ρ = (1, 0, 4, 3, 9)** — mixed leading + non-leading rows.
- **σ = (8, 1, 2, 0, 6)** — dead col = 1.
- **Diagonal primes:** `{19, 2, 43, 31, 97}`.
- **Below-diagonal composites (10):** `{9, 32, 33, 39, 42, 49, 91, 92,
  93, 99}`.

The triangulation has `max_diag = 97`, well within `decide`-checkable
range for `Nat.Prime` (primes ≥ 150 typically need `norm_num`).

## What this does not do

- Does **not** close the general `exists_invertible_submatrix` `sorry`
  in `Basic.lean`. That remains the only outstanding obligation.
- Does **not** unblock `W ∈ {7, 9, 11, 13, 14, 15, 16, ..., 66}`. The
  S144 DP-search is a complete decision procedure for the leading-row
  family, and it confirms these W's are structurally obstructed. They
  need either:
  - `Matrix.det_of_blockTriangular` for non-triangulizable matrices
    (multi-session API development);
  - Cofactor-expansion-based determinant proofs (e.g., W=9's 7×7 sub-
    matrix has det = 1 by row cofactor expansion).

## Next-action

Two paths forward at the Arc 2 frontier:

1. **Block-triangular API sub-arc.** Develop `Matrix.det_of_blockTriangular`
   wrappers for the `W = 9` corner (7×7 with `1+3+3` block structure
   conjectured by S128/S129; the DP confirms no upper-triangulation but
   numpy confirms `det(rows 0..6) = 1` for the `R × R` sub-matrix with
   leading rows). Multi-session: 1 session of Lean infrastructure +
   1 session per corner closure.
2. **Cofactor-expansion approach.** For the W=9 corner specifically,
   the 7×7 sub-matrix of leading rows `{0..6}` and cols `{0, 1, 2, 3, 4,
   6, 7}` has det = 1. This can be proven in Lean by recursive cofactor
   expansion using `Matrix.det_succ_row_zero` — ~250-300 lines but
   self-contained (no new infrastructure).

Either path is multi-session; the leading-row family is now empirically
exhausted at `R ≤ 22`. Single-session leading-row corners are no longer
available.

## Process notes for future agents

- The DP-based reachability check is fast: full scan `W ∈ [2, 72]` with
  `R ≤ 22` runs in under 60s on a single core. It is a complete decision
  procedure for the leading-row family. **Do not rely on prior session
  notes** ("structurally obstructed") without re-running the DP — S128
  and S129 missed W=10's triangulation because their search was
  restricted to leading-row prefixes.
- For new corners `R > 8`, mathlib lacks `Fin.prod_univ_X`; we add a
  local `prod_univ_X'` lemma (cf. S143's `prod_univ_nine'`).
- The default `maxHeartbeats 200000` is insufficient for `R ≥ 9`
  (simp blow-up scales as `R²`); `set_option maxHeartbeats 2000000 in`
  works for R = 9 (S143). Higher R may need more.

## Self-evaluation

- **Q1: What did I produce that was not in the project before?**
  Three new Lean theorems (W=10 corner closure, sorry-free); one new
  primality helper (chiP 97); a comprehensive DP-based enumeration
  script that maps the leading-row family's coverage exactly; a
  CLOSED_PATHS row documenting the family's exhaustion; correction of
  the S128/S129 W=10 obstruction claim.

- **Q2: What edges did my work compose or cite?**
  E2.1 (the parent edge, refined inline by adding W=10 to the closed
  set).

- **Q3: If my session produced only duplicate closures, why?**
  Not applicable — produced a genuine new Lean theorem plus a
  substantive negative-shape map.

- **Q4: What is the next-action for the next agent?**
  Either pivot to block-triangular API development (multi-session) or
  attempt the W=9 corner via cofactor expansion (single-session
  ambitious). The leading-row family is exhausted at R ≤ 22.
