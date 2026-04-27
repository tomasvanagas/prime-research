# Session 128 — Arc 2 (Lean L1) Route A''''''' W=8 corner closed

**Mode:** Arc continuation (Arc 2 — Lean Formalisation Track).
**Run:** #119.
**Self-grade:** **B-grade** (substantive refinement of the orthogonal-
corner family; seventh unconditional `mps_bond_dim` instance, first
extension of the BlockTriangular route to a wheel `W ≥ 6` where
`R < W`).

## What this session produced

A new sorry-free Lean 4 theorem and three supporting declarations:

```
theorem chiP_seventeen_eq_one : chiP 17 = 1
theorem chiP_thirty_seven_eq_one : chiP 37 = 1
theorem exists_invertible_submatrix_W_eq_8_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    ∃ (ρ : Fin 5 → Fin (8 ^ j))
      (σ : Fin 5 → Fin (8 ^ ((j + 1) - j))),
      IsUnit ((unfolding 8 (j + 1) j).submatrix ρ σ)
theorem mps_bond_dim_W_eq_8_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    (unfolding 8 (j + 1) j).rank = 5
```

`#print axioms` confirms only `[propext, Classical.choice, Quot.sound]`
for all four. Total addition: ~280 lines in
`experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`.
`lake build` succeeds (8315 jobs, 14s incremental).

## Mathematical content

For `W = 8`, `d = j + 1`, the formula gives
`R = min(8^j, φ(8) · 8^0 + 1) = min(8^j, 5) = 5` for every `j ≥ 1`.
Live columns are `k ∈ {0, 2, 4, 6}` (residues `1, 3, 5, 7 (mod 8)`,
where `gcd(k+1, 8) = 1`); pick `R = φ(8) + 1 = 5` columns by adding
the unique dead column `k = 1` with `chiP 2 = 1` at row 0.

The first 5 rows `{0, 1, 2, 3, 4}` restricted to the 5 chosen columns
`{0, 1, 2, 4, 6}` form a `5 × 5` matrix with the support pattern

```
     col 0  1  2  4  6
row 0  0  1  1  1  1
row 1  0  0  1  1  0
row 2  1  0  1  0  1
row 3  0  0  0  1  1
row 4  0  0  0  1  0
```

This is rank 5 (each row contains a 1 not present in any earlier row,
after suitable column re-ordering). The permutation
`ρ ↦ (2, 0, 1, 3, 4)` and `σ ↦ (0, 1, 2, 6, 4)` triangularises it:
diagonal primes are `{17, 2, 11, 31, 37}`, lower-triangle composites
are `{1, 9, 10, 25, 26, 27, 33, 34, 35, 39}`. The unit witness is the
product of the diagonal `1·1·1·1·1 = 1` over `ℚ`.

## Proof structure

The proof skeleton mirrors `mps_bond_dim_W_eq_5_d_eq_j_plus_1` exactly
(the BlockTriangular template introduced at S117):

1. **Pre-compute 5 diagonal entries:** each cell reduces by `change`
   to `chiP (i * 8 ^ ((j+1) - j) + k + 1)`, simplifies via `h_sub :
   (j+1) - j = 1`, then `norm_num` evaluates the integer to one of
   `{2, 11, 17, 31, 37}`, closed by the corresponding `chiP_..._eq_one`
   helper.
2. **Pre-compute 10 below-diagonal entries:** same pattern, integers
   evaluate to composites in `{1, 9, 10, 25, 26, 27, 33, 34, 35, 39}`,
   closed by `simp [chiP, h_not_prime_..]` after a `decide` step.
3. **Establish `BlockTriangular id`:** `intro i k h_lt; fin_cases i <;>
   fin_cases k`. The 15 vacuous (`k.val < i.val` false) cases close via
   `simp [id_eq, Fin.lt_def] at h_lt; exact absurd h_lt (by decide)`.
   The 10 below-diagonal cases reduce to the precomputed `hL..` zero
   facts after `simp only [Matrix.submatrix_apply, hMρ_def, hMσ_def,
   hne_..]`.
4. **Compute the determinant:** `rw [Matrix.det_of_upperTriangular
   h_blocktri, Fin.prod_univ_five]`, substitute the diagonal entries
   via `hD..`, finish with `norm_num` (goal: `IsUnit (1 * 1 * 1 * 1 * 1
   : ℚ)`).

The main theorem `mps_bond_dim_W_eq_8_d_eq_j_plus_1` follows the W=4
and W=6 pattern: upper bound via the general `upper_bound` lemma
(since `rank_le_width` gives only `rank ≤ 8`, not the sharp `rank ≤ 5`),
lower bound via `Matrix.rank_of_isUnit` and `Matrix.rank_submatrix_le`
applied to the corner-case prime exhibit.

## Falsification

The proof was checked by `lake build` (no `sorry`, no `axiom`
introductions; the only pre-existing `sorry` at line 467 inside
`exists_invertible_submatrix` is unaffected). `#print axioms`
confirms only the standard mathlib axioms.

## What edges this composes / cites

- **E2.1** — the MPS bond-dim identity. This session adds a seventh
  unconditional Lean instance.

## Why this is B-grade not A-grade

The session is a **substantive refinement** of the orthogonal-corner
family: it extends the BlockTriangular pattern from W=5 (where `R = W`)
to W=8 (where `R < W`, with one dead column added to the live-column
set). It is structurally the same proof technique, scaled to a slightly
different live/dead-column ratio.

It is **not A-grade** because:
- The mathematical content (BlockTriangular triangulation of the chiP
  slab) was already established at S117 for W=5.
- W=8 reuses the `det_of_upperTriangular` + `Fin.prod_univ_five` skeleton
  unchanged.
- No new structural fact about `exists_invertible_submatrix` (the
  general-case `sorry`) was discovered; the prime-density existential
  remains the only outstanding obligation.
- The "novelty" is purely empirical: confirming that a particular row/
  column permutation (found by hand) triangularises the 5×5 chiP slab
  at W=8.

## Why W=7 was deferred

I considered W=7 first since it was named in the next-action; for `W = 7`,
`R = min(7^j, φ(7) + 1) = 7` (since `7` is prime), so we'd need a
`7 × 7` invertible submatrix.

Empirical analysis of the first 7 rows of the W=7 slab shows the
support pattern admits exactly two perfect matchings of 1-entries:
- M1: `(0→6, 1→3, 2→2, 3→1, 4→0, 5→5, 6→4)`
- M2: `(0→6, 1→3, 2→4, 3→1, 4→2, 5→5, 6→0)`

Both matchings reduce to a residual `3 × 3` sub-issue on rows
`{2, 4, 6}` × cols `{0, 2, 4}` (or equivalent), where each remaining
column has exactly two 1's in our chosen row set — no column has a
single 1, so triangulation cannot proceed via the standard
"position 0 = column with single 1" greedy algorithm.

The 3×3 sub-issue submatrix
```
            col 0  2  4
   row 2     0  1  1       chiP 15, chiP 17, chiP 19
   row 4     1  1  0       chiP 29, chiP 31, chiP 33
   row 6     1  0  1       chiP 43, chiP 45, chiP 47
```
has determinant `0·(1-0) - 1·(1-0) + 1·(0-1) = -2 ≠ 0`, so the 7×7
matrix IS invertible — but no row/column permutation triangularises
it. Closing W=7 cleanly therefore needs either (a) `Matrix.det_fin_seven`
(absent from mathlib), (b) a manual Laplace expansion or row reduction,
or (c) the use of a non-leading row set for `j ≥ 2` (analogous to
W=6's row 5 trick) — multi-session work. W=8 fits in one session;
W=7 does not.

## Next-action for the next agent

1. **Route A^{(8)} (W = 9 corner).** `R = min(9^j, φ(9) + 1) = 7` for
   `j ≥ 1`. Live cols `{0, 1, 3, 4, 6, 7}` plus one dead col with
   `chiP` non-zero (e.g. col `2`, `chiP 3 = 1`). Need a `7 × 7`
   triangulation. Mathlib lacks `Fin.prod_univ_seven` directly but
   `Finset.prod_univ_succ` chains can build it; alternatively
   `Matrix.det_succ_column_zero` recursively expands. Single-session
   feasible if the triangulation exists.
2. **Route A^{(9)} (W = 7 corner, multi-session).** As analysed above,
   the leading-row triangulation fails. Either (a) accept a Laplace-
   expansion proof for the 7×7 j=1 matrix, or (b) for `j ≥ 2`, use a
   non-leading row that breaks the residual cycle (e.g. row 22 contains
   only the prime 157, contributing a clean column-2 entry).
3. **Route C (mathlib PNT).** Mathlib has `PrimeNumberTheorem` quantitative.
   For `R ≪ x / log x` regimes this yields enough primes; would close
   the low-density tail of the general `exists_invertible_submatrix`.
4. **Route A** (full Hoheisel-grade primes-in-short-intervals). Multi-
   session arc by itself.

## Self-evaluation per CLAUDE.md

1. **What did I produce that was not in the project before?** A
   sorry-free Lean 4 proof of `mps_bond_dim` for the corner case
   `(W = 8, d = j + 1)`, including two new prime helpers (chiP at
   17 and 37), an explicit `5 × 5` BlockTriangular construction, and
   a unit-witness of `IsUnit` for the resulting submatrix. ~280 Lean
   lines, all building.
2. **What edges did my work compose or cite?** E2.1 (the underlying
   identity).
3. **If duplicate-only, why?** Not duplicate; this is the first
   instance of `mps_bond_dim` Lean closure for `W = 8`. It is
   structurally a refinement of the W=5 template, but the specific
   permutation and the row/column choices are new.
4. **Next action:** see above. The cleanest single-session next step
   is Route A^{(8)} (W = 9), assuming a 7×7 triangulation can be found.
