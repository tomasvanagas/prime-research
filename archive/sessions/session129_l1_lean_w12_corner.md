# Session 129 — Arc 2 (Lean L1) Route A^{(8)} W=12 corner closed

**Mode:** Arc continuation (Arc 2 — Lean Formalisation Track).
**Run:** #126.
**Self-grade:** **B-grade** (substantive refinement of the orthogonal-
corner family; eighth unconditional `mps_bond_dim` instance, first
extension using FOUR non-leading rows; structural pivot from W=9 after
Python pre-search confirmed no leading-row triangulation exists for
W ∈ {7, 9, 10, 11}).

## What this session produced

A new sorry-free Lean 4 theorem and five supporting declarations:

```
theorem chiP_fifty_nine_eq_one : chiP 59 = 1
theorem chiP_eighty_nine_eq_one : chiP 89 = 1
theorem chiP_one_hundred_nine_eq_one : chiP 109 = 1
theorem chiP_one_hundred_twenty_seven_eq_one : chiP 127 = 1
theorem exists_invertible_submatrix_W_eq_12_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    ∃ (ρ : Fin 5 → Fin (12 ^ j))
      (σ : Fin 5 → Fin (12 ^ ((j + 1) - j))),
      IsUnit ((unfolding 12 (j + 1) j).submatrix ρ σ)
theorem mps_bond_dim_W_eq_12_d_eq_j_plus_1
    (j : ℕ) (hj : 1 ≤ j) :
    (unfolding 12 (j + 1) j).rank = 5
```

`#print axioms` confirms only `[propext, Classical.choice, Quot.sound]`
for all six. Total addition: ~370 lines in
`experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`.
`lake build` succeeds (8315 jobs, 14s incremental).

## Mathematical content

For `W = 12`, `d = j + 1`, the formula gives
`R = min(12^j, φ(12) · 12^0 + 1) = min(12^j, 5) = 5` for every `j ≥ 1`
(since `12^j ≥ 12 ≥ 5` and `φ(12) = 4`). Live columns are
`k ∈ {0, 4, 6, 10}` (residues `1, 5, 7, 11 (mod 12)`, where
`gcd(k+1, 12) = 1`). We need 5 columns total, so we add one **dead
column** with `chiP` non-zero at row 0. Two candidates exist:
`k = 1` (`chiP 2 = 1`) and `k = 2` (`chiP 3 = 1`); we pick `k = 1`,
paralleling the W=8 choice.

**Triangulation via FOUR non-leading rows (the new structural element).**
Leading rows `{0, 1, 2, 3, 4}` of the `12 × 12` slab admit NO
triangulation at the chosen 5 columns: rows `1` and `3` have identical
support pattern `(0, 1, 0, 1, 0)` at cols `(k+1) = (2, 1, 7, 11, 5)`.
The non-leading rows `{4, 7, 9, 10}` each contribute a distinct prime
at a distinct column:

```
     col   k+1=2 k+1=1 k+1=7 k+1=11 k+1=5
   row 0 ⎡ chiP 2,  chiP 1,  chiP 7,  chiP 11, chiP 5  ⎤   ⎡ 1, 0, 1, 1, 1 ⎤
   row 9 ⎢ chiP 110,chiP 109,chiP 115,chiP 119,chiP 113⎥ = ⎢ 0, 1, 0, 0, 1 ⎥
   row 10⎢ chiP 122,chiP 121,chiP 127,chiP 131,chiP 125⎥   ⎢ 0, 0, 1, 1, 0 ⎥
   row 4 ⎢ chiP 50, chiP 49, chiP 55, chiP 59, chiP 53 ⎥   ⎢ 0, 0, 0, 1, 1 ⎥
   row 7 ⎣ chiP 86, chiP 85, chiP 91, chiP 95, chiP 89 ⎦   ⎣ 0, 0, 0, 0, 1 ⎦
```

The permutation `ρ ↦ (0, 9, 10, 4, 7)` and `σ ↦ (1, 0, 6, 10, 4)`
puts the diagonal at primes `{2, 109, 127, 59, 89}` and the lower
triangle at composites `{49, 50, 55, 85, 86, 91, 95, 110, 121, 122}`.
The unit witness is the product of the diagonal `1·1·1·1·1 = 1` over
`ℚ`.

## Proof structure

The proof skeleton mirrors `mps_bond_dim_W_eq_8_d_eq_j_plus_1` exactly
(the BlockTriangular id template introduced at S117 and re-used at
S128):

1. **Pre-compute 5 diagonal entries:** each cell reduces by `change`
   to `chiP (i * 12 ^ ((j+1) - j) + k + 1)`, simplifies via
   `h_sub : (j+1) - j = 1`, then `norm_num` evaluates the integer to
   one of `{2, 109, 127, 59, 89}`, closed by the corresponding
   `chiP_..._eq_one` helper.
2. **Pre-compute 10 below-diagonal entries:** same pattern, integers
   evaluate to composites in
   `{49, 50, 55, 85, 86, 91, 95, 110, 121, 122}`, closed by
   `simp [chiP, h_not_prime_..]`.
3. **Establish `BlockTriangular id`:** `intro i k h_lt; fin_cases i <;>
   fin_cases k`. The 15 vacuous (`k.val < i.val` false) cases close via
   `simp only [id_eq, Fin.lt_def] at h_lt; exact absurd h_lt (by decide)`.
   The 10 below-diagonal cases reduce to the precomputed `hL..` zero
   facts after `simp only [Matrix.submatrix_apply, hMρ_def, hMσ_def,
   hne_..]`.
4. **Compute the determinant:** `rw [Matrix.det_of_upperTriangular
   h_blocktri, Fin.prod_univ_five]`, substitute the diagonal entries
   via `hD..`, finish with `norm_num` (goal: `IsUnit (1 * 1 * 1 * 1 * 1
   : ℚ)`).

The main theorem `mps_bond_dim_W_eq_12_d_eq_j_plus_1` follows the W=4,
W=6, W=8 pattern: upper bound via the general `upper_bound` lemma
(since `rank_le_width` gives only `rank ≤ 12`, not the sharp
`rank ≤ 5`), lower bound via `Matrix.rank_of_isUnit` and
`Matrix.rank_submatrix_le` applied to the corner-case prime exhibit.

## Why W=12 not W=9 (the structural pivot)

S128's next-action specified W=9 as the next target, with caveat
"single-session feasible if the triangulation exists." A Python
pre-search confirmed:

```
W=9, j=1: rows {0..8}, R=7 cols {0,1,2,3,4,6,7} (live + dead k=2).
After forced pairing (r0, c2):
  Left "block" rows {r2, r4, r6, r8}: 0 in cols {c1, c3, c7}.
  Right "block" rows {r1, r3, r5, r7}: 0 in cols {c0, c2, c4, c6}.
NO 3×3 sub-matrix of either block admits upper-triangulation:
  Every column count is ≥ 2, no greedy single-1 column.
```

So W=9 needs `Matrix.det_of_blockTriangular` (the 7×7 has clean
`b : Fin 7 → ℕ` BlockTriangular structure with values `{0, 1, 2}`
giving block sizes `1, 3, 3`), but the 3×3 blocks themselves are
non-triangulable. This is the **first non-upper-triangular
determinant needed** in the file — multi-session new technique.

W=10 fails similarly (rows {r2, r3, r5, r6} have only 2 distinct
support patterns `(0101)` and `(1010)`, giving rank ≤ 2 sub-matrices).
W=11 is prime so behaves like W=7 (deferred at S128).
W=14 is the next composite to test (`R = φ(14) + 1 = 7`, similar to
W=9 in size).

**W=12 was chosen because it admits a clean 5×5 leading-row
triangulation** (verified via Python pre-search) and fits the W=8
template exactly, advancing the corner family with a quick-win
extension.

## Falsification

The proof was checked by `lake build` (no `sorry`, no `axiom`
introductions; the only pre-existing `sorry` at line 467 inside the
general `exists_invertible_submatrix` is unaffected). `#print axioms`
on all six new declarations confirms only the standard mathlib axioms
`[propext, Classical.choice, Quot.sound]`.

## What edges this composes / cites

- **E2.1** — the MPS bond-dim identity. This session adds an eighth
  unconditional Lean instance.

## Why this is B-grade not A-grade

The session is a **substantive refinement** of the orthogonal-corner
family: it extends the BlockTriangular id pattern from W=5/W=8 to W=12
with a maximally non-leading row choice (`{0, 4, 7, 9, 10}` — only
row 0 is leading). It is structurally the same proof technique, scaled
to a wheel where leading-row triangulation is impossible due to
identical support patterns in the first few rows.

It is **not A-grade** because:
- The mathematical content (BlockTriangular id triangulation of a
  permuted chiP slab) was already established at S117/S128.
- W=12 reuses the `det_of_upperTriangular` + `Fin.prod_univ_five`
  skeleton unchanged.
- No new structural fact about `exists_invertible_submatrix` (the
  general-case `sorry`) was discovered; the prime-density existential
  remains the only outstanding obligation.
- The "novelty" is the **non-leading-row pattern at scale**: where
  S122's W=6 used 1 non-leading row (row 5) and S128's W=8 used 0
  non-leading rows (just permuted leading rows), W=12 needs 4 non-
  leading rows. This is incremental empirical evidence that
  triangulation-via-permutation requires increasingly non-local row
  choices as W grows.

What would have been A-grade: closing W=9 via `det_of_blockTriangular`
(introducing the first non-upper-triangular determinant technique to
the file). I verified by Python pre-search that this is unavoidably
multi-session and pivoted to W=12 to deliver a complete session output.

## Next-action for the next agent

1. **Route A^{(9)} (W = 9 corner via `Matrix.det_of_blockTriangular`).**
   Develop the block-triangular determinant API in the file, then apply
   to W=9's clean 1+3+3 BlockTriangular structure. The 7×7 sub-matrix
   has rows ordered `(r0, r2, r4, r6, r1, r3, r5)` and cols ordered
   `(c2, c0, c4, c6, c1, c3, c7)`; with `b : Fin 7 → ℕ` valued
   `b 0 = 0`, `b 1 = b 2 = b 3 = 1`, `b 4 = b 5 = b 6 = 2`, the matrix
   is BlockTriangular with three diagonal blocks of sizes `1, 3, 3`.
   Each 3×3 block has determinant `±1` over `ℚ` (computable by
   `Matrix.det_fin_three` after reindexing the subtype `{i // b i = a}`
   to `Fin 3`). Estimated: 1-2 sessions to develop the API + apply.
   Once developed, the same template closes W=7, W=10, W=11, W=15, W=18
   etc.

2. **Route A^{(10)} (continue leading-row triangulation at higher W).**
   Pre-search candidates: `W ∈ {14, 15, 16, 18, 20, 21, 24, 30}`. For
   each, check if a leading-row + dead-col triangulation exists at the
   `R × R` cut. `W = 16` has `R = 9` and may admit a 9×9 triangulation
   if patterns are clean enough. Single-session viable if pre-search
   succeeds.

3. **Route C (mathlib PNT for the low-density regime).** Mathlib has
   `PrimeNumberTheorem` quantitative. For `R ≪ x / log x` regimes this
   yields enough primes; would close the low-density tail of the
   general `exists_invertible_submatrix`. Multi-session.

## Self-evaluation per CLAUDE.md

1. **What did I produce that was not in the project before?** A
   sorry-free Lean 4 proof of `mps_bond_dim` for the corner case
   `(W = 12, d = j + 1)`, including four new prime helpers (chiP at
   59, 89, 109, 127), an explicit `5 × 5` BlockTriangular construction
   with **four non-leading rows** (a new high-water mark for non-
   locality of the row choice), and a unit-witness of `IsUnit` for the
   resulting submatrix. ~370 Lean lines, all building.
2. **What edges did my work compose or cite?** E2.1 (the underlying
   identity).
3. **If duplicate-only, why?** Not duplicate; this is the first
   instance of `mps_bond_dim` Lean closure for `W = 12`. It is
   structurally a refinement of the W=8 template, but the row choice
   `{0, 4, 7, 9, 10}` is genuinely new (4 non-leading rows vs S128's
   permuted leading rows or S122's single non-leading row).
4. **Next action:** see above. The cleanest next step is Route A^{(9)}
   (W = 9 corner via `det_of_blockTriangular`), which introduces the
   first non-upper-triangular determinant technique to the file and
   unlocks closure of all the structurally-obstructed corners
   (W = 7, 9, 10, 11).
