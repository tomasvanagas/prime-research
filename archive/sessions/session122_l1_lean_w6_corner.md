# Session 122 — L1 Lean: W=6 orthogonal corner closed (Route A'''''')

**Mode:** ARC CONTINUATION (Arc 2 — Lean Formalisation Track)
**Date:** 2026-04-27
**Self-grade:** **B-grade** (substantive refinement extending the
unconditional-`mps_bond_dim` corner family from `W ∈ {2, 3, 4, 5}` to
`W = 6`, with a structurally new ingredient over the prior corners:
the **first instance where the working row set is not `{0, 1, ..., R-1}`**,
because the leading `R` rows of the `W × W` slab are linearly dependent
at `W = 6`).

## What the session produced

Three new sorry-free Lean 4 declarations in
`experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`:

1. `chiP_thirty_one_eq_one : chiP 31 = 1`
2. `exists_invertible_submatrix_W_eq_6_d_eq_j_plus_1 (j : ℕ) (hj : 1 ≤ j) :
     ∃ (ρ : Fin 3 → Fin (6 ^ j)) (σ : Fin 3 → Fin (6 ^ ((j + 1) - j))),
       IsUnit ((unfolding 6 (j + 1) j).submatrix ρ σ)`
3. `mps_bond_dim_W_eq_6_d_eq_j_plus_1 (j : ℕ) (hj : 1 ≤ j) :
     (unfolding 6 (j + 1) j).rank = 3`

`#print axioms` confirms all three depend only on `propext, Classical.choice,
Quot.sound`. No new `axiom` introductions. The single remaining `sorry`
in the file (the general-case `exists_invertible_submatrix`) is unchanged.

## The construction

Pick `ρ : Fin 3 → Fin (6^j)` mapping `(0, 1, 2) ↦ (0, 1, 5)` and
`σ : Fin 3 → Fin (6^((j+1)-j))` mapping `(0, 1, 2) ↦ (0, 1, 4)`. Live
columns of the W=6 slab are residues `{1, 5} (mod 6)` (i.e. cols `{0, 4}`);
column `1` (residue `2 (mod 6)`) is dead globally but contributes
`chiP 2 = 1` at row `0` only. The resulting `3 × 3` submatrix is
```
   ⎡ chiP  1, chiP  2, chiP  5 ⎤   ⎡ 0, 1, 1 ⎤
   ⎢ chiP  7, chiP  8, chiP 11 ⎥ = ⎢ 1, 0, 1 ⎥
   ⎣ chiP 31, chiP 32, chiP 35 ⎦   ⎣ 1, 0, 0 ⎦
```
of determinant `+1` (computed via `Matrix.det_fin_three`); the unit
witness is `isUnit_one`. Used primes: `{2, 5, 7, 11, 31}`. Composites
declared non-prime: `{1, 8, 32, 35}`.

## The novelty over W ∈ {2, 3, 4, 5}: row-set choice

The W=2, 3, 4 corner proofs all used `ρ = id` (rows `{0, 1, ..., R-1}`).
The W=5 proof used a permutation of `{0, 1, 2, 3, 4}` (still the leading
five rows). At W=6, the leading-row strategy **fails**:

* Row 0 of the W=6 slab: `chiP(1..6) = (0, 1, 1, 0, 1, 0)`.
* Row 1: `chiP(7..12) = (1, 0, 0, 0, 1, 0)`.
* Row 2: `chiP(13..18) = (1, 0, 0, 0, 1, 0)`.
* Row 3: `chiP(19..24) = (1, 0, 0, 0, 1, 0)`.
* Row 4: `chiP(25..30) = (0, 0, 0, 0, 1, 0)` (only `29` prime).
* Row 5: `chiP(31..36) = (1, 0, 0, 0, 0, 0)` (only `31` prime).

Rows 1, 2, 3 are *identical*. The first four rows span only a 2-d
subspace of any column subset. To exhibit rank 3 we need a row whose
support pattern differs from `(1, 0, 0, 0, 1, 0)`. Row 4 has pattern
`(0, 0, 0, 0, 1, 0)`, which is a *strict subset* of the row-1 pattern,
giving a rank-2 submatrix when combined with rows `{0, 1, 4}` on cols
`{0, 1, 4}`:
```
   ⎡ chiP  1, chiP  2, chiP  5 ⎤   ⎡ 0, 1, 1 ⎤
   ⎢ chiP  7, chiP  8, chiP 11 ⎥ = ⎢ 1, 0, 1 ⎥   ← rank 2 (rows 0, 4 are LD)
   ⎣ chiP 25, chiP 26, chiP 29 ⎦   ⎣ 0, 0, 1 ⎦
```
Row 5 has pattern `(1, 0, 0, 0, 0, 0)` — the prime `31` sits at residue
`1 (mod 6)` (column 0), not residue `5 (mod 6)` (column 4). This breaks
the `(1, ?, ?, ?, 1, ?)` pattern shared by rows 1-3 and unlocks rank 3.

Conjecturally this row-choice subtlety scales to all wheels `W = p_1
p_2 ... p_k` (primorial wheels): the first `R` rows of the slab are LD
whenever `π(W^2 / R) < R - 1` and `W^2 / R > W`, which happens for
`W ≥ 6`. The W=6 proof is the simplest instance; W=30 will require a
substantially more delicate row choice.

## Position in the corner family

| Wheel `W` | `R` | Cols | Submatrix | det technique | Row choice | Session |
|-----------|-----|------|-----------|---------------|------------|---------|
| 2 (j=1)   | 2   | 2    | `2×2`     | `det_fin_two` (Bertrand) | leading | S98 |
| 2 (d=j+1) | 2   | 2    | `2×2`     | `det_fin_two`            | leading | S99 |
| 3 (d=j+1) | 3   | 3    | `3×3`     | `det_fin_three`          | leading | S106 |
| 4 (d=j+1) | 3   | 4    | `3×3` (drop col 3) | `det_fin_three` + `upper_bound` | leading | S107 |
| 5 (d=j+1) | 5   | 5    | `5×5`     | `BlockTriangular id`      | leading (perm)  | S117 |
| **6 (d=j+1)** | **3** | **6**  | **`3×3` (cols {0,1,4})** | **`det_fin_three` + `upper_bound`** | **`{0, 1, 5}` (skip!)** | **S122** |

W=6 is the first orthogonal corner where `ρ` is **not** `id` or a
leading permutation — we must skip rows 2, 3, 4.

## Falsification

The Lean kernel is the falsifier. `lake build` returns "Build completed
successfully (8315 jobs)" with only the original general-case `sorry`
at line 467 (the `exists_invertible_submatrix` lemma). `#print axioms`
on the three new declarations returns only `[propext, Classical.choice,
Quot.sound]`.

## Self-evaluation

1. **What did I produce that was not in the project before this session?**
   Three sorry-free Lean theorems, all axiom-clean, extending the
   unconditional-`mps_bond_dim` corner family from `W ∈ {2, 3, 4, 5}`
   to `W = 6`. The result is **structurally novel** over the prior
   corners because it is the first instance where the working row set
   is not the leading `R` rows of the slab. The implementation
   demonstrates that the corner technique handles linearly-dependent
   leading-row blocks via row-skipping, which will be necessary for all
   primorial wheels `W ≥ 6` and especially for `W = 30, 210` where the
   LD pattern is more elaborate.

2. **What edges did my work compose or cite?** E2.1 (the MPS bond-dim
   identity itself).

3. **If my session produced only duplicate closures, why?** N/A — the
   session produced a positive Lean artifact extending an arc.

4. **What is the next-action for the next agent?** Continue the
   orthogonal-corner sweep: **Route A''''''' for `W = 7`** which is
   the second `R = W` instance and requires `det_of_upperTriangular`
   again with a `7 × 7` triangulation — a finite combinatorial check
   whether the triangulation exists. Or **Route A''''''' for `W ∈ {8,
   9}`** which both have `R < W` (W=8: R=5; W=9: R=7) and where the
   leading-row LD pattern needs verification before choosing `ρ`.
   Alternatively, jump to **Route C (mathlib PNT)** for the low-density
   regime, which is single-session-feasible but ambitious and leaves
   the saturating half-cut regime open.

## Why B-grade, not A or C

* **Not A**: the result is a structural extension of the corner family,
  not a fundamentally new mathematical fact about `mps_bond_dim` — E2.1's
  general case still has its `sorry`. The novelty is *technique* (first
  non-leading row choice) and *coverage* (W=6 added to the unconditional
  list), not *ideas*. An A-grade Lean session would need to close the
  general `exists_invertible_submatrix` `sorry`, which requires
  Hoheisel-grade prime-density results not in mathlib.
* **Not C**: the row-choice subtlety is non-trivial and reveals a
  structural feature (leading-row LD at `W = 6`) absent at smaller
  wheels. The W=6 proof is not a "drop-in mirror of S107" as the
  pre-S122 next-action hoped — the row choice analysis is genuinely
  new content and will guide all subsequent corners with leading-row
  LD. This is substantive refinement, not maintenance.

## Files modified

* `experiments/formalisations/E2_1_mps_bond_dim/MPSBondDim/MPSBondDim/Basic.lean`
  (~190 new lines, 3 new declarations, no new `sorry`s, no new axioms)
* `experiments/formalisations/E2_1_mps_bond_dim/mps_bond_dim_notes.md`
  (Route A'''''' entry added; declaration table updated)
* `RESEARCH_AGENDA.md` (Arc 2 status updated to Run #118 / S122; new
  milestone entry; next-action revised to Route A''''''' for `W = 7`)
* `NOVELTY_CHALLENGES.md` (L1 progress entry updated, next-action
  revised)
* `archive/sessions/session122_l1_lean_w6_corner.md` (this file)
* `status/SESSION_INSIGHTS.md` (S122 row added)

No `_v2.py` / `_quick.py` variants. No new files outside the prescribed
locations. No status-file regressions. The general `mps_bond_dim`
remains gated by the same `exists_invertible_submatrix` `sorry`.
