# L1 — E2.1 MPS bond-dim identity (Lean 4 formalisation, in progress)

**Lean source:** `MPSBondDim/MPSBondDim/Basic.lean`
**Toolchain:** `leanprover/lean4:v4.30.0-rc2` + mathlib `v4.30.0-rc2` (lake project under `MPSBondDim/`).
**Build status (S152):** `lake build` succeeds. **1 `sorry` placeholder remains, isolated to a pure prime-density existential** (S83 reorganised the file: `lower_bound` is now closed by a structural reduction to `exists_invertible_submatrix`, which is the new home of the only `sorry`. S90 added the trivial floor `1 ≤ rank` as a separate lemma. **S98 closed the corner case `(W = 2, j = 1)` unconditionally via Bertrand** — `mps_bond_dim_W_eq_2_j_eq_1 : (unfolding 2 d 1).rank = 2` for every `d ≥ 2`, sorry-free. **S99 closed the orthogonal corner case `(W = 2, d = j + 1)` *without even needing Bertrand*** — `mps_bond_dim_W_eq_2_d_eq_j_plus_1 : (unfolding 2 (j+1) j).rank = 2` for every `j ≥ 1`, using only `Nat.prime_two` and `Nat.prime_three`. **S106 extended the orthogonal-corner argument to `W = 3`** — `mps_bond_dim_W_eq_3_d_eq_j_plus_1 : (unfolding 3 (j+1) j).rank = 3` for every `j ≥ 1`, using `Matrix.det_fin_three` and the explicit primes `2, 3, 5, 7`; the unit witness is `isUnit_one.neg : IsUnit (-1 : ℚ)`. **S107 extended the orthogonal-corner argument to `W = 4`** — `mps_bond_dim_W_eq_4_d_eq_j_plus_1 : (unfolding 4 (j+1) j).rank = 3` for every `j ≥ 1`, using `Matrix.det_fin_three`, the explicit primes `2, 3, 5, 7, 11`, and (for the upper bound) the general `upper_bound` lemma — the first orthogonal-corner instance where `rank_le_width` is not tight (it gives only `rank ≤ 4`, not the sharp `rank ≤ 3 = φ(4) · 4^0 + 1`). **S117 extended the orthogonal-corner argument to `W = 5`** — `mps_bond_dim_W_eq_5_d_eq_j_plus_1 : (unfolding 5 (j+1) j).rank = 5` for every `j ≥ 1`, using `Matrix.det_of_upperTriangular` (since mathlib has no `det_fin_four` or `det_fin_five`), `Fin.prod_univ_five`, the explicit primes `2, 3, 5, 7, 11, 19, 23`, and a row/column permutation that triangularises the `chiP 1..25` slab. **First instance with `R = W` (all `W` columns needed) and first instance not relying on `det_fin_three`.** **S122 extended the orthogonal-corner argument to `W = 6`** — `mps_bond_dim_W_eq_6_d_eq_j_plus_1 : (unfolding 6 (j+1) j).rank = 3` for every `j ≥ 1`, using `Matrix.det_fin_three`, the explicit primes `2, 3, 5, 7, 11, 31`, and (for the upper bound) the general `upper_bound` lemma. **First orthogonal-corner instance where the working row set is not `{0, 1, ..., R-1}`** — rows `{1, 2, 3}` of the `6 × 6` slab are linearly dependent (all three windows `chiP 7..12, 13..18, 19..24` have identical support pattern `(1, 0, 0, 0, 1, 0)`), forcing the choice `ρ ↦ (0, 1, 5)` with `chiP 31` providing the third linearly independent row. **S128 extended the orthogonal-corner argument to `W = 8`** — `mps_bond_dim_W_eq_8_d_eq_j_plus_1 : (unfolding 8 (j+1) j).rank = 5` for every `j ≥ 1`, using `Matrix.det_of_upperTriangular`, `Fin.prod_univ_five`, the explicit primes `2, 11, 17, 31, 37`, and (for the upper bound) the general `upper_bound` lemma (since `rank_le_width` gives only `rank ≤ 8`, not the sharp `rank ≤ 5 = φ(8) · 8^0 + 1`). The triangulation `ρ ↦ (2, 0, 1, 3, 4)` and `σ ↦ (0, 1, 2, 6, 4)` puts the diagonal at primes `{17, 2, 11, 31, 37}` and the lower triangle at composites `{1, 9, 10, 25, 26, 27, 33, 34, 35, 39}`. **Seventh unconditional instance; fifth instance over a wheel `W ≥ 3`; second instance using `det_of_upperTriangular`.** **S129 extended the orthogonal-corner argument to `W = 12`, skipping the structurally-obstructed `W ∈ {7, 9, 10, 11}` corners** — `mps_bond_dim_W_eq_12_d_eq_j_plus_1 : (unfolding 12 (j+1) j).rank = 5` for every `j ≥ 1`, using `Matrix.det_of_upperTriangular`, `Fin.prod_univ_five`, the explicit primes `2, 59, 89, 109, 127`, and (for the upper bound) the general `upper_bound` lemma. The triangulation `ρ ↦ (0, 9, 10, 4, 7)` and `σ ↦ (1, 0, 6, 10, 4)` puts the diagonal at primes `{2, 109, 127, 59, 89}` and the lower triangle at composites `{49, 50, 55, 85, 86, 91, 95, 110, 121, 122}`. **Eighth unconditional instance; sixth instance over a wheel `W ≥ 3`; third instance using `det_of_upperTriangular`.** **First instance using FOUR non-leading rows** (only row `0` is leading) — extends S122's W=6 single-non-leading-row trick to the maximally non-leading regime. **S137 extended the orthogonal-corner argument to `W = 18`** — `mps_bond_dim_W_eq_18_d_eq_j_plus_1 : (unfolding 18 (j+1) j).rank = 7` for every `j ≥ 1`, using `Matrix.det_of_upperTriangular`, `Fin.prod_univ_seven`, the explicit primes `2, 29, 43, 109, 179, 211, 293`, and (for the upper bound) the general `upper_bound` lemma (since `rank_le_width` gives only `rank ≤ 18`, not the sharp `rank ≤ 7 = φ(18) · 18^0 + 1`). The triangulation `ρ ↦ (0, 2, 9, 1, 11, 6, 16)` and `σ ↦ (1, 6, 16, 10, 12, 0, 4)` puts the diagonal at primes `{2, 43, 179, 29, 211, 109, 293}` and the lower triangle at composites `{20, 25, 35, 38, 110, 115, 119, 121, 125, 164, 169, 200, 205, 209, 215, 289, 290, 295, 299, 301, 305}`. **Ninth unconditional instance; seventh instance over a wheel `W ≥ 3`; fourth instance using `det_of_upperTriangular`.** **First instance with `R = 7`.** **W=14 (also `R = 7`) was tried and structurally obstructed**: the 14×14 j=1 slab admits no leading-row triangulation with `ρ < 14` (rows 2 and 5 of the W=14 slab have identical support pattern at the chosen 7 cols, and exhaustive search over rows in `[0, 14)` finds zero upper-triangulations). The general-case `sorry` is unchanged.)
**S152 closed the orthogonal corner `(W = 9, d = j + 1)` — the FIRST closure of an S128/S129/S144 "block-triangular-required" wheel.** `mps_bond_dim_W_eq_9_d_eq_j_plus_1 : (unfolding 9 (j+1) j).rank = 7` for every `j ≥ 1`, sorry-free. Uses `Matrix.det_fromBlocks_zero₂₁` in a NESTED 1+(3+3) decomposition (under `finSumFinEquiv : Fin 1 ⊕ Fin 6 ≃ Fin 7` outer + `finSumFinEquiv : Fin 3 ⊕ Fin 3 ≃ Fin 6` inner). The 1×1 outer block contributes `det = 1` via `det_fin_one`; the two 3×3 inner blocks each contribute `det = -1` via `det_fin_three`; total `det = 1 · (-1) · (-1) = 1`, hence `IsUnit`. **First instance using `det_fromBlocks_zero₂₁`** — orthogonal to the previous nine corner closures (all using `det_of_upperTriangular`). Permutation `ρ ↦ (0, 1, 3, 5, 2, 4, 6), σ ↦ (2, 1, 3, 7, 0, 4, 6)` from S151 pre-search; primes `{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61}` (4 new at S152 — `chiP_thirteen, _forty_one, _fifty_three, _sixty_one`). ~610 Lean lines (49 entry-lemmas + 49+36 = 85 fromBlocks reindex case checks + structural assembly). **Eleventh unconditional `mps_bond_dim` instance; tenth over a wheel `W ≥ 3`.**

**No `axiom` introductions.** All new declarations have `#print axioms` returning only `[propext, Classical.choice, Quot.sound]`.

## What this file formalises

The identity from `novel/mps_bond_dimension.md`:

> For `chiP : [1, W^d] → {0,1}` the prime indicator and `M^(j)` the
> `(W^j × W^(d-j))` unfolding via base-`W` reshape, for every `W ≥ 2`
> and every cut `1 ≤ j < d`:
>
>     rank M^(j)  =  min ( W^j ,  φ(W) · W^(d-j-1) + 1 ).

EDGES.md edge: **E2.1**. Empirically saturated in `novel/mps_bond_dimension.md`
for `W ∈ {2, 6, 30, 210}` and `d` up to 20.

## Lean layout

The proof is decomposed into 8 declarations:

| name                            | role                                           | status   |
|---------------------------------|------------------------------------------------|----------|
| `chiP`                          | prime indicator `ℕ → ℚ`                        | def      |
| `unfolding`                     | the `(W^j × W^(d-j))` matrix over `ℚ`          | def      |
| `rank_le_min_dim`               | trivial ceiling `rank ≤ W^j`                   | **done** |
| `row_support_coprime`           | nonzero entries imply `gcd(k+1, W) = 1`        | **done** |
| `live_columns_count`            | CRT count: live columns = `φ(W) · W^(d-j-1)`   | **done** |
| `upper_bound`                   | `rank ≤ φ(W) · W^(d-j-1) + 1`                  | **done** |
| `chiP_two_eq_one` *(S90)*       | `chiP 2 = 1` (uses `Nat.prime_two`)            | **done** |
| `entry_zero_one_eq_one` *(S90)* | unfolding entry at (0,1) is 1                  | **done** |
| `one_le_rank_unfolding` *(S90)* | trivial floor `1 ≤ rank`                       | **done** |
| `exists_invertible_submatrix`   | exhibit ∃ `R × R` IsUnit submatrix             | `sorry`  |
| `lower_bound`                   | `min(W^j, φ(W)·W^(d-j-1)+1) ≤ rank`            | **done** |
| `mps_bond_dim`                  | **main theorem** (term application)            | **done** |
| `exists_invertible_submatrix_W_eq_2_j_eq_1` *(S98)* | corner case Bertrand exhibit | **done** |
| `mps_bond_dim_W_eq_2_j_eq_1` *(S98)* | corner case `rank = 2` (sorry-free)       | **done** |
| `chiP_three_eq_one` *(S99)*     | `chiP 3 = 1`                                   | **done** |
| `exists_invertible_submatrix_W_eq_2_d_eq_j_plus_1` *(S99)* | second corner exhibit (no Bertrand) | **done** |
| `mps_bond_dim_W_eq_2_d_eq_j_plus_1` *(S99)* | second corner case `rank = 2` (sorry-free) | **done** |
| `chiP_five_eq_one` *(S106)*    | `chiP 5 = 1`                                   | **done** |
| `chiP_seven_eq_one` *(S106)*   | `chiP 7 = 1`                                   | **done** |
| `exists_invertible_submatrix_W_eq_3_d_eq_j_plus_1` *(S106)* | W=3 corner exhibit (no Bertrand) | **done** |
| `mps_bond_dim_W_eq_3_d_eq_j_plus_1` *(S106)* | W=3 corner case `rank = 3` (sorry-free) | **done** |
| `chiP_eleven_eq_one` *(S107)*  | `chiP 11 = 1`                                  | **done** |
| `exists_invertible_submatrix_W_eq_4_d_eq_j_plus_1` *(S107)* | W=4 corner exhibit (no Bertrand) | **done** |
| `mps_bond_dim_W_eq_4_d_eq_j_plus_1` *(S107)* | W=4 corner case `rank = 3` (sorry-free) | **done** |
| `chiP_nineteen_eq_one` *(S117)*       | `chiP 19 = 1`                            | **done** |
| `chiP_twenty_three_eq_one` *(S117)*   | `chiP 23 = 1`                            | **done** |
| `exists_invertible_submatrix_W_eq_5_d_eq_j_plus_1` *(S117)* | W=5 corner exhibit (no Bertrand, BlockTriangular route) | **done** |
| `mps_bond_dim_W_eq_5_d_eq_j_plus_1` *(S117)* | W=5 corner case `rank = 5` (sorry-free) | **done** |
| `chiP_thirty_one_eq_one` *(S122)*     | `chiP 31 = 1`                            | **done** |
| `exists_invertible_submatrix_W_eq_6_d_eq_j_plus_1` *(S122)* | W=6 corner exhibit (no Bertrand, det_fin_three with row choice {0,1,5}) | **done** |
| `mps_bond_dim_W_eq_6_d_eq_j_plus_1` *(S122)* | W=6 corner case `rank = 3` (sorry-free) | **done** |
| `chiP_seventeen_eq_one` *(S128)*       | `chiP 17 = 1`                            | **done** |
| `chiP_thirty_seven_eq_one` *(S128)*    | `chiP 37 = 1`                            | **done** |
| `exists_invertible_submatrix_W_eq_8_d_eq_j_plus_1` *(S128)* | W=8 corner exhibit (no Bertrand, BlockTriangular route) | **done** |
| `mps_bond_dim_W_eq_8_d_eq_j_plus_1` *(S128)* | W=8 corner case `rank = 5` (sorry-free) | **done** |
| `chiP_fifty_nine_eq_one` *(S129)*           | `chiP 59 = 1`                            | **done** |
| `chiP_eighty_nine_eq_one` *(S129)*          | `chiP 89 = 1`                            | **done** |
| `chiP_one_hundred_nine_eq_one` *(S129)*     | `chiP 109 = 1`                           | **done** |
| `chiP_one_hundred_twenty_seven_eq_one` *(S129)* | `chiP 127 = 1`                       | **done** |
| `exists_invertible_submatrix_W_eq_12_d_eq_j_plus_1` *(S129)* | W=12 corner exhibit (4 non-leading rows + BlockTriangular route) | **done** |
| `mps_bond_dim_W_eq_12_d_eq_j_plus_1` *(S129)* | W=12 corner case `rank = 5` (sorry-free) | **done** |
| `chiP_twenty_nine_eq_one` *(S137)*               | `chiP 29 = 1`                            | **done** |
| `chiP_forty_three_eq_one` *(S137)*               | `chiP 43 = 1`                            | **done** |
| `chiP_one_hundred_seventy_nine_eq_one` *(S137)*  | `chiP 179 = 1`                           | **done** |
| `chiP_two_hundred_eleven_eq_one` *(S137)*        | `chiP 211 = 1`                           | **done** |
| `chiP_two_hundred_ninety_three_eq_one` *(S137)*  | `chiP 293 = 1`                           | **done** |
| `exists_invertible_submatrix_W_eq_18_d_eq_j_plus_1` *(S137)* | W=18 corner exhibit (R=7, BlockTriangular route, mixed leading/non-leading rows) | **done** |
| `mps_bond_dim_W_eq_18_d_eq_j_plus_1` *(S137)* | W=18 corner case `rank = 7` (sorry-free) | **done** |
| `chiP_forty_seven_eq_one` *(S143)*                   | `chiP 47 = 1`                            | **done** |
| `chiP_one_hundred_forty_nine_eq_one` *(S143)*        | `chiP 149 = 1`                           | **done** |
| `chiP_one_hundred_ninety_nine_eq_one` *(S143)*       | `chiP 199 = 1`                           | **done** |
| `chiP_two_hundred_forty_one_eq_one` *(S143)*         | `chiP 241 = 1`                           | **done** |
| `chiP_three_hundred_thirty_seven_eq_one` *(S143)*    | `chiP 337 = 1`                           | **done** |
| `prod_univ_nine'` *(S143, private)*                  | local `Fin 9` product expansion          | **done** |
| `exists_invertible_submatrix_W_eq_20_d_eq_j_plus_1` *(S143)* | W=20 corner exhibit (R=9, BlockTriangular route, mixed leading/non-leading rows; `set_option maxHeartbeats 2000000`) | **done** |
| `mps_bond_dim_W_eq_20_d_eq_j_plus_1` *(S143)* | W=20 corner case `rank = 9` (sorry-free) | **done** |
| `chiP_ninety_seven_eq_one` *(S144)*                  | `chiP 97 = 1`                            | **done** |
| `exists_invertible_submatrix_W_eq_10_d_eq_j_plus_1` *(S144)* | W=10 corner exhibit (R=5, BlockTriangular route, ρ ↦ (1, 0, 4, 3, 9)) | **done** |
| `mps_bond_dim_W_eq_10_d_eq_j_plus_1` *(S144)* | W=10 corner case `rank = 5` (sorry-free; refutes S128/S129 obstruction claim) | **done** |
| `chiP_thirteen_eq_one` *(S152)*           | `chiP 13 = 1`                                  | **done** |
| `chiP_forty_one_eq_one` *(S152)*          | `chiP 41 = 1`                                  | **done** |
| `chiP_fifty_three_eq_one` *(S152)*        | `chiP 53 = 1`                                  | **done** |
| `chiP_sixty_one_eq_one` *(S152)*          | `chiP 61 = 1`                                  | **done** |
| `exists_invertible_submatrix_W_eq_9_d_eq_j_plus_1` *(S152)* | W=9 corner exhibit (R=7, **first `det_fromBlocks_zero₂₁` route** — `(1+6)→(3+3)` nested block-DIAGONAL decomposition; permutation from S151 pre-search) | **done** |
| `mps_bond_dim_W_eq_9_d_eq_j_plus_1` *(S152)* | W=9 corner case `rank = 7` (sorry-free; closes the first S128/S129/S144 "block-triangular-required" wheel) | **done** |

(Note: `mps_bond_dim` itself contains no `sorry`; it is a 3-line
`Nat.le_antisymm` of the auxiliary lemmas. As of S83, `lower_bound` is
also `sorry`-free: its proof is a 6-line structural reduction citing
`exists_invertible_submatrix`. The single remaining `sorry` lives
entirely inside `exists_invertible_submatrix` and captures the
prime-density content of the informal proof. Once that exhibit is in
hand, `lower_bound` and `mps_bond_dim` close automatically.)

The five completed proofs (plus the term-mode main theorem):

1. `rank_le_min_dim` is one line: a direct citation of mathlib's
   `Matrix.rank_le_height`.
2. `row_support_coprime` is a 30-line number-theoretic argument: nonzero
   entry ⇒ `n = i·W^(d-j) + k + 1` is prime; `i ≥ 1` and `j < d` give
   `n > W`; hence no prime factor of `W` divides `n`, so `gcd(n, W) = 1`;
   reducing mod `W` (using `gcd_add_mul_right_left` after rewriting
   `W^(d-j) = W^(d-j-1) · W`) yields `gcd(k+1, W) = 1`.
3. `live_columns_count` is a ~110-line two-stage argument (S75):
   - **Stage A (Fin → range):** `Finset.card_bij` with the value-projection
     `k ↦ k.val` shows the count over `Fin (W^(d-j))` equals the count over
     `Finset.range (W^(d-j))`.
   - **Stage B (multi-block induction):** for every `M : ℕ`,
     `|{n ∈ range(W·M) : gcd(n+1,W)=1}| = M · φ(W)` by induction on `M`.
     The successor step splits `range(W·(M+1)) = range(W·M) ∪ Ico(W·M)(W·M+W)`
     (disjoint), invokes the IH on the first piece, and reduces the second
     piece to `Nat.filter_coprime_Ico_eq_totient W 1` via the bijection
     `n ↔ n + 1 − W·M` (`Finset.card_bij'`). The bijection's gcd-preservation
     uses `Nat.gcd_add_mul_right_left` (the same identity that
     `row_support_coprime` uses).
   - Combine with `W^(d-j) = W · W^(d-j-1)` (single-occurrence `conv_lhs`
     rewrite to avoid the subterm-clash with `d-j-1`) and instantiate at
     `M := W^(d-j-1)`.
4. **`upper_bound` is a ~80-line column-span argument (S76).** Strategy:
   - Define `i₀ : Fin(W^j) := ⟨0, hWj_pos⟩` (the row-0 index) and
     `e0 : Fin(W^j) → ℚ := Pi.single i₀ 1` (the row-0 unit vector).
   - Define `GoodCols : Finset (Fin(W^(d-j)))` as the live columns —
     `|GoodCols| = φ(W)·W^(d-j-1)` by `live_columns_count`.
   - Define the generating Finset `S := insert e0 (GoodCols.image col)`.
     `|S| ≤ |GoodCols| + 1 = φ(W)·W^(d-j-1) + 1`.
   - Show every column lies in `Submodule.span ℚ (S : Set _)`:
     - **Good columns** (`gcd(k+1, W) = 1`): in `S` directly.
     - **Bad columns** (`gcd(k+1, W) ≠ 1`): all entries at rows `i ≥ 1`
       vanish (by `row_support_coprime`), so the column equals
       `(unfolding i₀ k) • e0`, hence in `Submodule.span ℚ {e0} ⊆ span S`.
   - Apply `Matrix.rank_eq_finrank_span_cols`, `Submodule.span_le`,
     `Submodule.finrank_mono`, then `finrank_span_finset_le_card S` to
     conclude `rank ≤ |S| ≤ φ(W)·W^(d-j-1) + 1`.
5. **`mps_bond_dim` (main theorem) is a 3-line term-mode proof (S76).**
   `Nat.le_antisymm` of:
   - `Nat.le_min.mpr ⟨rank_le_min_dim, upper_bound⟩` for the upper side
     (`rank ≤ min(W^j, φ(W)·W^(d-j-1) + 1)`).
   - `lower_bound` for the lower side.
   The main theorem itself contains no `sorry`. It currently still has
   an open obligation only because `lower_bound` is not yet closed —
   once `lower_bound` loses its `sorry`, `mps_bond_dim` is automatically
   closed without modification.

## What the next session needs to do

Only one open obligation remains:

**`exists_invertible_submatrix`** — the prime-density exhibit. State:
```
∃ (ρ : Fin R → Fin (W^j)) (σ : Fin R → Fin (W^(d-j))),
   IsUnit ((unfolding W d j).submatrix ρ σ)
```
where `R = min(W^j, φ(W)·W^(d-j-1) + 1)`.

### S90 update — Route A is structurally insufficient with current mathlib

After auditing mathlib's `NumberTheory.Bertrand` and `NumberTheory.PrimeCounting`,
**Route A as originally outlined cannot be completed without first formalising
Hoheisel-type density results that mathlib does not yet have.**

The argument: for the `i`-th row of the unfolding, we need a prime in
the half-open interval `(i · W^(d-j), (i+1) · W^(d-j)]`. The interval
length is `W^(d-j)`, the endpoint is `(i+1) · W^(d-j)`. Bertrand's
postulate (`Nat.bertrand` / `Nat.exists_prime_lt_and_le_two_mul`)
provides a prime in `(n, 2n]`, i.e., **interval-length equal to endpoint**.
For our row `i`, interval-length `W^(d-j)` equals endpoint
`(i+1) · W^(d-j)` only when `i = 0`. For `i ≥ 1`, we need primes in
intervals of length `1/(i+1)` of the endpoint — strictly shorter than
Bertrand provides.

In the worst case (`R = φ(W)·W^(d-j-1) + 1`, the typical "large `j`"
regime), `i` can be as large as `φ(W)·W^(d-j-1)`, so the required
interval-length-to-endpoint ratio is `1/φ(W) · W^(-(d-j-1))` —
polynomially small in `W^(d-j-1)`. This is a Hoheisel-type
("primes in `[n, n + n^θ]` for some `θ < 1`") density question, far
beyond mathlib's current toolbox.

### What S90 closed unconditionally

* **`chiP_two_eq_one : chiP 2 = 1`** (one line, uses `Nat.prime_two`).
* **`entry_zero_one_eq_one`**: matrix entry at row 0, column 1 equals 1.
* **`one_le_rank_unfolding : 1 ≤ (unfolding W d j).rank`** — the floor
  case `R = 1` of the lower bound. Witness: the `1 × 1` submatrix at
  row 0, column 1 is `chiP 2 = 1 ≠ 0`, hence `IsUnit`. Combined with
  `Matrix.rank_of_isUnit` and `Matrix.rank_submatrix_le`. ~25 lines.
  No prime-density beyond `Nat.prime_two`.

This is the only piece of the lower bound that is unconditional. It
closes nothing new in `mps_bond_dim` itself (since `R ≥ 2` whenever
`W ≥ 2`, `1 ≤ j`, `j < d`), but it isolates the trivial portion of
the prime-density argument from the deep portion.

### Routes for closing `exists_invertible_submatrix` (revised)

* **Route A (Hoheisel-grade).** Formalise primes-in-short-intervals at
  the polynomial scale. Beyond a single session — likely a multi-session
  arc by itself, possibly more than this entire L1 formalisation has
  been so far.

* **Route A' (Bertrand-only sub-cases). DONE S98.** Closed
  `exists_invertible_submatrix` for the corner case `(W = 2, j = 1)`
  via two new theorems: `exists_invertible_submatrix_W_eq_2_j_eq_1`
  (the prime exhibit) and `mps_bond_dim_W_eq_2_j_eq_1 : (unfolding 2 d 1).rank = 2`
  (the unconditional rank statement). Both are sorry-free. The
  construction picks σ as `(column 1, column p − 2^(d−1) − 1)` for the
  Bertrand prime `p ∈ (2^(d−1), 2 · 2^(d−1)]`; the resulting `2 × 2`
  submatrix is upper-triangular `[[1, ?], [0, 1]]` with `det = 1`,
  because `2^(d−1) + 2` is even and `> 2`, hence not prime. ~70 Lean
  lines. The general case (the original `exists_invertible_submatrix`)
  remains the only `sorry` in the file.

* **Route A'' (orthogonal corner, even simpler than A'). DONE S99.**
  Closed `exists_invertible_submatrix` for the orthogonal corner case
  `(W = 2, d = j + 1)` (i.e. `d - j = 1`) via two new theorems:
  `exists_invertible_submatrix_W_eq_2_d_eq_j_plus_1` (the prime
  exhibit) and `mps_bond_dim_W_eq_2_d_eq_j_plus_1 : (unfolding 2 (j+1) j).rank = 2`
  for every `j ≥ 1`. Both sorry-free. **No Bertrand required**: the
  matrix has only `2^(d-j) = 2^1 = 2` columns, so we take *both*; rows
  `{0, 1}` (available since `2^j ≥ 2` for `j ≥ 1`) plus the column-swap
  `σ = (col 1, col 0)` give the `2 × 2` submatrix
  ```
     ⎡ unfolding(0, 1),  unfolding(0, 0) ⎤   ⎡ chiP 2, chiP 1 ⎤   ⎡ 1, 0 ⎤
     ⎣ unfolding(1, 1),  unfolding(1, 0) ⎦ = ⎣ chiP 4, chiP 3 ⎦ = ⎣ 0, 1 ⎦
  ```
  of determinant `1`. Only `Nat.prime_two`, `Nat.prime_three`,
  `Nat.not_prime_one`, and decidability of `¬ Nat.Prime 4` are used.
  ~110 Lean lines (counting docstrings + the helper `chiP_three_eq_one`).
  Together with Route A', the file now has unconditional Lean proofs of
  `mps_bond_dim` whenever **either** `j = 1` **or** `d - j = 1` —
  i.e. on the entire boundary of the `(j, d − j)` integer grid. The
  case `(j = 1, d = 2)` is the unique overlap; `j ≥ 2` (with `d = j + 1`)
  is genuinely new content.

* **Route A''' (orthogonal corner, `W = 3`). DONE S106.**
  Closed `exists_invertible_submatrix` for `(W = 3, d = j + 1)` via
  four new declarations: `chiP_five_eq_one`, `chiP_seven_eq_one`,
  `exists_invertible_submatrix_W_eq_3_d_eq_j_plus_1`, and
  `mps_bond_dim_W_eq_3_d_eq_j_plus_1 : (unfolding 3 (j+1) j).rank = 3`
  for every `j ≥ 1`. All sorry-free; `#print axioms` confirms only
  `propext, Classical.choice, Quot.sound`. The matrix has only
  `3^(d-j) = 3^1 = 3` columns, so we take all three; rows `{0, 1, 2}`
  (available since `3^j ≥ 3` for `j ≥ 1`) and `σ = id` give the `3 × 3`
  submatrix
  ```
     ⎡ chiP 1, chiP 2, chiP 3 ⎤   ⎡ 0, 1, 1 ⎤
     ⎢ chiP 4, chiP 5, chiP 6 ⎥ = ⎢ 0, 1, 0 ⎥
     ⎣ chiP 7, chiP 8, chiP 9 ⎦   ⎣ 1, 0, 0 ⎦
  ```
  of determinant `−1` (computed via `Matrix.det_fin_three`). The unit
  witness is `isUnit_one.neg : IsUnit (-(1 : ℚ))`, simplified to
  `IsUnit (-1)` via `ring_nf`. Uses only `Nat.prime_two`,
  `Nat.prime_three`, `Nat.prime_five`, `Nat.prime_seven` (all
  `decide`-checkable) and decidability of `¬ Nat.Prime` for `1, 4, 6,
  8, 9`. ~150 Lean lines. **First unconditional `mps_bond_dim`
  instance over a wheel `W ≥ 3`.** Pattern is reusable: any `W` with
  an explicit invertible `(R(W) × R(W))` chiP-submatrix in the first
  `R(W)` rows of the `d − j = 1` slab admits the same closure.

* **Route A''''' (orthogonal corner, `W = 5`). DONE S117.**
  Closed `exists_invertible_submatrix` for `(W = 5, d = j + 1)` via four
  new declarations: `chiP_nineteen_eq_one`, `chiP_twenty_three_eq_one`,
  `exists_invertible_submatrix_W_eq_5_d_eq_j_plus_1`, and
  `mps_bond_dim_W_eq_5_d_eq_j_plus_1 : (unfolding 5 (j+1) j).rank = 5`
  for every `j ≥ 1`. All sorry-free; `#print axioms` confirms only
  `propext, Classical.choice, Quot.sound`. The matrix is `5^j × 5`; we
  pick **all five columns** (the dead column at `k = 4` contributes the
  diagonal entry `chiP 5 = 1` which bumps rank from `4 = φ(5)` to `5`)
  and rows `{0, 1, 2, 3, 4}` (available since `5^j ≥ 5` for `j ≥ 1`).

  **First non-`det_fin_three` proof.** Mathlib has `Matrix.det_fin_three`
  but no `det_fin_four` / `det_fin_five`. Instead, `ρ` and `σ` are chosen
  as *permutations* such that the `5 × 5` submatrix is upper triangular
  with `1` on the diagonal:
  ```
     ⎡ chiP  5, chiP  4, chiP  1, chiP  2, chiP  3 ⎤   ⎡ 1, 0, 0, 1, 1 ⎤
     ⎢ chiP 20, chiP 19, chiP 16, chiP 17, chiP 18 ⎥   ⎢ 0, 1, 0, 1, 0 ⎥
     ⎢ chiP 15, chiP 14, chiP 11, chiP 12, chiP 13 ⎥ = ⎢ 0, 0, 1, 0, 1 ⎥.
     ⎢ chiP 10, chiP  9, chiP  6, chiP  7, chiP  8 ⎥   ⎢ 0, 0, 0, 1, 0 ⎥
     ⎣ chiP 25, chiP 24, chiP 21, chiP 22, chiP 23 ⎦   ⎣ 0, 0, 0, 0, 1 ⎦
  ```
  Concretely `ρ : Fin 5 → Fin (5^j)` permutes original rows to
  `(0, 3, 2, 1, 4)` and `σ : Fin 5 → Fin (5^((j+1)-j))` permutes original
  columns to `(4, 3, 0, 1, 2)`. The diagonal hits the five primes
  `{5, 19, 11, 7, 23}` and the lower triangle hits the ten composites
  `{20, 15, 14, 10, 9, 6, 25, 24, 21, 22}` — all in `[1, 25]`.

  **Determinant via `det_of_upperTriangular`.** The proof:
  - Pre-compute the 5 diagonal entries (each = 1) and the 10 below-
    diagonal entries (each = 0) via 15 `change → rw [h_sub] → norm_num →
    chiP_..._eq_one / simp [chiP, ¬ Nat.Prime _]` steps.
  - Establish `(submatrix Mρ Mσ).BlockTriangular id` by `intro i k h_lt`
    and `fin_cases i <;> fin_cases k`. The 15 vacuous cases (where
    `id k < id i` is false) close via
    `simp [id_eq, Fin.lt_def] at h_lt` then
    `exact absurd h_lt (by decide)`. The 10 below-diagonal cases reduce
    to `exact hLij`.
  - Apply `Matrix.det_of_upperTriangular` to rewrite the determinant as
    `∏ i, M i i`, then `Fin.prod_univ_five` to expand to a 5-term product.
    Substitute the 5 diagonal entries and finish with `norm_num`.

  Uses `Nat.prime_two, prime_three, prime_five, prime_seven, prime_eleven,
  prime_nineteen, prime_twenty_three` (via the `chiP` helpers) and
  decidability of non-primality for `1, 4, 6, 9, 10, 14, 15, 20, 21, 22,
  24, 25`. ~250 Lean lines. **Third unconditional `mps_bond_dim` instance
  over a wheel `W ≥ 3`; the first instance with `R = W`** (so all `W`
  columns must be retained); **the first instance using
  `det_of_upperTriangular`** rather than a `det_fin_n` formula. The
  technique scales to every wheel `W` for which we can exhibit a
  permutation of `Fin W` that triangularises the slab `chiP 1 .. chiP W^2`
  — the existence of such a permutation is a finite check at each `W`,
  not a number-theoretic fact.

  **Why no `upper_bound` here.** Unlike W=4 (where `rank_le_width`
  gives only `rank ≤ 4`, not the sharp `≤ 3`), at W=5 the column count
  `5^((j+1)-j) = 5` exactly equals `R = min(5^j, 5) = 5`, so
  `rank_le_width` directly yields the sharp upper bound.

* **Route A'''''' (orthogonal corner, `W = 6`). DONE S122.**
  Closed `exists_invertible_submatrix` for `(W = 6, d = j + 1)` via three
  new declarations: `chiP_thirty_one_eq_one`,
  `exists_invertible_submatrix_W_eq_6_d_eq_j_plus_1`, and
  `mps_bond_dim_W_eq_6_d_eq_j_plus_1 : (unfolding 6 (j+1) j).rank = 3`
  for every `j ≥ 1`. All sorry-free; `#print axioms` confirms only
  `propext, Classical.choice, Quot.sound`. The matrix has `6^(d-j) =
  6` columns; live columns are `{0, 4}` (residues `1, 5 (mod 6)`), so
  `R = φ(6) + 1 = 3`. We pick the two live columns plus dead column
  `1` (which contributes `chiP 2 = 1` at row `0`).

  **Row choice subtlety (the key novelty over W ∈ {2, 3, 4, 5}).** The
  first three rows of the `6 × 6` slab are *not* linearly independent:
  rows `1, 2, 3` (corresponding to `chiP 7..12, 13..18, 19..24`) all
  have the identical support pattern `(1, 0, 0, 0, 1, 0)`, since each
  of those windows contains exactly two primes at residues `1, 5 (mod 6)`.
  We therefore pick rows `{0, 1, 5}` — the row-`5` window `chiP 31..36`
  is `(1, 0, 0, 0, 0, 0)` because `31` is the only prime in that window
  and it sits at residue `1 (mod 6)` (column `0`). With `ρ ↦ (0, 1, 5)`
  and `σ ↦ (0, 1, 4)` the `3 × 3` submatrix is
  ```
     ⎡ chiP  1, chiP  2, chiP  5 ⎤   ⎡ 0, 1, 1 ⎤
     ⎢ chiP  7, chiP  8, chiP 11 ⎥ = ⎢ 1, 0, 1 ⎥
     ⎣ chiP 31, chiP 32, chiP 35 ⎦   ⎣ 1, 0, 0 ⎦
  ```
  of determinant `+1` (computed via `Matrix.det_fin_three`); the unit
  witness is `isUnit_one`. **First orthogonal-corner instance where
  the working row set is not `{0, 1, ..., R-1}`** — every prior corner
  (W ∈ {2, 3, 4, 5}) used the leading rows.

  **Upper-bound subtlety:** `rank_le_width` gives only `rank ≤ 6`, so
  we cite the general `upper_bound`, evaluating to `φ(6) · 6^0 + 1 =
  3` (same pattern as W=4).

  Uses `Nat.prime_two, prime_three, prime_five, prime_seven,
  prime_eleven, prime_thirty_one` (last new at S122, all
  `decide`-checkable) and decidability of non-primality for `1, 8,
  32, 35`. ~190 Lean lines. **Sixth unconditional `mps_bond_dim`
  instance; fourth instance over a wheel `W ≥ 3`.** Sets the template
  for higher-`W` corners where the first `R` rows are LD: identify
  the row pattern of the `W × W` slab and pick a row-set spanning all
  three live-column-residues plus a single dead column.

* **Route A''''''' (orthogonal corner, `W = 8`). DONE S128.**
  Closed `exists_invertible_submatrix` for `(W = 8, d = j + 1)` via
  four new declarations: `chiP_seventeen_eq_one`,
  `chiP_thirty_seven_eq_one`,
  `exists_invertible_submatrix_W_eq_8_d_eq_j_plus_1`, and
  `mps_bond_dim_W_eq_8_d_eq_j_plus_1 : (unfolding 8 (j+1) j).rank = 5`
  for every `j ≥ 1`. All sorry-free; `#print axioms` confirms only
  `propext, Classical.choice, Quot.sound`. The matrix is `8^j × 8`;
  live columns are `{0, 2, 4, 6}` (residues `1, 3, 5, 7 (mod 8)`,
  i.e. `k + 1` odd), so `R = φ(8) + 1 = 5`. We pick the four live
  columns plus dead column `1` (the unique dead column whose `chiP`
  at row 0 is non-zero, since `chiP 2 = 1`).

  **Triangulation via permutation (BlockTriangular route, à la W=5).**
  The first 5 rows `{0, 1, 2, 3, 4}` of the `8 × 8` slab restricted to
  cols `{0, 1, 2, 4, 6}` are linearly independent. The permutation
  `ρ ↦ (2, 0, 1, 3, 4)` and `σ ↦ (0, 1, 2, 6, 4)` triangularises the
  `5 × 5` submatrix:
  ```
     ⎡ chiP 17, chiP 18, chiP 19, chiP 23, chiP 21 ⎤   ⎡ 1, 0, 1, 1, 0 ⎤
     ⎢ chiP  1, chiP  2, chiP  3, chiP  7, chiP  5 ⎥   ⎢ 0, 1, 1, 1, 1 ⎥
     ⎢ chiP  9, chiP 10, chiP 11, chiP 15, chiP 13 ⎥ = ⎢ 0, 0, 1, 0, 1 ⎥.
     ⎢ chiP 25, chiP 26, chiP 27, chiP 31, chiP 29 ⎥   ⎢ 0, 0, 0, 1, 1 ⎥
     ⎣ chiP 33, chiP 34, chiP 35, chiP 39, chiP 37 ⎦   ⎣ 0, 0, 0, 0, 1 ⎦
  ```
  with diagonal primes `{17, 2, 11, 31, 37}`.

  **Determinant via `det_of_upperTriangular`.** Same skeleton as W=5:
  pre-compute the 5 diagonal entries (each `= 1`) and the 10 below-
  diagonal entries (each `= 0`); establish `BlockTriangular id` by
  `fin_cases i <;> fin_cases k`; apply `det_of_upperTriangular`,
  `Fin.prod_univ_five`, then `norm_num`.

  **Upper-bound subtlety:** `rank_le_width` gives only `rank ≤ 8`,
  not the sharp `rank ≤ 5 = φ(8) + 1`. We cite the general
  `upper_bound`, which evaluates to `φ(8) · 8^0 + 1 = 4 · 1 + 1 = 5`
  in this corner (same pattern as W=4 and W=6).

  Uses `Nat.prime_two, prime_eleven, prime_thirty_one` (existing) and
  `Nat.prime_seventeen, prime_thirty_seven` (new at S128, all
  `decide`-checkable), plus decidability of non-primality for
  `1, 9, 10, 25, 26, 27, 33, 34, 35, 39`. ~280 Lean lines.
  **Seventh unconditional `mps_bond_dim` instance; fifth instance
  over a wheel `W ≥ 3`; second instance using `det_of_upperTriangular`
  (after W=5).** Confirms that the BlockTriangular pattern from W=5
  scales to wheels with composite `W` and `R < W` (`R = 5` here, vs.
  `R = W = 5` at W=5).

* **Route A^{(8)} (orthogonal corner, `W = 12`). DONE S129.**
  Closed `exists_invertible_submatrix` for `(W = 12, d = j + 1)` via
  six new declarations: `chiP_fifty_nine_eq_one`,
  `chiP_eighty_nine_eq_one`, `chiP_one_hundred_nine_eq_one`,
  `chiP_one_hundred_twenty_seven_eq_one`,
  `exists_invertible_submatrix_W_eq_12_d_eq_j_plus_1`, and
  `mps_bond_dim_W_eq_12_d_eq_j_plus_1 : (unfolding 12 (j+1) j).rank = 5`
  for every `j ≥ 1`. All sorry-free; `#print axioms` confirms only
  `propext, Classical.choice, Quot.sound`. The matrix is `12^j × 12`;
  live columns are `{0, 4, 6, 10}` (residues `1, 5, 7, 11 (mod 12)`),
  giving `φ(12) = 4` live columns. We need `R = φ(12) + 1 = 5`
  columns total, so we add one dead column with `chiP` non-zero at
  row 0 — two candidates exist (`k = 1` with `chiP 2 = 1` and `k = 2`
  with `chiP 3 = 1`); we pick `k = 1`.

  **Skipping `W ∈ {7, 9, 10, 11}` is structural, not bookkeeping.** All
  four W's hit the same obstruction: NO `R × R` sub-matrix of the
  `W × W` slab admits a triangulation with rows `< W` (which is the
  binding constraint at `j = 1`). For W=7, the 7-prime live-col matrix
  has all column-counts `≥ 2` (no greedy single-1 col exists). For
  W=9, after fixing `(r0, c2)` as the unique pairing for the dead col,
  the remaining `6 × 6` sub-matrix decomposes into two **block-
  diagonal** `3 × 3` blocks (left rows/cols disjoint from right),
  neither of which admits standalone upper-triangulation; closing W=9
  needs `Matrix.det_of_blockTriangular` (multi-session new technique).
  W=10 has all candidate row sets generating only rank 2 (the
  multiplicity-2 residue pattern). W=11 is prime so behaves like W=7.

  **Triangulation via NON-LEADING rows (W=12).** The first 5 rows
  `{0, 1, 2, 3, 4}` of the `12 × 12` slab give linearly dependent
  windows (row 1 and row 3 have identical support pattern at the chosen
  cols). Non-leading rows `{4, 7, 9, 10}` each contribute a distinct
  prime at a distinct column, breaking the linear-dependence. The
  permutation `ρ ↦ (0, 9, 10, 4, 7)` and `σ ↦ (1, 0, 6, 10, 4)`
  triangularises the `5 × 5` submatrix:
  ```
     ⎡ chiP   2, chiP   1, chiP   7, chiP  11, chiP   5 ⎤   ⎡ 1, 0, 1, 1, 1 ⎤
     ⎢ chiP 110, chiP 109, chiP 115, chiP 119, chiP 113 ⎥   ⎢ 0, 1, 0, 0, 1 ⎥
     ⎢ chiP 122, chiP 121, chiP 127, chiP 131, chiP 125 ⎥ = ⎢ 0, 0, 1, 1, 0 ⎥.
     ⎢ chiP  50, chiP  49, chiP  55, chiP  59, chiP  53 ⎥   ⎢ 0, 0, 0, 1, 1 ⎥
     ⎣ chiP  86, chiP  85, chiP  91, chiP  95, chiP  89 ⎦   ⎣ 0, 0, 0, 0, 1 ⎦
  ```
  with diagonal primes `{2, 109, 127, 59, 89}` and below-diagonal
  composites `{49, 50, 55, 85, 86, 91, 95, 110, 121, 122}`.

  **Determinant via `det_of_upperTriangular`.** Same skeleton as W=5
  and W=8: pre-compute the 5 diagonal entries (each `= 1`) and the 10
  below-diagonal entries (each `= 0`); establish `BlockTriangular id`
  by `fin_cases i <;> fin_cases k`; apply `det_of_upperTriangular`,
  `Fin.prod_univ_five`, then `norm_num`.

  **Upper-bound subtlety:** `rank_le_width` gives only `rank ≤ 12`,
  not the sharp `rank ≤ 5 = φ(12) + 1`. We cite the general
  `upper_bound`, which evaluates to `φ(12) · 12^0 + 1 = 4 · 1 + 1 = 5`
  (same pattern as W=4, W=6, W=8).

  Uses `Nat.prime_two` (existing) and four new `chiP_X_eq_one` helpers
  for `X ∈ {59, 89, 109, 127}` (all `decide`-checkable), plus
  decidability of non-primality for
  `49, 50, 55, 85, 86, 91, 95, 110, 121, 122`. ~370 Lean lines.
  **Eighth unconditional `mps_bond_dim` instance; sixth instance
  over a wheel `W ≥ 3`; third instance using `det_of_upperTriangular`
  (after W=5 and W=8).** **First instance using FOUR non-leading rows**
  — extends S122's W=6 single-non-leading-row trick (rows `{0, 1, 5}`)
  to the maximally non-leading regime where only row `0` is leading.

* **Route A^{(9)} (orthogonal corner, `W = 18`). DONE S137.**
  Closed `exists_invertible_submatrix` for `(W = 18, d = j + 1)` via
  six new declarations: `chiP_twenty_nine_eq_one`,
  `chiP_forty_three_eq_one`, `chiP_one_hundred_seventy_nine_eq_one`,
  `chiP_two_hundred_eleven_eq_one`, `chiP_two_hundred_ninety_three_eq_one`,
  `exists_invertible_submatrix_W_eq_18_d_eq_j_plus_1`, and
  `mps_bond_dim_W_eq_18_d_eq_j_plus_1 : (unfolding 18 (j+1) j).rank = 7`
  for every `j ≥ 1`. All sorry-free; `#print axioms` confirms only
  `propext, Classical.choice, Quot.sound`. The matrix is `18^j × 18`;
  live columns are `{0, 4, 6, 10, 12, 16}` (residues `1, 5, 7, 11, 13, 17
  (mod 18)`), giving `φ(18) = 6` live columns. We need `R = 7`, so we add
  one dead column with `chiP` non-zero — two candidates exist (`k = 1`
  with `chiP 2 = 1` and `k = 2` with `chiP 3 = 1`); we pick `k = 1`.

  **Triangulation via mixed leading/non-leading rows (W=18).**
  The leading rows `{0, 1, 2, 3, 4, 5, 6}` of the 18×18 slab admit no
  triangulation; rows `2` and `5` have identical support at the 7 chosen
  columns. Picking rows `{0, 1, 2, 6, 9, 11, 16}` works — in particular
  rows `9, 11, 16` contribute primes `179, 211, 293` at distinct columns
  while the leading rows `0, 1, 2, 6` cover the other four diagonal
  primes `{2, 29, 43, 109}`. Permutation
  `ρ ↦ (0, 2, 9, 1, 11, 6, 16)` and `σ ↦ (1, 6, 16, 10, 12, 0, 4)`
  triangularises the `7 × 7` submatrix:
  ```
     ⎡ chiP   2, chiP   7, chiP  17, chiP  11, chiP  13, chiP   1, chiP   5 ⎤   ⎡ 1, 1, 1, 1, 1, 0, 1 ⎤
     ⎢ chiP  38, chiP  43, chiP  53, chiP  47, chiP  49, chiP  37, chiP  41 ⎥   ⎢ 0, 1, 1, 1, 0, 1, 1 ⎥
     ⎢ chiP 164, chiP 169, chiP 179, chiP 173, chiP 175, chiP 163, chiP 167 ⎥   ⎢ 0, 0, 1, 1, 0, 1, 1 ⎥
     ⎢ chiP  20, chiP  25, chiP  35, chiP  29, chiP  31, chiP  19, chiP  23 ⎥ = ⎢ 0, 0, 0, 1, 1, 1, 1 ⎥.
     ⎢ chiP 200, chiP 205, chiP 215, chiP 209, chiP 211, chiP 199, chiP 203 ⎥   ⎢ 0, 0, 0, 0, 1, 1, 0 ⎥
     ⎢ chiP 110, chiP 115, chiP 125, chiP 119, chiP 121, chiP 109, chiP 113 ⎥   ⎢ 0, 0, 0, 0, 0, 1, 1 ⎥
     ⎣ chiP 290, chiP 295, chiP 305, chiP 299, chiP 301, chiP 289, chiP 293 ⎦   ⎣ 0, 0, 0, 0, 0, 0, 1 ⎦
  ```
  with diagonal primes `{2, 43, 179, 29, 211, 109, 293}` and below-
  diagonal composites `{20, 25, 35, 38, 110, 115, 119, 121, 125, 164,
  169, 200, 205, 209, 215, 289, 290, 295, 299, 301, 305}`.

  **Determinant via `det_of_upperTriangular`.** Same skeleton as W=5,
  W=8, W=12: pre-compute the 7 diagonal entries (each `= 1`) and the
  21 below-diagonal entries (each `= 0`); establish `BlockTriangular id`
  by `fin_cases i <;> fin_cases k`; apply `det_of_upperTriangular`,
  `Fin.prod_univ_seven`, then `norm_num`.

  **Upper-bound subtlety:** `rank_le_width` gives only `rank ≤ 18`,
  not the sharp `rank ≤ 7 = φ(18) + 1`. We cite the general
  `upper_bound`, which evaluates to `φ(18) · 18^0 + 1 = 6 · 1 + 1 = 7`
  (same pattern as W=4, W=6, W=8, W=12).

  **Recursion-depth lesson.** The five new prime helpers (29, 43, 179,
  211, 293) and the not-prime helpers for composites ≥ 150 cannot use
  `decide` — they hit Lean's `maxRecDepth`. Switching to `norm_num`
  resolves this; this is the file's first use of `norm_num` for
  primality (existing helpers up to W=12 stayed under the limit at
  `decide`).

  Uses `Nat.prime_two` (existing) and five new `chiP_X_eq_one` helpers
  for `X ∈ {29, 43, 179, 211, 293}` (all `norm_num`-checkable; `decide`
  hits `maxRecDepth` for the largest two), plus `chiP_one_hundred_nine_eq_one`
  (existing from S129). ~520 Lean lines (the largest single-corner block
  in the file). **Ninth unconditional `mps_bond_dim` instance; seventh
  instance over a wheel `W ≥ 3`; fourth instance using
  `det_of_upperTriangular` (after W=5, W=8, W=12).** **First instance
  with `R = 7`.**

  **Why W=14 is structurally obstructed (negative shape edge).**
  Like W ∈ {7, 9, 10, 11}, W=14 cannot be closed by leading-row
  triangulation. Concretely: the 14×14 W=14 j=1 slab has rank 7 (so
  `mps_bond_dim` holds), but exhaustive search over column permutations
  + greedy row selection within `[0, 14)` finds zero upper-triangulations.
  The obstruction: rows 2 and 5 of the W=14 slab have identical support
  pattern `(1, 1, 0, 1, 0, 1)` at cols `{0, 2, 4, 8, 10, 12}` (cf. the
  primes `29, 31, 37, 41` resp. `71, 73, 79, 83`); any 7-row set must
  include only one of them, and the remaining 13-row pool plus 1-row
  redundancy structure does not admit a triangulation. The closure
  needs `Matrix.det_of_blockTriangular` (the same multi-session technique
  W=9 needs).

* **Route A^{(10)} (orthogonal corner, `W = 20`). DONE S143.**
  Closed `exists_invertible_submatrix` for `(W = 20, d = j + 1)` via
  six new declarations: `chiP_forty_seven_eq_one`,
  `chiP_one_hundred_forty_nine_eq_one`,
  `chiP_one_hundred_ninety_nine_eq_one`,
  `chiP_two_hundred_forty_one_eq_one`,
  `chiP_three_hundred_thirty_seven_eq_one`,
  `exists_invertible_submatrix_W_eq_20_d_eq_j_plus_1`, and
  `mps_bond_dim_W_eq_20_d_eq_j_plus_1 : (unfolding 20 (j+1) j).rank = 9`
  for every `j ≥ 1`. All sorry-free; `#print axioms` confirms only
  `propext, Classical.choice, Quot.sound`. The matrix is `20^j × 20`;
  live columns are `{0, 2, 6, 8, 10, 12, 16, 18}` (residues `1, 3, 7, 9,
  11, 13, 17, 19 (mod 20)`), giving `φ(20) = 8` live columns. We need
  `R = 9`, so we add one dead column with `chiP` non-zero — two
  candidates exist (`k = 1` with `chiP 2 = 1` and `k = 4` with
  `chiP 5 = 1`); we pick `k = 1`.

  **Triangulation via mixed leading/non-leading rows (W=20).**
  The leading rows `{0..8}` of the 20×20 slab admit no triangulation
  (verified by Python pre-search). Picking rows
  `{0, 1, 2, 7, 9, 10, 12, 14, 16}` works — non-leading rows
  `{9, 10, 12, 14, 16}` contribute primes `{199, 211, 241, 293, 337}`
  and leading rows `{0, 1, 2, 7}` cover `{2, 23, 47, 149}`. Permutation
  `ρ ↦ (0, 2, 9, 14, 1, 7, 12, 16, 10)` and `σ ↦ (1, 6, 18, 12, 2, 8, 0,
  16, 10)` triangularises the `9 × 9` submatrix to upper triangular
  with diagonal primes `{2, 47, 199, 293, 23, 149, 241, 337, 211}` and
  36 below-diagonal composites `{22, 27, 33, 39, 42, 142, 143, 147,
  153, 159, 182, 187, 201, 202, 203, 207, 209, 213, 217, 219, 242, 243,
  247, 249, 253, 259, 282, 287, 299, 321, 322, 323, 327, 329, 333,
  339}`. Pre-search picked the triangulation that minimises the maximum
  diagonal prime; `max_diag = 337` is unavoidable across all 600
  triangulations (300 per dead-column choice).

  **Determinant via `det_of_upperTriangular`.** Same skeleton as W=5,
  W=8, W=12, W=18: pre-compute the 9 diagonal entries (each `= 1`) and
  the 36 below-diagonal entries (each `= 0`); establish `BlockTriangular
  id` by `fin_cases i <;> fin_cases k` (81 subgoals); apply
  `det_of_upperTriangular`. **New helper: `prod_univ_nine'`.** Mathlib
  provides `Fin.prod_univ_eight` but not `prod_univ_nine`; we add a
  private `prod_univ_nine'` lemma following mathlib's pattern verbatim
  (`rw [Fin.prod_univ_castSucc, Fin.prod_univ_eight]; rfl`).

  **Heartbeat scaling.** The default `maxHeartbeats 200000` is
  insufficient for `R = 9`. The `R²` simp blow-up (81 vs 49 fin_cases
  subgoals) plus the 9-deep if-then-else chain (vs 7-deep at W=18)
  pushes simp past the limit. We scope `set_option maxHeartbeats
  2000000 in` to just the `exists_invertible_submatrix_W_eq_20_d_eq_j_plus_1`
  declaration. This is the file's first heartbeat-bumped declaration
  and a forward-compatibility note: every higher-`R` corner will need
  similar treatment.

  **Upper-bound subtlety:** `rank_le_width` gives only `rank ≤ 20`,
  not the sharp `rank ≤ 9 = φ(20) + 1`. We cite the general
  `upper_bound`, which evaluates to `φ(20) · 20^0 + 1 = 8 · 1 + 1 = 9`.

  Uses `Nat.prime_two`, `chiP_twenty_three_eq_one` (existing from S117),
  `chiP_two_hundred_eleven_eq_one` (existing from S137),
  `chiP_two_hundred_ninety_three_eq_one` (existing from S137), and five
  new helpers for `{47, 149, 199, 241, 337}` (all `norm_num`-checkable).
  ~570 Lean lines (the largest single-corner block in the file, beating
  W=18's 520). **Tenth unconditional `mps_bond_dim` instance; eighth
  instance over a wheel `W ≥ 3`; fifth instance using
  `det_of_upperTriangular`.** **First instance with `R = 9`** and
  **first heartbeat-bumped declaration in the file**.

  **Why W ∈ {15, 16, 24, 30} are structurally obstructed (S143
  pre-search).** All four of these prime-power-free wheels with
  `R = φ(W) + 1 = 9` admit zero leading-row+dead-col upper-
  triangulations within rows `[0, W)`. They join `W ∈ {7, 9, 10, 11,
  14}` in the "`det_of_blockTriangular`-required" set. Concretely:
  - W=15: live = `{0, 1, 3, 6, 7, 10, 12, 13}`, dead candidates =
    `{2, 4}`; both column sets fail.
  - W=16: live = `{0, 2, 4, 6, 8, 10, 12, 14}`, dead = `{1}`; fails.
  - W=24: live = `{0, 4, 6, 10, 12, 16, 18, 22}`, dead = `{1, 2}`; both fail.
  - W=30: live = `{0, 6, 10, 12, 16, 18, 22, 28}`, dead = `{1, 2, 4}`; all three fail.

  The closed-W set is now `{2, 3, 4, 5, 6, 8, 12, 18, 20}`; the
  block-triangular-required set is `{7, 9, 10, 11, 14, 15, 16, 24, 30}`.
  The next single-session leading-row candidates are W ∈ {21, 22, 25,
  26, 27, 28, 33, 35, 36, ...} (untested). At W=21, R = 13 — substantially
  more work; at W=22, R = 11; at W=25 (= 5²), R = 21. Cleanest next
  candidates are W ∈ {21, 22, 28} where R ∈ {11, 12, 13}; ambitious
  candidates are W ∈ {25, 27, 35} where R is in the 19-25 range
  (mathlib has no `Fin.prod_univ_X` for X ≥ 9, so each new R needs a
  local `prod_univ_X` lemma — the W=20 closure adds `prod_univ_nine'`;
  successor closures will add the rest).

* **Route A^{(11)} (orthogonal corner, `W = 10`). DONE S144.**
  Closed `exists_invertible_submatrix` for `(W = 10, d = j + 1)` via
  three new declarations: `chiP_ninety_seven_eq_one`,
  `exists_invertible_submatrix_W_eq_10_d_eq_j_plus_1`, and
  `mps_bond_dim_W_eq_10_d_eq_j_plus_1 : (unfolding 10 (j+1) j).rank = 5`
  for every `j ≥ 1`. All sorry-free; `#print axioms` confirms only
  `propext, Classical.choice, Quot.sound`. The matrix is `10^j × 10`;
  live columns `{0, 2, 6, 8}` (residues `1, 3, 7, 9 (mod 10)`) plus the
  dead column `1` (`chiP 2 = 1` at row `0`) give `R = φ(10) + 1 = 5`.

  **Refutation of S128/S129 obstruction claim.** S128/S129 listed W=10
  as "structurally obstructed" in the same paragraph with W ∈ {7, 9, 11},
  citing a "multiplicity-2 residue pattern" argument. The S144 DP-based
  search refuted this. Permutation `ρ ↦ (1, 0, 4, 3, 9)`, `σ ↦ (8, 1, 2,
  0, 6)` triangularises the `5 × 5` submatrix to upper triangular with
  diagonal primes `{19, 2, 43, 31, 97}` and below-diagonal composites
  `{9, 32, 33, 39, 42, 49, 91, 92, 93, 99}`. **Row 9** is the key — the
  earlier search apparently restricted to row prefixes `{0, 1, 2, 3, 4}`,
  which is indeed insufficient.

  **Determinant via `det_of_upperTriangular`.** Same skeleton as W=8/W=12:
  pre-compute the 5 diagonal entries (each `= 1`) and the 10 below-
  diagonal entries (each `= 0`), establish `BlockTriangular id` by
  `fin_cases i <;> fin_cases k`, apply `det_of_upperTriangular`,
  expand via `Fin.prod_univ_five`, finish with `norm_num`.

  **Upper-bound subtlety:** as with W=4/6/8/12, `rank_le_width` gives
  only `rank ≤ 10`, not the sharp `rank ≤ 5 = φ(10) + 1`. We cite the
  general `upper_bound`, which evaluates to `φ(10) · 10^0 + 1 = 5`.

  Uses `Nat.prime_two` (existing), `chiP_nineteen_eq_one` (S117),
  `chiP_thirty_one_eq_one` (S122), `chiP_forty_three_eq_one` (S137),
  and one new helper `chiP_ninety_seven_eq_one` (S144,
  `decide`-checkable). ~310 Lean lines. **Eleventh unconditional
  `mps_bond_dim` instance; ninth instance over a wheel `W ≥ 3`; sixth
  instance using `det_of_upperTriangular`.** **First instance refuting
  an entry on the S128/S129 "structurally obstructed" list.**

  **S144 leading-row enumeration (definitive map).** A DP-based search
  over `W ∈ [2, 72]` with `R = φ(W) + 1 ≤ 22` (script:
  `leading_row_search.py`) shows the leading-row + dead-col upper-
  triangulation route closes **exactly**:
  - `W ∈ {2, 3, 4, 5, 6, 8, 10, 12, 18, 20}` — algorithmically reachable.

  All other tested wheels are structurally obstructed:
  - `{7, 9, 11, 13, 14, 15, 16, 17, 19, 21, 22, 24, 25, 26, 27, 28, 30,
    32, 33, 34, 36, 38, 40, 42, 44, 48, 50, 54, 60, 66}` — DP-confirmed
    no upper-triangulation in rows `[0, W)` with `R ≤ 22`.

  Wheels with `R > 22` (W ∈ {23, 29, 31, 35, 37, …}) remain untested
  but the obstruction pattern (every wheel between W=21 and W=66 with
  `R ≤ 22` is obstructed) suggests the leading-row family is essentially
  exhausted at the manageable-`R` parameter range. The **next single-
  session paths require either:**
  (a) `Matrix.det_of_blockTriangular` for non-triangulizable matrices
      (multi-session sub-arc; would unlock W ∈ {7, 9, 11, 14, 15, 16,
      24, 30, …} collectively);
  (b) Higher-`R` corners `(W ≥ 23, R ≥ 22)` — would need new local
      `prod_univ_X` lemmas plus heartbeat-bumping (cf. S143's
      `prod_univ_nine'` and `set_option maxHeartbeats 2000000`); but
      these are likely also structurally obstructed if the W ≤ 66 trend
      continues;
  (c) Cofactor-expansion-based determinant proofs for non-triangular
      sub-matrices (e.g., W=9's 7×7 sub-matrix has det = 1 via row
      cofactor expansion but no upper-triangular permutation; ~250
      Lean lines per such corner).

* **Route A^{(12)} (orthogonal corner, `W = 9`, `det_fromBlocks_zero₂₁` route). DONE S152.**
  Closed `exists_invertible_submatrix` for `(W = 9, d = j + 1)` via six
  new declarations: `chiP_thirteen_eq_one`, `chiP_forty_one_eq_one`,
  `chiP_fifty_three_eq_one`, `chiP_sixty_one_eq_one`,
  `exists_invertible_submatrix_W_eq_9_d_eq_j_plus_1`, and
  `mps_bond_dim_W_eq_9_d_eq_j_plus_1 : (unfolding 9 (j+1) j).rank = 7`
  for every `j ≥ 1`. All sorry-free; `#print axioms` confirms only
  `propext, Classical.choice, Quot.sound`. The matrix is `9^j × 9`;
  live columns are `{0, 2, 4, 6}` (residues `1, 3, 5, 7 (mod 9)` —
  φ(9) = 6, but only 4 are coprime to 9 because 9 = 3² and residues
  divisible by 3 are dead). Wait: φ(9) = 6, live residues are
  `{1, 2, 4, 5, 7, 8}` so live cols `{0, 1, 3, 4, 6, 7}`. We pick
  six live cols plus dead col 2 (`chiP 3 = 1` at row 0), totalling
  `R = 7`.

  **Triangulation IMPOSSIBLE; block-DIAGONAL `(1 + 3 + 3)` works.**
  The S144 DP-based search confirmed the 9×9 `j=1` slab admits NO
  leading-row + dead-col upper-triangulation in rows `[0, 9)`. The
  S151 pre-search (Python script `w9_blocktriangular_search.py`)
  enumerated `(1 + 3 + 3)` block-DIAGONAL decompositions and found 32
  valid candidates. The minimum-new-helpers choice picks
  `ρ ↦ (0, 1, 3, 5, 2, 4, 6)` and `σ ↦ (2, 1, 3, 7, 0, 4, 6)`, requiring
  only 4 new chiP helpers `{13, 41, 53, 61}`.

  After this permutation, the 7×7 submatrix has the structure:
  ```
     1 1 0 0 | 0 1 1     ρ=0  ← block-1 row (1×1: just chiP 3 = 1)
     0 1 1 1 | 0 0 0     ρ=1
     0 1 1 0 | 0 0 0     ρ=3  ← block-2 (3×3): primes {11, 13, 17, 29, 31, 47, 53}
     0 1 0 1 | 0 0 0     ρ=5
     ----------+------
     0 0 0 0 | 1 1 0     ρ=2
     0 0 0 0 | 1 1 1     ρ=4  ← block-3 (3×3): primes {19, 23, 37, 41, 43, 59, 61}
     0 0 0 0 | 0 1 1     ρ=6
  ```
  with off-diagonal blocks: row 0 → cols 1-3 = `(1,0,0)`, row 0 → cols
  4-6 = `(0,1,1)`, rows 1-3 → cols 4-6 = all zero, rows 4-6 → col 0 =
  all zero, rows 4-6 → cols 1-3 = all zero.

  **Determinant via NESTED `det_fromBlocks_zero₂₁`.** The Lean proof
  decomposes the 7×7 in two stages:
  1. **Layer 1** (1 + 6 split via `finSumFinEquiv : Fin 1 ⊕ Fin 6 ≃ Fin 7`):
     `Mexp = fromBlocks A B 0 D` with `A : Matrix (Fin 1) (Fin 1) ℚ`
     (just `[[1]]`), `B : Matrix (Fin 1) (Fin 6) ℚ` = `[[1, 0, 0, 0, 1, 1]]`,
     and `D : Matrix (Fin 6) (Fin 6) ℚ` (the 6×6 block-diagonal core).
     Lower-left zero ✓ (col 0 of rows 1-6 is all zero).
     `Mexp.det = A.det * D.det` via `det_fromBlocks_zero₂₁`.
     `A.det = 1` via `det_fin_one`.
  2. **Layer 2** (3 + 3 split via `finSumFinEquiv : Fin 3 ⊕ Fin 3 ≃ Fin 6`):
     `D = fromBlocks D1 0 0 D2` (BOTH off-diagonals zero, but we use
     `det_fromBlocks_zero₂₁` since lower-left is zero). `D.det = D1.det * D2.det`.
     `D1.det = -1` and `D2.det = -1` via `det_fin_three`.
  3. Combine: `Mexp.det = 1 · ((-1) · (-1)) = 1`, hence `IsUnit`.

  **Avoids the 4×4 det issue.** A first attempt used a 4 + 3 split for
  Layer 1, but the 4×4 block A had a non-trivial determinant (`det = -1`)
  requiring `det_succ_column_zero` + simp expansion. Mathlib has no
  `det_fin_four`, and `simp` blows up with `maxRecursionDepth` errors
  on the cofactor expansion + 4-fold sum. The 1+6 split sidesteps this
  by reducing all block-determinant computations to `det_fin_one`
  (1×1) and `det_fin_three` (3×3).

  **Layer-1 and Layer-2 equality proofs.** Each fromBlocks equality
  is proved by `Matrix.ext + rcases i / j with .. | .. <;> fin_cases ..
  <;> rfl` — 49 cases for layer 1 (Fin 7 = Fin 1 ⊕ Fin 6 reindex) and
  36 cases for layer 2 (Fin 6 = Fin 3 ⊕ Fin 3 reindex). Total 85 case
  checks, each closed by `rfl` since both sides are concrete ℚ literals.

  Uses `Nat.prime_two, prime_three, prime_five, prime_seven, prime_eleven,
  prime_seventeen, prime_nineteen, prime_twenty_three, prime_twenty_nine,
  prime_thirty_one, prime_thirty_seven, prime_forty_three, prime_forty_seven,
  prime_fifty_nine` (existing) and four new
  `chiP_X_eq_one` helpers for `X ∈ {13, 41, 53, 61}` (S152,
  `decide`-checkable). Plus `decide`-based non-primality for 31 composites
  in `[1, 62]`. ~600 Lean lines (49 entry-lemmas + 85 fromBlocks case
  checks + structural assembly).

  **Eleventh unconditional `mps_bond_dim` instance; tenth instance over a
  wheel `W ≥ 3`; FIRST instance using `Matrix.det_fromBlocks_zero₂₁`** —
  a determinant technique orthogonal to the previous nine corner closures
  (W ∈ {2, 3, 4, 5, 6, 8, 10, 12, 18, 20}, all using `det_of_upperTriangular`).
  **First closure of an S128/S129/S144 "block-triangular-required" wheel.**

  **Refutes the S151 "deferred" claim.** S151's note said the Lean
  implementation needs a separate session and ~600 lines. S152 confirms
  the line estimate but executes in one session via the nested-fromBlocks
  trick (avoiding the 4×4 det). The S151 LESSON about `rw [h_sub]`
  motive failures still holds — the entry-by-entry approach is the
  workaround.

  **Once-W=9 is closed (now actually true):** Route A^{(12)} becomes the
  template for W ∈ {7, 11, 14, 15, 16, 24, 30, ...} (all S144-obstructed
  cases). Each may require its own block-decomposition search (Python
  pre-search analogous to S151's). The Lean pattern (nested
  `det_fromBlocks_zero₂₁` with 1+(n) outer split + `det_fin_one` /
  `det_fin_three` per block) is reusable. Higher `R` values may require
  more nested splits.

* **Route A^{(13)} (orthogonal corner, `W = 7`, `det_fromBlocks_zero₂₁` route). DONE S159.**
  Closed `exists_invertible_submatrix` for `(W = 7, d = j + 1)` via two
  new declarations: `exists_invertible_submatrix_W_eq_7_d_eq_j_plus_1` and
  `mps_bond_dim_W_eq_7_d_eq_j_plus_1 : (unfolding 7 (j+1) j).rank = 7`
  for every `j ≥ 1`. **Sorry-free; ZERO new prime helpers** —
  every prime in the 7×7 submatrix
  (`{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47}`) was already
  declared from prior corner closures (S98, S99, S106, S117, S122,
  S128, S152). `#print axioms` confirms only `propext, Classical.choice,
  Quot.sound`.

  **Pre-search: 2 candidates, both block-DIAGONAL.** The S159 pre-search
  (`w7_blocktriangular_search.py`) on the W=7 7×7 slab found exactly two
  valid `1+(3+3)` BlockTriangular decompositions; both block-diagonal,
  both needing zero new prime helpers. The chosen permutation is
  `ρ ↦ (0, 1, 3, 5, 2, 4, 6)` and `σ ↦ (6, 1, 3, 5, 0, 2, 4)`. The matrix
  is `7^j × 7`; live cols are `{0, 1, 2, 3, 4, 5}` (W=7 is prime, so all
  six are coprime to 7) plus dead col `6` (the self-prime witness `chiP 7
  = 1` at row 0). `R = φ(7) + 1 = 7`.

  **Triangulation IMPOSSIBLE; block-DIAGONAL `(1 + 3 + 3)` works** with
  the same Lean assembly as W=9 (S152). Layer-1 split via
  `finSumFinEquiv : Fin 1 ⊕ Fin 6 ≃ Fin 7` gives `A : Matrix (Fin 1) (Fin
  1) ℚ = [[1]]`, `B : Matrix (Fin 1) (Fin 6) ℚ = [[1, 0, 0, 0, 1, 1]]`,
  and `D : Matrix (Fin 6) (Fin 6) ℚ` (the 6×6 block-diagonal core).
  Layer-2 split via `finSumFinEquiv : Fin 3 ⊕ Fin 3 ≃ Fin 6` gives `D₁
  = [[0,1,1],[1,0,0],[1,0,1]]` (det `-1`) and `D₂ =
  [[0,1,1],[1,1,0],[1,0,1]]` (det `-2`). Total
  `Mexp.det = 1 · ((-1) · (-2)) = 2`.

  **First `det ≠ ±1` corner; closes via `Ne.isUnit`.** Unlike all prior
  corners (W=2..6, 8..10, 12, 18, 20 give total `det = ±1`), W=7's total
  det is `2`. The closing step changes from `IsUnit 1` (handled by
  `norm_num` on `(1 : ℚ).IsUnit`) to `Ne.isUnit` from `(2 : ℚ) ≠ 0`.
  Demonstrates the BlockTriangular template is robust to nontrivial
  determinant magnitudes.

  **First corner with ZERO new helpers.** Helper-amortisation milestone:
  the `chiP_*_eq_one` library has reached steady-state reuse. Ten new
  proof-internal `h_not_prime_X` declarations (composites `{6, 9, 15, 18,
  24, 27, 33, 36, 42, 45}` not seen in W=9). ~600 Lean lines (49 entry-
  lemmas + 85 fromBlocks case checks + structural assembly + ten extra
  composite non-primality declarations).

  **Twelfth unconditional `mps_bond_dim` instance; eleventh instance
  over a wheel `W ≥ 3`; second instance using
  `Matrix.det_fromBlocks_zero₂₁`** (after S152 W=9). Closes the second
  S128/S129/S144 "block-triangular-required" wheel. Refutes S128's
  prediction that W=7 needed "a multi-session new technique" — actual
  closure is single-session via S152's nested-fromBlocks pattern with
  a fresh pre-search.

* **Route A'''' (orthogonal corner, `W = 4`). DONE S107.**
  Closed `exists_invertible_submatrix` for `(W = 4, d = j + 1)` via
  three new declarations: `chiP_eleven_eq_one`,
  `exists_invertible_submatrix_W_eq_4_d_eq_j_plus_1`, and
  `mps_bond_dim_W_eq_4_d_eq_j_plus_1 : (unfolding 4 (j+1) j).rank = 3`
  for every `j ≥ 1`. All sorry-free; `#print axioms` confirms only
  `propext, Classical.choice, Quot.sound`. The matrix has `4^(d-j) =
  4` columns; column `3` is `chiP` at multiples of `4` (all zeros), so
  we drop it and pick the three live columns `{0, 1, 2}`. With rows
  `{0, 1, 2}` (available since `4^j ≥ 4` for `j ≥ 1`) the `3 × 3`
  submatrix is
  ```
     ⎡ chiP 1, chiP 2, chiP 3  ⎤   ⎡ 0, 1, 1 ⎤
     ⎢ chiP 5, chiP 6, chiP 7  ⎥ = ⎢ 1, 0, 1 ⎥
     ⎣ chiP 9, chiP 10, chiP 11⎦   ⎣ 0, 0, 1 ⎦
  ```
  of determinant `−1` (computed via `Matrix.det_fin_three`). Uses only
  `Nat.prime_two`, `Nat.prime_three`, `Nat.prime_five`,
  `Nat.prime_seven`, `Nat.prime_eleven` (the last new at S107, all
  `decide`-checkable) and decidability of `¬ Nat.Prime` for `1, 4, 6,
  9, 10`. **Upper-bound subtlety:** `rank_le_width` gives only
  `rank ≤ 4`, not the sharp `rank ≤ 3`; we therefore cite the general
  `upper_bound` lemma, which evaluates to `φ(4) · 4^0 + 1 = 3` here.
  This is the **first orthogonal-corner instance where the live-column
  count strictly beats the column count**, so `upper_bound` becomes
  necessary; subsequent `W ∈ {6, 7, …}` corners will follow the same
  upper-bound pattern. ~190 Lean lines. **Second unconditional
  `mps_bond_dim` instance over a wheel `W ≥ 3`.**

* **Route B (Vandermonde / character).** Re-prove the lower bound by
  a non-density argument. The Arc 4 follow-on note conjectures that
  the spike eigenvectors of `M^T M` are Dirichlet character vectors
  restricted to the `chiP` support. If true, the rank lower bound
  would follow from the linear independence of distinct characters
  (a Vandermonde-style fact in cyclotomic fields). This sidesteps
  prime-density but introduces character theory infrastructure.
  Speculative; needs the conjecture to be formalised first.

* **Route C (PNT).** Mathlib does have a quantitative prime number
  theorem (`PrimeNumberTheorem`). For the small-density regime
  (`R ≪ x / log x`), PNT implies enough primes in the right intervals.
  For the high-density regime (e.g. half-cut `R = φ(W)·W^(d/2-1)`),
  we still need short-interval primes. Route C closes the easy regime
  but not the saturating one.

The structural reduction `lower_bound ← exists_invertible_submatrix`
(closed in S83 via `Matrix.rank_of_isUnit` and `Matrix.rank_submatrix_le`)
is mechanical and never needs to be revisited.

**File structure note (S76, refined S83):** the file is now ordered
auxiliary-lemmas-first (`rank_le_min_dim`, `row_support_coprime`,
`live_columns_count`, `upper_bound`, `exists_invertible_submatrix`,
`lower_bound`), followed by the main theorem `mps_bond_dim`. The
opening docstring/decomposition section makes the architecture
explicit. The `exists_invertible_submatrix` declaration sits between
`upper_bound` and `lower_bound` and is the only place where prime
density is invoked.

## Falsification

This formalisation is falsified if any session produces a `lake build`
failure on `MPSBondDim/Basic.lean`. The Lean kernel is the falsifier.

## Why this matters

E2.1 is the cleanest of the L1–L5 queue: no analytic number theory
(no PNT, no zeros), no circuit complexity, no Kt-randomness. The proof
is ~50 lines of informal mathematics. Completing it is a proof-of-concept
that the project's Lean track can produce verified artifacts in a
session-tractable timeframe.
