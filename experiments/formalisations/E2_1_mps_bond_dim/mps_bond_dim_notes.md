# L1 — E2.1 MPS bond-dim identity (Lean 4 formalisation, in progress)

**Lean source:** `MPSBondDim/MPSBondDim/Basic.lean`
**Toolchain:** `leanprover/lean4:v4.30.0-rc2` + mathlib `v4.30.0-rc2` (lake project under `MPSBondDim/`).
**Build status:** `lake build` succeeds. **1 `sorry` placeholder remains, isolated to a pure prime-density existential** (S83 reorganised the file: `lower_bound` is now closed by a structural reduction to `exists_invertible_submatrix`, which is the new home of the only `sorry`. S90 added the trivial floor `1 ≤ rank` as a separate lemma. **S98 closed the corner case `(W = 2, j = 1)` unconditionally via Bertrand** — `mps_bond_dim_W_eq_2_j_eq_1 : (unfolding 2 d 1).rank = 2` for every `d ≥ 2`, sorry-free. **S99 closed the orthogonal corner case `(W = 2, d = j + 1)` *without even needing Bertrand*** — `mps_bond_dim_W_eq_2_d_eq_j_plus_1 : (unfolding 2 (j+1) j).rank = 2` for every `j ≥ 1`, using only `Nat.prime_two` and `Nat.prime_three`. **S106 extended the orthogonal-corner argument to `W = 3`** — `mps_bond_dim_W_eq_3_d_eq_j_plus_1 : (unfolding 3 (j+1) j).rank = 3` for every `j ≥ 1`, using `Matrix.det_fin_three` and the explicit primes `2, 3, 5, 7`; the unit witness is `isUnit_one.neg : IsUnit (-1 : ℚ)`. **S107 extended the orthogonal-corner argument to `W = 4`** — `mps_bond_dim_W_eq_4_d_eq_j_plus_1 : (unfolding 4 (j+1) j).rank = 3` for every `j ≥ 1`, using `Matrix.det_fin_three`, the explicit primes `2, 3, 5, 7, 11`, and (for the upper bound) the general `upper_bound` lemma — the first orthogonal-corner instance where `rank_le_width` is not tight (it gives only `rank ≤ 4`, not the sharp `rank ≤ 3 = φ(4) · 4^0 + 1`). The general-case `sorry` is unchanged.)
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
