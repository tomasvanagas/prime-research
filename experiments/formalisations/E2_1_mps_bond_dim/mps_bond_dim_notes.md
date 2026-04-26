# L1 — E2.1 MPS bond-dim identity (Lean 4 formalisation, in progress)

**Lean source:** `MPSBondDim/MPSBondDim/Basic.lean`
**Toolchain:** `leanprover/lean4:v4.30.0-rc2` + mathlib `v4.30.0-rc2` (lake project under `MPSBondDim/`).
**Build status:** `lake build` succeeds. **1 `sorry` placeholder remains** (down from 3 after S76 closed `upper_bound` and reduced `mps_bond_dim` to a term application of the lemmas).
**No `axiom` introductions.**

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

The proof is decomposed into 6 declarations:

| name                   | role                                       | status   |
|------------------------|--------------------------------------------|----------|
| `chiP`                 | prime indicator `ℕ → ℚ`                    | def      |
| `unfolding`            | the `(W^j × W^(d-j))` matrix over `ℚ`      | def      |
| `rank_le_min_dim`      | trivial ceiling `rank ≤ W^j`               | **done** |
| `row_support_coprime`  | nonzero entries imply `gcd(k+1, W) = 1`    | **done** |
| `live_columns_count`   | CRT count: live columns = `φ(W) · W^(d-j-1)` | **done** |
| `upper_bound`          | `rank ≤ φ(W) · W^(d-j-1) + 1`              | **done** |
| `lower_bound`          | `min(W^j, φ(W) · W^(d-j-1) + 1) ≤ rank`    | `sorry`  |
| `mps_bond_dim`         | **main theorem** (term application)        | **done** |

(Note: `mps_bond_dim` itself contains no `sorry`; it is a 3-line
`Nat.le_antisymm` of the auxiliary lemmas. But it transitively depends
on `lower_bound`, which still has a `sorry`. Once `lower_bound` is
closed, the main theorem is automatically closed without modification.)

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

**`lower_bound`** — the harder side. Informally, exhibit a set of
`min(W^j, φ(W) · W^(d-j-1) + 1)` rows whose restrictions to the live
columns are linearly independent over ℚ. The informal argument
hand-waves over a prime-counting density; making this rigorous in Lean
will probably require either (i) a counting-via-PNT step or (ii) a
clever reformulation that reduces to a generic "many rows with distinct
non-zero patterns" combinatorial lemma.

**File structure note (S76):** the file was reorganised so the lemmas
(`rank_le_min_dim`, `row_support_coprime`, `live_columns_count`,
`upper_bound`, `lower_bound`) now appear before the main theorem
`mps_bond_dim`. This was needed for `mps_bond_dim`'s term-mode proof
to refer to the lemmas. The narrative-first ordering (statement at
the top) was traded for proof structure; the file's opening
docstring/decomposition section makes the architecture explicit.

## Falsification

This formalisation is falsified if any session produces a `lake build`
failure on `MPSBondDim/Basic.lean`. The Lean kernel is the falsifier.

## Why this matters

E2.1 is the cleanest of the L1–L5 queue: no analytic number theory
(no PNT, no zeros), no circuit complexity, no Kt-randomness. The proof
is ~50 lines of informal mathematics. Completing it is a proof-of-concept
that the project's Lean track can produce verified artifacts in a
session-tractable timeframe.
