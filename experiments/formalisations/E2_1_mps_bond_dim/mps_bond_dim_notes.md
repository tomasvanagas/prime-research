# L1 — E2.1 MPS bond-dim identity (Lean 4 formalisation, in progress)

**Lean source:** `MPSBondDim/MPSBondDim/Basic.lean`
**Toolchain:** `leanprover/lean4:v4.30.0-rc2` + mathlib `v4.30.0-rc2` (lake project under `MPSBondDim/`).
**Build status:** `lake build` succeeds. 4 `sorry` placeholders remain (down from 5 in the initial skeleton).
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
| `mps_bond_dim`         | **main theorem**                           | `sorry`  |
| `upper_bound`          | `rank ≤ φ(W) · W^(d-j-1) + 1`              | `sorry`  |
| `lower_bound`          | `min(W^j, φ(W) · W^(d-j-1) + 1) ≤ rank`    | `sorry`  |
| `rank_le_min_dim`      | trivial ceiling `rank ≤ W^j`               | **done** |
| `row_support_coprime`  | nonzero entries imply `gcd(k+1, W) = 1`    | **done** |
| `live_columns_count`   | CRT count: live columns = `φ(W) · W^(d-j-1)` | `sorry`  |

The two completed proofs:

1. `rank_le_min_dim` is one line: a direct citation of mathlib's
   `Matrix.rank_le_height`.
2. `row_support_coprime` is a 30-line number-theoretic argument: nonzero
   entry ⇒ `n = i·W^(d-j) + k + 1` is prime; `i ≥ 1` and `j < d` give
   `n > W`; hence no prime factor of `W` divides `n`, so `gcd(n, W) = 1`;
   reducing mod `W` (using `gcd_add_mul_right_left` after rewriting
   `W^(d-j) = W^(d-j-1) · W`) yields `gcd(k+1, W) = 1`.

## What the next session needs to do

In rough order of decreasing tractability:

1. **`live_columns_count`** — the CRT count
   `#{k ∈ [0, W^(d-j)) : gcd(k+1, W) = 1} = φ(W) · W^(d-j-1)`.
   Periodicity in `W`-blocks; each block contains exactly `φ(W)` live
   columns (a one-`W`-block lemma plus an induction on `d-j-1`).
   Mathlib has `Nat.totient_eq_card_lt_and_coprime`, `Nat.totient` lemmas,
   and `Finset.filter` + `Finset.card_eq_sum_ones` machinery.

2. **`upper_bound`** — combines `row_support_coprime` and
   `live_columns_count`. The structural step: rows `i ≥ 1` of `unfolding W d j`
   live in the subspace `V := {v : Fin(W^(d-j)) → ℚ | ∀ k, gcd(k+1, W) ≠ 1 → v k = 0}`,
   which has dimension `φ(W) · W^(d-j-1)`. Adding row `0` gives at most
   one further dimension. Use `Matrix.rank_le` plus
   `Submodule.finrank_le_finrank_of_le`.

3. **`lower_bound`** — the harder side. Informally, exhibit a set of
   `min(W^j, φ(W) · W^(d-j-1) + 1)` rows whose restrictions to the live
   columns are linearly independent. The informal argument hand-waves
   over a prime-counting density; making this rigorous in Lean will
   probably require either (i) a counting-via-PNT step or (ii) a clever
   reformulation that reduces to a generic "many rows with distinct
   non-zero patterns" combinatorial lemma.

4. **`mps_bond_dim`** — `Nat.le_antisymm` of the two bounds, plus a small
   case split between the two arguments to `min` (using
   `rank_le_min_dim` to handle the `W^j ≤ φ(W) · W^(d-j-1) + 1` case).

## Falsification

This formalisation is falsified if any session produces a `lake build`
failure on `MPSBondDim/Basic.lean`. The Lean kernel is the falsifier.

## Why this matters

E2.1 is the cleanest of the L1–L5 queue: no analytic number theory
(no PNT, no zeros), no circuit complexity, no Kt-randomness. The proof
is ~50 lines of informal mathematics. Completing it is a proof-of-concept
that the project's Lean track can produce verified artifacts in a
session-tractable timeframe.
