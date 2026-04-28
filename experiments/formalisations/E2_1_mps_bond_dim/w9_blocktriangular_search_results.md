# W=9 BlockTriangular pre-search (S151)

## Summary

The S144 leading-row enumeration confirmed that the `(W=9, d=j+1)`
orthogonal-corner case admits **no leading-row + dead-col upper-
triangulation** with rows in `[0, 9)` (the binding constraint at `j=1`).

This script confirms a **different** structural shape closes W=9: the
matrix is **block-DIAGONAL as 1 + 3 + 3** under a specific row/col
permutation. Each `3 × 3` diagonal block has `det = ±1`, computable in
Lean via `Matrix.det_fin_three`.

This unlocks **Route A^{(12)}** (a new technique in the file): nest two
applications of `Matrix.det_fromBlocks_zero_21` (1+6 then 3+3) and
compute the two `3 × 3` block dets via `det_fin_three`. The total is
`1 · (-1) · (-1) = 1`.

## Numerical artefacts

```
Found 32 valid (1+3+3) BlockTriangular decompositions for W=9 corner.
Minimum new chiP helpers needed: 4
Block-diagonal-only candidates: 32 (every viable decomposition has
   the upper-right 3x3 block also zero — so every decomposition is
   actually block-DIAGONAL, not just block-upper-triangular)

Top candidate (minimum new chiP helpers):
  rho = (0, 1, 3, 5, 2, 4, 6)
  sigma = (2, 1, 3, 7, 0, 4, 6)
  A-block det = -1, B-block det = -1
  new primes needed: {13, 41, 53, 61}
```

The full 7×7 submatrix (over chiP, i.e. {0, 1}-valued):
```
   1 1 0 0 | 0 1 1     ρ=0
   0 1 1 1 | 0 0 0     ρ=1
   0 1 1 0 | 0 0 0     ρ=3
   0 1 0 1 | 0 0 0     ρ=5
   ----------+------
   0 0 0 0 | 1 1 0     ρ=2
   0 0 0 0 | 1 1 1     ρ=4
   0 0 0 0 | 0 1 1     ρ=6
```
With determinant exactly `1` (verified via SymPy).

## Block decomposition

After applying `Matrix.det_succ_column_zero` to peel off the `(0, 0)`
entry (the dead-col witness `chiP 3 = 1`), the remaining `6 × 6` minor:
```
   1 1 1 | 0 0 0
   1 1 0 | 0 0 0
   1 0 1 | 0 0 0
   ------+------
   0 0 0 | 1 1 0
   0 0 0 | 1 1 1
   0 0 0 | 0 1 1
```
This is **block-diagonal** as `3 + 3`: both upper-right and lower-left
`3 × 3` blocks are zero. Apply `Matrix.det_fromBlocks_zero_21` (after
reindexing via `finSumFinEquiv : Fin 3 ⊕ Fin 3 ≃ Fin 6`) to get
`det(6×6) = det(B1) * det(B2)`.

Block 1 (rows 1, 3, 5 × cols 1, 3, 7):
```
B1 = [chiP 11, chiP 13, chiP 17;     [1, 1, 1;
      chiP 29, chiP 31, chiP 35;  =   1, 1, 0;
      chiP 47, chiP 49, chiP 53]      1, 0, 1]
det(B1) = -1 (via det_fin_three)
```

Block 2 (rows 2, 4, 6 × cols 0, 4, 6):
```
B2 = [chiP 19, chiP 23, chiP 25;     [1, 1, 0;
      chiP 37, chiP 41, chiP 43;  =   1, 1, 1;
      chiP 55, chiP 59, chiP 61]      0, 1, 1]
det(B2) = -1 (via det_fin_three)
```

## Lean implementation status (S151)

**Status: BLOCKED on Lean tactic complexity.**

The mathematical structure is fully validated in Python and
straightforward in principle: two `det_fromBlocks_zero_21`
applications + two `det_fin_three` calls. However, the Lean
implementation hit a **`rw [h_sub] motive not type-correct`**
issue when attempting to prove the submatrix equals an explicit
`!![...]` matrix via `Matrix.ext + fin_cases <;> fin_cases`.

The root cause: after `fin_cases`, the goal's σ Fin values carry
proof terms `_ < 9^((j+1)-j)` whose types depend on `(j+1)-j`.
Rewriting `(j+1)-j → 1` via `rw` requires the motive to be type-
correct over the dependent expression, which fails when the proof
terms appear in the goal's underlying term.

The fix: do **49 individual `have` lemmas** in the W=20 style — each
of form `unfolding 9 (j+1) j ⟨r, _⟩ ⟨s, _⟩ = const`, with each lemma
using the standard `change chiP (...) = const; rw [h_sub]; norm_num;
exact chiP_X_eq_one` pattern. Each lemma is ~6 lines. 49 lemmas =
~300 lines + the structural assembly = ~600 lines total.

This is achievable in a follow-up session.

## What this script enables for the next session

1. **Concrete permutation:** the next agent doesn't need to re-search.
   Use `ρ = (0, 1, 3, 5, 2, 4, 6)`, `σ = (2, 1, 3, 7, 0, 4, 6)`.
2. **chiP helper requirements:** new helpers needed are
   `chiP_thirteen_eq_one`, `chiP_forty_one_eq_one`,
   `chiP_fifty_three_eq_one`, `chiP_sixty_one_eq_one` — each a
   one-liner via `Nat.Prime` + `decide`.
3. **Composite decidability:** the 31 distinct composites used are
   `{1, 4, 8, 10, 12, 14, 16, 20, 21, 22, 25, 26, 28, 30, 32, 34,
   35, 38, 39, 40, 44, 46, 48, 49, 50, 52, 55, 56, 57, 58, 62}`.
4. **Block dets are -1, -1.** Total det is +1.

## What this would do for Arc 2

Closing W=9 unconditionally introduces the **first BlockTriangular
non-id technique** in the file (i.e. the first time we use
`det_fromBlocks_zero_21` rather than `det_of_upperTriangular`).

This unlocks the full `det_of_blockTriangular`-required set:
`{7, 9, 11, 13, 14, 15, 16, 17, 19, 21, ...}` — many of which can
be closed by similar 1+3+3 or 1+a+b BlockTriangular decompositions
with appropriate row/col permutations.

## Falsifying conditions

- If the permutation `ρ = (0, 1, 3, 5, 2, 4, 6), σ = (2, 1, 3, 7, 0, 4, 6)`
  does NOT give a BlockTriangular submatrix with the claimed det, this
  document is wrong. (Falsified by SymPy: det = 1.)
- If 4 new chiP helpers (13, 41, 53, 61) are NOT sufficient, this
  document is wrong. (Falsified by enumeration: every helper used in
  the 1+3+3 decomposition is in the existing set ∪ {13, 41, 53, 61}.)
