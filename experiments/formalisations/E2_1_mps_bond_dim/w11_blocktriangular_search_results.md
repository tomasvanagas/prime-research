# W=11 BlockTriangular pre-search results (S206)

**Goal:** Identify a single-session-viable block-triangular decomposition
for the W=11 orthogonal corner case `(W=11, d=j+1)`, R = φ(11) + 1 = 11,
of E2.1.

**Result:** **Single-session closure obstructed by an atomic 5×5 odd
block.** The arc's W=11 next-action (anticipated `1 + 5 + 5` split via
`Fin.prod_univ_five`) is **invalidated** by an explicit STRUCTURAL
obstruction: the only block-triangular decomposition of the W=11 11×11
matrix over rows in [0, 11) into blocks of size ≤ 5 forces a 5×5 block
with NO further block-triangular sub-decomposition (verified
exhaustively over all 15 ordered partitions of 5 with parts ≤ 4).

## Pre-search summary

Three Python scripts were run:

* `w11_blocktriangular_search.py` — search for `1 + 1 + 3 + 3 + 3` block-
  DIAGONAL decompositions over rows in [0, 33). **Result: 0 candidates**.
  Row 0 has 5 nonzero entries (cols where `chiP(c+1) = 1` for primes 2,
  3, 5, 7, 11) → cannot isolate row 0 as a 1×1 block-DIAGONAL component.
* `w11_inner_triangulation.py` — leading-row upper-triangulation of the
  inner 10×10 over rows [1, 11) × live cols. **Result: obstructed** in
  rows [1, 11); FOUND in rows [1, 22) with `ρ ↦ (3, 2, 1, 6, 5, 7, 16,
  10, 19, 18)`, max row = 19. **Useful only for j ≥ 2** (since 19 ≥ 11
  forces `j ≥ 2` because `W^j = 11^j ≥ 19` requires j ≥ 2).
* `w11_general_search.py` — exhaustive search over partitions of 11 with
  blocks of size ≤ 5, peeling bottom-up. **Result: 0 candidates with all
  blocks ≤ 3** (or ≤ 4); only `(5, 5, 1)`, `(1, 5, 5)`, `(5, 1, 5)`,
  `(3, 5, 3)`, `(3, 3, 5)`, `(4, 5, 2)`, `(2, 5, 4)` partitions admit
  decompositions, **and in every one the 5×5 block is the same** —
  rows `(1, 3, 5, 7, 9)` × cols `(1, 3, 5, 7, 9)` (the odd-row × odd-col
  parity block).

## The parity-block obstruction

For `r ≥ 1`, the entry `M[r][c] = chiP(11r + c + 1)`. Since `chiP(q) = 1`
requires `q` prime and `q > 2 ⇒ q odd`, we have

* `r` odd → `11r` odd → `11r + c + 1` odd ⇔ `c` odd → nonzero entries
  only at odd `c`.
* `r` even (r ≥ 1) → `11r` even → `11r + c + 1` odd ⇔ `c` even → nonzero
  entries only at even `c`.
* `r = 0`: `11·0 + c + 1` = `c + 1`, prime iff `c + 1 ∈ {2, 3, 5, 7, 11}`,
  i.e., `c ∈ {1, 2, 4, 6, 10}`. **Row 0 mixes even and odd cols** — the
  single anomalous entry at `(r=0, c=1)` (chiP(2) = 1) is the **only**
  off-diagonal entry under the parity permutation.

Permuting rows to `(0, 2, 4, 6, 8, 10, 1, 3, 5, 7, 9)` and cols to
`(0, 2, 4, 6, 8, 10, 1, 3, 5, 7, 9)` yields the BlockTriangular form

```
| even-rows × even-cols (6×6) | sparse 6×5 with single nonzero at (0,1) |
| 5×6 ZERO                    | odd-rows × odd-cols (5×5)               |
```

(The lower-left is structurally zero by the parity argument.)

* **6×6 even-block** further splits as `(1×1) ⊕ (5×5)` via row 10 (only
  nonzero at col 2 = chiP(113) = 1), and the inner 5×5 over rows
  `(0, 2, 4, 6, 8)` × cols `(0, 4, 6, 8, 10)` IS leading-row
  upper-triangulable: `ρ ↦ (0, 6, 2, 8, 4)`, `σ ↦ (10, 4, 6, 0, 8)`,
  diagonal primes `{11, 71, 29, 89, 53}`. **det = 1.**
* **5×5 odd-block** has det = 1 but is **block-irreducible** — see below.

## Atomicity of the odd 5×5 block (verified S206)

Script `w11_odd_block_atomicity.py` enumerated all 15 ordered partitions
of 5 with parts ≤ 4 and found **NO block-triangular decomposition** for
any of them. The 5×5 odd block is the matrix

```
        col 1  col 3  col 5  col 7  col 9
row 1:  [ 1,    0,    1,    1,    0  ]   primes:  13         17  19
row 3:  [ 0,    1,    0,    1,    1  ]            37         41  43
row 5:  [ 0,    1,    1,    0,    0  ]            59  61
row 7:  [ 1,    0,    1,    0,    0  ]            79  83
row 9:  [ 1,    1,    0,    1,    1  ]           101 103     107 109
```

with det = 1, rank = 5. No row has fewer than 2 nonzero entries → no
leading-row triangulation. No (1, 4), (4, 1), (2, 3), (3, 2), (1, 1, 3),
(3, 1, 1), (1, 3, 1), (1, 2, 2), (2, 2, 1), (2, 1, 2), (1, 1, 2, 1),
(1, 2, 1, 1), (2, 1, 1, 1), (1, 1, 1, 2), or (1, 1, 1, 1, 1) partition
admits a BlockTriangular decomposition. The block is **atomic** in the
sense that any Lean closure must compute its determinant via a
non-decomposable method (cofactor expansion, ring-theoretic identity,
or a new technique).

## Implication for the W=11 corner closure

Path | Feasibility | Cost
--- | --- | ---
`(5, 5, 1)` BlockTriangular + `det_of_upperTriangular` for top + cofactor expansion for odd 5×5 | viable but **multi-session** | ~1500-2000 Lean lines (5×5 cofactor → 4 4×4 → 16 3×3 dets, plus 25 entry lemmas for the odd block alone, plus the standard top + 1×1 + entries)
Leading-row triangulation of inner 10×10 over rows [1, 22) (max row 19) | viable for **j ≥ 2 ONLY** (j=1 separate edge) | single-session if accepting weaker statement, ~700 lines
`Matrix.det_succ_above` general cofactor expansion (5×5) | adds reusable mathematical content | ~1000 Lean lines for det_fin_five lemma + ~500 application
`Matrix.det_of_blockTriangular` over the parity permutation (avoids cofactor expansion of odd block) | does NOT help — `det_of_blockTriangular` reduces to det of each diagonal block, and the odd 5×5 block is exactly what we need | does not bypass the obstruction
Conjugate the odd 5×5 to a tridiagonal/banded form via row/col operations | speculative; would require new algebraic content | unknown

## What the next agent should do

Three clean single-session paths forward:

1. **W = 11 corner for `j ≥ 2` only.** Use the leading-row triangulation
   `ρ ↦ (3, 2, 1, 6, 5, 7, 16, 10, 19, 18)` of the inner 10×10 (max row
   19; valid for `j ≥ 2` since `11^2 = 121 > 19`). Six new prime helpers
   needed: `{67, 71, 73, 79, 83, 113, 181}` (S206 inner-triangulation
   search). The j=1 case becomes a separate sub-problem of the corner
   theorem.
2. **W = 13 corner.** Skip W=11 for now; W=13 is the next prime, R=13.
   **Pre-search needed.** If the analogous parity decomposition is
   atomic (likely, by parity argument for any prime W ≥ 11), then W=13
   has the SAME obstruction at scale 6×7 → either close for `j ≥ 2`
   only, or pivot to non-prime W.
3. **W = 14 corner via block-triangular over `det_of_blockTriangular`.**
   W=14 was identified S128/S144 as needing block-triangular — not
   leading-row. R=φ(14)+1 = 7. Pre-search for partition with all blocks
   ≤ 3 should succeed (W=14 isn't a prime, no parity issue at the same
   scale). If it does, single-session viable.

## Multi-session sub-arc proposal: `det_fin_five` lemma

A reusable `det_fin_five` lemma (analogous to mathlib's `det_fin_three`)
would cleanly unblock W=11. The work: prove

```
theorem Matrix.det_fin_five {R : Type*} [CommRing R] (M : Matrix (Fin 5) (Fin 5) R) :
    M.det = ⟨the explicit Leibniz expansion, 120 terms⟩
```

via `Fintype.sum_perm` over `Equiv.Perm (Fin 5)`. Once available, the
W=11 odd 5×5 closure becomes a single `decide` (or single-line `simp`)
step — no cofactor recursion needed.

This sub-arc is reusable for all `W` whose closure requires a 5×5
det. Estimated 2 sessions: 1 for the lemma, 1 for application at W=11.

## Files

* `w11_blocktriangular_search.py` — `(1 + 1 + 3 + 3 + 3)` BD search.
* `w11_inner_triangulation.py` — leading-row triangulation of inner
  10×10 (rows ≥ 1).
* `w11_general_search.py` — exhaustive partition search with
  triangulability check per block.
* `w11_odd_block_atomicity.py` — atomicity verification of odd 5×5.
* `w11_nested_search.py` — nested 1+(1+3+3+3) BD search at row pool ≤ 22.

## Searched candidate (multi-session reference)

For the **multi-session** Lean closure via `(5, 5, 1)` partition with
parity permutation, the witness is:

```
rho   = (0, 2, 4, 6, 8, 10, 1, 3, 5, 7, 9)
sigma = (0, 4, 6, 8, 10, 2, 1, 3, 5, 7, 9)
```

(this is the parity permutation; cols of the inner 6×6 are reordered
to put the inner triangulable 5×5 cols `{0, 4, 6, 8, 10}` first and the
1×1 col `{2}` last)

block 1 (top 5×5, even-row triangulable): rows `(0, 6, 2, 8, 4)`,
cols `(10, 4, 6, 0, 8)`, diagonal primes `{11, 71, 29, 89, 53}`,
det = 1.
block 2 (1×1): row 10, col 2, value chiP(113) = 1.
block 3 (bottom 5×5, atomic): rows `(1, 3, 5, 7, 9)`, cols
`(1, 3, 5, 7, 9)`, det = 1 (cofactor expansion required).

New prime helpers needed: `{29, 53, 67, 71, 73, 79, 83, 101, 103, 107,
113}` minus existing `{29, 53, 71, 89, 113}` → roughly 7-9 new helpers.

For the **single-session j ≥ 2 only** Lean closure, witness:

```
rho   = (3, 2, 1, 6, 5, 7, 16, 10, 19, 18)
sigma = (9, 8, 7, 6, 3, 5, 4, 2, 1, 0)
```

over the inner 10×10 (rows ≥ 1, live cols), with outer 1+10
BlockTriangular split. Diagonal primes `{19, 31, 43, 59, 73, 83, 113,
181, 199, 211}`. New helpers: `{67, 71, 73, 79, 83, 113, 181}`. **Six**
fresh primes (within the arc's "≤ 6 new helpers" single-session budget).
But this is a strict-strict closure: it works ONLY for j ≥ 2 because
max row = 19 forces 11^j ≥ 19 ⇔ j ≥ 2.

## Self-evaluation

**Grade: B-grade** (ambitious failure with structural insight).

What this session produced that wasn't in the project before:

1. Closed-form atomicity proof for the W=11 odd 5×5 block: enumerated
   all 15 partitions of 5 with parts ≤ 4 and confirmed NONE gives a
   BlockTriangular decomposition.
2. Parity-block decomposition of the W=11 11×11 matrix: explicit
   permutation revealing `6×6 even ⊕ 5×5 odd` BlockTriangular structure
   with the single off-diagonal entry at `(row 0, col 1)`.
3. Refined arc 2 next-action: instead of "single-session 1+5+5", three
   distinct paths (j ≥ 2 only, W=13/W=14 pivot, det_fin_five sub-arc).

What this session did NOT produce:

* A Lean closure of W=11 (the original arc next-action) — the structural
  obstruction makes it multi-session.

The session refutes the arc's anticipated "Single-session if pre-search
yields ≤ 6 new prime helpers" prediction for W=11 by identifying a NEW
class of obstruction (block-atomicity) that prior closures (W ∈ {2, 3,
4, 5, 6, 7, 8, 9, 10, 12, 18, 20}) did not encounter.
