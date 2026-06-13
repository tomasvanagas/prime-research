# w11_odd_block_atomicity — results

**Result (S206):** The W=11 odd 5×5 block (rows `{1, 3, 5, 7, 9}` ×
cols `{1, 3, 5, 7, 9}`) is **atomically block-irreducible**:
exhaustive enumeration over all 15 ordered partitions of 5 with parts
≤ 4 yields **zero** BlockTriangular decompositions.

Block matrix:
```
        col 1  col 3  col 5  col 7  col 9
row 1:  [ 1,    0,    1,    1,    0  ]
row 3:  [ 0,    1,    0,    1,    1  ]
row 5:  [ 0,    1,    1,    0,    0  ]
row 7:  [ 1,    0,    1,    0,    0  ]
row 9:  [ 1,    1,    0,    1,    1  ]
```

det = 1, rank = 5. No row has fewer than 2 nonzero entries → no
leading-row triangulation. No (1, 4), (4, 1), (2, 3), (3, 2), (1, 1, 3),
(3, 1, 1), (1, 3, 1), (1, 2, 2), (2, 2, 1), (2, 1, 2), (1, 1, 1, 2),
(1, 1, 2, 1), (1, 2, 1, 1), (2, 1, 1, 1), or (1, 1, 1, 1, 1) admits a
BlockTriangular form.

**Falsification statement:** if any partition of 5 with parts ≤ 4
yielded a BlockTriangular decomposition of the W=11 odd 5×5 block, the
single-session arc-anticipated `1 + 5 + 5` route would close W=11 by
combining `det_fromBlocks_zero_21` with the smaller-block dets. The
empirical exhaustion confirms the decomposition does not exist.

**Implication:** Lean closure of W=11 via the `1 + 5 + 5` parity-
permutation route requires non-decomposable determinant computation
of the 5×5 odd block — i.e., 5×5 cofactor expansion (mathlib's
`Matrix.det_succ_column_zero` recursively, ~16 `det_fin_three`
invocations) or a new `Matrix.det_fin_five` lemma.

See `w11_blocktriangular_search_results.md` for the multi-session
paths forward.
