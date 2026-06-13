# W=14 j=1 shape search — results

Results of `w14_j1_shape_search.py`. Efficient partition-based
enumeration of all block shapes (with parts ≤ 4) that admit a
BlockTriangular decomposition of the W=14 j=1 inner 6×6.

**Per-shape results (across all 587 invertible row sets):**

| Shape       | Hits | Lean closure feasible? |
|-------------|-----:|------------------------|
| (5, 1)      |   ≥3 | ✗ (no det_fin_five) |
| (1, 5)      |   ≥3 | ✗ (no det_fin_five) |
| (4, 2)      |    0 | — |
| (2, 4)      |   ≥3 | ✗ (no det_fin_four) |
| (3, 3)      |    0 | — |
| (4, 1, 1)   |    0 | — |
| (1, 4, 1)   |   ≥3 | ✗ (no det_fin_four) |
| (1, 1, 4)   |   ≥3 | ✗ (no det_fin_four) |
| (3, 2, 1)   |    0 | — |
| (3, 1, 2)   |    0 | — |
| **(2, 3, 1)** |  ≥3 | ✓ — **single-session viable** |
| (1, 3, 2)   |    0 | — |
| (2, 1, 3)   |    0 | — |
| (1, 2, 3)   |    0 | — |
| (2, 2, 2)   |    0 | — |

**Headline result.** `(2, 3, 1)` is the **unique** shape with all parts
≤ 3 that admits a BT decomposition. The asymmetry is forced by row 8
(unique support-size-1 row in the inner 6×6) being pinned as the
C-block. The 5×5 left after removing row 8 admits a (2, 3) inner BT
split but no other parts-≤-3 split.

**Refutation of S206 prediction.** S206 hypothesised that "composite W
avoids the W=11-style atomic-block obstruction". W=14 (composite)
shows the obstruction is NOT controlled by primality of W; instead
it is controlled by the distribution of row-support sizes — and W=14
narrowly escapes atomicity only via the unique (2, 3, 1) shape.

The Lean closure was implemented at S235 using this shape; see
`w14_blocktriangular_search_results.md` for the full proof structure.
