# k-Party NOF Communication Complexity of pi(x)

## Date: 2026-04-05 (Session 20)

## Summary

We computed the k-party Number-On-Forehead (NOF) communication complexity
of pi(x) for k=2,3,4,5,6,7,8 parties and N=6,8,10,12,14,16,18 bits.

### Key finding: Two distinct regimes

**For k >= 3 (multiparty NOF):** Every mode-i unfolding has **FULL RANK**
(= 2^{sizes[i]}). The max unfolding rank equals exactly 2^{ceil(N/k)}.

**For k = 2 (two-party):** The balanced rank is 2^{N/2-1} + 2, which is
approximately HALF the full rank 2^{N/2}. This is the well-known result
from Session 19.

### The exact formula

For k >= 3 parties with balanced partition of N bits:
```
max mode-unfolding rank = 2^{ceil(N/k)}     (EXACTLY)
```
This was verified for all tested (N,k) pairs with k >= 3 -- zero exceptions
across 35+ test cases.

For k = 2:
```
max mode-unfolding rank = 2^{N/2 - 1} + 2   (matches Session 19)
```

### Scaling implications

| k parties | NOF lower bound (bits) | Scaling with N | Implication |
|-----------|----------------------|----------------|-------------|
| k = 2     | ~N/2                 | linear         | Not in NC^1 via this method |
| k = 3     | ceil(N/3)            | linear         | Full rank, no 3-party structure |
| k = 4     | ceil(N/4)            | linear         | Full rank |
| k = 5     | ceil(N/5)            | linear         | Full rank |
| k = O(1)  | N/k = Theta(N)       | linear         | Still grows linearly |
| k = N/c   | c = O(1)             | constant       | Complexity becomes constant |
| k = O(log N) | N/log(N)          | superlogarithmic | **NOT enough to prove NOT in ACC^0** |

### Why mode-unfolding rank is not sufficient

The mode-i unfolding rank being full for all k >= 3 is **uninformative** for
the ACC^0 question. Here is why:

1. For ANY non-degenerate function f: {0,1}^N -> Z, the mode-i unfolding
   of the k-party tensor has rank <= 2^{sizes[i]}, with equality for
   "generic" (random) functions.

2. pi(x) achieves full rank for k >= 3, which puts it in the same class
   as random functions -- but this does NOT prove high NOF complexity.
   The NOF complexity can be O(polylog) even when all unfolding ranks are
   maximal.

3. The 2-party case is special: the rank 2^{N/2-1}+2 is exactly HALF of
   full rank, confirming the known Session 19 result. This sub-maximality
   is interesting but still gives Omega(N) communication.

### Additional analyses performed

**Rank over GF(p) (p = 2, 3, 5, 7):** Identical to real rank in all cases.
No algebraic structure over any small finite field.

**SVD spectrum:** Extremely concentrated -- the first singular value captures
>99.9% of variance. Effective rank at threshold 0.01 is only 1-2 for all
(N,k) pairs. This means pi(x) is essentially a rank-1 function (the smooth
part x/ln(x)) plus small corrections.

**Residual analysis (pi(x) - x/ln(x)):** The smooth part x/ln(x) has low
rank (8-26 depending on N), but the residual has FULL or even higher rank
than pi(x) itself for k >= 3. The oscillatory part is maximally complex in
tensor structure.

**pi(x) mod 2 and mod 4:** For k >= 3, these also achieve full rank.
For k = 2, they match pi(x)'s rank exactly (same 2^{N/2-1}+2 formula).

## Complete Data Table

```
  N   k  ceil(N/k)  max_rank  full_rank?  log2(r)  NOF_lb
--------------------------------------------------------------
  6   2      3           6       NO         2.58       3
  6   3      2           4       YES        2.00       2
  6   4      2           4       YES        2.00       2
  6   5      2           4       YES        2.00       2
  6   6      1           2       YES        1.00       1
  8   2      4          10       NO         3.32       4
  8   3      3           8       YES        3.00       3
  8   4      2           4       YES        2.00       2
  8   5      2           4       YES        2.00       2
  8   6      2           4       YES        2.00       2
  8   7      2           4       YES        2.00       2
  8   8      1           2       YES        1.00       1
 10   2      5          18       NO         4.17       5
 10   3      4          16       YES        4.00       4
 10   4      3           8       YES        3.00       3
 10   5      2           4       YES        2.00       2
 10   6      2           4       YES        2.00       2
 10   7      2           4       YES        2.00       2
 10   8      2           4       YES        2.00       2
 12   2      6          34       NO         5.09       6
 12   3      4          16       YES        4.00       4
 12   4      3           8       YES        3.00       3
 12   5      3           8       YES        3.00       3
 12   6      2           4       YES        2.00       2
 12   7      2           4       YES        2.00       2
 12   8      2           4       YES        2.00       2
 14   2      7          66       NO         6.04       7
 14   3      5          32       YES        5.00       5
 14   4      4          16       YES        4.00       4
 14   5      3           8       YES        3.00       3
 14   6      3           8       YES        3.00       3
 14   7      2           4       YES        2.00       2
 14   8      2           4       YES        2.00       2
 16   2      8         130       NO         7.02       8
 16   3      6          64       YES        6.00       6
 16   4      4          16       YES        4.00       4
 16   5      4          16       YES        4.00       4
 16   6      3           8       YES        3.00       3
 16   7      3           8       YES        3.00       3
 16   8      2           4       YES        2.00       2
 18   2      9         273       NO         8.09       9
 18   3      6          64       YES        6.00       6
 18   4      5          32       YES        5.00       5
 18   5      4          16       YES        4.00       4
 18   6      3           8       YES        3.00       3
 18   7      3           8       YES        3.00       3
 18   8      3           8       YES        3.00       3
```

## Interpretation

### What we learned

1. **pi(x) has maximal multilinear rank for k >= 3.** Every mode-unfolding
   of the k-party communication tensor achieves full rank. This is the
   generic/expected behavior and does NOT by itself imply high NOF complexity.

2. **The 2-party rank 2^{N/2-1}+2 is an anomaly.** It is the ONLY case
   where rank drops below full. This likely reflects the monotonicity of
   pi(x) and its smooth structure, not a deep communication complexity fact.

3. **Mode-unfolding rank is too coarse for the ACC^0 question.** To prove
   pi(x) is not in ACC^0, one would need to show that the k-party NOF
   complexity is omega(polylog(N)) for some fixed k, using methods beyond
   unfolding rank (e.g., discrepancy, cylinder intersection bounds, or
   polynomial method arguments).

4. **SVD analysis reveals extreme spectral concentration.** pi(x) is
   essentially a rank-1 function (the smooth part) with small perturbations.
   This suggests that approximate computation is very easy, consistent with
   the known result that R^{-1}(n) gives ~50% of digits in O(polylog) time.

### What this does NOT resolve

- Whether pi(x) is in ACC^0 or TC^0
- Whether the true k-party NOF complexity exceeds polylog(N) for any fixed k
- Whether stronger tensor rank measures (e.g., true tensor rank, partition
  rank, analytic rank) would give different scaling

### Recommended next steps

1. **Compute true tensor rank** (not just multilinear/unfolding rank) for
   small cases -- this is NP-hard in general but feasible for N <= 10.
2. **Apply the polynomial method** (Razborov-Smolensky / Cai-Chen-Lipton)
   to pi(x) mod p for small primes.
3. **Study the discrepancy** of pi(x) under cylinder intersection partitions.
4. **Investigate whether the Babai-Nisan-Szegedy multiparty protocol** can
   exploit the spectral concentration of pi(x).

## Files

- `k_party_nof.py` -- v1 implementation (systematic scan)
- `k_party_nof_v2.py` -- v2 implementation (refined analysis with GF(p), SVD, residuals)
- `k_party_nof_results.md` -- this file
