# Results: Persistent Homology of Detrended Prime Sequence (F1)

**Date:** 2026-04-25
**Script:** `persistent_homology_primes.py`
**Conjecture tested:** TCDP — the number of long-lived 0-dim persistence bars
in the detrended prime trajectory is o(N / log N).

## Setup
- 0-dim PH via single-linkage hierarchical clustering on 2D point cloud
  (n/N, normalized d_n) where d_n = p_n - n ln n - n (ln ln n - 1).
- "Long-lived" = bar length > 5x median.
- Compared to two baselines: Cramér model pseudoprimes and IID Gaussian.

## Key numbers

| N    | True primes | Cramér | IID Gaussian | N/log N | long/(N/log N) |
|------|-------------|--------|--------------|---------|----------------|
| 200  | 0           | 0      | 2            | 37.7    | 0.000          |
| 500  | 3           | 6      | 6            | 80.5    | 0.037          |
| 1000 | 11          | 3      | 8            | 144.8   | 0.076          |
| 2000 | 14          | 22     | 10           | 263.1   | 0.053          |

Scaling fit (N = 200, 400, 800, 1600 → long = 0, 2, 5, 8):
- long_lived ≈ 0.005 * N (linear, ~0.5% of N)
- All counts well below N/log N target

## Verdict: WEAKLY POSITIVE BUT INCONCLUSIVE

The long-lived bar count stays at ~0.5-1% of N across N ∈ [200, 2000]. This is
strictly sub-(N/log N) for this range — TCDP holds **as a quantitative observation**.

However:
- Absolute counts are too small (0-14) for reliable scaling conclusions.
- Comparison to Cramér model: true primes are sometimes BELOW (N=500, 2000),
  sometimes ABOVE (N=1000). Cramér model is itself sub-N/log N. This suggests the
  sublinearity is a generic property of bounded-detrended-trajectory PH, not
  specific to primes.

The 5x-median threshold is sensitive: with a stricter threshold (10x), most counts
go to zero; with a looser threshold (3x), all baselines explode together.

## Failure mode

Likely **E (Equivalence)**: the few long-lived bars are not enough to reconstruct
the trajectory uniquely. To recover p_n exactly from PH, one would need both the
*structure* of the persistence diagram AND the assignments of bars to specific
positions n — which is encoded by the ORIGINAL 1D order, i.e., already requires
all of pi(x).

## What would change the verdict

1. **Reconstruction test:** can a "barcode" with k bars reproduce p_n with error
   < 0.5? If yes for k = polylog(N), TCDP is meaningful.
2. **Larger N:** the sublinear scaling at small N might break at N >= 10^4.
3. **Baseline gap:** if true_primes(N) is *consistently* below Cramér(N) by a
   factor that grows with N, true primes have specific extra structure beyond
   the Cramér model.

## Conclusion

The PH approach detects fewer "topological events" than naive O(N) but does not
provide a constructive reconstruction. Without an inverse map (PH-bars -> p_n),
this gives at best a lower bound on the *information content* needed, not a
polylog algorithm. Closing this approach would require either:
- proving no inverse map exists (information-theoretic hard barrier), or
- exhibiting one (algorithmic breakthrough).

Status after test: **interesting topological measurement, no algorithmic
implication**. Not promoted to novel/.
