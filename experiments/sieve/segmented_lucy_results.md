# Segmented Lucy DP: Cache & Wheel Optimization — Results

**Date:** 2026-04-05 (Session 41)
**Script:** `segmented_lucy.py`

## What Was Tested

Three optimizations applied to the baseline Lucy DP from v10:
1. **Wheel-6** — skip p=2,3 iterations (initialize with coprime-to-6 count)
2. **Wheel-30** — skip p=2,3,5 iterations (initialize with coprime-to-30 count)
3. **Wheel-30 + prefetch** — add __builtin_prefetch hints for cache lines

All verified 100% correct against known pi(x) values (9/9 tests, x up to 10^9).

## Key Findings

### Marginal speedups only

| x | Baseline (ms) | Wheel-6 | Wheel-30 | W30+prefetch | Best speedup |
|---|---|---|---|---|---|
| 10^5 | 0.013 | 0.010 | 0.010 | 0.010 | 1.35x |
| 10^6 | 0.060 | 0.052 | 0.051 | 0.050 | 1.21x |
| 10^7 | 0.302 | 0.275 | 0.266 | 0.271 | 1.14x |
| 10^8 | 1.53 | 1.42 | 1.43 | 1.41 | 1.09x |
| 10^9 | 7.78 | 7.48 | 7.47 | 7.50 | 1.04x |
| 10^10 | 39.5 | 38.1 | 39.7 | 38.9 | 1.04x |
| 10^11 | 199 | 212 | 196 | 197 | 1.02x |

**Maximum speedup: 35% at x=10^5, dropping to 2-4% for x > 10^9.**

### Why optimizations don't help much

1. **Wheel skip saves few iterations:** For p=2, the inner loop has limit ≈ sqrt(x)/4.
   For p=3, limit ≈ sqrt(x)/9. These are small fractions of total work which is
   dominated by the O(sqrt(x)) primes up to sqrt(x).

2. **GCC -O3 already optimizes well:** The compiler generates efficient code for the
   simple inner loop. Prefetch hints don't improve on the hardware prefetcher.

3. **Memory bandwidth is not the bottleneck:** At x=10^9, arrays are ~250KB each,
   fitting in L2/L3 cache. The computation is ALU-bound, not memory-bound.

4. **Diminishing returns at large x:** More primes means more inner-loop iterations
   per prime. The p=2,3,5 savings become negligible relative to total work.

### Why primecount is 100-1000x faster

Wheel factorization and cache prefetching are **within-algorithm** optimizations.
The 100-1000x gap between v10 and primecount comes from:

1. **Different algorithm:** Gourdon decomposes pi(x) into ordinary/easy/hard special
   leaves, each computed with specialized methods. Lucy DP computes everything uniformly.
2. **Segmented sieve:** The hard special leaves use a cache-friendly segmented sieve
   that processes intervals fitting in L1 cache.
3. **SIMD vectorization:** AVX512 POPCNT processes 512 bits of sieve at once.
4. **OpenMP parallelism:** Multiple cores work on independent sieve segments.

**These are algorithmic and architectural changes, not parameter tweaks.**

## Verdict

**CLOSED** as significant optimization path.
**Failure Mode:** E (Equivalence) — within-Lucy-DP optimizations yield at most 35%
speedup, nowhere near the 100-1000x needed to match primecount. The gap is algorithmic.

## Engineering Recommendation

For v10:
- Wheel-30 could be integrated for a free 2-10% improvement, but it's not worth the
  code complexity for the project's purposes.
- Matching primecount requires reimplementing Gourdon's algorithm from scratch — a
  multi-thousand-line C++ effort, essentially forking primecount itself.
- For the project's goals (understanding barriers, not competing on speed), v10 with
  naive Lucy DP is adequate.

## One-Line Summary

Wheel-30 and cache prefetching give only 2-35% speedup on Lucy DP; the 100-1000x gap to primecount requires a completely different algorithm (Gourdon).

**Session:** 41
