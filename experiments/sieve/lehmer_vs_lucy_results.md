# Lehmer vs Lucy DP Benchmark — Results

**Date:** 2026-04-05 (Session 41)
**Script:** `lehmer_vs_lucy.py`

## What Was Tested

Three C implementations of pi(x) benchmarked head-to-head:
1. **Lucy DP** — the naive floor-value DP from v10 (baseline)
2. **Lehmer** — pi(x) = phi(x,a) + a - 1 - P2 where a = pi(x^{1/3})
3. **Meissel-Lehmer** — same but with a = pi(x^{1/4}) and explicit P3 term

All three verified 100% correct: exhaustive on [2, 10000], spot-checks to 10^10.

## Key Findings

### Lucy DP wins decisively

| x | Lucy (ms) | Lehmer (ms) | M-L (ms) | Lucy advantage |
|---|---|---|---|---|
| 10^6 | 0.063 | 0.190 | 0.233 | 3x faster |
| 10^8 | 1.53 | 8.25 | 8.79 | 5x faster |
| 10^9 | 7.80 | 62.1 | 56.7 | 7-8x faster |
| 10^10 | 31.6 | 405.8 | 303.0 | 10-13x faster |

**Lucy DP is 3-13x FASTER than both Lehmer variants.** The gap widens at larger x.

### Why Lucy DP wins

1. **Cache friendliness:** Lucy DP iterates sequentially over two arrays of size O(sqrt(x)). 
   Near-perfect L1/L2 cache utilization.

2. **No recursion:** Lehmer's phi function uses deep inclusion-exclusion recursion. 
   The 16K hash cache has ~30% collision rate, and the recursive call tree 
   for phi(x, a) with a = pi(x^{1/3}) has depth proportional to a ~ x^{1/3}/ln(x).

3. **Simple inner loop:** Lucy DP's inner loop is just array index arithmetic and subtraction. 
   The Lehmer P2 loop calls lucy_pi() for each semiprime, adding function call overhead.

4. **Theoretical advantage is constant-factor, not asymptotic for practical x:**
   Lehmer is O(x^{2/3}/log x) vs O(x^{2/3}), but the log x improvement 
   (factor ~7 at x=10^10) is overwhelmed by the ~50x worse constant.

### Why primecount is still faster than both

primecount (Gourdon variant) achieves its 100-1000x speedup over our v10 NOT by 
using Lehmer-style decomposition with recursive phi. Instead:

1. **Segmented sieve for phi:** Replaces recursive inclusion-exclusion with a 
   cache-friendly segmented sieve, eliminating recursion entirely.
2. **Special leaves decomposition:** Splits P2 into ordinary, easy, and hard 
   special leaves with different algorithms for each.
3. **SIMD/AVX512:** Vectorized bit counting for sieve operations.
4. **OpenMP:** Parallel computation across cores.

The lesson: **the decomposition (Lehmer/Gourdon) only helps when phi is computed 
via segmented sieve, not via recursion.** A recursive phi implementation is worse 
than Lucy DP.

## Verdict

**CLOSED** as an improvement path for v10.
**Failure Mode:** E (Equivalence) — naive Lehmer is slower than Lucy DP; the 
Deleglise-Rivat/Gourdon speedup comes from segmented sieve infrastructure, 
not from the Lehmer decomposition alone.

## Implication for v10 Engineering

To match primecount, v10 would need:
1. Segmented sieve (replacing recursive phi) — ~10x improvement
2. Special leaves decomposition — ~5x improvement  
3. SIMD + parallelism — ~10x improvement
4. Total potential: ~500x, matching the observed gap

This is essentially reimplementing primecount. For the project's purposes, 
v10 with Lucy DP (p(10^9) in 7.8ms for pi(x), 175ms for full p(n)) is 
adequate. The barrier is asymptotic (O(x^{2/3}) for any variant), not 
constant-factor.

## One-Line Summary

Naive Lehmer is 5-13x SLOWER than Lucy DP due to recursive phi; primecount's speed comes from segmented sieve infrastructure, not decomposition alone.

**Session:** 41
