# Schoenfeld-Cramer Experiments: Results

**Task #4, Experiments 5 & 6**
**Date:** 2026-04-05
**Script:** `schoenfeld_cramer.py`

## Experiment 5: Schoenfeld's Explicit Bounds under RH

### Part A: Schoenfeld Bound Size

Under RH, |pi(x) - li(x)| < sqrt(x)*ln(x)/(8*pi) for x >= 2657 (Schoenfeld 1976).

| x | Schoenfeld bound | x^(2/3) (ML) | Ratio |
|---|---|---|---|
| 10^4 | 36.7 | 464 | 0.079 |
| 10^8 | 7,329 | 215,443 | 0.034 |
| 10^12 | 1.1M | 100M | 0.011 |
| 10^16 | 147M | 46.4B | 0.003 |
| 10^20 | 18.3B | 21.5T | 0.0009 |

The Schoenfeld bound is always smaller than x^(2/3), and the ratio shrinks as x grows. However, the sieve cost over the Schoenfeld interval is O(sqrt(x)*ln(x)*ln(ln(x))), which is still O(x^{1/2+epsilon}) -- better asymptotically than O(x^{2/3}) but the constant factors and the need for RH make this impractical.

### Part B: Sieve Cost Comparison

| x | Interval width 2E | Sieve cost | ML cost O(x^{2/3}) | Ratio |
|---|---|---|---|---|
| 10^4 | 73 | 163 | 464 | 0.35 |
| 10^8 | 14,700 | 42,700 | 215,443 | 0.20 |
| 10^12 | 2.2M | 7.3M | 100M | 0.073 |
| 10^20 | ~10^10.3 | ~10^11.1 | ~10^13.3 | 0.0065 |
| 10^50 | ~10^26 | ~10^26.6 | ~10^33.3 | 2e-7 |
| 10^100 | ~10^51.3 | ~10^52 | ~10^66.7 | 2e-15 |

**Finding:** Sieve-in-Schoenfeld-interval is O(x^{1/2+epsilon}), technically better than Meissel-Lehmer O(x^{2/3}). But this requires RH and still nowhere near polylog.

### Part C: Error Reduction with K Zeta Zeros

Tested correction using explicit formula with K zeros from our 1000-zero dataset.

**x = 10^6 (pi(x) = 78,498):**

| K zeros | Corrected estimate | Actual error | Heuristic bound |
|---|---|---|---|
| 0 (li only) | 78,627.5 | 129.5 | -- |
| 1 | 78,622.3 | 124.3 | 95,434 |
| 10 | 78,616.6 | 118.6 | 17,352 |
| 100 | 78,601.5 | 103.5 | 1,890 |
| 500 | 78,602.4 | 104.4 | 381 |
| 1000 | 78,600.4 | 102.4 | 191 |

**Key observation:** The correction converges VERY slowly. Even with 1000 zeros, the error only drops from 129.5 to 102.4 -- a mere 21% improvement. The oscillatory terms from higher zeros partially cancel but leave a substantial residual.

**x = 10^8 (pi(x) = 5,761,455):**

| K zeros | Actual error |
|---|---|
| 0 | 754.4 |
| 100 | 612.4 |
| 1000 | 638.3 |

At x = 10^8, the convergence is even worse -- the error actually fluctuates rather than monotonically decreasing, because higher zeros have quasi-random phases.

### Part D: Zeros Needed for Error < 1

| x | K needed | Comment |
|---|---|---|
| 10^6 | ~10^5.3 | 200,000 zeros |
| 10^10 | ~10^7.7 | 50 million zeros |
| 10^20 | ~10^13.3 | 20 trillion zeros |
| 10^50 | ~10^29.1 | Absurd |
| 10^100 | ~10^54.7 | Utterly impossible |

**Verdict:** To achieve error < 1 (needed for exact p(n)), the number of zeros required scales as O(sqrt(x) * log^2(x)), which is worse than Meissel-Lehmer.

### Part E: Convergence Profile at x = 10^6

The correction from K = 1 to K = 1000 zeros showed:
- K=1: error 124.3 (from 129.5 baseline)
- K=10: error 118.6
- K=100: error 103.5
- K=500: error 104.4
- K=1000: error 102.4

The convergence is logarithmic at best, with significant oscillation. The error does not approach zero systematically because we need O(sqrt(x)*log^2(x)) ~ 200,000 zeros for x = 10^6.

## Experiment 6: Cramer's Conjecture Search Algorithm

### Part A: Algorithm Implementation and Phase Costs

Implemented the full Cramer search: (1) R^{-1}(n) -> x0, (2) pi(x0), (3-4) walk to p(n).

| n | p(n) | Offset x0-p(n) | Walk steps | t_count | t_walk | %counting |
|---|---|---|---|---|---|---|
| 100 | 541 | -5.1 | 1 | 0.0001s | 0.0001s | 0.2% |
| 1,000 | 7,919 | +2.9 | 0 | 0.0002s | 0.004s | 0.4% |
| 10,000 | 104,729 | +38.0 | 3 | 0.0003s | 0.005s | 0.9% |
| 100,000 | 1,299,709 | +23.6 | 1 | 0.0015s | 0.007s | 3.4% |
| 1,000,000 | 15,485,863 | -1,824 | 114 | 0.0075s | 0.017s | 12.8% |
| 5,000,000 | 86,028,121 | -1,223 | 65 | 0.024s | 0.037s | 24.7% |

All results verified correct. The counting phase (pi(x0)) grows as the dominant cost as n increases, reaching 25% at n=5M and will dominate for larger n.

### Part B: Gap Statistics (Cramer Verification)

| n | p(n) | ln^2(p) | Max gap (20 neighbors) | Max/ln^2 | Avg gap | Avg/ln^2 |
|---|---|---|---|---|---|---|
| 1,000 | 7,919 | 80.6 | 30 | 0.37 | 9.6 | 0.12 |
| 10,000 | 104,729 | 133.6 | 24 | 0.18 | 11.2 | 0.08 |
| 100,000 | 1,299,709 | 198.2 | 28 | 0.14 | 13.5 | 0.07 |
| 1,000,000 | 15,485,863 | 274.1 | 50 | 0.18 | 17.3 | 0.06 |

Gaps are well within the Cramer bound (max ratio 0.37, well under Granville's C ~ 1.12). This confirms the walk phase is extremely cheap.

### Part C: Theoretical Bottleneck Analysis

| x | Count phase O(x^{2/3}) | Walk phase O(ln^4 x) | Ratio |
|---|---|---|---|
| 10^8 | 10^5.3 | 10^5.1 | ~2 |
| 10^12 | 10^8.0 | 10^5.8 | ~160 |
| 10^20 | 10^13.3 | 10^6.7 | 10^6.7 |
| 10^50 | 10^33.3 | 10^8.2 | 10^25.1 |
| 10^100 | 10^66.7 | 10^9.5 | 10^57.2 |

At x = 10^100, the counting phase is 10^57 times more expensive than the walk phase.

## Summary and Verdict

### Experiment 5 (Schoenfeld): PATH CLOSED
- Schoenfeld bounds give a sieve-based algorithm of complexity O(x^{1/2+epsilon}), which is better than Meissel-Lehmer O(x^{2/3}) asymptotically
- BUT still exponentially far from polylog
- Iterative refinement with zeta zeros converges too slowly: need O(sqrt(x)*log^2(x)) zeros for exactness
- Numerically verified: 1000 zeros barely dent the error at x = 10^6 (21% reduction)
- The information bottleneck is clear: each zero contributes O(1/K) correction, but there are O(sqrt(x)) zeros that matter

### Experiment 6 (Cramer): PATH CLOSED
- Cramer's conjecture makes the SEARCH phase trivial: O(ln^4 x)
- But the COUNTING phase (computing pi(x0)) remains the bottleneck at O(x^{2/3})
- At x = 10^100, counting is 10^57x more expensive than searching
- Key insight: **Cramer solves the wrong problem** -- the hard part of computing p(n) is not finding a prime near x0, it's determining which prime is the nth one

### Combined Insight
Both approaches reveal the same fundamental structure:
1. **Approximation is easy:** R^{-1}(n) gives ~50% of digits in polylog time
2. **Locating primes is easy:** Under Cramer, once you know where to look, gaps are O(ln^2 x)
3. **Counting primes is hard:** Computing pi(x) exactly is the irreducible bottleneck
4. The counting bottleneck connects to the information-theoretic barrier: the exact value of pi(x) encodes information from O(sqrt(x)) zeta zeros with GUE-random phases
