# Elliott-Halberstam Conjecture & Gap Structure Experiments

**Task #4, Experiments 3 & 4**
**Date:** 2026-04-05

## Experiment 3: Elliott-Halberstam Conjecture and Counting

### Setup
- Tested primorial moduli q = 6, 30, 210
- For each q, computed pi(x; q, a) for all coprime residues a
- Summed to reconstruct pi(x) and compared with li(x)/phi(q) approximation
- Tested x = 10^3, 10^4, 10^5, 10^6

### Results

**Reconstruction is exact (tautology):** Summing pi(x;q,a) over all coprime residues + primes dividing q always gives pi(x) exactly. This is just partitioning primes by residue class -- no information gained.

**Error cancellation in the approximation:**

| q | x | sum of errors (li/phi approx) | RMS per-class error | cancellation ratio |
|---|---|---:|---:|---:|
| 6 | 10^3 | -11.6 | 6.5 | 0.89 |
| 6 | 10^6 | -131.6 | 67.9 | 0.97 |
| 30 | 10^3 | -12.6 | 2.5 | 0.63 |
| 30 | 10^6 | -132.6 | 21.8 | 0.76 |
| 210 | 10^3 | -13.6 | 1.0 | 0.30 |
| 210 | 10^6 | -133.6 | 13.2 | 0.21 |

The cancellation ratio = |sum of errors| / (phi(q) * RMS). Lower means MORE cancellation. With larger q (210), individual errors are smaller but they do NOT cancel well enough -- the sum of errors converges to li(x) - pi(x), which is the same error as using li(x) directly.

**Timing: residue-class summation is vastly slower:**

| x | Direct primepi | Residue sum (q=30) | Slowdown |
|---|---:|---:|---:|
| 10^4 | 0.00003s | 0.02s | ~1800x |
| 10^5 | 0.000004s | 0.36s | ~90000x |
| 10^6 | 0.000004s | 5.5s | ~1.3M x |

### Verdict

**EH does not help with counting.** The Elliott-Halberstam conjecture controls the error in pi(x;q,a) - li(x)/phi(q) for individual residue classes. When you sum over all coprime residue classes:
- The exact sum gives pi(x) (tautology, no speedup)
- The approximate sum (using li(x)/phi(q) for each class) gives li(x), which has the same error as computing li(x) directly
- EH provides information about **distribution** of primes across residue classes, not about **total count**
- The residue-class approach is orders of magnitude slower than direct computation

**Path closed.** EH is fundamentally about equidistribution, not enumeration.

---

## Experiment 4A: Cramer's Model for Prime Gaps

### Setup
- Generated all primes up to 10^7 (664,579 primes)
- Binned by decade, computed gap statistics
- Tested Cramer's model: normalized gaps g_n/ln(p_n) should be ~Exp(1)

### Results

| Range | Mean gap | Expected (ln p) | Max gap | Cramer bound (ln^2 p) | Max/Cramer | Norm. mean | Norm. std |
|-------|---:|---:|---:|---:|---:|---:|---:|
| [10, 100) | 4.29 | 3.73 | 8 | 21.2 | 0.38 | 1.15 | 0.44 |
| [100, 1K) | 6.35 | 6.10 | 20 | 47.7 | 0.42 | 1.04 | 0.60 |
| [1K, 10K) | 8.48 | 8.43 | 36 | 84.8 | 0.42 | 1.01 | 0.72 |
| [10K, 100K) | 10.76 | 10.74 | 72 | 132.5 | 0.54 | 1.00 | 0.76 |
| [100K, 1M) | 13.06 | 13.04 | 114 | 190.9 | 0.60 | 1.00 | 0.80 |
| [1M, 10M) | 15.36 | 15.35 | 154 | 259.8 | 0.59 | 1.00 | 0.83 |

Cramer's model fits well:
- Normalized mean gap converges to 1.0 (Exp(1) mean)
- Normalized std converges toward 1.0 (Exp(1) std)
- Max gaps stay well below Cramer's ln^2(p) bound (ratio 0.4--0.6)

---

## Experiment 4B: Search Cost with R^{-1}(n)

### Setup
- Computed R^{-1}(n) via Newton's method on li(x) = n
- Search interval: [R^{-1}(n) - 2*ln^2(R^{-1}(n)), R^{-1}(n) + 2*ln^2(R^{-1}(n))]
- Counted primes in interval, checked if p(n) is present

### Results

| n | p(n) | R^{-1}(n) | Error | Rel. error | Interval size | Primes in interval | p(n) found? |
|---:|---:|---:|---:|---:|---:|---:|:---:|
| 10 | 29 | 20.2 | -8.8 | -3.0e-1 | 37 | 12 | Yes |
| 100 | 541 | 488.7 | -52.3 | -9.7e-2 | 154 | 23 | Yes |
| 1,000 | 7,919 | 7,763.0 | -156.0 | -2.0e-2 | 322 | 35 | Yes |
| 10,000 | 104,729 | 104,269.8 | -459.2 | -4.4e-3 | 535 | 45 | No |
| 50,000 | 611,953 | 611,012.8 | -940.2 | -1.5e-3 | 711 | 55 | No |
| 100,000 | 1,299,709 | 1,298,172.4 | -1,536.6 | -1.2e-3 | 793 | 43 | No |
| 500,000 | 7,368,787 | 7,365,837.3 | -2,949.7 | -4.0e-4 | 1,001 | 63 | No |

Key observations:
- R^{-1}(n) underestimates p(n) systematically (error is negative)
- The error |R^{-1}(n) - p(n)| grows faster than C*ln^2 for larger n
- For n >= 10,000, p(n) falls OUTSIDE the 2*ln^2 search interval
- The search interval contains O(ln(x)) primes as expected
- A larger safety factor (or better R^{-1} approximation) would be needed

**Even if p(n) is in the interval:** we still need to know pi(boundary) to identify which prime is the nth one. This is the counting bottleneck.

---

## Experiment 4C: The Counting Bottleneck

### Cost Comparison

| n | p(n) | Sequential O(n) | Meissel-Lehmer O(x^{2/3}) | EH optimistic O(x^{1/2} ln) | Polylog O(ln^3) |
|---:|---:|---:|---:|---:|---:|
| 1,000 | 7,919 | 1,000 | 397 | 799 | 723 |
| 10,000 | 104,729 | 10,000 | 2,222 | 3,741 | 1,545 |
| 100,000 | 1,299,709 | 100,000 | 11,910 | 16,049 | 2,790 |

Key insight: Even under the most optimistic reading of EH, the counting cost is O(x^{1/2+eps}), which is WORSE than Meissel-Lehmer's O(x^{2/3}) at these scales but better asymptotically. Neither comes close to polylog.

---

## Overall Verdict

**Both paths closed for polylog p(n):**

1. **Elliott-Halberstam does not help counting.** EH is about equidistribution of primes in residue classes. Summing the approximations recovers li(x), not pi(x). The distribution information does not translate to a faster way to compute the total count.

2. **Gap structure does not bypass counting.** Even with:
   - Perfect knowledge of gap distribution (Cramer's model)
   - R^{-1}(n) giving O(ln^2(x))-sized search interval
   - Sieving that interval cheaply
   
   ...we still need pi(x) at the interval boundary to determine WHICH prime is p(n). Computing pi(x) exactly costs O(x^{2/3}) (Meissel-Lehmer) or O(x^{1/2+eps}) (analytic), and no gap conjecture reduces this.

**The counting bottleneck is fundamental:** To find p(n), you must either count primes up to some point (O(x^{2/3}) at best) or accumulate all gaps from p(1) (O(n) = O(x/ln x)). Neither is polylog. The gap between these approaches and polylog is enormous and no known conjecture bridges it.
