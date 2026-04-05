# Modular Counting Attack: Results

**Date:** 2026-04-05 (Session 28)
**Experiment:** `experiments/circuit_complexity/modular_counting_attack.py`

## Question

Does decomposing pi(x) by residue class mod M (primorial) give SMALLER circuits
for each per-class counting function pi_r(x)?

This is distinct from CRT approaches (Sessions 3,7,9,13,20-22,24) -- we are asking
about the circuit complexity of the COUNTING problem, not computing p(n) mod q.

## Results Summary

### Experiment 1: Mixed-Radix Conditional Entropy

**Finding:** Conditioning on residue class DOES reduce entropy, but the reduction
**SHRINKS as N grows** (for fixed M).

| N bits | M=6    | M=30   | M=210  | M=2310 |
|--------|--------|--------|--------|--------|
| 10     | 8.0%   | 14.1%  | 50.6%  | --     |
| 12     | 6.3%   | 14.6%  | 24.0%  | --     |
| 14     | 3.0%   | 12.8%  | 26.5%  | 52.2%  |
| 16     | 1.1%   | 7.6%   | 25.9%  | 40.9%  |

**Interpretation:** The entropy reduction comes from two sources:
1. Removing small-prime divisibility (a constant contribution)
2. Knowing the residue class (reduces uncertainty about local prime density)

The first source is O(1) -- it cannot grow with x. The second source also appears
to shrink: as x grows, the within-class distribution approaches the unconditional
distribution. This is consistent with Bombieri-Vinogradov: for small M, the
per-class distribution is asymptotically li(x)/phi(M) with the SAME relative
fluctuations as pi(x) itself.

**Verdict:** Entropy reduction is a FINITE-SIZE EFFECT, not a structural advantage.

### Experiment 2: Per-Class Circuit Complexity

**Finding:** Per-class functions are HARDER, not simpler.

| N bits | M=6 ratio | M=30 ratio | M=210 ratio |
|--------|-----------|------------|-------------|
| 10     | 2.97      | 3.56       | --          |
| 12     | 2.99      | 3.70       | 4.13        |
| 14     | 3.00      | 3.73       | 4.27        |
| 16     | 3.00      | 3.75       | 4.35        |

The transition ratio (per-class normalized transitions / full normalized transitions)
is consistently **3x-4x LARGER** for per-class functions.

**Why:** Full pi(x) mod 2 has strong autocorrelation (~0.80) because consecutive
integers are rarely both prime. But in the per-class function, consecutive q values
map to x values that are M apart -- these are essentially independent for primality.
The per-class function looks like a random binary string with density ~1/ln(x),
while the full function has exploitable regularity from the guaranteed gaps between
primes.

The per-class function has **LOST** the sequential structure that makes pi(x) partially
regular. Each class is a sparse, nearly-random subsequence.

**Verdict:** ANTI-WIN. Decomposition INCREASES per-function complexity.

### Experiment 3: Cross-Class Mutual Information

**Finding:** Classes are nearly independent (I/H < 0.01 for most cases).

As N grows, the mutual information ratio DECREASES:
- N=12, M=210: I/H = 0.042 (mildly dependent)
- N=14, M=210: I/H = 0.009 (near-independent)
- N=16, M=210: I/H = 0.003 (near-independent)

**Interpretation:** This confirms the Bombieri-Vinogradov prediction: for moduli
M << sqrt(x), the per-class counts are nearly independent. The prime distribution
across different coprime residue classes is governed by different Dirichlet L-functions,
and under GRH these are independent.

**However:** Independence means pi(x) = sum_r pi_r(x) is a sum of phi(M) independent
terms, each of which is as complex as computing a random-looking function.
Independence does NOT help -- it means there is no shortcut from one class to another.

### Experiment 4: Divide-and-Conquer Scaling

**Finding:** Total transitions across all classes is LESS than full (ratio 0.95-0.99),
but the "win" is purely from removing composite-class contributions.

The win converges to ratio ~1.0 as N grows:
- N=10, M=30: 0.947
- N=12, M=30: 0.984
- N=14, M=30: 0.995
- N=16, M=30: 0.999

**Interpretation:** The only savings come from not needing to count composites
divisible by small primes -- this is exactly what wheel factorization gives:
a constant factor improvement. As x grows, the savings become negligible relative
to the total cost.

### Experiment 5: Per-Class Entropy Scaling

**Finding:** Total entropy phi(M) * H_per_class GROWS with M.

| M    | phi(M) | Total entropy |
|------|--------|---------------|
| 2    | 1      | 0.72          |
| 6    | 2      | 1.76          |
| 30   | 8      | 7.63          |
| 210  | 48     | 47.39         |
| 2310 | 480    | 471.50        |

Total entropy grows as ~phi(M), while per-class entropy stays ~1 bit.
This is EXACTLY what we'd expect for independent subproblems of roughly equal
difficulty: decomposing does not reduce total work.

## Conclusions

### This approach is CLOSED.

**The decomposition fails for three interconnected reasons:**

1. **Per-class functions are MORE complex, not less.** The sequential regularity
   of pi(x) (consecutive integers, guaranteed gaps between primes) is destroyed
   by the M-spacing of each class. Each class function looks nearly random.

2. **Classes are independent, which is BAD.** Independence means there is no
   shortcut -- we must solve each class separately. The total work is
   phi(M) * (per-class cost), and per-class cost is NOT smaller.

3. **Entropy reduction is a finite-size effect.** The only entropy saved comes
   from removing small-prime divisibility, which is O(1) information regardless
   of x. The remaining per-class entropy is asymptotically identical to the
   full problem's entropy per surviving integer.

**The fundamental issue:** Wheel decomposition trades one hard problem pi(x) for
phi(M) independent hard problems pi_r(x), each of which is individually HARDER
(per input bit) than the original. The only saving is from not computing the
trivially-composite entries, which is the well-known constant-factor improvement
of wheel sieving. There is no complexity-class improvement.

**Connection to known barriers:** This is an instance of the "Information Loss"
failure mode (see `novel/failure_taxonomy.md`). The decomposition preserves total
information content but distributes it across independent subproblems, destroying
the sequential structure that provides the only known computational regularity.

## Status

**Add to CLOSED_PATHS.md:**
| Wheel decomposition circuit complexity | FAIL | I (Info Loss) | Per-class circuits 3-4x MORE complex; independence prevents shortcuts; total entropy grows as phi(M) | 28 |
