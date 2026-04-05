# Batch Mobius Sieve: Recursive Halving for phi(x,a)

**Date:** 2026-04-05  
**Verdict:** CLOSED -- exponential blowup makes this strictly worse than standard methods.

## Idea

Replace the standard Meissel-Lehmer recursion (remove one prime at a time) with recursive halving:

```
phi(x, a) = sum_{d | Q} mu(d) * phi(x/d, a/2)
```

where Q = product of the second-half primes {p_{a/2+1},...,p_a}. This reduces the recursion tree depth from pi(sqrt(x)) to O(log(pi(sqrt(x)))), but at the cost of exponential branching at each level (2^{half_size} squarefree divisors).

## Results

### Correctness
- **x = 10^4:** All three methods (Lucy_Hedgehog, Batch Mobius, Standard recursive) agree: pi(10^4) = 1229. Verified against sympy.

### Performance Comparison

| x | pi(sqrt(x)) | Lucy ops | Lucy time | BMS calls | BMS time | Std calls |
|---|-------------|----------|-----------|-----------|----------|-----------|
| 10^4 | 26 | 994 | 0.000146s | 682 | 0.004610s | 7,959 |
| 10^5 | 66 | 5,002 | 0.000641s | SKIP (2^33) | -- | -- |
| 10^6 | 168 | 24,652 | 0.003585s | SKIP (2^84) | -- | -- |
| 10^7 | 447 | 122,586 | 0.017167s | SKIP (2^224) | -- | -- |

### Branching Factor Explosion

The top-level Mobius sum has 2^{pi(sqrt(x))/2} terms:

| x | pi(sqrt(x)) | Top-level divisors |
|---|-------------|-------------------|
| 10^4 | 26 | 2^13 = 8,192 |
| 10^5 | 66 | 2^33 = 8.6 billion |
| 10^6 | 168 | 2^84 ~ 10^25 |
| 10^7 | 447 | 2^224 ~ 10^67 |

For x = 10^4, the batch method did work (682 phi calls, thanks to heavy pruning -- 10,287 terms pruned where x/d < 1). But it was still 30x slower than Lucy_Hedgehog due to the overhead of generating and iterating over 8,192 divisors at the top level.

### Computation Tree Structure (x = 10^4)

```
Level 0: up to 8,192 divisors (13 primes in second half)
Level 1: up to 128 divisors (7 primes in second half)
Level 2: up to 8 divisors (3 primes in second half)
Level 3: up to 4 divisors (2 primes in second half)
```

Tree depth: 4 (vs. depth 26 for standard recursion). But the breadth explodes.

## Analysis

### Why This Doesn't Work

The fundamental issue: **recursive halving trades depth for exponential breadth, but the total work is unchanged or worse.**

- **Standard recursion:** depth = pi(sqrt(x)), branching = 2 per node, heavily pruned by x/p < 1 bounds and caching. Effective work ~ O(x^{2/3}) with Meissel-Lehmer optimizations.
- **Batch Mobius halving:** depth = O(log(pi(sqrt(x)))), but branching = 2^{half_size} per level. Total leaves = product of branching factors = 2^{pi(sqrt(x))}, same as the unpruned standard tree.

The Mobius sum `sum_{d|Q} mu(d) * phi(x/d, a/2)` is mathematically identical to the inclusion-exclusion that the standard recursion computes -- it just reorganizes the order of operations. No cancellation or compression occurs; the individual terms must still be evaluated.

### Pruning Helps But Not Enough

At x = 10^4, pruning eliminated 10,287 of ~16,000 total terms (64%). But the surviving terms still required 682 phi evaluations, and the divisor generation overhead dominated. As x grows, the number of surviving terms grows exponentially, making the method completely infeasible beyond x ~ 10^4.

### Connection to Known Barriers

This is a specific instance of the general inclusion-exclusion barrier:
- The Legendre sieve phi(x, a) fundamentally requires tracking which integers survive removal of each prime.
- Any reorganization of the inclusion-exclusion (whether one-at-a-time, batch halving, or any other grouping) computes the same combinatorial quantity.
- The O(x^{2/3}) barrier of Meissel-Lehmer comes from clever avoidance of most terms via the P2/P3 decomposition, not from reorganizing the Mobius sum.

## Verdict

**CLOSED.** The recursive halving approach to phi(x, a) offers no computational advantage over standard methods. It trades logarithmic depth for exponential breadth, with the total work remaining 2^{pi(sqrt(x))} in the worst case. Lucy_Hedgehog's O(x^{2/3}) approach is fundamentally different (it works in the "dual" space of floor quotients) and cannot be improved by rearranging the Mobius inversion.
