# Explicit Circuit Synthesis for pi(x) -- Results

## Experiment (Session 28)

For N-bit inputs x in [0, 2^N - 1], compute pi(x) and analyze Boolean circuit complexity
via BDD (Binary Decision Diagram) construction with multiple variable orderings.

**This is the first explicit BDD/circuit measurement for pi(x) in this project.**

## Key Finding

**BDD(LSB of pi) ~ 2^(0.73*N).** This is WORSE than sqrt(x) = 2^(0.5*N), approaching
the trivial 2^N bound. BDDs are more restricted than general circuits, so this is an
upper bound on BDD complexity and does NOT directly bound circuit complexity.

However, combined with OBDD ~ 2^(0.79*N) from Session 20, the evidence suggests that
branching-program representations of pi(x) are genuinely hard -- harder than sqrt(x).

## Summary Table: BDD Nodes per Output Bit

| N | pi(2^N-1) | Output bits | LSB BDD | log2(LSB BDD) | Total BDD | BDD/N^2 | BDD/N^3 | BDD/2^(N/2) |
|---|-----------|-------------|---------|---------------|-----------|---------|---------|-------------|
| 4 | 6 | 3 | 7 | 2.81 | 18 | 0.44 | 0.109 | 1.75 |
| 5 | 11 | 4 | 13 | 3.70 | 35 | 0.52 | 0.104 | 2.30 |
| 6 | 18 | 5 | 19 | 4.25 | 62 | 0.53 | 0.088 | 2.38 |
| 7 | 31 | 5 | 35 | 5.13 | 108 | 0.71 | 0.102 | 3.09 |
| 8 | 54 | 6 | 53 | 5.73 | 196 | 0.83 | 0.104 | 3.31 |
| 9 | 97 | 7 | 82 | 6.36 | 312 | 1.01 | 0.113 | 3.62 |
| 10 | 172 | 8 | 147 | 7.20 | 534 | 1.47 | 0.147 | 4.59 |
| 11 | 309 | 9 | 247 | 7.95 | 876 | 2.04 | 0.186 | 5.46 |
| 12 | 564 | 10 | 421 | 8.72 | 1436 | 2.92 | 0.244 | 6.58 |
| 13 | 1028 | 11 | 698 | 9.45 | 2402 | 4.13 | 0.318 | 7.71 |
| 14 | 1900 | 11 | 1207 | 10.24 | 3941 | 6.16 | 0.440 | 9.43 |

## Per-Bit BDD Sizes (all output bits)

Higher bits (MSB) have smaller BDDs; LSB is always hardest.

```
N=4:  [7, 7, 4]
N=5:  [13, 9, 8, 5]
N=6:  [19, 14, 12, 11, 6]
N=7:  [35, 28, 23, 15, 7]
N=8:  [53, 45, 39, 31, 20, 8]
N=9:  [82, 67, 55, 48, 32, 19, 9]
N=10: [147, 110, 84, 74, 54, 36, 19, 10]
N=11: [247, 180, 132, 105, 85, 60, 35, 21, 11]
N=12: [421, 292, 206, 159, 126, 97, 61, 40, 22, 12]
N=13: [698, 484, 331, 248, 215, 156, 115, 72, 45, 25, 13]
N=14: [1207, 803, 549, 384, 318, 251, 188, 117, 73, 37, 14]
```

**Pattern:** MSB BDD ~ N+1 (trivially small, just encodes magnitude).
LSB BDD ~ 2^(0.73*N). Middle bits interpolate.

## Scaling Analysis

### Exponential Fit

Fitting log2(LSB BDD) = beta * N + gamma:

- **beta = 0.7332** (exponent)
- **gamma = -0.10** (constant)

This means: **BDD(LSB) ~ 2^(0.73*N) = x^0.73**

Comparison:
- Trivial bound: 2^N (beta = 1.0)
- OBDD (Session 20): 2^(0.79*N) (beta = 0.79)
- **BDD multi-ordering (this): 2^(0.73*N) (beta = 0.73)**
- sqrt(x) barrier: 2^(N/2) (beta = 0.5)
- Polynomial: N^k (beta ~ 0)

### Power Law Fit

Fitting BDD = c * N^alpha:

- **alpha = 4.07** (exponent)
- **c = 0.02** (constant)

So BDD ~ 0.02 * N^4. But this fit is misleading: the BDD/N^3 ratio is GROWING
(0.109 -> 0.440), so N^4 is not stable either. The exponential model is a better fit.

### Key Ratios

| N | BDD/2^(0.8N) | BDD/2^(0.73N) | BDD/N^3 |
|---|-------------|---------------|---------|
| 4 | 0.762 | 0.917 | 0.109 |
| 6 | 0.682 | 0.808 | 0.088 |
| 8 | 0.628 | 0.739 | 0.104 |
| 10 | 0.574 | 0.717 | 0.147 |
| 12 | 0.542 | 0.669 | 0.244 |
| 14 | 0.513 | 0.631 | 0.440 |

- BDD/2^(0.73N) is slowly decreasing -- the fit slightly overestimates
- BDD/N^3 is clearly growing -- polynomial model breaks down
- True scaling is likely between 2^(0.65N) and 2^(0.75N)

## Influence Analysis (Threshold Gate Proxy)

Total influence of LSB(pi(x)) over N-bit inputs:

| N | Total influence | Influence/N |
|---|----------------|-------------|
| 4 | 2.25 | 0.56 |
| 6 | 3.12 | 0.52 |
| 8 | 3.91 | 0.49 |
| 10 | 4.88 | 0.49 |
| 12 | 5.54 | 0.46 |
| 14 | 6.31 | 0.45 |

**Influence grows as ~0.5*N.** This is moderate -- parity has influence = N,
and random functions have influence ~N/2. The LSB of pi(x) has influence ~N/2,
consistent with pseudorandom behavior.

A single threshold gate can realize functions with influence <= N. Since
influence(LSB) ~ N/2 < N, a single threshold gate MIGHT suffice for each N
(though this is only a necessary condition, not sufficient).

## LSB Minimum Circuit (Exhaustive, N=4)

For N=4, the minimum AND/OR/XOR circuit (with free NOT) for pi(x) mod 2 has
**3 gates**. For N >= 5, the exhaustive search exceeded memory limits.

## Comparison to Previous Models

| Model | LSB Scaling | Exponent beta | Source |
|-------|------------|---------------|--------|
| ANF (GF2) | Full degree N, ~50% coeffs | N/A (dense) | Session 13 |
| OBDD (fixed order) | 2^(0.79*N) | 0.79 | Session 20 |
| BDD (best of k orders) | 2^(0.73*N) | 0.73 | **This experiment** |
| Communication rank | 2^(N/2-1)+2 | 0.50 | Session 17 |

**Hierarchy:** Comm rank (0.50) < BDD (0.73) < OBDD (0.79) < trivial (1.0).

The BDD with multiple orderings is better than single-order OBDD (0.73 vs 0.79),
as expected. The gap is small, suggesting variable ordering optimization helps
somewhat but doesn't fundamentally change the scaling.

The communication rank (beta = 0.50) corresponds to the sqrt barrier and is a
*lower bound* on communication complexity. The BDD is an *upper bound* on BDD
complexity. These are different models measuring different things.

## Implications for Circuit Complexity

1. **BDD complexity does NOT equal circuit complexity.** There exist functions with
   exponential BDD size but polynomial circuit size (e.g., integer multiplication).
   So BDD ~ 2^(0.73N) does NOT prove circuit complexity is super-polynomial.

2. **The BDD gives a circuit UPPER BOUND.** Any BDD of size S can be converted to a
   circuit of size O(S). So the circuit for LSB(pi(x)) is at most ~2^(0.73N).
   But this is not useful -- we already know O(x^{2/3}) algorithms exist.

3. **The interesting question remains open:** does pi(x) have polynomial-size general
   Boolean circuits? The BDD evidence says "BDDs can't do it in poly size" but
   doesn't answer whether unrestricted circuits can.

4. **Next steps would be:** SAT-based exact circuit minimization for N <= 10, or
   ABC/AIGER synthesis tools for larger N. These could find circuits much smaller
   than BDDs.
