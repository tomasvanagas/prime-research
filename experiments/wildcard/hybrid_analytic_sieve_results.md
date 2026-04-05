# Hybrid Analytic-Sieve Method for p(n)

## Concept
Combine analytic approximation (R⁻¹(n) + zero corrections) with local sieving to find p(n). Only need π(x) precise enough to locate p(n) in a small interval, then sieve locally.

## Key Results

### R⁻¹(n) approximation quality
| n | p(n) | |R⁻¹(n) - p(n)| | |error|/√p(n) |
|---|------|-----------------|---------------|
| 100 | 541 | 4.10 | 0.176 |
| 1000 | 7919 | 4.05 | 0.046 |
| 10000 | 104729 | 39.31 | 0.122 |
| 50000 | 611953 | 165.24 | 0.211 |

Error scales as O(√p(n)) with bounded ratio 0.05-0.25.

### Zero corrections make things WORSE for partial sums!
For n=50000 (p(n)=611953):
- K=0 (R function only): |error in π̂| = 12.4
- K=100 zeros: |error| = 205.6
- K=500 zeros: |error| = 207.6

The truncated zero sum x^ρ/ρ at 1000 zeros introduces MORE error than it removes, because:
1. The series converges conditionally (not absolutely)
2. The simplified x^ρ/ρ doesn't match the full R(x^ρ) Gram series
3. With 1000 zeros, γ_max ≈ 1419, but for x=611953, we need γ_max ~ √x ≈ 782

### Iterative refinement
- Converges in 3-7 iterations to a fixed point
- BUT the fixed point has error ~2000 for n=50000 (worse than initial estimate!)
- Bottleneck: π(lo) estimation needs O(√x) zeros for O(1) accuracy

### Complexity analysis
| Method | Cost | p(10^100) ops |
|--------|------|---------------|
| Eratosthenes | O(x log log x) | 10^102 |
| Meissel-Lehmer | O(x^{2/3}) | 10^66 |
| Lagarias-Odlyzko | O(x^{1/2+ε}) | 10^50 |
| Hybrid (optimal) | O(x^{1/2+ε}) | 10^50 |
| **Polylog target** | O(log^k x) | **10^2 - 10^3** |

The hybrid method matches but cannot beat Lagarias-Odlyzko.

## Why Polylog is Blocked
For polylog(x) zeros: error tail ~ √x/polylog(x) → unbounded.
The error is dominated by the LOWEST unfixed zero, which contributes O(x^{1/2}/γ_1). Even one zero contributes O(√x/14.13), which for x=10^100 is O(10^{49}).

## Verdict
DEAD END for polylog. The hybrid approach cannot circumvent the O(√x) information requirement. The zero sum's conditional convergence makes partial sums unreliable — sometimes NO zeros is better than some zeros.
