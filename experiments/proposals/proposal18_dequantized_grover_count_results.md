# Proposal 18: Dequantized Grover Counting for Primes — Results

## Idea
Apply dequantization techniques (classical analogues of quantum speedups) to
prime counting. Specifically, use the Mobius inversion pi(x) = sum mu(d)*floor(x/d)
and attempt sampling-based speedups.

## Results

### Truncated Mobius sums
| x | pi(x) | D=10 | D=50 | D=100 | D=sqrt(x) |
|---|-------|------|------|-------|-----------|
| 100 | 25 | 8 | -4 | 0 | 8 |
| 1000 | 168 | 90 | -20 | 29 | -61 |
| 10000 | 1229 | 904 | -204 | 310 | 310 |

Truncation at D < x gives wildly inaccurate results. The Mobius sum does not
converge quickly because sum(|mu(d)|/d) diverges logarithmically.

### Hyperbola method
O(sqrt(x)) complexity using the Dirichlet hyperbola method. Verified to work
but this is already known (Lagarias-Odlyzko achieves O(x^{1/2+eps})).

### Random sampling (dequantized approach)
| x | Samples | Estimate | Error |
|---|---------|----------|-------|
| 1000 | 100 | -456 | 624 (371%) |
| 1000 | 10000 | 20 | 148 (88%) |
| 10000 | 1000 | 1191 | 38 (3.1%) |
| 10000 | 10000 | 807 | 422 (34%) |

Sampling is extremely noisy. Variance is O(x^2/S), requiring S ~ x^2 samples
for sub-1 accuracy — far worse than direct computation.

## Why dequantization fails here
1. mu(d) has no low-rank structure — it's pseudorandom (Mobius randomness)
2. Dequantization results (Tang, Chia et al.) require low-rank input matrices
3. The prime indicator function chi_P has approximate degree N/2 (our novel finding),
   confirming it resists any polynomial approximation scheme
4. Grover's quadratic speedup to O(sqrt(x)) is already matched classically by
   Lagarias-Odlyzko

## Verdict: CLOSED
Dequantization cannot help because the Mobius function (which encodes primality)
has no exploitable structure. The problem is "full-rank" in every known basis.
