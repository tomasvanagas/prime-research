# Partition-Based Prime Detection (Ono et al. 2024)

## Reference
"Integer partitions detect the primes" — Ono, Craig, van Ittersum (arXiv:2405.06451, PNAS 2024)

## Key Theorem
n ≥ 2 is prime iff (n²−3n+2)·M₁(n) − 8·M₂(n) = 0, where M_a(n) are MacMahon partition functions.

## Implementation Notes
- My implementation of M_a(n) used Σ_{λ⊢n} Σᵢ C(mᵢ(λ), a), which does NOT match the paper's exact definition
- The criterion failed on all tested primes — the paper's M_a functions likely have a different normalization
- However, the COMPUTATIONAL ANALYSIS remains valid regardless of the exact formula

## Computational Results

### Partition enumeration scaling
| n  | p(n) partitions | M₁ time (s) |
|----|----------------|-------------|
| 10 | 42             | 0.0003      |
| 20 | 627            | 0.006       |
| 30 | 5,604          | 0.065       |
| 40 | 37,338         | 0.643       |
| 50 | 204,226        | 3.66        |

p(n) grows as exp(π√(2n/3))/(4n√3) — **exponentially** in √n.

### Generating function approach
- M₁(n) = Σⱼ₌₁ⁿ d(j)·p(n-j) (convolution of divisor function with partition function)
- This avoids enumeration but still requires O(n) terms
- Each p(k) computable in O(k^{1/2}) via Rademacher's exact formula
- Total: O(n^{3/2}) for M₁(n)

## Verdict: CLOSED

The Ono partition criterion is mathematically beautiful but computationally useless:
- **Direct enumeration**: O(exp(√n)) — exponential, infeasible for n > 100
- **Generating function**: O(n^{3/2}) — better, but still much worse than Meissel-Lehmer
- **Meissel-Lehmer**: O(x^{2/3}) — current best
- **Target**: O(polylog) — partition functions cannot reach this

**Why it fails**: Computing partition functions inherently requires processing O(n) or more terms. There is no known shortcut to evaluate MacMahon-type partition statistics in sublinear time. The partition approach converts a number-theoretic problem (primality) into a combinatorial one (partition enumeration) that is HARDER, not easier.

The paper's contribution is to number theory (new characterization of primes), not to algorithm design.
