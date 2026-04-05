# Transfer Matrix Sieve — Results

**Experiment:** Can the sieve of Eratosthenes be encoded as a matrix product with repeated squaring for O(polylog) computation?

**Verdict:** CLOSED — FAILS. Equivalence barrier.

## Key Findings

1. **LCM growth is the fundamental barrier**: The transfer matrix dimension equals lcm(primes up to sqrt(x)). By Chebyshev's theorem, log(lcm) ~ sqrt(x), so the dimension is exp(sqrt(x)) — worse than Meissel-Lehmer's x^{2/3}.

2. **Blocking doesn't help enough**: Grouping primes into manageable blocks still yields total work >> x^{2/3} for all tested scales (10^4 to 10^100). At x=10^100: blocked work ~ 10^100 vs Meissel-Lehmer ~ 10^67.

3. **Euler product approximation converges to Mertens' ratio**: x * prod(1-1/p) / pi(x) → 2*e^{-gamma} ≈ 1.123. This is the classical Mertens theorem — no surprise.

4. **Depth-k sieve reduces dimension but increases correction complexity**: Using only primes up to x^{1/k} for k>2 reduces the matrix dimension but requires handling unsieved composites separately — exactly what Meissel-Lehmer already does.

## Numbers
| x | Primes to sieve | Total work (blocked) | Meissel-Lehmer | Speedup |
|---|---|---|---|---|
| 10^10 | ~9000 | 1.45e+10 | 4.64e+06 | 3.2e-4x (SLOWER) |
| 10^50 | ~2e+23 | 1.44e+50 | 2.15e+33 | 1.5e-17x (SLOWER) |
| 10^100 | ~9e+47 | 1.44e+100 | 4.64e+66 | 3.2e-34x (SLOWER) |

## Failure Mode
**EQUIVALENCE**: The transfer matrix approach is equivalent to the full inclusion-exclusion sieve in matrix form. The multiplicative independence of primes forces exponential state space. No compression is possible because each prime contributes an independent factor to the state.

**Session:** 40 (Fresh Perspective)
