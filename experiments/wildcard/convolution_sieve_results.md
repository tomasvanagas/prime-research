# Convolution Sieve Experiment Results

## Idea

Represent the Sieve of Eratosthenes as a polynomial product
`P(z) = prod_{p <= sqrt(x)} (1 - z^p)` and investigate whether FFT-based
polynomial arithmetic or evaluation at roots of unity can compute `pi(x)`
faster than direct methods.

## Key Discovery: Structural Mismatch

**The polynomial product P(z) does NOT encode the Legendre sieve.**

- P(z) has nonzero coefficients at **subset-SUM** positions of primes (additive structure)
- The Legendre sieve formula needs **subset-PRODUCT** divisors of the primorial (multiplicative structure)
- FFT operates on power series (additive exponents) and cannot accelerate Dirichlet convolution (multiplicative structure)

This is the fundamental reason the approach fails.

## Experiment Results

### Experiment 1: Coefficient Structure of P(z)

| x | B | #primes | Nonzero coeffs | Subset sums <= x | All at subset-sum positions? | Max |coeff| |
|---|---|---------|----------------|-------------------|------------------------------|-------------|
| 100 | 10 | 4 | 8 | 12 | Yes | 1 |
| 500 | 22 | 8 | 32 | 72 | Yes | 1 |
| 1000 | 31 | 11 | 92 | 155 | Yes | 2 |

- Nonzero coefficients occur at positions that are subset sums of sieving primes
- Some subset sums have zero coefficient due to cancellation (different subsets of same parity summing to same value)
- P(1) = 0 always (since factor (1 - z^2) evaluates to 0 at z=1)
- For a=11 primes, max |coeff| = 2, showing cancellations already occur

### Experiment 2: Legendre Sieve Correctness

| x | true pi(x) | Legendre pi(x) | Match | t_Legendre | t_sieve | Ratio |
|---|------------|----------------|-------|------------|---------|-------|
| 100 | 25 | 25 | OK | 0.000029s | 0.000011s | 2.5x |
| 500 | 95 | 95 | OK | 0.000171s | 0.000012s | 13.8x |
| 1000 | 168 | 168 | OK | 0.001862s | 0.000020s | 94x |
| 5000 | 669 | 669 | OK | 0.739943s | 0.000155s | 4782x |

- Legendre inclusion-exclusion is correct but **exponentially slower** than direct sieve
- At x=5000, a=19 primes, 2^19 = 524K divisors to enumerate. Takes 0.74s vs 0.0002s for sieve.
- x=10000 would require 2^25 ~ 33M enumerations -- already impractical vs O(x) sieve

### Experiment 3: Sequential vs FFT Polynomial Product Timing

| x | #primes | t_sequential | t_FFT | t_direct_sieve | FFT/seq speedup |
|---|---------|-------------|-------|----------------|-----------------|
| 10,000 | 25 | 0.0006s | 0.0076s | 0.0001s | 0.07x (FFT slower) |
| 50,000 | 48 | 0.0020s | 0.0023s | 0.0004s | 0.89x |
| 100,000 | 65 | 0.0047s | 0.0041s | 0.0008s | 1.13x |

- FFT divide-and-conquer only beats sequential for x >= 100K (marginal)
- Both polynomial methods are **5-50x slower** than direct sieve
- The polynomial product is O(a * x) sequential or O(x log^2 x) FFT -- both worse than O(x log log x) sieve

### Experiment 4: Roots of Unity

P(omega^j) = 0 whenever any sieving prime p divides j*p modulo m. Since 2 is always a sieving prime, P(omega^j) = 0 whenever m | 2j. For small primes {2,3,5}, P is zero at most roots of unity for small m.

Key observation: the non-zero values at roots of unity encode **additive** structure of subset sums, not the **multiplicative** sieve structure needed for pi(x). No useful prime-counting information can be extracted.

### Experiment 5: Dirichlet Polynomial

| x | a = pi(B) | 2^a total divisors | Divisors <= x | Exceeding x | Match |
|---|-----------|-------------------|---------------|-------------|-------|
| 100 | 4 | 16 | 14 | 2 | OK |
| 500 | 8 | 256 | 73 | 183 | OK |
| 1000 | 11 | 2048 | 153 | 1895 | OK |
| 5000 | 19 | 524288 | 690 | 523598 | OK |

- The vast majority of squarefree divisors of the primorial **exceed x** and are wasted work
- For x=5000: only 690 of 524K divisors contribute (0.13%)
- The correct algebraic object is the Dirichlet polynomial F(s) = prod(1 - p^{-s}), not the power series P(z)

### Experiment 6: Complexity Crossover

| x | pi(sqrt(x)) | 2^{pi(sqrt(x))} | x^{2/3} | Winner |
|---|------------|-----------------|---------|--------|
| 10^2 | 4 | 16 | 22 | Legendre |
| 10^3 | 11 | 2048 | 100 | Meissel-Lehmer |
| 10^4 | 25 | 3.4 x 10^7 | 464 | Meissel-Lehmer |
| 10^6 | 168 | 3.7 x 10^50 | 10^4 | Meissel-Lehmer |
| 10^10 | 9592 | inf | 4.6 x 10^6 | Meissel-Lehmer |
| 10^100 | ~10^48 | 2^{10^48} | 10^67 | Meissel-Lehmer |

- Legendre enumeration wins ONLY for x < ~200
- For all practically interesting x, Meissel-Lehmer dominates decisively
- The exponential blowup 2^{pi(sqrt(x))} makes the divisor approach utterly hopeless for large x

### Experiment 7: Polynomial vs Sieve Structure

Verified that P(z) encodes **signed subset-sum counts** (additive), not the **multiplicative sieve indicator**. The nonzero positions of P(z) for primes {2,3,5,7} are {0,2,3,8,9,14,15,17} -- these are subset sums like 2+3=5... wait, 8=3+5, 9=2+7, 14=2+5+7, etc. These have nothing to do with composite numbers or the sieve output.

## Verdict: CLOSED

The polynomial convolution approach to the prime sieve fails due to a **fundamental structural mismatch**:

1. **Polynomial multiplication** (acceleratable by FFT) operates on power series where exponents **add**
2. **The prime sieve** requires multiplicative number theory where divisors **multiply**
3. The correct algebraic framework is **Dirichlet series/convolution**, which has no known sublinear algorithm for partial sums
4. The brute-force Legendre approach enumerates 2^{pi(sqrt(x))} divisors, which is **exponentially worse** than Meissel-Lehmer's O(x^{2/3}/log x) for any x > 10^3

No amount of clever polynomial arithmetic can bridge the gap between additive and multiplicative structure. The approach reduces to the classical Legendre sieve, known since the 1800s, and offers no improvement.

## Relation to Known Work

- The Legendre sieve (1808) is exactly what the divisor enumeration computes
- Meissel (1870) and Lehmer (1959) improved this to O(x^{2/3}) by clever recursive decomposition that avoids full divisor enumeration
- The polynomial product P(z) = prod(1-z^p) appears in partition theory (related to subset-sum problems) but is not useful for prime counting
- The impossibility of fast Dirichlet convolution partial sums is related to open questions in analytic number theory
