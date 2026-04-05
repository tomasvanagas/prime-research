# Additive Combinatorics & Arithmetic Progressions: Results

**Script:** additive_combinatorics.py

## What Was Tested
Five additive combinatorial approaches to computing p(n): Green-Tao AP interpolation, Hardy-Littlewood circle method inversion, Selberg/GPY sieve weight formulas, Goldbach representation extraction, and Erdos-Kac prime detection.

## Key Findings
- Green-Tao APs exist among primes but interpolating from AP structure requires knowing the primes first (circularity)
- Circle method inversion requires evaluating exponential sums over primes, which is equivalent to knowing pi(x)
- Sieve weight formulas (Selberg/GPY) give density estimates but not exact counts without O(x^{2/3}) work
- Goldbach representations encode primes but extracting individual primes requires the full sum
- Erdos-Kac gives distributional information only, not individual prime values

## Verdict
**CLOSED** -- Failure Mode: C (Circularity) / E (Equivalence)

## One-Line Summary
All five additive combinatorial angles require knowing primes to compute primes, or reduce to the explicit formula.
