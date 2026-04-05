# Additive Number Theory & Combinatorics Approaches: Results

**Script:** additive_theory.py

## What Was Tested
Seven approaches from additive number theory: partition-based formulas (Ono-Craig inversion), Goldbach additive basis inversion, Mobius function cumulative patterns, Mertens function shortcuts, Legendre sieve acceleration, Dirichlet series / prime zeta inversion, and smooth number counts. Tested for n=1..10000.

## Key Findings
- Partition-based formulas (Ono-Craig) provide a primality criterion but computing the partition functions requires O(sqrt(n)) divisor sums per candidate
- Goldbach inversion is circular: needs primes to count representations
- Mobius/Mertens cumulative sums have the same complexity as pi(x) via inclusion-exclusion
- Legendre sieve acceleration still requires O(x^{2/3}) work at best
- Dirichlet series inversion reduces to the explicit formula with zeta zeros
- Smooth number counts via Dickman's function give approximations but not exact values

## Verdict
**CLOSED** -- Failure Mode: C (Circularity) / E (Equivalence)

## One-Line Summary
All seven additive number theory angles either require primes as input or reduce to the explicit formula.
