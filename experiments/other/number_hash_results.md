# Number-Theoretic Hash / Perfect Mapping Approach: Results

**Script:** number_hash.py

## What Was Tested
Five approaches to constructing a perfect mapping n -> p(n): CRT with precomputed tables and residue pattern analysis, minimal perfect hash via polynomial mod M, Pisano-like periodicity of p(n) mod m, Stern-Brocot / Calkin-Wilf tree for primes, and Conway FRACTRAN optimization.

## Key Findings
- CRT residue patterns: p(n) mod m shows no periodicity as a function of n for any tested modulus m
- Polynomial hash p(n) = (a*n^2 + b*n + c) mod M fails catastrophically; primes are not polynomial in n
- Pisano-like periodicity: p(n) mod m is NOT periodic in n (unlike Fibonacci mod m)
- Stern-Brocot / Calkin-Wilf trees enumerate all rationals but provide no shortcut for locating primes
- Conway FRACTRAN: the 14-fraction program generates primes but requires exponentially many iterations

## Verdict
**CLOSED** -- Failure Mode: C (Circularity) / I (Information Loss)

## One-Line Summary
No perfect hash or periodic pattern maps n to p(n); prime residues mod m have no periodic dependence on n.
