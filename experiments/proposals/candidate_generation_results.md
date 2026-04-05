# Candidate Generation + AKS/Miller-Rabin Verification -- Results

**Script:** candidate_generation.py

## What Was Tested
Generate a small candidate set guaranteed to contain p(n) using R^{-1}(n) plus progressive filters (mod 30 coprime residues, zeta zero oscillatory correction, small-prime sieving), then verify each candidate for primality. Measure how small the candidate set can get.

## Key Findings
- R^{-1}(n) gives initial interval of width O(sqrt(x) * log(x)) under RH.
- Mod 30 filter (coprime residues to 2*3*5): reduces candidates by factor ~30/8 = 3.75x.
- Small-prime sieving further reduces by a constant factor (product of (1-1/p) for primes up to sieve bound).
- Zeta zero correction (with K zeros) narrows interval width by factor ~K, but computing each zero contribution costs O(polylog).
- After all cheap filters, candidate set size is still O(sqrt(x) / log^2(x)) -- far from O(1).
- Primality testing each candidate is fast (O(polylog) per test), but the bottleneck is WHICH candidate is the nth prime, not whether candidates are prime.
- Determining the correct index among surviving candidates requires counting, which is O(sqrt(x)).

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (E) -- candidate set remains O(sqrt(x)/log^2(x)) after all cheap filters; index determination costs O(sqrt(x)).

## One-Line Summary
Candidate generation + verification: filters reduce candidates by constant factors; O(sqrt(x)/log^2(x)) candidates remain, index still hard.
