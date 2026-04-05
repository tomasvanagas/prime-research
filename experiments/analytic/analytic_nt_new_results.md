# Seven Lesser-Known Analytic NT Approaches: Results

**Script:** analytic_nt_new.py

## What Was Tested
Eight analytic number theory approaches not previously explored: (1) Nyman-Beurling-Baez-Duarte criterion, (2) Li's criterion / Keiper-Li sequence, (3) Weil's explicit formula with optimized test functions, (4) Guinand's formula (Bessel variant), (5) Deuring-Heilbronn zero repulsion, (6) Turan's power sum method, (7) Brun's exact sieve with error terms, (8) Dusart/Rosser-Schoenfeld tight bounds.

## Key Findings
- Nyman-Beurling: characterizes RH equivalence but does not compute pi(x)
- Li's criterion / Keiper-Li: sequence growth rate encodes zero locations but extracting pi(x) requires all zeros
- Weil explicit formula with optimized test functions: reduces constant factor but not asymptotic O(sqrt(x)) zero requirement
- Guinand's formula: Bessel function variant of explicit formula; same zero sum barrier
- Deuring-Heilbronn: zero repulsion gives local improvements but cannot eliminate the global zero sum
- Turan's power sum: requires O(x) evaluations to detect primes
- Brun's sieve: exact with error terms, but error terms require zeta zeros to bound
- Dusart bounds: tightest proven bounds still have O(sqrt(x)/log(x)) gap, insufficient for exactness

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence -- all reduce to explicit formula / zeta zero information)

## One-Line Summary
Eight lesser-known analytic NT methods (Nyman-Beurling, Li, Guinand, Turan, etc.) all reduce to the explicit formula barrier.
