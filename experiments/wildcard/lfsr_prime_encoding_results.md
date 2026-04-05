# LFSR Prime Encoding: Results

**Script:** lfsr_prime_encoding.py

## What Was Tested
Three approaches: (1) linear complexity of the prime indicator sequence over GF(2) via Berlekamp-Massey, (2) Mills' constant computability, (3) Matiyasevich-Robinson polynomial for prime generation.

## Key Findings
- Berlekamp-Massey on prime indicator over GF(2): linear complexity = N/2 for all tested sequence lengths -- maximal, confirming the sequence is pseudorandom over GF(2)
- This matches Session 24 results: L/N = 0.5 for all finite fields tested
- Mills' constant: circular -- computing A requires knowing all primes; no independent characterization exists
- Matiyasevich-Robinson polynomial: the 26-variable polynomial exists but finding the "right" inputs for the nth prime requires search over exponentially many candidates
- No shortcut to evaluate the polynomial at inputs producing p(n) without already knowing p(n)

## Verdict
**CLOSED**
**Failure Mode:** Circularity (C) for Mills and Matiyasevich; Information Loss (I) for LFSR -- linear complexity N/2 means the sequence is incompressible over GF(2).

## One-Line Summary
LFSR prime encoding: linear complexity = N/2 (maximal); Mills' constant circular; Matiyasevich polynomial requires exponential search.
