# Error-Correcting Code Perspective on Primes: Results

**Script:** ecc_approach.py

## What Was Tested

Whether treating the approximate formula (Cipolla/R^{-1}) as a "noisy channel" and the exact prime as the "codeword" enables error correction. Four approaches: (1) modular syndrome decoding via small primes, (2) nearest prime decoding from approximation, (3) prime constellation matching, (4) wheel factorization + position recovery with ensemble voting.

## Key Findings

- Wheel position prediction: primes mod 30 have 8 valid residues (1,7,11,13,17,19,23,29) but residue is unpredictable from n alone
- Nearest prime decoding: works for small n but hit rate decreases as primes grow (gaps grow relative to approximation error)
- Prime gap distribution: follows known Cramer-type models but is not predictive for individual gaps
- Ensemble voting: combines multiple heuristics but accuracy degrades for larger n
- The "syndrome" (error pattern) encodes the same irreducible information as the prime itself
- No error-correcting structure exists because the "noise" delta(n) = p(n) - approx(n) is information-theoretically incompressible

## Verdict

**CLOSED** -- Failure Mode: Information Loss (I). The error term delta(n) carries ~50% of the bits of p(n) and has no exploitable structure for error correction.

## One-Line Summary

ECC perspective fails -- the approximation error delta(n) = p(n) - R^{-1}(n) is incompressible noise, not a correctable error pattern.
