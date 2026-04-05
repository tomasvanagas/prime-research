# Prime Race Correction: Results

**Date:** 2026-04-04 (Session 6)
**Script:** prime_race_correction.py

## What Was Tested
Whether multiple prime race biases (Chebyshev bias mod 3, 4, 5, 7, 8, 11, 12) can triangulate the exact error in R^{-1}(n) and help select the correct prime among nearby candidates.

## Key Findings
- Residue class biases measured for first 5000 primes across 7 moduli
- Base approach (nearest prime to Cipolla approximation) vs bias-corrected: negligible improvement
- Chebyshev bias is O(1/ln(p)) per prime -- far too weak to distinguish candidates separated by O(1)
- For p(10^100): bias ~1/235, need to distinguish ~10^48 candidates
- Bias provides ~0.5 bits of information out of ~170 bits needed

## Verdict
**CLOSED**
**Failure Mode:** I (Information Loss) -- Chebyshev bias is O(1/ln(p)), providing negligible information vs the ~170 correction bits needed

## One-Line Summary
Chebyshev prime race bias is O(1/ln(p)) -- provides ~0.5 bits out of ~170 needed; cannot meaningfully correct R^{-1}(n).
