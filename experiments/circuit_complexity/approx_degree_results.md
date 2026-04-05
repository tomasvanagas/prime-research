# Approximate Degree over GF(p): Results

**Script:** `approx_degree.py`
**Session:** 23
**See also:** `approx_degree_prime_results.md` for the main L_inf approximate degree result

## What Was Tested
Approximate polynomial degree of the prime indicator chi_P over GF(p) for small primes p, using Mobius inversion on the Boolean lattice to compute ANF coefficients mod p. Tests whether chi_P could be in ACC^0 via Razborov-Smolensky.

## Key Findings
- Over GF(p), the multilinear ANF of chi_P has high-degree nonzero coefficients for all primes p tested
- Agreement fraction with low-degree approximations is close to 1/p (trivial random baseline)
- No GF(p) captures useful structure that the real-valued approximate degree misses

## Verdict
**CLOSED** -- subsumed by the main approximate degree result (adeg ~ N/2).
**Failure Mode:** Equivalence (GF(p) analysis confirms random-like behavior seen over R)

## One-Line Summary
GF(p) approximate degree of chi_P shows no ACC^0 structure; agreement with low-degree polynomials matches random functions.
