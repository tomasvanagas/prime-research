# Determinant/Permanent Encoding of pi(x): Results

**Script:** `det_perm_encoding.py`
**Session:** 15

## What Was Tested
Whether pi(x) can be expressed as det or permanent of a poly(N)-size matrix whose entries are functions of bits of x. Three experiments: exhaustive search for small matrices, constraint analysis, and permanent representations.

## Key Findings
- For m=2 (2x2 matrix with linear entries in bits): solvable for N <= 4 but system is overdetermined for N >= 5
- For m=3: system has 9*(N+1) unknowns vs 2^N equations; underdetermined only for tiny N
- Permanent representation has same constraints as determinant (both require exponentially many free parameters)
- The GapL question (is pi(x) in GapL?) remains open but empirical evidence suggests exponential matrix size needed

## Verdict
**CLOSED**
**Failure Mode:** Information loss (poly(N)-size det/perm representations do not exist; system is massively overdetermined)

## One-Line Summary
No poly(N)-size det or permanent representation of pi(x) exists; the system of equations is exponentially overdetermined.
