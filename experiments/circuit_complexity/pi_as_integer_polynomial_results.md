# pi(x) as Integer Polynomial: Results

**Script:** `pi_as_integer_polynomial.py`
**Session:** 13

## What Was Tested
Degree and coefficient structure of the unique multilinear polynomial over Z representing the prime indicator chi_P: {0,1}^N -> {0,1}. Compared with the GF(2) ANF to see if integer arithmetic provides lower degree.

## Key Findings
- The multilinear representation over Z has degree N (same as over GF(2))
- Integer coefficients grow exponentially: max |c_S| ~ 2^{N/2}
- The number of nonzero coefficients is ~50% of all 2^N monomials (random-like sparsity)
- Working over Z does NOT reduce degree compared to GF(2); the representations are equivalent in structure
- High-degree terms cannot be eliminated by integer cancellations

## Verdict
**CLOSED**
**Failure Mode:** Information loss (integer polynomial has same degree N and random-like structure as GF(2) ANF)

## One-Line Summary
Multilinear polynomial for chi_P over Z has degree N with exponential coefficients; no improvement over GF(2).
