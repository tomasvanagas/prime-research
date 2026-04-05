# Nonlinear Comparisons: Results

**Date:** 2026-04-04 (Session 14)
**Script:** nonlinear_comparisons.py

## What Was Tested
Comparisons and thresholding of floor values: parity of floor(x/k), gap functions, modular floor values, and bitwise operations (XOR/AND/OR) as potential prime-counting features.

## Key Findings
- Parity of floor(x/k): correlation with primality ~0 (0.009 at x=100)
- Gap function g(k,x): 96%+ of primes are invisible (gap=0 at prime k)
- Modular floor values: correlation < 0.12
- Bitwise operations carry no useful prime signal

## Verdict
**CLOSED** -- see nonlinear_sieve_summary.md for detailed results.
**Failure Mode:** E (Equivalence)

## One-Line Summary
Comparisons/thresholding of floor values carry near-zero correlation with primality; 96%+ of primes invisible to gap function.
