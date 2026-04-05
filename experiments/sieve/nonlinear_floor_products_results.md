# Nonlinear Floor Products: Results

**Date:** 2026-04-04 (Session 14)
**Script:** nonlinear_floor_products.py

## What Was Tested
Products of floor values floor(x/a)*floor(x/b) as nonlinear features for computing pi(x), including correction terms, polynomial fits, and information content analysis.

## Key Findings
- Correction term delta(a,b,x) does NOT separate primes from composites
- Polynomial fit (K=10, quadratic): max_err=1.34 on training, ~10 on test -- no generalization
- x=10000 has only 199 distinct floor values; products expand to ~6882 but O(x) info still needed

## Verdict
**CLOSED** -- see nonlinear_sieve_summary.md for detailed results.
**Failure Mode:** E (Equivalence) -- nonlinear combinations of floor values still require O(sqrt(x)) distinct values, equivalent to Meissel-Lehmer.

## One-Line Summary
Nonlinear floor products expand feature space but cannot compress O(sqrt(x)) floor values needed for exact pi(x).
