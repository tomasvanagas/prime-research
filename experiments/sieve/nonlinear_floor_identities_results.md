# Nonlinear Floor Identities: Results

**Date:** 2026-04-04 (Session 14)
**Script:** nonlinear_floor_identities.py

## What Was Tested
Recursive identities for pi(x) using floor values, integer coefficient search, and higher-order differences of floor functions.

## Key Findings
- Recursive identities: max_err ~2.6 with 3-5 recursive terms; correction comparable to prime gaps
- No simple rational formula exists; best-fit coefficients are irrational
- Higher-order differences: products give lcm-divisibility, not primality
- Overfitting scales with degree: K=20, deg=2 gives 285/290 exact on train, max_err=10.49 on test

## Verdict
**CLOSED** -- see nonlinear_sieve_summary.md for detailed results.
**Failure Mode:** E (Equivalence)

## One-Line Summary
Floor value identities and polynomial fits memorize training data but fail outside; no closed-form integer relation exists.
