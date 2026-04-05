# Nonlinear Multiplicative Structure: Results

**Date:** 2026-04-04 (Session 14)
**Script:** nonlinear_multiplicative_structure.py

## What Was Tested
GCD structure, floor value collisions, quadratic sieve analogs, and nonlinear Mobius combinations for computing pi(x).

## Key Findings
- GCD(floor(x/k), k) does not distinguish primes from composites
- Only ~30% of collision groups contain primes -- not a useful discriminator
- Nonlinear Mobius (quadratic in M(x/k) values): max_err=14.6 train, 44.4 test -- complete failure

## Verdict
**CLOSED** -- see nonlinear_sieve_summary.md for detailed results.
**Failure Mode:** E (Equivalence)

## One-Line Summary
Multiplicative structure of floor values (GCD, collisions, quadratic Mobius) fails to reduce complexity below O(x^{2/3}).
