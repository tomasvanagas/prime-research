# Nonlinear Selberg Parity: Results

**Date:** 2026-04-04 (Session 14)
**Script:** nonlinear_selberg_parity.py

## What Was Tested
Whether nonlinear operations can break the Selberg parity barrier and reduce sieve complexity: Legendre sieve, split sieve, prime vs semiprime distinguishing.

## Key Findings
- Parity barrier confirmed: Legendre sieve finds all primes above sqrt(x) but cannot count them efficiently (2^{pi(sqrt(x))} terms)
- Split sieve = conjunction = full Legendre sieve -- no improvement
- Nonlinear ops CAN distinguish primes from semiprimes, but this requires O(sqrt(n)) divisibility tests per number giving O(x^{3/2}) total -- WORSE
- Key: parity barrier is about correctness; the efficiency barrier comes from O(sqrt(x)) distinct floor values

## Verdict
**CLOSED** -- see nonlinear_sieve_summary.md for detailed results.
**Failure Mode:** E (Equivalence)

## One-Line Summary
Nonlinear ops break parity barrier but not efficiency barrier; O(sqrt(x)) floor values remain the bottleneck.
