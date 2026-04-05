# Monotone Complexity of pi(x): Results

**Script:** `monotone_complexity.py`
**Session:** 15

## What Was Tested
Whether monotone circuit lower bound techniques (Razborov approximation method, Alon-Boppana) apply to the threshold function f_k(x) = [pi(x) >= k], which is monotone in the bits of x.

## Key Findings
- f_k(x) IS monotone: flipping a bit from 0 to 1 increases x, hence pi(x) can only increase
- However, monotone complexity lower bounds do NOT transfer to general circuits efficiently: C(f) >= M(f) / 2^N is trivial
- The slice function transfer (tighter for specific Hamming weights) does not apply because f_k is not a slice function
- Razborov's approximation method requires identifying "easy" sub-functions that approximate f_k -- for primality-based thresholds, no such easy sub-functions exist
- Communication complexity of f_k is also exponential (full rank)

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (monotone lower bound techniques do not transfer usefully to general circuits for f_k)

## One-Line Summary
Monotone complexity of [pi(x) >= k] has exponential monotone lower bounds, but these do not transfer to useful general circuit lower bounds.
