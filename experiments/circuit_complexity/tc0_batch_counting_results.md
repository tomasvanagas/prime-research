# TC^0 Batch Counting Analysis: Results

**Script:** `tc0_batch_counting.py`
**Session:** 16

## What Was Tested
Whether PRIMES in TC^0 (via BPSW) enables counting primes in [1,x] without enumerating all x candidates. Five escape routes: TC^0 MAJORITY fan-in reduction, divide-and-conquer in TC^0, batch structure in modular exponentiation, algebraic structure of f(k) = 2^{k-1} mod k, and #TC^0 counting.

## Key Findings
- MAJORITY fan-in: can sum x bits in constant depth, but the x input bits must each be COMPUTED (exponential circuit size)
- Divide-and-conquer: splitting [1,x] into halves does not help because each half still requires enumerating candidates
- Batch modular exponentiation: the coupling of exponent and modulus through k makes f(k) = 2^{k-1} mod k cryptographically pseudorandom across consecutive k
- Algebraic structure: f(k) has no polynomial or algebraic relationship across different k values
- #TC^0: whether #TC^0 is in NC is THE key open question; no resolution found

## Verdict
**CLOSED**
**Failure Mode:** Information loss (modular exponentiation is pseudorandom across consecutive k; #TC^0 vs NC is open)

## One-Line Summary
All five TC^0 batch counting routes fail: modular exponentiation is pseudorandom across k; #TC^0 in NC is the key open question.
