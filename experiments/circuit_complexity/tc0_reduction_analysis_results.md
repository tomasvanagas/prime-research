# TC^0 Reduction Analysis for pi(x): Results

**Script:** `tc0_reduction_analysis.py`
**Session:** 13

## What Was Tested
Whether pi(x) can be reduced to functions known to be in TC^0 (addition, multiplication, division, powering, sorting, MAJORITY). Combined theoretical analysis of what is/is not in TC^0 with computational experiments.

## Key Findings
- TC^0 contains: addition, multiplication, division, iterated multiplication, powering x^n mod m (HAB 2002)
- TC^0 does NOT contain (or unknown): factoring, PRIMES, discrete log
- pi(x) = sum of primality indicators; if PRIMES in TC^0, then the sum is TC^0 BUT requires 2^N parallel evaluations
- The reduction pi(x) -> "count satisfying inputs of TC^0 circuit" is #TC^0, not TC^0 itself
- No known TC^0 reduction avoids enumerating all candidates

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (pi(x) reduces to #TC^0 counting, which is not known to be in NC)

## One-Line Summary
pi(x) reduces to #TC^0 (counting TC^0 circuit satisfying inputs); no TC^0 reduction avoids candidate enumeration.
