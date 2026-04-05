# Approximate Degree of pi(x) mod 2: Results

**Script:** `approx_degree_counting.py`
**Session:** 12
**See also:** `approx_degree_prime_results.md`

## What Was Tested
Compared approximate degree of the prime indicator chi_P vs pi(x) mod 2 (counting parity) using LP-based L_inf approximation on {0,1}^N for N=4..10.

## Key Findings
- adeg(pi(x) mod 2) is comparable to adeg(chi_P), both scaling as ~N/2
- Session 12 proved pi(x) mod 2 is as hard as pi(x), so this confirms the counting parity is not easier
- No separation between indicator and counting parity in approximate degree

## Verdict
**CLOSED** -- counting parity is no easier than the indicator.
**Failure Mode:** Equivalence (pi(x) mod 2 reduces to chi_P, same approximate degree)

## One-Line Summary
Approximate degree of pi(x) mod 2 matches chi_P (~N/2); counting parity is as hard as detection.
