# Arithmetic Term Construction for p(n): Results

**Script:** arithmetic_term.py

## What Was Tested
Whether a practical fixed-length arithmetic term (using +, -, *, //, ^, floor, mod) can compute p(n) with bounded intermediate values. Four approaches: optimized Willans formula, Miller-Rabin as arithmetic term, R^{-1}(n) with rounding, and sieve encoding via matrix operations.

## Key Findings
- Willans formula via Wilson's theorem requires O(n*ln(n)) multiplications of numbers up to k ~ 10^102; totally infeasible
- Miller-Rabin can be expressed as an arithmetic term but encoding "find the nth number passing MR" still requires a loop/sum over all candidates
- R^{-1}(n) rounding fails because error ~ sqrt(p(n)) >> gap ~ ln(p(n)) by factor ~10^47
- Matrix encoding of sieve requires matrix size proportional to the sieve range, not polylog

## Verdict
**CLOSED** -- Failure Mode: E (Equivalence) / I (Information Loss)

## One-Line Summary
Arithmetic terms for p(n) exist theoretically (Prunescu-Shunia) but all practical constructions have superexponential intermediate values.
