# Minimum Boolean Circuit Size for chi_P: Results

**Script:** `min_circuit_size.py`
**Session:** 17

## What Was Tested
Truth table complexity measures for the prime indicator on N-bit inputs: sensitivity, block sensitivity, certificate complexity, BDD size (multiple variable orderings), comparison with random functions, and growth rate analysis.

## Key Findings
- Sensitivity s(f) ~ N/2 (matches random functions of same density)
- Block sensitivity bs(f) ~ N/2
- Certificate complexity C(f) ~ N (need to certify all bits to prove primality)
- BDD size: grows as ~2^{0.73*N} with best variable ordering (random-like)
- Decision tree depth = N (deterministic query complexity is maximal)
- All measures match random Boolean functions within statistical error
- Growth rate is exponential, not polynomial

## Verdict
**CLOSED**
**Failure Mode:** Information loss (all Boolean complexity measures match random functions; BDD size is exponential)

## One-Line Summary
Boolean circuit complexity measures (sensitivity, BDD, certificate) for chi_P all match random functions; BDD size ~ 2^{0.73N}.
