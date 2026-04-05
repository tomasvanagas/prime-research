# Additive Construction Approach: Results

**Script:** additive_construction.py

## What Was Tested
Whether p(n) can be decomposed additively as p(n) = A(n) + B(n) where A(n) is easily computable and B(n) is small enough to determine by brute search. Tested modular computation of the correction delta(n) = p(n) - R^{-1}(n).

## Key Findings
- delta(n) has 5.04 bits/prime entropy and autocorrelation r(1) = 0.996 (smooth random walk)
- B(n) = p(n) - A(n) is as hard to compute as p(n) itself for any smooth A(n)
- Modular decomposition (delta mod 2, mod 3, etc.) does not yield independently computable residues
- The correction term carries ~50% of the information content, which is incompressible

## Verdict
**CLOSED** -- Failure Mode: I (Information Loss)

## One-Line Summary
No additive decomposition separates p(n) into a cheap computable part plus a small brute-forceable correction.
