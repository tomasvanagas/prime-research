# Optimal Zero Refinement (Sparse Zero Selection): Results

**Script:** optimal_zero_refinement.py

## What Was Tested
Whether selecting optimally-chosen zeta zeros (rather than consecutive) can reduce the number needed for exact pi(x). Also: using GUE/Montgomery pair correlation statistics to approximate the full zero sum without individual zeros.

## Key Findings
- Each zero contributes ~x^{1/2} / |rho| * cos(gamma*ln(x) + phase) to pi(x)
- For a given x, the "most important" zeros are those where cos(gamma*ln(x)) ~ 1 (constructive)
- Greedy selection: choosing highest-contribution zeros gives ~2-5x improvement over consecutive
- BUT: for different x values, different zeros are most important -- no universal sparse subset works
- GUE statistical approximation: replaces individual zero contributions with their expected value, but the variance (which IS the information) is O(1) per zero pair
- Montgomery pair correlation: tells us about spacing statistics but not individual zero positions
- The GUE "replacement" gives the SMOOTH part (R^{-1}) only -- the random residual is the barrier
- Sparse zero selection is x-dependent and cannot be precomputed for arbitrary targets

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence -- sparse selection is x-dependent; GUE statistics give only the smooth part)

## One-Line Summary
Optimal/sparse zero selection: gives 2-5x improvement but is x-dependent; GUE statistics recover only R^{-1}(n).
