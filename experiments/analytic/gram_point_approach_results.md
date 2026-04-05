# Gram Points + Riemann-Siegel Approach: Results

**Script:** gram_point_approach.py

## What Was Tested
Computing pi(x) using Gram points and the Riemann-Siegel formula: (1) exploiting the sign pattern of Z(t) at Gram points, (2) Odlyzko-Schonhage batching, (3) selective zero contributions (most-significant zeros), (4) optimal test functions in explicit formula.

## Key Findings
- Riemann-Siegel Z(t) can be evaluated in O(sqrt(t)) time per point
- Gram points g_n (where theta(g_n) = n*pi) allow sign detection of Z(t) but Gram's law violations (~27% of zeros) prevent simple zero counting
- Odlyzko-Schonhage algorithm: computes N zeros near height T in O(N * T^eps) amortized -- the state of the art but still O(T^{1/2}) per zero
- Selective zero contribution (using largest cos(gamma*ln(x)) terms): constant-factor improvement, same asymptotics
- Optimal test functions (smooth cutoffs): reduce error constant but not scaling
- The Gram point structure does not bypass the need to account for individual zeros

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence -- Gram points and Riemann-Siegel are efficient tools for the explicit formula but do not change its O(sqrt(x)) requirement)

## One-Line Summary
Gram points + Riemann-Siegel give efficient zero computation but do not reduce the O(sqrt(x)) zero requirement for exact pi(x).
