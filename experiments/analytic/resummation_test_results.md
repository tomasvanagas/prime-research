# Resummation Techniques vs Summation Barrier: Results

**Script:** resummation_test.py

## What Was Tested
Five resummation methods applied to the oscillatory zero sum in the explicit formula, tested at pi(7919)=1000, pi(104729)=10000, pi(1299709)=100000: (1) Borel summation / Borel-Pade, (2) Pade approximants of partial sums, (3) Euler-Maclaurin with zero density integration, (4) sequence acceleration (Richardson, Aitken, Shanks), (5) Mellin-Barnes / saddle-point.

## Key Findings
- Borel summation: the zero sum is NOT Borel summable -- the Borel transform has singularities on the positive real axis
- Borel-Pade: approximates the Borel transform but inherits the singularity problem
- Pade approximants of partial sums: converge to incorrect value due to essential singularity
- Euler-Maclaurin: replaces discrete sum with integral + corrections; the integral IS the smooth part (R^{-1})
- Shanks transformation: best practical method, ~5x acceleration, but hits a wall at O(sqrt(x)) zeros
- Mellin-Barnes: reduces to the Perron integral, no savings
- All methods confirm: the zero sum is fundamentally non-summable to fewer than O(sqrt(x)) terms

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence -- the zero sum is not amenable to resummation below O(sqrt(x)) terms)

## One-Line Summary
Five resummation methods (Borel, Pade, Euler-Maclaurin, Shanks, Mellin-Barnes): zero sum is not resummable below O(sqrt(x)) terms.
