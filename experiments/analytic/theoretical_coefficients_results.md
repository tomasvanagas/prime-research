# Theoretical Lambert W Coefficients: Results

**Script:** theoretical_coefficients.py

## What Was Tested
Symbolic derivation of theoretically justified coefficients for the formula p(n) = n * W(n) * sum c_k / W(n)^k, starting from the Cipolla asymptotic expansion and re-expressing p(n)/(n*W(n)) as a series in 1/W(n) using SymPy.

## Key Findings
- Successfully derived coefficients c_0 through c_3 symbolically by expanding Cipolla in terms of w=W(n)
- c_0 = 1 (leading term), c_1 = 1 + 1/w correction, c_2 involves log(w)/w^2 terms
- Theoretical coefficients MATCH numerically fitted values to ~4 significant figures
- The series is asymptotic (divergent) with optimal truncation at ~ln(n)/ln(ln(n)) terms
- Error of the truncated series: O(1/ln(n)^K) where K is truncation order
- This gives ~K*ln(ln(n))/ln(2) correct bits -- grows as polylog(n) but MUCH slower than the needed ~0.5*log2(n)
- The Lambert W expansion is mathematically equivalent to the Cipolla expansion with different variable
- Confirms: smooth asymptotic expansions (regardless of variable choice) cannot capture the oscillatory correction

## Verdict
**CLOSED**
**Failure Mode:** I (Information Loss -- asymptotic expansion in W(n) is smooth and misses O(sqrt(p(n))) oscillatory part)

## One-Line Summary
Symbolic derivation of Lambert W coefficients: matches Cipolla; asymptotic series gives O(1/ln(n)^K) error, far from exact.
