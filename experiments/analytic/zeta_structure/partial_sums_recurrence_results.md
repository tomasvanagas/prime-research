# Partial Sums Recurrence Analysis: Results

**Script:** partial_sums_recurrence.py

## What Was Tested
Whether partial sums S_K(x) = sum_{k=1}^K cos(gamma_k*log(x))/sqrt(1/4+gamma_k^2) satisfy any recurrence relation of order 1-20, including: (1) linear recurrences, (2) nonlinear (polynomial degree 2,3) recurrences, (3) recurrences on differences d_K = S_K - S_{K-1}, (4) convergence rate analysis.

## Key Findings
- Linear recurrence of order r: residual decreases slowly with r but never reaches machine precision
- Order 20 residual: relative error ~10^{-2} -- NOT a recurrence (true recurrence would give ~10^{-15})
- Nonlinear recurrences: slightly better fit but still relative error ~10^{-3} at degree 3 -- overfitting
- Difference recurrence: d_K = cos(gamma_K*log(x))/sqrt(...) has no simpler recurrence than the defining formula
- Convergence rate: |S_K - S_1000| decays as K^{-0.5} (algebraic, not geometric) -- consistent with random phases
- The partial sums do NOT satisfy any low-order recurrence at ANY of the tested x values
- This is consistent with the zeta zeros being algebraically independent (Schanuel's conjecture)

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence -- no recurrence found; partial sums must be computed term-by-term)

## One-Line Summary
Partial sums of explicit formula satisfy no recurrence of order <=20; convergence is algebraic (K^{-0.5}), not geometric.
