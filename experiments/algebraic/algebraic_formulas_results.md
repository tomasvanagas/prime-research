# Algebraic and Closed-Form Formulas: Results

**Script:** algebraic_formulas.py

## What Was Tested

Five algebraic/closed-form approaches to computing p(n): (1) Mills' constant, (2) Willans' formula (Wilson's theorem), (3) prime-counting via floor/trig sums, (4) interpolation / curve-fitting of the prime sequence (main focus), (5) modular exponentiation / AKS-adjacent ideas.

## Key Findings

- Mills' constant: computing A requires knowing the primes; digits of A needed grow with n
- Willans' formula: requires O(x^2) factorial computations -- no advantage over sieving
- Floor/trig sums: theoretically exact but computationally equivalent to sieving
- Interpolation (main findings):
  - Asymptotic correction model: ~0.01-0.05% relative error
  - Polynomial interpolation: exact on training set, catastrophic extrapolation
  - Rational functions: comparable to asymptotic model, no advantage
  - Cipolla expansion + residual fitting: best overall but CANNOT predict exact primes
  - Error is always a significant fraction of the prime gap
  - KEY INSIGHT: no smooth function can predict exact primes; best smooth approximation has error O(sqrt(p)/ln(p)) which exceeds gap O(ln(p)) for large p
- AKS/modular: no algebraic shortcut; binary search + Meissel-Lehmer remains best

## Verdict

**CLOSED** -- Failure Mode: Information Loss (I) / Equivalence (E). Exact formulas (Mills, Willans) are circular or exponential; interpolation hits the fundamental irregularity of primes.

## One-Line Summary

Mills' constant, Willans' formula, floor/trig sums, and interpolation all fail -- exact formulas are circular and smooth approximations cannot bridge the prime gap.
