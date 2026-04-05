# Novel Exact Formula Attempts: Results

**Script:** novel_exact_formula.py

## What Was Tested
Five formula ideas: exact floor rounding of R^{-1}(n), binary digit extraction via pi(x) at specific points, implicit definition via fixed-point contraction, number-theoretic transform (NTT/DFT) approach, and compositional formula p(n) = f1(f2(...fk(n)...)).

## Key Findings
- Exact rounding: error |R^{-1}(n) - p(n)| ~ O(sqrt(p(n))) vastly exceeds gap ~ O(log^2(p(n))); error/gap ratio ~10^47 for p(10^100)
- Even under RH, explicit formula error O(sqrt(x)*log(x)) >> gap by factor ~10^47
- Binary digit extraction: each bit depends on pi(x) at a specific point, which is O(x^{2/3}) to compute
- Fixed-point iteration T(x) = n*ln(x) converges to ~n*ln(n) (PNT approximation) but NOT to p(n) exactly
- NTT/DFT: the prime indicator has high-frequency components (from zeta zeros) that cannot be cheaply isolated
- Compositional formulas reduce to nested approximations that don't improve beyond R^{-1}(n)

## Verdict
**CLOSED** -- Failure Mode: I (Information Loss)

## One-Line Summary
All exact formula attempts fail because the rounding/correction step requires information equivalent to O(sqrt(x)) zeta zeros.
