# Iterative Zero Refinement: Results

**Script:** iterative_zero_refinement.py

## What Was Tested
An iterative scheme for p(n): start with Li^{-1}(n), compute partial zero sum with T1 << x zeros, refine, iterate with more zeros. Tests whether self-correction (weak dependence of correction on x) reduces total zero count needed.

## Key Findings
- The partial zero sum converges: each batch of zeros improves accuracy
- However, the correction depends on x only through x^{i*gamma} terms -- sensitivity is O(1) in x, not O(x^{-1/2}) as hoped
- Self-correction provides at most a constant-factor improvement (2-3x fewer zeros per iteration), not an asymptotic improvement
- After iteration, still need K = O(x^{1/2}/ln(x)) total zeros to reach error < 0.5
- The convergence rate matches the known explicit formula truncation error: E(K) ~ x^{1/2} * K^{-1/2} * (log K)

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (E) -- iterative refinement converges at the same rate as direct truncation of the explicit formula; no reduction in total zero count.

## One-Line Summary
Iterative zero refinement: self-correction gives constant-factor savings only; total zeros needed remains O(x^{1/2}/ln(x)).
