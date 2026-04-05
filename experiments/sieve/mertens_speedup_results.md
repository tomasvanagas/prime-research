# Mertens Function Speedup: Results

**Date:** 2026-04-04 (Session context)
**Script:** mertens_speedup.py

## What Was Tested
Why D(x) = sum d(n) admits O(x^{1/3}) via hyperbola method but M(x) = sum mu(n) is stuck at O(x^{2/3}). Implemented and benchmarked: brute-force M(x), Lucy DP M(x) in O(x^{2/3}), optimized Lucy with memoization, D(x) via hyperbola in O(sqrt(x)), D(x) in O(x^{1/3}) via block decomposition, and analysis of why M(x) resists speedup.

## Key Findings
- D(x) = sum d(n) has a closed-form sum of floor functions (no recursion), enabling O(x^{1/3})
- M(x) is self-referential: M(x) = 1 - sum_{n>=2} M(floor(x/n)), requiring O(sqrt(x)) recursive evaluations
- No known factorization mu = f*g (Dirichlet convolution) exists where both sum_f and sum_g are easy
- H-T achieves O(x^{3/5}) for M(x) via elementary methods but this does NOT help pi(x)
- The Liouville function L(x) is computationally equivalent to M(x)
- M(x) speedup is a genuine theoretical advance but does NOT break the O(x^{2/3}) barrier for pi(x)

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence) -- M(x) speedup to O(x^{3/5}) does not propagate to pi(x); the conversion reintroduces O(x^{2/3})

## One-Line Summary
M(x) admits O(x^{3/5}) via H-T but D(x) admits O(x^{1/3}) due to closed-form structure; neither helps pi(x) below O(x^{2/3}).
