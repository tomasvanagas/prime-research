# V9 Local Formula: Results

**Date:** 2026-04-04 (Session context)
**Script:** v9_local_formula.py

## What Was Tested
Local prime counting approach: use R^{-1}(n) as initial estimate, sieve a small window around it, and use R(x) as a proxy for pi(x) to determine the offset. Falls back to Lucy DP when local formula fails.

## Key Findings
- Local formula accuracy degrades with n: 77% exact at n<=100, 32.4% at n<=1000, 12.4% at n<=10000
- First failure at n=49: predicted 229 instead of 227
- R(x) as pi(x) proxy: error grows from +0.77 at x=100 to +29.42 at x=10^6
- R(x) - pi(x) oscillates with amplitude O(sqrt(x)*ln(x)/(8*pi)), exceeding 0.5 around x~100-200
- Local formula works only when |R(x) - pi(x)| < 0.5 -- true for tiny x only
- V9 with Lucy DP fallback: correct but costs O(x^{2/3}) for the fallback

## Verdict
**CLOSED**
**Failure Mode:** I (Information Loss) -- R(x) approximation error exceeds 0.5 for moderate x; exact pi(x) requires O(x^{2/3})

## One-Line Summary
Local formula using R(x) as pi(x) proxy works only for tiny n (first failure at n=49); fundamentally limited by R(x) oscillatory error.
