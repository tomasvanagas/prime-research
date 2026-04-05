# Lambert W / R^{-1} Error Analysis: Results

**Script:** lambert_w_error_analysis.py

## What Was Tested
Comprehensive analysis of delta(n) = p(n) - R^{-1}(n): exact error at p(10^k), error correlation with prime gaps and local density, error mod small primes, hybrid search range narrowing, binary search with cheaper oracles, and proven vs observed error bounds.

## Key Findings
- At p(10^k), |delta| grows roughly as sqrt(p(n)): delta(10^6) ~ 100, delta(10^9) ~ 3000
- Error correlates weakly (~0.2) with local prime gap; no strong predictive signal
- delta mod small primes: uniform distribution, no exploitable pattern
- Hybrid approach: R^{-1}(n) narrows search to interval of width ~2*sqrt(p(n)); then sieve within that interval
- Binary search with pi(x) oracle: log2(search_interval) ~ 0.5*log2(p(n)) steps, each costing O(p(n)^{2/3})
- RH-conditional bound: |pi(x) - li(x)| < sqrt(x)*ln(x)/(8*pi), matching observed errors
- Unconditional bound (Trudgian): weaker by factor ~100 but same scaling

## Verdict
**CLOSED**
**Failure Mode:** I (Information Loss -- R^{-1}(n) error is O(sqrt(p(n))), containing ~0.5*log2(n) irreducible bits)

## One-Line Summary
R^{-1}(n) error analysis: |delta| ~ O(sqrt(p(n))), no exploitable structure in residuals; confirms ~50% digit barrier.
