# Random Matrix Theory + Zeta Zero Proxy: Results

**Script:** rmt_zeta_proxy.py

## What Was Tested
Using GUE random matrix eigenvalues as statistical proxies for zeta zeros to compute the correction term sum_rho R(x^rho) without individual zeros. Also: whether bounding the correction tightly enough enables exact p(n) identification from a bounded window.

## Key Findings
- GUE proxy: replacing individual zeros with GUE-distributed eigenvalues gives EXPECTED value of correction = 0 (it cancels in expectation)
- The variance of the GUE correction proxy is O(ln(ln(x))) -- matches actual zeta zero variance (consistent with GUE conjecture)
- But the ACTUAL correction for specific x is a specific realization, not the expectation
- Bounding: GUE statistics give |correction| < C*sqrt(x)*sqrt(ln(ln(x))) with high probability, but this bound is too loose for exact rounding
- Tight bound requires knowing individual zero contributions, defeating the proxy approach
- The correction term contains ~0.5*log2(x) bits of actual information that cannot be replaced by statistical estimates
- RMT gives excellent statistical predictions but cannot determine specific prime positions

## Verdict
**CLOSED**
**Failure Mode:** I (Information Loss -- GUE proxy captures statistics but loses the specific O(sqrt(x)) correction needed for exactness)

## One-Line Summary
GUE random matrix proxy for zeta zeros: captures statistics correctly but loses specific correction value; cannot achieve exactness.
