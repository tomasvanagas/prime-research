# CRT Reconstruction: Results

**Script:** crt_reconstruction.py

## What Was Tested
Comprehensive empirical study quantifying why CRT-based reconstruction of p(n) fails: number of moduli needed, accuracy of smooth approximation li(x)/phi(q) for pi(x; q, a), information rate per modulus, and candidate set size after partial CRT.

## Key Findings
- Need O(log(p(n))/log(q)) CRT moduli to uniquely determine p(n) -- roughly 100+ moduli for p(10^9)
- Smooth approximation li(x)/phi(q) for pi(x; q, a) has error O(x^{1/2}/phi(q)) -- too large for exact rounding
- Each modulus provides O(log q) bits of information about p(n)
- L-function zeros for non-trivial characters are ADDITIONAL zeros beyond zeta zeros -- cost increases, not decreases
- Partial CRT narrows candidates but still requires O(x^{2/3}) work per modulus to resolve

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (E) -- CRT reconstruction reduces to computing pi(x; q, a) for many (q, a), which collectively requires more work than computing pi(x) directly.

## One-Line Summary
CRT reconstruction quantified: ~100+ moduli needed, each costing O(x^{2/3}); L-function zeros add cost; strictly worse than direct pi(x).
