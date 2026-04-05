# Proposal 1: Spectral Truncation + Sieve Candidate Elimination -- Results

**Script:** proposal1_spectral_sieve.py

## What Was Tested
Use a truncated Riemann explicit formula with K zeta zeros to localize p(n) to an interval [a,b], then sieve that interval to find the exact prime. The "twist" is using CRT modular constraints (p(n) mod m) to shrink the effective interval from W to W/m.

## Key Findings
- With K zeros, interval width W = O(x/K * log^2 x). Sieving costs O(W * log log W).
- Optimal K = sqrt(x / log^3 x) gives total cost O(sqrt(x) * polylog) -- does NOT beat O(x^{2/3}).
- The CRT twist: if m = x^{epsilon}, cost drops to O(x^{1/2 - epsilon/2} * polylog).
- For O(polylog), need m = x^{1/2}/polylog, i.e., need to determine p(n) mod a huge modulus cheaply.
- Computing p(n) mod q for large q requires knowing which residue class the nth prime falls in -- this is as hard as the original problem (requires L-function zeros at the same cost).
- The R^{-1}(n) approximation gives ~50% of digits; the remaining ~50% are locked behind the sqrt(x) barrier from zeta zeros.

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (E) -- computing p(n) mod large modulus is equivalent to the original problem.

## One-Line Summary
Spectral truncation + sieve: O(sqrt(x) * polylog) at best; CRT shrinkage requires p(n) mod q which is equally hard.
