# Few Zeros + Primality Search: Results

**Script:** few_zeros_search.py

## What Was Tested
Using R^{-1}(n) plus correction from K zeta zeros to narrow the search interval for p(n), then checking candidates with Miller-Rabin. Also explores FMM, FFT-based zero evaluation, and grouped zero approximation.

## Key Findings
- With K zeros, error ~ sqrt(x) / (K * ln(x)); candidates to check ~ error / ln(x)
- For p(10^100): K=10^6 zeros still leaves ~10^40 candidates -- infeasible
- Error does NOT converge faster than 1/K -- no unexpected cancellation found
- FMM for zero sum: reduces evaluation cost per point but not the NUMBER of zeros needed
- FFT-based batch evaluation: useful for computing pi(x) at many x values simultaneously, but same total work
- Grouped zero approximation: loses precision in the oscillatory phases, introducing O(1) errors per group
- The approach is viable as a practical algorithm but cannot achieve polylog complexity

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence -- narrowing search interval still requires O(sqrt(x)) zeros to get interval width < prime gap)

## One-Line Summary
Few-zeros + primality search: K zeros leave ~sqrt(x)/K candidates; need K~sqrt(x)/ln(x) for unique identification.
