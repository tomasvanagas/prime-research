# Sieve Function Compression: Results

**Script:** sieve_function_compression.py

## What Was Tested
Whether the Meissel-Lehmer sieve function S(v, p) as a function of v is compressible: piecewise polynomial structure, Fourier sparsity, low-rank approximation via SVD, and whether compression could reduce per-prime update cost from O(sqrt(x)) to O(polylog).

## Key Findings
- S(v, p) is piecewise linear but the number of breakpoints grows with each sieve step -- no bounded-complexity representation
- Fourier analysis: S(v, p) has broadband spectrum; no sparse representation in frequency domain
- SVD of the S(v, p) matrix (indexed by v and sieve step p): singular values decay slowly; 99% energy requires O(sqrt(x)) singular values
- The sieve update S(v,p) = S(v,p-1) - [S(v/p,p-1) - S(p-1,p-1)] introduces new breakpoints at every v = k*p, preventing compression
- Total breakpoints after all sieve steps: O(x/ln(x)) -- same order as the number of primes

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (E) -- sieve function complexity grows with each prime; breakpoints accumulate to O(x/ln(x)); no compression below O(x^{2/3}).

## One-Line Summary
Sieve function S(v,p) compression: breakpoints grow O(x/ln(x)), SVs decay slowly, broadband spectrum -- no sublinear representation.
