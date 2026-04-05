# NTT Sieve: Results

**Script:** ntt_sieve.py

## What Was Tested
Whether the sieve of Eratosthenes can be accelerated via NTT (Number-Theoretic Transform) by viewing it as a Dirichlet convolution: chi_P = mu_P * 1, computed via Dirichlet series multiplication in the s-domain (Mellin transform), with inversion via Perron's formula.

## Key Findings
- Dirichlet convolution is correctly expressed as pointwise multiplication in the s-domain (Dirichlet series)
- Perron formula inversion: need O(x^{1/2+epsilon}) sampling points on the critical line for error < 0.5
- Each sampling point requires evaluating the Dirichlet series at s = sigma + it, which costs O(x) naively or O(x^{1/2}) via Euler product truncation
- Total cost via Perron: O(x^{1/2+epsilon}) evaluations * O(x^{1/2}) each = O(x^{1+epsilon}) -- worse than direct sieving
- The "NTT speedup" analogy breaks: FFT works because convolution over Z/NZ has O(N) "frequency points"; Dirichlet convolution over Z needs O(x^{1/2+epsilon}) "s-points"

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (E) -- Dirichlet convolution via Perron formula requires O(x^{1/2+epsilon}) critical-line evaluations, reducing to the analytic method (Lagarias-Odlyzko).

## One-Line Summary
NTT sieve via Dirichlet convolution + Perron inversion: needs O(x^{1/2+epsilon}) s-evaluations -- reduces to Lagarias-Odlyzko, no FFT-like speedup.
