# Prime Zeta Oracle / Prony's Method: Results

**Script:** prime_zeta_oracle.py

## What Was Tested
Extracting individual primes from the prime zeta function P(s) = sum p^{-s} using: (1) evaluating P(s) at N values and solving for p(1)..p(N) (Vandermonde-like system), (2) exponential signal E(t) = sum e^{-t*p(k)} with Prony's method, (3) computing E(t) from -zeta'(s)/zeta(s) without knowing primes.

## Key Findings
- P_N(s) at N values of s: the system is a Vandermonde-like system in p(n)^{-s_k}; solvable in principle but condition number grows EXPONENTIALLY with N -- numerically impossible for N > ~20
- Prony's method: recovers frequencies (primes) from signal samples, but needs O(N) samples to recover N primes
- Computing E(t) from zeta: E(t) = sum_p e^{-tp} relates to integral of -zeta'(s)/zeta(s) via Mellin; evaluating this integral requires the explicit formula
- The approach is fundamentally **circular**: computing E(t) or P(s) efficiently already requires the Euler product (i.e., knowing primes)
- Zeta can be computed from its Dirichlet series (not needing primes), but extracting individual primes back from zeta IS the hard problem
- Prony/ESPRIT: works for recovering O(K) exponentials from O(K) samples but cannot skip to the n-th one in polylog time

## Verdict
**CLOSED**
**Failure Mode:** C (Circularity -- evaluating E(t) or extracting primes from P(s) requires the Euler product or explicit formula)

## One-Line Summary
Prime zeta oracle / Prony: extracting primes from P(s) has exponential condition number; E(t) computation reduces to explicit formula.
