# Generating Function Coefficient Extraction: Results

**Script:** generating_function_extraction.py

## What Was Tested
Extracting p(n) from generating functions: (1) OGF F(z) = sum p(n)*z^n via Cauchy integral / FFT, (2) Dirichlet generating function G(s) = sum p(n)/n^s, (3) EGF E(z) = sum p(n)*z^n/n!, (4) Prime zeta function P(s) = sum 1/p^s.

## Key Findings
- OGF F(z) has **natural boundary** at |z|=1 (Fabry gap theorem) -- cannot evaluate on contour for Cauchy extraction
- EGF E(z) converges for all z but extracting n-th coefficient requires O(n) evaluations (no fast method)
- Dirichlet series G(s): no Euler product structure, cannot be evaluated without knowing primes -- **circular**
- P(s) = sum mu(k)/k * ln(zeta(ks)): computable from zeta but has natural boundary at Re(s)=0
- Extracting individual primes from P(s) requires inverting the sum, which is as hard as the original problem
- FFT-based coefficient extraction requires N function evaluations for N-th coefficient -- O(n), not polylog

## Verdict
**CLOSED**
**Failure Mode:** C (Circularity -- evaluating prime-encoding GFs requires knowing primes, or extraction costs O(n))

## One-Line Summary
Generating function coefficient extraction fails: OGF has natural boundary; EGF/DGF extraction costs O(n); P(s) inversion is circular.
