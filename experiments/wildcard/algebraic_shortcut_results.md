# Algebraic Shortcut: Results

**Script:** algebraic_shortcut.py

## What Was Tested
Six algebraic/arithmetic approaches to bypass analytic number theory: Wilson's theorem for prime detection, Bernoulli number/zeta value extraction, cyclotomic polynomial structure, polynomial evaluation mod n, matrix permanent/determinant formulation, and the factorial prime counting formula.

## Key Findings
- Wilson's theorem approach: computing (k-1)! mod k for each k is O(k) per evaluation; summing over k gives O(x^2) total -- worse than sieving
- Costa-Gerbicz-Harvey (2014) computes n! mod p in O(sqrt(p) polylog p), but this still requires per-prime work
- Cyclotomic polynomials: degree phi(n) requires knowing factorization; no shortcut to prime counting
- Matrix formulation: expressing pi(x) as a determinant would require a matrix of size O(x/ln(x))
- Factorial prime counting formula is exact but inherently O(n) with expensive per-term factorials

## Verdict
**CLOSED**
**Failure Mode:** Circularity (C) -- all algebraic reformulations either require knowing primes or reduce to computations at least as expensive as direct counting.

## One-Line Summary
Algebraic shortcuts (Wilson, cyclotomic, determinant, factorial formula) all require O(x) or worse per-element computation -- no polylog path.
