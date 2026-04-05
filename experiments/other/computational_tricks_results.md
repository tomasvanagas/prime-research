# Computational Tricks for p(n): Results

**Script:** computational_tricks.py

## What Was Tested
Seven approaches from deep number theory: Kim's non-abelian Chabauty / anabelian geometry, Deninger's dynamical system for zeta zeros, Haran's F_1 geometry, Zagier-style formulas, integer complexity and addition chains, Conway-Guy / Stohr sequences, and inverse sieve via CRT construction.

## Key Findings
- Kim's anabelian geometry gives Diophantine finiteness results but no explicit prime formula
- Deninger's dynamical system is conjectural; even if it exists, extracting zeros requires spectral theory equivalent to the explicit formula
- F_1 geometry is a framework, not a computational method; no known algorithm emerges
- Zagier-style formulas are asymptotic, not exact
- Integer complexity of p(n) grows as O(log p(n)), not providing a shortcut
- Conway-Guy sequences are not related to primes in a useful way
- CRT inverse sieve requires knowing p(n) mod small primes, which is as hard as knowing p(n)

## Verdict
**CLOSED** -- Failure Mode: E (Equivalence) / C (Circularity)

## One-Line Summary
Deep number-theoretic frameworks (anabelian, F_1, Deninger) remain conjectural and provide no computational shortcut.
