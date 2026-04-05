# Hypergeometric and Modular Special Functions: Results

**Script:** hypergeometric_modular.py

## What Was Tested

Systematic exploration of hypergeometric functions, modular forms, q-series, theta functions, and elliptic functions for an O(polylog n) formula for p(n). Eight major sub-approaches (~20 sub-approaches total): hypergeometric 2F1 grid search, modular form coefficient encoding, CM form prime/composite distinction, theta function r_k(n), partition congruences, q-Pochhammer, Hecke L-function explicit formula, and j-invariant / Heegner numbers.

## Key Findings

- Hypergeometric 2F1 grid search: best match 6/20 primes (cubic polynomial fit), diverges by n=30
- Modular form coefficients: multiplicativity + growth constraints prevent encoding p(n)
- CM forms: distinguishing prime/composite eigenvalues requires factoring (circular)
- Theta function r_k(n) = divisor sums: need O(sqrt(n)) to compute
- Partition congruences: primes appear as moduli, not as outputs
- q-Pochhammer at fixed q: smooth in n, no prime discontinuity
- Hecke L-function: same zero-counting cost as Riemann explicit formula
- j-invariant: only 9 Heegner numbers, finite resource

## Verdict

**CLOSED** -- Failure Mode: Circularity (C) / Equivalence (E) / Information Loss (I). All special function approaches either require primes as input, reduce to the explicit formula, or cannot encode the irregular prime sequence.

## One-Line Summary

Systematic search across hypergeometric, modular, q-series, theta, and elliptic functions finds no polylog formula -- all ~20 sub-approaches hit circularity, equivalence, or information barriers.
