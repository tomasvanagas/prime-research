# Prime Race Narrowing: Results

**Date:** 2026-04-04 (Session 7)
**Script:** prime_race_narrowing.py

## What Was Tested
Multi-residue constraint propagation: whether combining wheel factorization, Chebyshev biases, and lattice constraints can narrow candidates for p(n) to O(1) without computing pi(x).

## Key Findings
- Wheel density: after sieving by first k primes, density ~ e^{-gamma}/ln(p_k); for ~1 candidate per gap need P ~ N^{0.56}
- For N=10^102: need primes up to 10^57, requiring 10^54 primes just for the wheel -- far worse than Lucy DP
- Constraint propagation cannot distinguish primes from each other: all primes pass all coprimality tests
- Chebyshev bias: accuracy barely above random baseline (1/phi(q)) for all tested moduli
- Lattice constraint: range from R^{-1} is ~10^51 wide; after wheel mod 2310 still ~10^{47.5} candidates; need ~160 bits to select
- All constraints equivalent to sieving -- no lattice/constraint shortcut exists

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence) -- constraint propagation reduces to sieving; counting step remains O(x^{2/3})

## One-Line Summary
Wheel/constraint narrowing cannot avoid the pi(x) counting step; all constraint propagation is equivalent to sieving.
