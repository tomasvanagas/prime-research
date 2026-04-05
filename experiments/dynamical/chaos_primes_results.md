# Chaos Theory and Dynamical Systems: Results

**Script:** chaos_primes.py

## What Was Tested
Six approaches connecting primes to chaos/dynamical systems: (1) iterative map x_{n+1} = x_n * ln(x_n) / ln(ln(x_n * n)), (2) finding map T(x,n) such that T(p(n),n) = p(n+1) via least-squares fitting, (3) Ulam spiral angular analysis, (4) zeta function inversions, (5) Collatz-like iterations, (6) Lambert W function formulas. Tested on first 22000 primes.

## Key Findings
- Iterative map trajectories diverge from prime sequence after ~5 steps; no stable attractor at primes
- Fitted map T(x,n) achieves low residual on training data but has no predictive power (overfits gap noise)
- Ulam spiral angles show modular structure (primes avoid multiples of 2,3,5) but no new information beyond wheel factorization
- Zeta function inversion recovers the explicit formula -- same summation barrier
- Collatz-like maps either diverge, cycle, or require pi(x) as subroutine
- Lambert W formulas are equivalent to known asymptotic expansions (n*W(n) ~ n*ln(n))

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence -- all approaches reduce to known methods: PNT asymptotics, explicit formula, or sieve)

## One-Line Summary
Chaos/dynamical approaches (iterative maps, Ulam spiral, Collatz-like, Lambert W) all reduce to known asymptotics or require pi(x) as subroutine.
