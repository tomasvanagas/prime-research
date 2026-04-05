# Dynamical Systems for Primes (Session 8): Results

**Script:** dynamical_primes.py

## What Was Tested
Five families of dynamical/ergodic approaches: (1) Collatz-like prime-generating maps and Rowland acceleration, (2) Furstenberg topology computational content extraction, (3) ergodic averages along primes -- reversal attempts, (4) symbolic dynamics on prime indicator sequences, (5) fixed points of iterated maps at primes (Newton on Chebyshev psi). Tested against reference primes up to 100000.

## Key Findings
- Rowland sequence a(n) = a(n-1) + gcd(n, a(n-1)) produces primes but at O(p^2) cost per prime -- worse than trial division
- Furstenberg topology quantification recovers Mertens' theorem; no computational shortcut
- Ergodic averages converge to predicted densities but cannot be inverted to extract individual p(n)
- Symbolic dynamics: prime indicator has near-maximal block complexity; no finite automaton generates it
- Newton iteration on Chebyshev psi converges to R_inv(n), not p(n) -- same smooth approximation

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence -- Rowland is O(p^2), Furstenberg recovers sieve theory, ergodic averages give density not individuals, symbolic dynamics has maximal complexity)

## One-Line Summary
Collatz/Rowland maps are O(p^2) per prime; Furstenberg topology, ergodic averages, and symbolic dynamics provide no computational shortcut.
