# Session 10 Dynamical Systems: Results

**Script:** s10_dynamical_systems.py

## What Was Tested
Six dynamical systems strategies tested for O(polylog n) prime computation: (1) Conway's PRIMEGAME (FRACTRAN), (2) Collatz-like maps T(x) = x + pi(x), (3) symbolic dynamics -- prime indicator shift space, (4) transfer operator / Ruelle zeta (Mayer), (5) Furstenberg topology quantitative analysis, (6) cellular automata for prime generation.

## Key Findings
- FRACTRAN: step count grows ~O(p^2.5) with register values of size O(2^p); inherently sequential, no skip possible
- Collatz-like T(x) = x + pi(x): requires pi(x) as subroutine (O(x^{2/3})); gap p(n+1) - T(p(n)) does not converge
- Symbolic dynamics: prime indicator has topological entropy ~1, all binary blocks appear; not sofic or SFT; no finite automaton generates it
- Transfer operator (Mayer): can find individual zeta zeros in O(polylog) but summing O(sqrt(x)) zeros still required for pi(x)
- Furstenberg topology: quantification recovers classical sieve = O(x log log x)
- Cellular automata: best rule matches chi_P at ~random baseline; sieve-as-CA requires O(sqrt(N)) steps and O(N) space
- Compression ratio of chi_P[2:1002]: ~0.15 (slightly better than random, due to density ~1/ln(n))

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence -- all six reduce to known barriers: sequential iteration, summation barrier, or sieve)

## One-Line Summary
All 6 dynamical systems approaches (FRACTRAN, Collatz, symbolic dynamics, transfer operator, Furstenberg, CA) reduce to known O(x^{2/3}) or worse barriers.
