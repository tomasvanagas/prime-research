# Dynamical Systems Approaches (Session 4): Results

**Script:** dynamical_systems.py

## What Was Tested
Seven dynamical systems approaches: (1) Conway's PRIMEGAME (FRACTRAN) acceleration, (2) Bost-Connes system KMS states, (3) symbolic dynamics complexity of prime indicator, (4) Selberg trace / hyperbolic geodesics, (5) Collatz-like prime-detecting maps, (6) Ulam spiral / lattice dynamics, (7) ergodic averages on multiplicative functions.

## Key Findings
- FRACTRAN/PRIMEGAME: register grows as 2^p, requiring O(2^p) steps to reach prime p; no acceleration possible
- Bost-Connes system: partition function Z(beta) encodes primes but evaluation at beta=1 requires summing over all primes
- Symbolic dynamics: prime indicator shift space has topological entropy ~1 (near full shift); not sofic or SFT
- Selberg trace formula: eigenvalue-to-prime correspondence is exact but computing eigenvalues requires O(sqrt(x)) terms
- Collatz-like maps: orbits are unpredictable without running them; no shortcut to step n
- Ergodic averages: converge to Mertens-type constants; cannot extract individual primes

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence -- all reduce to known barriers: FRACTRAN is exponential, Bost-Connes requires infinite sums, Selberg trace = explicit formula)

## One-Line Summary
Seven dynamical systems frameworks (FRACTRAN, Bost-Connes, symbolic dynamics, Selberg trace, Collatz, lattice, ergodic) all reduce to known barriers.
