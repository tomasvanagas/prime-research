# P-recursive / Holonomic / D-finite Exploration: Results

**Script:** holonomic_prime.py

## What Was Tested
Whether p(n) or transforms of it satisfy linear recurrences with polynomial coefficients (P-recursive / holonomic sequences). Five experiments: direct P-recursiveness test, transforms (Euler, binomial, Mobius, EGF/OGF), matrix power formulation search, Cloitre-type recurrence, and prime zeta / contour integration feasibility.

## Key Findings
- Direct test: p(n) does NOT satisfy any linear recurrence with polynomial coefficients of order up to 6 and degree up to 4
- No transform (Euler, binomial, Mobius) of p(n) is P-recursive either
- If p(n) were holonomic, binary splitting would give O(polylog n) computation; the failure confirms the information barrier
- Cloitre-type recurrences produce primes but are not P-recursive (they use floor/gcd, not polynomial coefficients)
- The prime sequence is provably not D-finite (well-known result: its generating function has natural boundary at |z|=1)

## Verdict
**CLOSED** -- Failure Mode: I (Information Loss)

## One-Line Summary
p(n) is provably not P-recursive/holonomic; binary splitting via linear recurrence is impossible.
