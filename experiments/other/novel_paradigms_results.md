# Genuinely Novel Paradigms for Computing nth Prime: Results

**Script:** novel_paradigms.py

## What Was Tested
Eight unconventional directions: Levin Universal Search / algorithmic information theory, cellular automata (Rules 30, 110, 150), Busy Beaver connections, Curry-Howard / type-theoretic structure, topos-theoretic analysis, hypercomputation / real-number oracles, CRT + function-of-p(n) encoding, and reverse mathematics / axiom strength.

## Key Findings
- Levin Universal Search: theoretical optimality but constant factor is 2^{K(program)}, astronomically large
- Cellular automata: Rules 30/110 are Turing-complete but no CA rule generates primes faster than O(n) steps
- Busy Beaver: BB(n) grows faster than any computable function; no useful connection to primes
- Curry-Howard: the proof that p(n) is prime encodes a Pratt certificate, not a shortcut to finding p(n)
- Topos theory: internal logic of Zariski topos relates to primes but provides no computational algorithm
- Hypercomputation: real-number oracles could solve it but don't exist in standard computation models
- CRT encoding: requires knowing p(n) mod small primes, which is circular
- Reverse mathematics: computing p(n) requires at least ACA_0, confirming it is not computationally trivial

## Verdict
**CLOSED** -- Failure Mode: C (Circularity) / E (Equivalence)

## One-Line Summary
All eight unconventional paradigms (cellular automata, Busy Beaver, topos, hypercomputation, etc.) fail to provide a polylog algorithm.
