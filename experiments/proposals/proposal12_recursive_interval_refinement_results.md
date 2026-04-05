# Proposal 12: Recursive Interval Refinement via Pratt Certificates -- Results

**Script:** proposal12_recursive_interval_refinement.py

## What Was Tested
Hierarchical approach: (1) R^{-1}(n) gives x_approx in O(polylog); (2) interval [x-E, x+E] contains O(sqrt(x)/log^2(x)) primes; (3) use CRT modular constraints, Beatty sequence approximation, and iterated logarithmic expansion to narrow candidates.

## Key Findings
- Asymptotic expansion p(n) ~ n*ln(n) + n*ln(ln(n)) - n + ... converges, but each term only adds ~1 bit of precision. After O(log n) terms, error is still O(sqrt(n)).
- CRT constraints (p(n) mod q for small q) each give ~1 bit, but computing them requires L-function zeros at O(sqrt(x)) cost (same barrier as Proposals 1-2).
- Beatty sequence approximation: p(n) is NOT a Beatty sequence; prime gaps are irregular and their distribution encodes zeta-zero information.
- Pratt certificates prove primality of a candidate but don't help locate p(n) by index.
- The interval after all cheap refinements still contains O(sqrt(x)/log^2(x)) candidate primes -- far from O(1).
- Resolving which candidate is the nth prime requires counting primes in the interval, costing O(sqrt(x)) at minimum.

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (E) -- all refinement steps either give O(1) bits cheaply or require O(sqrt(x)) work; no intermediate regime.

## One-Line Summary
Recursive interval refinement: cheap steps give O(1) bits each; narrowing to O(1) candidates still requires O(sqrt(x)) counting.
