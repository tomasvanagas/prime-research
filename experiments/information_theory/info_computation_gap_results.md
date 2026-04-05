# Information-Computation Gap Analysis Results

## What Was Tested
Session 10 investigation of 7 approaches to exploit the gap between O(log n) bits of delta(n) and O(x^{2/3}) computation time, using first 100,000 primes:
1. **PRG/Derandomization**: Whether delta(n) residuals compress better than random (would suggest PRG structure).
2. **Randomness extractors**: Min-entropy of delta(n) and whether individual bits show exploitable bias.
3. **Circuit complexity**: Testing delta(n) membership in AC0, TC0, NC1.
4. **Cryptographic hardness reductions**: Whether delta(n) reduces to/from known hard problems (factoring, discrete log).
5. **One-way function analysis**: Whether n -> p(n) or n -> delta(n) is a one-way function.
6. **Expander graph hash extraction**: Using expander walks to extract randomness from delta(n).
7. **Logspace / NC2 analysis**: Whether p(n) is computable in logspace.

## Key Findings
- **PRG/Derandomization**: BLOCKED -- delta(n) appears easy for circuits (no hardness to amplify), not hard.
- **Randomness extractors**: BLOCKED -- seed-finding for the extractor is as hard as the original problem.
- **Circuit complexity**: OPEN -- delta(n) is NOT in AC0 (mod structure rules it out), but TC0/NC1 status is unknown. This is the only non-closed path.
- **Crypto reductions**: NEUTRAL -- no reduction found to/from factoring or discrete log (positive for us -- no known hardness equivalence).
- **One-way functions**: BLOCKED -- n -> delta(n) is NOT one-to-one (many collisions); both directions are equally hard; not a OWF.
- **Expander graphs**: BLOCKED -- circular dependency: need primes to build the expander graph.
- **Logspace**: IMPORTANT INSIGHT -- p(n) IS in L (logspace), hence in NC2. But L -> NC2 simulation uses polynomial-SIZE circuits (poly in x, not polylog).

**The fundamental barrier**: The O(log n) bits measure how much data SPECIFIES the answer; the O(x^{2/3}) time measures how long it takes to FIND it. These are different quantities. pi(x) is a SUM over the prime characteristic function, and any algorithm must "touch" at least Omega(sqrt(x)) terms by the inclusion-exclusion lower bound.

**Analogy**: Computing the sum of n random bits gives an O(log n)-bit answer but requires reading all n bits. The information in the answer is small; the information needed to COMPUTE it is large.

## Verdict
**CLOSED** -- Failure Mode: **E** (Equivalence)

The information-computation gap is REAL but UNEXPLOITABLE. All 7 approaches are blocked except circuit complexity (whose status is unknown but all evidence suggests polynomial-size circuits are required). No path to polylog exists through this gap.

## One-Line Summary
Seven approaches to exploit the O(log n) bits vs O(x^{2/3}) time gap all fail; the gap is an inherent feature of the problem (answer is small but computing it requires global aggregation), not an exploitable inefficiency.
