# Proposal 10: CRT Reconstruction of pi(x) via Modular Prime Counting -- Results

**Script:** proposal10_crt_modular_pi.py

## What Was Tested
Compute pi(x) mod m for many small moduli m via modular Legendre sieve, then reconstruct pi(x) exactly via CRT. The key question: is pi(x) mod m cheaper to compute than pi(x) exactly?

## Key Findings
- Modular Legendre sieve: phi(x, a) mod m can be computed by reducing floor(x/d) mod m, which is trivial per term.
- However, the Legendre sieve has 2^{pi(sqrt(x))} terms. Even mod m, grouping by d mod m only reduces by factor m, still exponential.
- Meissel-Lehmer computes pi(x) in O(x^{2/3}). Computing pi(x) mod m via Meissel-Lehmer still costs O(x^{2/3}) -- the modular reduction doesn't save intermediate steps.
- Need O(log x / log log x) moduli with product > x, each at O(x^{2/3}) cost.
- Total: O(x^{2/3} * log x / log log x) -- WORSE than a single exact computation.
- The Schoof analogy fails: for elliptic curves, |E(F_p)| mod l is determined by Frobenius on l-torsion (low-dimensional). For pi(x), the analogous structure involves all zeta/L-function zeros (infinite-dimensional).

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (E) -- pi(x) mod m costs the same as pi(x) exactly; no modular shortcut for Meissel-Lehmer.

## One-Line Summary
CRT on pi(x) mod m: each modular computation costs O(x^{2/3}) same as exact; Schoof analogy fails (infinite cohomology).
