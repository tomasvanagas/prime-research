# Proposal 2: CRT Modular Reconstruction of delta(n) -- Results

**Script:** proposal2_crt_modular.py

## What Was Tested
Compute delta(n) = p(n) - R^{-1}(n) by finding delta(n) mod q for several small primes q via Chebyshev bias / prime race analysis, then reconstruct via CRT.

## Key Findings
- Trivial cases work: p(n) mod 2 = 1 for n >= 2; p(n) mod 6 is 1 or 5 for n >= 3.
- For non-trivial moduli, determining p(n) mod q requires computing pi(x; q, a) for each residue class a, which involves L-function zeros.
- The L-function explicit formula for pi(x; q, a) has the same sqrt(x)-scale oscillatory structure as the zeta zero sum for pi(x).
- CRT needs O(log n / log log n) moduli to reconstruct delta(n), and each modulus requires its own L-function zero sum at cost O(sqrt(x)).
- Total cost: O(sqrt(x) * log(n) / log(log(n))) -- worse than direct computation.

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (E) -- p(n) mod q requires L-function zero sums of the same complexity as the original problem.

## One-Line Summary
CRT on delta(n) mod small primes: each modular residue requires L-function zeros at O(sqrt(x)) cost, no savings.
