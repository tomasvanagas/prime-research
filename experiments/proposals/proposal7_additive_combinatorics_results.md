# Proposal 7: Additive Combinatorics / Green-Tao Machinery -- Results

**Script:** proposal7_additive_combinatorics.py

## What Was Tested
Explore whether Green-Tao transference principle, GPY sieve weights, or Gowers-norm uniformity can be used to extract p(n) as a weighted sum or threshold function. Test prime distribution in residue classes mod primorials.

## Key Findings
- Primes distribute nearly uniformly across coprime residue classes mod W = primorial(p), confirming Dirichlet's theorem.
- The Green-Tao machinery proves EXISTENCE of patterns (arithmetic progressions) but doesn't provide efficient SELECTION of the nth prime.
- GPY weights are designed to detect small prime gaps, not to locate specific primes by index.
- p(n) = min{m : pi(m) >= n} is a THRESHOLD/SELECTION problem, not a counting problem. Additive combinatorics tools address counting and density, not selection.
- Extracting the nth term of a Dirichlet series (prime zeta function) without computing all previous terms is equivalent to the original problem.
- No generating-function trick avoids the pi(x) = n inversion.

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (E) -- additive combinatorics gives counting/density results, not selection; extracting the nth prime from a generating function is equivalent to computing pi(x).

## One-Line Summary
Additive combinatorics / Green-Tao: addresses counting and density, not nth-element selection; no shortcut for p(n).
