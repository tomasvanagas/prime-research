# Randomized Algorithms for pi(x): Results

**Script:** `randomized_pi.py`
**Session:** 15

## What Was Tested
Whether randomness helps compute pi(x) exactly in sublinear time. Five approaches: random sampling of zeta zeros, randomized inclusion-exclusion (probabilistic sieve), hash-based counting, RP/BPP oracle exploitation, and quantum counting.

## Key Findings
- Zeta zero sampling: variance analysis shows ALL K_min ~ x^{0.27} terms are needed; random subsampling gives O(1) variance only if M = K_min (no savings)
- Randomized I-E: massive cancellation in the sieve means subsampling introduces catastrophic errors
- Hash-based counting: hashing preserves prime density but does not help count without testing each candidate
- BPP oracle: even with a perfect BPP primality oracle, counting still requires enumerating candidates
- Quantum counting: Grover gives O(sqrt(x)) queries, which is worse than the classical O(x^{2/3}) work

## Verdict
**CLOSED**
**Failure Mode:** Information loss (randomness does not help: zeta zero variance requires all terms; sieve cancellation prevents subsampling)

## One-Line Summary
Randomized approaches for exact pi(x) all fail: zeta zero subsampling has prohibitive variance, probabilistic sieve has catastrophic cancellation.
