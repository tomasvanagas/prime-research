# Probabilistic / Randomized Approaches for Exact p(n): Results

**Script:** probabilistic_exact.py

## What Was Tested
Eight randomized approaches: Monte Carlo pi(x), Las Vegas candidate generation, random walk on primes, importance sampling, randomized binary search with probabilistic pi(x), birthday-paradox prime isolation, random sieving, and probabilistic certification.

## Key Findings
- Monte Carlo pi(x): needs O(x/ln(x)) random samples for error < 0.5, which is O(x) work
- Las Vegas candidate generation: R^{-1}(n) + random correction finds primes but cannot verify the INDEX without pi(x)
- Random walk convergence: O(sqrt(x)) steps to converge, no better than deterministic
- Importance sampling: optimal proposal distribution is 1/ln(n), but variance still requires O(x^2/ln(x)) samples
- Randomized binary search: each step needs pi(x) estimation with error < 0.5, which costs O(x^{2/3})
- Birthday-paradox: no useful collision structure in primes
- Random sieving: partial sieve reduces candidates but still O(x^{2/3}) asymptotically
- Probabilistic certificate: can verify primality in polylog but not the prime's index

## Verdict
**CLOSED** -- Failure Mode: E (Equivalence) / I (Information Loss)

## One-Line Summary
Randomness does not help: all probabilistic approaches require O(x^{2/3}) work or more for exact p(n).
