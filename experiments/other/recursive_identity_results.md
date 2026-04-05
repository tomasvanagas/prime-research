# Recursive Identity Exploration: Results

**Script:** recursive_identity.py

## What Was Tested
Whether a recursive identity connects pi(x) to pi(x/2), pi(x/4), etc. yielding an O(polylog) divide-and-conquer algorithm. Six approaches: doubling formula pi(2x) vs 2*pi(x), squaring formula pi(x^2), Buchstab identity inversions, Legendre generalization, Meissel-Lehmer phi recursion analysis, and binary splitting of prime counting.

## Key Findings
- Doubling: pi(2x) - 2*pi(x) = pi(2x, x) (primes in (x, 2x]), which requires sieving the interval
- Squaring: pi(x^2) involves all primes up to x^2; the correction from pi(x) is not efficiently computable
- Buchstab identity: phi(x,a) = phi(x,a-1) - phi(x/p_a, a-1) gives a binary tree of depth pi(sqrt(x)) with O(x^{2/3}) leaves
- Legendre generalization: recursion depth is O(log x) but branching factor creates exponentially many nodes
- Meissel-Lehmer: the phi recursion has O(x^{2/3}) non-trivial leaves regardless of recursion strategy
- Binary splitting: reduces to the same O(x^{2/3}) Meissel-Lehmer complexity

## Verdict
**CLOSED** -- Failure Mode: E (Equivalence)

## One-Line Summary
All recursive identities for pi(x) have O(x^{2/3}) non-trivial leaves; no divide-and-conquer achieves polylog.
