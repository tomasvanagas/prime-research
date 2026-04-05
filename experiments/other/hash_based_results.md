# Hash/Lookup-Based Approaches: Results

**Script:** hash_based.py

## What Was Tested
Six radical ideas: perfect hash function h(n)=p(n), compressed oracle / succinct data structure, algebraic integer approach, self-correcting formula (R^{-1}(n) + correction), multi-precision approximation cascade, and probabilistic generation with proof.

## Key Findings
- Perfect hash: no known hash function maps n to p(n) without encoding the full prime table
- Compressed oracle: random access to compressed prime sequences costs at least O(sqrt(n)) by information-theoretic arguments
- Algebraic integer: primes are integers, not usefully approximated by algebraic numbers of low degree
- Self-correcting: R^{-1}(n) correction requires pi(x) evaluation to determine which candidate is the nth prime
- Multi-approximation cascade: multiple approximations with different error characteristics still cannot cancel errors below O(sqrt(x))
- Probabilistic with proof: generating candidates is easy but proving the index requires pi(x)

## Verdict
**CLOSED** -- Failure Mode: C (Circularity) / I (Information Loss)

## One-Line Summary
Hash, lookup, and cascade approaches all fail because verifying the prime index requires computing pi(x).
