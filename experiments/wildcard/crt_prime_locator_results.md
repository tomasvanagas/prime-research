# CRT Prime Locator: Results

**Script:** crt_prime_locator.py

## What Was Tested
Computing p(n) mod q for many small primes q via counting primes in arithmetic progressions pi(x; q, a), then reconstructing p(n) via CRT. Uses R^{-1}(n) as initial approximation.

## Key Findings
- Computing pi(x; q, a) for each residue class a mod q requires enumerating primes up to x -- cost O(x/ln(x)) per modulus
- No known way to compute pi(x; q, a) faster than pi(x) itself
- Each CRT modulus multiplies the total cost by a constant factor
- The number of moduli needed grows as O(log p(n)) to cover all bits
- Total cost: O(log(p(n)) * x^{2/3}) -- strictly worse than computing p(n) directly

## Verdict
**CLOSED**
**Failure Mode:** Circularity (C) -- computing p(n) mod q requires pi(x; q, a) which costs at least O(x^{2/3}), the same as computing pi(x).

## One-Line Summary
CRT prime locator via arithmetic progressions: each modulus costs O(x^{2/3}), multiplying total work -- no savings over direct computation.
