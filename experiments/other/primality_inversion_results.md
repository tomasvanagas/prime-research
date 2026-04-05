# Deterministic Primality Testing Inversion: Results

**Script:** primality_inversion.py

## What Was Tested
Whether primality tests can be "inverted" to find the nth prime without brute-force enumeration. Five approaches: Miller-Rabin witness inversion, Fermat quotient arithmetic properties, Wilson's theorem inversion with fast factorial, constructive wheel sieve acceleration via CRT, and Legendre's formula fast evaluation.

## Key Findings
- Miller-Rabin witness inversion: the witness pattern a^d mod n encodes n's primality but cannot be inverted to find the nth prime
- Fermat quotients q_p(a) = (a^{p-1} - 1)/p have arithmetic properties but no pattern useful for prime enumeration
- Wilson's theorem: (n-1)! mod n = n-1 iff n is prime, but computing (n-1)! mod n costs O(n) multiplications
- Wheel sieve + CRT: accelerates constant factors but doesn't change O(x^{2/3}) complexity class
- Legendre's formula: pi(x) = phi(x, a) + a - 1 - P2 where phi requires O(x^{2/3}) evaluations
- All primality test inversions reduce to testing candidates one by one

## Verdict
**CLOSED** -- Failure Mode: C (Circularity)

## One-Line Summary
Primality tests cannot be inverted; finding the nth prime still requires testing or counting candidates.
