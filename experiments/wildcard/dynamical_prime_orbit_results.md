# Dynamical Prime Orbit: Results

**Script:** dynamical_prime_orbit.py

## What Was Tested
Four novel dynamical approaches: (1) linear dynamical systems over Z/MZ whose orbits encode p(n), (2) sieve as GF(2) matrix product with fast-forwarding, (3) polynomial recurrence hunting via Berlekamp-Massey on p(n) mod m, (4) continued fraction structure of p(n)/n and Stern-Brocot navigation.

## Key Findings
- Linear systems over Z/MZ: p(n) mod M does not satisfy any linear recurrence of dimension <= 6 for any tested M; verification on held-out data fails for all (M, d) pairs
- Sieve as GF(2) matrix product: the product of per-prime projection matrices has no exploitable structure; matrix size grows as primorial
- Berlekamp-Massey on p(n) mod m: linear complexity is N/2 for all tested fields -- maximal, confirming LFSR results from Session 24
- Continued fraction of p(n)/n: partial quotients show no pattern; CF length grows linearly with n

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (E) -- linear recurrences over Z/MZ have maximal complexity; GF(2) matrix products and CF structure offer no compression.

## One-Line Summary
Dynamical orbit approaches (linear mod M, GF(2) product, Berlekamp-Massey, CF): all show maximal complexity, no fast-forwarding possible.
