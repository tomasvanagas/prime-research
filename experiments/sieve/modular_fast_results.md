# Modular Fast (Modular Arithmetic Approaches): Results

**Date:** 2026-04-04 (Session 5)
**Script:** modular_fast.py

## What Was Tested
Six modular arithmetic approaches to computing p(n): (1) Legendre's formula mod small primes + CRT, (2) Meissel-Mertens constant approach, (3) Euler product exact value, (4) Lucas sequences / Fibonacci entry points, (5) quadratic residue pattern reversal, (6) Bertrand postulate iterations + fast primality.

## Key Findings
- Approach 1 (Legendre + CRT): Correct but requires 2^{pi(sqrt(x))} work -- exponential in pi(sqrt(x))
- Approach 2 (Meissel-Mertens): Exact S(x) encodes primes but requires knowing them; approximate S(x) cannot recover pi(x) -- circular
- Approach 3 (Euler product): Encodes primes but computing it requires knowing them -- circular
- Approach 4 (Lucas/Fibonacci): Entry points create mapping primes->integers but inverting requires factoring large Fibonacci numbers
- Approach 5 (QR patterns): Characterize primes mod M but don't help compute p(n) without testing candidates
- Approach 6 (Bertrand + MR): Sequential Miller-Rabin is 10-100x slower than sieving; recursive halving gives O(n ln n)

## Verdict
**CLOSED**
**Failure Mode:** C (Circularity) for approaches 2-3; E (Equivalence) for approach 1; I (Information Loss) for approaches 4-6

## One-Line Summary
All 6 modular arithmetic approaches fail: Legendre+CRT is exponential, Euler product is circular, Lucas requires factoring, MR is slower than sieving.
