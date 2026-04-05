# Modular Pi (CRT for pi(x)): Results

**Date:** 2026-04-04 (Session 4)
**Script:** modular_pi.py

## What Was Tested
Whether pi(x) mod m can be computed cheaper than pi(x) exactly, then reconstructed via CRT. Four angles: (1) parity of pi(x), (2) character sums over primes, (3) p-adic / Kubota-Leopoldt connection, (4) power sums S_k(x) = sum_{p<=x} p^k.

## Key Findings
- Angle 1 (parity): pi(x) mod 2 is AS HARD as pi(x) itself; a parity oracle gives primality testing; li(x) predicts parity at ~50% (random chance)
- Angle 2 (character sums): Cannot avoid computing pi(x); trivial character IS pi(x); non-trivial ones give O(sqrt(x)) noise
- Angle 3 (p-adic): p-adic methods encode algebraic/class-number info, not counting info; no shortcut to pi(x) mod p
- Angle 4 (power sums): S_k(x) mod m offers no shortcut; computational bottleneck (sieving) is independent of accumulation function
- CRT reconstruction adds ~330x overhead for x=10^100 even if modular computation were cheaper -- and it is not
- pi(x) mod m costs the same as pi(x) exactly: the Meissel-Lehmer recursion tree is identical

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence) -- computing pi(x) mod m visits the same recursion tree as pi(x) exact

## One-Line Summary
pi(x) mod m costs the same as pi(x) exactly; parity oracle implies primality test; CRT adds overhead without reducing complexity.
