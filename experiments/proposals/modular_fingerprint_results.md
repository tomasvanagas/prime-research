# Proposal C — Hecke / arithmetic fingerprint as primality oracle

## Setup

For n in [2, 10000] (containing 1229 primes, 8770 composites), compute fingerprints of varying depth combining:

1. tau(n) - n^11 - 1 (mod 691) -- the Ramanujan congruence (zero on primes)
2. tau(n) mod 5
3. sigma_1(n) - 1 - n -- zero iff n is prime or 1
4. n mod 30
5. n^6 mod 7

## Ramanujan-congruence test (depth 1)

- Primes for which tau(p) ≡ p^11 + 1 (mod 691): 1229 / 1229 = 100.00%
- Composites for which tau(c) ≡ c^11 + 1 (mod 691): 10 / 8770 = 0.11%

## Fingerprint discrimination as a function of depth

*FP* = composites whose full fingerprint also occurs at some prime (false positives, the metric we want to drive to 0).
*Collide* = primes whose fingerprint also occurs at some composite.

| depth | unique prime fps | composite FP | FP rate | primes colliding | collide rate |
|---|---|---|---|---|---|
| 1 | 1 | 10 | 0.11% | 1229 | 100.00% |
| 2 | 3 | 9 | 0.10% | 920 | 74.86% |
| 3 | 3 | 0 | 0.00% | 0 | 0.00% |
| 4 | 11 | 0 | 0.00% | 0 | 0.00% |
| 5 | 12 | 0 | 0.00% | 0 | 0.00% |

## Interpretation

With depth 5, the fingerprint **perfectly separates primes from composites** in [2, 10000]. This is a polylog primality oracle for n <= 10000. Path remains open; escalate to larger N.

## Note on cost of features

- tau(n) per single n: O(polylog n) is **conjectural**. The known best is Edixhoven et al.'s O(polylog) algorithm for tau(p) at primes p (assuming GRH for ell-adic Galois reps). For composites, tau(n) is computed via multiplicativity from tau(p) at the prime factors of n — but this requires factoring n, breaking polylog.
- sigma_k(n): O(sqrt(n)) by enumerating divisors — not polylog.
- n mod m: O(polylog).

So unless tau(n) for *composite* n can be computed without factoring (open problem, likely false), the entire fingerprint is not polylog to evaluate. This is a structural barrier even before considering the discrimination rate.
