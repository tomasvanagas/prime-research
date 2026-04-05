# Ono Partition Characterization + p-adic Lifting: Results

**Script:** ono_padic_lift.py

## What Was Tested

Whether the Ono-Craig-van Ittersum (2024) partition characterization of primes can be made efficient via p-adic methods. Tested: (1) implementing the Ono characterization for small n, (2) whether it preserves congruences mod small primes, (3) modular vs exact computation cost of MacMahon partition functions M_k(n).

## Key Findings

- Partition Mobius transform T(n) = sum_{d|n} mu(n/d)*M_2(d) is an approximation, not the exact Ono criterion
- Modular computation of M_k has same operation count as exact -- just bounded values, no asymptotic improvement
- Ono characterization requires M_k(d) for ALL divisors d of n; computing M_k(d) costs O(d^2) per divisor, total O(n^{1+epsilon})
- This is MUCH WORSE than Meissel-Lehmer's O(n^{2/3})
- Ramanujan-type congruences give M_k(n) mod l for special (n,l) pairs but not general polylog computation
- p-adic lifting does not help because the partition DP has O(n^2) inherent complexity

## Verdict

**CLOSED** -- Failure Mode: Equivalence (E). The Ono characterization is structurally beautiful but computationally inferior to all known methods (O(n^{1+eps}) vs O(n^{2/3})).

## One-Line Summary

Ono partition primality characterization + p-adic lifting is computationally worse than Meissel-Lehmer -- O(n^{1+eps}) vs O(n^{2/3}).
