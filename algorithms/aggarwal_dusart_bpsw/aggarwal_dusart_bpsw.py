#!/usr/bin/env python3
"""
aggarwal_dusart_bpsw.py — C4 unified p_n library

Composes E6.6 (Aggarwal binary search optimality) + E6.8 (Dusart bracket of
width n) + E5.1 (BPSW primality oracle, TC^0 conditional) into one
algorithm and benchmark suite for computing p_n exactly.

Modes:
  --mode aggarwal     Pure Aggarwal binary search using pi_lucy(x), no BPSW.
  --mode bpsw_walk    Pure BPSW-walk from 2 (no Dusart, no Aggarwal).
  --mode hybrid       Aggarwal narrowing to width K, then BPSW walk in [L,R].

Falsification (pre-stated): Wall-clock timings at p(10^k) for k in {4,5,6,7}.
All three modes must agree with sympy.prime(n). The HYBRID mode is expected
to dominate AGGARWAL_PURE in wall-clock for every n (BPSW-walk over O(log n)
candidates beats one extra Lucy-DP pi(x) call). It must also dominate
BPSW_WALK for n >= 10^4 (the Aggarwal narrowing avoids walking O(n) primes).
If either domination fails the composition is non-improving and the
composition's claim of being the natural unified library is FALSIFIED.
"""

import argparse
import math
import sys
import time
from typing import Callable, Tuple


# ============================================================================
# E6.8 — Dusart bracket
# ============================================================================
#
# For n >= 6:
#   n (log n + log log n - 1)  <=  p_n  <=  n (log n + log log n).
# Width = n exactly (logarithmically narrow on the scale of p_n ~ n log n).
# Reference: Dusart 2010, 2018; literature/aggarwal_2025_analysis.md.

_SMALL_PRIMES = [None, 2, 3, 5, 7, 11]


def dusart_bounds(n: int) -> Tuple[int, int]:
    """Return integer (L, R) with p_n in [L, R] for n >= 6.
    For n < 6, return a tight 2-element bracket containing the small prime."""
    if n < 1:
        raise ValueError("n must be >= 1")
    if n < 6:
        p = _SMALL_PRIMES[n]
        return (max(p - 1, 1), p + 1)
    ln_n = math.log(n)
    ln_ln_n = math.log(ln_n)
    # Floor for L (slightly looser than Dusart, safe).
    L = int(math.floor(n * (ln_n + ln_ln_n - 1)))
    R = int(math.ceil(n * (ln_n + ln_ln_n)))
    return (L, R)


# ============================================================================
# E5.1 — BPSW primality (TC^0 conditional, verified to 2^64)
# ============================================================================
#
# BPSW = Strong-MR base 2 AND Strong-Lucas (Selfridge parameters).
# Pre-filter by trial division on small primes; reject perfect squares
# before strong Lucas (Selfridge D-search would otherwise loop forever).

_SMALL_PRIMES_TRIAL = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)


def jacobi(a: int, n: int) -> int:
    """Jacobi symbol (a/n) for odd n > 0."""
    if n <= 0 or n % 2 == 0:
        raise ValueError("Jacobi requires odd positive n")
    a %= n
    result = 1
    while a != 0:
        while a % 2 == 0:
            a //= 2
            r = n & 7
            if r == 3 or r == 5:
                result = -result
        a, n = n, a
        if (a & 3) == 3 and (n & 3) == 3:
            result = -result
        a %= n
    return result if n == 1 else 0


def miller_rabin_strong(n: int, a: int) -> bool:
    """Strong Miller-Rabin probable-prime test for given base a; n odd > 2.
    Returns True iff n passes (could be prime or strong-pseudoprime to base a)."""
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    x = pow(a, d, n)
    if x == 1 or x == n - 1:
        return True
    for _ in range(s - 1):
        x = (x * x) % n
        if x == n - 1:
            return True
    return False


def is_perfect_square(n: int) -> bool:
    if n < 0:
        return False
    r = math.isqrt(n)
    return r * r == n


def strong_lucas_prp(n: int) -> bool:
    """Strong Lucas probable-prime test (Selfridge parameters).
    Caller must reject perfect squares first."""
    # Selfridge D-search: D in {5, -7, 9, -11, 13, ...} until Jacobi(D|n) = -1.
    D = 5
    sign = 1
    while True:
        a = D % n
        j = jacobi(a, n)
        if j == -1:
            break
        if j == 0:
            # gcd(D, n) > 1; either |D| == n (so n is small prime) or n is composite.
            return abs(D) == n
        # Next candidate: alternate sign and bump magnitude by 2.
        sign = -sign
        D = sign * (abs(D) + 2)

    P = 1
    Q = (1 - D) // 4

    # n + 1 = d * 2^s with d odd.
    d = n + 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1

    # Compute (U_d, V_d, Q^d) mod n by binary expansion of d.
    inv2 = pow(2, -1, n)
    U, V = 1, P
    Qk = Q  # Q^k where k=1
    bits = bin(d)[3:]  # skip "0b" and leading "1"
    # current k = 1 (the leading bit handled implicitly)
    for bit in bits:
        # Double: (U, V, Q^k) -> (U_{2k}, V_{2k}, Q^{2k})
        U = (U * V) % n
        V = (V * V - 2 * Qk) % n
        Qk = (Qk * Qk) % n
        if bit == '1':
            # Add: (U, V) -> ((P*U + V)/2 mod n, (D*U + P*V)/2 mod n)
            U_new = ((P * U + V) * inv2) % n
            V_new = ((D * U + P * V) * inv2) % n
            U, V = U_new, V_new
            Qk = (Qk * Q) % n  # Q^{k+1}

    # U_d == 0 mod n  =>  Lucas pseudoprime.
    if U == 0:
        return True
    # Strong condition: V_{d * 2^r} == 0 mod n for some 0 <= r < s.
    if V == 0:
        return True
    for _ in range(s - 1):
        V = (V * V - 2 * Qk) % n
        Qk = (Qk * Qk) % n
        if V == 0:
            return True
    return False


def is_bpsw(n: int) -> bool:
    """Baillie-PSW primality test.
    Verified deterministic for n < 2^64 (E5.1, conditional above)."""
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0:
        return False
    for p in _SMALL_PRIMES_TRIAL:
        if n == p:
            return True
        if n % p == 0:
            return False
    if is_perfect_square(n):
        return False
    if not miller_rabin_strong(n, 2):
        return False
    return strong_lucas_prp(n)


# ============================================================================
# pi(x) oracle — Lucy_Hedgehog DP, O(x^{2/3}) time, O(sqrt x) space
# ============================================================================
#
# Used as the practical pi(x) oracle for Aggarwal binary search.
# In Aggarwal 2025 the asymptotic claim O(sqrt(n) (log n)^4) requires HKM
# (O~(sqrt x)) which is heavyweight to implement; Lucy DP is the practical
# state-of-the-art for x reachable in this benchmark (x up to ~10^9).
#
# See literature/primecount_analysis.md and v10_c_accelerated.py for the
# C-accelerated version of essentially this code.

def pi_lucy(x: int) -> int:
    """Lucy_Hedgehog DP for pi(x). Pure Python."""
    if x < 2:
        return 0
    sqrtx = math.isqrt(x)
    while (sqrtx + 1) * (sqrtx + 1) <= x:
        sqrtx += 1
    while sqrtx * sqrtx > x:
        sqrtx -= 1
    # small[v] = (after sieve) pi(v) for v in {1, .., sqrtx}.
    small = [v - 1 if v >= 1 else 0 for v in range(sqrtx + 1)]
    # large[i] = (after sieve) pi(x // i) for i in {1, .., sqrtx}.
    large = [0] * (sqrtx + 1)
    for i in range(1, sqrtx + 1):
        large[i] = (x // i) - 1
    for p in range(2, sqrtx + 1):
        if small[p] == small[p - 1]:
            continue  # p composite
        sp = small[p - 1]
        p2 = p * p
        end = min(sqrtx, x // p2)
        for i in range(1, end + 1):
            d = i * p
            if d <= sqrtx:
                large[i] -= large[d] - sp
            else:
                large[i] -= small[x // d] - sp
        for v in range(sqrtx, p2 - 1, -1):
            small[v] -= small[v // p] - sp
    return large[1]


# ============================================================================
# E6.6 — Aggarwal binary search
# ============================================================================
#
# Find smallest x in [L, R] with pi(x) >= n. By definition that x is p_n.
# Maintains:  pi(R)     >= n
#             pi(L - 1) <  n
# with initial bracket from Dusart, both invariants hold.

def aggarwal_search(n: int,
                    pi_oracle: Callable[[int], int] = pi_lucy,
                    bracket: Tuple[int, int] = None,
                    counter: dict = None) -> int:
    if bracket is None:
        L, R = dusart_bounds(n)
    else:
        L, R = bracket
    if counter is None:
        counter = {}
    counter.setdefault('pi_calls', 0)
    while L < R:
        mid = (L + R) // 2
        c = pi_oracle(mid)
        counter['pi_calls'] += 1
        if c >= n:
            R = mid
        else:
            L = mid + 1
    return L


# ============================================================================
# BPSW walk
# ============================================================================
#
# Walk forward, count primes via BPSW.  Used in two ways:
#   (a) bpsw_walk_only(n): start at 2, count to n. Pure baseline.
#   (b) hybrid: start at L (final Aggarwal bracket lower endpoint), count
#       until we hit the residual = n - pi(L - 1)-th prime in [L, R].

def bpsw_walk_from(start: int, count_target: int,
                   counter: dict = None) -> int:
    """Return the count_target-th prime >= start, found via BPSW."""
    if counter is None:
        counter = {}
    counter.setdefault('bpsw_calls', 0)
    if count_target <= 0:
        raise ValueError("count_target must be positive")
    n = start
    found = 0
    if n <= 2:
        counter['bpsw_calls'] += 1
        if is_bpsw(2):
            found = 1
            if found == count_target:
                return 2
        n = 3
    if n % 2 == 0:
        n += 1
    while True:
        counter['bpsw_calls'] += 1
        if is_bpsw(n):
            found += 1
            if found == count_target:
                return n
        n += 2


def bpsw_walk_only(n: int, counter: dict = None) -> int:
    """Pure BPSW baseline: count primes from 2."""
    return bpsw_walk_from(2, n, counter)


# ============================================================================
# Hybrid (the C4 composition)
# ============================================================================
#
# 1. Dusart bounds (E6.8) -> initial bracket [L0, R0] of width n.
# 2. Aggarwal binary search (E6.6) using pi_lucy until R - L <= K.
# 3. One pi(L - 1) call to get residual r = n - pi(L - 1).
# 4. BPSW walk (E5.1) from L for the r-th prime in [L, ...].
#
# Conditional propagation note: BPSW correctness is unconditional for
# n < 2^64 (S64; verified empirically). For p_n < 2^64 (n < ~4.4e17),
# this algorithm is unconditional. For larger n, a single BPSW pseudoprime
# in the final bracket would shift the answer by ONE prime; conditionality
# is 1-to-1, no amplification through Aggarwal's wrapper.

def hybrid_pn(n: int, K: int = 64,
              pi_oracle: Callable[[int], int] = pi_lucy,
              counter: dict = None) -> int:
    if counter is None:
        counter = {}
    counter.setdefault('pi_calls', 0)
    counter.setdefault('bpsw_calls', 0)

    L, R = dusart_bounds(n)

    # Step 2: Aggarwal narrowing.
    while R - L > K:
        mid = (L + R) // 2
        c = pi_oracle(mid)
        counter['pi_calls'] += 1
        if c >= n:
            R = mid
        else:
            L = mid + 1
    if L == R:
        return L  # Aggarwal already converged before BPSW kicks in.

    # Step 3: residual.
    if L <= 1:
        pi_left = 0
    else:
        pi_left = pi_oracle(L - 1)
        counter['pi_calls'] += 1
    residual = n - pi_left
    if residual <= 0:
        # Should not happen: invariant pi(L-1) < n is maintained by binary
        # search. Fall back to sweeping from start of bracket.
        raise RuntimeError(f"hybrid invariant violation: pi(L-1)={pi_left} >= n={n}")

    # Step 4: BPSW walk in [L, R].
    return bpsw_walk_from(L, residual, counter)


# ============================================================================
# Verification + benchmarks
# ============================================================================

def _format_counter(counter: dict) -> str:
    parts = []
    for k in ('pi_calls', 'bpsw_calls'):
        if k in counter:
            parts.append(f"{k}={counter[k]}")
    return " ".join(parts)


def run_correctness_smoke():
    """Quick correctness smoke test against sympy.prime for small n."""
    try:
        from sympy import prime as sym_prime
    except ImportError:
        print("[WARN] sympy unavailable; skipping correctness smoke test.")
        return
    test_cases = [1, 2, 3, 5, 6, 10, 25, 50, 100, 168, 169, 500, 1000]
    print("Correctness smoke test (vs sympy.prime):")
    for n in test_cases:
        truth = sym_prime(n)
        results = {
            'agg':    aggarwal_search(n),
            'bpsw':   bpsw_walk_only(n),
            'hybrid': hybrid_pn(n),
        }
        ok = all(v == truth for v in results.values())
        flag = "OK " if ok else "FAIL"
        print(f"  n={n:6d}  truth={truth:8d}  agg={results['agg']:8d}  "
              f"bpsw={results['bpsw']:8d}  hybrid={results['hybrid']:8d}  [{flag}]")
        if not ok:
            print(f"    *** MISMATCH: {results} vs {truth}")
            sys.exit(1)


def run_benchmark(ns, modes=('agg', 'bpsw', 'hybrid'), K=64,
                  hybrid_K_sweep=None):
    """Wall-clock benchmarks at the given n-values.
    Verifies all results agree with sympy.prime when available."""
    try:
        from sympy import prime as sym_prime
        have_sympy = True
    except ImportError:
        have_sympy = False
        print("[WARN] sympy unavailable; agreement check skipped.")

    print(f"\nBenchmark (K={K}):")
    header = f"{'n':>10s} {'truth':>14s}"
    for m in modes:
        header += f"  {m:>10s}_t  {m:>10s}_n"
    header += "  agree"
    print(header)

    rows = []
    for n in ns:
        truth = sym_prime(n) if have_sympy else None
        result = {'n': n, 'truth': truth}
        agree = True
        for m in modes:
            counter = {}
            t0 = time.perf_counter()
            if m == 'agg':
                pn = aggarwal_search(n, counter=counter)
            elif m == 'bpsw':
                pn = bpsw_walk_only(n, counter=counter)
            elif m == 'hybrid':
                pn = hybrid_pn(n, K=K, counter=counter)
            else:
                raise ValueError(f"unknown mode {m!r}")
            t1 = time.perf_counter()
            result[m + '_pn'] = pn
            result[m + '_time'] = t1 - t0
            result[m + '_counter'] = counter
            if have_sympy and pn != truth:
                agree = False
        result['agree'] = agree
        rows.append(result)

        # Print row
        line = f"{n:10d}"
        line += f" {truth if truth is not None else 'NA':>14}"
        for m in modes:
            line += f"  {result[m + '_time']:9.3f}s  {result[m + '_pn']:>10d}"
        line += f"  {'YES' if agree else '*** NO ***'}"
        print(line)
        # operation counts
        ops_line = " " * 25
        for m in modes:
            ops_line += f"  {_format_counter(result[m + '_counter']):>22s}"
        print(ops_line)

    # Optional K-sweep: hybrid only
    if hybrid_K_sweep:
        print("\nHybrid K-sweep:")
        n_target = hybrid_K_sweep['n']
        Ks = hybrid_K_sweep['Ks']
        truth = sym_prime(n_target) if have_sympy else None
        print(f"  n = {n_target}, truth = {truth}")
        for K_ in Ks:
            counter = {}
            t0 = time.perf_counter()
            pn = hybrid_pn(n_target, K=K_, counter=counter)
            t1 = time.perf_counter()
            ok = "OK" if (truth is None or pn == truth) else "FAIL"
            print(f"  K={K_:8d}  time={t1-t0:7.3f}s  pi_calls={counter.get('pi_calls', 0):3d}  "
                  f"bpsw_calls={counter.get('bpsw_calls', 0):8d}  pn={pn}  [{ok}]")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--smoke', action='store_true',
                    help='Run small-n correctness smoke test only.')
    ap.add_argument('--bench', action='store_true',
                    help='Run wall-clock benchmark.')
    ap.add_argument('--ns', type=int, nargs='*', default=None,
                    help='Override n values for --bench.')
    ap.add_argument('--K', type=int, default=64,
                    help='Bracket-width threshold for hybrid (default 64).')
    ap.add_argument('--Ksweep', type=int, default=None,
                    help='If set, sweep K over {4, 16, 64, 256, 1024, 4096} '
                         'at this n.')
    ap.add_argument('--mode', choices=('agg', 'bpsw', 'hybrid'), default=None,
                    help='Run a single mode at --n (no benchmarking framework).')
    ap.add_argument('--n', type=int, default=None)
    args = ap.parse_args()

    if args.mode is not None:
        if args.n is None:
            print("--mode requires --n", file=sys.stderr)
            sys.exit(2)
        counter = {}
        t0 = time.perf_counter()
        if args.mode == 'agg':
            pn = aggarwal_search(args.n, counter=counter)
        elif args.mode == 'bpsw':
            pn = bpsw_walk_only(args.n, counter=counter)
        elif args.mode == 'hybrid':
            pn = hybrid_pn(args.n, K=args.K, counter=counter)
        t1 = time.perf_counter()
        print(f"mode={args.mode} n={args.n}  pn={pn}  time={t1-t0:.3f}s  "
              f"{_format_counter(counter)}")
        return

    if args.smoke or not (args.bench or args.Ksweep):
        run_correctness_smoke()

    if args.bench:
        if args.ns:
            ns = args.ns
        else:
            ns = [10, 100, 1000, 10_000, 100_000, 1_000_000]
        sweep = None
        if args.Ksweep:
            sweep = {'n': args.Ksweep,
                     'Ks': [4, 16, 64, 256, 1024, 4096, 16384, 65536]}
        run_benchmark(ns, K=args.K, hybrid_K_sweep=sweep)


if __name__ == '__main__':
    main()
