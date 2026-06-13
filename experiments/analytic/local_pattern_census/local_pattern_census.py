#!/usr/bin/env python3
"""
local_pattern_census.py — S491 continuation: the local pattern census
D_k(x) of the primes, against the exact Hardy-Littlewood admissible set.

Object (from the S491 width-spectrum cliff): D_k(x) = number of DISTINCT
exact occupancy patterns among aligned k-windows of the prime indicator,
chi_P[mk .. (m+1)k), m < x/k. The width-dichotomy measurement found
D_16(2^24) = 107 vs 4707 for matched random — a 44x rigidity. This
script explains the 107 exactly:

  A pattern S subset {0..k-1} is ADMISSIBLE for aligned k-windows iff
    (i)  for every prime q | k: no o in S with q | o
         (else km+o is divisible by q for every window m), and
    (ii) for every prime q <= k with q !| k: the residues
         {o mod q : o in S} do not cover all of Z_q
         (else every window has a q-divisible element at an occupied
         offset).
  Under the prime k-tuple conjecture every admissible pattern occurs as
  an exact pattern (prime at S, composite off S) for infinitely many m;
  inadmissible patterns occur at most for finitely many small-m windows
  (where km+o IS the obstructing prime q itself). Prediction:

      D_k(x) -> A_k + E_k   (A_k = #admissible, E_k = #small-x
                              exceptional patterns, a finite set)

This is a falsifiable finite-x census of the k-tuple conjecture's
qualitative content, and it pins the asymptotic height of the
width-profile cliff at exactly A_k + E_k.

What would falsify: D_k(x) exceeding A_k + (verified exceptional count)
at large x (an inadmissible pattern occurring at large m — would
CONTRADICT plain divisibility, so really a bug check); or persistent
non-convergence (admissible patterns of small weight |S| never
appearing — would be evidence against k-tuple uniformity).

Usage: python3 local_pattern_census.py [--k 16] [--nmax 28]
"""

import argparse

import numpy as np


def sieve(limit):
    s = np.ones(limit, dtype=bool)
    s[:2] = False
    for p in range(2, int((limit - 1) ** 0.5) + 1):
        if s[p]:
            s[p * p::p] = False
    return s


def primes_upto(m):
    return [p for p in range(2, m + 1) if all(p % d for d in
                                              range(2, int(p ** 0.5) + 1))]


def admissible_patterns(k):
    """Exact enumeration of admissible aligned-k patterns."""
    div_primes = [q for q in primes_upto(k) if k % q == 0]
    other_primes = [q for q in primes_upto(k) if k % q != 0]
    # condition (i) restricts S to offsets coprime-compatible:
    # o with q|o for q|k can never be occupied
    allowed = [o for o in range(k)
               if all(o % q != 0 for q in div_primes)]
    adm = []
    for mask in range(1 << len(allowed)):
        S = [allowed[i] for i in range(len(allowed)) if (mask >> i) & 1]
        ok = True
        for q in other_primes:
            if len({o % q for o in S}) == q:
                ok = False
                break
        if ok:
            adm.append(frozenset(S))
    return adm, allowed


def realized_patterns(chi, k):
    """Distinct aligned-k patterns as frozensets, with first occurrence m."""
    nwin = len(chi) // k
    win = chi[:nwin * k].reshape(nwin, k)
    codes = win.dot(1 << np.arange(k, dtype=np.int64))
    first = {}
    uniq, idx = np.unique(codes, return_index=True)
    for c, i in zip(uniq, idx):
        S = frozenset(int(o) for o in range(k) if (int(c) >> o) & 1)
        first[S] = int(i)  # window index m of first occurrence
    return first


def count_admissible_dfs(k):
    """A_k by DFS with pruning. Admissibility is downward-closed
    (a subset of an admissible set is admissible), so a partial subset
    that already covers all residues mod some prime kills its whole
    supertree. Node count ~ O(A_k * k) — reaches k ~ 60."""
    div_primes = [q for q in primes_upto(k) if k % q == 0]
    other_primes = [q for q in primes_upto(k) if k % q != 0]
    allowed = [o for o in range(k)
               if all(o % q != 0 for q in div_primes)]
    full = {q: (1 << q) - 1 for q in other_primes}
    obits = [{q: 1 << (o % q) for q in other_primes} for o in allowed]

    count = 0
    # iterative DFS over (index, coverage) — coverage as tuple of masks
    qlist = other_primes
    stack = [(0, tuple(0 for _ in qlist))]
    while stack:
        idx, cov = stack.pop()
        if idx == len(allowed):
            count += 1
            continue
        # exclude offset idx
        stack.append((idx + 1, cov))
        # include offset idx (prune if some prime becomes fully covered)
        newcov = []
        dead = False
        for qi, q in enumerate(qlist):
            c = cov[qi] | obits[idx][q]
            if c == full[q]:
                dead = True
                break
            newcov.append(c)
        if not dead:
            stack.append((idx + 1, tuple(newcov)))
    return count


def entropy_scan(kmax):
    """A_k for k = 8, 12, ..., kmax; test rate models
    log2(A_k)/k -> const  vs  log2(A_k)/k ~ c/ln k (Mertens shape)."""
    import math
    print(f"{'k':>4} {'A_k':>12} {'log2A/k':>9} {'(log2A/k)*ln(k)':>16}")
    for k in range(8, kmax + 1, 4):
        a = count_admissible_dfs(k)
        rate = math.log2(a) / k
        print(f"{k:>4} {a:>12} {rate:>9.4f} {rate*math.log(k):>16.4f}")


def main(k, nmax):
    adm, allowed = admissible_patterns(k)
    adm_set = set(adm)
    print(f"k = {k}: allowed offsets (cond i): {allowed}")
    print(f"A_{k} (admissible patterns, exact): {len(adm)}")
    wmax = max((len(s) for s in adm), default=0)
    from collections import Counter
    cnt = Counter(len(s) for s in adm)
    print(f"admissible by weight |S|: "
          f"{dict(sorted(cnt.items()))} (max weight {wmax})")

    print(f"\n{'x':>12} {'D_k(x)':>8} {'adm.':>6} {'except.':>8} "
          f"{'missing adm.':>13}")
    chi = sieve(1 << nmax)
    for n in range(16, nmax + 1, 2):
        first = realized_patterns(chi[:1 << n], k)
        realized = set(first)
        n_adm = len(realized & adm_set)
        exc = realized - adm_set
        missing = adm_set - realized
        print(f"{1 << n:>12} {len(realized):>8} {n_adm:>6} "
              f"{len(exc):>8} {len(missing):>13}")
        if n == nmax:
            print(f"\nexceptional (inadmissible, small-x) patterns at "
                  f"x = 2^{nmax}: {len(exc)}")
            for S in sorted(exc, key=lambda s: (first[s], sorted(s))):
                m = first[S]
                print(f"  m={m:>4} window [{m*k},{(m+1)*k}): "
                      f"offsets {sorted(S)}")
            print(f"\nadmissible not yet realized at x = 2^{nmax}: "
                  f"{len(missing)} (by weight: "
                  f"{dict(sorted(Counter(len(s) for s in missing).items()))})")
            for S in sorted(missing, key=lambda s: (len(s), sorted(s)))[:10]:
                print(f"  weight {len(S)}: offsets {sorted(S)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=16)
    ap.add_argument("--nmax", type=int, default=28)
    ap.add_argument("--entropy-scan", type=int, default=0,
                    help="compute A_k for k = 8..N step 4 and fit rate models")
    args = ap.parse_args()
    if args.entropy_scan:
        # regression guard: DFS must reproduce the enumerative counts
        assert count_admissible_dfs(8) == 13
        assert count_admissible_dfs(16) == 106
        assert count_admissible_dfs(32) == 3573
        entropy_scan(args.entropy_scan)
    else:
        main(args.k, args.nmax)
