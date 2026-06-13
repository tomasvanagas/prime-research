#!/usr/bin/env python3
"""
automaton_width_dichotomy.py — S491 brainstorm build: the width spectrum
of the prime-counting pipeline as a VERIFICATION-cost measure.

New object. For a Boolean predicate f on n bits read MSB-first, the
per-cut Nerode width W_j(f) = #distinct subfunctions after j bits. The
S491 transfer-DP primitive shows: the multilinear extension f~ is
evaluable at an ARBITRARY field point in O(sum_j W_j) field ops GIVEN an
explicit read-once branching program. So the width profile is not just a
computation lower bound (S17/S20/S23/S28 used it that way) — it is the
EXACT price a sum-check verifier pays per point, and what separates
verifiable from unverifiable routes is whether the minimal-width BP is
CONSTRUCTIBLE without already knowing the function.

Measured here:
  1. Width profile W_j(chi_P) for n = 16..24, vs matched-density random
     control. Yields the verification exponent c in sum_j W_j ~ x^c —
     the cost of ONE-SHOT sum-check verification of pi(x) = sum chi_P(w)
     if an oracle handed the verifier the minimal OBDD.
  2. Width profile of the division relation [u = floor(v/p)] (interleaved
     bits) — the wiring of the S491 Lucy protocol. Expect max width
     p + O(1) and TIGHT (Myhill-Nerode: remainders pairwise
     distinguishable), proving the protocol's O(n*p) wiring check is
     optimal within the explicit-automaton paradigm; improving it
     requires prover-supplied carry witnesses (GKR-style), not a better
     automaton.
  3. Comparator [v >= M] width (= 3), the other wiring ingredient.

Falsification: width profile of chi_P found polylog at some cut order
(would give a polylog one-shot verifier AND contradict E2.1); division
relation width o(p) (would beat the S491 wiring without aux witnesses);
chi_P profile indistinguishable from random at all cuts (would refute
the 2^{0.79N} OBDD-size edge S20/S28).

Usage: python3 automaton_width_dichotomy.py [--nmax 24]
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


def width_profile(table_bool):
    """W_j = #distinct rows of the 2^j x 2^(n-j) reshape, j = 1..n-1."""
    n = int(np.log2(len(table_bool)))
    widths = []
    for j in range(1, n):
        rows = np.packbits(table_bool.reshape(1 << j, 1 << (n - j)), axis=1)
        widths.append(len(np.unique(rows, axis=0)))
    return widths


def division_relation_profile(n, p):
    """Width profile of [u = floor(v/p)] over interleaved bits
    (v_0 u_0 v_1 u_1 ... MSB-first)."""
    N = 1 << n
    v = np.arange(N, dtype=np.int64)
    acc = np.zeros((N, N), dtype=bool)
    acc[v, v // p] = True
    # interleave bits: index = sum_j (v_j*2^(2(n-1-j)+1) + u_j*2^(2(n-1-j)))
    table = np.zeros(1 << (2 * n), dtype=bool)
    vv, uu = np.meshgrid(v, v, indexing="ij")
    idx = np.zeros_like(vv)
    for j in range(n):
        vbit = (vv >> (n - 1 - j)) & 1
        ubit = (uu >> (n - 1 - j)) & 1
        idx |= (vbit << (2 * (n - 1 - j) + 1)) | (ubit << (2 * (n - 1 - j)))
    table[idx.ravel()] = acc.ravel()
    return width_profile(table)


def comparator_profile(n, M):
    v = np.arange(1 << n)
    return width_profile(v >= M)


def main(nmax):
    rng = np.random.default_rng(0)

    print("=== 1. chi_P width profile vs matched random ===")
    print(f"{'n':>3} {'peak_P':>9} {'argpeak':>8} {'sumW_P':>10} "
          f"{'c=log(sumW)/n':>14} {'peak_rand':>10} {'sumW_rand':>10}")
    for n in range(16, nmax + 1, 2):
        chi = sieve(1 << n)
        wp = width_profile(chi)
        rand = np.zeros(1 << n, dtype=bool)
        rand[rng.choice(np.arange(2, 1 << n), size=int(chi.sum()),
                        replace=False)] = True
        wr = width_profile(rand)
        c = np.log2(sum(wp)) / n
        print(f"{n:>3} {max(wp):>9} {int(np.argmax(wp))+1:>8} "
              f"{sum(wp):>10} {c:>14.4f} {max(wr):>10} {sum(wr):>10}")
        if n == nmax:
            print(f"    full profile  W_j(chi_P), n={n}:")
            print("    " + " ".join(str(w) for w in wp))
            print(f"    full profile  W_j(random), n={n}:")
            print("    " + " ".join(str(w) for w in wr))
            half = n // 2
            print(f"    middle-cut width: chi_P {wp[half-1]} vs "
                  f"2^(n/2-1)+2 = {(1 << (half - 1)) + 2} (E2.1 rank bound)")

    print("\n=== 2. division relation [u = floor(v/p)] width (n=8) ===")
    print(f"{'p':>4} {'max_width':>10} {'max/p':>8}  profile(first 10 cuts)")
    for p in (3, 5, 7, 11, 13):
        prof = division_relation_profile(8, p)
        print(f"{p:>4} {max(prof):>10} {max(prof)/p:>8.2f}  "
              f"{prof[:10]}")

    print("\n=== 3. comparator [v >= M] width (n=12) ===")
    for M in (37, 1000, 4091):
        prof = comparator_profile(12, M)
        print(f"  M={M:>5}: max width = {max(prof)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=24)
    args = ap.parse_args()
    main(args.nmax)
