"""
Thread 11 / Slot 5 — Hardy-Littlewood singular series for Ulam slope-+/-1
diagonals, and structural identification of the LP-cover constant
c = lim_{N -> infty} LP_p(N) / LP_r(N).

Slot 4 (S487) measured LP_p / LP_r = 0.7774 +/- 0.0055 across N in
[10^3, 10^6] (six orders of magnitude). Slot 5 task: provide the
structural characterisation of the constant via HL Conjecture F on
the slope-+/-1 quadratic-prime sequences carrying 69% of LP weight.

Empirical calculus:
- LP_r(N) = sqrt(N) (random points: integer-tight axis cover)
- LP_p(N) = c * sqrt(N), c ≈ 0.776
- Ratio kappa_p / kappa_r ≈ 1.288 = "average HL boost"
  on LP-active slope-+/-1 diagonals

Each slope-+/-1 Ulam diagonal is the union of 2 quadratic sequences
of form f(k) = 4k^2 + b*k + c, one per "side" of the spiral ring.
HL Conjecture F gives the prime-counting density:
  #{k <= K : f(k) prime} ~ (C_f / 2) * K / log(f(K))
with HL singular series C_f = prod_p (1 - sigma_f(p)/p) / (1 - 1/p),
sigma_f(p) = #{k mod p : p | f(k)}.

This script:
1. enumerates slope-+/-1 Ulam diagonals at intercepts |c| <= C_max;
2. for each, finds the underlying quadratic pair (4k^2 + b*k + c);
3. computes HL singular series via prod_p truncated at p <= P_max;
4. compares the empirical prime count on the diagonal to HL prediction;
5. computes the LP-weighted average of HL singular series and compares
   to the empirical 1/0.776 = 1.288.

CLI:
  python p11_hl_singular_series.py --N 100000 --K 5 --P_max 1000

Output: per-diagonal table + LP-weighted average + theoretical c bound.
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from math import gcd

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p11_ulam_line_cover import sieve_primes, ulam_coords_array


def legendre(a: int, p: int) -> int:
    """Legendre symbol (a/p). p odd prime."""
    a = a % p
    if a == 0:
        return 0
    r = pow(a, (p - 1) // 2, p)
    return 1 if r == 1 else -1


def hl_singular_series(a: int, b: int, c: int, P_max: int = 10000) -> float:
    """HL singular series C_f for f(n) = a*n^2 + b*n + c.

    C_f = prod_p (1 - sigma_f(p) / p) / (1 - 1/p)
    where sigma_f(p) = #{n in {0, ..., p-1} : f(n) ≡ 0 (mod p)}.

    For f irreducible (D = b^2 - 4ac NOT a perfect square):
    - p odd, p ∤ a: sigma_f(p) = 1 + (D/p)
        - if (D/p) = -1: sigma = 0, factor = (1 - 0)/(1 - 1/p) = p/(p-1)
        - if (D/p) = +1: sigma = 2, factor = (1 - 2/p)/(1 - 1/p) = (p-2)/(p-1)
        - if (D/p) =  0: sigma = 1, factor = (1 - 1/p)/(1 - 1/p) = 1
    - p odd, p | a: separate case (handled by direct enumeration)
    - p = 2: enumerate directly.

    Truncated product up to P_max.
    """
    D = b * b - 4 * a * c
    # D = 0 or perfect square => reducible => infinite or trivial; skip.
    if D == 0:
        return 0.0
    if D > 0:
        s = int(math.isqrt(D))
        if s * s == D:
            return 0.0  # reducible
    sieve = bytearray(b"\x01") * (P_max + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(P_max ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = bytearray(len(sieve[i * i :: i]))
    primes = [i for i in range(P_max + 1) if sieve[i]]
    log_C = 0.0
    for p in primes:
        # Direct enumeration of sigma_f(p) for safety (handles p | a).
        sigma = 0
        for n in range(p):
            if (a * n * n + b * n + c) % p == 0:
                sigma += 1
        # factor = (1 - sigma/p) / (1 - 1/p) = (p - sigma) / (p - 1)
        if p - 1 == 0:
            continue
        factor = (p - sigma) / (p - 1)
        if factor <= 0:
            return 0.0
        log_C += math.log(factor)
    return math.exp(log_C)


def diagonal_quadratic_decomp(N: int, ulam_coords: list, line_a: int, line_b: int, line_c: int) -> list[int]:
    """Find all integers v in [1, N] on the line b*x - a*y = c."""
    on_line = []
    for v in range(1, N + 1):
        x, y = ulam_coords[v]
        if line_b * x - line_a * y == line_c:
            on_line.append(v)
    return on_line


def fit_quadratic_from_consecutive(seq: list[int]) -> tuple[int, int, int] | None:
    """Given a sorted list of integers from one quadratic side of a
    slope-+/-1 Ulam line, fit f(k) = a*k^2 + b*k + c.

    Uses 3 consecutive points starting from k=0:
      f(0) = c
      f(1) = a + b + c
      f(2) = 4a + 2b + c
    => a = (f(2) - 2*f(1) + f(0)) / 2
       b = f(1) - f(0) - a
       c = f(0)
    """
    if len(seq) < 3:
        return None
    f0, f1, f2 = seq[0], seq[1], seq[2]
    second_diff = f2 - 2 * f1 + f0
    if second_diff % 2 != 0:
        return None
    a = second_diff // 2
    b = f1 - f0 - a
    c = f0
    # Verify on a fourth point if available.
    if len(seq) >= 4:
        if a * 9 + b * 3 + c != seq[3]:
            return None
    return (a, b, c)


def split_into_quadratics(on_line: list[int]) -> list[tuple[int, int, int]]:
    """Split a sorted line-occupancy list into quadratic side-sequences.

    Heuristic (verified empirically in slot 1): on a slope-+/-1 line in
    Ulam spiral, integers come from at most 2 sides per ring, and the
    'same-side' subsequence has 2nd-difference = 8 (= 2 * 4).

    Greedy: extract longest 2nd-diff = 8 run; remove; repeat.
    """
    quadratics = []
    seq = sorted(on_line)
    while len(seq) >= 3:
        # Find first run of length >= 3 with consistent 2nd-diff = 8.
        for start in range(len(seq) - 2):
            run = [seq[start]]
            for i in range(start + 1, len(seq)):
                # Test if seq[i] continues the run with 2nd-diff = 8.
                if len(run) >= 2:
                    expected = 2 * run[-1] - run[-2] + 8
                    if seq[i] == expected:
                        run.append(seq[i])
                else:
                    run.append(seq[i])
            if len(run) >= 3:
                # Check 2nd-diff = 8 for the run.
                second_diffs = [run[k + 2] - 2 * run[k + 1] + run[k] for k in range(len(run) - 2)]
                if all(d == 8 for d in second_diffs):
                    q = fit_quadratic_from_consecutive(run)
                    if q is not None:
                        quadratics.append(q)
                        for v in run:
                            seq.remove(v)
                        break
        else:
            break
    return quadratics


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=100000,
                        help="upper bound for primes")
    parser.add_argument("--C_max", type=int, default=300,
                        help="max |intercept| for slope-+/-1 lines to enumerate")
    parser.add_argument("--P_max", type=int, default=1000,
                        help="max prime in HL singular series truncation")
    parser.add_argument("--out", type=str,
                        default="hl_singular_series_results.csv")
    args = parser.parse_args()

    print(f"# Slope-+/-1 Ulam diagonal HL analysis at N = {args.N}, |c| <= {args.C_max}")
    print(f"# HL singular series truncated at p <= {args.P_max}")

    coords = ulam_coords_array(args.N)
    primes_arr = sieve_primes(args.N)
    primes_set = set(primes_arr)
    pi_N = len(primes_arr)
    print(f"# pi(N) = {pi_N}, sqrt(N) = {math.sqrt(args.N):.2f}")

    rows = []
    for direction in [(1, -1), (1, 1)]:
        a_dir, b_dir = direction
        for c in range(-args.C_max, args.C_max + 1):
            on_line = diagonal_quadratic_decomp(args.N, coords, a_dir, b_dir, c)
            if len(on_line) < 6:
                continue
            n_total = len(on_line)
            n_prime = sum(1 for v in on_line if v in primes_set)
            quadratics = split_into_quadratics(on_line)
            if not quadratics:
                continue
            # Compute HL constants for each quadratic, then average for the line.
            hl_constants = [hl_singular_series(qa, qb, qc, args.P_max)
                            for (qa, qb, qc) in quadratics]
            hl_constants = [h for h in hl_constants if h > 0]
            if not hl_constants:
                continue
            hl_mean = sum(hl_constants) / len(hl_constants)
            hl_max = max(hl_constants)
            # HL prediction for primes on line up to N:
            # density per quadratic = C_f / 2 / log(N)
            # primes on quadratic = (C_f / 2) * K / log(N)
            # where K = number of integers from this quadratic <= N
            # roughly K ~ sqrt(N/4) = sqrt(N)/2
            # primes on line ≈ sum over quads (C_q / 2) * K_q / log(K_q^2/4)
            density_observed = n_prime / n_total if n_total else 0.0
            rows.append({
                "direction": f"({a_dir}, {b_dir})",
                "c": c,
                "n_on_line": n_total,
                "n_primes": n_prime,
                "density": density_observed,
                "n_quadratics": len(quadratics),
                "quadratics": str(quadratics),
                "hl_constants": str([round(h, 4) for h in hl_constants]),
                "hl_mean": hl_mean,
                "hl_max": hl_max,
                "predicted_density_x_logN": hl_mean / 2,
            })

    # Sort by primes-on-line descending.
    rows.sort(key=lambda r: -r["n_primes"])
    print(f"\n# {len(rows)} slope-+/-1 lines with quadratic decomposition found")
    print("\n## Top-30 by prime count")
    print(f"{'dir':>10} {'c':>5} {'n_on':>5} {'n_pr':>5} {'dens':>6} {'#q':>3} {'hl_mean':>8} {'hl_max':>8}")
    for r in rows[:30]:
        print(f"{r['direction']:>10} {r['c']:>5} {r['n_on_line']:>5} {r['n_primes']:>5} "
              f"{r['density']:>6.3f} {r['n_quadratics']:>3} "
              f"{r['hl_mean']:>8.4f} {r['hl_max']:>8.4f}")

    # Total primes on the LP-active subset (slope-+/-1 lines).
    n_primes_on_diags = sum(r["n_primes"] for r in rows)
    n_lines_with_decomp = len(rows)
    avg_hl = sum(r["hl_mean"] * r["n_primes"] for r in rows) / max(n_primes_on_diags, 1)
    avg_hl_unweighted = sum(r["hl_mean"] for r in rows) / max(n_lines_with_decomp, 1)
    log_N = math.log(args.N)
    print(f"\n## Aggregate")
    print(f"   n primes on slope-+/-1 lines (counted with multiplicity): {n_primes_on_diags}")
    print(f"   n distinct slope-+/-1 lines with quadratic decomposition: {n_lines_with_decomp}")
    print(f"   Prime-weighted mean HL singular series: {avg_hl:.4f}")
    print(f"   Line-uniform mean HL singular series:   {avg_hl_unweighted:.4f}")
    print(f"   Empirical 1/0.776 (target):             1.2887")
    print(f"   log(N):                                  {log_N:.4f}")

    out_path = os.path.join(os.path.dirname(__file__), args.out)
    with open(out_path, "w", newline="") as f:
        if rows:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
    print(f"\nCSV written: {out_path}")

    return rows, avg_hl, avg_hl_unweighted


if __name__ == "__main__":
    main()
