#!/usr/bin/env python3
"""
C3 — Brandt obstructions × per-bit difficulty.

Composes E5.8 (Brandt-class barrier on π(x) mod 2) with E1.3 (per-bit
difficulty boundary). Defines

    π_J(x) := bit J of π(x)        (binary representation; J=0 is LSB)
    s_J^(N) := ( π_J(0), π_J(1), ..., π_J(2^N - 1) ) ∈ {0,1}^(2^N)

and computes the bounded-Kt complexity Kt_b(s_J^(N)) for J = 0,...,N-1
and N = 3,...,7 using the 3-bit-per-op stack VM imported from
brandt_mktp/ (L_MAX=12, T_MAX=4096; saturation = 61).

Pre-stated falsifiers (see brandt_per_bit_results.md):
  F1 — uniform saturation: bounded-VM doesn't distinguish J → DUPLICATE.
  F2 — J-stratified easy/hard cut matching E1.3's 0.5N boundary → novel
       per-bit Kt-stratified family.
  F3 — anomalous transition (cut differs from E1.3 by >1 bit, or
       non-monotonic).
"""

import math
import os
import sys
import time
from typing import Dict, List, Tuple

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
BRANDT_MKTP_DIR = os.path.normpath(
    os.path.join(THIS_DIR, "..", "brandt_mktp")
)
sys.path.insert(0, BRANDT_MKTP_DIR)
from brandt_mktp import Kt_bounded, sieve_primes, L_MAX, T_MAX  # noqa: E402

INF = 4 * L_MAX + int(math.ceil(math.log2(T_MAX + 1)))  # = 61


def pi_table(N_bits: int) -> List[int]:
    """π(x) for x = 0, 1, ..., 2^N_bits - 1."""
    M = (1 << N_bits) - 1
    sieve = sieve_primes(M)
    out = [0] * (M + 1)
    count = 0
    for x in range(M + 1):
        if sieve[x]:
            count += 1
        out[x] = count
    return out


def per_bit_string(pi_vals: List[int], J: int) -> str:
    """Extract bit J of π(x) over x ∈ [0, len(pi_vals))."""
    return "".join("1" if ((v >> J) & 1) else "0" for v in pi_vals)


def baseline_kts() -> Dict[str, int]:
    """Reference Kt_b values for calibration patterns at length 16."""
    n = 16
    targets = {
        "0^16": "0" * n,
        "1^16": "1" * n,
        "(01)^8": "01" * (n // 2),
        "(10)^8": "10" * (n // 2),
    }
    return {k: Kt_bounded(s) for k, s in targets.items()}


def run_table(N_values: List[int]) -> Dict[Tuple[int, int], dict]:
    """For each N and each J ∈ [0, N-1], compute Kt_b(s_J^(N)).

    Returns dict keyed by (N, J) with fields:
        kt: bounded Kt of s_J^(N)
        saturated: kt == INF
        ones: count of 1s in s_J^(N)
        zeros: count of 0s
        len: 2^N
        identically_zero: bool (no x in [0, 2^N) has bit J set in π(x))
    """
    results: Dict[Tuple[int, int], dict] = {}
    for N in N_values:
        pivals = pi_table(N)
        n = 1 << N
        max_pi = max(pivals)
        for J in range(N):
            s = per_bit_string(pivals, J)
            ones = s.count("1")
            zeros = n - ones
            identically_zero = (max_pi >> J) == 0
            t0 = time.time()
            kt = Kt_bounded(s)
            dt = time.time() - t0
            results[(N, J)] = {
                "N": N,
                "J": J,
                "kt": kt,
                "saturated": kt == INF,
                "ones": ones,
                "zeros": zeros,
                "len": n,
                "identically_zero": identically_zero,
                "max_pi": max_pi,
                "string": s if n <= 32 else None,
                "runtime_s": dt,
            }
    return results


def threshold_J_E13(N: int) -> float:
    """E1.3's predicted easy/hard boundary at 0.5N."""
    return 0.5 * N


def classify_outcome(results: Dict[Tuple[int, int], dict]) -> str:
    """Verdict logic per the pre-stated falsifiers F1/F2/F3."""
    # Group by N.
    by_N: Dict[int, List[Tuple[int, int, bool]]] = {}
    for (N, J), info in results.items():
        by_N.setdefault(N, []).append((J, info["kt"], info["saturated"]))
    # F1 if every (N, J) saturates.
    if all(info["saturated"] for info in results.values()):
        return "F1 (collapse): every Kt_b saturates at INF=61"
    # Else check J-stratification.
    f2_satisfied_per_N = {}
    for N, rows in by_N.items():
        rows_sorted = sorted(rows)
        # 'Easy zone' J >= 0.5N (i.e., J >= ceil(N/2)).
        # 'Hard zone' J < 0.5N (J < ceil(N/2)).
        boundary = math.ceil(N / 2)
        hard_J = [(j, kt, sat) for (j, kt, sat) in rows_sorted if j < boundary]
        easy_J = [(j, kt, sat) for (j, kt, sat) in rows_sorted if j >= boundary]
        all_hard_saturated = all(sat for (_, _, sat) in hard_J)
        any_easy_compressible = any(not sat for (_, _, sat) in easy_J)
        f2_satisfied_per_N[N] = (all_hard_saturated, any_easy_compressible)
    if all(a and b for (a, b) in f2_satisfied_per_N.values()):
        return "F2 (J-stratified): hard-J all saturate, easy-J at least one compresses"
    # Otherwise, F3.
    return "F3 (anomalous): pattern doesn't match either F1 or F2"


def main():
    print("=" * 76)
    print("C3 — Brandt obstructions × per-bit difficulty (S105)")
    print("=" * 76)
    print()
    print(f"Bounded-Kt simulator: L_MAX = {L_MAX}, T_MAX = {T_MAX}, INF = {INF}")
    print()

    # Calibration.
    print("CALIBRATION — Kt_b for reference 16-bit patterns:")
    cal = baseline_kts()
    for name, kt in cal.items():
        tag = "saturated" if kt == INF else "compressed"
        print(f"  {name:8s} -> Kt_b = {kt:3d}  [{tag}]")
    print()

    # Main table.
    N_values = [3, 4, 5, 6, 7]
    print(f"PER-BIT TABLE for N ∈ {N_values}, J = 0..N-1")
    print("    N   J  | length |  ones  zeros |  Kt_b  | sat? | id-zero? | E1.3-zone")
    print("    --- ---+--------+--------------+--------+------+----------+----------")
    t0 = time.time()
    results = run_table(N_values)
    total_runtime = time.time() - t0
    for N in N_values:
        boundary = math.ceil(N / 2)
        for J in range(N):
            r = results[(N, J)]
            zone = "easy" if J >= boundary else "hard"
            sat = "Y" if r["saturated"] else " "
            iz = "Y" if r["identically_zero"] else " "
            print(
                f"    {N:3d} {J:3d}  | {r['len']:6d} |  {r['ones']:4d}   {r['zeros']:4d} |"
                f"   {r['kt']:3d}  |  {sat}   |    {iz}     | {zone}"
            )
        print()

    print(f"Total runtime: {total_runtime:.2f} s")
    print()

    # Outcome verdict.
    verdict = classify_outcome(results)
    print(f"VERDICT: {verdict}")
    print()

    # Brief structural commentary on O1-O4 for the per-bit family.
    print("STRUCTURAL ANALYSIS — O1..O4 obstructions on the per-bit family:")
    print()
    print("  O1 (oracle-dependent vs fixed). Each π_J is a fixed total")
    print("     Boolean function. Per-bit decomposition does NOT relax O1.")
    print()
    print("  O2 (self-referential Kt). π_J answers 'is bit J of π(x) set',")
    print("     not 'is z complicated'. No self-referential measure exists.")
    print()
    print("  O3 (density of hard witnesses). The per-bit FAMILY {π_J : J}")
    print("     gives N strings per N, of which the high-J ones are")
    print("     identically-zero or near-zero. The hard-J 'witnesses' are")
    print("     finite in number (≤ ⌈log₂ π(2^N)⌉) — no infinite supply of")
    print("     growing-length Kt-bounded-random prefixes.")
    print()
    print("  O4 (uniform vs circuit). Same as for π_0 — Brandt's bound")
    print("     does not produce circuit lower bounds.")
    print()
    print("Conclusion: O1, O2, O4 unaffected by per-bit decomposition;")
    print("O3 weakly restated (finite family vs single function) but the")
    print("'infinite supply of Kt-random prefixes' requirement is still")
    print("not met. The minimal weakening of O1 the spec asked about does")
    print("NOT exist within the per-bit family. C3 closes structurally")
    print("regardless of the empirical Kt_b table.")
    print()
    print("=" * 76)


if __name__ == "__main__":
    main()
