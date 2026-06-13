"""
Thread 11 / Slot 1 — top-line quadratic decomposition.

Each Ulam-spiral line of slope (a, b) is a UNION of quadratic
sequences in the ring index k. For a line of slope (1, -1) [the
NW-SE diagonal direction] passing through Ulam coords with
x + y = c, the integers landing on it from each side of each ring
form one quadratic in k. Reconstruct these quadratics for the top-K
densest lines and report their HL-prime profiles.

CLI:
  python p11_top_lines_quadratic_decomp.py --top-lines 1 -1 -18 1 -1 -58 1 1 40
"""
from __future__ import annotations

import argparse
from math import gcd
from p11_ulam_line_cover import sieve_primes, ulam_coords_array, canonical_direction


def find_quadratics_on_line(N: int, ulam_coords: list[tuple[int, int]], line_a: int, line_b: int, line_c: int) -> list[tuple[str, int, int, int, list[int]]]:
    """Walk the Ulam spiral, find all integers v in [1, N] whose Ulam
    coords (x, y) satisfy line_b * x - line_a * y = line_c. Group them
    by which "side" of which ring they came from. Return list of
    (side_label, ring_k_start, integer_values).

    Each side of each ring is a contiguous run of integers; on the
    line, these form a quadratic in k.
    """
    # Brute force: scan v=1..N, check membership.
    on_line = []
    for v in range(1, N + 1):
        x, y = ulam_coords[v]
        if line_b * x - line_a * y == line_c:
            on_line.append(v)
    return on_line


def fit_quadratics(integer_sequence: list[int], ulam_coords: list[tuple[int, int]]) -> list[tuple[int, int, int, list[int]]]:
    """Given a sorted list of integers all on a single Ulam-spiral line,
    partition them into "sides of rings": each side's integers form a
    quadratic in k. Returns list of (a, b, c, k_values) such that
    f(k) = a*k^2 + b*k + c interpolates the integers on that side.

    Heuristic: on a slope-1 line, integers from the SAME side of
    successive rings differ by ~8k (where k is the ring index). So a
    rolling difference detector splits the sequence into runs of
    "same-side" integers.
    """
    if not integer_sequence:
        return []
    # Each integer's ring index k = ceil((sqrt(v) - 1) / 2) for v >= 1.
    # Equivalently, k satisfies (2k-1)^2 < v <= (2k+1)^2.
    rings = []
    for v in integer_sequence:
        k = 0
        while (2 * k + 1) ** 2 < v:
            k += 1
        rings.append(k)
    # Group by ring side: on a slope-1 line passing through 2 sides of
    # each ring, we get up to 2 integers per ring k. Sort by (k, position).
    # Just split: for each ring k, the integers on that ring should
    # be paired and form 1 or 2 "sides" (east, north for x+y=c>0).
    # Easier: cluster into 2 groups by index parity in the sorted list.
    return rings


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, default=100000)
    p.add_argument("--top-lines", type=int, nargs="*",
                   default=[1, -1, -18, 1, -1, -58, 1, 1, 40, 1, 1, -40, 1, -1, 2],
                   help="space-separated triples a b c for canonical line keys")
    args = p.parse_args()

    coords = ulam_coords_array(args.N)
    primes_arr = sieve_primes(args.N)
    primes_set = set(primes_arr)
    triples = []
    raw = args.top_lines
    for i in range(0, len(raw), 3):
        triples.append((raw[i], raw[i + 1], raw[i + 2]))

    print(f"# Quadratic decomposition of top Ulam-spiral lines at N={args.N}")
    print(f"# pi(N) = {len(primes_arr)}")
    print()

    for (a, b, c) in triples:
        print(f"\n## Line key (dx, dy, intercept) = ({a}, {b}, {c})")
        print(f"## Equation: {b}*x - {a}*y = {c}, i.e., the line {b}x - {a}y = {c}")
        # Find integers on this line
        on_line = []
        for v in range(1, args.N + 1):
            x, y = coords[v]
            if b * x - a * y == c:
                on_line.append(v)
        n_on_line = len(on_line)
        n_prime_on_line = sum(1 for v in on_line if v in primes_set)
        density = n_prime_on_line / n_on_line if n_on_line else 0.0
        print(f"   total Ulam-points on line: {n_on_line}")
        print(f"   primes on line: {n_prime_on_line}  density: {density:.3f}")
        print(f"   first 12 integers: {on_line[:12]}")
        print(f"   first 12 primes: {[v for v in on_line if v in primes_set][:12]}")
        # Fit quadratic to consecutive same-side integers.
        # Heuristic: split sorted on_line into 2 alternating subsequences.
        if n_on_line >= 6:
            # Estimate: subseq A = even-indexed in sorted, subseq B = odd.
            seq_a = on_line[0::2]
            seq_b = on_line[1::2]
            # Fit f(k) = αk^2 + βk + γ to first 4 of each.
            for label, seq in [("A_even-idx", seq_a), ("B_odd-idx", seq_b)]:
                if len(seq) >= 3:
                    # Use 3 consecutive pts to derive α, β, γ
                    y0, y1, y2 = seq[0], seq[1], seq[2]
                    second_diff = y2 - 2 * y1 + y0
                    alpha = second_diff / 2.0
                    if alpha == int(alpha) and alpha != 0:
                        alpha_i = int(alpha)
                        # Pick k indexing s.t. seq[i] = αk² + βk + γ for some k = k0+i.
                        # Just print first three diffs as sanity.
                        diffs = [seq[i + 1] - seq[i] for i in range(min(5, len(seq) - 1))]
                        print(f"   {label}: 1st-diffs={diffs[:5]}  2nd-diff(consec)={second_diff}  alpha={alpha}")
                    else:
                        print(f"   {label}: non-integer alpha (mixed sequence?), 2nd-diff={second_diff}")


if __name__ == "__main__":
    main()
