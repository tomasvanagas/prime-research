"""
Thread 11 / Slot 1 — minimum line cover of primes on the Ulam spiral.

Build Ulam-spiral coordinates for n = 1..N. Restrict to primes <= N.
Compute (an upper bound for) the minimum number of straight lines
covering all prime points, via a greedy "remove-the-densest-line"
procedure. Compare to a matched-baseline: pi(N) random integers from
[1, N] mapped through the same Ulam spiral. Identify the top-K
densest lines and check whether they correspond to Hardy-Littlewood
quadratic-prime sequences.

CLI:
  python p11_ulam_line_cover.py --N 10000 --baseline-trials 5
  python p11_ulam_line_cover.py --N 100000 --bounded-K 200

Output: prints scaling table; writes per-N CSVs and top-line tables.
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import random
import sys
import time
from collections import defaultdict
from math import gcd


def sieve_primes(N: int) -> list[int]:
    if N < 2:
        return []
    sieve = bytearray(b"\x01") * (N + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(N ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = bytearray(len(sieve[i * i :: i]))
    return [i for i in range(N + 1) if sieve[i]]


def ulam_coords_array(N: int) -> list[tuple[int, int]]:
    """Return list `coords` of length N+1 where coords[n] = (x, y) for n >= 1.
    coords[0] is unused (set to (0, 0)).

    Ulam spiral convention: 1 at origin, 2 to the right, 3 above 2, then
    leftward across the top, etc. Leg lengths follow 1, 1, 2, 2, 3, 3, 4, 4, ...
    Directions cycle R, U, L, D.
    """
    coords = [(0, 0)] * (N + 1)
    if N < 1:
        return coords
    coords[1] = (0, 0)
    if N == 1:
        return coords
    directions = [(1, 0), (0, 1), (-1, 0), (0, -1)]  # R, U, L, D
    dir_idx = 0
    leg_length = 1
    n = 1
    x, y = 0, 0
    while n < N:
        for _ in range(2):
            dx, dy = directions[dir_idx]
            for _ in range(leg_length):
                x += dx
                y += dy
                n += 1
                coords[n] = (x, y)
                if n >= N:
                    return coords
            dir_idx = (dir_idx + 1) % 4
        leg_length += 1
    return coords


def canonical_direction(dx: int, dy: int) -> tuple[int, int]:
    """Canonicalise direction so it's unique per line: gcd-reduce, fix sign."""
    g = gcd(abs(dx), abs(dy))
    if g == 0:
        return (0, 0)
    dx //= g
    dy //= g
    # Canonical: positive dx, or (dx=0, dy=1)
    if dx < 0 or (dx == 0 and dy < 0):
        dx, dy = -dx, -dy
    return (dx, dy)


def all_lines_through_pairs(points: list[tuple[int, int]]) -> dict[tuple[int, int, int], list[int]]:
    """For each pair of points compute the canonical line they define,
    accumulate the set of point-indices on each line.

    Returns: line_key -> sorted list of indices on that line.
    line_key = (dx_canonical, dy_canonical, intercept) where
    intercept c = dy * x - dx * y.

    Cost: O(M^2) where M = len(points).
    """
    n_points = len(points)
    line_indices: dict[tuple[int, int, int], set] = defaultdict(set)
    for i in range(n_points):
        x1, y1 = points[i]
        for j in range(i + 1, n_points):
            x2, y2 = points[j]
            dx = x2 - x1
            dy = y2 - y1
            cdx, cdy = canonical_direction(dx, dy)
            c = cdy * x1 - cdx * y1
            key = (cdx, cdy, c)
            line_indices[key].add(i)
            line_indices[key].add(j)
    return {k: sorted(v) for k, v in line_indices.items()}


def greedy_line_cover_from_lines(line_indices: dict, n_points: int) -> tuple[int, list]:
    """Greedy: pick the line covering the most uncovered points; remove
    those points; repeat. Singletons (uncovered points after no line of
    size >= 2 remains) each take a line of their own.

    Uses heap-based max-find with lazy deletion for O(L log L + I) cost
    where L is # lines and I is total point-line incidences.

    Returns: (total_lines_used, list_of_chosen_lines).
    """
    import heapq
    # Build inverse map: for each point, set of lines containing it.
    point_to_lines: dict[int, set] = defaultdict(set)
    for line_key, idxs in line_indices.items():
        for i in idxs:
            point_to_lines[i].add(line_key)
    # Mutable line membership.
    line_remaining: dict[tuple, set] = {k: set(v) for k, v in line_indices.items()}
    # Max-heap by negative size.
    heap: list = []
    for line_key, members in line_remaining.items():
        if len(members) >= 2:
            heapq.heappush(heap, (-len(members), line_key))
    covered = set()
    chosen = []
    while heap:
        neg_sz, line_key = heapq.heappop(heap)
        if line_key not in line_remaining:
            continue
        actual_sz = len(line_remaining[line_key])
        if actual_sz != -neg_sz:
            # stale entry — re-push if still big enough
            if actual_sz >= 2:
                heapq.heappush(heap, (-actual_sz, line_key))
            continue
        if actual_sz < 2:
            continue
        members = line_remaining.pop(line_key)
        chosen.append((line_key, sorted(members)))
        covered.update(members)
        for i in members:
            for other in list(point_to_lines[i]):
                if other == line_key:
                    continue
                if other in line_remaining:
                    line_remaining[other].discard(i)
                    new_sz = len(line_remaining[other])
                    if new_sz < 2:
                        del line_remaining[other]
                    else:
                        heapq.heappush(heap, (-new_sz, other))
    singleton_count = n_points - len(covered)
    total_lines = len(chosen) + singleton_count
    return total_lines, chosen


def random_baseline_indices(N: int, k: int, rng: random.Random) -> list[int]:
    """Random sample of k integers from [1, N] without replacement, sorted."""
    if k >= N:
        return list(range(1, N + 1))
    return sorted(rng.sample(range(1, N + 1), k))


# ---------------------------------------------------------------------------
# Hardy-Littlewood quadratic-prime sequences (well-known prime-rich quadratics)
# ---------------------------------------------------------------------------
HL_QUADRATICS = [
    # (label, a, b, c) for f(n) = a*n^2 + b*n + c
    # The canonical class-number-1 Heegner polynomials (n^2 + n + p):
    ("Euler-41",     1,  1, 41),    # Heegner d=-163; prime for n=0..39
    ("Heegner-17",   1,  1, 17),    # Heegner d=-67;  prime for n=0..15
    ("Heegner-11",   1,  1, 11),    # Heegner d=-43;  prime for n=0..9
    ("Heegner-5",    1,  1,  5),    # Heegner d=-19;  prime for n=0..3
    # Legendre prime-rich:
    ("Legendre-29",  2,  0, 29),    # 2n^2 + 29 prime for n=0..28
    # Other variants:
    ("Plus-41",      1,  0, 41),    # n^2 + 41
    ("Euler-41-NEG", 1, -1, 41),    # n^2 - n + 41 (same primes as Euler-41 but reversed)
    # Controls:
    ("SqrTriv",      1,  0,  0),    # n^2 (corner sequence; not prime-rich)
    ("OddSqr",       4,  4,  1),    # (2n+1)^2 = 4n^2 + 4n + 1
]


def quadratic_to_ulam_line(a: int, b: int, c: int, N: int, ulam: list[tuple[int, int]]) -> tuple[tuple[int, int, int] | None, int, int, list[tuple[int, int]]]:
    """Plot first few values f(n) = an^2 + bn + c on the Ulam spiral; check
    if they are collinear. Return canonical line key + count of values <=N
    + collinear count + list of points.

    Skips duplicate values (in case f(n_1) == f(n_2)).
    """
    seen_values = set()
    pts = []
    for n_val in range(0, max(500, int(N ** 0.5) + 50)):
        v = a * n_val * n_val + b * n_val + c
        if v < 1 or v > N:
            continue
        if v in seen_values:
            continue
        seen_values.add(v)
        pts.append(ulam[v])
    if len(pts) < 2:
        return None, len(pts), 0, pts
    # Find the first pair of DISTINCT points to determine the candidate line.
    line_key = None
    intercept = None
    cdx = cdy = None
    for k in range(1, len(pts)):
        x0, y0 = pts[0]
        x1, y1 = pts[k]
        if (x0, y0) != (x1, y1):
            dx = x1 - x0
            dy = y1 - y0
            cdx, cdy = canonical_direction(dx, dy)
            if (cdx, cdy) == (0, 0):
                continue
            intercept = cdy * x0 - cdx * y0
            break
    if cdx is None:
        return None, len(pts), 0, pts
    collinear_count = 0
    for x, y in pts:
        if cdy * x - cdx * y == intercept:
            collinear_count += 1
    final_key = (cdx, cdy, intercept) if collinear_count == len(pts) else None
    return final_key, len(pts), collinear_count, pts


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, default=10000)
    p.add_argument("--baseline-trials", type=int, default=3)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--top-k-lines", type=int, default=20)
    p.add_argument("--out-dir", type=str, default=".")
    args = p.parse_args()

    rng = random.Random(args.seed)
    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    print(f"# Thread 11 / Slot 1 — Ulam-spiral line cover, N={args.N}")
    sys.stdout.flush()

    # Step 1: Ulam spiral coords + primes <= N
    t0 = time.time()
    coords = ulam_coords_array(args.N)
    primes = sieve_primes(args.N)
    pi_N = len(primes)
    spread_x = max(coords[i][0] for i in range(1, args.N + 1)) - min(coords[i][0] for i in range(1, args.N + 1))
    spread_y = max(coords[i][1] for i in range(1, args.N + 1)) - min(coords[i][1] for i in range(1, args.N + 1))
    print(f"# pi(N)={pi_N}, spiral spread = {spread_x+1} x {spread_y+1}, time {time.time()-t0:.2f}s")

    # Step 2: prime points
    t0 = time.time()
    prime_points = [coords[p] for p in primes]
    line_dict = all_lines_through_pairs(prime_points)
    print(f"# distinct lines through >=2 primes: {len(line_dict)}, time {time.time()-t0:.2f}s")

    # Step 3: greedy line cover
    t0 = time.time()
    L_primes, chosen_primes = greedy_line_cover_from_lines(line_dict, len(prime_points))
    print(f"# L_Ulam_greedy(primes, N={args.N}) = {L_primes}, time {time.time()-t0:.2f}s")
    print(f"#   chosen multi-lines: {len(chosen_primes)}, singletons: {L_primes - len(chosen_primes)}")

    # Step 4: top-K densest lines
    chosen_sorted = sorted(chosen_primes, key=lambda x: -len(x[1]))[: args.top_k_lines]
    print(f"\n# Top-{args.top_k_lines} densest greedy-chosen lines on primes (Ulam):")
    print(f"# {'rank':>4} {'dx':>4} {'dy':>4} {'intercept':>10} {'count':>6}")
    top_lines_for_csv = []
    for r, ((dx, dy, c), members) in enumerate(chosen_sorted, 1):
        print(f"# {r:>4} {dx:>4} {dy:>4} {c:>10} {len(members):>6}")
        top_lines_for_csv.append((r, dx, dy, c, len(members)))

    # Step 4b: HL quadratic correspondence
    print(f"\n# Hardy-Littlewood quadratic-prime sequences vs greedy lines:")
    for name, a, b, cc in HL_QUADRATICS:
        line_key, n_values_in_range, n_collinear, _pts = quadratic_to_ulam_line(a, b, cc, args.N, coords)
        # Find this line in chosen_primes (if any).
        found_count = 0
        if line_key is not None:
            for k, members in chosen_primes:
                if k == line_key:
                    found_count = len(members)
                    break
        # Also: search all_lines (line_dict) — useful if the greedy
        # algorithm picked a different (overlapping) line; report
        # raw prime-count on the same line.
        raw_prime_count = len(line_dict.get(line_key, [])) if line_key else 0
        print(
            f"#   {name:<14} f(n)={a}n^2{b:+d}n{cc:+d}  vals<=N: {n_values_in_range:4d}  "
            f"collinear-on-spiral: {n_collinear:4d}/{n_values_in_range:<4d}  "
            f"primes-on-this-line: {raw_prime_count:3d}"
        )

    # Step 5: matched-baseline random
    print(f"\n# Matched-baseline (random integers in [1,N] of size pi(N)={pi_N})")
    L_baseline_list = []
    for trial in range(args.baseline_trials):
        t0 = time.time()
        random_ints = random_baseline_indices(args.N, pi_N, rng)
        random_pts = [coords[i] for i in random_ints]
        rl = all_lines_through_pairs(random_pts)
        L_rand, _ = greedy_line_cover_from_lines(rl, len(random_pts))
        L_baseline_list.append(L_rand)
        print(f"#   trial {trial+1}: L_random_greedy = {L_rand}, time {time.time()-t0:.2f}s")
    L_random_mean = sum(L_baseline_list) / len(L_baseline_list)
    L_random_var = sum((x - L_random_mean) ** 2 for x in L_baseline_list) / max(1, len(L_baseline_list) - 1)
    L_random_std = math.sqrt(L_random_var) if len(L_baseline_list) > 1 else 0.0
    print(f"#   L_random_mean = {L_random_mean:.1f} +/- {L_random_std:.1f}")

    # Z-score: prime line cover vs random
    if L_random_std > 0:
        z = (L_random_mean - L_primes) / L_random_std
        print(f"#   z-score (random - prime) / std = {z:+.2f}sigma  (positive = primes have FEWER lines)")
    else:
        z = float("nan")

    # Save summary CSV
    summary_path = os.path.join(out_dir, f"summary_N{args.N}.csv")
    with open(summary_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["N", "pi_N", "L_primes_greedy", "L_random_mean", "L_random_std",
                    "n_baseline_trials", "z_score", "n_distinct_lines"])
        w.writerow([args.N, pi_N, L_primes, f"{L_random_mean:.2f}", f"{L_random_std:.2f}",
                    args.baseline_trials, f"{z:.3f}" if not math.isnan(z) else "NA",
                    len(line_dict)])

    # Top lines CSV
    toplines_path = os.path.join(out_dir, f"top_lines_N{args.N}.csv")
    with open(toplines_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "dx", "dy", "intercept", "prime_count"])
        for row in top_lines_for_csv:
            w.writerow(row)

    print(f"\n# wrote {summary_path}, {toplines_path}")
    print(f"# DONE N={args.N}")


if __name__ == "__main__":
    main()
