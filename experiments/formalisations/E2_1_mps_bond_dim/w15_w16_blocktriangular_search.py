"""
S245: BlockTriangular pre-search for W=15 and W=16 orthogonal-corner cases
of E2.1 (MPS bond dimension).

Both W=15 = 3·5 and W=16 = 2^4 have R = phi(W) + 1 = 9. They are on the
S128/S143/S144 "block-triangular-required" set (leading-row + dead-col
upper-triangulation route is exhausted).

The W=14 closure (S235, R=7) used a 1+(2+3+1) three-level fromBlocks
nest. For R=9 the analogous routes are partitions of 8 with parts ≤ 4,
since mathlib provides det_fin_one/two/three/four (det_fin_four via
Matrix.det_succ_row + det_fin_three) but NO det_fin_five at the time of
writing.

This implementation uses BITMASK-based zero-block checks for speed:
the W-row support pattern over the chosen 8 cols becomes a 8-bit mask,
and "lower-left zero" reduces to a bitwise AND test. Only candidates
that pass the bitwise filter pay for a numeric rank/det check.

USAGE:
    python3 w15_w16_blocktriangular_search.py
"""
from itertools import combinations
from math import gcd
import numpy as np
from sympy import isprime, totient


def ordered_partitions(n, max_part):
    """All ordered partitions of n into parts of size in [1, max_part]."""
    if n == 0:
        yield ()
        return
    for k in range(1, min(n, max_part) + 1):
        for rest in ordered_partitions(n - k, max_part):
            yield (k,) + rest


def find_bt_decomp_bitmask(M, sigs, sizes, want_block_diagonal=False):
    """
    Find BlockTriangular decomposition of size pattern `sizes` using
    bitmask zero-block checks for speed.

    M[i, k] in {0, 1}, sigs[i] = bitmask of M[i, :] over the 8 cols.
    Returns the FIRST hit; not exhaustive.
    """
    n = M.shape[0]
    assert n == 8 and sum(sizes) == 8
    rows = list(range(n))

    def recurse(avail_rows_mask, avail_cols_mask, idx):
        if idx == len(sizes):
            yield [], []
            return
        s = sizes[idx]
        avail_r = [r for r in rows if avail_rows_mask & (1 << r)]
        avail_c = [c for c in rows if avail_cols_mask & (1 << c)]
        for A_rows in combinations(avail_r, s):
            A_rows_mask = 0
            for r in A_rows:
                A_rows_mask |= 1 << r
            rest_rows_mask = avail_rows_mask & ~A_rows_mask
            rest_rows = [r for r in avail_r if not A_rows_mask & (1 << r)]
            for A_cols in combinations(avail_c, s):
                A_cols_mask = 0
                for c in A_cols:
                    A_cols_mask |= 1 << c
                # Lower-left zero: rest_rows × A_cols all zero
                bad = False
                for r in rest_rows:
                    if sigs[r] & A_cols_mask:
                        bad = True
                        break
                if bad:
                    continue
                # Optional: block-diagonal — also require A_rows × rest_cols zero
                rest_cols_mask = avail_cols_mask & ~A_cols_mask
                if want_block_diagonal:
                    bad_bd = False
                    for r in A_rows:
                        if sigs[r] & rest_cols_mask:
                            bad_bd = True
                            break
                    if bad_bd:
                        continue
                A = M[np.ix_(list(A_rows), list(A_cols))]
                # Use det check (faster + already needed)
                detA = int(round(np.linalg.det(A)))
                if detA == 0:
                    continue
                for sub_r, sub_c in recurse(rest_rows_mask, rest_cols_mask,
                                            idx + 1):
                    yield ([list(A_rows)] + sub_r,
                           [list(A_cols)] + sub_c)

    full_rows_mask = (1 << n) - 1
    full_cols_mask = (1 << n) - 1
    for block_rows, block_cols in recurse(full_rows_mask, full_cols_mask, 0):
        block_dets = []
        for s, A_rows, A_cols in zip(sizes, block_rows, block_cols):
            A = M[np.ix_(A_rows, A_cols)]
            block_dets.append(int(round(np.linalg.det(A))))
        full_det = 1
        for d in block_dets:
            full_det *= d
        return (block_rows, block_cols, block_dets, full_det)
    return None


def primes_in_submatrix(W, rho, sigma):
    primes = set()
    R = len(rho)
    for i in range(R):
        for k in range(R):
            n = rho[i] * W + sigma[k] + 1
            if isprime(n):
                primes.add(n)
    return primes


def composites_in_submatrix(W, rho, sigma):
    composites = set()
    R = len(rho)
    for i in range(R):
        for k in range(R):
            n = rho[i] * W + sigma[k] + 1
            if not isprime(n) and n >= 1:
                composites.add(n)
    return composites


def search_W(W, EXISTING_HELPERS):
    print(f"=== W={W} BlockTriangular pre-search ===", flush=True)
    R = int(totient(W)) + 1
    print(f"R = phi({W}) + 1 = {R}", flush=True)
    live = [k for k in range(W) if gcd(k + 1, W) == 1]
    print(f"Live cols (k+1 coprime to {W}): {live}  (count = {len(live)})",
          flush=True)
    dead_candidates = [k for k in range(W)
                       if gcd(k + 1, W) != 1 and isprime(k + 1)]
    print(f"Dead-col candidates (k+1 prime, gcd != 1): {dead_candidates}",
          flush=True)
    print(flush=True)

    inner_size = R - 1  # 8 for W=15, W=16

    # Iterate shapes from "fewest blocks" to "most blocks" — Lean ergonomics
    # prefer shallow nests. Sort by (block count, max part desc).
    shapes = sorted(ordered_partitions(inner_size, 4),
                    key=lambda s: (len(s), -max(s)))
    print(f"Inner-block shapes (ordered partitions of {inner_size}, "
          f"parts ≤ 4): {len(shapes)} shapes", flush=True)
    print(flush=True)

    all_results = []

    for dead in dead_candidates:
        print(f"--- Dead col k={dead} (chiP({dead + 1}) = 1) ---", flush=True)
        # Dead col must be zero on rows ≥ 1 for the outer (1 + inner) BT.
        bad_rows = [r for r in range(1, W) if isprime(r * W + dead + 1)]
        if bad_rows:
            print(f"  [skip] dead col k={dead} nonzero on rows {bad_rows}",
                  flush=True)
            continue

        # Build M_inner_full: (W-1) × 8 over rows [1, W) and live cols.
        M_inner_full = np.array(
            [[1 if isprime(r * W + c + 1) else 0 for c in live]
             for r in range(1, W)], dtype=int)
        # Bitmask sigs.
        sigs_full = []
        for r in range(W - 1):
            s = 0
            for k in range(inner_size):
                if M_inner_full[r, k]:
                    s |= 1 << k
            sigs_full.append(s)

        rows_avail = list(range(W - 1))
        invertible_inner = []
        for rows in combinations(rows_avail, inner_size):
            sub = M_inner_full[list(rows)]
            if np.linalg.matrix_rank(sub) == inner_size:
                invertible_inner.append(rows)
        print(f"  Invertible inner {inner_size}x{inner_size} subsets: "
              f"{len(invertible_inner)}", flush=True)

        for shape in shapes:
            best_for_shape = None
            for rows in invertible_inner:
                sub = M_inner_full[list(rows)]
                sub_sigs = [sigs_full[r] for r in rows]
                res = find_bt_decomp_bitmask(sub, sub_sigs, shape)
                if res is None:
                    continue
                block_rows, block_cols, block_dets, full_det = res
                if full_det == 0:
                    continue

                rho = [0]
                sigma = [dead]
                for AR, AC in zip(block_rows, block_cols):
                    for r_idx in AR:
                        rho.append(rows[r_idx] + 1)
                    for c_idx in AC:
                        sigma.append(live[c_idx])

                primes = primes_in_submatrix(W, rho, sigma)
                new_primes = primes - EXISTING_HELPERS

                cand = {
                    'shape': shape,
                    'rows_set': rows,
                    'block_rows': block_rows,
                    'block_cols': block_cols,
                    'block_dets': block_dets,
                    'full_det': full_det,
                    'rho': rho,
                    'sigma': sigma,
                    'new_primes': sorted(new_primes),
                    'primes': sorted(primes),
                    'max_row': max(rho),
                    'dead': dead,
                }
                if (best_for_shape is None
                        or len(new_primes) < len(best_for_shape['new_primes'])
                        or (len(new_primes) == len(best_for_shape['new_primes'])
                            and max(rho) < best_for_shape['max_row'])):
                    best_for_shape = cand
                # Early exit if we already have a 0-new-helpers candidate.
                if len(new_primes) == 0:
                    break

            if best_for_shape is None:
                print(f"  shape {shape}: NO HIT", flush=True)
            else:
                print(f"  shape {shape}: HIT, det={best_for_shape['full_det']:+d}, "
                      f"max_row={best_for_shape['max_row']}, "
                      f"new helpers={len(best_for_shape['new_primes'])}: "
                      f"{best_for_shape['new_primes']}", flush=True)
                all_results.append(best_for_shape)
        print(flush=True)

    if not all_results:
        print(f"=== W={W}: NO BlockTriangular decomposition found ===",
              flush=True)
        return None

    # Pick best by (n_new_primes, num_blocks, max_row).
    all_results.sort(key=lambda x: (len(x['new_primes']), len(x['shape']),
                                    x['max_row']))
    best = all_results[0]
    print(f"=== Best candidate for W={W} ===", flush=True)
    print(f"  shape       = {best['shape']}", flush=True)
    print(f"  dead col    = {best['dead']}", flush=True)
    print(f"  rho         = {best['rho']}", flush=True)
    print(f"  sigma       = {best['sigma']}", flush=True)
    print(f"  block_dets  = {best['block_dets']}", flush=True)
    print(f"  full det    = {best['full_det']}", flush=True)
    print(f"  max_row     = {best['max_row']} (need {W}^j ≥ "
          f"{best['max_row'] + 1}; j ≥ 1 OK if max_row < {W})", flush=True)
    print(f"  primes      = {best['primes']}", flush=True)
    print(f"  new helpers = {best['new_primes']}", flush=True)

    rho = best['rho']
    sigma = best['sigma']
    full = np.array([[1 if isprime(rho[i] * W + sigma[k] + 1) else 0
                      for k in range(R)] for i in range(R)], dtype=int)
    print(f"\n  Full {R}x{R} submatrix (rows = rho, cols = sigma):",
          flush=True)
    for row in full.tolist():
        print(f"    {row}", flush=True)
    det_val = int(round(np.linalg.det(full)))
    print(f"  Verified det = {det_val}", flush=True)

    composites = composites_in_submatrix(W, rho, sigma)
    composites.add(1)
    print(f"  Composites to witness: {sorted(composites)}", flush=True)

    return best


if __name__ == "__main__":
    EXISTING_HELPERS = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
        67, 89, 97, 109, 113, 127, 149, 173, 179, 181, 199, 211, 241, 293,
        337
    }
    print(f"Existing chiP_X_eq_one helpers: {sorted(EXISTING_HELPERS)}",
          flush=True)
    print(flush=True)

    for W in (15, 16):
        search_W(W, EXISTING_HELPERS)
        print(flush=True)
        print("=" * 60, flush=True)
        print(flush=True)
