"""
S235: Diagnose what's possible for the W=14 j=1 slab specifically.
The (1+3+3) BT search returned 0 candidates in rows [0, 14) but 38254
in rows [0, 28). The (1+3+3) shape is too restrictive for j=1.

This script tests an alternative shape: (1 + 5 + 1), i.e., outer
1×1 (row 0 + dead col), then inner 6×6 split as (5×5 + 1×1) with
the 1×1 being row 8 at col 0 (chiP(113) = 1, the only nonzero entry
in row 8 at the chosen 7 cols).

The (5×5 + 1×1) inner BT works iff:
- Some 5 rows from {r : M[r][col=0] = 0, r ≠ 8} × cols {2, 4, 8, 10, 12}
  give an invertible 5×5.
- The 5×5 ITSELF can then be (1+3+3) or (3+2) or (2+3) BT for further
  decomposition, OR closed via cofactor expansion (mathlib has no
  det_fin_five).

The script:
1. Identifies all 5-row subsets (from rows ∈ {1..13}\{8} with col-0 = 0)
   such that the 5×5 over cols {2, 4, 8, 10, 12} is invertible.
2. For each, attempts (3+2) or (2+3) BT decomposition of the 5×5 — these
   are single-session viable since we have det_fin_three for the 3×3 and
   det_fin_two (mathlib has Matrix.det_fin_two via Matrix.det_succ_above
   or similar) for the 2×2.
3. Reports the simplest decomposition.

USAGE:
    python3 w14_j1_analysis.py
"""
from itertools import combinations, permutations
from math import gcd
import numpy as np
from sympy import isprime, totient


def main():
    W = 14

    print(f"=== W={W} j=1 slab analysis ===")
    live = [k for k in range(W) if gcd(k + 1, W) == 1]
    print(f"Live cols (residues coprime to {W}): {live}")
    dead_candidates = [k for k in range(W)
                       if gcd(k + 1, W) != 1 and isprime(k + 1)]
    print(f"Dead-col candidates: {dead_candidates}")

    # Build the j=1 slab in numpy
    M = np.array([[1 if isprime(r * W + c + 1) else 0 for c in range(W)]
                  for r in range(W)], dtype=int)

    # Inner 13×6 (rows 1..13 × live cols)
    M_inner = M[1:W, :][:, live]
    print(f"\nInner matrix (rows 1..{W-1}, live cols):")
    for i, r in enumerate(range(1, W)):
        print(f"  r={r:2d}: {M_inner[i].tolist()}")

    # Find all 6-row subsets that give an invertible 6x6.
    invertible_6 = []
    rows_avail = list(range(W - 1))  # indices into M_inner
    for rows6 in combinations(rows_avail, 6):
        sub = M_inner[list(rows6)]
        if np.linalg.matrix_rank(sub) == 6:
            d = round(np.linalg.det(sub))
            invertible_6.append((rows6, d))

    print(f"\nInvertible 6×6 inner submatrices: {len(invertible_6)}")
    if not invertible_6:
        print("  None — even rank-6 inner is impossible. Cannot close W=14 j=1.")
        return

    # For each invertible 6×6, look for (5+1) BT decomposition where the
    # 1×1 block is some row × col with all other entries in that col zero.
    # I.e. find a (row, col) such that M_inner[rows6, col] is zero everywhere
    # EXCEPT at one row position. That row goes to the bottom; col goes
    # to the right; the upper 5×5 (other 5 rows × other 5 cols) needs to
    # be invertible.
    candidates_5plus1 = []
    for rows6, d6 in invertible_6:
        sub = M_inner[list(rows6)]  # 6×6
        for col_idx in range(6):
            col_data = sub[:, col_idx]
            nonzero_positions = [i for i, v in enumerate(col_data) if v != 0]
            if len(nonzero_positions) == 1:
                row_idx = nonzero_positions[0]
                # 1×1 block: row_idx × col_idx, value = M_inner[rows6[row_idx], col_idx]
                # 5×5 upper block: other 5 rows × other 5 cols
                other_rows = [i for i in range(6) if i != row_idx]
                other_cols = [j for j in range(6) if j != col_idx]
                upper = sub[np.ix_(other_rows, other_cols)]
                if np.linalg.matrix_rank(upper) == 5:
                    d5 = round(np.linalg.det(upper))
                    candidates_5plus1.append({
                        'rows6': rows6,
                        'd6': d6,
                        # Convert M_inner index to slab row index by +1
                        'corner_row': rows6[row_idx] + 1,
                        'corner_col': live[col_idx],
                        'upper_rows': [rows6[i] + 1 for i in other_rows],
                        'upper_cols': [live[j] for j in other_cols],
                        'd5': d5,
                    })

    print(f"\n(5+1) BT candidates (corner row has only 1 nonzero col within rows6): "
          f"{len(candidates_5plus1)}")
    if not candidates_5plus1:
        print("  No (5+1) BT structure for the inner 6×6.")
        return

    # For each (5+1) candidate, check if the 5×5 admits (3+2) or (2+3) BT.
    # Try all (row_perm, col_perm) of the 5×5 with the requested partition.
    def find_bt_5(M5, partition):
        """partition is e.g. (3, 2). Returns (row_perm, col_perm) or None."""
        a, b = partition
        for row_perm in permutations(range(5)):
            for col_perm in permutations(range(5)):
                P = M5[np.ix_(list(row_perm), list(col_perm))]
                # block lower-left (rows a..4, cols 0..a-1) zero
                if np.any(P[a:, :a] != 0):
                    continue
                # diagonal blocks invertible
                A = P[:a, :a]
                B = P[a:, a:]
                if np.linalg.matrix_rank(A) != a:
                    continue
                if np.linalg.matrix_rank(B) != b:
                    continue
                return (row_perm, col_perm)
        return None

    # We also try BD (block-diagonal) for cleanness — upper-right zero too.
    def find_bd_5(M5, partition):
        a, b = partition
        for row_perm in permutations(range(5)):
            for col_perm in permutations(range(5)):
                P = M5[np.ix_(list(row_perm), list(col_perm))]
                if np.any(P[a:, :a] != 0):
                    continue
                if np.any(P[:a, a:] != 0):
                    continue
                A = P[:a, :a]
                B = P[a:, a:]
                if np.linalg.matrix_rank(A) != a:
                    continue
                if np.linalg.matrix_rank(B) != b:
                    continue
                return (row_perm, col_perm)
        return None

    # Take the first 20 (5+1) candidates and check
    print(f"\nTesting all (5+1) candidates for inner (3+2) or (2+3) BT:")
    found = []
    for cand in candidates_5plus1:
        upper_rows_idx = [r - 1 for r in cand['upper_rows']]  # slab → M_inner
        upper_cols_idx = [live.index(c) for c in cand['upper_cols']]
        M5 = M_inner[np.ix_(upper_rows_idx, upper_cols_idx)]
        for partition in [(3, 2), (2, 3)]:
            res = find_bt_5(M5, partition)
            if res is not None:
                bd_res = find_bd_5(M5, partition)
                found.append({
                    **cand,
                    'partition': partition,
                    'row_perm': res[0],
                    'col_perm': res[1],
                    'is_bd': bd_res is not None,
                })
                break

    print(f"\n=== (5+1) → (3+2)/(2+3) BT decompositions found: {len(found)} ===")
    if not found:
        print("None. Try (1+5+1) or other shape.")
        # As a fallback, look for any (1+1+5+1+1) shape — single-row blocks
        # corresponding to rows with single-entry support.
        return

    # Score by min new helpers
    EXISTING_HELPERS = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
                        47, 53, 59, 61, 89, 97, 109, 127, 149, 179, 199, 211,
                        241, 293, 337}

    scored = []
    for entry in found:
        # Compute primes and composites in the full 7×7 submatrix.
        # rho = (0, corner_row, *upper_rows reordered by row_perm)
        # sigma = (dead_col, corner_col, *upper_cols reordered by col_perm)
        # But we're flexible about whether corner_col goes first or last.
        # For a clean (1 + 5 + 1) outer with row 0 first, corner_row last:
        #   actual_rho = [0] + reordered upper rows + [corner_row]
        #   actual_sigma = [dead] + reordered upper cols + [corner_col]
        actual_upper_rows = [entry['upper_rows'][i] for i in entry['row_perm']]
        actual_upper_cols = [entry['upper_cols'][j] for j in entry['col_perm']]
        # Test both dead choices
        for dead in dead_candidates:
            rho = [0] + actual_upper_rows + [entry['corner_row']]
            sigma = [dead] + actual_upper_cols + [entry['corner_col']]
            primes_used = set()
            for i in range(7):
                for k in range(7):
                    n = rho[i] * W + sigma[k] + 1
                    if isprime(n):
                        primes_used.add(n)
            new_primes = primes_used - EXISTING_HELPERS

            # Verify outer (1+6) is BT (dead col zero on rows 1..6 of rho).
            # That's automatically true by dead-col choice.

            # Verify (5+1) inner is BT.
            # Verify the 5×5 block is BT or BD as found.
            # For Lean assembly, we'd use nested fromBlocks.
            scored.append({
                'rho': tuple(rho),
                'sigma': tuple(sigma),
                'partition': entry['partition'],
                'is_bd': entry['is_bd'],
                'd5': entry['d5'],
                'corner_value': 1,
                'dead': dead,
                'new_primes': sorted(new_primes),
                'n_new': len(new_primes),
            })

    scored.sort(key=lambda e: (e['n_new'], 0 if e['is_bd'] else 1))
    print(f"\nTop candidates (min new helpers):")
    min_new = scored[0]['n_new']
    seen = 0
    for s in scored:
        if s['n_new'] > min_new:
            break
        print(f"\n  rho   = {s['rho']}")
        print(f"  sigma = {s['sigma']}")
        print(f"    inner partition = (5+1), 5×5 sub-partition = {s['partition']}")
        print(f"    inner 5×5 BD-able: {s['is_bd']}")
        print(f"    new prime helpers: {s['new_primes']} ({s['n_new']} total)")
        seen += 1
        if seen >= 5:
            break

    # Verify the top candidate
    best = scored[0]
    full = np.array([[1 if isprime(best['rho'][i] * W + best['sigma'][k] + 1) else 0
                      for k in range(7)] for i in range(7)])
    print(f"\n=== Verification of top candidate ===")
    print(f"Full 7×7 submatrix:")
    for row in full:
        print(f"  {row.tolist()}")
    print(f"det = {round(np.linalg.det(full))}")

    primes_used = set()
    composites = set()
    for i in range(7):
        for k in range(7):
            n = best['rho'][i] * W + best['sigma'][k] + 1
            if isprime(n):
                primes_used.add(n)
            elif n > 1:
                composites.add(n)
    composites.add(1)
    print(f"\nPrimes appearing: {sorted(primes_used)}")
    print(f"Composites to witness: {sorted(composites)}")
    new_primes = primes_used - EXISTING_HELPERS
    print(f"New chiP helpers needed: {sorted(new_primes)}")


if __name__ == "__main__":
    main()
