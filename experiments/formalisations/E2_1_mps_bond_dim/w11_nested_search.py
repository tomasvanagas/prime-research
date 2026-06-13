"""
S206: Search for a nested block-triangular structure for W=11.

Outer split: 1 + 10 (BlockTriangular with row 0 + dead col at top-left
1×1 with chiP(11) = 1, and 10×10 inner block over rows {1..10} × live cols
{0..9}).

Inner 10×10 split: try block-triangular structures with blocks of size
≤ 3 (since mathlib has det_fin_one / det_fin_two / det_fin_three but not
det_fin_four). The cleanest target is 1 + 3 + 3 + 3.

In a block-TRIANGULAR (1+3+3+3) inner structure (block ordering
B_1 < B_2 < B_3 < B_4):

  - B_1 × any other block: free (upper-triangular, no constraint on B_1 row)
  - B_2 × B_1: zero  (3 rows × 1 col)
  - B_3 × (B_1 ∪ B_2): zero  (3 rows × 4 cols)
  - B_4 × (B_1 ∪ B_2 ∪ B_3): zero  (3 rows × 7 cols)
  - Each diagonal block must be invertible.

The outer 1+10 split's bottom-left zero requires: rows 1..10 have 0 at the
dead col (k=10). chiP(r*11 + 10 + 1) = chiP(11r + 11) = chiP(11(r+1)) = 0
for r ≥ 1 (since 11(r+1) is composite for r ≥ 1). Verified.

USAGE:
    python3 w11_nested_search.py
"""
from itertools import combinations, permutations
from math import gcd
from sympy import isprime, totient, Matrix


def find_nested_1_3_3_3(W=11, row_pool_max=11, max_to_keep=200):
    """Search for outer 1+10, inner 1+3+3+3 BlockTriangular over rows
    [1, row_pool_max) × live cols of W=11."""
    live = sorted(k for k in range(W) if gcd(k + 1, W) == 1)
    R_inner = 10
    assert len(live) == R_inner
    cols_set = list(live)  # 10 live cols
    rows_avail = list(range(1, row_pool_max))

    # Build the inner matrix
    M = [[1 if isprime(r * W + c + 1) else 0 for c in cols_set]
         for r in range(row_pool_max)]

    candidates = []
    # Sweep: pick (r1, c1) for the inner 1×1 block; then partition the
    # remaining 9 cols into A_cols, B_cols, C_cols and the 9 remaining
    # rows from rows_avail \ {r1} into A_rows, B_rows, C_rows such that
    #   - col c1 zero on (A_rows ∪ B_rows ∪ C_rows)
    #   - cols A zero on (B_rows ∪ C_rows)
    #   - cols B zero on (C_rows)
    #   - each block (A, B, C) invertible (det_fin_three)
    #   - row r1's entry at c1 is 1 (for det_fin_one block = [[1]])
    for r1 in rows_avail:
        # row r1 must have a "dominant" 1 at some col c1 such that no
        # row in rows_avail \ {r1} placed in B or C blocks has nonzero at c1.
        # Try every c1 where M[r1][c1] = 1.
        for c1_idx, c1 in enumerate(cols_set):
            if M[r1][c1_idx] != 1:
                continue
            # Define candidate ABC rows: rows_avail \ {r1}
            ABC_rows_pool = [r for r in rows_avail if r != r1]
            # Filter: rows that have 0 at col c1 (otherwise can't be in B/C)
            ABC_rows_filtered = [r for r in ABC_rows_pool
                                 if M[r][c1_idx] == 0]
            if len(ABC_rows_filtered) < 9:
                continue
            ABC_cols_idx = [i for i in range(R_inner) if i != c1_idx]
            assert len(ABC_cols_idx) == 9

            # For each (3, 3, 3) ordered partition of cols
            for A_cols_idx in combinations(ABC_cols_idx, 3):
                rest_after_A = [i for i in ABC_cols_idx if i not in A_cols_idx]
                for B_cols_idx in combinations(rest_after_A, 3):
                    C_cols_idx = tuple(i for i in rest_after_A
                                       if i not in B_cols_idx)

                    # B-block rows must be 0 on A_cols
                    B_eligible = [r for r in ABC_rows_filtered
                                  if all(M[r][i] == 0 for i in A_cols_idx)]
                    # C-block rows must be 0 on A_cols ∪ B_cols
                    C_eligible = [r for r in ABC_rows_filtered
                                  if all(M[r][i] == 0
                                         for i in A_cols_idx + B_cols_idx)]
                    if len(C_eligible) < 3 or len(B_eligible) < 3:
                        continue

                    # First, find C_rows ⊆ C_eligible with detC ≠ 0
                    for C_rows in combinations(C_eligible, 3):
                        MC = Matrix([[M[r][i] for i in C_cols_idx]
                                     for r in C_rows])
                        dC = MC.det()
                        if dC == 0:
                            continue
                        # B_rows ⊆ B_eligible \ C_rows with detB ≠ 0
                        B_pool = [r for r in B_eligible if r not in C_rows]
                        for B_rows in combinations(B_pool, 3):
                            MB = Matrix([[M[r][i] for i in B_cols_idx]
                                         for r in B_rows])
                            dB = MB.det()
                            if dB == 0:
                                continue
                            # A_rows ⊆ ABC_rows_filtered \ (C_rows ∪ B_rows)
                            A_pool = [r for r in ABC_rows_filtered
                                      if r not in C_rows and r not in B_rows]
                            for A_rows in combinations(A_pool, 3):
                                MA = Matrix([[M[r][i] for i in A_cols_idx]
                                             for r in A_rows])
                                dA = MA.det()
                                if dA == 0:
                                    continue
                                candidates.append({
                                    'r1': r1, 'c1': c1, 'c1_idx': c1_idx,
                                    'A_rows': A_rows, 'A_cols_idx': A_cols_idx,
                                    'A_kvals': tuple(cols_set[i] for i in A_cols_idx),
                                    'B_rows': B_rows, 'B_cols_idx': B_cols_idx,
                                    'B_kvals': tuple(cols_set[i] for i in B_cols_idx),
                                    'C_rows': C_rows, 'C_cols_idx': C_cols_idx,
                                    'C_kvals': tuple(cols_set[i] for i in C_cols_idx),
                                    'detA': dA, 'detB': dB, 'detC': dC,
                                })
                                if len(candidates) >= max_to_keep:
                                    return candidates, cols_set
    return candidates, cols_set


def primes_used(cand, W=11):
    """All primes appearing as 1-entries in the 11×11 submatrix.

    Recall: outer 1×1 is at row 0, dead col k=10 → chiP(11) = 1.
    Inner 1×1 is at row r1, col c1 → chiP(r1*11 + c1 + 1) = 1.
    Inner 3×3 blocks A, B, C contribute their entries.
    Plus row 0 may have 1s on inner cols (allowed, upper-right of outer).
    """
    primes = set()
    # Outer 1×1: chiP(11) = 1 (W = 11 itself is prime)
    primes.add(11)
    # Inner 1×1: chiP at (r1, c1)
    primes.add(cand['r1'] * W + cand['c1'] + 1)
    # ABC blocks
    for r in cand['A_rows']:
        for k in cand['A_kvals']:
            n = r * W + k + 1
            if isprime(n):
                primes.add(n)
    for r in cand['B_rows']:
        for k in cand['B_kvals']:
            n = r * W + k + 1
            if isprime(n):
                primes.add(n)
    for r in cand['C_rows']:
        for k in cand['C_kvals']:
            n = r * W + k + 1
            if isprime(n):
                primes.add(n)
    # Upper-right of outer (row 0 × inner cols): allowed but contribute primes
    # for the chiP entry computations. These are entries not affecting det_X
    # but DO appear in the submatrix definition for hE0i etc.
    for k in cols_set_global:
        if k == cand['c1']:
            continue
        n = 0 * W + k + 1
        if isprime(n):
            primes.add(n)
    # Upper-right entries from row 0 to inner cols
    return primes


cols_set_global = None  # will be set below


if __name__ == "__main__":
    EXISTING_HELPERS = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
                        47, 53, 59, 61, 89, 97, 109, 127, 149, 179, 199,
                        211, 241, 293, 337}

    # First try rows in [1, 11)
    print("=== W=11 nested 1+(1+3+3+3) BlockTriangular search ===")
    print()
    cands, cols_set = find_nested_1_3_3_3(W=11, row_pool_max=11)
    cols_set_global = cols_set
    print(f"rows in [1, 11): found {len(cands)} candidates.")
    if not cands:
        print("\nNo candidates in [1, 11); retrying with [1, 22).")
        cands, cols_set = find_nested_1_3_3_3(W=11, row_pool_max=22)
        cols_set_global = cols_set
        print(f"rows in [1, 22): found {len(cands)} candidates.")

    if not cands:
        print("\nNo candidates found.")
        exit(0)

    # Score
    scored = []
    for c in cands:
        primes = primes_used(c)
        new_primes = primes - EXISTING_HELPERS
        # Prefer rows in [1, 11) (j=1 case viable) -- compute max row
        max_row = max(c['r1'], max(c['A_rows']), max(c['B_rows']), max(c['C_rows']))
        scored.append((max_row, len(new_primes), sorted(new_primes), c))
    # Sort by max_row first (j=1 viable), then by new helpers
    scored.sort(key=lambda x: (x[0], x[1]))

    # Take the lowest max_row candidates first
    min_max_row = scored[0][0]
    print(f"\nBest max_row achievable: {min_max_row} (j ≥ ⌈log_11({min_max_row}+1)⌉ needed)")
    if min_max_row < 11:
        print("  → j ≥ 1 sufficient (W^1 = 11 > all rows)")
    elif min_max_row < 121:
        print("  → j ≥ 2 needed (W^2 = 121 > all rows)")

    # Filter to top 5 with min_max_row
    top = [s for s in scored if s[0] == min_max_row][:5]
    if len(top) < 5:
        # Add more
        others = [s for s in scored if s[0] > min_max_row][:5 - len(top)]
        top.extend(others)

    print()
    print("Top candidates (sorted by max_row, then new_primes):")
    for max_row, n_new, new_primes, c in top[:8]:
        rho_full = (0, c['r1']) + tuple(c['A_rows']) + tuple(c['B_rows']) + tuple(c['C_rows'])
        sigma_full = (10, c['c1']) + c['A_kvals'] + c['B_kvals'] + c['C_kvals']
        print(f"  max_row = {max_row}")
        print(f"  rho   = {rho_full}")
        print(f"  sigma = {sigma_full}")
        print(f"    detA = {c['detA']}, detB = {c['detB']}, detC = {c['detC']}")
        diag_primes_inner = (c['r1']*11 + c['c1'] + 1,)
        for i in range(3):
            diag_primes_inner += (c['A_rows'][i]*11 + c['A_kvals'][i] + 1,)
        for i in range(3):
            diag_primes_inner += (c['B_rows'][i]*11 + c['B_kvals'][i] + 1,)
        for i in range(3):
            diag_primes_inner += (c['C_rows'][i]*11 + c['C_kvals'][i] + 1,)
        print(f"    'diagonal' primes inner blocks: (varies, blocks are 3×3)")
        print(f"    new primes: {new_primes} (count={n_new})")
        print()

    # Verify the top candidate is actually block-triangular
    cand = scored[0][3]
    rho_full = [0, cand['r1']] + list(cand['A_rows']) + list(cand['B_rows']) + list(cand['C_rows'])
    sigma_full = [10, cand['c1']] + list(cand['A_kvals']) + list(cand['B_kvals']) + list(cand['C_kvals'])
    M_full = Matrix([[1 if isprime(rho_full[i] * 11 + sigma_full[k] + 1) else 0
                      for k in range(11)] for i in range(11)])
    print("=== Verification of top candidate ===")
    print("Full 11×11 submatrix M[i][k] = chiP(rho[i]*11 + sigma[k] + 1):")
    for i, row in enumerate(M_full.tolist()):
        print(f"  row {i} (rho={rho_full[i]}): {row}")
    print(f"det(full) = {M_full.det()}")
    print(f"Block determinants: 1 * 1 * detA * detB * detC = "
          f"{cand['detA'] * cand['detB'] * cand['detC']}")
