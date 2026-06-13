"""
S235: Find the (2, 3, 1) BT decomposition for W=14 j=1 inner 6×6 with
minimum new chiP-prime helpers needed.

Build full 7×7 BT structure: (1 outer) + (2, 3, 1) inner.
"""
from itertools import combinations
from math import gcd
import numpy as np
from sympy import isprime


EXISTING_HELPERS = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
                    47, 53, 59, 61, 89, 97, 109, 127, 149, 179, 199, 211,
                    241, 293, 337}


def admits_bt_231(M):
    """For shape (2, 3, 1), enumerate ALL valid BT decompositions."""
    n = M.shape[0]
    rows = list(range(n))
    cols = list(range(n))
    results = []
    for A_rows in combinations(rows, 2):
        rest_r1 = [r for r in rows if r not in A_rows]
        for B_rows in combinations(rest_r1, 3):
            C_rows = [r for r in rest_r1 if r not in B_rows]
            for A_cols in combinations(cols, 2):
                rest_c1 = [c for c in cols if c not in A_cols]
                for B_cols in combinations(rest_c1, 3):
                    C_cols = [c for c in rest_c1 if c not in B_cols]
                    if any(M[r][c] != 0 for r in B_rows for c in A_cols):
                        continue
                    if any(M[r][c] != 0 for r in C_rows for c in A_cols):
                        continue
                    if any(M[r][c] != 0 for r in C_rows for c in B_cols):
                        continue
                    A = M[np.ix_(list(A_rows), list(A_cols))]
                    if np.linalg.matrix_rank(A) != 2:
                        continue
                    B = M[np.ix_(list(B_rows), list(B_cols))]
                    if np.linalg.matrix_rank(B) != 3:
                        continue
                    C = M[np.ix_(list(C_rows), list(C_cols))]
                    if np.linalg.matrix_rank(C) != 1:
                        continue
                    results.append((A_rows, list(B_rows), C_rows,
                                    A_cols, list(B_cols), C_cols))
    return results


def main():
    W = 14
    live = [k for k in range(W) if gcd(k + 1, W) == 1]
    M = np.array([[1 if isprime(r * W + c + 1) else 0 for c in range(W)]
                  for r in range(W)], dtype=int)
    M_inner = M[1:W, :][:, live]

    invertible_6 = []
    rows_avail = list(range(W - 1))
    for rows6 in combinations(rows_avail, 6):
        sub = M_inner[list(rows6)]
        if np.linalg.matrix_rank(sub) == 6:
            invertible_6.append(rows6)

    print(f"Invertible inner 6×6: {len(invertible_6)}")

    # Enumerate all (2, 3, 1) BT decompositions across all invertible row sets.
    all_candidates = []
    for rows6 in invertible_6:
        sub = M_inner[list(rows6)]
        decomps = admits_bt_231(sub)
        for d in decomps:
            all_candidates.append((rows6, d))

    print(f"Total (2, 3, 1) BT decompositions: {len(all_candidates)}")

    # For each candidate and each dead-col choice, score by new helpers needed.
    dead_candidates = [k for k in range(W)
                       if gcd(k + 1, W) != 1 and isprime(k + 1)]

    scored = []
    for rows6, (Ar, Br, Cr, Ac, Bc, Cc) in all_candidates:
        # Convert M_inner indices to slab indices and live-col positions to col ks.
        actual_A_rows = [rows6[i] + 1 for i in Ar]
        actual_B_rows = [rows6[i] + 1 for i in Br]
        actual_C_rows = [rows6[i] + 1 for i in Cr]
        actual_A_cols = [live[j] for j in Ac]
        actual_B_cols = [live[j] for j in Bc]
        actual_C_cols = [live[j] for j in Cc]
        all_rows = [0] + actual_A_rows + actual_B_rows + actual_C_rows
        for dead in dead_candidates:
            all_cols = [dead] + actual_A_cols + actual_B_cols + actual_C_cols
            # Compute primes in 7×7 sub
            primes_used = set()
            for r in all_rows:
                for c in all_cols:
                    n = r * W + c + 1
                    if isprime(n):
                        primes_used.add(n)
            new_primes = primes_used - EXISTING_HELPERS
            # Compute composites that need decide-witnesses
            composites = set()
            for r in all_rows:
                for c in all_cols:
                    n = r * W + c + 1
                    if not isprime(n) and n > 1:
                        composites.add(n)
            composites.add(1)
            scored.append({
                'all_rows': all_rows,
                'all_cols': all_cols,
                'A_rows': actual_A_rows,
                'B_rows': actual_B_rows,
                'C_rows': actual_C_rows,
                'A_cols': actual_A_cols,
                'B_cols': actual_B_cols,
                'C_cols': actual_C_cols,
                'dead': dead,
                'new_primes': sorted(new_primes),
                'n_new': len(new_primes),
                'n_composites': len(composites),
                'max_row': max(all_rows),
            })

    # Sort by (n_new, max_row, n_composites)
    scored.sort(key=lambda x: (x['n_new'], x['max_row'], x['n_composites']))

    print(f"\nMin new helpers across all candidates: {scored[0]['n_new']}")
    print(f"Min max_row at min new helpers: "
          f"{min(s['max_row'] for s in scored if s['n_new'] == scored[0]['n_new'])}")
    print()

    # Show the top 5 candidates
    print("=== Top 5 candidates (min new helpers, then min max_row) ===")
    for i, s in enumerate(scored[:5]):
        print(f"\n--- Candidate #{i+1} ---")
        print(f"  rho   = (0, {s['A_rows'][0]}, {s['A_rows'][1]}, "
              f"{s['B_rows'][0]}, {s['B_rows'][1]}, {s['B_rows'][2]}, "
              f"{s['C_rows'][0]})")
        print(f"  sigma = ({s['dead']}, {s['A_cols'][0]}, {s['A_cols'][1]}, "
              f"{s['B_cols'][0]}, {s['B_cols'][1]}, {s['B_cols'][2]}, "
              f"{s['C_cols'][0]})")
        print(f"  max_row = {s['max_row']}")
        print(f"  new prime helpers: {s['new_primes']} ({s['n_new']} total)")
        print(f"  composites to witness: {s['n_composites']}")
        # Verify
        rho = [0] + s['A_rows'] + s['B_rows'] + s['C_rows']
        sigma = [s['dead']] + s['A_cols'] + s['B_cols'] + s['C_cols']
        full = np.array([[1 if isprime(rho[i] * W + sigma[k] + 1) else 0
                          for k in range(7)] for i in range(7)])
        det = round(np.linalg.det(full))
        print(f"  full 7×7 det: {det}")
        # Print the matrix
        for row in full:
            print(f"    {row.tolist()}")


if __name__ == "__main__":
    main()
