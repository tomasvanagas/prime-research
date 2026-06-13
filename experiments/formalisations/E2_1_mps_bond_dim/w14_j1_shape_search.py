"""
S235: Efficient block-shape search for W=14 j=1 inner 6×6.

Tests specific shapes with bitmask-based partition enumeration.
"""
from itertools import combinations
from math import gcd
import numpy as np
from sympy import isprime


def admits_bt_decomp_2(M, sizes):
    """For shape (a, b), test if BT exists (lower-left zero, both dets)."""
    n = M.shape[0]
    a, b = sizes
    rows = list(range(n))
    cols = list(range(n))
    for A_rows in combinations(rows, a):
        B_rows = [r for r in rows if r not in A_rows]
        for A_cols in combinations(cols, a):
            B_cols = [c for c in cols if c not in A_cols]
            # Lower-left: B_rows × A_cols zero
            if any(M[r][c] != 0 for r in B_rows for c in A_cols):
                continue
            A = M[np.ix_(list(A_rows), list(A_cols))]
            if np.linalg.matrix_rank(A) != a:
                continue
            B = M[np.ix_(list(B_rows), list(B_cols))]
            if np.linalg.matrix_rank(B) != b:
                continue
            return (A_rows, list(B_rows), A_cols, list(B_cols))
    return None


def admits_bt_decomp_3(M, sizes):
    """For shape (a, b, c), test BT existence."""
    n = M.shape[0]
    a, b, c = sizes
    rows = list(range(n))
    cols = list(range(n))
    for A_rows in combinations(rows, a):
        rest_r1 = [r for r in rows if r not in A_rows]
        for B_rows in combinations(rest_r1, b):
            C_rows = [r for r in rest_r1 if r not in B_rows]
            for A_cols in combinations(cols, a):
                rest_c1 = [c for c in cols if c not in A_cols]
                for B_cols in combinations(rest_c1, b):
                    C_cols = [c for c in rest_c1 if c not in B_cols]
                    # Check zero blocks
                    if any(M[r][c] != 0 for r in B_rows for c in A_cols):
                        continue
                    if any(M[r][c] != 0 for r in C_rows for c in A_cols):
                        continue
                    if any(M[r][c] != 0 for r in C_rows for c in B_cols):
                        continue
                    A = M[np.ix_(list(A_rows), list(A_cols))]
                    if np.linalg.matrix_rank(A) != a:
                        continue
                    B = M[np.ix_(list(B_rows), list(B_cols))]
                    if np.linalg.matrix_rank(B) != b:
                        continue
                    C = M[np.ix_(list(C_rows), list(C_cols))]
                    if np.linalg.matrix_rank(C) != c:
                        continue
                    return (A_rows, list(B_rows), C_rows,
                            A_cols, list(B_cols), C_cols)
    return None


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
    print()

    # Test each 2-block shape
    shapes_2 = [(5, 1), (1, 5), (4, 2), (2, 4), (3, 3)]
    for shape in shapes_2:
        hits = 0
        first_hit = None
        for rows6 in invertible_6:
            sub = M_inner[list(rows6)]
            res = admits_bt_decomp_2(sub, shape)
            if res is not None:
                hits += 1
                if first_hit is None:
                    first_hit = (rows6, res)
                if hits >= 3:
                    break
        if first_hit:
            rows6, (Ar, Br, Ac, Bc) = first_hit
            print(f"  Shape {shape}: HITS≥{hits}, e.g.")
            print(f"    A_rows = {[rows6[i]+1 for i in Ar]}, B_rows = {[rows6[i]+1 for i in Br]}")
            print(f"    A_cols = {[live[j] for j in Ac]}, B_cols = {[live[j] for j in Bc]}")
        else:
            print(f"  Shape {shape}: 0 hits")

    # Test 3-block shapes with parts ≤ 4
    shapes_3 = [(4, 1, 1), (1, 4, 1), (1, 1, 4),
                (3, 2, 1), (3, 1, 2), (2, 3, 1), (1, 3, 2), (2, 1, 3), (1, 2, 3),
                (2, 2, 2)]
    print()
    for shape in shapes_3:
        hits = 0
        first_hit = None
        for rows6 in invertible_6:
            sub = M_inner[list(rows6)]
            res = admits_bt_decomp_3(sub, shape)
            if res is not None:
                hits += 1
                if first_hit is None:
                    first_hit = (rows6, res)
                if hits >= 3:
                    break
        if first_hit:
            rows6, (Ar, Br, Cr, Ac, Bc, Cc) = first_hit
            print(f"  Shape {shape}: HITS≥{hits}, e.g.")
            print(f"    A_rows = {[rows6[i]+1 for i in Ar]}, "
                  f"B_rows = {[rows6[i]+1 for i in Br]}, "
                  f"C_rows = {[rows6[i]+1 for i in Cr]}")
            print(f"    A_cols = {[live[j] for j in Ac]}, "
                  f"B_cols = {[live[j] for j in Bc]}, "
                  f"C_cols = {[live[j] for j in Cc]}")
        else:
            print(f"  Shape {shape}: 0 hits")


if __name__ == "__main__":
    main()
