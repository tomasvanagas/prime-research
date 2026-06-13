"""
S245: Verify the 4×4 inner blocks of the W=15 and W=16 best candidates are
"part-3 atomic" — they admit no further BT decomposition with parts ≤ 3.

This isolates the precise structural reason BOTH W=15 and W=16 require
det_fin_four (and not just det_fin_three): the 4×4 sub-blocks themselves
do not split.

USAGE:
    python3 w15_w16_inner_4x4_atomicity.py
"""
from itertools import combinations
import numpy as np
from sympy import isprime


def part3_atomic(M):
    """Return True iff M admits no BT decomposition with parts ≤ 3.

    For 4×4: shapes (1, 3), (3, 1), (1, 1, 2), (1, 2, 1), (2, 1, 1),
    (2, 2), (1, 1, 1, 1).
    """
    n = M.shape[0]
    rows = list(range(n))

    def admits_bt(sizes):
        """Return True iff M admits BT with the given sizes."""
        def rec(prev_r, prev_c, idx):
            if idx == len(sizes):
                return True
            s = sizes[idx]
            avail_r = [r for r in rows if r not in prev_r]
            avail_c = [c for c in rows if c not in prev_c]
            for AR in combinations(avail_r, s):
                rest_r = [r for r in avail_r if r not in AR]
                for AC in combinations(avail_c, s):
                    if any(M[r, c] != 0 for r in rest_r for c in AC):
                        continue
                    A = M[np.ix_(list(AR), list(AC))]
                    if np.linalg.matrix_rank(A) != s:
                        continue
                    if rec(prev_r + list(AR), prev_c + list(AC), idx + 1):
                        return True
            return False
        return rec([], [], 0)

    shapes = [(1, 3), (3, 1), (2, 2), (1, 1, 2), (1, 2, 1), (2, 1, 1),
              (1, 1, 1, 1)]
    for shape in shapes:
        if admits_bt(shape):
            return (False, shape)
    return (True, None)


def main():
    # W=15 best (S245): A-block on rows {1, 3, 7, 13} × cols {1, 3, 7, 13}
    W = 15
    A_rows = [1, 3, 7, 13]
    A_cols = [1, 3, 7, 13]
    A = np.array([[1 if isprime(r * W + c + 1) else 0 for c in A_cols]
                  for r in A_rows], dtype=int)
    print(f"W={W} A-block (rows {A_rows} × cols {A_cols}):")
    for row in A.tolist():
        print(f"  {row}")
    print(f"  det = {int(round(np.linalg.det(A)))}")
    atomic, shape = part3_atomic(A)
    print(f"  parts-≤3 atomic: {atomic}"
          + (f" (admits {shape})" if not atomic else ""))
    print()

    # W=15 best (S245): B-block on rows {2, 6, 8, 12} × cols {0, 6, 10, 12}
    B_rows = [2, 6, 8, 12]
    B_cols = [0, 6, 10, 12]
    B = np.array([[1 if isprime(r * W + c + 1) else 0 for c in B_cols]
                  for r in B_rows], dtype=int)
    print(f"W={W} B-block (rows {B_rows} × cols {B_cols}):")
    for row in B.tolist():
        print(f"  {row}")
    print(f"  det = {int(round(np.linalg.det(B)))}")
    atomic, shape = part3_atomic(B)
    print(f"  parts-≤3 atomic: {atomic}"
          + (f" (admits {shape})" if not atomic else ""))
    print()

    # W=16 best (S245): A on rows {1, 2, 3, 7} × cols {0, 6, 10, 12}
    W = 16
    A_rows = [1, 2, 3, 7]
    A_cols = [0, 6, 10, 12]
    A = np.array([[1 if isprime(r * W + c + 1) else 0 for c in A_cols]
                  for r in A_rows], dtype=int)
    print(f"W={W} A-block (rows {A_rows} × cols {A_cols}):")
    for row in A.tolist():
        print(f"  {row}")
    print(f"  det = {int(round(np.linalg.det(A)))}")
    atomic, shape = part3_atomic(A)
    print(f"  parts-≤3 atomic: {atomic}"
          + (f" (admits {shape})" if not atomic else ""))
    print()

    # W=16 best (S245): B on rows {5, 11, 13, 14} × cols {2, 4, 8, 14}
    B_rows = [5, 11, 13, 14]
    B_cols = [2, 4, 8, 14]
    B = np.array([[1 if isprime(r * W + c + 1) else 0 for c in B_cols]
                  for r in B_rows], dtype=int)
    print(f"W={W} B-block (rows {B_rows} × cols {B_cols}):")
    for row in B.tolist():
        print(f"  {row}")
    print(f"  det = {int(round(np.linalg.det(B)))}")
    atomic, shape = part3_atomic(B)
    print(f"  parts-≤3 atomic: {atomic}"
          + (f" (admits {shape})" if not atomic else ""))


if __name__ == "__main__":
    main()
