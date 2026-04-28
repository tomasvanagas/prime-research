"""
S151: Identify a `1 + 3 + 3` BlockTriangular permutation for the W=9
orthogonal-corner case (`(W=9, d=j+1)`, R = phi(9) + 1 = 7) of E2.1.

The S144 leading_row_search exhausted the leading-row + dead-col upper-
triangulation route and confirmed W=9 is structurally obstructed for
*upper triangular* under any row/col permutation (rows in [0, 9)).

This script searches for a coarser BlockTriangular structure: a
permutation that places the dead-col witness `chiP 3` at position
(0, 0) (a 1x1 top-left block), and partitions the remaining 6 rows
and 6 live columns into two blocks of size 3 such that the matrix
becomes BLOCK-DIAGONAL (both off-diagonal 3x3 blocks are zero).

Each diagonal 3x3 block having nonzero det then yields:
    det(7x7) = 1 * det(B1) * det(B2)

via two nested `Matrix.det_fromBlocks_zero_21` reductions in Lean.

USAGE:
    python3 w9_blocktriangular_search.py

S151 OUTPUT: 32 candidates exist; the minimum-new-prime-helpers
candidate uses 4 new primes {13, 41, 53, 61} and gives
    det(B1) = -1, det(B2) = -1, total det = 1.

The chosen permutation:
    rho ↦ (0, 1, 3, 5, 2, 4, 6)
    sigma ↦ (2, 1, 3, 7, 0, 4, 6)
"""
from itertools import combinations
from math import gcd
from sympy import isprime, totient, Matrix


def find_block_diagonal_decompositions(W=9):
    """
    For the (W, d=j+1) orthogonal corner, find all permutations of
    {0..W-1} rows and the 7 chosen cols (live ∪ {dead_2}) that yield a
    BlockTriangular structure (rho_0 = 0, dead-col = first col, then
    a 3+3 block-diagonal split of the remaining 6 rows × 6 live cols).
    """
    live = [k for k in range(W) if gcd(k + 1, W) == 1]
    dead_candidates = [k for k in range(W)
                       if gcd(k + 1, W) != 1 and isprime(k + 1)]
    if not dead_candidates:
        return []
    dead = dead_candidates[0]  # k=2 for W=9
    cols = sorted(live + [dead])
    dead_col_idx = cols.index(dead)
    live_col_idxs = [i for i in range(len(cols)) if i != dead_col_idx]

    M = [[1 if isprime(r * W + c + 1) else 0 for c in cols] for r in range(W)]
    R = int(totient(W)) + 1

    candidates = []
    rows_avail = list(range(1, W))  # rows 1..W-1 (row 0 is forced to dead col)
    for rows6 in combinations(rows_avail, R - 1):
        for A_rows in combinations(rows6, 3):
            B_rows = tuple(r for r in rows6 if r not in A_rows)
            for A_cols in combinations(live_col_idxs, 3):
                B_cols = tuple(c for c in live_col_idxs if c not in A_cols)
                # Lower-block constraint: M[B_rows, A_cols] = 0
                if any(M[r][c] != 0 for r in B_rows for c in A_cols):
                    continue
                # Block A determinant nonzero
                MA = Matrix([[M[r][c] for c in A_cols] for r in A_rows])
                dA = MA.det()
                if dA == 0:
                    continue
                # Block B determinant nonzero
                MB = Matrix([[M[r][c] for c in B_cols] for r in B_rows])
                dB = MB.det()
                if dB == 0:
                    continue
                # Also check upper-right block can be zero (block-DIAGONAL,
                # stronger than just upper-block-triangular):
                ur_zero = all(M[r][c] == 0 for r in A_rows for c in B_cols)
                candidates.append({
                    'A_rows': A_rows, 'A_cols': A_cols,
                    'B_rows': B_rows, 'B_cols': B_cols,
                    'detA': dA, 'detB': dB, 'block_diagonal': ur_zero,
                    'A_kvals': [cols[i] for i in A_cols],
                    'B_kvals': [cols[i] for i in B_cols],
                })
    return candidates, cols, M


def primes_used_by_candidate(cand, cols, M, W=9):
    """Compute the set of primes that appear in the (1+3+3) submatrix."""
    primes = set()
    # Row 0 entries (all 7 chosen cols)
    for c_idx in range(len(cols)):
        n = 0 * W + cols[c_idx] + 1
        if isprime(n):
            primes.add(n)
    # A block + A row entries to B cols (above-diagonal of 4x4 top-left
    # in 1+3+3 = both upper-block-triangular case)
    for r in cand['A_rows']:
        for c_idx in cand['A_cols']:
            n = r * W + cols[c_idx] + 1
            if isprime(n):
                primes.add(n)
        for c_idx in cand['B_cols']:
            n = r * W + cols[c_idx] + 1
            if isprime(n):
                primes.add(n)
    # B block entries (lower-right 3x3)
    for r in cand['B_rows']:
        for c_idx in cand['B_cols']:
            n = r * W + cols[c_idx] + 1
            if isprime(n):
                primes.add(n)
    return primes


if __name__ == "__main__":
    EXISTING_HELPERS = {2, 3, 5, 7, 11, 17, 19, 23, 29, 31, 37, 43, 47, 59,
                        89, 97, 109, 127, 149, 179, 199, 211, 241, 293, 337}

    candidates, cols, M = find_block_diagonal_decompositions(W=9)
    print(f"Found {len(candidates)} valid (1+3+3) BlockTriangular "
          f"decompositions for W=9 corner.")
    print()

    # Score by minimum new chiP-helpers needed
    scored = []
    for cand in candidates:
        primes = primes_used_by_candidate(cand, cols, M)
        new_primes = primes - EXISTING_HELPERS
        scored.append((len(new_primes), sorted(new_primes), cand))
    scored.sort(key=lambda x: x[0])

    min_new = scored[0][0]
    print(f"Minimum new chiP helpers needed: {min_new}")
    print(f"Block-diagonal-only candidates: "
          f"{sum(1 for _,_,c in scored if c['block_diagonal'])}")
    print()

    # Show the top candidates
    print("Top candidates (minimum new helpers):")
    for n_new, new_primes, cand in scored:
        if n_new > min_new:
            break
        rho = (0,) + cand['A_rows'] + cand['B_rows']
        sigma = (cols[2],) + tuple(cand['A_kvals']) + tuple(cand['B_kvals'])
        print(f"  rho = {rho}, sigma = {sigma}")
        print(f"    A-block dets: dA={cand['detA']}, dB={cand['detB']}")
        print(f"    block_diagonal={cand['block_diagonal']}")
        print(f"    new primes: {new_primes}")
        print()

    # Verify the chosen permutation gives det = 1
    print("=== Verification of chosen permutation ===")
    rho = [0, 1, 3, 5, 2, 4, 6]
    sigma = [2, 1, 3, 7, 0, 4, 6]
    full = Matrix([[1 if isprime(rho[i] * 9 + sigma[k] + 1) else 0
                    for k in range(7)] for i in range(7)])
    print("Full 7x7 submatrix:")
    for row in full.tolist():
        print(f"  {row}")
    print(f"det = {full.det()}")
