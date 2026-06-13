"""
S206: Identify a `1 + 1 + 3 + 3 + 3` BlockTriangular permutation for the
W=11 orthogonal-corner case (`(W=11, d=j+1)`, R = phi(11) + 1 = 11) of E2.1.

W=11 is prime, so live cols are {0..9} (every k+1 in {1..10} is coprime to
11) and dead col k=10 (k+1 = 11, prime, gcd = 11). chiP(11) = 1 at row 0.
R = 11.

S144 leading_row_search confirmed W=11 has NO leading-row + dead-col
upper-triangulation in rows [0, 11). This pre-search asks for the next
best thing: a 1 + 1 + 3 + 3 + 3 BlockTriangular structure, which is the
cleanest single-session Lean target -- it composes only `det_fin_one`
and `det_fin_three` and nests `det_fromBlocks_zero_21` four times.

Partition shape:
  block X:  1×1  (row 0, dead col)        — chiP(11) = 1
  block Y:  1×1  (some row r1, some col)  — chiP(p) = 1 for some prime p
  block A:  3×3  (rows A_rows, cols A_cols)
  block B:  3×3  (rows B_rows, cols B_cols)
  block C:  3×3  (rows C_rows, cols C_cols)
Block-DIAGONAL means all off-diagonal blocks zero.

USAGE:
    python3 w11_blocktriangular_search.py [--row-pool 11|22|33] \\
                                          [--shape 1+1+3+3+3|1+3+3+1+3]
"""
import sys
from itertools import combinations, permutations
from math import gcd
from sympy import isprime, totient, Matrix


def find_block_diagonal_11335(W=11, row_pool=11, max_to_keep=200):
    """
    Search for `1 + 1 + 3 + 3 + 3` block-DIAGONAL decompositions.

    Returns a list of candidates {rho, sigma, det, primes_used}, sorted
    by (number of new prime helpers needed, total det magnitude).
    """
    live = [k for k in range(W) if gcd(k + 1, W) == 1]
    dead_candidates = [k for k in range(W)
                       if gcd(k + 1, W) != 1 and isprime(k + 1)]
    if not dead_candidates:
        # W prime: dead col is k = W - 1 (k+1 = W).
        dead = W - 1
    else:
        dead = dead_candidates[0]

    # All R = phi(W) + 1 columns we will use (sorted ascending by k).
    cols_set = sorted(live + [dead])
    R = len(cols_set)
    assert R == int(totient(W)) + 1

    rows_avail = list(range(0, row_pool))
    M = [[1 if isprime(r * W + c + 1) else 0 for c in cols_set]
         for r in range(row_pool)]

    # row 0 has the dead-col witness chiP(W) = 1; we always place row 0 at
    # block X (1×1 block).
    X_row = 0
    X_col_idx = cols_set.index(dead)

    # row 0 must have ZERO entries on every other column (for block-diagonal).
    if any(M[0][c] != 0 for c in range(R) if c != X_col_idx):
        # row 0 has chiP(prime+1) at multiple cols (e.g. for W=11, primes
        # 2,3,5,7 are at cols 1,2,4,6 too) -- so row 0's row in the matrix
        # has multiple 1s. Block-DIAGONAL requires row 0 to be 0 outside
        # X-col. This means we cannot use row 0 for the dead-col 1×1 block
        # in a block-DIAGONAL split. Pivot to a BlockTriangular argument:
        # row 0 is in block X, but UPPER off-diagonal entries in row 0 are
        # allowed (since det_fromBlocks_zero_21 only requires LOWER-LEFT
        # zero). The "dead col is zero on rows other than 0" requirement
        # still must hold for block-triangular too.
        pass

    # To ensure each row used outside row 0 has zero in the dead column
    # (so that the dead-col contributes only the row-0 1×1 block), we
    # restrict to rows where M[r][X_col_idx] == 0.
    candidate_rows = [r for r in rows_avail
                      if r != X_row and M[r][X_col_idx] == 0]

    # Also, we want to find a SECOND 1×1 block (block Y) at some row r1
    # and col c1 with M[r1][c1] = 1, AND for block-DIAGONAL we need:
    #   - row r1 has 0s on all other cols (in cols_set), OR we accept
    #     block-TRIANGULAR which only needs LOWER-LEFT off-diagonal zero.
    #
    # We'll search for block-DIAGONAL first (cleaner Lean assembly).

    print(f"W={W}, row_pool={row_pool}, R={R}, dead={dead}, "
          f"candidate_rows={len(candidate_rows)}")
    print(f"row 0 nonzero col indices in cols_set: "
          f"{[c for c in range(R) if M[0][c] != 0]}")
    print()

    # Enumerate Y = (r1, c1) singleton pairs. For BD, M[r1][c1] = 1, and
    # ideally row r1 is zero on every other col (in our cols_set), and
    # col c1 is zero on every other row in candidate_rows.
    bd_candidates = []
    other_col_idxs = [c for c in range(R) if c != X_col_idx]
    for r1 in candidate_rows:
        for c1 in other_col_idxs:
            if M[r1][c1] != 1:
                continue
            # row r1 must be 0 on cols_set \ {c1, X_col_idx} for BD
            row_others = [c for c in other_col_idxs if c != c1]
            row_r1_other_nz = [c for c in row_others if M[r1][c] != 0]
            if row_r1_other_nz:
                continue
            # Now we have rows for blocks A, B, C from candidate_rows \ {r1}
            # and cols for blocks A, B, C from other_col_idxs \ {c1}.
            ABC_rows_pool = [r for r in candidate_rows if r != r1]
            ABC_cols_pool = [c for c in other_col_idxs if c != c1]
            # We need 9 rows out of len(ABC_rows_pool), and partition the
            # 9 cols (= len(ABC_cols_pool) since R - 2 = 9 for R=11) into
            # three 3-sets.
            if len(ABC_cols_pool) != 9:
                continue
            # Also, col c1 must be 0 on all rows in ABC_rows_pool (for BD)
            col_c1_other_nz = [r for r in ABC_rows_pool if M[r][c1] != 0]
            ABC_rows_pool = [r for r in ABC_rows_pool
                             if M[r][c1] == 0]

            # Finally for true BD we need col X_col_idx zero on all
            # ABC rows (already true by candidate_rows filter).
            # And row X_row (=0) zero on all ABC cols (need to enforce).
            row0_ABC_nz = [c for c in ABC_cols_pool if M[0][c] != 0]
            if row0_ABC_nz:
                # Not block-diagonal at row 0; would need block-triangular.
                continue

            # Search 3+3+3 partition of 9 cols × 9 rows for block-DIAGONAL
            # of three 3×3 blocks (each invertible).
            for A_cols in combinations(ABC_cols_pool, 3):
                rest_cols_after_A = [c for c in ABC_cols_pool if c not in A_cols]
                for B_cols in combinations(rest_cols_after_A, 3):
                    C_cols = tuple(c for c in rest_cols_after_A
                                   if c not in B_cols)
                    # For block-DIAGONAL with three 3×3 blocks, we partition
                    # rows ABC_rows_pool into three 3-sets too.
                    if len(ABC_rows_pool) < 9:
                        continue
                    # We're searching all 9-element subsets and all 3+3+3
                    # partitions; for efficiency use the "rows where col
                    # set X has any 1" filter heuristic below.
                    # For BD, A_rows × (B_cols + C_cols) must be all zero,
                    # B_rows × (A_cols + C_cols) all zero, C_rows × (A_cols
                    # + B_cols) all zero.
                    # Rows that are zero on (B_cols + C_cols) are A-eligible.
                    A_eligible = [r for r in ABC_rows_pool
                                  if all(M[r][c] == 0 for c in B_cols + C_cols)]
                    B_eligible = [r for r in ABC_rows_pool
                                  if all(M[r][c] == 0 for c in A_cols + C_cols)]
                    C_eligible = [r for r in ABC_rows_pool
                                  if all(M[r][c] == 0 for c in A_cols + B_cols)]
                    if len(A_eligible) < 3 or len(B_eligible) < 3 or len(C_eligible) < 3:
                        continue

                    # Brute force: try every disjoint triple from each
                    # eligibility set.
                    for A_rows in combinations(A_eligible, 3):
                        # Check A block invertible
                        MA = Matrix([[M[r][c] for c in A_cols] for r in A_rows])
                        dA = MA.det()
                        if dA == 0:
                            continue
                        rest_after_A = [r for r in B_eligible if r not in A_rows]
                        for B_rows in combinations(rest_after_A, 3):
                            MB = Matrix([[M[r][c] for c in B_cols] for r in B_rows])
                            dB = MB.det()
                            if dB == 0:
                                continue
                            rest_after_AB = [r for r in C_eligible
                                             if r not in A_rows and r not in B_rows]
                            for C_rows in combinations(rest_after_AB, 3):
                                MC = Matrix([[M[r][c] for c in C_cols] for r in C_rows])
                                dC = MC.det()
                                if dC == 0:
                                    continue
                                # Verify Y col is 0 on (A_rows ∪ B_rows ∪ C_rows)
                                ABC_rows = A_rows + B_rows + C_rows
                                if any(M[r][c1] != 0 for r in ABC_rows):
                                    continue
                                # Verify X col is 0 on (A_rows ∪ B_rows ∪ C_rows)
                                if any(M[r][X_col_idx] != 0 for r in ABC_rows):
                                    continue
                                bd_candidates.append({
                                    'rho': (X_row, r1) + tuple(A_rows) + tuple(B_rows) + tuple(C_rows),
                                    'sigma_idx': (X_col_idx, c1) + tuple(A_cols) + tuple(B_cols) + tuple(C_cols),
                                    'sigma_kvals': (cols_set[X_col_idx], cols_set[c1])
                                                   + tuple(cols_set[c] for c in A_cols)
                                                   + tuple(cols_set[c] for c in B_cols)
                                                   + tuple(cols_set[c] for c in C_cols),
                                    'detX': 1, 'detY': 1,
                                    'detA': dA, 'detB': dB, 'detC': dC,
                                    'X_row': X_row, 'X_col_idx': X_col_idx,
                                    'Y_row': r1, 'Y_col_idx': c1,
                                    'A_rows': A_rows, 'A_cols': A_cols,
                                    'B_rows': B_rows, 'B_cols': B_cols,
                                    'C_rows': C_rows, 'C_cols': C_cols,
                                    'partition': '1+1+3+3+3',
                                    'cols_set': cols_set,
                                })
                                if len(bd_candidates) >= max_to_keep:
                                    return bd_candidates
    return bd_candidates


def primes_used_by_candidate(cand, W=11):
    """Compute the set of primes appearing as `1` entries in the submatrix."""
    cols_set = cand['cols_set']
    rho = cand['rho']
    sigma_idx = cand['sigma_idx']
    primes = set()
    for i, r in enumerate(rho):
        for k, c_idx in enumerate(sigma_idx):
            n = r * W + cols_set[c_idx] + 1
            if isprime(n):
                primes.add(n)
    return primes


def verify_candidate(cand, W=11):
    """Construct the full 11×11 submatrix and verify it is block-diagonal
    with the expected determinant."""
    rho = cand['rho']
    sigma_idx = cand['sigma_idx']
    cols_set = cand['cols_set']
    M = Matrix([[1 if isprime(rho[i] * W + cols_set[sigma_idx[k]] + 1) else 0
                 for k in range(11)] for i in range(11)])
    expected = (cand['detX'] * cand['detY'] * cand['detA']
                * cand['detB'] * cand['detC'])
    return M, M.det(), expected


if __name__ == "__main__":
    EXISTING_HELPERS = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
                        47, 53, 59, 61, 89, 97, 109, 127, 149, 179, 199,
                        211, 241, 293, 337}

    row_pool = 11
    print(f"=== W=11 block-DIAGONAL (1+1+3+3+3) search, row_pool={row_pool} ===")
    print()
    cands = find_block_diagonal_11335(W=11, row_pool=row_pool, max_to_keep=2000)
    print(f"Found {len(cands)} candidates.")
    if not cands:
        print("\nNo candidates in [0, 11); retrying with [0, 22).")
        cands = find_block_diagonal_11335(W=11, row_pool=22, max_to_keep=2000)
        print(f"Found {len(cands)} candidates in [0, 22).")
    if not cands:
        print("\nNo candidates in [0, 22); retrying with [0, 33).")
        cands = find_block_diagonal_11335(W=11, row_pool=33, max_to_keep=2000)
        print(f"Found {len(cands)} candidates in [0, 33).")

    if not cands:
        print("\nNo (1+1+3+3+3) BD found.")
        sys.exit(0)

    # Score by number of new prime helpers needed.
    scored = []
    for c in cands:
        primes = primes_used_by_candidate(c)
        new_primes = primes - EXISTING_HELPERS
        scored.append((len(new_primes), sorted(new_primes), c))
    scored.sort(key=lambda x: (x[0], abs(x[2]['detA'] * x[2]['detB'] * x[2]['detC'])))

    min_new = scored[0][0]
    print(f"\nMinimum new chiP helpers needed: {min_new}")
    print()
    print("Top candidates:")
    seen = 0
    for n_new, new_primes, c in scored:
        if seen >= 8:
            break
        if n_new > min_new + 1:
            break
        print(f"  rho   = {c['rho']}")
        print(f"  sigma = {c['sigma_kvals']}")
        print(f"    detX={c['detX']}, detY={c['detY']}, "
              f"detA={c['detA']}, detB={c['detB']}, detC={c['detC']}")
        print(f"    total det = {c['detX']*c['detY']*c['detA']*c['detB']*c['detC']}")
        print(f"    new primes: {new_primes} (count={n_new})")
        primes_used = primes_used_by_candidate(c)
        print(f"    all primes used: {sorted(primes_used)}")
        seen += 1
        print()

    # Verify the top candidate.
    cand = scored[0][2]
    M, computed_det, expected_det = verify_candidate(cand)
    print("=== Verification of top candidate ===")
    print("Full 11×11 submatrix:")
    for row in M.tolist():
        print(f"  {row}")
    print(f"computed det = {computed_det}")
    print(f"expected det (from blocks) = {expected_det}")
    assert computed_det == expected_det, "block-product != actual det!"
    print("MATCH.")
