"""
S159: Identify a `1 + 3 + 3` BlockTriangular permutation for the W=7
orthogonal-corner case (`(W=7, d=j+1)`, R = phi(7) + 1 = 7) of E2.1.

W=7's 7×7 slab has live cols {0..5} (since 7 is prime, every k+1 ∈
{1..6} is coprime to 7) and dead col k=6 (k+1 = 7, prime). R = 7.

The S144 leading_row_search showed W=7 has NO leading-row + dead-col
upper-triangulation in rows [0, 7). This pre-search asks for the next
best thing: a 1+3+3 BlockTriangular structure (block-diagonal with the
dead-col witness chiP(7) at top-left, 1×1 block) where each diagonal
3×3 block has nonzero determinant.

USAGE:
    python3 w7_blocktriangular_search.py
"""
from itertools import combinations
from math import gcd
from sympy import isprime, totient, Matrix


def find_block_diagonal_decompositions(W=7, allow_extended_rows=False):
    """
    For the (W=7, d=j+1) corner, search 1+3+3 (block-DIAGONAL) splits.
    By default rows are in [1, W). Optionally allow rows up to 2*W-1
    if no decomposition is found in the basic regime.
    """
    live = [k for k in range(W) if gcd(k + 1, W) == 1]
    dead_candidates = [k for k in range(W)
                       if gcd(k + 1, W) != 1 and isprime(k + 1)]
    if not dead_candidates:
        # For W = prime (e.g. 7), all k in [0, W) have gcd(k+1, W) = 1
        # except k = W - 1 (where k+1 = W). The "dead col" is k = W-1
        # (a self-prime witness contributing chiP(W) = 1 at row 0).
        dead = W - 1  # k+1 = W (prime, but gcd = W ≠ 1)
    else:
        dead = dead_candidates[0]

    cols = sorted(live + [dead])
    dead_col_idx = cols.index(dead)
    live_col_idxs = [i for i in range(len(cols)) if i != dead_col_idx]

    R = int(totient(W)) + 1

    rows_avail = list(range(1, W if not allow_extended_rows else 2 * W))

    # Build the matrix on the available rows × chosen cols.
    M = [[1 if isprime(r * W + c + 1) else 0 for c in cols]
         for r in range(max(rows_avail) + 1)]

    candidates = []
    for rows6 in combinations(rows_avail, R - 1):
        for A_rows in combinations(rows6, 3):
            B_rows = tuple(r for r in rows6 if r not in A_rows)
            for A_cols in combinations(live_col_idxs, 3):
                B_cols = tuple(c for c in live_col_idxs if c not in A_cols)
                # Lower-left block (B_rows × A_cols) zero (block triangular)
                if any(M[r][c] != 0 for r in B_rows for c in A_cols):
                    continue
                MA = Matrix([[M[r][c] for c in A_cols] for r in A_rows])
                dA = MA.det()
                if dA == 0:
                    continue
                MB = Matrix([[M[r][c] for c in B_cols] for r in B_rows])
                dB = MB.det()
                if dB == 0:
                    continue
                # Upper-right block (A_rows × B_cols) zero too: block-diagonal
                ur_zero = all(M[r][c] == 0 for r in A_rows for c in B_cols)
                # Also check: dead-col rows must be zero except at row 0
                # (since the (1,0..6) and (0..6, 0) entries should match
                # the 1+6 fromBlocks structure via the dead col at row 0).
                # We already enforce row 0 → dead col by construction in
                # the rho/sigma assembly below; just verify dead col is
                # zero on rows in {A_rows ∪ B_rows}.
                dead_col_zero_on_rest = all(M[r][dead_col_idx] == 0
                                            for r in A_rows + B_rows)
                if not dead_col_zero_on_rest:
                    continue
                candidates.append({
                    'A_rows': A_rows, 'A_cols': A_cols,
                    'B_rows': B_rows, 'B_cols': B_cols,
                    'detA': dA, 'detB': dB, 'block_diagonal': ur_zero,
                    'A_kvals': [cols[i] for i in A_cols],
                    'B_kvals': [cols[i] for i in B_cols],
                    'dead': dead,
                })
    return candidates, cols, M


def primes_used_by_candidate(cand, cols, M, W=7):
    """Compute the set of primes in the (1+3+3) submatrix entries."""
    primes = set()
    # Row 0 entries (all 7 chosen cols)
    for c_idx in range(len(cols)):
        n = 0 * W + cols[c_idx] + 1
        if isprime(n):
            primes.add(n)
    # A block + A row-to-B-cols (will be zero in BD case but we still
    # compute the prime cells for any nonzero entries)
    for r in cand['A_rows']:
        for c_idx in cand['A_cols'] + cand['B_cols']:
            n = r * W + cols[c_idx] + 1
            if isprime(n):
                primes.add(n)
    # B block (B rows × B cols)
    for r in cand['B_rows']:
        for c_idx in cand['B_cols']:
            n = r * W + cols[c_idx] + 1
            if isprime(n):
                primes.add(n)
    return primes


if __name__ == "__main__":
    EXISTING_HELPERS = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
                        47, 53, 59, 61, 89, 97, 109, 127, 149, 179, 199, 211,
                        241, 293, 337}

    # First try rows in [0, W) only.
    candidates, cols, M = find_block_diagonal_decompositions(W=7, allow_extended_rows=False)
    print(f"Found {len(candidates)} valid (1+3+3) BlockTriangular "
          f"decompositions for W=7 corner (rows in [0, 7)).")
    if not candidates:
        print("\nNone in [0, 7); retrying with rows in [0, 14).")
        candidates, cols, M = find_block_diagonal_decompositions(W=7, allow_extended_rows=True)
        print(f"Found {len(candidates)} valid (1+3+3) decompositions in [0, 14).")
    print()

    if not candidates:
        print("No (1+3+3) BlockTriangular decomposition exists with rows in [0, 14).")
        print("Need to expand search further or try different block shapes.")
        exit(0)

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
    seen = 0
    for n_new, new_primes, cand in scored:
        if n_new > min_new:
            break
        rho = (0,) + cand['A_rows'] + cand['B_rows']
        sigma = (cand['dead'],) + tuple(cand['A_kvals']) + tuple(cand['B_kvals'])
        print(f"  rho   = {rho}")
        print(f"  sigma = {sigma}")
        print(f"    A-block dets: dA={cand['detA']}, dB={cand['detB']}")
        print(f"    block_diagonal={cand['block_diagonal']}")
        print(f"    new primes: {new_primes}")
        print()
        seen += 1
        if seen >= 5:
            print(f"  (... total {sum(1 for n,_,_ in scored if n == min_new)} at min={min_new}; stopping at 5)")
            break

    # Verify the top candidate gives det = ±1
    if scored:
        cand = scored[0][2]
        rho = [0] + list(cand['A_rows']) + list(cand['B_rows'])
        sigma = [cand['dead']] + list(cand['A_kvals']) + list(cand['B_kvals'])
        full = Matrix([[1 if isprime(rho[i] * 7 + sigma[k] + 1) else 0
                        for k in range(7)] for i in range(7)])
        print("=== Verification of top candidate ===")
        print("Full 7x7 submatrix:")
        for row in full.tolist():
            print(f"  {row}")
        print(f"det = {full.det()}")
