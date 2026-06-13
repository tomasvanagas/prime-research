"""
S235: Identify a `1 + 3 + 3` BlockTriangular permutation for the W=14
orthogonal-corner case `(W=14, d=j+1)`, R = phi(14) + 1 = 7, of E2.1.

Per S128/S137/S144 leading-row search, W=14 was ruled out for leading-
row + dead-col upper-triangulation (rows 2 and 5 of the 14×14 j=1 slab
have identical support pattern at the chosen 7 cols). It joins the
"block-triangular-required" set.

W=14 = 2·7 is composite. By the parity argument: for r >= 1,
M[r][c] = chiP(14r + c + 1). Since 14r is even, 14r + c + 1 is even iff
c is odd. Primes > 2 are odd, so for r >= 1 the only nonzero column
positions are even c. Live cols are k where gcd(k+1, 14) = 1:
    k+1 in {1, 3, 5, 9, 11, 13} -> k in {0, 2, 4, 8, 10, 12}  (all even)
Dead-col candidates (k+1 prime, gcd != 1):
    k+1 = 2 -> k = 1 (odd; row >=1 entries all zero by parity)
    k+1 = 7 -> k = 6 (even; row >=1 entries possibly nonzero)

This is a striking simplification vs W=7 / W=11: for W=14 the live
cols are ALL EVEN, and rows >= 1 contribute nonzero only at EVEN
cols. So the parity argument that broke W=11 (forcing an irreducible
5x5 odd block) does NOT bite here — there is no odd live col, so
no "odd block" can form.

Per S206 results doc, W=14 was suggested as the next single-session
target precisely because composite W avoids the W=11-style parity
atomicity. This pre-search verifies that and finds a viable
(1 + 3 + 3) decomposition.

USAGE:
    python3 w14_blocktriangular_search.py
"""
from itertools import combinations
from math import gcd
from sympy import isprime, totient, Matrix


def find_block_diagonal_decompositions(W=14, dead_choice=1, max_row=None):
    """
    For the (W=14, d=j+1) corner, search 1+3+3 BlockTriangular splits.
    `dead_choice` selects which dead-col witness to use (k+1 = 2 -> k=1,
    or k+1 = 7 -> k=6).
    `max_row` is the upper bound on rows considered (default W).
    """
    live = [k for k in range(W) if gcd(k + 1, W) == 1]
    dead_candidates = [k for k in range(W)
                       if gcd(k + 1, W) != 1 and isprime(k + 1)]
    assert dead_choice in dead_candidates, (
        f"dead_choice={dead_choice} not in {dead_candidates}")
    dead = dead_choice

    cols = sorted(live + [dead])
    dead_col_idx = cols.index(dead)
    live_col_idxs = [i for i in range(len(cols)) if i != dead_col_idx]

    R = int(totient(W)) + 1
    assert R == 7, f"Expected R=7 for W=14, got R={R}"

    if max_row is None:
        max_row = W
    rows_avail = list(range(1, max_row))

    M = [[1 if isprime(r * W + c + 1) else 0 for c in cols]
         for r in range(max_row)]

    # Enforce: dead col must be zero on rows 1..max_row-1 (else outer
    # 1+6 BlockTriangular fails). Verify upfront.
    for r in rows_avail:
        if M[r][dead_col_idx] != 0:
            print(f"  [warn] dead col k={dead} nonzero at row r={r} "
                  f"(value = chiP({r*W+dead+1})). Skipping this dead.")
            return [], cols, M

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
                ur_zero = all(M[r][c] == 0 for r in A_rows for c in B_cols)
                candidates.append({
                    'A_rows': A_rows, 'A_cols': A_cols,
                    'B_rows': B_rows, 'B_cols': B_cols,
                    'detA': dA, 'detB': dB, 'block_diagonal': ur_zero,
                    'A_kvals': [cols[i] for i in A_cols],
                    'B_kvals': [cols[i] for i in B_cols],
                    'dead': dead,
                })
    return candidates, cols, M


def primes_used_by_candidate(cand, cols, M, W=14):
    """Compute the set of primes in the (1+3+3) submatrix entries."""
    primes = set()
    for c_idx in range(len(cols)):
        n = 0 * W + cols[c_idx] + 1
        if isprime(n):
            primes.add(n)
    for r in cand['A_rows']:
        for c_idx in cand['A_cols'] + cand['B_cols']:
            n = r * W + cols[c_idx] + 1
            if isprime(n):
                primes.add(n)
    for r in cand['B_rows']:
        for c_idx in cand['B_cols'] + cand['A_cols']:
            n = r * W + cols[c_idx] + 1
            if isprime(n):
                primes.add(n)
    return primes


if __name__ == "__main__":
    W = 14
    EXISTING_HELPERS = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
                        47, 53, 59, 61, 89, 97, 109, 127, 149, 179, 199, 211,
                        241, 293, 337}

    print(f"=== W={W} BlockTriangular pre-search ===")
    print(f"R = phi({W}) + 1 = {int(totient(W)) + 1}")
    live = [k for k in range(W) if gcd(k + 1, W) == 1]
    print(f"Live cols (residues coprime to {W}): {live}")
    dead_candidates = [k for k in range(W)
                       if gcd(k + 1, W) != 1 and isprime(k + 1)]
    print(f"Dead-col candidates (k+1 prime, gcd != 1): {dead_candidates}")
    print()

    all_results = {}
    for dead in dead_candidates:
        print(f"--- Dead col k={dead} (chiP({dead + 1}) = 1) ---")
        candidates, cols, M = find_block_diagonal_decompositions(
            W=W, dead_choice=dead, max_row=W)
        print(f"  rows in [0, {W}): {len(candidates)} valid (1+3+3) "
              f"BlockTriangular decompositions")
        if not candidates:
            candidates, cols, M = find_block_diagonal_decompositions(
                W=W, dead_choice=dead, max_row=2 * W)
            print(f"  rows in [0, {2 * W}): {len(candidates)} valid "
                  f"(1+3+3) decompositions")
        all_results[dead] = (candidates, cols, M)
        print()

    # Pick best candidate across both dead choices, scored by min new helpers,
    # then by block_diagonal preference.
    scored = []
    for dead, (candidates, cols, M) in all_results.items():
        for cand in candidates:
            primes = primes_used_by_candidate(cand, cols, M)
            new_primes = primes - EXISTING_HELPERS
            scored.append((len(new_primes),
                           0 if cand['block_diagonal'] else 1,
                           sorted(new_primes), cand, cols, M))
    scored.sort(key=lambda x: (x[0], x[1]))

    if not scored:
        print("=== No (1+3+3) BlockTriangular decomposition found ===")
        print("Need to widen search or try a different block shape.")
        exit(0)

    min_new = scored[0][0]
    bd_count = sum(1 for s in scored if s[1] == 0)
    print(f"Total candidates: {len(scored)}")
    print(f"  block-diagonal (preferred):     {bd_count}")
    print(f"  block-triangular only:           {len(scored) - bd_count}")
    print(f"Minimum new chiP helpers needed: {min_new}")
    print()

    print("Top candidates (minimum new helpers):")
    seen = 0
    for n_new, bd_score, new_primes, cand, cols, M in scored:
        if n_new > min_new:
            break
        rho = (0,) + cand['A_rows'] + cand['B_rows']
        sigma = (cand['dead'],) + tuple(cand['A_kvals']) + tuple(cand['B_kvals'])
        bd_label = 'block-DIAGONAL' if cand['block_diagonal'] else 'block-triangular'
        print(f"  rho   = {rho}")
        print(f"  sigma = {sigma}")
        print(f"    type: {bd_label}")
        print(f"    A-block det = {cand['detA']}, B-block det = {cand['detB']}")
        print(f"    new primes: {new_primes}")
        print()
        seen += 1
        if seen >= 5:
            print(f"  (... showing first 5 of {sum(1 for s in scored if s[0] == min_new)} at min={min_new})")
            break

    # Verify the top candidate gives nonzero det
    if scored:
        _, _, _, cand, cols, M = scored[0]
        rho = [0] + list(cand['A_rows']) + list(cand['B_rows'])
        sigma = [cand['dead']] + list(cand['A_kvals']) + list(cand['B_kvals'])
        full = Matrix([[1 if isprime(rho[i] * W + sigma[k] + 1) else 0
                        for k in range(7)] for i in range(7)])
        print("=== Verification of top candidate ===")
        print("Full 7x7 submatrix:")
        for row in full.tolist():
            print(f"  {row}")
        print(f"det = {full.det()}")
        print(f"rho   = {rho}")
        print(f"sigma = {sigma}")

        # Also list composites that need h_not_prime witnesses.
        composites = set()
        for i in range(7):
            for k in range(7):
                n = rho[i] * W + sigma[k] + 1
                if not isprime(n) and n > 1:
                    composites.add(n)
        composites.add(1)
        print(f"\nComposites to witness (entries with chiP = 0): "
              f"{sorted(composites)}")
        primes = set()
        for i in range(7):
            for k in range(7):
                n = rho[i] * W + sigma[k] + 1
                if isprime(n):
                    primes.add(n)
        print(f"Primes appearing in submatrix: {sorted(primes)}")
        new_primes = primes - EXISTING_HELPERS
        print(f"New chiP helpers needed: {sorted(new_primes)}")
