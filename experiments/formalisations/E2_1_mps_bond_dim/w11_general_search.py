"""
S206: General block-triangular search for W=11 with rows in [0, 11).

Tries all partitions of 11 into parts of size ≤ 3 (so we use only
det_fin_one, det_fin_two, det_fin_three from mathlib).

Strategy: greedy block-peeling DP.
For each partition shape (n_1, ..., n_K) (sizes ≤ 3):
  - Place a "first block" of size n_1: pick n_1 rows R_1 and n_1 cols C_1
    such that the n_1 × n_1 submatrix M[R_1, C_1] is invertible.
  - For block-TRIANGULAR ordering, the FIRST block is upper-left in the
    permuted matrix; the rows in subsequent blocks must have 0 on the cols
    of this block (so the lower-left of the next "larger" sub-matrix is
    zero relative to C_1).

Equivalently and more efficiently: enumerate all chains
  C_1 ⊊ C_1 ∪ C_2 ⊊ ... ⊊ {0..10}
of col subsets where successive C_k are k-step extensions, and each
"layer" admits enough rows with support entirely in C_1 ∪ ... ∪ C_k.

We actually do a backward search: pick the LAST block first (bottom-right
of the permuted matrix). Rows in the last block have all their nonzeros
in the last block's cols. Then peel.

USAGE:
    python3 w11_general_search.py
"""
from itertools import combinations
from math import gcd
from sympy import isprime, totient, Matrix


def build_matrix(W, row_pool):
    """Return the row_pool × W matrix M[r][c] = chiP(r*W + c + 1)."""
    return [[1 if isprime(r * W + c + 1) else 0 for c in range(W)]
            for r in range(row_pool)]


def row_support(M, r, cols):
    """Set of col-indices in `cols` where M[r][c] != 0."""
    return frozenset(c for c in cols if M[r][c] != 0)


def search_block_triangular(M, all_rows, all_cols, partition,
                            existing_helpers=None, max_to_keep=200):
    """Search for a block-triangular decomposition of M[all_rows, all_cols]
    with block sizes given by `partition` (an ordered tuple, top-to-bottom
    of the block-triangular decomposition).

    Block-triangular structure:
      Permuted matrix has blocks B_1 (top-left), B_2, ..., B_K (bottom-right)
      Block B_i has size partition[i-1] × partition[i-1].
      Lower-left blocks (B_i × B_j with i > j) are zero.

    Equivalently:
      - For block_K (bottom-right): rows R_K have support ⊆ C_K.
      - For block_{K-1}: rows R_{K-1} have support ⊆ C_{K-1} ∪ C_K.
      - ...
      - For block_1: rows R_1 have support ⊆ all cols.

    We search by peeling from the BOTTOM (last block first):
      - Pick C_K (size partition[-1]) and R_K (size partition[-1]) such
        that for r ∈ R_K, support(r) ⊆ C_K, and M[R_K, C_K] is invertible.
      - Recurse on remaining rows × remaining cols with partition[:-1].
    """
    K = len(partition)
    n_total = sum(partition)
    if len(all_rows) < n_total or len(all_cols) < n_total:
        return []

    candidates = []

    def recurse(rows_remaining, cols_remaining, partition_remaining,
                R_chain_so_far, C_chain_so_far):
        if len(candidates) >= max_to_keep:
            return
        if not partition_remaining:
            candidates.append((tuple(R_chain_so_far), tuple(C_chain_so_far)))
            return
        # Pick the LAST block of the remaining partition.
        n_last = partition_remaining[-1]
        # Eligible rows: those with support (within cols_remaining)
        # of size ≤ n_last that lies entirely in cols_remaining (trivially
        # true) AND we'll pick C_last as a superset of these supports.
        # Strategy: enumerate C_last (n_last cols), then find R_last rows
        # whose support is ⊆ C_last.
        for C_last in combinations(cols_remaining, n_last):
            C_last_set = set(C_last)
            # Eligible rows for R_last: support ⊆ C_last_set, AND
            # support intersects C_last_set (else row is all-zero on
            # C_last → no contribution).
            eligible_rows = [r for r in rows_remaining
                             if all(M[r][c] == 0 for c in cols_remaining
                                    if c not in C_last_set)]
            if len(eligible_rows) < n_last:
                continue
            for R_last in combinations(eligible_rows, n_last):
                # Check submatrix M[R_last, C_last] is invertible
                MM = Matrix([[M[r][c] for c in C_last] for r in R_last])
                if MM.det() == 0:
                    continue
                # Recurse
                new_rows = [r for r in rows_remaining if r not in R_last]
                new_cols = [c for c in cols_remaining if c not in C_last_set]
                recurse(new_rows, new_cols, partition_remaining[:-1],
                        R_chain_so_far + [R_last],
                        C_chain_so_far + [tuple(C_last)])
                if len(candidates) >= max_to_keep:
                    return

    recurse(list(all_rows), list(all_cols), list(partition), [], [])
    return candidates


def primes_used_by_blocks(M, R_chain, C_chain, W):
    """All primes that contribute as 1-entries in the blocks (NOT in upper-
    right zones; only the chosen blocks)."""
    primes = set()
    for R_block, C_block in zip(R_chain, C_chain):
        for r in R_block:
            for c in C_block:
                if M[r][c]:
                    primes.add(r * W + c + 1)
    return primes


def all_primes_in_chain(M, R_chain, C_chain, W):
    """All primes (as 1-entries) anywhere in the union of all rows × all
    cols (so block-triangular UPPER-RIGHT entries are included; needed
    for Lean closure since we still must compute those entries' values)."""
    primes = set()
    all_R = sum([list(R) for R in R_chain], [])
    all_C = sum([list(C) for C in C_chain], [])
    for r in all_R:
        for c in all_C:
            if M[r][c]:
                primes.add(r * W + c + 1)
    return primes


if __name__ == "__main__":
    EXISTING_HELPERS = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
                        47, 53, 59, 61, 89, 97, 109, 127, 149, 179, 199,
                        211, 241, 293, 337}

    W = 11
    row_pool = 11  # j = 1: rows must be in [0, W^j) = [0, 11)
    M = build_matrix(W, row_pool)
    all_rows = list(range(row_pool))
    all_cols = list(range(W))

    # Try several partition shapes (descending preference order):
    PARTITIONS = [
        (5, 5, 1),
        (1, 5, 5),
        (5, 1, 5),
        (5, 3, 3),
        (3, 5, 3),
        (3, 3, 5),
        (5, 4, 2),
        (5, 2, 4),
        (4, 5, 2),
        (2, 5, 4),
        (4, 4, 3),
        (4, 3, 4),
        (3, 4, 4),
    ]

    best_per_partition = {}
    for partition in PARTITIONS:
        cands = search_block_triangular(M, all_rows, all_cols, partition,
                                        existing_helpers=EXISTING_HELPERS,
                                        max_to_keep=2000)
        if cands:
            # Score by min new primes
            scored = []
            for R_chain, C_chain in cands:
                primes = all_primes_in_chain(M, R_chain, C_chain, W)
                new_primes = primes - EXISTING_HELPERS
                scored.append((len(new_primes), sorted(new_primes),
                               R_chain, C_chain))
            scored.sort(key=lambda x: x[0])
            best_per_partition[partition] = scored[0]
        print(f"  partition {partition}: {len(cands)} candidates "
              f"{'(best new primes: ' + str(best_per_partition.get(partition, (0,))[0]) + ')' if cands else ''}")

    print()
    print("=== Triangulability of each block in best per-partition ===")
    for partition, (n_new, new_primes, R_chain, C_chain) in sorted(best_per_partition.items(),
                                                                   key=lambda kv: (kv[1][0], len(kv[0])))[:6]:
        print(f"\nPartition {partition} (R_chain bottom-first):")
        for i, (R, C) in enumerate(zip(R_chain, C_chain)):
            n = len(R)
            print(f"  block {i} (size {n}): rows={R}, cols={C}")
            if n >= 4:
                # Check if this block admits leading-row upper triangulation
                sigs = []
                for r in R:
                    s = 0
                    for k, c in enumerate(C):
                        if M[r][c] != 0:
                            s |= (1 << k)
                    sigs.append((r, s))
                # DP search for triangulation
                full_mask = (1 << n) - 1
                reach = [False] * (1 << n)
                pred = [-1] * (1 << n)
                reach[0] = True
                for T in range(1, 1 << n):
                    for c_bit in range(n):
                        if not (T & (1 << c_bit)):
                            continue
                        T_prev = T & ~(1 << c_bit)
                        if not reach[T_prev]:
                            continue
                        if any((sg & T_prev) == 0 and (sg & (1 << c_bit))
                               for _, sg in sigs):
                            reach[T] = True
                            pred[T] = c_bit
                            break
                if reach[full_mask]:
                    # Reconstruct
                    sigma_order = []
                    T = full_mask
                    while T != 0:
                        cb = pred[T]
                        sigma_order.append(cb)
                        T &= ~(1 << cb)
                    sigma_order.reverse()
                    rho_used = []
                    used = set()
                    placed = 0
                    for cb in sigma_order:
                        for ri, (r, sg) in enumerate(sigs):
                            if ri in used:
                                continue
                            if (sg & placed) == 0 and (sg & (1 << cb)):
                                rho_used.append(r)
                                used.add(ri)
                                placed |= (1 << cb)
                                break
                    rho = rho_used
                    sigma = [C[cb] for cb in sigma_order]
                    print(f"    triangulable: row-order={rho}, col-order={sigma}")
                    diag_primes = [rho[i] * 11 + sigma[i] + 1 for i in range(n)]
                    print(f"    diagonal primes: {diag_primes}")
                else:
                    print(f"    NOT triangulable: needs cofactor expansion")
            elif n <= 3:
                MM = Matrix([[M[r][c] for c in C] for r in R])
                print(f"    det = {MM.det()} (size {n} → mathlib det_fin_{n})")

    print()
    print("=== Summary by partition ===")
    sorted_partitions = sorted(best_per_partition.items(),
                               key=lambda kv: (kv[1][0], len(kv[0])))
    for partition, (n_new, new_primes, R_chain, C_chain) in sorted_partitions[:5]:
        print(f"\nPartition: {partition}")
        print(f"  R_chain (top → bottom): {R_chain}")
        print(f"  C_chain (top → bottom): {C_chain}")
        print(f"  new primes needed: {new_primes} (count={n_new})")
        # Verify
        all_R = sum([list(R) for R in R_chain], [])
        all_C = sum([list(C) for C in C_chain], [])
        full = Matrix([[M[r][c] for c in all_C] for r in all_R])
        block_dets = []
        for R, C in zip(R_chain, C_chain):
            MM = Matrix([[M[r][c] for c in C] for r in R])
            block_dets.append(MM.det())
        print(f"  block dets: {block_dets}, product = {Matrix.det(full)}")
        # Check block-triangular
        for i_idx, (Ri, Ci) in enumerate(zip(R_chain, C_chain)):
            for j_idx, (Rj, Cj) in enumerate(zip(R_chain, C_chain)):
                if i_idx > j_idx:
                    nonzero = [(r, c) for r in Ri for c in Cj if M[r][c]]
                    if nonzero:
                        print(f"  WARNING: block ({i_idx}, {j_idx}) "
                              f"has nonzero entries: {nonzero}")
