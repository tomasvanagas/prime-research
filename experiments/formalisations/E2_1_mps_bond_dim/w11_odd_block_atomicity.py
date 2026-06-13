"""
S206: Verify the W=11 odd 5×5 block (rows {1, 3, 5, 7, 9} × cols
{1, 3, 5, 7, 9}) is BLOCK-IRREDUCIBLE — no block-triangular decomposition
exists for ANY ordered partition with blocks of size ≤ 4.

This is the STRUCTURAL OBSTRUCTION to single-session Lean closure of the
W=11 corner case via the standard `(1 + W) + det_fromBlocks_zero_21` route.
"""
from itertools import combinations
from sympy import isprime, Matrix


def build_odd_block():
    """The 5×5 odd block of W=11."""
    W = 11
    odd_rows = (1, 3, 5, 7, 9)
    odd_cols = (1, 3, 5, 7, 9)
    return Matrix([[1 if isprime(r * W + c + 1) else 0 for c in odd_cols]
                   for r in odd_rows]), odd_rows, odd_cols


def check_blocktriangular(M, partition):
    """Search for block-triangular decomposition of M with given partition.
    Returns (R_chain, C_chain) if found, else None.

    R_chain[0] is the BOTTOM block (peeled first).
    """
    n = M.shape[0]
    assert sum(partition) == n
    all_rows = list(range(n))
    all_cols = list(range(n))
    found = []

    def recurse(rows_remaining, cols_remaining, partition_remaining,
                R_chain, C_chain):
        if found:
            return
        if not partition_remaining:
            found.append((tuple(R_chain), tuple(C_chain)))
            return
        n_last = partition_remaining[-1]
        for C_last in combinations(cols_remaining, n_last):
            C_last_set = set(C_last)
            other_cols = [c for c in cols_remaining if c not in C_last_set]
            eligible_rows = [r for r in rows_remaining
                             if all(M[r, c] == 0 for c in other_cols)]
            if len(eligible_rows) < n_last:
                continue
            for R_last in combinations(eligible_rows, n_last):
                MM = Matrix([[M[r, c] for c in C_last] for r in R_last])
                if MM.det() == 0:
                    continue
                new_rows = [r for r in rows_remaining if r not in R_last]
                new_cols = [c for c in cols_remaining if c not in C_last_set]
                recurse(new_rows, new_cols, partition_remaining[:-1],
                        R_chain + [R_last], C_chain + [C_last])
                if found:
                    return

    recurse(all_rows, all_cols, list(partition), [], [])
    return found[0] if found else None


def all_partitions(n, max_part=4):
    """All ordered partitions (compositions) of n with parts ≤ max_part."""
    if n == 0:
        yield ()
        return
    for first in range(1, min(n, max_part) + 1):
        for rest in all_partitions(n - first, max_part):
            yield (first,) + rest


if __name__ == "__main__":
    M, odd_rows, odd_cols = build_odd_block()
    print(f"Odd 5×5 block of W=11 (rows={odd_rows}, cols={odd_cols}):")
    for row in M.tolist():
        print(f"  {row}")
    print(f"det = {M.det()}, rank = {M.rank()}")
    print()

    print("=== Checking ALL ordered partitions of 5 with parts ≤ 4 ===")
    print()
    found_any = False
    for partition in all_partitions(5, max_part=4):
        if len(partition) == 1:
            continue  # trivial single-block "decomposition"
        result = check_blocktriangular(M, partition)
        status = "FOUND" if result else "obstructed"
        if result:
            found_any = True
        print(f"  partition {partition}: {status}")
        if result:
            print(f"    R_chain (bottom-first): {result[0]}")
            print(f"    C_chain (bottom-first): {result[1]}")
    print()
    if not found_any:
        print("ATOMIC: the W=11 odd 5×5 block admits NO non-trivial")
        print("block-triangular decomposition with blocks of size ≤ 4.")
    else:
        print("Some decomposition exists — see above.")
