"""
S206: Test whether the 10×10 inner block of W=11 admits upper triangulation.

Strategy: take the outer 1+10 BlockTriangular split with row 0 + dead-col
witness for the 1×1 block, and ask whether rows [1, 11) × live-cols form
a 10×10 matrix that upper-triangulates.

If YES → the W=11 corner closes via outer fromBlocks_zero_21 (1×1 + 10×10
upper-triangular) — single-session Lean assembly via det_fin_one + a local
prod_univ_ten' lemma + det_of_upperTriangular.

If NO → fall back to nested block-triangular search.
"""
from math import gcd
from sympy import isprime, totient


def upper_triangulate(W=11, row_pool_max=22, exclude_row0=True, exclude_dead_col=True):
    """DP-based upper-triangulation search on rows [a, row_pool_max) × cols.

    a = 1 if exclude_row0 else 0; cols = live cols if exclude_dead_col else
    live + dead.
    """
    live = [k for k in range(W) if gcd(k + 1, W) == 1]
    dead_candidates = [k for k in range(W)
                       if gcd(k + 1, W) != 1 and isprime(k + 1)]
    if not dead_candidates:
        dead = W - 1
    else:
        dead = dead_candidates[0]
    cols_set = sorted(live) if exclude_dead_col else sorted(live + [dead])
    R = len(cols_set)
    a = 1 if exclude_row0 else 0
    rows_avail = list(range(a, row_pool_max))

    sig = []
    row_idx = []
    for r in rows_avail:
        s = 0
        for k, c in enumerate(cols_set):
            if isprime(r * W + c + 1):
                s |= (1 << k)
        if s != 0:
            sig.append(s)
            row_idx.append(r)

    full = (1 << R) - 1
    reach = [False] * (1 << R)
    pred = [-1] * (1 << R)
    reach[0] = True

    for T in range(1, 1 << R):
        for c in range(R):
            if not (T & (1 << c)):
                continue
            T_prev = T & ~(1 << c)
            if not reach[T_prev]:
                continue
            if any((sig[r] & T_prev) == 0 and (sig[r] & (1 << c))
                   for r in range(len(sig))):
                reach[T] = True
                pred[T] = c
                break

    if not reach[full]:
        return None

    sigma_order = []
    T = full
    while T != 0:
        c = pred[T]
        sigma_order.append(c)
        T &= ~(1 << c)
    sigma_order.reverse()

    rho_rows, used = [], set()
    placed = 0
    for c_idx in sigma_order:
        picked = -1
        for ri in range(len(sig)):
            if ri in used:
                continue
            if (sig[ri] & placed) == 0 and (sig[ri] & (1 << c_idx)):
                picked = ri
                break
        if picked < 0:
            return None
        rho_rows.append(row_idx[picked])
        used.add(picked)
        placed |= (1 << c_idx)

    return (dead, rho_rows, [cols_set[ci] for ci in sigma_order], cols_set)


def primes_in_submatrix(rho_rows, sigma_kvals, W):
    """All primes in the submatrix entries (where chiP = 1)."""
    primes = set()
    for i, r in enumerate(rho_rows):
        for k, c in enumerate(sigma_kvals):
            n = r * W + c + 1
            if isprime(n):
                primes.add(n)
    return primes


if __name__ == "__main__":
    print("=== W=11: inner 10×10 (rows ≥ 1, live cols) upper-triangulation ===")
    print()
    # Try with rows in [1, 11)
    result = upper_triangulate(W=11, row_pool_max=11,
                               exclude_row0=True, exclude_dead_col=True)
    if result is None:
        print("[1, 11) × live-cols (10×10): OBSTRUCTED")
    else:
        dead, rho, sigma, _ = result
        print(f"[1, 11) × live-cols (10×10): FOUND")
        print(f"  rho   = {rho}")
        print(f"  sigma = {sigma}")
        print(f"  diagonal entries (chiP applied to): {[rho[i]*11 + sigma[i] + 1 for i in range(len(rho))]}")
        diag = [rho[i]*11 + sigma[i] + 1 for i in range(len(rho))]
        print(f"  all diagonal primes: {sorted(diag)}")
        primes = primes_in_submatrix(rho, sigma, 11)
        EXISTING_HELPERS = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
                            47, 53, 59, 61, 89, 97, 109, 127, 149, 179, 199,
                            211, 241, 293, 337}
        new_primes = primes - EXISTING_HELPERS
        print(f"  primes used: {sorted(primes)}")
        print(f"  new primes needed: {sorted(new_primes)} (count={len(new_primes)})")
    print()

    # If [1, 11) fails, try [1, 22)
    if result is None:
        result = upper_triangulate(W=11, row_pool_max=22,
                                   exclude_row0=True, exclude_dead_col=True)
        if result is None:
            print("[1, 22) × live-cols: OBSTRUCTED")
        else:
            dead, rho, sigma, _ = result
            print(f"[1, 22) × live-cols: FOUND")
            print(f"  rho   = {rho}")
            print(f"  sigma = {sigma}")
            diag = [rho[i]*11 + sigma[i] + 1 for i in range(len(rho))]
            print(f"  all diagonal primes: {sorted(diag)}")
            primes = primes_in_submatrix(rho, sigma, 11)
            EXISTING_HELPERS = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
                                47, 53, 59, 61, 89, 97, 109, 127, 149, 179, 199,
                                211, 241, 293, 337}
            new_primes = primes - EXISTING_HELPERS
            print(f"  primes used: {sorted(primes)}")
            print(f"  new primes needed: {sorted(new_primes)} (count={len(new_primes)})")
        print()

    # Also test the simpler full triangulation hypothesis
    print("=== Sanity check: full [0, 11) × all-cols (11×11) (S144 said NO) ===")
    result = upper_triangulate(W=11, row_pool_max=11,
                               exclude_row0=False, exclude_dead_col=False)
    print("Result:", "FOUND" if result else "OBSTRUCTED")
