"""
S144: Comprehensive enumeration of the leading-row triangulation family for
the orthogonal-corner case `(W, d = j + 1)` of E2.1 (MPS bond dimension).

For each W in [2, 72], we ask: does there exist a row permutation `rho` and a
column permutation `sigma` (with rows in [0, W) and columns in [0, W) chosen
from `live_cols(W) ∪ {dead}` for some dead column with `chiP(dead+1) = 1`)
such that `(unfolding W (j+1) j).submatrix rho sigma` is upper triangular
with `1` on the diagonal?

This is the precondition of the existing W ∈ {2, 3, 4, 5, 6, 8, 12, 18, 20}
Lean closures: each closure exhibits such (rho, sigma).

The search is dynamic-programming based: signatures of rows over the chosen
R cols form bitmasks, and we DP over subsets to decide whether a permutation
exists. The DP runs in O(2^R · R · W) time per dead-col candidate.

Parameters scanned: R = phi(W) + 1 up to 22 (i.e., 2^R ≤ 4 · 10^6 states).

USAGE:
    python3 leading_row_search.py
"""
from math import gcd
from sympy import isprime, totient


def find_triangulation_dp(W, R_cap=22):
    """
    Return (dead, rho, sigma) if a triangulation exists, otherwise None.
    `rho` is a list of row indices in [0, W); `sigma` is a list of col
    indices in [0, W).  The R x R submatrix `unfolding(rho_i, sigma_k)`
    (with rho_i in [0, W^j) embedded for any j >= 1 and entry =
    chiP(rho_i * W + sigma_k + 1)) is upper triangular with 1's on the
    diagonal under the canonical Fin-ordering.
    """
    live_cols = [k for k in range(W) if gcd(k + 1, W) == 1]
    dead_candidates = [k for k in range(W)
                       if gcd(k + 1, W) != 1 and isprime(k + 1)]
    R = int(totient(W)) + 1

    if not dead_candidates or R > R_cap:
        return None

    for dead in dead_candidates:
        cols_set = sorted(live_cols + [dead])
        sig = []
        for r in range(W):
            s = 0
            for k, c in enumerate(cols_set):
                if isprime(r * W + c + 1):
                    s |= (1 << k)
            sig.append(s)

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
                       for r in range(W)):
                    reach[T] = True
                    pred[T] = c
                    break

        if not reach[full]:
            continue

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
            for r in range(W):
                if r in used:
                    continue
                if (sig[r] & placed) == 0 and (sig[r] & (1 << c_idx)):
                    picked = r
                    break
            if picked < 0:
                rho_rows = None
                break
            rho_rows.append(picked)
            used.add(picked)
            placed |= (1 << c_idx)

        if rho_rows is None:
            continue

        return (dead, rho_rows, [cols_set[ci] for ci in sigma_order])

    return None


if __name__ == "__main__":
    closed = {2, 3, 4, 5, 6, 8, 12, 18, 20}
    R_cap = 22

    print("Leading-row triangulation enumeration for E2.1 orthogonal corner")
    print(f"Range scanned: W in [2, 72]; R_cap = {R_cap}")
    print()
    print(f"{'W':>3} | {'phi':>4} | {'R':>3} | result")
    print("-" * 50)

    closed_found, obstructed_found, skipped = [], [], []
    for W in range(2, 73):
        phi = int(totient(W))
        R = phi + 1
        if R > R_cap:
            print(f"{W:>3} | {phi:>4} | {R:>3} | SKIP (R > {R_cap})")
            skipped.append((W, R))
            continue

        result = find_triangulation_dp(W, R_cap=R_cap)
        if result is None:
            tag = "OBSTRUCTED"
            obstructed_found.append(W)
        else:
            dead, rho, sigma = result
            diag = [rho[i] * W + sigma[i] + 1 for i in range(R)]
            tag = f"FOUND (dead={dead}, max_diag={max(diag)})"
            closed_found.append(W)
        marker = " *" if W in closed else ""
        print(f"{W:>3} | {phi:>4} | {R:>3} | {tag}{marker}")

    print()
    print(f"Closed in Lean (S98..S143): {sorted(closed)}")
    print(f"Algorithmically reachable: {sorted(closed_found)}")
    print(f"Structurally obstructed:   {sorted(obstructed_found)}")
    print(f"Skipped (R > {R_cap}):          {sorted(W for W, _ in skipped)}")
    print()
    new_obstructions = sorted(set(obstructed_found) -
                              {7, 9, 10, 11, 14, 15, 16, 24, 30})
    print(f"NEW obstructions from S144 search: {new_obstructions}")
