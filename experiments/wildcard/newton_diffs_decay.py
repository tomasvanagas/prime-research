"""
Newton forward-difference series for pi(x).

Idea: Any function f: Z>=0 -> Z admits the Newton series expansion
    f(x) = sum_{k>=0} C(x, k) * Delta^k f(0)
where Delta^k f(0) is the k-th forward difference at 0.
If the difference table Delta^k pi(0) has fast decay, sparsity, or a closed
form, then pi(x) is evaluable in O(polylog) time by truncating the series.

Truncation error is bounded by max |Delta^k pi(.)| for k > K. So we test:
- How fast does max_x |Delta^k pi(x)| decay (or grow) with k?
- Is the difference table sparse in any structured basis?
- For partial recovery: how many leading terms are needed for pi(x) +-1?

Published context: it is well known that for completely multiplicative or
periodic sequences, Newton series may converge fast; for primes this is
expected to *not* converge nicely, but we measure whether there is any
"effective" partial fast formula.
"""

import numpy as np
from sympy import primerange


def pi_table(N):
    """Return [pi(0), pi(1), ..., pi(N)]."""
    primes = set(primerange(2, N + 1))
    arr = np.zeros(N + 1, dtype=np.int64)
    c = 0
    for n in range(N + 1):
        if n in primes:
            c += 1
        arr[n] = c
    return arr


def difference_table(seq):
    """Return Delta^k seq for k=0..len(seq)-1, as a 2D array (rows = k)."""
    N = len(seq)
    table = np.zeros((N, N), dtype=object)  # use python ints to avoid overflow
    for j in range(N):
        table[0, j] = int(seq[j])
    for k in range(1, N):
        for j in range(N - k):
            table[k, j] = table[k - 1, j + 1] - table[k - 1, j]
    return table


def main():
    N = 256  # difference table is N x N; bigger N exponentiates Delta-magnitude
    seq = pi_table(N)
    print(f"pi(0..{N}) computed; pi({N}) = {seq[N]}")

    # Build difference table at x = 0
    primes_idx = [int(p) for p in primerange(2, N + 1)]

    # We only need Delta^k pi(0) for k = 0..N
    # Compute via the formula Delta^k f(0) = sum_{j=0}^k (-1)^(k-j) C(k,j) f(j)
    from math import comb

    K = min(N, 200)
    delta0 = []
    for k in range(K + 1):
        s = 0
        for j in range(k + 1):
            s += ((-1) ** (k - j)) * comb(k, j) * int(seq[j])
        delta0.append(s)

    abs_d = [abs(int(d)) for d in delta0]
    print(f"\nFirst 16 Delta^k pi(0):")
    for k in range(16):
        print(f"  k={k:3d}  Delta^k pi(0) = {delta0[k]}")

    # Compute the magnitude growth.
    print(f"\nMagnitude profile (log2 |Delta^k pi(0)|):")
    growth = []
    for k in range(0, K + 1, 8):
        g = 0.0 if abs_d[k] == 0 else float(np.log2(float(abs_d[k]))) if abs_d[k] < 1e300 else int(abs_d[k]).bit_length() - 1
        growth.append((k, abs_d[k], g))
        print(f"  k={k:3d}  |Delta^k|={abs_d[k]:>40d}  log2={float(g):7.2f}")

    # Truncation: compute pi(x) via partial Newton sum and see error
    print(f"\nNewton-truncation reconstruction error (truncate at K terms):")
    from math import comb as C

    test_xs = [10, 50, 100, 200]
    for K_use in [5, 10, 20, 40, 80, 160]:
        if K_use > K:
            continue
        errs = []
        for x in test_xs:
            if x > N:
                continue
            est = sum(int(delta0[k]) * C(x, k) for k in range(K_use + 1))
            errs.append(abs(est - int(seq[x])))
        print(f"  K_use={K_use:3d}  errors at x={test_xs}: {errs}")

    # Sparsity check: count of nonzero, log-magnitude statistics
    nz = sum(1 for d in delta0 if d != 0)
    print(f"\nNonzero count among first {K+1} differences: {nz}")
    print(f"Largest |Delta^k|: {max(abs_d):e}")
    print(f"Smallest nonzero: {min(d for d in abs_d if d != 0)}")

    # Effective decay rate via linear fit on log magnitudes
    ks = np.array([k for k in range(K + 1) if abs_d[k] > 0])
    ys = np.array([float(int(abs_d[k]).bit_length() - 1) for k in ks])
    if len(ks) > 5:
        slope, intercept = np.polyfit(ks, ys, 1)
        print(f"\nLinear fit log2|Delta^k| ~ {slope:.4f} * k + {intercept:.4f}")
        if slope > 0:
            print(f"  -> magnitudes GROW exponentially in k (factor {2**slope:.4f}/step).")
            print(f"  -> Newton truncation does NOT converge; we cannot get poly-log.")
        else:
            print(f"  -> magnitudes decay; partial-formula fast eval may be possible.")


if __name__ == "__main__":
    main()
