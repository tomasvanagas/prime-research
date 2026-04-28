#!/usr/bin/env python3
"""S182 — Independent check of A_W = lim (1/log Q) sum_{sqf q <= Q} 1/phi(q).

S168 claimed A_W = 1 by Selberg-Delange, with empirical 1.04 at Q=5000.
Verify this directly at Q up to 10^7. If A_W is NOT converging to 1
(or converging to a value materially different from 1), the entire
'21%' chain falls apart.

Also: directly compute asymptote by fitting

    sum_{sqf q <= Q} 1/phi(q) = A * log(Q) + B + epsilon(Q)

and ask: is the slope of (sum - log Q) vs (something) consistent with A=1?
"""

import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def linear_sieve_squarefree_inv_phi(Q):
    """Return arrays:
       sqf[q] = is q squarefree (bool)
       phi[q] = euler totient
       cumulative_sum_inv_phi over sqf q.
    """
    smallest_prime = np.zeros(Q + 1, dtype=np.int64)
    for i in range(2, Q + 1):
        if smallest_prime[i] == 0:
            for j in range(i, Q + 1, i):
                if smallest_prime[j] == 0:
                    smallest_prime[j] = i

    phi = np.zeros(Q + 1, dtype=np.int64)
    sqf = np.ones(Q + 1, dtype=bool)
    phi[1] = 1
    for q in range(2, Q + 1):
        n = q
        result = q
        while n > 1:
            p = int(smallest_prime[n])
            result = result // p * (p - 1)
            mult = 0
            while n % p == 0:
                n //= p
                mult += 1
            if mult > 1:
                sqf[q] = False
        phi[q] = result
    return phi, sqf


def main():
    Qs = [100, 500, 1000, 5000, 10_000, 50_000, 100_000, 500_000, 1_000_000, 5_000_000]
    print(f"{'Q':>10} {'sum 1/phi over sqf':>22} {'A_W = sum/log(Q)':>18} {'sum - log Q (= B if A=1)':>26}")

    Q_max = max(Qs)
    print(f"Sieving phi/sqf up to {Q_max}...")
    phi, sqf = linear_sieve_squarefree_inv_phi(Q_max)
    print("Done.")

    inv_phi = np.where(sqf[1:], 1.0 / phi[1:], 0.0)
    cum = np.cumsum(inv_phi)  # cum[i] = sum over q in [1, i+1]

    rows = []
    for Q in Qs:
        s = float(cum[Q - 1])  # cum index Q-1 = sum over q in [1, Q]
        A = s / math.log(Q)
        diff = s - math.log(Q)
        rows.append((Q, s, A, diff))
        print(f"{Q:>10} {s:>22.6f} {A:>18.6f} {diff:>26.6f}")

    # Fit sum = A log Q + B over the high-Q tail
    log_Qs = np.array([math.log(Q) for Q in Qs])
    sums = np.array([r[1] for r in rows])
    # Use last 5 points
    idx = slice(5, None)
    A_fit, B_fit = np.polyfit(log_Qs[idx], sums[idx], 1)
    print()
    print(f"Linear fit on Q in {Qs[5]}..{Qs[-1]}: sum = {A_fit:.6f} * log(Q) + {B_fit:.6f}")
    print(f"  S168 claim: A = 1.0 (Selberg-Delange).")
    print(f"  This fit:    A = {A_fit:.6f}.")
    if abs(A_fit - 1.0) < 0.05:
        print(f"  CONSISTENT with A=1 (within 5%).")
    else:
        print(f"  REFUTE: |A - 1| = {abs(A_fit - 1.0):.4f} > 0.05.")

    # Also fit using all 10 points
    A_all, B_all = np.polyfit(log_Qs, sums, 1)
    print(f"Linear fit on all 10 Qs: A = {A_all:.6f}, B = {B_all:.6f}")

    # And: if we forced through (large-Q-asymptote, B), compute residual at small Q
    print()
    print("Predicted vs observed under A=1:")
    print(f"{'Q':>10} {'observed sum':>14} {'1*log(Q)':>10} {'observed - log(Q)':>20}")
    for Q, s, A, diff in rows:
        print(f"{Q:>10} {s:>14.6f} {math.log(Q):>10.4f} {diff:>20.6f}")
    # If A really is 1, the differences should converge to a constant B as Q -> inf.
    # Compute their second differences to see if they're stabilizing
    diffs = np.array([r[3] for r in rows])
    print()
    print("If A_W = 1 strictly, (sum - log Q) should converge to B (Wirsing's constant).")
    print(f"  observed differences: {diffs}")
    print(f"  successive deltas:    {np.diff(diffs)}")
    final_delta = diffs[-1] - diffs[-2]
    print(f"  last delta (Q={Qs[-2]} -> {Qs[-1]}): {final_delta:.6f}")
    if abs(final_delta) < 0.05:
        print(f"  CONSISTENT with stabilizing (last delta < 0.05).")
    else:
        print(f"  Still drifting — A may not equal 1 exactly.")

    import json
    out_path = HERE / "wirsing_results.json"
    out_path.write_text(json.dumps({
        "rows": [{"Q": Q, "sum": s, "A_W_partial": A, "B_implied": d} for Q, s, A, d in rows],
        "fit_high_tail": {"A": float(A_fit), "B": float(B_fit)},
        "fit_all": {"A": float(A_all), "B": float(B_all)},
    }, indent=2))
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
