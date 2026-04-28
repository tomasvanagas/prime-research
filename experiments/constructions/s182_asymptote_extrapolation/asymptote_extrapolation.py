#!/usr/bin/env python3
"""S182 verify probe — extend analytic side of S169's 21% claim to d=26..30
and ask: does cum(Q*=N^{0.21}) / (0.21 * pi(N)) actually trend to 1.0?

S169's existing data: ratio at d=14..24 = {1.330, 1.266, 1.260, 1.193,
1.172, 1.167}. Slow finite-N convergence claimed toward 1.0 (Wirsing-A
asymptote). Six points only; the asymptote is NOT directly observed.

This probe extends to d=26, 28, 30 using the *closed-form* main-term
(no Fourier sums needed). If the asymptote is genuinely 1.0, the ratio
should keep decreasing. If it stalls or asymptotes to a value > 1.0,
S168's "0.21 * pi(N)" is an OVERSTATEMENT — the true asymptote is
0.21 * (asymptotic ratio).

Closed-form main term (from S169.py):
    E_main(q, N) = mu(q)^2 * (pi(N) - omega(q))^2 / (phi(q) * N)
    cum(Q) = sum_{sqf q in [2, Q]} E_main(q, N)

For pi(N) at d in {26, 28, 30}, use a segmented sieve OR Cipolla.
We use a segmented sieve for d <= 28; the run is bounded.

Falsification test: fit ratio = a + b/log(N) + c/log(N)^2 to all 9
points. If best-fit a is significantly different from 1.0, S168's
asymptote claim is wrong.
"""

import math
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def factor(q):
    f = []
    n = q
    p = 2
    while p * p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            f.append((p, e))
        p += 1
    if n > 1:
        f.append((n, 1))
    return f


def is_squarefree(q):
    return q == 1 or all(e == 1 for _, e in factor(q))


def euler_phi(q):
    if q == 1:
        return 1
    res = q
    for p, _ in factor(q):
        res = res // p * (p - 1)
    return res


def omega(q):
    return len(factor(q))


def pi_N_via_sieve(N):
    """Direct sieve. OK up to N ~ 10^9 with bytearray memory."""
    sieve = bytearray(b"\x01") * (N + 1)
    sieve[0] = 0
    sieve[1] = 0
    for i in range(2, int(N ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = bytearray(len(sieve[i * i :: i]))
    return sum(sieve)


def cum_main(Q, pi_N, N):
    """Analytic main-term sum used by S169.

    cum(Q) = sum_{sqf q in [2, Q]} mu(q)^2 (pi(N) - omega(q))^2 / (phi(q) * N).
    """
    s = 0.0
    for q in range(2, Q + 1):
        if not is_squarefree(q):
            continue
        r = omega(q)
        s += (pi_N - r) ** 2 / (euler_phi(q) * N)
    return s


def wirsing_partial(Q):
    """Partial Wirsing sum: sum_{sqf q <= Q} 1/phi(q)."""
    s = 0.0
    for q in range(1, Q + 1):
        if is_squarefree(q):
            s += 1.0 / euler_phi(q)
    return s


def main():
    print("S182 — extend S169's analytic 21% test to high d")
    print("=" * 64)
    out = {}

    # Recompute S169's d=14..24 baseline (analytic only, fast)
    ds = [14, 16, 18, 20, 22, 24, 26, 28]
    rows = []
    print(f"{'d':>3} {'N':>12} {'pi(N)':>10} {'Q*':>5} {'cum/(0.21pi)':>15} {'A_W partial':>15} {'PNT factor':>12}")
    for d in ds:
        N = 2 ** d
        t0 = time.time()
        pi_N = pi_N_via_sieve(N)
        t_pi = time.time() - t0
        Q_star = max(2, int(round(N ** 0.21)))
        t0 = time.time()
        cum = cum_main(Q_star, pi_N, N)
        A_W = wirsing_partial(Q_star) / math.log(Q_star) if Q_star > 1 else float("nan")
        t_cum = time.time() - t0
        target = 0.21 * pi_N
        ratio = cum / target
        # PNT factor: pi(N) * log(N) / N
        pnt_factor = pi_N * math.log(N) / N
        rows.append({
            "d": d, "N": N, "pi_N": pi_N, "Q_star": Q_star,
            "cum": cum, "ratio": ratio, "A_W": A_W,
            "pnt_factor": pnt_factor,
        })
        print(f"{d:>3} {N:>12} {pi_N:>10} {Q_star:>5} {ratio:>15.4f} {A_W:>15.4f} {pnt_factor:>12.4f}  (pi:{t_pi:.1f}s cum:{t_cum:.1f}s)")
        out[d] = rows[-1]

    # Fit ratio = a + b/log(N) + c/log(N)^2
    log_Ns = np.array([math.log(r["N"]) for r in rows])
    ratios = np.array([r["ratio"] for r in rows])
    X = np.vstack([np.ones_like(log_Ns), 1.0 / log_Ns, 1.0 / log_Ns ** 2]).T
    coef, *_ = np.linalg.lstsq(X, ratios, rcond=None)
    print()
    print(f"Fit ratio = a + b/log(N) + c/log(N)^2")
    print(f"  a (asymptote) = {coef[0]:.5f}")
    print(f"  b              = {coef[1]:.5f}")
    print(f"  c              = {coef[2]:.5f}")
    pred = X @ coef
    print()
    print(f"{'d':>3} {'observed':>10} {'predicted':>10} {'residual':>12}")
    for r, p in zip(rows, pred):
        print(f"{r['d']:>3} {r['ratio']:>10.4f} {p:>10.4f} {r['ratio']-p:>+12.5f}")

    # Bootstrap-style: refit with each point dropped, see how stable a is
    print()
    print("Leave-one-out asymptote stability:")
    for i in range(len(rows)):
        idx = [j for j in range(len(rows)) if j != i]
        X_loo = X[idx]
        y_loo = ratios[idx]
        c_loo, *_ = np.linalg.lstsq(X_loo, y_loo, rcond=None)
        print(f"  drop d={rows[i]['d']:>2}: a = {c_loo[0]:.5f}")

    # Also fit with c=0 (linear)
    X2 = np.vstack([np.ones_like(log_Ns), 1.0 / log_Ns]).T
    coef2, *_ = np.linalg.lstsq(X2, ratios, rcond=None)
    print()
    print(f"Linear fit (c=0): a = {coef2[0]:.5f}, b = {coef2[1]:.5f}")

    # Forced a=1 (S168 prediction)
    X3 = (1.0 / log_Ns).reshape(-1, 1)
    coef3, *_ = np.linalg.lstsq(X3, ratios - 1.0, rcond=None)
    pred3 = 1.0 + X3 @ coef3
    rss3 = float(np.sum((ratios - pred3) ** 2))
    rss_unrestricted = float(np.sum((ratios - X @ coef) ** 2))
    print(f"Forced asymptote = 1.0:  b = {coef3[0]:.5f}, RSS = {rss3:.6f}")
    print(f"Unrestricted fit:        RSS = {rss_unrestricted:.6f}")

    # Save
    import json
    out_path = HERE / "asymptote_results.json"
    out_path.write_text(json.dumps({
        "rows": rows,
        "fit_full": {"a": float(coef[0]), "b": float(coef[1]), "c": float(coef[2])},
        "fit_linear": {"a": float(coef2[0]), "b": float(coef2[1])},
        "rss_forced_1": rss3,
        "rss_unrestricted": rss_unrestricted,
    }, indent=2))
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
