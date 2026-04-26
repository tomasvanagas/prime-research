"""
Proposal 2 — PSLQ on delta(n) against an analytic feature dictionary.

Compute delta(n) = p(n) - round(R^{-1}(n)) for n = 1..n_max.
Build feature matrix F[n, :] of polylog-computable analytic quantities.
Run PSLQ on sliding windows of (delta, F) rows; check if a stable, small-norm
integer relation appears across windows.

VERDICT: >=50% of windows yield the same low-norm integer relation ⇒ escalate.

(Fast version: scipy.special.expi for Li, no mpmath bottleneck.)
"""
from __future__ import annotations
import math
from typing import List

import numpy as np
import mpmath as mp
from scipy.special import expi
from scipy.special import digamma as sp_digamma
from sympy import primerange
from sympy.functions.combinatorial.numbers import mobius


def li(x):
    """Logarithmic integral via Ei."""
    return expi(np.log(x))


# Mobius coefficients up to k=30 (enough for R(x))
MU = np.array([int(mobius(k)) for k in range(1, 31)], dtype=np.int64)


def riemann_R(x: float) -> float:
    s = 0.0
    for k in range(1, 31):
        m = MU[k - 1]
        if m == 0:
            continue
        s += m / k * li(x ** (1.0 / k))
    return float(s)


def riemann_R_array(xs: np.ndarray) -> np.ndarray:
    out = np.zeros_like(xs, dtype=np.float64)
    for k in range(1, 31):
        m = MU[k - 1]
        if m == 0:
            continue
        out += m / k * expi(np.log(xs) / k)
    return out


def riemann_R_inv_bisect(n_target: float, lo: float, hi: float) -> float:
    """Bisect to find R(x) = n_target."""
    f_lo = riemann_R(lo) - n_target
    f_hi = riemann_R(hi) - n_target
    while f_hi < 0:
        hi *= 2
        f_hi = riemann_R(hi) - n_target
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        f_mid = riemann_R(mid) - n_target
        if f_mid < 0:
            lo, f_lo = mid, f_mid
        else:
            hi, f_hi = mid, f_mid
    return 0.5 * (lo + hi)


def first_n_primes(n_max: int) -> List[int]:
    upper = int(n_max * (math.log(n_max) + math.log(math.log(n_max)) + 2)) + 100
    return list(primerange(2, upper))[:n_max]


def liouville_partial(N: int) -> List[int]:
    spf = [0] * (N + 1)
    for i in range(2, N + 1):
        if spf[i] == 0:
            for j in range(i, N + 1, i):
                if spf[j] == 0:
                    spf[j] = i
    out = [0] * (N + 1)
    for k in range(1, N + 1):
        m = k
        cnt = 0
        while m > 1:
            p = spf[m]
            while m % p == 0:
                m //= p
                cnt += 1
        out[k] = (-1) ** cnt
    return list(np.cumsum(out[1:]))


def pslq_relation(rows: np.ndarray, max_coeff: int = 200):
    cols = rows.shape[1]
    avg = rows.mean(axis=0)
    # Use mpmath PSLQ on column averages.
    mp.mp.dps = 30
    try:
        vec = mp.pslq([mp.mpf(float(x)) for x in avg], tol=1e-12, maxcoeff=max_coeff)
        if vec is None:
            return None, None
        v = np.array([float(x) for x in vec])
        residual = float(np.max(np.abs(rows @ v)))
        return v, residual
    except Exception:
        return None, None


def main(n_max: int = 2000):
    print(f"# PSLQ on delta(n), n=1..{n_max}")
    primes = first_n_primes(n_max)
    print(f"# p({n_max}) = {primes[-1]}")

    print("# Computing R^{-1}(n) for all n ...")
    rinv = np.zeros(n_max, dtype=np.float64)
    lo = 2.0
    for i, n in enumerate(range(1, n_max + 1)):
        # R(x) ~ x/log x; start near n*log(n+2)
        guess = max(2.0, n * math.log(n + 2))
        x = riemann_R_inv_bisect(n, max(lo, 2.0), guess * 4)
        rinv[i] = x
        lo = max(2.0, x * 0.9)
    print(f"# R^{{-1}}(1)={rinv[0]:.4f}, R^{{-1}}({n_max})={rinv[-1]:.4f}")

    delta = np.array([p - round(r) for p, r in zip(primes, rinv)], dtype=np.float64)
    print(f"# delta stats: min={delta.min()}, max={delta.max()}, "
          f"mean={delta.mean():.3f}, std={delta.std():.3f}")

    # Liouville partial sums up to max prime
    print("# Building Liouville partial sums ...")
    L_max = primes[-1] + 10
    Lpart = liouville_partial(L_max)

    # Feature dictionary
    print("# Building feature dictionary ...")
    F = np.zeros((n_max, 7), dtype=np.float64)
    for i, n in enumerate(range(1, n_max + 1)):
        p = primes[i]
        F[i, 0] = math.log(p)
        F[i, 1] = math.log(math.log(max(p, 3)))
        F[i, 2] = float(sp_digamma(n))
        F[i, 3] = 1.0 / math.log(p)
        F[i, 4] = 1.0
        F[i, 5] = float(Lpart[p - 1]) if p - 1 < len(Lpart) else 0.0
        F[i, 6] = float(n)

    win = 30
    step = 50
    print(f"# Sliding PSLQ: window={win}, step={step}")
    relations = []
    n_windows = 0
    for start in range(0, n_max - win, step):
        rows = np.column_stack([delta[start:start + win], F[start:start + win]])
        v, res = pslq_relation(rows, max_coeff=200)
        n_windows += 1
        if v is not None:
            relations.append((start, v, res))

    print(f"\n## PSLQ summary across {n_windows} windows: "
          f"{len(relations)} non-null returns")

    if not relations:
        print("VERDICT: NO relations found. PROPOSAL 2 FAILS at n<=", n_max)
        return

    # Group by sign-canonicalized integer vector
    sigs = {}
    for start, v, res in relations:
        nz = [x for x in v if x != 0]
        if not nz:
            continue
        sign = 1.0 if nz[0] > 0 else -1.0
        key = tuple(np.round(v * sign).astype(int))
        sigs.setdefault(key, []).append((start, res))

    print(f"\n## Distinct relation signatures: {len(sigs)}")
    items = sorted(sigs.items(), key=lambda kv: -len(kv[1]))
    print("  Top 5 signatures (vec  -> count, mean residual):")
    for key, hits in items[:5]:
        mean_res = float(np.mean([r for _, r in hits]))
        print(f"    {key} -> count={len(hits)}, mean_res={mean_res:.3e}")

    most_common, hits = items[0]
    frac = len(hits) / n_windows
    print(f"\n## Most common signature appears in {frac*100:.1f}% of windows")
    if frac >= 0.5:
        print("VERDICT: PSLQ relation stable across windows ⇒ "
              "PROPOSAL 2 PROMISING (escalate).")
    else:
        print("VERDICT: relations unstable across windows ⇒ PROPOSAL 2 FAILS.")
    print("\n# DONE")


if __name__ == "__main__":
    main(n_max=2000)
