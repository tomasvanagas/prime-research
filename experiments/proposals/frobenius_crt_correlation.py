"""Proposal D: Frobenius-trace ensemble correlation.

For several elliptic curves E and primes q <= N, compute
    a_q(E) = q + 1 - |E(F_q)|       (Frobenius trace, |a_q| <= 2 sqrt(q))

Build cumulative sums S_x(E) = sum_{q <= x prime} a_q(E) and a few
nonlinear features (squares, cubes). Then test whether `pi(x) mod m`
for small m can be recovered from these features by *any* linear
combination via least squares.

A finding "pi(x) mod m can be predicted at >> chance from a small set
of L-coefficient cumulative sums" would suggest a CRT-based modular
algorithm. A negative finding closes the path.

Curves are computed naively (no Schoof needed for x <= 10000) — the
test asks whether the *information* is there, regardless of how it's
extracted.

Run: python3 frobenius_crt_correlation.py
"""
from __future__ import annotations

import math
import sys
import time

import numpy as np


def sieve_primes(N: int) -> list[int]:
    sv = [True] * (N + 1)
    sv[:2] = [False, False]
    for k in range(2, int(math.isqrt(N)) + 1):
        if sv[k]:
            for j in range(k * k, N + 1, k):
                sv[j] = False
    return [i for i, b in enumerate(sv) if b]


def count_points_on_curve(a: int, b: int, p: int) -> int:
    """Naive |E(F_p)| where E: y^2 = x^3 + a x + b mod p, p odd prime."""
    if p == 2:
        # Different equation form needed; skip and return 0 (caller handles)
        return p + 1  # not meaningful for our tests
    # Build quadratic-residue indicator
    qr = [False] * p
    for k in range(p):
        qr[(k * k) % p] = True
    count = 1  # point at infinity
    for x in range(p):
        rhs = (x * x % p * x + a * x + b) % p
        if rhs == 0:
            count += 1
        elif qr[rhs]:
            count += 2
    return count


def ap(a: int, b: int, p: int) -> int:
    return p + 1 - count_points_on_curve(a, b, p)


def main() -> None:
    t0 = time.time()
    print("# Frobenius-trace ensemble correlation with pi(x) mod m")
    N = 5000
    primes = sieve_primes(N)
    print(f"# Total primes <= {N}: {len(primes)}")

    # Curves to use; avoid singular ones (need 4a^3 + 27 b^2 != 0)
    curves: list[tuple[int, int]] = []
    for a, b in [(0, 1), (1, 0), (-1, 1), (1, -1), (2, 3), (-3, 4), (5, 7)]:
        if (4 * a**3 + 27 * b**2) != 0:
            curves.append((a, b))
    print(f"# Curves used: {curves}")
    print()

    # Compute a_q for each curve at each odd prime
    odd_primes = [q for q in primes if q > 2]
    a_q_table = np.zeros((len(curves), len(odd_primes)), dtype=np.int64)
    for ci, (a, b) in enumerate(curves):
        for qi, q in enumerate(odd_primes):
            try:
                a_q_table[ci, qi] = ap(a, b, q)
            except Exception:
                a_q_table[ci, qi] = 0
        if ci == 0:
            print(f"# Curve {ci} = y^2 = x^3 + {a}x + {b}: a_q[:10] = {a_q_table[ci, :10].tolist()}")
    print(f"# a_q computed in {time.time()-t0:.1f}s")
    print()

    # Build feature matrix at "checkpoints" x = odd_primes themselves.
    # Row i corresponds to evaluation at odd_primes[i].
    # Features: cumulative sums of a_q, a_q^2, a_q (mod m) etc.
    n_check = len(odd_primes)
    feats: list[np.ndarray] = []
    feat_names: list[str] = []

    # Cumulative sums per curve
    for ci, (a, b) in enumerate(curves):
        cs = np.cumsum(a_q_table[ci])
        feats.append(cs.astype(np.float64))
        feat_names.append(f"S_aq_E{ci}")
        cs2 = np.cumsum(a_q_table[ci] ** 2)
        feats.append(cs2.astype(np.float64))
        feat_names.append(f"S_aq2_E{ci}")
        # Cumulative sum of (a_q mod m) for various m
        for m in (2, 3, 5):
            cs_m = np.cumsum((a_q_table[ci] % m)).astype(np.float64)
            feats.append(cs_m)
            feat_names.append(f"S_aq_mod{m}_E{ci}")

    # Universal features
    feats.append(np.array([math.log(max(2, q)) for q in odd_primes]))
    feat_names.append("log_q")
    feats.append(np.arange(1, n_check + 1, dtype=np.float64))
    feat_names.append("idx")

    F = np.array(feats).T  # shape (n_check, F_dim)
    print(f"# Feature matrix shape: {F.shape}")

    # Targets: pi(odd_primes[i]) mod m  =  (i+2) mod m  (since odd_primes[i]
    # is the (i+2)-nd prime if we count 2 as first prime... let me redo)
    # pi(odd_primes[i]) = i + 2 because primes[0]=2, primes[1]=3=odd_primes[0],
    # so pi(odd_primes[i]) = i + 2 with our indexing.
    pi_at = np.array([i + 2 for i in range(n_check)])
    print()

    for m in (2, 3, 5, 7):
        target = pi_at % m
        target_centered = target.astype(np.float64) - target.mean()
        # Train on first 70%, test on last 30%
        split = int(0.7 * n_check)
        F_tr, F_te = F[:split], F[split:]
        y_tr, y_te = target[:split], target[split:]

        # Standardize features
        mu = F_tr.mean(axis=0)
        sigma = F_tr.std(axis=0)
        sigma[sigma == 0] = 1.0
        Fs_tr = (F_tr - mu) / sigma
        Fs_te = (F_te - mu) / sigma

        # Solve least squares for real-valued target (then round)
        # Use ridge regression to avoid overfitting
        ridge = 1e-3
        XtX = Fs_tr.T @ Fs_tr + ridge * np.eye(Fs_tr.shape[1])
        Xty = Fs_tr.T @ y_tr.astype(np.float64)
        coef = np.linalg.solve(XtX, Xty)
        pred_real = Fs_te @ coef
        pred = np.mod(np.round(pred_real).astype(np.int64), m)

        accuracy = float(np.mean(pred == y_te))
        # Chance baseline: most-frequent class
        from collections import Counter

        train_counts = Counter(y_tr.tolist())
        majority_class, _ = train_counts.most_common(1)[0]
        chance = float(np.mean(y_te == majority_class))
        uniform = 1.0 / m

        print(
            f"## pi(x) mod {m}: train n={split}, test n={n_check-split}; "
            f"accuracy={accuracy:.3f}  majority={chance:.3f}  uniform={uniform:.3f}"
        )

    print(f"\n# Total elapsed: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
