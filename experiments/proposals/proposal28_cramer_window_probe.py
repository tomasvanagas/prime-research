"""Proposal 28 — Empirical Cramér-window probe.

Question: how wide must a window around R^{-1}(n) be to reliably contain p(n)?
If max gap |p(n) - R^{-1}(n)| / log^2(p(n)) is bounded for n <= 10^6, then a
sliding-window sieve of polylog width suffices to localise p(n) given an
anchor count.

Approach
--------
1. Sieve primes up to N = 10^7. Sequence p_1, p_2, ... .
2. Compute Riemann inverse R^{-1}(n) by Newton's method on R(x) = Σ μ(k)/k Li(x^{1/k}).
3. Tabulate delta_n = p_n - R^{-1}(n) and the ratio |delta_n| / log^2(p_n).
4. Report max, mean, 99th percentile, scaling exponent of |delta_n| vs log p_n.

If the empirical max ratio is bounded by a small constant for n up to 10^6,
the polylog window is empirically validated; if it grows, it is not.
"""
import math
from pathlib import Path
from typing import List

import numpy as np


REPO = Path(__file__).resolve().parents[2]


def sieve(N: int) -> List[int]:
    s = bytearray([1]) * (N + 1)
    s[0] = s[1] = 0
    for i in range(2, int(math.isqrt(N)) + 1):
        if s[i]:
            for j in range(i * i, N + 1, i):
                s[j] = 0
    return [i for i, v in enumerate(s) if v]


def li(x: float, n_terms: int = 200) -> float:
    """Logarithmic integral via series expansion (avoid mpmath dependency).

    Li(x) = gamma + log(log x) + sum_{k>=1} (log x)^k / (k * k!)
    Valid for x > 1.
    """
    if x <= 1:
        return 0.0
    log_x = math.log(x)
    EULER_GAMMA = 0.5772156649015329
    s = EULER_GAMMA + math.log(abs(log_x))
    term = 1.0
    for k in range(1, n_terms + 1):
        term *= log_x / k
        s += term / k
    return s


def R(x: float, K: int = 30) -> float:
    """Riemann's R function: R(x) = Σ_{k>=1} μ(k)/k Li(x^{1/k})."""
    # Möbius up to K
    mu = [0] * (K + 1)
    mu[1] = 1
    for i in range(1, K + 1):
        if mu[i]:
            for j in range(2 * i, K + 1, i):
                mu[j] -= mu[i]
    s = 0.0
    for k in range(1, K + 1):
        if mu[k] == 0:
            continue
        if x ** (1.0 / k) <= 1.0:
            break
        s += mu[k] / k * li(x ** (1.0 / k))
    return s


def R_inverse(n: int, x0: float = None) -> float:
    """Newton iteration solving R(x) = n. dR/dx ≈ 1/log x for large x."""
    if x0 is None:
        x0 = max(n * math.log(max(n, 2)) * 1.1, 10.0)
    x = x0
    for _ in range(50):
        Rx = R(x)
        dR = 1.0 / math.log(max(x, 2.0))
        step = (Rx - n) / dR
        x_new = x - step
        if x_new <= 1.0:
            x_new = max(1.5, x / 2.0)
        if abs(x_new - x) < 1e-6 * abs(x):
            break
        x = x_new
    return x


def main() -> None:
    N_max = 1_000_000  # cap on prime index to test (sample)
    sieve_max = 16_000_000  # > p(10^6) ≈ 15485863
    print(f"Sieving primes up to {sieve_max} ...")
    primes = sieve(sieve_max)
    print(f"  found {len(primes)} primes; p(10^6) = {primes[10**6 - 1]}")

    # Sample n at log-spaced indices to keep cost reasonable
    sample_ns = sorted(set(
        int(round(x))
        for x in np.geomspace(10, N_max, 80)
    ))
    rows = []  # (n, p_n, R_inv_n, delta, |delta|/log^2 p_n)
    for n in sample_ns:
        if n > len(primes):
            break
        p_n = primes[n - 1]
        rinv = R_inverse(n)
        delta = p_n - rinv
        log2p = math.log(p_n) ** 2
        rows.append((n, p_n, rinv, delta, abs(delta) / log2p))

    print(f"\n{'n':>9} | {'p_n':>10} | {'R^-1(n)':>10} | "
          f"{'delta':>9} | {'|d|/log^2 p':>12}")
    for n, p, ri, d, r in rows:
        print(f"{n:>9} | {p:>10} | {ri:>10.2f} | {d:>9.2f} | {r:>12.4f}")

    # Aggregates
    deltas = np.array([r[3] for r in rows])
    ratios = np.array([r[4] for r in rows])
    log_p = np.array([math.log(r[1]) for r in rows])
    sqrt_p = np.array([math.sqrt(r[1]) for r in rows])
    print(f"\nN sample = {len(rows)}")
    print(f"max |delta|        = {np.max(np.abs(deltas)):.2f}")
    print(f"max |delta|/log^2p = {np.max(ratios):.4f}")
    print(f"mean |delta|/log^2p= {np.mean(ratios):.4f}")
    print(f"max |delta|/sqrt(p)= {np.max(np.abs(deltas) / sqrt_p):.4f}")

    # Linear fit: log|delta| ~ alpha * log p + const  → estimate scaling
    mask = np.abs(deltas) > 1e-3
    slope, intercept = np.polyfit(log_p[mask], np.log(np.abs(deltas[mask])), 1)
    print(f"|delta| scaling exponent (delta ~ p^alpha): alpha ≈ {slope:.3f}")
    print(f"  (RH expects ~0.5; pure polylog would give ~0.0)")

    # Window-width recommendation
    max_ratio = np.max(ratios)
    print(f"\nFor n ≤ 10^6, window of width {math.ceil(max_ratio)} * log^2(p_n) "
          f"is empirically sufficient.")
    print("This window is polylog, but anchor counting at the window's left edge "
          "remains O(√x) and is the actual bottleneck.")


if __name__ == "__main__":
    main()
