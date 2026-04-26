"""Proposal 27 — Hermite-mollified reverse explicit formula.

Question: does Gaussian/Hermite-mollification of the truncated explicit formula
psi(x) ≈ x - sum_{|gamma|<T} x^rho/rho
converge faster (in number K of zeros used) than the unmollified truncation,
for x in {100, 1000, 10000}?

Approach
--------
- Load Riemann zeros gamma_n.
- For target x, compute three estimators of psi(x):
    (a) unmollified truncated sum,
    (b) Gaussian-mollified truncated sum (kernel exp{-(gamma - log x)^2 / 2*sigma^2}),
    (c) Riesz-mean (1 - |gamma|/T)^k truncated sum, as a baseline mollification.
- Compare |estimator - psi_exact(x)| as a function of K = number of zeros used.

If (b) decays faster than (a) by more than a constant factor, mollification
helps. Otherwise the square-root-cancellation barrier dominates.

Run on small N. No iterative tuning - one shot.
"""
import math
import os
import sys
from pathlib import Path

import numpy as np


REPO = Path(__file__).resolve().parents[2]
ZERO_FILE_CANDIDATES = [
    REPO / "data" / "zeta_zeros_2000.txt",
    REPO / "data" / "zeta_zeros_1000.txt",
    REPO / "data" / "zeta_zeros_500.txt",
    REPO / "data" / "zeta_zeros_300.txt",
    REPO / "data" / "zeta_zeros_200.txt",
]


def load_zeros(max_zeros: int = 1000) -> np.ndarray:
    for f in ZERO_FILE_CANDIDATES:
        if f.exists():
            zs = []
            with open(f) as fh:
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        zs.append(float(line.split()[0]))
                    except ValueError:
                        pass
                    if len(zs) >= max_zeros:
                        break
            if zs:
                return np.array(zs)
    raise FileNotFoundError("No Riemann-zeros data file found in data/.")


def psi_exact(x: int) -> float:
    """psi(x) = sum_{p^k <= x} log p, computed via sieve up to x."""
    if x < 2:
        return 0.0
    sieve = bytearray([1]) * (x + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(math.isqrt(x)) + 1):
        if sieve[i]:
            for j in range(i * i, x + 1, i):
                sieve[j] = 0
    primes = [i for i, v in enumerate(sieve) if v]
    s = 0.0
    for p in primes:
        pk = p
        while pk <= x:
            s += math.log(p)
            pk *= p
    return s


def estimator_a_unmollified(x: float, gammas: np.ndarray) -> np.ndarray:
    """psi_K(x) for K = 1..len(gammas). Returns cumulative array."""
    log_x = math.log(x)
    sqrt_x = math.sqrt(x)
    contribs = np.empty(len(gammas))
    for i, g in enumerate(gammas):
        rho = 0.5 + 1j * g
        term = (x ** rho) / rho
        contribs[i] = -2.0 * term.real
    cum = np.cumsum(contribs)
    main = x - math.log(2 * math.pi) - 0.5 * math.log(max(1.0 - x ** -2, 1e-300))
    return main + cum


def estimator_b_gaussian(x: float, gammas: np.ndarray, sigma: float) -> np.ndarray:
    """Gaussian-mollified: weight each zero by exp(-(gamma - log x)^2 / (2 sigma^2))."""
    log_x = math.log(x)
    contribs = np.empty(len(gammas))
    for i, g in enumerate(gammas):
        rho = 0.5 + 1j * g
        kernel = math.exp(-((g - log_x) ** 2) / (2 * sigma * sigma))
        term = (x ** rho) / rho
        contribs[i] = -2.0 * term.real * kernel
    cum = np.cumsum(contribs)
    main = x - math.log(2 * math.pi) - 0.5 * math.log(max(1.0 - x ** -2, 1e-300))
    return main + cum


def estimator_c_riesz(x: float, gammas: np.ndarray, T: float, k: int = 2) -> np.ndarray:
    """Riesz-mean truncation: weight (1 - |gamma|/T)^k for |gamma| < T."""
    contribs = np.empty(len(gammas))
    for i, g in enumerate(gammas):
        if abs(g) >= T:
            contribs[i] = 0.0
            continue
        w = (1.0 - abs(g) / T) ** k
        rho = 0.5 + 1j * g
        term = (x ** rho) / rho
        contribs[i] = -2.0 * term.real * w
    cum = np.cumsum(contribs)
    main = x - math.log(2 * math.pi) - 0.5 * math.log(max(1.0 - x ** -2, 1e-300))
    return main + cum


def main() -> None:
    gammas = load_zeros(max_zeros=1000)
    print(f"Loaded {len(gammas)} zeros; first three gammas: {gammas[:3]}")

    summary_rows = []
    for x in [100, 1000, 10000]:
        psi_x = psi_exact(int(x))
        log_x = math.log(x)
        sqrt_x = math.sqrt(x)

        ks = [10, 25, 50, 100, 200, 400, 800]
        ks = [k for k in ks if k <= len(gammas)]

        # Unmollified
        cum_a = estimator_a_unmollified(x, gammas)
        # Gaussian-mollified at three sigma values
        sigmas = [0.5, 1.0, 2.0]
        gaussian_cums = {s: estimator_b_gaussian(x, gammas, s) for s in sigmas}
        # Riesz-mean at three T values
        Ts = [50.0, 200.0, 800.0]
        riesz_cums = {T: estimator_c_riesz(x, gammas, T, k=2) for T in Ts}

        print(f"\n=== x = {x},  psi(x) = {psi_x:.4f},  sqrt(x)/log(x) = {sqrt_x/log_x:.3f}")
        print(f"  K  | unmoll  |  G(0.5)  |  G(1)   |  G(2)   |  R(T=50) | R(T=200)| R(T=800)")
        for K in ks:
            row = [
                abs(cum_a[K - 1] - psi_x),
                abs(gaussian_cums[0.5][K - 1] - psi_x),
                abs(gaussian_cums[1.0][K - 1] - psi_x),
                abs(gaussian_cums[2.0][K - 1] - psi_x),
                abs(riesz_cums[50.0][K - 1] - psi_x),
                abs(riesz_cums[200.0][K - 1] - psi_x),
                abs(riesz_cums[800.0][K - 1] - psi_x),
            ]
            print(f"  {K:>3} | "
                  f"{row[0]:7.3f} |  {row[1]:7.3f} | {row[2]:7.3f} | {row[3]:7.3f} | "
                  f"{row[4]:7.3f} | {row[5]:7.3f} | {row[6]:7.3f}")
            summary_rows.append((x, K, *row))

    # Overall verdict: at fixed K, does any estimator beat unmollified by > 2x?
    # We measure ratio at maximum K used per x.
    print("\nVerdict probe: ratio (unmoll / best_mollified) at K = K_max per x")
    by_x = {}
    for x, K, a, *rest in summary_rows:
        if x not in by_x or K > by_x[x][0]:
            by_x[x] = (K, a, *rest)
    for x, vals in by_x.items():
        K, a, g05, g1, g2, r50, r200, r800 = vals
        best_moll = min(g05, g1, g2, r50, r200, r800)
        ratio = a / best_moll if best_moll > 0 else float('inf')
        print(f"  x={x}, K={K}: unmoll={a:.3f}, best_moll={best_moll:.3f}, "
              f"ratio={ratio:.2f}")


if __name__ == "__main__":
    main()
