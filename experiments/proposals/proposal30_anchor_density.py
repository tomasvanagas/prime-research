"""Proposal 30 — Cancellation-anchor density probe.

Question: how often does the truncated explicit-formula error
   E_K(y) = sum_{|gamma_k| < gamma_K-th, k <= K} 2 Re(y^rho / rho)
fall to polylog smallness inside a fixed integer window?

If anchors with |E_K(y)| <= sqrt(log y) have density bounded below in
windows of width polylog, we could navigate to one and use it. Predicted
GUE density: probability ≈ polylog/√y per integer y → 0.

Approach
--------
- Load Riemann zeros.
- For x in {10^3, 10^4, 10^5}, K in {10, 50, 200}:
  evaluate E_K(y) for integer y in [x, x + 10^4].
  Bucket |E_K| / sqrt(y); compute fraction with |E_K(y)| <= sqrt(log y),
  fraction with |E_K(y)| <= log^2(y), and the empirical CDF.
- Compare to GUE predicted density.
"""
import math
from pathlib import Path

import numpy as np


REPO = Path(__file__).resolve().parents[2]
ZERO_FILE = REPO / "data" / "zeta_zeros_1000.txt"


def load_zeros(max_zeros: int = 1000) -> np.ndarray:
    zs = []
    with open(ZERO_FILE) as fh:
        for line in fh:
            line = line.strip()
            if line:
                try:
                    zs.append(float(line.split()[0]))
                except ValueError:
                    pass
            if len(zs) >= max_zeros:
                break
    return np.array(zs)


def E_K(y: float, gammas: np.ndarray) -> float:
    """sum_{k=1..K} 2 Re(y^rho / rho), rho = 1/2 + i gamma_k."""
    sqrt_y = math.sqrt(y)
    log_y = math.log(y)
    s = 0.0
    for g in gammas:
        rho_real = 0.5
        rho_imag = g
        # y^rho = sqrt(y) * (cos(g log y) + i sin(g log y))
        c = math.cos(g * log_y)
        si = math.sin(g * log_y)
        # 1/rho = (1/2 - i g) / (1/4 + g^2)
        denom = 0.25 + g * g
        inv_re = 0.5 / denom
        inv_im = -g / denom
        # 2 Re(y^rho / rho) = 2 sqrt(y) (c * inv_re - si * inv_im)
        s += 2 * sqrt_y * (c * inv_re - si * inv_im)
    return s


def main() -> None:
    gammas = load_zeros(max_zeros=1000)
    print(f"Loaded {len(gammas)} zeros.")

    for x in [1000, 10000, 100000]:
        log_x = math.log(x)
        sqrt_x = math.sqrt(x)
        for K in [10, 50, 200]:
            ys = np.arange(x, x + 10_001, dtype=np.int64)
            es = np.empty(len(ys), dtype=float)
            for i, y in enumerate(ys):
                es[i] = E_K(float(y), gammas[:K])
            abs_e = np.abs(es)
            # thresholds
            t1 = math.sqrt(log_x)            # polylog: very tight
            t2 = log_x ** 2                  # Cramér-scale
            t3 = sqrt_x / log_x              # standard explicit-formula bound
            f1 = float(np.mean(abs_e <= t1))
            f2 = float(np.mean(abs_e <= t2))
            f3 = float(np.mean(abs_e <= t3))
            print(f"x={x}, K={K}: |E_K| stats over 10^4 ys: "
                  f"max={abs_e.max():.2f}, mean={abs_e.mean():.2f}, "
                  f"frac<=sqrt(log x)={f1:.4f} (thr={t1:.3f}), "
                  f"frac<=log^2 x={f2:.4f} (thr={t2:.3f}), "
                  f"frac<=sqrt(x)/log x={f3:.4f} (thr={t3:.3f})")

    # Expectation: as x grows, the fraction with |E_K| <= sqrt(log x) shrinks
    # toward zero, confirming exponentially sparse anchors. The fraction with
    # |E_K| <= sqrt(x)/log x stays near 1.0, since that is the typical scale.


if __name__ == "__main__":
    main()
