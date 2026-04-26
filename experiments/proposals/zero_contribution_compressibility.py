"""
Compressibility of the zeta-zero contribution matrix.

Build M[i, j] = 2 Re(Li(x_i^{rho_j})) for x_i in window, rho_j = 1/2 + i*t_j.
Test:
  1. Singular value decay (numerical low-rank).
  2. Sparsity in 2D Fourier basis (number of large coefficients).
  3. Sparsity in 2D Haar wavelet basis.

A polylog-recoverable signal needs polylog effective rank or polylog
coefficient count.
"""
import numpy as np
from mpmath import mp, mpc, li, mpf
from pathlib import Path


def load_zeros(path, K):
    with open(path) as f:
        ts = [float(line.strip()) for line in f.readlines()[:K]]
    return np.array(ts)


def build_contribution_matrix(x_lo, x_hi, ts):
    """M[i, j] = 2 Re(li(x_i^rho_j)) approximated via numpy complex pow + mpmath li.
       For speed: use the asymptotic approximation li(x^rho) ~ x^rho / (rho * log x)."""
    xs = np.arange(x_lo, x_hi, dtype=np.float64)
    log_xs = np.log(xs)
    M = np.zeros((len(xs), len(ts)), dtype=np.float64)
    for j, t in enumerate(ts):
        # rho = 1/2 + i*t
        # x^rho = sqrt(x) * exp(i * t * log x)
        # contribution to pi(x): -2 Re(li(x^rho))
        # asymptotic: li(x^rho) ~ x^rho / (rho * log x)
        # 2 Re(li(x^rho)) ~ 2 sqrt(x) * Re(exp(i*t*log x) / ((1/2+it) * log x))
        denom = (0.5 + 1j * t)
        contrib = 2 * np.sqrt(xs) * np.real(np.exp(1j * t * log_xs) / (denom * log_xs))
        M[:, j] = contrib
    return M, xs


def svd_decay(M, label):
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    S_norm = S / S[0]
    print(f"\n{label}")
    print(f"  shape={M.shape}, ||M||_F={np.linalg.norm(M):.3e}")
    print(f"  singular values (normalized) at k = 1, 2, 4, 8, …, 256:")
    ks = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
    for k in ks:
        if k <= len(S):
            print(f"    sigma_{k:4d}/sigma_0 = {S_norm[k-1]:.3e}")
    # Effective rank: number of sigmas above 1% relative
    eff_1pct = np.sum(S_norm > 0.01)
    eff_01pct = np.sum(S_norm > 0.001)
    print(f"  effective rank @1%   = {eff_1pct}")
    print(f"  effective rank @0.1% = {eff_01pct}")
    return S


def fourier_sparsity(M, label):
    F = np.fft.fft2(M)
    F_norm = np.abs(F) / np.abs(F).max()
    nonzero = np.sum(F_norm > 0.01)
    nonzero_strict = np.sum(F_norm > 0.001)
    print(f"  Fourier 2D coefficients > 1%   : {nonzero} / {M.size} = {nonzero/M.size*100:.2f}%")
    print(f"  Fourier 2D coefficients > 0.1% : {nonzero_strict} / {M.size} = {nonzero_strict/M.size*100:.2f}%")


def main():
    zeros_path = Path("data/zeta_zeros_8000.txt")
    K_max = 4000
    ts = load_zeros(zeros_path, K_max)

    # Window 1: small x for direct comparison
    print("=== Window x in [1000, 2000), K = 1000 zeros ===")
    M1, _ = build_contribution_matrix(1000, 2000, ts[:1000])
    svd_decay(M1, "M1 (1000 x 1000)")
    fourier_sparsity(M1, "M1")

    # Window 2: larger
    print("\n=== Window x in [5000, 6000), K = 2000 zeros ===")
    M2, _ = build_contribution_matrix(5000, 6000, ts[:2000])
    svd_decay(M2, "M2 (1000 x 2000)")
    fourier_sparsity(M2, "M2")

    # Comparison: random Gaussian matrix
    print("\n=== Random Gaussian baseline (1000 x 1000) ===")
    rng = np.random.default_rng(0)
    R = rng.standard_normal((1000, 1000))
    svd_decay(R, "Gaussian baseline")

    # Comparison: low-rank-plus-noise (truth for what we'd need)
    print("\n=== Synthetic rank-10 + small noise ===")
    A = rng.standard_normal((1000, 10))
    B = rng.standard_normal((10, 1000))
    LR = A @ B + 0.001 * rng.standard_normal((1000, 1000))
    svd_decay(LR, "rank-10 + noise")


if __name__ == "__main__":
    main()
