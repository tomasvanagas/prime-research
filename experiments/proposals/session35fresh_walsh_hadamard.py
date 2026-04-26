"""
Proposal 1 — Sparse Walsh–Hadamard spectrum of chi_P.

Test: does chi_P[0..N) have effective WH-sparsity polylog(N)?
- N = 16384.
- Compute exact FWHT.
- Sort coefficients by magnitude; report decay.
- Reconstruct using top-k for several k; measure max |pi_hat(x) - pi(x)|.

VERDICT condition: top-k with k <= 32 reconstructs pi(x) within ±1 for all x.
"""
import numpy as np
from sympy import isprime


def build_prime_indicator(N):
    chi = np.zeros(N, dtype=np.float64)
    for k in range(2, N):
        if isprime(k):
            chi[k] = 1.0
    return chi


def fwht(a):
    """In-place Fast Walsh–Hadamard Transform (sequency / Hadamard ordering)."""
    a = a.copy()
    h = 1
    n = len(a)
    while h < n:
        for i in range(0, n, h * 2):
            for j in range(i, i + h):
                x = a[j]
                y = a[j + h]
                a[j] = x + y
                a[j + h] = x - y
        h *= 2
    return a


def main():
    L = 14
    N = 1 << L  # 16384
    print(f"# Walsh-Hadamard spectrum of chi_P, N = 2^{L} = {N}")
    chi = build_prime_indicator(N)
    pi_true = np.cumsum(chi)

    W = fwht(chi)
    mag = np.abs(W)

    # Sort coefficients by magnitude
    order = np.argsort(-mag)
    sorted_mag = mag[order]

    print("\n## Top-20 WH coefficient magnitudes")
    for i in range(20):
        idx = order[i]
        print(f"  rank {i:>2}  index {idx:>5}  |W|={sorted_mag[i]:.4f}")

    # Decay diagnostics
    print("\n## Cumulative L2 mass captured by top-k")
    total_l2 = float(np.sum(W ** 2))
    cum = np.cumsum(sorted_mag ** 2)
    for k in [1, 4, 16, 64, 256, 1024, 4096, N // 2, N - 1]:
        frac = cum[min(k, N) - 1] / total_l2 if total_l2 > 0 else 0.0
        print(f"  top-{k:>5}: {frac*100:.4f}% of L2 mass")

    # Reconstruction quality: keep top-k Walsh coefficients,
    # invert (FWHT is involutive up to scale 1/N), compare prefix sums.
    print("\n## Reconstruction error in pi(x)")
    print("  k    max|pi_hat - pi|   mean|pi_hat - pi|   any-exact-x?")
    for k in [16, 64, 256, 1024, 4096, N // 2, N]:
        Wk = np.zeros_like(W)
        keep_idx = order[:k]
        Wk[keep_idx] = W[keep_idx]
        chi_hat = fwht(Wk) / N  # inverse FWHT (Hadamard is self-inverse / N)
        chi_hat_round = np.round(chi_hat)
        pi_hat = np.cumsum(chi_hat_round)
        diff = pi_hat - pi_true
        max_err = float(np.max(np.abs(diff)))
        mean_err = float(np.mean(np.abs(diff)))
        # Also try without rounding, for raw L_inf
        pi_hat_raw = np.cumsum(chi_hat)
        diff_raw = pi_hat_raw - pi_true
        max_err_raw = float(np.max(np.abs(diff_raw)))
        print(f"  {k:>5}  max={max_err:>8.2f}  mean={mean_err:>8.4f}   "
              f"max_raw={max_err_raw:.4f}")

    # Entropy of normalized squared spectrum
    p = sorted_mag ** 2
    p = p / p.sum()
    nz = p[p > 0]
    H = float(-(nz * np.log2(nz)).sum())
    print(f"\n## Spectral Shannon entropy: H = {H:.3f} bits "
          f"(uniform on {N} bins would be {np.log2(N):.3f})")

    # Test polylog claim explicitly: how many coefficients to recover pi exactly?
    print("\n## Smallest k that gives pi_hat == pi_true everywhere")
    for k_try in [N // 8, N // 4, N // 2, int(N * 0.75), N - 1, N]:
        Wk = np.zeros_like(W)
        Wk[order[:k_try]] = W[order[:k_try]]
        chi_hat = np.round(fwht(Wk) / N)
        pi_hat = np.cumsum(chi_hat)
        ok = np.array_equal(pi_hat, pi_true)
        print(f"  k = {k_try:>5}  exact_match = {ok}")

    print("\n# DONE")


if __name__ == "__main__":
    main()
