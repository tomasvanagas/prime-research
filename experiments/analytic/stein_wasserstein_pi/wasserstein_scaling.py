"""
Diagnostic: how does W_1(D, fitted Gaussian) scale with K?

Stein-CLT prediction for genuinely-Gaussian sample: W_1 ~ 1/sqrt(K).
If D's W_1 plateaus at a positive constant c > 0, this is the
quantitative non-Gaussianity signal.

Run a sequence of K values; compare D's W_1 to bootstrapped
Gaussian-control W_1 at each K.
"""

import json
import time
import numpy as np
import sympy
from scipy import stats, special
from stein_wasserstein_pi import (
    primepi_array, li_array, D_signal, wasserstein1_to_normal
)


def main():
    Ks = [200, 500, 1000, 2000, 5000, 10000]
    rng = np.random.default_rng(101)

    # Use a SINGLE precomputed pi array to avoid recomputing.
    # Sample at the largest K, then subsample for smaller K.
    K_max = max(Ks)
    print(f"Pre-computing pi(x_k) at K_max = {K_max}")
    xs_full = np.unique(np.geomspace(1e6, 1e7, K_max).astype(np.int64))
    K_max = len(xs_full)
    t0 = time.time()
    pis_full = primepi_array(xs_full)
    print(f"  done in {time.time()-t0:.1f}s")
    lis_full = li_array(xs_full.astype(np.float64))
    D_full = D_signal(xs_full.astype(np.float64),
                      pis_full.astype(np.float64), lis_full)

    # Use the K_max-fitted (mu, sigma); scaling regime requires
    # comparing same population-fitted Gaussian, not K-dependent fit.
    mu_pop = float(np.mean(D_full))
    sigma_pop = float(np.std(D_full, ddof=0))
    print(f"Population (K={K_max}): mean={mu_pop:.6f} std={sigma_pop:.6f}")

    table = []
    for K in Ks:
        if K > K_max:
            continue
        # Subsample: take every (K_max // K) anchors
        step = K_max // K
        idx = np.arange(0, K_max, step)[:K]
        D_sub = D_full[idx]
        K_eff = len(D_sub)

        # W_1 to fitted Gaussian (use SAMPLE mean/std for fair comparison
        # to Gaussian-control with K samples)
        mu_s = float(np.mean(D_sub))
        sigma_s = float(np.std(D_sub, ddof=0))
        W1_D = wasserstein1_to_normal(D_sub, mu_s, sigma_s)

        # Gaussian control: each trial uses its OWN sample-fitted (mu, sigma)
        # for both sample and fitted Gaussian, mirroring D̂ procedure.
        n_ctrl = 50
        W1_gauss = np.empty(n_ctrl)
        for b in range(n_ctrl):
            g = rng.normal(mu_s, sigma_s, size=K_eff)
            m_g, s_g = float(np.mean(g)), float(np.std(g, ddof=0))
            W1_gauss[b] = wasserstein1_to_normal(g, m_g, s_g)
        gauss_mean = float(W1_gauss.mean())
        gauss_std = float(W1_gauss.std())
        z = (W1_D - gauss_mean) / gauss_std if gauss_std > 0 else 0.0

        # Skew/kurt at this subsample
        skew = float(stats.skew(D_sub))
        kurt = float(stats.kurtosis(D_sub))

        row = {
            "K": K_eff,
            "W1_D": W1_D,
            "W1_gauss_mean": gauss_mean,
            "W1_gauss_std": gauss_std,
            "z_score": z,
            "skewness": skew,
            "excess_kurtosis": kurt,
            "stein_clt_rate": 1.0 / np.sqrt(K_eff),
            "W1_to_rate_ratio": W1_D * np.sqrt(K_eff),
        }
        table.append(row)
        print(f"K={K_eff:6d}: W_1(D)={W1_D:.6f}  gauss_ctrl="
              f"{gauss_mean:.6f}+-{gauss_std:.6f}  z={z:+6.2f}  "
              f"W_1*sqrt(K)={row['W1_to_rate_ratio']:.4f}  "
              f"skew={skew:+.3f} kurt={kurt:+.3f}", flush=True)
        # Save partial results in case of interruption
        with open("scaling_results.json", "w") as f:
            json.dump({
                "K_max": K_max,
                "mu_pop": mu_pop,
                "sigma_pop": sigma_pop,
                "table": table,
            }, f, indent=2)

    with open("scaling_results.json", "w") as f:
        json.dump({
            "K_max": K_max,
            "mu_pop": mu_pop,
            "sigma_pop": sigma_pop,
            "table": table,
        }, f, indent=2)
    print(f"\nWrote scaling_results.json")
    print()
    print("Interpretation:")
    print("  W_1(D)*sqrt(K) → constant  =>  W_1 ~ 1/sqrt(K) (Gaussian rate)")
    print("  W_1(D)*sqrt(K) → ∞         =>  W_1 plateaus, NON-GAUSSIAN")
    print("  W_1 / Gauss_ctrl → 1       =>  Gaussian-indistinguishable")


if __name__ == "__main__":
    main()
