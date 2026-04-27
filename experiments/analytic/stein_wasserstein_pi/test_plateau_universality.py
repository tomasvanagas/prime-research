"""S113 verification: is the c(X) ≈ 0.0083 plateau a generic feature of
ANY non-Gaussian distribution sampled at K log-uniform anchors and
compared (W_1) to the sample-fitted Gaussian?

S112 showed the plateau is consistent with random-phase oscillatory sums
of arbitrary (non-Riemann) frequencies in [10, 145]. This goes one step
further: do non-arithmetic, non-oscillatory distributions also produce a
plateau at ~0.04σ when matched to D's variance and (approximately)
kurtosis?

If yes: the W_1 plateau is just a generic non-Gaussianity quantity, with
magnitude determined by sample distribution shape (kurtosis primarily) +
sample size, NOT by anything Riemann/arithmetic. The "novel finite-x
Wasserstein bound" then reduces to "kurt(D) ≈ -0.41" — a measurement,
not an edge.

Independent of S108/S109/S110/S111/S112 code; uses only numpy + scipy.
"""

import numpy as np
from scipy import stats


def W1_to_normal(samples, mu, sigma):
    """W_1 between empirical distribution and N(mu, sigma^2).
    Uses K equispaced quantile points - mid-rank empirical."""
    K = len(samples)
    s = np.sort(samples)
    us = (np.arange(K) + 0.5) / K
    inv = stats.norm.ppf(us, loc=mu, scale=sigma)
    return float(np.mean(np.abs(s - inv)))


def run_protocol(sample_fn, K, n_trials, sigma_target, rng):
    """Apply same protocol as S108: K samples, fit Gaussian, compute W_1."""
    Ws = []
    kurts = []
    for _ in range(n_trials):
        x = sample_fn(K, rng)
        # standardize to unit variance, then scale to sigma_target
        x = (x - np.mean(x)) / np.std(x) * sigma_target
        mu_s = float(np.mean(x))
        sig_s = float(np.std(x))
        W = W1_to_normal(x, mu_s, sig_s)
        Ws.append(W)
        kurts.append(float(stats.kurtosis(x)))
    return {
        "W1_mean": float(np.mean(Ws)),
        "W1_std": float(np.std(Ws)),
        "kurt_mean": float(np.mean(kurts)),
    }


# ---- distribution generators ----

def gaussian(K, rng):
    return rng.normal(0, 1, size=K)


def uniform_sym(K, rng):
    """U[-sqrt(3), sqrt(3)] — variance 1, kurt = -1.2."""
    return rng.uniform(-np.sqrt(3), np.sqrt(3), size=K)


def arcsine_one_mode(K, rng):
    """cos(2π U) — single arcsine, kurt = -1.5."""
    return np.cos(2 * np.pi * rng.uniform(0, 1, size=K))


def arcsine_sum_50_unit_weights(K, rng):
    """Sum of 50 i.i.d. cos(U_k) with unit weights — central-limit shape."""
    n_modes = 50
    phases = rng.uniform(0, 2 * np.pi, size=(K, n_modes))
    return np.sum(np.cos(phases), axis=1)


def arcsine_sum_50_powerlaw_weights(K, rng):
    """Sum of 50 cos(U_k) with weights 1/k^0.5 — frequency-decay shape
    similar to explicit-formula 2/sqrt(¼+γ²) low-zero weights."""
    n_modes = 50
    weights = 1.0 / np.sqrt(np.arange(1, n_modes + 1))
    phases = rng.uniform(0, 2 * np.pi, size=(K, n_modes))
    return np.sum(weights * np.cos(phases), axis=1)


def two_gauss_mixture(K, rng):
    """0.5 N(-1, 0.3²) + 0.5 N(1, 0.3²) — bimodal, kurt ~ -1.78."""
    z = rng.binomial(1, 0.5, size=K)
    return rng.normal(2 * z - 1, 0.3, size=K)


def t_distribution_df10(K, rng):
    """t with df=10 — kurt = +1.0 (super-Gaussian)."""
    return rng.standard_t(df=10, size=K)


def laplace(K, rng):
    """Laplace (double exponential) — kurt = +3.0."""
    return rng.laplace(0, 1, size=K)


def low_zero_sum_log_uniform(K, rng):
    """Riemann low-zero explicit-formula sum, evaluated on
    K log-uniform anchors in [10⁶, 10⁷]. This is the analytic
    model of D̂."""
    GAMMAS = np.array([
        14.134725141734693, 21.022039638771555, 25.010857580145688,
        30.424876125859513, 32.935061587739189, 37.586178158825671,
        40.918719012147495, 43.327073280914999, 48.005150881167160,
        49.773832477672302, 52.970321477714460, 56.446247697063394,
        59.347044002602353, 60.831778524609810, 65.112544048081606,
        67.079810529494173, 69.546401711173979, 72.067157674481907,
        75.704690699083933, 77.144840068874805,
    ])
    xs = np.geomspace(1e6, 1e7, K)
    L = np.log(xs)
    weights = 2.0 / np.sqrt(0.25 + GAMMAS ** 2)
    phis = np.arctan(2 * GAMMAS)
    contribs = np.cos(np.outer(L, GAMMAS) - phis)  # (K, n_zeros)
    return -np.sum(weights * contribs, axis=1)  # match sign convention


# ---- experiment ----

def main():
    seed = 2026
    rng = np.random.default_rng(seed)
    sigma_target = 1.6  # approximate D's std at X=10⁶

    K_values = [200, 500, 1000, 2000, 5000, 10000]
    n_trials = 30

    distributions = {
        "Gaussian (control)": gaussian,
        "Uniform symmetric": uniform_sym,
        "Single arcsine cos(U)": arcsine_one_mode,
        "Sum 50 arcsines unit weights": arcsine_sum_50_unit_weights,
        "Sum 50 arcsines 1/sqrt(k) weights": arcsine_sum_50_powerlaw_weights,
        "Two-Gaussian mixture (bimodal)": two_gauss_mixture,
        "t-distribution df=10": t_distribution_df10,
        "Laplace": laplace,
        "Low-zero sum log-uniform (analytic D)": low_zero_sum_log_uniform,
    }

    results = {}
    for name, fn in distributions.items():
        print(f"\n=== {name} ===")
        # Use a fresh rng per distribution for clean comparison
        sub_rng = np.random.default_rng(seed)
        results[name] = {}
        for K in K_values:
            r = run_protocol(fn, K, n_trials, sigma_target, sub_rng)
            results[name][K] = r
            print(f"  K={K:5d}  W_1 = {r['W1_mean']:.5f} ± {r['W1_std']:.5f}  "
                  f"kurt = {r['kurt_mean']:+.3f}  W_1/σ = {r['W1_mean']/sigma_target:.5f}")

    # Summary: does W_1 PLATEAU (saturate) or DECAY (~ 1/sqrt(K))?
    print("\n\n=== Plateau diagnostic ===")
    print("Compare W_1 at K=10000 vs the 1/sqrt(K)-extrapolated value:")
    print(f"{'Distribution':50s}  {'K=200':>8s}  {'K=10000':>8s}  "
          f"{'1/sqrt(K) ratio':>15s}  {'Plateau?':>10s}")
    for name in distributions:
        W_low = results[name][200]["W1_mean"]
        W_high = results[name][10000]["W1_mean"]
        # If pure 1/sqrt(K), ratio should be sqrt(10000/200) = sqrt(50) ≈ 7.07
        ratio = W_low / W_high if W_high > 0 else float('inf')
        plateau = ratio < 3.0  # plateau if W_1 didn't decay as much as 1/sqrt(K)
        print(f"{name:50s}  {W_low:.5f}  {W_high:.5f}  {ratio:>15.2f}  {'YES' if plateau else ' no'}")
    # 1/sqrt(K) prediction: ratio = 7.07. Significant plateau = ratio << 7.07.

    # save
    import json
    with open("/tmp/verify_S113_universality_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print("\nWrote /tmp/verify_S113_universality_results.json")


if __name__ == "__main__":
    main()
