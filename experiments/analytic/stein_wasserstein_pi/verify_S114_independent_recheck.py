"""S114 verify: independent reproduction of S108's W_1 = 0.00829 at K=10000
using a fully different W_1 implementation, plus a focused falsification of
S113's "kurtosis fully predicts W_1/sigma" claim using a Beta(alpha,alpha)
distribution targeted at kurt = -0.41 (D_emp's kurtosis).

If S113's prediction holds, Beta(alpha,alpha) at kurt=-0.41 should give
W_1/sigma in [0.034, 0.046] (within ~10% of D_emp's 0.0376).
If Beta gives a value outside that band, S113's "kurtosis-only" claim is
too strong (but doesn't unsave A-grade — only makes S113 partial).

Independent of every existing project file: uses scipy emd via cdf
trapezoid integration AND scipy.stats.wasserstein_distance.
"""

import numpy as np
import sympy
from scipy import stats, special


# ---- D_emp computation (independent of stein_wasserstein_pi.py) ----

def D_emp(K, x_lo=1e6, x_hi=1e7):
    """Compute D̂_K = (pi(x) - Li(x)) log(x) / sqrt(x) on K log-uniform anchors
    in [x_lo, x_hi]. We round x to nearest int (sympy.primepi requires int)."""
    xs = np.geomspace(x_lo, x_hi, K)
    xs = np.unique(np.round(xs).astype(np.int64))
    pis = np.array([int(sympy.primepi(int(x))) for x in xs], dtype=np.int64)
    lis = special.expi(np.log(xs.astype(float)))
    D = (pis - lis) * np.log(xs.astype(float)) / np.sqrt(xs.astype(float))
    return xs, D


# ---- Two independent W_1 implementations ----

def W1_method_A_quantile(samples, mu, sigma):
    """Mid-rank quantile estimator (matches S113's protocol)."""
    K = len(samples)
    s = np.sort(samples)
    us = (np.arange(K) + 0.5) / K
    inv = stats.norm.ppf(us, loc=mu, scale=sigma)
    return float(np.mean(np.abs(s - inv)))


def W1_method_B_scipy(samples, mu, sigma, n_ref=200_000, seed=42):
    """Use scipy.stats.wasserstein_distance with a large-sample Gaussian
    reference (Monte Carlo Gaussian sample drawn at n_ref points)."""
    rng = np.random.default_rng(seed)
    ref = rng.normal(mu, sigma, size=n_ref)
    return float(stats.wasserstein_distance(samples, ref))


def W1_method_C_cdf_integral(samples, mu, sigma):
    """W_1(P, Q) = ∫ |F_P(x) - F_Q(x)| dx for 1D distributions.
    Compute this via fine grid integration."""
    s = np.sort(samples)
    K = len(s)
    # Grid covering range
    x_lo = float(min(s.min(), mu - 7 * sigma))
    x_hi = float(max(s.max(), mu + 7 * sigma))
    grid = np.linspace(x_lo, x_hi, 200_000)
    # F_emp on grid
    F_emp = np.searchsorted(s, grid, side='right') / K
    # F_gauss
    F_gauss = stats.norm.cdf(grid, loc=mu, scale=sigma)
    return float(np.trapezoid(np.abs(F_emp - F_gauss), grid))


# ---- run independent reproduction ----

def reproduce_D_emp():
    print("=" * 70)
    print("S114 verify-1: independent reproduction of S108 W_1 at K=10000")
    print("=" * 70)
    K = 10000
    print(f"Computing pi(x), Li(x) on K={K} anchors in [10^6, 10^7]...")
    xs, D = D_emp(K)
    K_eff = len(D)
    mu = float(np.mean(D))
    sigma = float(np.std(D))
    kurt = float(stats.kurtosis(D))
    print(f"  K_effective = {K_eff}")
    print(f"  mean(D)  = {mu:.6f}")
    print(f"  std(D)   = {sigma:.6f}")
    print(f"  kurt(D)  = {kurt:.4f}")

    W_A = W1_method_A_quantile(D, mu, sigma)
    W_B = W1_method_B_scipy(D, mu, sigma)
    W_C = W1_method_C_cdf_integral(D, mu, sigma)
    print(f"\n  W_1 (method A: mid-rank quantile)   = {W_A:.6f}  W_1/σ = {W_A/sigma:.5f}")
    print(f"  W_1 (method B: scipy MC Gaussian)   = {W_B:.6f}  W_1/σ = {W_B/sigma:.5f}")
    print(f"  W_1 (method C: CDF-integral)        = {W_C:.6f}  W_1/σ = {W_C/sigma:.5f}")
    print(f"\n  S108 reported (results_K10000_fixed.json): W_1 = 0.00829, W_1/σ = 0.03763")

    # Gaussian-control baseline z-score
    n_null = 200
    rng_null = np.random.default_rng(2026)
    Ws_null = np.array([
        W1_method_A_quantile(rng_null.normal(mu, sigma, K_eff), mu, sigma)
        for _ in range(n_null)
    ])
    z = (W_A - Ws_null.mean()) / Ws_null.std()
    print(f"\n  Gaussian-control W_1 mean = {Ws_null.mean():.6f}, std = {Ws_null.std():.6f}")
    print(f"  z-score (D vs sample-fitted Gaussian fluctuation) = {z:.2f}")
    return {"W_A": W_A, "W_B": W_B, "W_C": W_C, "sigma": sigma, "kurt": kurt, "z": float(z)}


# ---- S113 kurtosis-only claim falsification: Beta(α,α) at matched kurt ----

def beta_alpha_for_kurt(kurt_target):
    """For X ~ Beta(α,α) on (0,1), excess kurt of (X - 1/2) is
    kurt = -6 / (2α + 3).
    Solve for α at kurt_target."""
    alpha = (-6 / kurt_target - 3) / 2
    return alpha


def test_kurtosis_only_prediction():
    print("\n" + "=" * 70)
    print("S114 verify-2: falsify S113 'kurt fully predicts W_1/σ' via Beta")
    print("=" * 70)
    kurt_target = -0.41
    alpha = beta_alpha_for_kurt(kurt_target)
    print(f"  Building Beta({alpha:.3f}, {alpha:.3f}) -- predicted excess kurt = {kurt_target}")

    K = 10000
    n_trials = 30
    sigma_target = 1.6  # match D's std at X=10⁶
    seed = 2026
    rng = np.random.default_rng(seed)

    Ws = []
    kurts = []
    skews = []
    for t in range(n_trials):
        x = rng.beta(alpha, alpha, size=K) - 0.5
        x = (x - x.mean()) / x.std() * sigma_target
        mu_s = float(x.mean())
        sg_s = float(x.std())
        W = W1_method_A_quantile(x, mu_s, sg_s)
        Ws.append(W)
        kurts.append(float(stats.kurtosis(x)))
        skews.append(float(stats.skew(x)))
    W_mean = float(np.mean(Ws))
    W_std = float(np.std(Ws))
    kurt_mean = float(np.mean(kurts))
    skew_mean = float(np.mean(skews))
    W_over_sigma = W_mean / sigma_target

    print(f"\n  Beta(α,α) results, K={K}, n_trials={n_trials}:")
    print(f"    sample kurt mean = {kurt_mean:.4f} (target = {kurt_target})")
    print(f"    sample skew mean = {skew_mean:.4f}")
    print(f"    W_1 mean         = {W_mean:.6f} ± {W_std:.6f}")
    print(f"    W_1/σ            = {W_over_sigma:.5f}")

    print(f"\n  S113 predicted W_1/σ at kurt=-0.41: 0.0426")
    print(f"  D_emp observed W_1/σ at kurt=-0.41: 0.0376")
    print(f"  Beta(α,α) observed W_1/σ:          {W_over_sigma:.5f}")

    pred_band_lo = 0.0376 * 0.9
    pred_band_hi = 0.0376 * 1.1  # ±10% band around D_emp
    inside = pred_band_lo <= W_over_sigma <= pred_band_hi
    print(f"\n  Within ±10% band of D_emp [{pred_band_lo:.4f}, {pred_band_hi:.4f}]? {inside}")
    pred_band_lo2 = 0.0426 * 0.9
    pred_band_hi2 = 0.0426 * 1.1
    inside2 = pred_band_lo2 <= W_over_sigma <= pred_band_hi2
    print(f"  Within ±10% band of S113 prediction [{pred_band_lo2:.4f}, {pred_band_hi2:.4f}]? {inside2}")

    # also test a SECOND distribution at same kurt: scaled triangular?
    # Triangular (-a, 0, a): variance = a²/6, kurt = -0.6 (fixed). Not at -0.41.
    # Use a mixture instead: 0.5 N(0,σ_0²) + 0.5 N(±d, σ_d²) tuned to kurt=-0.41
    # Skip for time — Beta is the cleanest second-distribution check.

    return {"W_over_sigma_beta": W_over_sigma, "kurt_obtained": kurt_mean,
            "S113_pred": 0.0426, "D_emp_obs": 0.0376}


if __name__ == "__main__":
    r1 = reproduce_D_emp()
    r2 = test_kurtosis_only_prediction()
    print("\n" + "=" * 70)
    print("S114 SUMMARY")
    print("=" * 70)
    print(f"  W_A (mid-rank quantile)       = {r1['W_A']:.6f}  /σ = {r1['W_A']/r1['sigma']:.5f}")
    print(f"  W_C (CDF integral, indep)     = {r1['W_C']:.6f}  /σ = {r1['W_C']/r1['sigma']:.5f}")
    print(f"  S108 reported                 = 0.008290  /σ = 0.03763")
    diff_A = abs(r1['W_A'] - 0.00829) / 0.00829
    diff_C = abs(r1['W_C'] - 0.00829) / 0.00829
    print(f"  |W_A − S108|/S108 = {diff_A*100:.2f}%  ;  |W_C − S108|/S108 = {diff_C*100:.2f}%")
    print(f"  Beta(α,α) at kurt=-0.41 W_1/σ = {r2['W_over_sigma_beta']:.5f}")
    print(f"  Compared to D_emp = 0.0376 ; S113 prediction = 0.0426")
