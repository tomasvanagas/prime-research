"""
Stein's-method quantitative finite-x Gaussianity test for
D(x) = (pi(x) - Li(x)) * log(x) / sqrt(x).

Question (ATTACK_VECTORS.md §C5): does the empirical Wasserstein-1
distance W_1(D_hat_X, N(0, sigma_hat^2)) scale as the Stein-CLT
prediction O(1/sqrt(K)) (genuine Gaussianity), or plateau at a
positive constant (non-Gaussian arithmetic structure)?

If it plateaus quantitatively above the Gaussian-control baseline,
this is the FIRST quantitative non-Gaussianity result for
pi(x) - Li(x) at finite x.

Cross-domain technique: Stein's method (Chen-Goldstein-Shao 2011),
exchangeable-pair Wasserstein bound (Ross 2011, Probability Surveys).
"""

from __future__ import annotations
import argparse
import json
import time
import numpy as np
import sympy
from scipy import stats, special


# --------- core arithmetic ---------

def primepi_array(xs: np.ndarray) -> np.ndarray:
    """Compute pi(x) for each x in xs via sympy.primepi."""
    return np.array([int(sympy.primepi(int(x))) for x in xs], dtype=np.int64)


def li(x: float) -> float:
    """Logarithmic integral Li(x) = real(Ei(log x))."""
    return float(special.expi(np.log(x)))


def li_array(xs: np.ndarray) -> np.ndarray:
    return special.expi(np.log(xs))


def D_signal(xs: np.ndarray, pis: np.ndarray, lis: np.ndarray) -> np.ndarray:
    """Renormalised pi(x) - Li(x) deviation:
        D(x) = (pi(x) - Li(x)) * log(x) / sqrt(x)
    Under RH, |D(x)| = O(log x); D(x) -> centered weak limit."""
    return (pis - lis) * np.log(xs) / np.sqrt(xs)


# --------- W_1 (Wasserstein-1) ---------

def wasserstein1_to_normal(samples: np.ndarray,
                           mu: float | None = None,
                           sigma: float | None = None) -> float:
    """
    Compute W_1 between empirical distribution of `samples` and
    N(mu, sigma^2), using sorted-CDF integration in closed form.

    For sorted x_(1) <= ... <= x_(K),
        W_1 = sum_k integral_{(k-1)/K}^{k/K} |x_(k) - F^{-1}(u)| du
    where F is the Gaussian CDF.  We use the substitution u = F(t)
    and split each [x_(k-?), F^{-1}(k/K)] piece accordingly.
    """
    K = len(samples)
    if mu is None:
        mu = float(np.mean(samples))
    if sigma is None:
        sigma = float(np.std(samples, ddof=0))
    if sigma == 0:
        sigma = 1e-30
    s = np.sort(samples)
    # quantile boundaries u_k = k/K, k = 0..K
    us = np.arange(K + 1) / K
    # F^{-1}(u) for u in (0, 1); endpoints are infinite -> use eps clamp
    eps = 1e-12
    us_clipped = np.clip(us, eps, 1 - eps)
    inv = stats.norm.ppf(us_clipped, loc=mu, scale=sigma)
    inv[0] = mu - 6 * sigma  # numerical infinity sentinel
    inv[-1] = mu + 6 * sigma
    # On interval (u_{k-1}, u_k), empirical inverse is s[k-1] (constant).
    # |s[k-1] - F^{-1}(u)| du.  Use trapezoidal w/ 4 subdivisions per slab.
    M = 8
    total = 0.0
    for k in range(K):
        u_lo, u_hi = us_clipped[k], us_clipped[k + 1]
        # subdivide
        sub_us = np.linspace(u_lo, u_hi, M + 1)
        sub_inv = stats.norm.ppf(sub_us, loc=mu, scale=sigma)
        integrand = np.abs(s[k] - sub_inv)
        # trapezoidal in u-space
        total += np.trapezoid(integrand, sub_us)
    return float(total)


def wasserstein1_empirical_pair(a: np.ndarray, b: np.ndarray) -> float:
    """W_1 between two empirical distributions of equal size = mean
    of |sorted(a) - sorted(b)|."""
    if len(a) != len(b):
        # Resample / interpolate for unequal sizes
        K = min(len(a), len(b))
        a_s = np.sort(a)
        b_s = np.sort(b)
        ua = (np.arange(len(a_s)) + 0.5) / len(a_s)
        ub = (np.arange(len(b_s)) + 0.5) / len(b_s)
        u = np.linspace(0.5 / K, 1 - 0.5 / K, K)
        return float(np.mean(np.abs(np.interp(u, ua, a_s) - np.interp(u, ub, b_s))))
    return float(np.mean(np.abs(np.sort(a) - np.sort(b))))


# --------- Stein bounds (exchangeable-pair / linear) ---------

def stein_skew_kurt_bounds(samples: np.ndarray) -> dict:
    """Skewness/kurtosis-based "soft" Stein-CLT bounds.
    Returns dict with values + 95% bootstrap CIs."""
    z = (samples - np.mean(samples)) / np.std(samples)
    skew = float(stats.skew(z))
    kurt_excess = float(stats.kurtosis(z))  # excess (Fisher)
    # Bootstrap CI
    rng = np.random.default_rng(7)
    B = 2000
    skews = np.empty(B)
    kurts = np.empty(B)
    for i in range(B):
        bs = rng.choice(z, size=len(z), replace=True)
        skews[i] = stats.skew(bs)
        kurts[i] = stats.kurtosis(bs)
    return {
        "skewness": skew,
        "kurtosis_excess": kurt_excess,
        "skew_CI95": [float(np.quantile(skews, 0.025)),
                      float(np.quantile(skews, 0.975))],
        "kurt_CI95": [float(np.quantile(kurts, 0.025)),
                      float(np.quantile(kurts, 0.975))],
    }


def stein_exchangeable_pair_bound(samples: np.ndarray,
                                  rng: np.random.Generator) -> dict:
    """
    Construct an exchangeable pair (W, W') and compute the
    Stein-Wasserstein bound (Chatterjee-Shao / Ross 2011 Eq 4.1):

        W_1(W, Z) <= sqrt(Var(E[W' - W | W]))/lambda
                     + (1/(2 lambda)) * E|(W'-W)^2 - 2 lambda Var(W)|
                     + (1/(2 lambda)) * E|W' - W|^3

    Standard exchangeable pair on a finite sample: pick i uniformly,
    set W = w_i; pick j uniformly distinct from i, set W' = w_j.

    Then E[W' - W | W = w_i] = (1/(K-1)) sum_{j != i} (w_j - w_i)
                              = (K * mean - w_i) / (K-1) - w_i if mean=0
                              -> -lambda * w_i with lambda = K/(K-1)
    after centring; we centre then standardise to unit variance so
    that lambda is unambiguous.
    """
    z = (samples - np.mean(samples)) / np.std(samples, ddof=0)
    K = len(z)
    # E[W'-W | W = z_i] for centred z is -K/(K-1) * z_i
    # so the regression coefficient lambda = K/(K-1) ~ 1.
    lam = K / (K - 1)

    # Compute moments of (W' - W) by exact enumeration if K small,
    # else Monte Carlo.
    n_mc = 200_000
    i = rng.integers(0, K, size=n_mc)
    j = rng.integers(0, K - 1, size=n_mc)
    j[j >= i] += 1  # ensure j != i
    diff = z[j] - z[i]

    E_diff_sq = float(np.mean(diff ** 2))
    E_abs_diff_cu = float(np.mean(np.abs(diff) ** 3))

    # Conditional E[(W'-W)^2 | W = z_i] for centred z, var=1:
    # = (1/(K-1)) sum_{j != i} (z_j - z_i)^2
    cond_sq_sum = (np.sum(z ** 2) + K * z ** 2 - 2 * z * np.sum(z)) / (K - 1)
    cond_sq = cond_sq_sum  # array of size K
    # 2 lambda Var(W) = 2 lambda  (W standardised => Var=1)
    target = 2 * lam
    bias_term = float(np.mean(np.abs(cond_sq - target)))
    var_cond = float(np.mean((cond_sq - target) ** 2))

    bound_term1 = np.sqrt(var_cond) / lam  # variance of regression residual
    bound_term2 = bias_term / (2 * lam)
    bound_term3 = E_abs_diff_cu / (2 * lam)

    bound_total = bound_term1 + bound_term2 + bound_term3
    return {
        "lambda": lam,
        "E_diff_sq": E_diff_sq,
        "E_abs_diff_cu": E_abs_diff_cu,
        "stein_bound_term1_residual_var": bound_term1,
        "stein_bound_term2_bias": bound_term2,
        "stein_bound_term3_third_moment": bound_term3,
        "stein_bound_total": bound_total,
    }


# --------- decomposition / structural diagnostics ---------

def low_zero_modes(xs: np.ndarray, n_zeros: int = 10) -> np.ndarray:
    """
    Riemann's explicit formula's leading oscillatory contribution:
        D(x) = -(2 log x / sqrt(x)) sum_{rho} Li(x^rho)
             ~ -2 sum_rho cos(gamma_rho log x) / |rho|
    after the standard partial-fraction approximation.

    Build the matrix of low-zero contributions [cos(gamma_k log x_i),
    sin(gamma_k log x_i)] and return basis matrix.  The first n_zeros
    Riemann zeros (imag parts) are hard-coded -- standard reference
    values from Odlyzko's tables.
    """
    GAMMAS = np.array([
        14.134725141734693, 21.022039638771555, 25.010857580145688,
        30.424876125859513, 32.935061587739189, 37.586178158825671,
        40.918719012147495, 43.327073280914999, 48.005150881167160,
        49.773832477672302, 52.970321477714460, 56.446247697063394,
        59.347044002602353, 60.831778524609810, 65.112544048081606,
        67.079810529494173, 69.546401711173979, 72.067157674481907,
        75.704690699083933, 77.144840068874805,
    ])
    g = GAMMAS[:n_zeros]
    L = np.log(xs)
    # n_zeros pairs of basis functions [cos(g_k L), sin(g_k L)]
    cos_part = np.cos(np.outer(L, g))
    sin_part = np.sin(np.outer(L, g))
    return np.hstack([cos_part, sin_part])  # shape (K, 2*n_zeros)


def project_low_zeros(D: np.ndarray, basis: np.ndarray) -> dict:
    """Least-squares projection of D onto low-zero modes.
    Returns coefficients and residual."""
    coefs, *_ = np.linalg.lstsq(basis, D, rcond=None)
    fitted = basis @ coefs
    residual = D - fitted
    var_total = float(np.var(D))
    var_resid = float(np.var(residual))
    var_explained = var_total - var_resid
    return {
        "coefs": coefs.tolist(),
        "var_total": var_total,
        "var_explained_by_low_zeros": var_explained,
        "var_residual": var_resid,
        "frac_explained": var_explained / var_total if var_total > 0 else 0.0,
        "residual": residual,
        "fitted": fitted,
    }


# --------- main experiment ---------

def run(K: int = 1000, x_lo: float = 1e6, x_hi: float = 1e7,
        seed: int = 42, n_zeros: int = 20,
        n_gauss_controls: int = 200) -> dict:
    rng = np.random.default_rng(seed)
    print(f"K={K} log-uniform anchors in [{x_lo:.0e}, {x_hi:.0e}]")

    # Anchors: log-uniform
    xs = np.unique(np.geomspace(x_lo, x_hi, K).astype(np.int64))
    K = len(xs)  # may shrink slightly due to dedup

    print(f"Computing pi(x_k) for K={K} anchors ...")
    t0 = time.time()
    pis = primepi_array(xs)
    print(f"  pi computed in {time.time()-t0:.2f}s")

    lis = li_array(xs.astype(np.float64))
    D = D_signal(xs.astype(np.float64), pis.astype(np.float64), lis)

    mu_hat = float(np.mean(D))
    sigma_hat = float(np.std(D, ddof=0))
    print(f"D moments: mean={mu_hat:.6f}, std={sigma_hat:.6f}")

    # Empirical W_1 to fitted Gaussian
    W1_D = wasserstein1_to_normal(D, mu_hat, sigma_hat)
    print(f"W_1(D_hat, N({mu_hat:.4f}, {sigma_hat**2:.4f})) = {W1_D:.6f}")

    # Gaussian control: B = n_gauss_controls bootstrap trials.
    # PROPER COMPARISON: each trial uses its OWN sample-fitted (mu, sigma)
    # for both sample and fitted Gaussian, mirroring what we do for D.
    # This excludes the trial-vs-fitted dispersion (which would be a
    # different test, comparing absolute means/vars rather than shape).
    print(f"Computing Gaussian control: {n_gauss_controls} trials "
          f"(sample-fitted)...")
    W1_gauss_samples = []
    for b in range(n_gauss_controls):
        gsamp = rng.normal(mu_hat, sigma_hat, size=K)
        m_s, s_s = float(np.mean(gsamp)), float(np.std(gsamp, ddof=0))
        W1_gauss_samples.append(
            wasserstein1_to_normal(gsamp, m_s, s_s)
        )
    W1_gauss_samples = np.array(W1_gauss_samples)
    W1_gauss_mean = float(np.mean(W1_gauss_samples))
    W1_gauss_std = float(np.std(W1_gauss_samples))
    W1_gauss_q975 = float(np.quantile(W1_gauss_samples, 0.975))
    print(f"  Gaussian control W_1: {W1_gauss_mean:.6f} +- {W1_gauss_std:.6f}")
    print(f"  97.5% quantile: {W1_gauss_q975:.6f}")
    z_score = (W1_D - W1_gauss_mean) / W1_gauss_std
    print(f"  W_1(D) z-score vs Gaussian-control: {z_score:.2f}")

    # Bernoulli null: K Bernoulli(1/log(median x)) samples scaled to match D
    print("Computing Bernoulli null control ...")
    p_med = 1.0 / np.log(np.median(xs))
    bern_W1s = []
    for b in range(min(50, n_gauss_controls)):
        bs = rng.binomial(1, p_med, size=K).astype(float)
        bs = (bs - bs.mean()) / (bs.std() + 1e-30) * sigma_hat + mu_hat
        bern_W1s.append(wasserstein1_to_normal(bs, mu_hat, sigma_hat))
    bern_W1_mean = float(np.mean(bern_W1s))
    print(f"  Bernoulli null W_1: {bern_W1_mean:.6f}")

    # Stein's method: skew/kurt and exchangeable-pair bound
    sk = stein_skew_kurt_bounds(D)
    sep = stein_exchangeable_pair_bound(D, rng)
    print(f"Skewness: {sk['skewness']:.4f} (95% CI {sk['skew_CI95']})")
    print(f"Excess kurtosis: {sk['kurtosis_excess']:.4f} (95% CI {sk['kurt_CI95']})")
    print(f"Stein exchangeable-pair bound: {sep['stein_bound_total']:.4f}")

    # Low-zero structure decomposition
    print(f"Decomposing D against {n_zeros} lowest Riemann zeros ...")
    basis = low_zero_modes(xs.astype(np.float64), n_zeros=n_zeros)
    proj = project_low_zeros(D, basis)
    print(f"  Variance explained by {n_zeros} low zeros: "
          f"{proj['frac_explained']*100:.1f}%")

    # Residual: how Gaussian is the low-zero residual?
    resid = proj['residual']
    if np.std(resid) > 0:
        z_resid = (resid - resid.mean()) / resid.std()
        sk_resid = float(stats.skew(z_resid))
        kt_resid = float(stats.kurtosis(z_resid))
        W1_resid = wasserstein1_to_normal(resid, resid.mean(), resid.std())
    else:
        sk_resid = kt_resid = W1_resid = float('nan')
    print(f"  Residual skew={sk_resid:.4f} kurt_excess={kt_resid:.4f}")
    print(f"  Residual W_1: {W1_resid:.6f}")

    # Final A-grade / B-grade decision
    print()
    print("=" * 70)
    print("  A-GRADE TEST: W_1(D) >> 1/sqrt(K) and >> Gaussian-control")
    rate = 1.0 / np.sqrt(K)
    inflate = W1_D / W1_gauss_mean
    print(f"  1/sqrt(K) = {rate:.6f}")
    print(f"  W_1(D) / W_1(Gauss control) = {inflate:.4f}")
    print(f"  W_1(D) / sigma_hat = {W1_D / sigma_hat:.6f} (relative)")
    if z_score > 3:
        print(f"  *** Z-SCORE {z_score:.2f} > 3 — quantitative non-Gaussianity")
    else:
        print(f"  Z-SCORE {z_score:.2f} -- W_1 within Gaussian control noise")

    return {
        "K_effective": K,
        "x_range": [float(xs.min()), float(xs.max())],
        "D_mean": mu_hat,
        "D_std": sigma_hat,
        "W1_D_to_normal": W1_D,
        "W1_gauss_control_mean": W1_gauss_mean,
        "W1_gauss_control_std": W1_gauss_std,
        "W1_gauss_q975": W1_gauss_q975,
        "W1_D_z_score": z_score,
        "W1_bernoulli_null": bern_W1_mean,
        "stein_skew_kurt": sk,
        "stein_exchangeable_pair_bound": {
            k: v for k, v in sep.items()
        },
        "low_zero_decomposition": {
            "n_zeros": n_zeros,
            "var_total": proj["var_total"],
            "var_explained": proj["var_explained_by_low_zeros"],
            "var_residual": proj["var_residual"],
            "frac_explained": proj["frac_explained"],
            "residual_skew": sk_resid,
            "residual_kurt_excess": kt_resid,
            "residual_W1_to_normal": W1_resid,
        },
        "rate_1_over_sqrtK": rate,
        "W1_inflation_factor": inflate,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--K", type=int, default=1000)
    p.add_argument("--x-lo", type=float, default=1e6)
    p.add_argument("--x-hi", type=float, default=1e7)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--n-zeros", type=int, default=20)
    p.add_argument("--n-controls", type=int, default=200)
    p.add_argument("--out", type=str, default="results.json")
    args = p.parse_args()
    res = run(K=args.K, x_lo=args.x_lo, x_hi=args.x_hi, seed=args.seed,
              n_zeros=args.n_zeros, n_gauss_controls=args.n_controls)
    with open(args.out, "w") as f:
        json.dump(res, f, indent=2)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
