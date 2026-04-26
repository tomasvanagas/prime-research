"""
Prime Theta Function Modular Test
==================================

Hypothesis: If a "prime theta function" of the form
    Theta_P(t) = sum_{p prime, p <= N} f(p, t)
admits a Jacobi-like functional equation Theta_P(1/t) = g(t) * Theta_P(t)
(possibly up to lower-order terms), then we can evaluate Theta_P(t) at
small t (where many primes contribute) using its value at large t (where
only a few do). This would yield a polylog evaluator and -- after Mellin
inversion -- a polylog pi(x) algorithm.

Theoretical background: Riemann's classical theta theta(t) = sum_n exp(-pi n^2 t)
satisfies theta(1/t) = sqrt(t) theta(t) due to Poisson summation over Z.
Poisson summation generalizes to lattices, but PRIMES are not a lattice.
The prime zeta function P(s) = sum p^{-s} has a NATURAL BOUNDARY at Re(s)=0
(Landau 1903), which structurally PREDICTS no Jacobi-like functional equation.

This experiment EMPIRICALLY tests several theta-like prime kernels for
functional-equation-like behavior, using high-precision arithmetic.

Test kernels:
  K1: Theta_P(t) = sum_{p prime, p<=N} exp(-pi p^2 t)
  K2: Theta_P(t) = sum_{p prime, p<=N} exp(-pi p t)
  K3: Theta_P(t) = sum_{p prime, p<=N} (-1)^k * exp(-pi p^2 t)  (alternating)
  K4: prime-counting variant: sum_{p prime, p<=N} 1/(1+p^2 t)

For each, we test:
  Test A: Power law Theta_P(1/t) ~ t^alpha * Theta_P(t)
          Fit alpha by least-squares in log-space.
  Test B: Goodness of fit (R^2). For modular structure, R^2 ~ 1.
          For pseudo-random behavior, R^2 << 1 or alpha varies wildly.
  Test C: Compare to RANDOM control with same density.

If any kernel passes Test B with R^2 > 0.95 across multiple t-ranges,
that's a hit deserving deeper investigation.
"""
import numpy as np
from sympy import sieve
import time


def primes_up_to(N):
    return np.array(list(sieve.primerange(2, N + 1)), dtype=np.float64)


def theta_K1(primes, t):
    # Use stable evaluation: subtract max exponent
    a = -np.pi * primes**2 * t
    if len(a) == 0:
        return 0.0
    return float(np.exp(a).sum())


def theta_K2(primes, t):
    a = -np.pi * primes * t
    if len(a) == 0:
        return 0.0
    return float(np.exp(a).sum())


def theta_K3(primes, t):
    # Alternating sign by index
    a = -np.pi * primes**2 * t
    if len(a) == 0:
        return 0.0
    signs = (-1.0) ** np.arange(len(primes))
    return float((signs * np.exp(a)).sum())


def theta_K4(primes, t):
    # Rational kernel
    return float(np.sum(1.0 / (1.0 + primes**2 * t)))


KERNELS = {
    "K1_gauss_sq": theta_K1,
    "K2_exp_lin":  theta_K2,
    "K3_alt_sq":   theta_K3,
    "K4_rational": theta_K4,
}


def random_control(N, density, seed=0):
    """Return a random subset of [2, N] with same density as primes <= N."""
    rng = np.random.default_rng(seed)
    n_prime = max(1, int(round(density * (N - 1))))
    return rng.choice(np.arange(2, N + 1), size=n_prime, replace=False).astype(np.float64)


def fit_power_law(t_vals, theta_vals, theta_inv_vals):
    """
    Fit log(Theta_P(1/t)) = alpha * log(t) + log(C) + log(Theta_P(t))
    -> log(Theta_P(1/t)/Theta_P(t)) = alpha * log(t) + log(C)
    Returns (alpha, C, R^2) of the linear fit.
    """
    eps = 1e-300
    safe_inv = np.maximum(np.abs(theta_inv_vals), eps)
    safe = np.maximum(np.abs(theta_vals), eps)
    y = np.log(safe_inv / safe)
    x = np.log(t_vals)
    if len(x) < 3:
        return np.nan, np.nan, np.nan
    A = np.vstack([x, np.ones_like(x)]).T
    sol, residuals, rank, sv = np.linalg.lstsq(A, y, rcond=None)
    alpha, log_C = sol
    y_pred = A @ sol
    ss_res = float(np.sum((y - y_pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(alpha), float(np.exp(log_C)), r2


def run_kernel(name, kernel, points, t_vals):
    """Compute Theta_P(t) and Theta_P(1/t) for each t."""
    pairs = []
    for t in t_vals:
        theta_t = kernel(points, t)
        theta_inv = kernel(points, 1.0 / t)
        if not (np.isfinite(theta_t) and np.isfinite(theta_inv)):
            continue
        if theta_t == 0 or theta_inv == 0:
            continue
        pairs.append((t, theta_t, theta_inv))
    if len(pairs) < 4:
        return None
    pairs = np.array(pairs)
    return fit_power_law(pairs[:, 0], pairs[:, 1], pairs[:, 2])


def main():
    print("=" * 76)
    print("Prime Theta Modular Test")
    print("=" * 76)

    Ns = [1000, 10000, 100000]
    # t-range chosen so that exp(-pi p^2 t) doesn't underflow / overflow
    # Need 1/t in similar range; geometric grid around t=1.

    for N in Ns:
        primes = primes_up_to(N)
        density = len(primes) / (N - 1)
        ctrl = random_control(N, density, seed=42)
        # Avoid t too small/large for floating point.
        # For K1 with N=10^5: pi*p^2*t ~ pi*10^10*t ; need < 700 for exp.
        # So t > 700 / (pi * N^2) and t < (pi * N^2) / 700.
        t_min = max(1e-12, 800.0 / (np.pi * N**2))
        t_max = min(1e12, (np.pi * N**2) / 800.0)
        # Use 30 log-spaced t values within (t_min, t_max), avoiding t=1.
        t_vals = np.geomspace(t_min, t_max, 30)
        # Drop t very close to 1 (would give Theta_P(t) == Theta_P(1/t) trivially)
        t_vals = t_vals[np.abs(np.log(t_vals)) > 0.5]

        print(f"\n--- N = {N}, density = {density:.5f}, primes = {len(primes)} ---")
        print(f"    t range: {t_vals.min():.3e} -> {t_vals.max():.3e}")
        print(f"{'kernel':14}  {'alpha (P)':>11}  {'R^2 (P)':>9}  {'alpha (R)':>11}  {'R^2 (R)':>9}  verdict")

        for name, kernel in KERNELS.items():
            res_p = run_kernel(name, kernel, primes, t_vals)
            res_c = run_kernel(name, kernel, ctrl, t_vals)
            if res_p is None or res_c is None:
                print(f"{name:14}  numerically degenerate")
                continue
            alpha_p, _, r2_p = res_p
            alpha_c, _, r2_c = res_c

            verdict = "no fit"
            if r2_p > 0.95:
                if r2_c > 0.95 and abs(alpha_p - alpha_c) < 0.05:
                    verdict = "MATCHES random"  # not prime-specific
                elif abs(r2_p - r2_c) > 0.1:
                    verdict = "PRIME-SPECIFIC fit"
                else:
                    verdict = "fits, but ~ random"
            elif r2_p > 0.7:
                verdict = "weak fit"
            print(f"{name:14}  {alpha_p:>11.4f}  {r2_p:>9.4f}  {alpha_c:>11.4f}  {r2_c:>9.4f}  {verdict}")

    # --- Sanity check: classical theta(t) = sum_n exp(-pi n^2 t)
    # should give R^2 ~ 1 with alpha = 1/2 (as theta(1/t) = sqrt(t) theta(t))
    print("\n--- SANITY CHECK: classical theta over ALL integers ---")
    N = 10000
    all_n = np.arange(1, N + 1, dtype=np.float64)

    def classical_theta(t):
        a = -np.pi * all_n**2 * t
        return 1.0 + 2 * float(np.exp(a).sum())  # n=0 contributes 1, full theta is sym

    t_vals = np.geomspace(1e-3, 1e3, 30)
    t_vals = t_vals[np.abs(np.log(t_vals)) > 0.5]
    pairs = []
    for t in t_vals:
        a = classical_theta(t)
        b = classical_theta(1.0 / t)
        if a > 0 and b > 0:
            pairs.append((t, a, b))
    pairs = np.array(pairs)
    if len(pairs) >= 4:
        alpha, C, r2 = fit_power_law(pairs[:, 0], pairs[:, 1], pairs[:, 2])
        print(f"classical theta:  alpha = {alpha:.4f}  (expect 0.5)   R^2 = {r2:.6f}  (expect ~1)")
        print(f"                  prefactor C = {C:.4f}  (expect 1.0)")

    print("\nDone.")


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"\nTotal time: {time.time() - t0:.2f}s")
