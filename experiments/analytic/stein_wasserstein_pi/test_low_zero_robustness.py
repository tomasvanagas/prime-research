"""
Robustness test: is the 'n=3 zeros match' coincidence or real?

Approach 1: average W_1 over many random sub-anchors
Approach 2: empirical W_1 vs theoretical W_1 at varying x ranges
Approach 3: compute exact W_1 plateau as the *expected* W_1 when D
            converges to its explicit-formula limit.
"""

import numpy as np
import sympy
from scipy import stats
from stein_wasserstein_pi import (
    primepi_array, li_array, D_signal, wasserstein1_to_normal
)
from structural_explanation import explicit_formula_D, ODLYZKO_GAMMAS_50


def main():
    # Use a finer K=20000 sampling but restricted to the empirical x range
    # so we can subsample various K values from a pre-computed pi table.
    K_ref = 10000
    print(f"Computing reference D for K={K_ref} ...")
    xs = np.unique(np.geomspace(1e6, 1e7, K_ref).astype(np.int64))
    K_ref = len(xs)
    pis = primepi_array(xs)
    lis = li_array(xs.astype(np.float64))
    D_emp = D_signal(xs.astype(np.float64), pis.astype(np.float64), lis)

    # 1) Plateau at K=K_ref
    mu_emp, sigma_emp = D_emp.mean(), D_emp.std()
    W1_emp = wasserstein1_to_normal(D_emp, mu_emp, sigma_emp)
    print(f"\nReference: K={K_ref}, mean={mu_emp:.4f}, std={sigma_emp:.4f}")
    print(f"W_1(D_emp) = {W1_emp:.6f}")
    print(f"Skew = {stats.skew(D_emp):.4f}, Kurt_excess = {stats.kurtosis(D_emp):.4f}")

    # 2) Theoretical at MUCH FINER sampling — does W_1(D_th(n)) converge?
    print(f"\n--- Theoretical W_1 at K_ref = {K_ref}, varying n_zeros ---")
    n_list = list(range(1, 21)) + [25, 30, 40, 50]
    print(f"  {'n':>3} | {'std_th':>8} {'skew_th':>9} {'kurt_th':>9} "
          f"{'W_1(D_th)':>11}")
    for n in n_list:
        D_th = explicit_formula_D(xs.astype(np.float64), n_zeros=n)
        m_th, s_th = D_th.mean(), D_th.std()
        sk, kt = stats.skew(D_th), stats.kurtosis(D_th)
        W1 = wasserstein1_to_normal(D_th, m_th, s_th)
        marker = "  <--" if abs(W1 - W1_emp) < 0.001 else ""
        print(f"  {n:>3} | {s_th:>8.4f} {sk:>9.4f} {kt:>9.4f} {W1:>11.6f}{marker}")

    # 3) Average W_1 over MANY random x ranges of width log(10)
    print(f"\n--- W_1 of D_emp restricted to varying log-uniform sub-windows ---")
    # K=10000 spans log range [6, 7] in log10, width 1.0.
    # Try sub-windows of width 0.5 (log10): [6.0, 6.5], [6.25, 6.75], ...
    rng = np.random.default_rng(31)
    n_windows = 10
    W1s = []
    starts = np.linspace(6.0, 6.5, n_windows)
    log_xs = np.log10(xs)
    for s0 in starts:
        mask = (log_xs >= s0) & (log_xs <= s0 + 0.5)
        D_sub = D_emp[mask]
        if len(D_sub) < 100:
            continue
        m, sg = D_sub.mean(), D_sub.std()
        W1 = wasserstein1_to_normal(D_sub, m, sg)
        sk = stats.skew(D_sub)
        kt = stats.kurtosis(D_sub)
        W1s.append(W1)
        print(f"  log10(x) in [{s0:.2f}, {s0+0.5:.2f}]: K={len(D_sub):4d}  "
              f"W_1={W1:.5f}  skew={sk:+.3f} kurt={kt:+.3f}")
    print(f"  Mean W_1 over windows: {np.mean(W1s):.5f}")
    print(f"  Std W_1 over windows:  {np.std(W1s):.5f}")

    # 4) Theoretical at SAME windows: does the structural prediction track?
    print(f"\n--- Theoretical D_th(50) at SAME sub-windows ---")
    D_th50 = explicit_formula_D(xs.astype(np.float64), n_zeros=50)
    W1s_th = []
    for s0 in starts:
        mask = (log_xs >= s0) & (log_xs <= s0 + 0.5)
        D_sub = D_th50[mask]
        if len(D_sub) < 100:
            continue
        m, sg = D_sub.mean(), D_sub.std()
        W1 = wasserstein1_to_normal(D_sub, m, sg)
        sk = stats.skew(D_sub)
        kt = stats.kurtosis(D_sub)
        W1s_th.append(W1)
        print(f"  log10(x) in [{s0:.2f}, {s0+0.5:.2f}]: "
              f"W_1={W1:.5f}  skew={sk:+.3f} kurt={kt:+.3f}")
    print(f"  Mean W_1 (theoretical): {np.mean(W1s_th):.5f}")

    # Correlation of W_1 across windows: theoretical vs empirical
    if len(W1s) == len(W1s_th):
        corr = np.corrcoef(W1s, W1s_th)[0, 1]
        print(f"  Corr of W_1 (emp vs theoretical) across windows: {corr:.4f}")

    # 5) "Random-phase" control: replace gamma_k with uniform random phases
    #    in the explicit formula and compute W_1.  This isolates the
    #    arithmetic-zero effect from generic-cosine-sum.
    print(f"\n--- Random-phase control: same gamma weights, random phases ---")
    g = ODLYZKO_GAMMAS_50[:50]
    rho_abs = np.sqrt(0.25 + g ** 2)
    weights = 2.0 / rho_abs
    L = np.log(xs.astype(np.float64))
    W1s_rand = []
    sk_rand = []
    kt_rand = []
    for trial in range(50):
        phases = rng.uniform(0, 2 * np.pi, size=50)
        D_rand = -np.sum(weights * np.cos(np.outer(L, g) - phases), axis=1)
        m, sg = D_rand.mean(), D_rand.std()
        W1 = wasserstein1_to_normal(D_rand, m, sg)
        W1s_rand.append(W1)
        sk_rand.append(stats.skew(D_rand))
        kt_rand.append(stats.kurtosis(D_rand))
    print(f"  Random-phase W_1: mean={np.mean(W1s_rand):.5f} +- {np.std(W1s_rand):.5f}")
    print(f"  Random-phase skew: mean={np.mean(sk_rand):+.3f} +- {np.std(sk_rand):.3f}")
    print(f"  Random-phase kurt: mean={np.mean(kt_rand):+.3f} +- {np.std(kt_rand):.3f}")
    z_emp_vs_rand = (W1_emp - np.mean(W1s_rand)) / np.std(W1s_rand)
    print(f"  Z(W_1 emp vs random-phase): {z_emp_vs_rand:+.2f}")


if __name__ == "__main__":
    main()
