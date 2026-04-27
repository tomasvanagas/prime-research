"""
Structural explanation: derive the Wasserstein-1 plateau from the
explicit formula's low-zero contribution.

Riemann's explicit formula (under RH) gives:
    pi(x) - Li(x) = -sum_rho Li(x^rho) - log(2) + small terms.
For real-part = 1/2 zeros rho = 1/2 + i*gamma,
    Li(x^rho) = exp(rho log x) / (rho log x)  ~  x^(1/2) e^{i gamma log x} / (rho log x)
So the leading oscillatory contribution after multiplying by log(x)/sqrt(x):
    D(x) = -(log x / sqrt(x)) * sum_rho Li(x^rho) - log(2) log(x)/sqrt(x) + ...
        ~ -2 sum_{gamma>0} cos(gamma * log x - arctan(2 gamma)) / sqrt(1/4 + gamma^2)
                                                                   ~ 2/gamma for large gamma
Sum over the FIRST few gammas dominates the structural non-Gaussianity.

Compare:
- Empirical D(x_k) for K log-uniform anchors
- Synthetic D_n(x) = -2 sum_{rho=1..n} ... computed from Odlyzko zeros
- Wasserstein distance of EACH to fitted Gaussian
"""

import numpy as np
import sympy
from scipy import stats, special
from stein_wasserstein_pi import (
    primepi_array, li_array, D_signal, wasserstein1_to_normal
)

# First 50 imag parts of nontrivial Riemann zeros (Odlyzko)
ODLYZKO_GAMMAS_50 = np.array([
    14.134725141734693, 21.022039638771555, 25.010857580145688,
    30.424876125859513, 32.935061587739189, 37.586178158825671,
    40.918719012147495, 43.327073280914999, 48.005150881167160,
    49.773832477672302, 52.970321477714460, 56.446247697063394,
    59.347044002602353, 60.831778524609810, 65.112544048081606,
    67.079810529494173, 69.546401711173979, 72.067157674481907,
    75.704690699083933, 77.144840068874805, 79.337375020249366,
    82.910380854086030, 84.735492980517050, 87.425274613125229,
    88.809111207634465, 92.491899270558484, 94.651344040519887,
    95.870634228245309, 98.831194218193198, 101.317851005731391,
    103.725538040478339, 105.446623052147350, 107.168611184276645,
    111.029535543171830, 111.874659177002281, 114.320220915452712,
    116.226680321519285, 118.790782866264250, 121.370125001581826,
    122.946829293551812, 124.256818554343445, 127.516683879596495,
    129.578704199956021, 131.087688531366664, 133.497737203717945,
    134.756509753373328, 138.116042055126725, 139.736208952121477,
    141.123707404038466, 143.111845808910235,
])


def explicit_formula_D(x: np.ndarray, n_zeros: int = 50) -> np.ndarray:
    """
    Approximate D(x) from the leading explicit-formula expansion:
        D(x) = -2 sum_{k=1..n} cos(gamma_k log x - phi_k) / sqrt(1/4 + gamma_k^2)
    where phi_k = arctan(2 gamma_k) accounts for the real-part-1/2
    correction in 1/(1/2 + i gamma_k).

    The constant -log(2) log(x)/sqrt(x) is included as the bias.
    """
    g = ODLYZKO_GAMMAS_50[:n_zeros]
    L = np.log(x)
    rho_abs = np.sqrt(0.25 + g ** 2)  # |rho|
    phi = np.arctan(2 * g)  # arg(1/2 + i*gamma) = arctan(2 gamma)
    # Leading oscillatory: -2 sum cos(gamma L - phi) / |rho|  (factor of 2
    # for both rho and rho-conjugate pair)
    # For each k: contribution = -2 cos(g_k * L - phi_k) / |rho_k|
    cos_part = np.cos(np.outer(L, g) - phi)  # shape (K, n_zeros)
    weights = 2.0 / rho_abs  # shape (n_zeros,)
    D_osc = -np.sum(cos_part * weights, axis=1)
    bias = -np.log(2) * L / np.sqrt(x)
    return D_osc + bias


def analyse():
    K = 5000
    print(f"Computing empirical D for K={K} ...")
    xs = np.unique(np.geomspace(1e6, 1e7, K).astype(np.int64))
    K = len(xs)
    pis = primepi_array(xs)
    lis = li_array(xs.astype(np.float64))
    D_emp = D_signal(xs.astype(np.float64), pis.astype(np.float64), lis)

    mu_emp = D_emp.mean()
    sigma_emp = D_emp.std()
    skew_emp = stats.skew(D_emp)
    kurt_emp = stats.kurtosis(D_emp)
    W1_emp = wasserstein1_to_normal(D_emp, mu_emp, sigma_emp)
    print(f"\nEMPIRICAL D:")
    print(f"  mean={mu_emp:.4f}  std={sigma_emp:.4f}  skew={skew_emp:.4f}  "
          f"kurt={kurt_emp:.4f}")
    print(f"  W_1(D, N(mu_emp, sigma_emp^2)) = {W1_emp:.6f}")

    # Theoretical D from n=1, 5, 10, 20, 50 zeros
    print(f"\nEXPLICIT-FORMULA-BASED D (truncated to n zeros):")
    print(f"  {'n':>3} | {'mean':>10} {'std':>10} {'skew':>10} {'kurt':>10} "
          f"{'W_1':>10} {'corr_emp':>10}")
    for n in [1, 2, 3, 5, 10, 20, 50]:
        D_th = explicit_formula_D(xs.astype(np.float64), n_zeros=n)
        m_th = D_th.mean()
        s_th = D_th.std()
        sk = stats.skew(D_th)
        kt = stats.kurtosis(D_th)
        W1 = wasserstein1_to_normal(D_th, m_th, s_th)
        # Correlation with empirical
        corr = float(np.corrcoef(D_emp, D_th)[0, 1])
        print(f"  {n:>3} | {m_th:>10.4f} {s_th:>10.4f} {sk:>10.4f} {kt:>10.4f} "
              f"{W1:>10.6f} {corr:>10.4f}")

    # Direct comparison: residual after subtracting D_th
    print(f"\nRESIDUAL ANALYSIS (D_emp minus D_th(n=20)):")
    D_th20 = explicit_formula_D(xs.astype(np.float64), n_zeros=20)
    # Center each
    D_emp_c = D_emp - mu_emp
    D_th20_c = D_th20 - D_th20.mean()
    # Best-fit scalar
    alpha = np.dot(D_emp_c, D_th20_c) / np.dot(D_th20_c, D_th20_c)
    print(f"  Best-fit alpha (scalar match): {alpha:.4f}")
    resid = D_emp_c - alpha * D_th20_c
    print(f"  Residual: std={resid.std():.4f} skew={stats.skew(resid):.4f} "
          f"kurt={stats.kurtosis(resid):.4f}")
    print(f"  Variance of D_emp explained: "
          f"{1.0 - resid.var()/D_emp_c.var():.4f}")

    print(f"\nKEY OBSERVATION:")
    print(f"  W_1 plateau value c_emp = {W1_emp:.6f}")
    # Predict from low zeros: W_1 of D_th(n=20) to fitted Gaussian
    D_th20_W1 = wasserstein1_to_normal(D_th20, D_th20.mean(), D_th20.std())
    print(f"  W_1 of explicit-formula D_th(20 zeros) to fit Gaussian = "
          f"{D_th20_W1:.6f}")
    print(f"  Ratio empirical/theoretical = {W1_emp / D_th20_W1:.4f}")


if __name__ == "__main__":
    analyse()
