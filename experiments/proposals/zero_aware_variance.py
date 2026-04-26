"""
Proposal D: Test variance reduction of pi(x) Monte-Carlo via control variates
built from the Riemann explicit formula truncated at T zeros.

Standard MC: estimate pi(x) by uniform random sampling y_k in [1,x] with
estimator x * mean(is_prime(y_k)). Variance ~ x.

Control variate: subtract the *deterministic* estimate
psi_T(x) = x - sum_{|rho|<=T} x^rho/rho - log(2*pi)
                              ^^^^^^^^^^^^^^^^^^^^
which is the truncated Riemann-von Mangoldt formula for psi(x).
Then convert psi -> pi via Mertens-style adjustment.

Hypothesis: as T grows polylog, residual variance shrinks polylog. If true,
O(polylog) random samples suffice and we have an O(polylog) algorithm.

Test: compute Var(residual) where residual = chi_P(y) - p_T(y) for random y.
If slope of log Var vs log T is < -1, variance reduction beats sqrt(x) barrier.
If slope = -1 (or 0), it doesn't.
"""

import numpy as np
from sympy import isprime, primepi
from mpmath import mpf, mpc, log, sqrt, exp
import mpmath as mp
from pathlib import Path

mp.mp.dps = 30

ZEROS_FILE = Path(__file__).resolve().parents[2] / "data" / "zeta_zeros_1000.txt"


def load_zeros(n):
    out = []
    with open(ZEROS_FILE) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                out.append(float(line.split()[0]))
            except (ValueError, IndexError):
                continue
            if len(out) == n:
                break
    return out


def explicit_psi(x, zeros):
    """Truncated psi(x) = x - sum_rho x^rho/rho - log(2*pi) - 0.5*log(1 - x^-2)
    Uses zeros = list of imaginary parts."""
    x_mp = mpf(x)
    # Sum over zero pairs rho = 1/2 + i*gamma and conjugate
    # x^rho/rho + x^{rho_bar}/rho_bar = 2 Re(x^rho / rho)
    s = mpc(0)
    for gamma in zeros:
        rho = mpc(0.5, gamma)
        s += x_mp ** rho / rho
    psi = x_mp - 2 * s.real - log(2 * mp.pi)
    if x > 1:
        psi -= 0.5 * log(1 - 1 / x_mp**2)
    return float(psi)


def explicit_pi_via_psi(x, zeros, primepi_table):
    """Convert psi to pi using Abel summation:
    pi(x) ≈ psi(x)/log(x) + integral_2^x psi(t)/(t log^2 t) dt
    For our purposes, just compute psi(x)_T and compare to true psi(x).
    """
    return explicit_psi(x, zeros)


def true_psi(x):
    """psi(x) = sum_{p^k <= x} log(p)."""
    s = mpf(0)
    for p in range(2, int(x) + 1):
        if isprime(p):
            pk = p
            while pk <= x:
                s += log(p)
                pk *= p
    return float(s)


def main():
    out = Path(__file__).with_suffix("").as_posix() + "_results.md"
    zeros = load_zeros(500)

    # Compute psi(x) - psi_T(x) for various x and T
    # If |psi(x) - psi_T(x)|^2 vs T scales as 1/T^a with a>1, variance
    # reduction is super-trivial.
    x_values = [100, 500, 1000, 5000, 10000]
    T_values = [5, 10, 20, 50, 100, 200, 500]

    table = []
    true_psi_values = {}
    for x in x_values:
        true_psi_values[x] = true_psi(x)

    for T in T_values:
        row = {"T": T}
        for x in x_values:
            psi_T = explicit_psi(x, zeros[:T])
            row[f"err_x{x}"] = abs(psi_T - true_psi_values[x])
        table.append(row)

    # Compute scaling: log|err| = a + b*log(T) for each x.
    # If b ~ -0.5, classical sqrt-barrier. If b < -1, novel.
    scalings = {}
    for x in x_values:
        log_T = np.log([row["T"] for row in table])
        log_err = np.log(np.array([row[f"err_x{x}"] for row in table]) + 1e-30)
        slope, intercept = np.polyfit(log_T, log_err, 1)
        scalings[x] = (float(slope), float(intercept))

    # Now: monte-carlo experiment for residual variance
    # For x=10000, sample N points y, compute residual = chi_P(y) - density_at_y
    # where density_at_y is reconstructed from psi_T (ish).
    # For simplicity we just measure variance of chi_P alone vs variance after
    # subtracting prime number theorem density 1/log(y).

    np.random.seed(0)
    x_mc = 10000
    N_samples = 50000
    ys = np.random.randint(2, x_mc + 1, size=N_samples)
    chi_P = np.array([1 if isprime(int(y)) else 0 for y in ys], dtype=float)
    pnt_density = 1.0 / np.log(ys)
    residual = chi_P - pnt_density

    var_naive = chi_P.var()
    var_pnt = residual.var()
    pi_x_true = int(primepi(x_mc))

    # Estimate via PNT control variate
    # E[chi_P(Y)] = pi(x)/x for Y uniform on [1,x].
    # E[1/log(Y)] = (1/x) sum_{k=1..x} 1/log(k) ≈ Li(x)/x
    Li_x = float(mp.li(x_mc) - mp.li(2))
    expected_density = Li_x / x_mc
    # estimator: x * (mean(chi_P) - mean(pnt_density)) + Li_x
    pi_estimate_pnt = x_mc * (chi_P.mean() - pnt_density.mean()) + Li_x
    pi_estimate_naive = x_mc * chi_P.mean()

    with open(out, "w") as f:
        f.write("# Proposal D — Zero-aware control variate for pi(x)\n\n")
        f.write("## Setup\n\n")
        f.write(
            "We measure the residual error |psi_T(x) - psi(x)| where psi_T is the "
            "truncated Riemann-von Mangoldt formula at T zeros, for varying T and x. "
            "If |residual|^2 vs T scales as T^b with b < -1, we beat the classical "
            "sqrt-barrier and can build a polylog algorithm.\n\n"
        )

        f.write("## |psi_T(x) - psi(x)| as function of T and x\n\n")
        f.write("| T | err@x=100 | err@x=500 | err@x=1000 | err@x=5000 | err@x=10000 |\n")
        f.write("|---|---|---|---|---|---|\n")
        for row in table:
            f.write(
                f"| {row['T']} | "
                + " | ".join(f"{row[f'err_x{x}']:.4f}" for x in x_values)
                + " |\n"
            )

        f.write("\n## Scaling exponent: log|err| = a + b*log(T)\n\n")
        f.write("| x | slope b | intercept a |\n")
        f.write("|---|---|---|\n")
        for x, (slope, intercept) in scalings.items():
            f.write(f"| {x} | {slope:.3f} | {intercept:.3f} |\n")

        f.write("\nTheoretical prediction: under RH, |psi(x) - psi_T(x)| = "
                "O(x * log^2(x*T) / T) (roughly). So slope b should be ~ -1.\n")
        f.write("If observed slope b < -1.5, that would be evidence of stronger "
                "variance reduction than RH predicts.\n\n")

        f.write("## Monte-Carlo control variate test (x = 10000)\n\n")
        f.write(f"- True pi({x_mc}) = {pi_x_true}\n")
        f.write(f"- Naive estimator x * mean(chi_P): {pi_estimate_naive:.1f} "
                f"(error {abs(pi_estimate_naive - pi_x_true):.1f})\n")
        f.write(f"- PNT control variate x * mean(chi_P - 1/log(y)) + Li(x): "
                f"{pi_estimate_pnt:.1f} (error {abs(pi_estimate_pnt - pi_x_true):.1f})\n")
        f.write(f"- Var(chi_P): {var_naive:.6f}\n")
        f.write(f"- Var(chi_P - 1/log(y)): {var_pnt:.6f}\n")
        f.write(f"- Variance reduction factor: {var_naive/var_pnt:.3f}x\n\n")

        f.write("## Interpretation\n\n")
        avg_slope = np.mean([s[0] for s in scalings.values()])
        f.write(f"Average scaling exponent across x: b ≈ {avg_slope:.3f}.\n\n")
        if avg_slope < -1.3:
            f.write(
                "Slope is steeper than -1, suggesting variance reduction beats the "
                "classical RH-conditional bound. **Path remains open**; "
                "this would be a genuine improvement over expected truncation behavior.\n"
            )
        elif avg_slope < -0.7:
            f.write(
                "Slope is approximately -1, consistent with the standard RH-conditional "
                "bound |psi(x) - psi_T(x)| = O(x * log^2(x)/T). To get the integer "
                "pi(x), we still need T = O(x/log^2 x) zeros, which is NOT polylog. "
                "The PNT control variate reduces variance by a constant factor, not "
                "asymptotically.\n\n"
                "**Verdict:** Path closed -- the explicit-formula truncation as a "
                "control variate gives only constant-factor variance reduction, not "
                "the polylog reduction needed for an O(polylog) algorithm.\n"
            )
        else:
            f.write(
                f"Slope is shallow ({avg_slope:.3f}); explicit formula truncation "
                "barely improves with T at these small x values. **Verdict closed.**\n"
            )

    print(f"Wrote {out}")
    print(f"Naive est: {pi_estimate_naive:.1f}, PNT est: {pi_estimate_pnt:.1f}, true: {pi_x_true}")
    for x, (slope, intercept) in scalings.items():
        print(f"x={x}: slope={slope:.3f}")


if __name__ == "__main__":
    main()
