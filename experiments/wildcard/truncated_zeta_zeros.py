"""
Truncated zeta polynomial zeros as a substitute for true Riemann zeros.

Idea: in the explicit formula psi(x) = x - sum_rho x^rho/rho - ...,
true zeta zeros rho require Riemann-Siegel evaluation O(sqrt(t)) per zero.
The partial sum zeta_N(s) = sum_{n=1..N} n^{-s} is a polynomial in {n^{-s}},
trivially evaluable in O(N) per point.

QUESTION: Do zeros of zeta_N(s) on the critical line approximate Riemann zeros?

If yes -> we get fast approximations to zeros, plug into explicit formula,
get pi(x) approximation.
If no -> truncated zeta has ZEROS UNRELATED to analytic continuation's zeros,
and this approach fails fundamentally.

Test:
- For N = 10, 50, 100, 500, 1000:
  - Find local minima of |zeta_N(1/2 + it)| for t in (0, 60).
  - Compare to first 10 true Riemann zeros.
- Plot the error vs N.
"""

from mpmath import mp, mpc, mpf, zetazero
mp.dps = 15  # lower precision for speed


def truncated_zeta_value(s, N):
    """zeta_N(s) = sum_{n=1..N} n^{-s} as Dirichlet polynomial."""
    return sum(mpf(n) ** (-s) for n in range(1, N + 1))


def find_truncated_zero_candidates(N, t_max=60.0, samples=1500):
    """Find local minima of |zeta_N(1/2 + it)| (candidate zero locations)."""
    half = mpf("0.5")
    ts = [t_max * (k + 1) / samples for k in range(samples)]
    abs_vals = []
    for t in ts:
        s = mpc(half, t)
        z = truncated_zeta_value(s, N)
        abs_vals.append((t, abs(z)))

    minima = []
    for i in range(1, len(abs_vals) - 1):
        a = abs_vals[i - 1][1]
        b = abs_vals[i][1]
        c = abs_vals[i + 1][1]
        if b < a and b < c and float(b) < 0.2:
            minima.append((float(abs_vals[i][0]), float(b)))
    return minima


def main():
    # First 10 true Riemann zeros' imaginary parts
    print("=" * 60)
    print("Are zeros of partial sum zeta_N(s) close to Riemann zeros?")
    print("=" * 60)
    true_zeros = [float(zetazero(k).imag) for k in range(1, 11)]
    print(f"True Riemann zeros: {[round(t, 2) for t in true_zeros]}\n")

    print(f"{'N':>6} {'#minima':>10} {'avg err':>12} {'min err':>10}")
    print("-" * 50)
    for N in [5, 10, 20, 50, 100, 200, 500]:
        minima = find_truncated_zero_candidates(N, t_max=60, samples=1200)
        approx_t = [t for t, _ in minima]
        if approx_t:
            errs = []
            for tz in true_zeros:
                closest = min(approx_t, key=lambda a: abs(a - tz))
                errs.append(abs(tz - closest))
            avg_err = sum(errs) / len(errs)
            min_err = min(errs)
        else:
            avg_err = float("inf")
            min_err = float("inf")
        print(f"{N:>6} {len(approx_t):>10} {avg_err:>12.4f} {min_err:>10.4f}")

    print()
    print("Interpretation:")
    print("- If truncated zeros approximated Riemann zeros, error -> 0 as N -> infinity.")
    print("- A non-decreasing or growing error means truncation kills the relevant zeros.")


if __name__ == "__main__":
    main()
