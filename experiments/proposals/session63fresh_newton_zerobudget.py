"""
P4: Iterated Newton with progressive zero-budget for solving pi(x) = n.

p(n) is the inverse of pi: pi(p(n)) = n. Newton iteration for solving
pi(x) - n = 0 converges quadratically (digits double each iteration).

The key idea: an APPROXIMATE pi using only K zeta zeros has error
~ x / (sqrt(x) log x * K) (heuristic) — decreases roughly with 1/K. So:
  - Iteration 0:  use few zeros (low cost)
  - Iteration 1:  double the zero budget; the Newton step uses better pi
  - Iteration k:  use K_k = 2^(k+1) zeros

By doubling, total zero work is geometric: K_total = sum 2^k = 2 * K_final.
If K_final = polylog(n), total is polylog.

Test:
  For n in {100, 1000, 5000, 10000}, run Newton starting from x_0 = R^{-1}(n),
  using K_k zeros at iteration k. Report:
    - zero budget per iteration
    - residual |x_k - p_n_truth| per iteration
    - whether residual halves quadratically
    - final K_total to reach |x - p_n| < 1
"""

import sys
import numpy as np
from mpmath import mp, mpf, mpc, log as mlog, ei, riemannr
from sympy.ntheory import primerange

mp.prec = 80

ZEROS_PATH = "data/zeta_zeros_1000.txt"


def load_zeros():
    with open(ZEROS_PATH) as f:
        return [mpf(line.strip()) for line in f if line.strip()]


def true_p_n(n_max):
    """Returns sorted list of first n_max primes."""
    upper = max(20, int(n_max * (np.log(n_max) + np.log(np.log(n_max + 2)) + 2)))
    while True:
        primes = list(primerange(2, upper + 1))
        if len(primes) >= n_max:
            return primes[:n_max]
        upper = int(upper * 1.5)


def riemann_R_inverse(n_target):
    """Numerical inverse of R(x) = n_target via Newton (no zero corrections)."""
    n = float(n_target)
    if n < 5:
        x = mpf(2 * (n + 1))
    else:
        x = mpf(n) * mlog(mpf(n)) + mpf(n) * mlog(mlog(mpf(max(2.0, float(n)))))
    for _ in range(60):
        Rx = riemannr(x)
        Rp = mpf(1) / mlog(x)
        delta = (Rx - mpf(n_target)) / Rp
        x = x - delta
        if abs(delta) < mpf("1e-30"):
            break
    return x


def pi_approx_zeros(x, zeros, K):
    """Approximate pi(x) = R(x) - sum_{j<=K} 2 Re Ei(rho_j log x)."""
    log_x = mlog(mpf(x))
    R = riemannr(mpf(x))
    total = mpc(0)
    for j in range(min(K, len(zeros))):
        g = zeros[j]
        rho = mpc(mpf("0.5"), g)
        rho_bar = mpc(mpf("0.5"), -g)
        total += ei(rho * log_x) + ei(rho_bar * log_x)
    return float(R - total.real)


def pi_prime_approx_zeros(x, zeros, K):
    """Derivative d pi_K(x) / dx where pi_K = R - sum_{j<=K} 2 Re Ei(rho log x).
       d/dx Ei(rho log x) = e^(rho log x) / (rho log x) * (rho / x)
                          = x^rho / (x log x)
       So d/dx [Ei(rho log x) + Ei(rho_bar log x)] = (x^rho + x^{rho_bar})/(x log x)
                                                   = 2 Re(x^rho) / (x log x).
       d/dx R(x) is approximately 1 / log(x).
    """
    log_x = mlog(mpf(x))
    R_prime = mpf(1) / log_x  # leading term
    total = mpf(0)
    for j in range(min(K, len(zeros))):
        g = zeros[j]
        rho = mpc(mpf("0.5"), g)
        x_rho = mpc(x) ** rho
        # d/dx [Ei(rho log x)] = x^rho / (x log x)
        # add x^rho + x^{rho_bar} = 2 Re(x^rho)
        total += 2 * x_rho.real / (mpf(x) * log_x)
    return float(R_prime - total)


def newton_with_progressive_zeros(n_target, zeros, p_truth, max_iter=8):
    """Newton iteration for pi(x) = n_target, with K_k = 2^(k+1) zeros at step k."""
    x = riemann_R_inverse(n_target)
    history = []
    K_used_total = 0
    for k in range(max_iter):
        K_k = 2 ** (k + 1)
        K_k = min(K_k, len(zeros))
        pi_x = pi_approx_zeros(float(x), zeros, K_k)
        pi_p = pi_prime_approx_zeros(float(x), zeros, K_k)
        if abs(pi_p) < 1e-15:
            break
        x_new = mpf(x) - mpf(pi_x - n_target) / mpf(pi_p)
        residual = float(x_new - p_truth)  # how far from true p(n)
        history.append({
            "iter": k,
            "K": K_k,
            "x_before": float(x),
            "pi_x_K_zeros": pi_x,
            "x_after": float(x_new),
            "residual": residual,
        })
        K_used_total += K_k
        x = x_new
        # stop if residual is essentially zero
        if abs(residual) < 0.5:
            break
    return history, K_used_total, float(x)


def main():
    print("=" * 70)
    print("P4: Iterated Newton with progressive zero-budget")
    print("=" * 70)

    zeros = load_zeros()
    print(f"\nLoaded {len(zeros)} zeta zeros.")

    test_ns = [100, 1000, 5000, 10000]
    primes = true_p_n(max(test_ns))

    print()
    for n in test_ns:
        p_n_truth = mpf(primes[n - 1])
        print(f"n = {n:>5}, p(n) = {primes[n-1]} (truth)")
        history, K_total, x_final = newton_with_progressive_zeros(
            n, zeros, p_n_truth, max_iter=8
        )
        print(f"  start: x_0 = R^(-1)({n}) = {history[0]['x_before']:.4f}, "
              f"residual_0 = {history[0]['x_before'] - float(p_n_truth):+.4f}")
        for h in history:
            print(f"  iter {h['iter']:>1}: K={h['K']:>3}  "
                  f"x={h['x_after']:>13.5f}  pi_K(x)={h['pi_x_K_zeros']:>9.4f}  "
                  f"residual={h['residual']:+.4f}")
        print(f"  final x = {x_final:.4f}, true p(n) = {primes[n-1]}, "
              f"err = {x_final - primes[n-1]:+.4f}")
        print(f"  total zeros used = {K_total}")
        print()

    print("Interpretation:")
    print("  - If residual halves quadratically per iteration AND K_total = O(log^c n),")
    print("    polylog evaluation may be feasible.")
    print("  - If residual stagnates above some floor, the K-zero approximation is")
    print("    insufficient regardless of Newton.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
