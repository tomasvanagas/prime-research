"""
Borel resummation of the Riemann explicit formula.

Idea: The truncated Riemann sum
    psi_T(x) = x  -  sum_{|gamma|<T} x^rho/rho  -  log(2*pi)  -  (1/2)log(1-x^{-2})
converges to psi(x) only as T -> infinity, at rate ~1/T (or worse since terms
have size O(x/gamma)). Borel-resum the *partial-sum sequence* a_K = S_K - S_{K-1}
in K (an integer index over zeros sorted by ordinate). If the Borel transform
B(z) = sum_K a_K z^K / K! analytically continues to a tractable function on R+,
the Borel sum  int_0^infty e^{-z} B(z) dz  may converge to the true sum.

Test:
- Use 200, 500, 1000 zeros from data/.
- Compute psi_K(x) for K = 1..1000 and several test x.
- Compare:
    (a) raw partial sum
    (b) Cesaro mean of partial sums (well-known)
    (c) Borel-Pade summation:  fit Pade [M/N] to B(z) = sum a_K z^K / K!,
        then evaluate int_0^infty e^{-z} Pade(z) dz numerically.
- Verdict: does Borel-Pade beat Cesaro?
"""

import os
import sys
import math
import mpmath
from mpmath import mp, mpf, mpc, exp, log, quad, polyval

mp.dps = 40


def load_zeros(path, K):
    """Load first K positive imaginary parts of zeta zeros."""
    with open(path) as f:
        gammas = [mpf(line.strip()) for line in f if line.strip()]
    return gammas[:K]


def true_psi(x):
    """Reference Chebyshev psi(x) via direct sum over prime powers (small x only)."""
    # psi(x) = sum_{p^k <= x} log p
    from sympy import primerange
    s = mpf(0)
    for p in primerange(2, int(x) + 1):
        pk = p
        while pk <= x:
            s += mpf(int(p)).__class__(math.log(p))  # convert via float OK at this dps
            pk *= p
    # rebuild precisely
    s = mpf(0)
    for p in primerange(2, int(x) + 1):
        lp = mp.log(p)
        pk = p
        while pk <= x:
            s += lp
            pk *= p
    return s


def explicit_psi_partial(x, gammas, K):
    """psi_K(x) approximation from explicit formula using K zeros (each contributing rho and conj rho)."""
    x_mp = mpf(x)
    main = x_mp - mp.log(2 * mp.pi) - mpf("0.5") * mp.log(1 - 1 / (x_mp * x_mp))
    # zero contributions: each zero rho = 1/2 + i*gamma contributes  x^rho/rho + x^{rho_bar}/rho_bar
    s = mpf(0)
    half = mpf("0.5")
    log_x = mp.log(x_mp)
    sqrt_x = mp.sqrt(x_mp)
    for k in range(K):
        g = gammas[k]
        # x^rho = x^(1/2) * exp(i*gamma*log x) = sqrt(x) * (cos + i sin)
        c = mp.cos(g * log_x)
        s_ = mp.sin(g * log_x)
        # x^rho/rho + x^rhobar/rhobar = 2*Re(x^rho/rho)
        # rho = 1/2 + i g => 1/rho = (1/2 - i g)/(1/4 + g^2)
        denom = mpf("0.25") + g * g
        re_inv = half / denom
        im_inv = -g / denom
        # Re(x^rho * 1/rho) = sqrt(x)*(c*re_inv - s_*im_inv)
        contrib = 2 * sqrt_x * (c * re_inv - s_ * im_inv)
        s -= contrib  # explicit formula has MINUS sign on zero sum
    return main + s


def borel_pade_sum(coefs, M, N):
    """Borel-Pade resummation: form B(z) = sum_K coefs[K]*z^K/K!, fit Pade [M/N], integrate."""
    # B(z) is power series with coefficients coefs[K]/K!
    bcoefs = []
    fact = mpf(1)
    for K, a in enumerate(coefs):
        if K > 0:
            fact *= K
        bcoefs.append(a / fact)
    # Pade approximation
    try:
        p, q = mp.pade(bcoefs, M, N)
    except Exception as e:
        return None
    # Integrate e^{-z} * P(z)/Q(z) on [0, inf)
    def integrand(z):
        try:
            num = polyval(p[::-1], z)
            den = polyval(q[::-1], z)
            if abs(den) < mpf("1e-30"):
                return mpf(0)
            return mp.exp(-z) * num / den
        except Exception:
            return mpf(0)
    try:
        return quad(integrand, [0, mp.inf])
    except Exception:
        return None


def cesaro_mean(partials):
    """Cesaro mean of partial sums."""
    n = len(partials)
    cesums = []
    csum = mpf(0)
    for i, p in enumerate(partials):
        csum += p
        cesums.append(csum / (i + 1))
    return cesums


def main():
    data = "/apps/aplikacijos/prime-research/data/zeta_zeros_1000.txt"
    K_max = 200  # use first 200 zeros (each is an "amplitude tier")
    gammas = load_zeros(data, K_max)
    print(f"Loaded {len(gammas)} zero ordinates.")

    test_xs = [50, 100, 500, 1000]

    for x in test_xs:
        print("\n" + "=" * 60)
        print(f"x = {x}")
        true_p = true_psi(x)
        print(f"  true psi(x)        = {mp.nstr(true_p, 12)}")

        # cumulative partial sums of explicit-formula approximation
        partials = []
        for K in range(1, K_max + 1):
            psi_K = explicit_psi_partial(x, gammas, K)
            partials.append(psi_K)

        raw_K = partials[-1]
        print(f"  raw psi_{K_max}(x)        = {mp.nstr(raw_K, 12)}, err = {mp.nstr(raw_K - true_p, 6)}")

        # Cesaro
        cesums = cesaro_mean(partials)
        print(f"  Cesaro mean        = {mp.nstr(cesums[-1], 12)}, err = {mp.nstr(cesums[-1] - true_p, 6)}")

        # Borel-Pade on the increment sequence a_K = partials[K] - partials[K-1]
        increments = [partials[0]] + [partials[i] - partials[i - 1] for i in range(1, len(partials))]
        # try a few (M, N) Pade orders
        best_err = None
        best_MN = None
        best_val = None
        for M in [5, 10, 15, 20]:
            for N in [5, 10, 15, 20]:
                if M + N + 1 > len(increments):
                    continue
                val = borel_pade_sum(increments, M, N)
                if val is None or not mp.isfinite(val):
                    continue
                # Borel-Pade output is *sum* of increments, so it equals partials[K_max] in the limit
                err = abs(val - true_p)
                if best_err is None or err < best_err:
                    best_err = err
                    best_MN = (M, N)
                    best_val = val
        if best_val is not None:
            print(f"  Borel-Pade best    = {mp.nstr(best_val, 12)}, err = {mp.nstr(best_err, 6)}, (M,N)={best_MN}")
        else:
            print("  Borel-Pade FAILED (all orders unstable)")

        # Convergence rate of partials
        # rate ~ slope of log|err| vs log K
        import math as _m
        Ks = list(range(20, K_max + 1, 20))
        errs = [abs(partials[k - 1] - true_p) for k in Ks]
        if all(e > 0 for e in errs):
            xs = [_m.log(k) for k in Ks]
            ys = [float(mp.log(e)) for e in errs]
            n = len(xs)
            mx = sum(xs) / n
            my = sum(ys) / n
            num = sum((xs[i] - mx) * (ys[i] - my) for i in range(n))
            den = sum((xs[i] - mx) ** 2 for i in range(n))
            slope = num / den if den != 0 else float("nan")
            print(f"  raw partial conv rate: |err| ~ K^{slope:.3f}")


if __name__ == "__main__":
    main()
