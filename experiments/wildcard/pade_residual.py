"""
Padé approximants of the residual delta(n) = p(n) - R^{-1}(n).

The residual is "random-like" by Fourier, wavelet, autocorrelation tests
(per project's pseudorandomness_of_pi.md). But Padé approximants exploit
*pole structure* not detected by those bases. If delta(n), viewed as samples
of a continuous function, has hidden meromorphic structure, a Padé [M/N]
should crush a polynomial of degree M+N.

Test:
- Compute p(1)..p(N) for N = 5000 via sympy (fast).
- Compute R^{-1}(n) via Newton inversion of the Riemann R function for each n.
  R(x) = sum_{k=1}^infty mu(k)/k * li(x^{1/k})   (truncated at k=20)
- Build delta(n) = p(n) - R^{-1}(n).
- Construct Pade [M/N] from delta(n_0), ..., delta(n_0 + L) at consecutive integers.
  Compare extrapolation RMSE on next L points to a polynomial of degree M+N.
- If Pade dramatically wins for some (M,N), there's exploitable structure.
"""

import math
import numpy as np
from sympy import prime, mobius
from mpmath import mp, mpf, mpc, li, polyval, exp, log

mp.dps = 25


def riemann_R(x, kmax=20):
    """R(x) = sum_{k>=1} mu(k)/k * li(x^{1/k})."""
    x = mpf(x)
    s = mpf(0)
    for k in range(1, kmax + 1):
        mu = mobius(k)
        if mu == 0:
            continue
        s += mpf(mu) / k * li(x ** (mpf(1) / k))
    return s


def R_inv(n, x_guess=None, iters=60):
    """Solve R(x) = n by bisection (robust)."""
    n_mp = mpf(n)
    # bracket: x in [2, max(2, 4 n log(n+2))]
    lo = mpf(2)
    hi = mpf(max(10, 8 * n * float(mp.log(n + 2))))
    # ensure bracket
    while riemann_R(hi) < n_mp:
        hi *= 2
    for _ in range(iters):
        mid = (lo + hi) / 2
        if riemann_R(mid) < n_mp:
            lo = mid
        else:
            hi = mid
        if (hi - lo) < mpf("1e-12"):
            break
    return (lo + hi) / 2


def pade_extrapolate(values, M, N, n_predict):
    """Given values v[0]..v[M+N], build Pade [M/N] in n, predict v[n_predict]."""
    L = M + N
    if len(values) < L + 1:
        return None
    coefs = [mpf(v) for v in values[: L + 1]]
    try:
        p, q = mp.pade(coefs, M, N)
    except Exception:
        return None
    # evaluate p(n)/q(n) at n = n_predict
    n_mp = mpf(n_predict)
    num = polyval(p[::-1], n_mp)
    den = polyval(q[::-1], n_mp)
    if abs(den) < mpf("1e-30"):
        return None
    return num / den


def poly_extrapolate(values, deg, n_predict):
    """Polynomial fit of degree deg, predict at n_predict."""
    L = deg + 1
    if len(values) < L:
        return None
    xs = np.arange(L, dtype=float)
    ys = np.array([float(v) for v in values[:L]])
    coefs = np.polyfit(xs, ys, deg)
    return float(np.polyval(coefs, n_predict))


def main():
    print("Computing residuals delta(n) = p(n) - R^{-1}(n) for n=1..200...")
    N = 200  # keep small — R_inv is expensive
    primes = []
    for n in range(1, N + 1):
        primes.append(int(prime(n)))
    print(f"  primes: p(1)={primes[0]}, p({N})={primes[-1]}")

    # warm Newton with previous solution
    deltas = []
    Rinvs = []
    x_guess = None
    for n in range(1, N + 1):
        x = R_inv(n, x_guess=x_guess)
        x_guess = float(x) + 0.5
        Rinvs.append(x)
        deltas.append(float(primes[n - 1]) - float(x))
        if n % 50 == 0:
            print(f"  n={n}: p(n)={primes[n-1]}, R^-1={float(x):.4f}, delta={deltas[-1]:.4f}")

    print("\nDelta statistics:")
    print(f"  mean = {np.mean(deltas):.4f}")
    print(f"  std  = {np.std(deltas):.4f}")
    print(f"  range = [{min(deltas):.3f}, {max(deltas):.3f}]")

    # Pade vs polynomial extrapolation
    print("\nExtrapolation test:")
    print("  Use deltas[0..L-1] to predict deltas[L..L+W-1].")
    L = 60
    W = 20
    if N < L + W:
        print(f"  Need N >= {L + W}, have {N}. Adjusting.")
        L = N // 2
        W = min(20, N - L)

    train = deltas[:L]
    test_xs = list(range(L, L + W))
    test_ys = deltas[L : L + W]

    print(f"\n  {'method':<25} {'RMSE':>12} {'max err':>12}")
    print("  " + "-" * 51)

    # polynomials
    for deg in [3, 6, 10, 15, 20]:
        if deg + 1 > L:
            continue
        preds = [poly_extrapolate(train, deg, x) for x in test_xs]
        if any(p is None for p in preds):
            continue
        errs = [abs(preds[i] - test_ys[i]) for i in range(W)]
        rmse = math.sqrt(sum(e * e for e in errs) / W)
        print(f"  poly deg={deg:<3}             {rmse:>12.4f} {max(errs):>12.4f}")

    # Pade
    for M in [3, 5, 8, 10, 15]:
        for Ndeg in [3, 5, 8, 10, 15]:
            if M + Ndeg + 1 > L:
                continue
            preds = []
            for x in test_xs:
                v = pade_extrapolate(train, M, Ndeg, x)
                preds.append(float(v) if v is not None else None)
            if any(p is None or not np.isfinite(p) for p in preds):
                continue
            errs = [abs(preds[i] - test_ys[i]) for i in range(W)]
            rmse = math.sqrt(sum(e * e for e in errs) / W)
            tag = f"Pade [{M}/{Ndeg}]"
            print(f"  {tag:<25} {rmse:>12.4f} {max(errs):>12.4f}")


if __name__ == "__main__":
    main()
