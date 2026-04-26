"""
Free-probability moment test on delta(x) = pi(x) - Li(x).

Hypothesis (Item 3 of the 2026-04-26 fresh-perspective brainstorm):
treat the contribution of each non-trivial zeta zero rho to delta(x)
at log-scale log x as a free random variable. The R-transform composes
additively under free convolution. If the empirical distribution of
delta(x)/sqrt(x) over a dyadic window x in [N, 2N] matches a known
free-convolution form (free Poisson, free Wigner / semicircle,
Marchenko-Pastur), it suggests an analytic shortcut to sums-over-zeros
via R-transform inversion -- a path no other angle has explored.

Concretely we compute, on each dyadic window:

    z_n = (pi(n) - Li(n)) * sqrt(log n) / sqrt(n)        (Selberg scaling)

and report the empirical moments E[z_n^k] for k = 1..6, contrasting
against:

    Gaussian              : even k -> (k-1)!! ; odd -> 0
    Free Wigner (semic.)  : k=2 m -> Catalan(m) ; odd -> 0
    Free Poisson(1)       : k=2 m -> Catalan(m) ; k=2 m+1 -> Catalan(m)
    Marchenko-Pastur(c=1) : same Catalan-shape as Wigner

(Free Wigner moments assume variance 1; we rescale.)

Falsifier:
- If empirical moments match free-Wigner (even Catalan, odd 0) within
  a few percent at large window, *that* is the structural shortcut.
- If they match Gaussian, classic CLT applies (Selberg-style) and no
  freeness shortcut is available.
- If they match neither, document the residual.
"""

import math
import time

import numpy as np
import sympy


# --------------------------------------------------------------------------- #
# helpers                                                                     #
# --------------------------------------------------------------------------- #


def sieve_pi(n_max: int) -> np.ndarray:
    """pi(k) for k = 0..n_max (cumulative count)."""
    sieve = np.ones(n_max + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(n_max**0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = False
    return np.cumsum(sieve.astype(np.int64))


def li_array(n_max: int) -> np.ndarray:
    """
    Logarithmic integral Li(n) for integer n = 0..n_max.

    Use the relation Li(n) ≈ Ei(log n), and for stability we evaluate
    via mpmath at n = 2..n_max in chunks. Set Li(0) = Li(1) = 0 by
    convention (Li(2) is the customary base point).
    """
    out = np.zeros(n_max + 1, dtype=np.float64)
    # use math.lgamma-free: just use mpmath for numerical accuracy
    from mpmath import li as mp_li  # type: ignore[import-not-found]

    for k in range(2, n_max + 1):
        out[k] = float(mp_li(k))
    return out


def catalan(m: int) -> int:
    return math.comb(2 * m, m) // (m + 1)


def double_factorial(k: int) -> int:
    res = 1
    while k > 1:
        res *= k
        k -= 2
    return res


def gaussian_moment(k: int) -> float:
    if k % 2:
        return 0.0
    return float(double_factorial(k - 1))


def wigner_moment(k: int) -> float:
    """Variance-1 semicircle: even k = 2m -> Catalan(m), odd 0."""
    if k % 2:
        return 0.0
    m = k // 2
    return float(catalan(m))


def free_poisson_moment(k: int, lam: float = 1.0) -> float:
    """
    Compound free Poisson with rate lam, mean 1: moments are
    sum_{r=1..k} N(k,r) * lam^r where N(k,r) is the Narayana number.
    For lam=1 these reduce to total Catalan(k).
    """
    if k == 0:
        return 1.0
    if abs(lam - 1.0) < 1e-12:
        # sum_r N(k,r) = Catalan(k)
        return float(catalan(k))
    total = 0.0
    for r in range(1, k + 1):
        narayana = math.comb(k, r) * math.comb(k, r - 1) // k
        total += narayana * (lam**r)
    return total


def marchenko_pastur_moment(k: int, c: float = 1.0) -> float:
    """
    M-P distribution with shape c (ratio): moments are
        m_k = sum_{r=1..k} (1/k) * C(k,r) * C(k,r-1) * c^r = N_{k,r}*c^r
    For c=1 reduces to Catalan(k).
    """
    if k == 0:
        return 1.0
    total = 0.0
    for r in range(1, k + 1):
        narayana = math.comb(k, r) * math.comb(k, r - 1) // k
        total += narayana * (c**r)
    return total


# --------------------------------------------------------------------------- #
# main experiment                                                             #
# --------------------------------------------------------------------------- #


def empirical_moments(z: np.ndarray, kmax: int = 6) -> list[float]:
    """Sample moments E[z^k], k=1..kmax."""
    out = []
    for k in range(1, kmax + 1):
        out.append(float(np.mean(z**k)))
    return out


def run() -> None:
    print("=== Free-probability moment test on delta(x) = pi(x) - Li(x) ===\n")
    t0 = time.time()
    n_max = 1 << 18  # 262144
    print(f"Sieving pi up to {n_max:,} ...")
    pi_arr = sieve_pi(n_max)
    print(f"  done in {time.time() - t0:.2f} s; pi({n_max}) = {pi_arr[-1]}")

    t1 = time.time()
    print("Computing Li(n) on integer grid ...")
    li_arr = li_array(n_max)
    print(f"  done in {time.time() - t1:.2f} s")

    delta = pi_arr.astype(np.float64) - li_arr  # signed fluctuation
    # standard normalization: variance of delta over (N, 2N) ~ N/log(N)^2
    # Selberg-style scaling: z_n = delta(n) * log(n) / sqrt(n)
    # (this rescales so var(z_n) ~ const independent of N, asymptotically)
    n_idx = np.arange(n_max + 1, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        scale = np.where(n_idx >= 2, np.sqrt(n_idx) / np.log(n_idx), 1.0)
        z = delta / scale  # so var(z_n) ~ O(1) by Selberg

    # report moments on dyadic windows
    print("\n-- Empirical moments of delta(n) * log(n) / sqrt(n) --")
    print("(this scaling gives variance ~ O(1) -- Cramer's heuristic)\n")
    print(
        f"{'window':>20}  {'count':>8}  {'m1':>8}  {'m2':>8}  "
        f"{'m3':>8}  {'m4':>8}  {'m5':>9}  {'m6':>10}"
    )
    rows = []
    for j in range(8, 18):
        a, b = 1 << j, 1 << (j + 1)
        if b > n_max:
            break
        zj = z[a:b]
        zj = zj[np.isfinite(zj)]
        # standardize to mean 0, variance 1 within window so we can
        # compare cleanly to standardized free distributions
        mu, sig = float(np.mean(zj)), float(np.std(zj))
        if sig < 1e-12:
            continue
        zs = (zj - mu) / sig
        m = empirical_moments(zs, kmax=6)
        rows.append((j, a, b, len(zj), mu, sig, m))
        print(
            f"  [{a:>7},{b:>7})  {len(zj):>8}  "
            f"{m[0]:>+8.4f}  {m[1]:>8.4f}  "
            f"{m[2]:>+8.4f}  {m[3]:>8.4f}  "
            f"{m[4]:>+9.4f}  {m[5]:>10.4f}"
        )

    # show the reference moment sequences (variance-1)
    print("\n-- Reference moments (variance-1, mean 0) --\n")
    print(f"{'distribution':>22}  {'m1':>8}  {'m2':>8}  {'m3':>8}  "
          f"{'m4':>8}  {'m5':>9}  {'m6':>10}")
    for name, fn in [
        ("Gaussian (CLT)", lambda k: gaussian_moment(k)),
        ("Wigner / semicircle", lambda k: wigner_moment(k)),
        ("Free Poisson (lam=1)", lambda k: free_poisson_moment(k, 1.0)),
        ("Marchenko-Pastur c=1", lambda k: marchenko_pastur_moment(k, 1.0)),
    ]:
        # all 4 reference distributions have raw moments above; convert
        # to "standardized" (variance 1, mean 0) for fair comparison
        raw = [fn(k) for k in range(1, 7)]
        # WARNING: free Poisson(1) has mean 1, var 1 -> shift+rescale below
        if name == "Free Poisson (lam=1)":
            mu_ref = 1.0
            var_ref = 1.0
            sig_ref = math.sqrt(var_ref)
        elif name == "Marchenko-Pastur c=1":
            mu_ref = 1.0
            var_ref = 1.0
            sig_ref = math.sqrt(var_ref)
        else:
            mu_ref = 0.0
            var_ref = raw[1]
            sig_ref = math.sqrt(var_ref)

        # standardize raw moments via central-moment recurrence
        # (mu_k of (X-mu)/sig in terms of raw moments of X)
        def std_moment(k: int) -> float:
            # E[((X-mu)/sig)^k] = (1/sig^k) sum_j C(k,j) (-mu)^(k-j) E[X^j]
            tot = 0.0
            for jj in range(k + 1):
                ej = 1.0 if jj == 0 else raw[jj - 1]
                tot += math.comb(k, jj) * ((-mu_ref) ** (k - jj)) * ej
            return tot / (sig_ref**k)

        std = [std_moment(k) for k in range(1, 7)]
        print(
            f"  {name:>22}  {std[0]:>+8.4f}  {std[1]:>8.4f}  "
            f"{std[2]:>+8.4f}  {std[3]:>8.4f}  "
            f"{std[4]:>+9.4f}  {std[5]:>10.4f}"
        )

    # build a *closeness* table: L_inf distance between the empirical
    # standardized moment vector at the largest dyadic window and each
    # reference, on (m3, m4, m5, m6) only (m1, m2 are forced by
    # standardization).
    if rows:
        last = rows[-1]
        m_emp = last[6]
        cmps = {}
        for name, mu_ref, var_ref, raw_fn in [
            ("Gaussian", 0.0, 1.0, gaussian_moment),
            ("Wigner", 0.0, 1.0, wigner_moment),
            ("FreePoisson(1)", 1.0, 1.0, lambda k: free_poisson_moment(k, 1.0)),
            ("MP(c=1)", 1.0, 1.0, lambda k: marchenko_pastur_moment(k, 1.0)),
        ]:
            sig_ref = math.sqrt(var_ref)
            raw = [raw_fn(k) for k in range(1, 7)]
            std = []
            for k in range(1, 7):
                tot = 0.0
                for jj in range(k + 1):
                    ej = 1.0 if jj == 0 else raw[jj - 1]
                    tot += math.comb(k, jj) * ((-mu_ref) ** (k - jj)) * ej
                std.append(tot / (sig_ref**k))
            # compare m3..m6
            d = max(abs(m_emp[k] - std[k]) for k in (2, 3, 4, 5))
            cmps[name] = d
        print("\n-- Closeness (L_inf on m3..m6) of empirical vs reference --")
        print(f"  largest-window empirical: {[round(x, 4) for x in m_emp]}")
        for name, d in sorted(cmps.items(), key=lambda kv: kv[1]):
            print(f"    {name:>16}: L_inf = {d:.4f}")

        winner = min(cmps, key=cmps.get)
        print(f"\n  Best fit: {winner}  (L_inf = {cmps[winner]:.4f})")
        print("  -- if Gaussian wins, classical CLT (no freeness shortcut)")
        print("  -- if Wigner wins, asymptotic freeness exists -> R-transform path")
        print("  -- if FreePoisson / MP wins, more exotic structure")

    print(f"\nTotal time: {time.time() - t0:.2f} s")


if __name__ == "__main__":
    run()
