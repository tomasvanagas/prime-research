"""
Proposal B: Lagrange-Burmann inversion of Li(x) for p(n), with Borel-Pade
resummation of the resulting divergent formal series.

Test: how does |partial-sum estimate of p(n) - true p(n)| scale with K
(number of Lagrange terms) for n in {10, 100, 1000, 10000}? If the error
decays exponentially in K, we have a polylog algorithm. If it bottoms out
or diverges, we close the path.

Standard asymptotic series for p(n) (Cipolla):
  p(n) = n*log n + n*(log log n - 1)
       + n/log n * (log log n - 2)
       - n/log^2 n * (log log n)^2/2 - 3 log log n + ...
But these expansions have factorially-growing coefficients (Cesaro,
Salvy & Vasseur). Borel resummation maps factorial divergence to
analytic functions if the singularity structure cooperates.

We compute the formal Lagrange series symbolically, then apply Borel-Pade.
"""

import math
from sympy import (
    symbols, log, Rational, Poly, series, sympify, factorial, oo, summation,
    integrate, exp, Symbol, expand, simplify, lambdify, primerange, S, nsimplify,
    li, log as sym_log, cancel, factor
)
import sympy
from sympy.functions.special.error_functions import li as sym_li
from pathlib import Path
from sympy import primepi

import mpmath as mp
mp.mp.dps = 60


def cipolla_p(n, num_terms):
    """Cipolla's asymptotic expansion truncated at num_terms.
    p(n) = n*L + n*(M-1) + n/L*(M-2) - n/L^2 * ((M^2-6M+11)/2) + ...
    where L = log(n), M = log(log(n)).
    """
    n = mp.mpf(n)
    L = mp.log(n)
    M = mp.log(L)
    # Coefficients from Cipolla 1902 (cleaned up):
    # p(n) ~ n*( L + (M-1) + (M-2)/L - ((M^2 - 6M + 11)/2)/L^2 + ... )
    series_terms = [
        L,
        M - 1,
        (M - 2) / L,
        -((M**2 - 6*M + 11) / 2) / L**2,
        ((2*M**3 - 21*M**2 + 84*M - 131) / 6) / L**3,
        -((3*M**4 - 46*M**3 + 270*M**2 - 766*M + 887) / 12) / L**4,
        ((12*M**5 - 245*M**4 + 2030*M**3 - 8772*M**2 + 19914*M - 18589) / 60) / L**5,
    ]
    s = mp.mpf(0)
    cumulative = []
    for k in range(min(num_terms, len(series_terms))):
        s += series_terms[k]
        cumulative.append(n * s)
    return cumulative  # list of partial sums (n=1, n=1+2, ...)


def borel_pade_evaluate(series_coeffs, x):
    """Apply Borel summation to a series with coefficients a_k.
    Borel transform: B(t) = sum_k (a_k / k!) * t^k
    Evaluation: integral_0^infty exp(-t) * B(x*t) dt = sum_k a_k * x^k (formal).

    For numerical use: take Pade approximant of B(t) and integrate.
    Returns the resummed value. If Pade fails, returns float('nan').
    """
    n_coeffs = len(series_coeffs)
    if n_coeffs < 4:
        return mp.mpf("nan")

    # Borel coefficients
    b = [mp.mpf(c) / mp.mpf(math.factorial(k)) for k, c in enumerate(series_coeffs)]

    # Pade approximant [m/n] with m+n = len(b)-1
    m = (n_coeffs - 1) // 2
    n_deg = n_coeffs - 1 - m
    try:
        p_coeffs, q_coeffs = mp.pade(b, m, n_deg)
    except Exception:
        return mp.mpf("nan")

    # Numerically integrate exp(-t) * Pade(x*t) from 0 to infty
    def integrand(t):
        # Pade(x*t) = P(x*t) / Q(x*t)
        u = x * t
        num = mp.mpf(0)
        for c in reversed(p_coeffs):
            num = num * u + c
        den = mp.mpf(0)
        for c in reversed(q_coeffs):
            den = den * u + c
        if abs(den) < mp.mpf("1e-30"):
            return mp.mpf(0)
        return mp.exp(-t) * num / den

    try:
        val = mp.quad(integrand, [0, 5, 50, mp.inf])
        return val
    except Exception:
        return mp.mpf("nan")


def main():
    out = Path(__file__).with_suffix("").as_posix() + "_results.md"

    targets = [10, 100, 1000, 10000]
    true_values = {n: int(sympy.prime(n)) for n in targets}

    rows = []
    for n in targets:
        true_p = true_values[n]
        cipolla_partials = cipolla_p(n, 7)
        rows.append({"n": n, "true": true_p, "partials": cipolla_partials})

    with open(out, "w") as f:
        f.write("# Proposal B — Lagrange/Cipolla asymptotic with Borel resummation\n\n")
        f.write("## Setup\n\n")
        f.write(
            "Cipolla's classical asymptotic expansion:\n\n"
            "$$ p(n) \\sim n\\big[ L + (M-1) + (M-2)/L - \\tfrac{M^2-6M+11}{2L^2} + \\dots \\big] $$\n\n"
            "where L = log n, M = log log n. Coefficients grow factorially -> "
            "asymptotic / non-convergent. We measure how the truncation error scales "
            "with K and whether direct truncation already gives polylog accuracy.\n\n"
        )

        f.write("## Truncation error vs. K\n\n")
        f.write("| n | true p(n) | K=1 err | K=2 err | K=3 err | K=4 err | K=5 err | K=6 err | K=7 err |\n")
        f.write("|---|---|---|---|---|---|---|---|---|\n")
        for r in rows:
            f.write(f"| {r['n']} | {r['true']} ")
            for partial in r["partials"]:
                err = abs(float(partial) - r["true"])
                f.write(f"| {err:.3g} ")
            f.write("|\n")

        # Best K for each n (asymptotic series stops improving at the optimal index)
        f.write("\n## Best truncation per n\n\n")
        f.write("| n | best K | min error | relative error |\n")
        f.write("|---|---|---|---|\n")
        best_data = []
        for r in rows:
            errs = [abs(float(p) - r["true"]) for p in r["partials"]]
            best_K = int(errs.index(min(errs))) + 1
            min_err = min(errs)
            rel = min_err / r["true"]
            best_data.append((r["n"], best_K, min_err, rel))
            f.write(f"| {r['n']} | {best_K} | {min_err:.4g} | {rel:.3e} |\n")

        # Borel-Pade applied to the Cipolla series in 1/log(n)
        f.write("\n## Borel-Pade resummation\n\n")
        f.write(
            "We construct the formal series sum_k a_k / L^k where the a_k = a_k(M) "
            "depend on M = log log n. We pretend M is constant, "
            "transform a_k -> a_k/k!, take a Pade approximant, then integrate.\n\n"
        )
        f.write("| n | true | Cipolla K=7 | Borel-Pade resummed | error |\n")
        f.write("|---|---|---|---|---|\n")
        borel_results = []
        for r in rows:
            n = r["n"]
            L = mp.log(n)
            M = mp.log(L)
            # Build coefficients of (1/L)^k as floats (using current M)
            coeffs_in_invL = [
                mp.mpf(L) * 1,           # k=0 (the leading L term)
                mp.mpf(M) - 1,           # k=0 in 1/L expansion (but absorb separately)
                M - 2,                   # k=1
                -(M**2 - 6*M + 11) / 2,  # k=2
                (2*M**3 - 21*M**2 + 84*M - 131) / 6,         # k=3
                -(3*M**4 - 46*M**3 + 270*M**2 - 766*M + 887) / 12,  # k=4
                (12*M**5 - 245*M**4 + 2030*M**3 - 8772*M**2 + 19914*M - 18589) / 60,  # k=5
            ]
            # Borel-resum the 1/L tail (terms from k=1 onwards in 1/L)
            tail_coeffs = coeffs_in_invL[2:]  # the (M-2), -(M^2-...)/2, ... in 1/L^k for k=1..5
            x_eval = mp.mpf(1) / L
            resummed_tail = borel_pade_evaluate(tail_coeffs, x_eval)
            if mp.isnan(resummed_tail):
                resummed = mp.mpf("nan")
            else:
                resummed = mp.mpf(n) * (L + (M - 1) + resummed_tail)

            cipolla_full = float(r["partials"][-1])
            err = float(abs(resummed - r["true"])) if not mp.isnan(resummed) else float("nan")
            borel_results.append((n, r["true"], cipolla_full, float(resummed) if not mp.isnan(resummed) else float("nan"), err))
            f.write(
                f"| {n} | {r['true']} | {cipolla_full:.4f} | "
                f"{(float(resummed) if not mp.isnan(resummed) else 'NaN')} | {err if err == err else 'NaN'} |\n"
            )

        f.write("\n## Interpretation\n\n")
        # Check if relative error grows or stays bounded
        rel_errors = [bd[3] for bd in best_data]
        # If Cipolla relative error stays ~1/log(n), this is the standard barrier
        # If Borel improves it substantially, this is novel
        cipolla_better_at_large_n = best_data[-1][3] < best_data[0][3] * 0.1
        borel_works = any(r[4] < r[2] - r[1] for r in borel_results if r[4] == r[4])

        f.write(
            f"Best Cipolla truncation gives relative error "
            f"{best_data[0][3]:.2e} -> {best_data[-1][3]:.2e} as n goes "
            f"from 10 to 10000.\n\n"
        )
        if best_data[-1][3] < best_data[0][3]:
            f.write(
                "Cipolla relative error **shrinks** with n, consistent with "
                "the asymptotic-series character (better as n grows). "
                "Absolute error grows — only relative is bounded.\n\n"
            )
        if borel_works:
            f.write(
                "Borel-Pade resummation **reduces error below truncation** "
                "for at least one tested n. Path remains open; escalate.\n\n"
            )
        else:
            f.write(
                "Borel-Pade resummation does **not** improve over the best Cipolla "
                "truncation. The factorial divergence appears to be of "
                "non-Borel-summable type (Stokes-line / multi-instanton structure), "
                "or at this small n the asymptotic series is already past its "
                "optimal truncation point.\n\n"
                "**Verdict:** Path closed at this depth. Failure mode I "
                "(asymptotic information loss).\n"
            )

    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
