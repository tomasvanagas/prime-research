#!/usr/bin/env python3
"""
Extended-basis identity search for f(x) = pi(x) - R(x).

Session 29 closed the Novel Identity Search direction with 7 experiments
covering elementary functions, two zeta-zero oscillations (gamma_1, gamma_2),
PSLQ functional/shift recurrences, WZ certificates, LLL minimal polynomials,
ODEs, and Volterra kernels. Three sub-bases were NOT covered:

  A. Many zeta-zero oscillations (gamma_3..gamma_10), eight extra modes.
  B. Arithmetic-function partial sums M(x)=sum mu(n), L_lambda(x)=sum lambda(n),
     Phi(x)=sum phi(n), psi(x)=sum Lambda(n) -- mainstream summatory bases
     besides psi (which Session 29 partially summed for the Chebyshev link).
  C. Mahler-type functional relations f(x), f(x^2), f(x^3) -- distinct from
     the f(2x), f(3x), f(4x) tested earlier; these probe self-similarity
     under exponential rescaling rather than linear rescaling.

Method: recompute f(x) = pi(x) - R(x) at chosen sample points using mpmath
at 60-digit precision, then run PSLQ on each basis with cross-validation at
distinct points. A "valid" relation must (i) have nonzero coefficient on f
and (ii) survive cross-check residual < 1e-6 at a different point.

Date: 2026-04-25
"""

from __future__ import annotations
import time
import numpy as np
from sympy import primepi, mobius
from mpmath import mp, mpf, log, sqrt, li, sin, cos, pslq

mp.dps = 60

# ---------------------------------------------------------------------------
# High-precision f(x) = pi(x) - R(x), R(x) = sum mu(k)/k * li(x^{1/k})
# ---------------------------------------------------------------------------

def riemann_R(x_mpf, terms: int = 60) -> mpf:
    """Riemann's R-function via the Mobius series; terms<=60 is plenty for x<=10^6."""
    s = mpf(0)
    for k in range(1, terms + 1):
        mk = int(mobius(k))
        if mk == 0:
            continue
        s += mpf(mk) / k * li(x_mpf ** (mpf(1) / k))
    return s


def f_of(x_int: int) -> mpf:
    """Compute f(x) = pi(x) - R(x) at high precision."""
    return mpf(int(primepi(x_int))) - riemann_R(mpf(x_int))


# ---------------------------------------------------------------------------
# Arithmetic-function partial sums via a single Mobius/Liouville sieve
# ---------------------------------------------------------------------------

def make_arith_sums(N: int):
    """Return arrays M[k] = sum_{n<=k} mu(n), L[k] = sum_{n<=k} lambda(n),
    PHI[k] = sum_{n<=k} phi(n), PSI[k] = sum_{n<=k} Lambda(n) for k=0..N.

    Liouville lambda(n) = (-1)^Omega(n) where Omega counts prime factors with
    multiplicity. Lambda(n) = log(p) if n=p^k, else 0.
    """
    # Smallest prime factor sieve
    spf = list(range(N + 1))
    for i in range(2, int(N ** 0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, N + 1, i):
                if spf[j] == j:
                    spf[j] = i

    mu = [0] * (N + 1)
    lam = [0] * (N + 1)
    phi = [0] * (N + 1)
    Lam = [mpf(0)] * (N + 1)
    mu[1] = 1
    lam[1] = 1
    phi[1] = 1

    for n in range(2, N + 1):
        # factor n via spf
        m = n
        primes = {}
        while m > 1:
            p = spf[m]
            c = 0
            while m % p == 0:
                m //= p
                c += 1
            primes[p] = c
        # mu(n)
        if any(c > 1 for c in primes.values()):
            mu[n] = 0
        else:
            mu[n] = (-1) ** len(primes)
        # lambda(n) = (-1)^Omega
        lam[n] = (-1) ** sum(primes.values())
        # phi(n)
        ph = n
        for p in primes:
            ph = ph // p * (p - 1)
        phi[n] = ph
        # Lambda(n)
        if len(primes) == 1:
            (p,) = primes.keys()
            Lam[n] = mpf(int(p)).__class__(int(p))  # placeholder
            Lam[n] = mp.log(mpf(int(p)))

    M = [0] * (N + 1)
    L = [0] * (N + 1)
    PHI = [0] * (N + 1)
    PSI = [mpf(0)] * (N + 1)
    for n in range(1, N + 1):
        M[n] = M[n - 1] + mu[n]
        L[n] = L[n - 1] + lam[n]
        PHI[n] = PHI[n - 1] + phi[n]
        PSI[n] = PSI[n - 1] + Lam[n]
    return M, L, PHI, PSI


# ---------------------------------------------------------------------------
# Zeta zero loader (real parts of imaginary parts gamma_k)
# ---------------------------------------------------------------------------

def load_zeros(n: int = 20):
    with open("data/zeta_zeros_2000.txt") as fh:
        out = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            out.append(mpf(line))
            if len(out) >= n:
                break
    return out


# ---------------------------------------------------------------------------
# PSLQ helper with cross-validation
# ---------------------------------------------------------------------------

def cross_check(rel, vec_at_x2, eps_thresh: float = 1e-6) -> tuple[bool, mpf]:
    res2 = abs(sum(mpf(c) * v for c, v in zip(rel, vec_at_x2)))
    return float(res2) < eps_thresh, res2


def run_pslq(label: str, basis_names: list[str], vecs: dict[int, list], maxcoeff=10**8):
    """For each x in vecs, run PSLQ on vecs[x]. Cross-validate any relation with
    nonzero f-coefficient at every other x in vecs.
    """
    print(f"\n  ----- {label} -----")
    print(f"    Basis ({len(basis_names)} terms): {basis_names}")
    points = sorted(vecs.keys())
    survivors = []
    for x in points:
        vec = vecs[x]
        t0 = time.time()
        try:
            rel = pslq(vec, maxcoeff=maxcoeff, maxsteps=10000)
        except Exception as e:
            print(f"    x={x}: PSLQ failed -- {e}")
            continue
        elapsed = time.time() - t0
        if rel is None:
            print(f"    x={x}: no relation ({elapsed:.1f}s)")
            continue
        residual = abs(sum(mpf(c) * v for c, v in zip(rel, vec)))
        coeff_f = rel[0]
        terms = [f"{rel[i]}*{basis_names[i]}" for i in range(len(rel)) if rel[i] != 0]
        print(f"    x={x}: residual={float(residual):.2e}  coeff(f)={coeff_f}")
        print(f"          relation: {' + '.join(terms)[:160]}")
        if coeff_f == 0 or float(residual) > 1e-10:
            print(f"          -> spurious or non-tight; skip cross-check")
            continue
        # Cross-validate at all other x
        all_pass = True
        for x2 in points:
            if x2 == x:
                continue
            ok, res2 = cross_check(rel, vecs[x2])
            print(f"          cross-check x={x2}: residual={float(res2):.2e}  {'PASS' if ok else 'FAIL'}")
            if not ok:
                all_pass = False
        if all_pass:
            survivors.append((x, rel, residual))
            print(f"          *** SURVIVES CROSS-VALIDATION at all {len(points)-1} points ***")
    return survivors


# ===========================================================================
# Main
# ===========================================================================

def main():
    print("=" * 78)
    print("EXTENDED-BASIS IDENTITY SEARCH for f(x) = pi(x) - R(x)")
    print("Sub-bases not covered by Session 29's seven experiments")
    print("=" * 78)

    # ---- Choose sample points ----
    # Section A test points (extended zeta basis): need only f(x) and basis terms
    section_A_points = [1000, 5000, 20000, 50000]
    # Section B (arithmetic partial sums): use moderate x to keep sieve cheap
    section_B_points = [1000, 2000, 5000, 10000]
    # Section C (Mahler): need f(x), f(x^2), f(x^3); cap at x^3 <= 100000
    section_C_points = [10, 20, 30, 40, 46]   # 46^3 = 97336 < 100000

    all_xs_needed = sorted(set(
        section_A_points + section_B_points
        + [x for x in section_C_points]
        + [x ** 2 for x in section_C_points]
        + [x ** 3 for x in section_C_points]
    ))

    print(f"\nComputing f(x) at {len(all_xs_needed)} points (mp.dps={mp.dps})...")
    fx_cache: dict[int, mpf] = {}
    for x in all_xs_needed:
        t0 = time.time()
        fx_cache[x] = f_of(x)
        if time.time() - t0 > 0.5:
            print(f"    f({x}) = {float(fx_cache[x]):+.6f}  ({time.time()-t0:.2f}s)")
    print(f"  done. cached {len(fx_cache)} values.")

    # ---- Load extended zeta zeros ----
    zeros = load_zeros(10)
    print(f"\nLoaded {len(zeros)} zeta zeros (gamma_1..gamma_10).")

    # ====================================================================
    # SECTION A: Extended zeta-zero oscillation basis
    # ====================================================================
    print("\n" + "=" * 78)
    print("SECTION A: PSLQ over basis with gamma_1..gamma_10 oscillations")
    print("=" * 78)

    A_basis_names = (
        ["f(x)", "1", "log(x)", "sqrt(x)", "li(x)/x", "li(sqrt(x))"]
        + [f"sin(g{k+1}*log(x))" for k in range(10)]
        + [f"cos(g{k+1}*log(x))" for k in range(10)]
    )

    A_vecs: dict[int, list] = {}
    for x in section_A_points:
        xm = mpf(x)
        logx = log(xm)
        vec = [
            fx_cache[x],
            mpf(1),
            logx,
            sqrt(xm),
            li(xm) / xm,
            li(sqrt(xm)),
        ]
        for g in zeros:
            vec.append(sin(g * logx))
        for g in zeros:
            vec.append(cos(g * logx))
        A_vecs[x] = vec

    sectionA_survivors = run_pslq("A: 26-term basis (6 elementary + 20 zeta-osc)",
                                   A_basis_names, A_vecs, maxcoeff=10**10)

    # ====================================================================
    # SECTION B: Arithmetic-function partial sums
    # ====================================================================
    print("\n" + "=" * 78)
    print("SECTION B: PSLQ vs arithmetic partial sums M(x), L_lambda(x), Phi(x), psi(x)")
    print("=" * 78)
    print("  (sieving Mobius/Liouville/totient/Lambda up to N=10000)")

    t0 = time.time()
    M, L_lam, PHI, PSI = make_arith_sums(max(section_B_points))
    print(f"  sieve done in {time.time()-t0:.1f}s")

    # Normalize all integer-valued sums by sqrt(x) or x so PSLQ cannot find a
    # trivial integer combination among them; force any short relation to
    # involve f or be irrational.
    B_basis_names = [
        "f(x)", "log(x)", "1/sqrt(x)",
        "M(x)/sqrt(x)", "L_lambda(x)/sqrt(x)",
        "(Phi(x) - 3x^2/pi^2)/x", "(psi(x)-x)/sqrt(x)",
        "M(x)*log(x)/sqrt(x)", "L_lambda(x)*log(x)/sqrt(x)",
    ]
    B_vecs: dict[int, list] = {}
    pi_const = mp.pi
    for x in section_B_points:
        xm = mpf(x)
        logx = log(xm)
        sqx = sqrt(xm)
        vec = [
            fx_cache[x],
            logx,
            1 / sqx,
            mpf(M[x]) / sqx,
            mpf(L_lam[x]) / sqx,
            (mpf(PHI[x]) - 3 * xm**2 / pi_const**2) / xm,
            (PSI[x] - xm) / sqx,
            mpf(M[x]) * logx / sqx,
            mpf(L_lam[x]) * logx / sqx,
        ]
        B_vecs[x] = vec

    sectionB_survivors = run_pslq("B: 11-term basis (arithmetic partial sums)",
                                   B_basis_names, B_vecs, maxcoeff=10**10)

    # ====================================================================
    # SECTION C: Mahler-type functional relation f(x), f(x^2), f(x^3)
    # ====================================================================
    print("\n" + "=" * 78)
    print("SECTION C: Mahler functional test  f(x), f(x^2), f(x^3) + elementary")
    print("=" * 78)

    # Drop integer 'x' from basis to avoid trivial x = k*1 hits. Use only
    # transcendental basis elements so any short integer combination must
    # express f via the others.
    C_basis_names = ["f(x)", "f(x^2)", "f(x^3)", "log(x)", "sqrt(x)/log(x)", "1/log(x)", "log(x)^2"]
    C_vecs: dict[int, list] = {}
    for x in section_C_points:
        xm = mpf(x)
        logx = log(xm)
        vec = [
            fx_cache[x],
            fx_cache[x ** 2],
            fx_cache[x ** 3],
            logx,
            sqrt(xm) / logx,
            1 / logx,
            logx ** 2,
        ]
        C_vecs[x] = vec

    sectionC_survivors = run_pslq("C: 8-term Mahler basis", C_basis_names, C_vecs,
                                   maxcoeff=10**10)

    # ====================================================================
    # Summary
    # ====================================================================
    print("\n" + "=" * 78)
    print("SUMMARY")
    print("=" * 78)
    total_survivors = (
        len(sectionA_survivors) + len(sectionB_survivors) + len(sectionC_survivors)
    )
    print(f"  Section A survivors (extended zeta basis):       {len(sectionA_survivors)}")
    print(f"  Section B survivors (arithmetic partial sums):   {len(sectionB_survivors)}")
    print(f"  Section C survivors (Mahler functional):         {len(sectionC_survivors)}")
    print(f"  TOTAL surviving cross-validated relations:       {total_survivors}")
    if total_survivors == 0:
        print("\n  *** No identity in any extended sub-basis. ***")
        print("  Strengthens Session 29 closure of Novel Identity Search.")
    else:
        print("\n  *** RELATION SURVIVED -- INVESTIGATE ***")
    print("Done.")


if __name__ == "__main__":
    main()
