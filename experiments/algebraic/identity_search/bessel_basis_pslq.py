#!/usr/bin/env python3
"""
Bessel-function-basis PSLQ identity search for f(x) = pi(x) - R(x).

S29 (Session 29) thoroughly tested PSLQ in elementary + li + zeta-oscillation
bases. The strict closure "f is algebraically independent of all TESTED bases"
does not logically preclude an identity in untested bases. Bessel functions
appear naturally in the Selberg trace formula (K_{i*r}-kernels on the spectral
side) and Mellin-Barnes representations of L-functions, so a Bessel-basis test
is a strict extension of S29.

Basis tested (10 elements):
    f(x), 1, log(x), sqrt(x), li(x),
    K_0(log x), I_0(log x),                  # modified Bessel
    J_0(g1 * log x), Y_0(g1 * log x),        # oscillatory Bessel @ zeta zero
    K_0(2*pi*sqrt(log x))                    # partition-asymptotic kernel

This basis is disjoint from S29's basis (which used sin/cos(g_k * log x), not
Bessel J/Y/I/K). PSLQ is run at four x-values with cross-validation at nearby
points; a relation is "valid" only if coeff(f) != 0, residual < 1e-10, AND
cross-validation residual stays < 1e-5.

A negative outcome (no surviving relation) adds Bessel basis to the closed set,
strengthening the pseudorandomness battery for f(x).

Date: 2026-04-26 (Session 68 deep focus revisit of Task 3)
"""

import numpy as np
from mpmath import mp, mpf, log, sqrt, li, pslq, pi as mpi
from mpmath import besselk, besseli, besselj, bessely
import os, time

mp.dps = 50

GAMMA1 = mpf('14.134725141734693790457251983562470270784257115699')

DATA = np.load(os.path.join(os.path.dirname(__file__), "fx_data.npz"))
X_ARR = DATA['x']
F_ARR = DATA['f']

def get_f(x_val):
    idx = int(x_val) - 2
    if 0 <= idx < len(F_ARR):
        return mpf(str(F_ARR[idx]))
    return None

BASIS_NAMES = [
    'f(x)', '1', 'log(x)', 'sqrt(x)', 'li(x)',
    'K0(logx)', 'I0(logx)',
    'J0(g1*logx)', 'Y0(g1*logx)',
    'K0(2pi*sqrt(logx))',
]

def build_vec(xv, fx):
    xm = mpf(xv)
    L = log(xm)
    return [
        fx, mpf(1), L, sqrt(xm), li(xm),
        besselk(0, L), besseli(0, L),
        besselj(0, GAMMA1 * L), bessely(0, GAMMA1 * L),
        besselk(0, 2 * mpi * sqrt(L)),
    ]

def fmt_relation(rel):
    parts = []
    for c, n in zip(rel, BASIS_NAMES):
        if c != 0:
            parts.append(f"{c}*{n}")
    return " + ".join(parts) if parts else "(none)"

def main():
    test_points = [5000, 10000, 50000, 100000]
    print("=" * 72)
    print("Bessel-basis PSLQ identity search for f(x) = pi(x) - R(x)")
    print(f"Basis: {BASIS_NAMES}")
    print("=" * 72)

    findings = []

    for xv in test_points:
        fx = get_f(xv)
        if fx is None:
            print(f"  x={xv}: data unavailable")
            continue
        vec = build_vec(xv, fx)
        print(f"\nx = {xv}, f(x) = {float(fx):+.6f}")
        t0 = time.time()
        try:
            rel = pslq(vec, maxcoeff=10**6, maxsteps=8000)
        except Exception as e:
            print(f"  PSLQ failed: {e}")
            findings.append((xv, None, None, None, None))
            continue
        dt = time.time() - t0
        if rel is None:
            print(f"  No relation ({dt:.1f}s)")
            findings.append((xv, None, None, None, None))
            continue

        residual = abs(sum(mpf(c) * v for c, v in zip(rel, vec)))
        coeff_f = rel[0]
        rs = fmt_relation(rel)
        print(f"  relation ({dt:.1f}s): {rs}")
        print(f"  residual = {residual}")
        print(f"  coeff(f) = {coeff_f}")

        # Cross-validate at xv + 1000 (or +200 if near upper bound)
        delta = 1000 if xv + 1000 <= 100000 else 200
        xv2 = xv + delta
        fx2 = get_f(xv2)
        if fx2 is None:
            print(f"  cross-check: no data at x={xv2}")
            findings.append((xv, rs, residual, coeff_f, None))
            continue
        vec2 = build_vec(xv2, fx2)
        res2 = abs(sum(mpf(c) * v for c, v in zip(rel, vec2)))
        print(f"  cross-check at x={xv2}: residual = {res2}")
        valid = (coeff_f != 0) and (float(residual) < 1e-10) and (float(res2) < 1e-5)
        if valid:
            print("  *** SURVIVES CROSS-VALIDATION ***")
        else:
            why = []
            if coeff_f == 0:
                why.append("coeff(f)=0 (no f-bearing identity)")
            if float(residual) >= 1e-10:
                why.append(f"residual {float(residual):.2e} too large")
            if float(res2) >= 1e-5:
                why.append(f"cross-check residual {float(res2):.2e} too large")
            print(f"  FAILS: {', '.join(why)}")
        findings.append((xv, rs, residual, coeff_f, res2))

    # Random-control sanity check: replace f with same-variance gaussian noise
    print("\n" + "=" * 72)
    print("Random control: replace f with iid Gaussian, same std as f at x=10^5")
    print("=" * 72)
    rng = np.random.default_rng(20260426)
    f_std = float(np.std(F_ARR))
    for trial in range(3):
        xv = 50000
        fake_f = mpf(str(rng.normal(0.0, f_std)))
        vec = build_vec(xv, fake_f)
        try:
            rel = pslq(vec, maxcoeff=10**6, maxsteps=8000)
        except Exception:
            rel = None
        if rel is None:
            print(f"  trial {trial+1}: no relation (good — random has no identity)")
            continue
        residual = abs(sum(mpf(c) * v for c, v in zip(rel, vec)))
        # Cross-check with another random draw at xv2
        xv2 = xv + 1000
        fake_f2 = mpf(str(rng.normal(0.0, f_std)))
        vec2 = build_vec(xv2, fake_f2)
        res2 = abs(sum(mpf(c) * v for c, v in zip(rel, vec2)))
        print(f"  trial {trial+1}: PSLQ found {fmt_relation(rel)}")
        print(f"     residual at fit point = {residual}")
        print(f"     residual at cross-check = {res2} (>>1 expected if spurious)")

    # Summary
    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)
    any_valid = False
    for (xv, rs, residual, coeff_f, res2) in findings:
        if rs is None:
            print(f"  x={xv}: no relation")
            continue
        valid = coeff_f is not None and coeff_f != 0 and \
                residual is not None and float(residual) < 1e-10 and \
                res2 is not None and float(res2) < 1e-5
        any_valid = any_valid or valid
        flag = " <<<VALID>>>" if valid else ""
        print(f"  x={xv}: residual={float(residual):.2e}, coeff(f)={coeff_f}, "
              f"cross-check={float(res2) if res2 is not None else 'N/A'}{flag}")

    if any_valid:
        print("\n*** AT LEAST ONE RELATION SURVIVES — INVESTIGATE FURTHER ***")
    else:
        print("\nNo Bessel-basis identity for f(x) survives cross-validation.")
        print("Bessel/J/Y/I/K basis joins the closed set for novel-identity search.")

if __name__ == "__main__":
    main()
