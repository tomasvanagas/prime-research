"""
Proposal 4 — Truncated explicit formula on Dirichlet characters mod q.

Tests pi(x; q, a) reconstruction from finite zero truncation of L(s, chi mod q).
Uses *known* zeros stored on disk (data/zeta_zeros_2000.txt) for the principal
character and *computes* zeros of L(s, chi mod 12) for non-principal chars
via mpmath's mp.zeros — but limited to a small number to keep runtime short.

For q=12, phi(12)=4 characters total (principal chi_0 plus three real chars).

Question: at fixed small T, does rounded pi_AP(x;12,a) match the true sieve
count?

(Fast version: limit non-principal zero search to T_max=40 and at most 30
zeros per character.)
"""
from __future__ import annotations
import math

import numpy as np
import mpmath as mp
from scipy.special import expi
from sympy import primerange


mp.mp.dps = 25


def Li(x):
    """Real-valued logarithmic integral."""
    return float(expi(math.log(x)))


def Li_complex(z):
    """Logarithmic integral for complex z, via mpmath.li."""
    return mp.li(z)


def true_pi_AP(x, q, a):
    return sum(1 for p in primerange(2, x + 1) if p % q == a)


def char_mod_12_table():
    """4 real characters mod 12 (each value in {-1, 0, 1})."""
    coprimes = [1, 5, 7, 11]
    chars = [
        {1: 1, 5: 1, 7: 1, 11: 1},     # principal chi_0
        {1: 1, 5: -1, 7: -1, 11: 1},
        {1: 1, 5: 1, 7: -1, 11: -1},
        {1: 1, 5: -1, 7: 1, 11: -1},
    ]
    return chars, coprimes


def L_chi_mod_12(s_complex, chi):
    """L(s, chi) for chi mod 12 via Hurwitz-zeta sum."""
    q = 12
    coprimes = [1, 5, 7, 11]
    s = mp.mpc(s_complex)
    total = mp.mpc(0)
    for a in coprimes:
        total += chi[a] * mp.zeta(s, mp.mpf(a) / q)
    return total * mp.power(q, -s)


def find_L_zeros(chi, T_max=40, max_zeros=30):
    """Find imaginary parts of L(1/2+it, chi) zeros for 0 < t < T_max."""
    def f(t):
        z = L_chi_mod_12(mp.mpc('0.5', float(t)), chi)
        return float(z.real), float(z.imag)

    zeros = []
    t = 0.5
    step = 0.5
    prev_re, prev_im = f(t)
    while t < T_max and len(zeros) < max_zeros:
        t += step
        re, im = f(t)
        # zero of L(1/2+it, chi) is t* such that both Re,Im pass through 0.
        # For real characters, L(1/2+it) is conjugate of L(1/2-it), but on the
        # critical line, function values are complex. Track sign change of
        # the squared magnitude's local minima.
        if re * prev_re < 0:
            try:
                root = mp.findroot(
                    lambda u: float(L_chi_mod_12(mp.mpc('0.5', float(u)),
                                                 chi).real),
                    [t - step, t], solver='anderson'
                )
                t_root = float(root)
                # Verify the imaginary part is also small there
                _, im_at_root = f(t_root)
                if abs(im_at_root) < 0.1 and (not zeros or t_root > zeros[-1] + 0.05):
                    zeros.append(t_root)
            except Exception:
                pass
        prev_re, prev_im = re, im
    return zeros


def truncated_pi_AP(x, a, char_idx, char_table, zero_lists):
    """
    Truncated explicit formula for pi(x; q, a):
      pi(x;q,a) ≈ Li(x)/phi(q) - (1/phi(q)) sum_chi chibar(a) * S_chi(x)
    S_chi(x) = sum over zeros 1/2 + i*gamma_k:
                Li(x^{1/2 + i gamma_k}) + Li(x^{1/2 - i gamma_k})
    (For real characters only — non-real chars need conjugation but mod 12
     all chars are real.)
    """
    phi_q = 4
    smooth = Li(x) / phi_q
    osc = mp.mpc(0)
    for ci in range(4):
        chi_t = char_table[ci]
        chi_a = chi_t[a]  # real for chars mod 12
        gammas = zero_lists[ci]
        S = mp.mpc(0)
        for g in gammas:
            rho = mp.mpc(0.5, g)
            S += mp.li(mp.power(x, rho))
            S += mp.li(mp.power(x, mp.mpc(0.5, -g)))
        osc += chi_a * S
    val = smooth - osc / phi_q
    return float(val.real)


def main():
    chars, coprimes = char_mod_12_table()
    print("# Truncated explicit formula on Dirichlet APs, q = 12, phi(q) = 4")

    # Load principal-character zeros (= zeta zeros)
    zeta_zeros = []
    with open("/apps/aplikacijos/prime-research/data/zeta_zeros_2000.txt") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                zeta_zeros.append(float(line))
    print(f"# Loaded {len(zeta_zeros)} zeta zeros (for principal chi_0)")

    # For non-principal characters, find zeros up to T_max=40
    print("# Finding zeros of non-principal L(s, chi mod 12) ...")
    L_zeros_nonprinc = []
    for ci in range(1, 4):
        zs = find_L_zeros(chars[ci], T_max=40, max_zeros=20)
        print(f"#   character {ci}: found {len(zs)} zeros, "
              f"first 5 = {[round(z,2) for z in zs[:5]]}")
        L_zeros_nonprinc.append(zs)

    # Build per-T zero-list slices
    Ts = [5, 10, 20, 40]
    xs = [1000, 3000, 5000, 7000, 10000]

    print("\n## Approximate vs true pi(x; 12, a)")
    print("# x      a    T    approx        truth   |diff|  exact?")

    n_total = 0
    n_exact = {T: 0 for T in Ts}
    for T in Ts:
        zeta_T = [g for g in zeta_zeros if g < T][:T] if T < len(zeta_zeros) else zeta_zeros[:T]
        non_T = [[g for g in zs if g < T][:T] for zs in L_zeros_nonprinc]
        zero_lists = [zeta_T] + non_T
        for x in xs:
            for a in coprimes:
                approx = truncated_pi_AP(x, a, None, chars, zero_lists)
                truth = true_pi_AP(x, 12, a)
                rounded = round(approx)
                diff = abs(approx - truth)
                exact = (rounded == truth)
                if exact:
                    n_exact[T] += 1
                n_total_T = len(xs) * len(coprimes)
                print(f"  x={x:>5} a={a:>2} T={T:>2}  approx={approx:>9.3f}  "
                      f"truth={truth:>5}  |diff|={diff:>7.3f}  "
                      f"exact={exact}")
        print()

    print("\n## Summary by T:")
    n_total_per_T = len(xs) * len(coprimes)
    for T in Ts:
        print(f"  T={T:>2}: {n_exact[T]}/{n_total_per_T} exact "
              f"({n_exact[T]/n_total_per_T*100:.1f}%)")

    print("\n## VERDICT")
    if n_exact[20] == n_total_per_T:
        print("VERDICT: T=20 zeros suffice for ALL (x, a) — PROPOSAL 4 PROMISING.")
    elif n_exact[40] == n_total_per_T:
        print("VERDICT: T=40 needed — promising but T must be measured vs x scaling.")
    else:
        print(f"VERDICT: T=40 still gives {n_exact[40]}/{n_total_per_T} exact — "
              "PROPOSAL 4 NEEDS more zeros than polylog at this x range. FAILS.")
    print("\n# DONE")


if __name__ == "__main__":
    main()
