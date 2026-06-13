"""C9.a — Divisor-only restriction of the Mertens-Ramanujan spike intensity.

Build M_W^{div}(n) := sum_{q | W, q sqf} mu(gcd(q,n)) / phi(q/gcd(q,n))
and verify the closed-form identity:

    M_W^{div}(n)  =  [gcd(n, W) = 1] · W / phi(W)         for W squarefree
    M_W^{div}(n)  =  [gcd(n, rad(W)) = 1] · rad(W)/phi(rad(W))  generally

Plus comparison to full T_Q (S191), Pearson with chi_P, L^2 deviation.

Run:
    python3 spike_divisor_only.py --N 100000
    python3 spike_divisor_only.py --N 1000000  --no-fullcompare  # skip full T_W
"""

from __future__ import annotations

import argparse
import json
import math
import time
from functools import lru_cache
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# arithmetic primitives
# ---------------------------------------------------------------------------

def sieve_chi_p(N: int) -> np.ndarray:
    """Return chi_P[1..N] as a {0,1} array indexed 1..N (length N+1)."""
    s = np.ones(N + 1, dtype=np.int8)
    s[0] = 0
    s[1] = 0
    for p in range(2, int(math.isqrt(N)) + 1):
        if s[p]:
            s[p * p :: p] = 0
    return s


def squarefree_divisors(W: int) -> list[int]:
    primes = []
    n = W
    p = 2
    while p * p <= n:
        if n % p == 0:
            primes.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        primes.append(n)
    out = [1]
    for p in primes:
        out = out + [d * p for d in out]
    return out


def radical(W: int) -> int:
    r = 1
    n = W
    p = 2
    while p * p <= n:
        if n % p == 0:
            r *= p
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        r *= n
    return r


@lru_cache(maxsize=None)
def euler_phi(n: int) -> int:
    if n == 1:
        return 1
    result = n
    nn = n
    p = 2
    while p * p <= nn:
        if nn % p == 0:
            while nn % p == 0:
                nn //= p
            result -= result // p
        p += 1
    if nn > 1:
        result -= result // nn
    return result


@lru_cache(maxsize=None)
def mobius(n: int) -> int:
    if n == 1:
        return 1
    nn = n
    primes_seen = 0
    p = 2
    while p * p <= nn:
        if nn % p == 0:
            nn //= p
            if nn % p == 0:
                return 0
            primes_seen += 1
        p += 1
    if nn > 1:
        primes_seen += 1
    return -1 if primes_seen % 2 else 1


# ---------------------------------------------------------------------------
# divisor-only spike intensity M_W^{div}
# ---------------------------------------------------------------------------

def m_div_W_array(W: int, N: int) -> np.ndarray:
    """M_W^{div}(n) for n=1..N, computed from the explicit definition."""
    sqf_divs = squarefree_divisors(W)
    out = np.zeros(N + 1, dtype=np.float64)
    n_arr = np.arange(N + 1, dtype=np.int64)
    for q in sqf_divs:
        # gcd(q, n) for each n
        gcd_vals = np.gcd(n_arr, q)
        # weight per gcd value
        for d in squarefree_divisors(q):
            mask = gcd_vals == d
            if not mask.any():
                continue
            qd = q // d
            term = mobius(d) / euler_phi(qd) if qd >= 1 else 0.0
            out[mask] += term
    return out[1:]  # n=1..N


def m_div_W_closed_form(W: int, N: int) -> np.ndarray:
    """[gcd(n, rad(W)) = 1] · rad(W) / phi(rad(W)), for n=1..N."""
    rW = radical(W)
    factor = rW / euler_phi(rW)
    n_arr = np.arange(1, N + 1, dtype=np.int64)
    coprime = np.gcd(n_arr, rW) == 1
    return coprime.astype(np.float64) * factor


# ---------------------------------------------------------------------------
# full M_Q (S191), for comparison: sum over ALL squarefree q <= Q
# ---------------------------------------------------------------------------

def squarefree_upto(Q: int) -> list[int]:
    out = [1]
    is_sqf = np.ones(Q + 1, dtype=bool)
    is_sqf[0] = False
    for p in range(2, int(math.isqrt(Q)) + 1):
        is_sqf[p * p :: p * p] = False
    for q in range(2, Q + 1):
        if is_sqf[q]:
            out.append(q)
    return out


def m_full_Q_array(Q: int, N: int) -> np.ndarray:
    """Full M_Q(n) (S191): sum over all squarefree q <= Q.  n=1..N."""
    sqf = squarefree_upto(Q)
    out = np.zeros(N + 1, dtype=np.float64)
    n_arr = np.arange(N + 1, dtype=np.int64)
    for q in sqf:
        if q == 1:
            out += 1.0
            continue
        gcd_vals = np.gcd(n_arr, q)
        # for each squarefree divisor d of q (d = gcd is automatically sqf since q is sqf)
        # but we only iterate over distinct values that actually appear
        unique_d = np.unique(gcd_vals)
        for d in unique_d:
            d = int(d)
            qd = q // d
            term = mobius(d) / euler_phi(qd) if qd >= 1 else 0.0
            out[gcd_vals == d] += term
    return out[1:]


# ---------------------------------------------------------------------------
# main experiment
# ---------------------------------------------------------------------------

def run(N: int, full_compare: bool) -> dict:
    print(f"N = {N}")
    t0 = time.time()
    chi_p = sieve_chi_p(N)[1:]  # length N
    n_primes = int(chi_p.sum())
    rho = n_primes / N
    print(f"  pi(N) = {n_primes}  rho = {rho:.6f}  (sieve {time.time() - t0:.2f}s)")

    out: dict = {"N": N, "pi(N)": n_primes, "rho": rho, "cells": []}

    primorials = [(2, 2), (6, "2·3"), (30, "2·3·5"), (210, "2·3·5·7"), (2310, "2·3·5·7·11")]
    sqf_nonprimorial = [(15, "3·5"), (105, "3·5·7")]
    nonsqf = [(12, "2²·3 ; rad=6"), (60, "2²·3·5 ; rad=30"), (420, "2²·3·5·7 ; rad=210")]

    all_W = primorials + sqf_nonprimorial + nonsqf

    for W, factorisation in all_W:
        cell: dict = {"W": W, "factorisation": factorisation}
        rW = radical(W)
        phi_rW = euler_phi(rW)
        ratio = rW / phi_rW
        cell["rad(W)"] = rW
        cell["phi(rad(W))"] = phi_rW
        cell["W/phi(rad(W))_factor"] = ratio
        cell["squarefree"] = (W == rW)

        t1 = time.time()
        m_emp = m_div_W_array(W, N)
        m_pred = m_div_W_closed_form(W, N)
        t_div = time.time() - t1

        # F1 / F2 / F3: pointwise identity
        abs_err = np.max(np.abs(m_emp - m_pred))
        cell["pointwise_max_abs_err"] = float(abs_err)
        cell["pointwise_identity_passes"] = bool(abs_err < 1e-12)

        # F4: mean = 1
        mean_M = float(m_emp.mean())
        cell["mean_M_div"] = mean_M
        cell["mean_unit_passes"] = bool(abs(mean_M - 1.0) < 1e-3)

        # F5: L^2 deviation closed form
        T_div = rho * m_emp
        n_arr = np.arange(1, N + 1, dtype=np.int64)
        coprime_mask = np.gcd(n_arr, rW) == 1
        # variance closed form ignoring boundary primes that divide W:
        # composites coprime: T = rho · ratio, chi = 0
        # primes coprime:     T = rho · ratio, chi = 1
        # not coprime:        T = 0,            chi = 1 if n is one of the small primes p|W and p<=N, else 0
        l2_emp = float(np.mean((T_div - chi_p) ** 2))
        # closed-form prediction (Bernoulli-on-coprime, ignoring p|W primes):
        phi_w = phi_rW / rW
        l2_pred_main = rho * (1.0 - rho / phi_w)
        # boundary correction: primes p | W contribute (chi_p=1, T=0) at those n
        boundary_primes = [p for p in [2, 3, 5, 7, 11, 13] if p <= N and rW % p == 0]
        n_boundary = len(boundary_primes)
        l2_pred_corrected = (rho * (1.0 - rho / phi_w) * 1.0 - 0.0)
        # the main-term predicts assuming primes p|W are negligible mass; for the
        # exact mean-square use the explicit count:
        # E[(T-chi)^2] = (1/N) [ |coprime, prime| (T_v - 1)^2 + |coprime, comp| T_v^2 + |non-coprime, prime: p|W| (0-1)^2 + |non-coprime, comp| 0 ]
        T_val = rho * ratio
        n_coprime = int(coprime_mask.sum())
        n_coprime_primes = int(chi_p[coprime_mask].sum())
        n_coprime_comp = n_coprime - n_coprime_primes
        n_noncoprime = N - n_coprime
        n_noncoprime_primes = n_primes - n_coprime_primes
        l2_explicit = (
            n_coprime_primes * (T_val - 1.0) ** 2
            + n_coprime_comp * T_val ** 2
            + n_noncoprime_primes * 1.0  # (0-1)^2 = 1
        ) / N
        cell["l2_emp"] = l2_emp
        cell["l2_pred_main"] = l2_pred_main
        cell["l2_pred_explicit"] = l2_explicit
        cell["l2_rel_err_explicit"] = float(abs(l2_emp - l2_explicit) / l2_explicit) if l2_explicit > 0 else 0.0
        cell["l2_rel_err_main"] = float(abs(l2_emp - l2_pred_main) / l2_pred_main) if l2_pred_main > 0 else 0.0
        cell["n_boundary_primes_in_W"] = n_noncoprime_primes
        cell["l2_passes"] = bool(cell["l2_rel_err_explicit"] < 1e-9)

        # F6: Pearson(chi_P, T_div) — should be effectively zero (T_div constant on coprimes)
        cov_div = float(np.cov(chi_p.astype(np.float64), T_div, ddof=0)[0, 1])
        var_chi = float(np.var(chi_p.astype(np.float64)))
        var_T = float(np.var(T_div))
        if var_T > 0 and var_chi > 0:
            pearson_div = cov_div / math.sqrt(var_chi * var_T)
        else:
            pearson_div = 0.0
        cell["pearson_chi_T_div"] = pearson_div

        # full T_W comparison (only for tractable sizes; full T_W cost is O(W · omega(n)·N) ~ O(W·N))
        if full_compare and W <= 2310:
            t2 = time.time()
            m_full = m_full_Q_array(W, N)
            T_full = rho * m_full
            t_full = time.time() - t2

            cov_full = float(np.cov(chi_p.astype(np.float64), T_full, ddof=0)[0, 1])
            var_T_full = float(np.var(T_full))
            pearson_full = cov_full / math.sqrt(var_chi * var_T_full) if var_T_full > 0 else 0.0
            cell["pearson_chi_T_full"] = pearson_full
            cell["full_minus_div_pearson_gap"] = pearson_full - pearson_div
            cell["F6_passes"] = bool(pearson_full > pearson_div + 1e-6)
            cell["t_full_sec"] = round(t_full, 2)
            l2_full = float(np.mean((T_full - chi_p) ** 2))
            cell["l2_T_full_minus_chi_p"] = l2_full
        else:
            cell["pearson_chi_T_full"] = None
            cell["F6_passes"] = None

        cell["t_div_sec"] = round(t_div, 2)
        out["cells"].append(cell)

        passes_summary = {
            "F1/F2/F3 pointwise": cell["pointwise_identity_passes"],
            "F4 mean": cell["mean_unit_passes"],
            "F5 L²": cell["l2_passes"],
            "F6 Pearson_full > Pearson_div": cell["F6_passes"],
        }
        print(f"  W={W} ({factorisation}): rad={rW}, ratio={ratio:.4f}")
        print(f"    abs_err={abs_err:.2e}  mean={mean_M:.6f}  pearson_div={pearson_div:+.4e}")
        if cell["pearson_chi_T_full"] is not None:
            print(f"    pearson_full={cell['pearson_chi_T_full']:+.4e}  l2_emp={l2_emp:.4e}  l2_pred={cell['l2_pred_explicit']:.4e}")
        else:
            print(f"    pearson_full=skipped  l2_emp={l2_emp:.4e}  l2_pred={cell['l2_pred_explicit']:.4e}")
        print(f"    {passes_summary}")

    # global pass/fail
    all_pass_F1_F4 = all(c["pointwise_identity_passes"] and c["mean_unit_passes"] for c in out["cells"])
    all_pass_F5 = all(c["l2_passes"] for c in out["cells"])
    all_pass_F6 = all(c["F6_passes"] for c in out["cells"] if c["F6_passes"] is not None)
    out["F1_F4_all_pass"] = all_pass_F1_F4
    out["F5_all_pass"] = all_pass_F5
    out["F6_all_pass"] = all_pass_F6

    print(
        f"\n=== Verdict: F1/F2/F3/F4 pass = {all_pass_F1_F4} ; F5 pass = {all_pass_F5} ; F6 pass = {all_pass_F6} ==="
    )
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=100_000)
    parser.add_argument("--no-fullcompare", action="store_true",
                        help="skip the full T_W comparison (faster)")
    parser.add_argument("--out", type=str, default="spike_divisor_only_results.json")
    args = parser.parse_args()
    out = run(N=args.N, full_compare=not args.no_fullcompare)
    out_path = Path(__file__).parent / args.out
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nresults written to {out_path}")


if __name__ == "__main__":
    main()
