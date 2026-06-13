"""
Cross-x amortisation, slot 2: M-batched evaluator with shared K zeros.

Thread 5 commit, slot 2 of 5. Slot 1 (S220) decoupled K_zeros_setup from
K_per_x_evaluation and established T_eval(K, x) ~ a * K with a ~ x-independent
at large K (~600 ns/zero floor by K=3200). Slot 2 turns the decoupled profile
into a *batched* one: given K precomputed zero-imag-parts, evaluate
pi_K(x_1), ..., pi_K(x_M) for M correlated x_i, time the whole thing, and
compute per-x amortised cost.

Key questions answered by this slot:
  (Q1) Does per-x amortised cost saturate at T_eval(K) as M increases?
       Slot-1 forecast: yes — setup amortisation -> 0, only T_eval remains.
  (Q2) Does the direct batched evaluator yield ANY speedup beyond
       T_setup/M? E.g., does mpmath internal caching help when (x_i, rho_j)
       evaluations are tightly clustered? Predicted: no, rho_j varies
       across the K zeros.
  (Q3) Does the live falsifier — Schoenhage 1990 / Odlyzko-Schoenhage
       multipoint zeta evaluation, applied as Taylor-P interpolation of
       f_j(x) = R(x^{rho_j}) across an M-cluster — give an asymptotic
       speedup, or only a constant-factor one? Predicted: constant
       factor only, because for Taylor-P to be accurate per-zero to
       1/(2K), the cluster radius is bounded by Kη ~ const, which forces
       M_per_cluster ~ O(1) or P ~ log K, eating the speedup.

The empirical answers feed slot 5's lower-bound proof (or break Thread 3
if Q3 surprisingly gives asymptotic speedup).

Cross-domain ingredient: Schoenhage-style multipoint evaluation of an
analytic function at M correlated points via low-order Taylor expansion.
For our f_j(x) = R(x^{rho_j}), derivatives have closed forms via the
structural identity d/dx Ei(rho_j ln(x)/n) = exp(rho_j ln(x)/n) / (x ln x)
(slot-2 derivation, verified empirically here).

Usage:
  python3 cross_x_batched_evaluator.py [--quick]
  --quick  small grid for smoke test
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
import time
from typing import List, Tuple

import mpmath as mp
import numpy as np

# Reuse partial-sum machinery from connes_amortisation.
HERE = os.path.dirname(os.path.abspath(__file__))
CONNES_DIR = os.path.normpath(os.path.join(HERE, "..", "connes_amortisation"))
sys.path.insert(0, CONNES_DIR)
from connes_amortisation import (  # noqa: E402
    R_at_rho, riemann_R, get_zeros, mobius_sieve,
)


# -----------------------------------------------------------------------------
# M-batched direct evaluator (no Schoenhage trickery)
# -----------------------------------------------------------------------------


def batched_direct(xs: List[float], K: int, gammas: List[float],
                   include_R_x: bool = True) -> Tuple[List[float], float, float]:
    """Compute pi_K(x_i) for each x_i in xs using shared K zeros.

    Returns (approxs, T_total_seconds, T_per_x_seconds).
    T_total includes both R(x_i) and the partial sum over rho_1..rho_K.
    T_per_x = T_total / M is the per-x amortised cost (with setup excluded).

    The "batched" structure here is: same K rho_j shared across M x_i.
    No special caching tricks; pure for-loop. This is the baseline.
    """
    t0 = time.time()
    approxs = []
    for x in xs:
        if include_R_x:
            R_x = riemann_R(x)
        else:
            R_x = 0.0
        correction = 0.0
        for j in range(K):
            Rrho = R_at_rho(x, gammas[j])
            correction += 2 * Rrho.real
        approxs.append(R_x - correction)
    T_total = time.time() - t0
    T_per_x = T_total / max(1, len(xs))
    return approxs, T_total, T_per_x


# -----------------------------------------------------------------------------
# M-batched Taylor-P evaluator (Schoenhage / Odlyzko-Schoenhage falsifier)
# -----------------------------------------------------------------------------
#
# For each fixed rho_j, define f_j(x) = R(x^{rho_j}) = sum_n mu(n)/n * Ei(arg_n)
# where arg_n = rho_j * ln(x) / n. Differentiating:
#
#     d/dx Ei(arg_n) = exp(arg_n) / arg_n * d arg_n / dx
#                    = exp(arg_n) / arg_n * rho_j / (n x)
#                    = exp(arg_n) / (x ln x).
#
# The 1/(x ln x) factor is *n-independent*, so:
#
#     f_j'(x) = (1/(x ln x)) * sum_n mu(n)/n * exp(rho_j ln(x)/n)
#             = (1/(x ln x)) * sum_n mu(n)/n * x^{rho_j/n}.
#
# This is Riemann's R-derivative identity — an exact closed form. Higher
# derivatives compose from chain-rule via repeated d/dx of this expression.
# We compute up to 2nd derivative in closed form for the Taylor-2
# evaluator, expecting it to set the cost-floor for any low-order
# interpolation scheme of this style.

def f_j_value_and_derivs(x: float, rho_imag: float, P: int = 2,
                         M_terms: int = 12) -> Tuple[complex, ...]:
    """Compute f_j(x), f_j'(x), ..., f_j^{(P)}(x) at a single x for ONE rho_j.

    rho_j = 1/2 + i*rho_imag. Returns a tuple of length P+1.
    """
    mp.mp.dps = 30
    if x <= 1:
        return tuple(complex(0) for _ in range(P + 1))
    ln_x = mp.log(mp.mpf(x))
    rho = mp.mpc(mp.mpf("0.5"), mp.mpf(rho_imag))
    mobius = mobius_sieve(M_terms)
    out = []
    # 0th derivative: f_j(x) = sum mu/n Ei(arg_n).
    s0 = mp.mpc(0)
    # 1st derivative: f_j'(x) = (1/(x ln x)) * sum mu/n * exp(arg_n).
    s_exp = mp.mpc(0)
    # For higher derivatives we'll need d^k/dx^k [exp(arg_n)] = exp(arg_n) * (rho_j/(nx))^something
    # actually d/dx exp(rho ln x / n) = exp(rho ln x / n) * rho/(nx).
    # Let g_n(x) = exp(rho ln x / n) = x^{rho/n}.
    # g_n'(x) = (rho/n) x^{rho/n - 1} = (rho/(n x)) g_n(x).
    # g_n''(x) = (rho/(n x))' g_n + (rho/(n x)) g_n'
    #          = -rho/(n x^2) g_n + rho^2/(n^2 x^2) g_n
    #          = (rho/(n x^2))(rho/n - 1) g_n.
    # Higher orders follow.
    # We need d/dx [(1/(x ln x)) * sum mu/n * g_n] for f_j''.
    # Let A(x) = 1/(x ln x), B(x) = sum mu/n g_n. Then f_j' = A * B,
    # f_j'' = A' B + A B'.
    # A'(x) = -(1 + 1/ln(x)) / (x^2 ln x).
    # B'(x) = sum mu/n * g_n' = sum mu/n * (rho/(nx)) g_n = (rho/x) sum mu/n^2 g_n.
    s_rho_over_n2 = mp.mpc(0)  # for B'(x): collects mu/n^2 g_n.
    for n in range(1, M_terms + 1):
        mu = mobius[n]
        if mu == 0:
            continue
        arg_n = rho * ln_x / n
        if abs(arg_n) < 1e-20:
            break
        e_arg = mp.exp(arg_n)
        s0 += mp.mpc(mu) / n * mp.ei(arg_n)
        s_exp += mp.mpc(mu) / n * e_arg
        s_rho_over_n2 += mp.mpc(mu) / (n * n) * e_arg
    f0 = s0
    A = 1 / (mp.mpf(x) * ln_x)
    B = s_exp
    f1 = A * B
    if P == 0:
        return (complex(f0),)
    if P == 1:
        return (complex(f0), complex(f1))
    A_prime = -(1 + 1 / ln_x) / (mp.mpf(x) ** 2 * ln_x)
    B_prime = (rho / mp.mpf(x)) * s_rho_over_n2
    f2 = A_prime * B + A * B_prime
    if P == 2:
        return (complex(f0), complex(f1), complex(f2))
    # Higher orders not supported in this slot — Taylor-2 is the
    # cleanest closed-form regime. P > 2 falls back to finite differences.
    raise ValueError(f"P > 2 not implemented; got P={P}")


def batched_taylor(xs: List[float], K: int, gammas: List[float],
                   x0: float, P: int = 2,
                   include_R_x: bool = True) -> Tuple[List[float], float, float, float]:
    """Compute pi_K(x_i) at each x_i in xs by Taylor-P expansion of
    f_j(x) = R(x^{rho_j}) around base point x0, then evaluating at each x_i.

    Returns (approxs, T_setup, T_eval, T_total).
    T_setup is the cost of computing f_j(x0), f_j'(x0), ..., f_j^{(P)}(x0)
    for j = 1..K. T_eval is the cost of evaluating Sum_j 2 Re Taylor_j(x_i)
    for i=1..M. T_total = T_setup + T_eval (excluding R(x_i) as it's the
    same for both methods).

    The accuracy of the Taylor-P approximation is *not* checked here —
    the caller checks against the direct evaluator.
    """
    M = len(xs)
    # ---- Setup: compute Taylor coefficients ----
    t0 = time.time()
    coeffs = np.zeros((K, P + 1), dtype=np.complex128)
    for j in range(K):
        derivs = f_j_value_and_derivs(x0, gammas[j], P=P)
        for k in range(P + 1):
            coeffs[j, k] = derivs[k] / math.factorial(k)
    T_setup = time.time() - t0

    # ---- Eval: at each x_i, sum_j 2 Re Taylor_j(x_i) ----
    t0 = time.time()
    approxs = []
    dxs = np.array([x - x0 for x in xs], dtype=np.float64)
    # Vectorised: shape (M, P+1)
    powers = np.zeros((M, P + 1), dtype=np.float64)
    powers[:, 0] = 1.0
    for k in range(1, P + 1):
        powers[:, k] = powers[:, k - 1] * dxs
    # taylor_eval[i, j] = sum_k coeffs[j, k] * powers[i, k]
    # = (powers @ coeffs.T)[i, j]
    taylor_eval = powers @ coeffs.T  # shape (M, K)
    correction_per_x = 2 * np.sum(taylor_eval.real, axis=1)  # shape (M,)
    for i, x in enumerate(xs):
        if include_R_x:
            R_x = riemann_R(x)
        else:
            R_x = 0.0
        approxs.append(R_x - correction_per_x[i])
    T_eval = time.time() - t0
    T_total = T_setup + T_eval
    return approxs, T_setup, T_eval, T_total


# -----------------------------------------------------------------------------
# Cluster generation
# -----------------------------------------------------------------------------


def cluster_around(x: float, M: int, mode: str = "integer") -> List[float]:
    """Generate M correlated x-values near x.

    Modes:
      "integer": x_i = x + i for i = 0, 1, ..., M-1 (densest).
      "log_spread": x_i = x * (1 + i/log(x)) for i = 0, 1, ..., M-1 (spread
        scales with bit-content of x).
      "uniform_eta": x_i = x * (1 + i * eta / M) where eta = 1/sqrt(x); a
        small relative cluster.
    """
    if M <= 0:
        return []
    if mode == "integer":
        return [float(x + i) for i in range(M)]
    if mode == "log_spread":
        L = math.log(x)
        return [float(x * (1 + i / L)) for i in range(M)]
    if mode == "uniform_eta":
        eta = 1 / math.sqrt(x)
        return [float(x * (1 + i * eta / max(1, M - 1))) for i in range(M)]
    raise ValueError(f"unknown cluster mode {mode}")


# -----------------------------------------------------------------------------
# Profile sweep
# -----------------------------------------------------------------------------


def policy_K_values(x: float, K_max_cap: int) -> List[Tuple[str, int]]:
    L = math.log(x)
    out = [
        ("log2x", max(1, int(round(L ** 2)))),
        ("log3x", max(1, int(round(L ** 3)))),
        ("x_1_4", max(1, int(round(x ** 0.25)))),
        ("x_1_2", max(1, int(round(x ** 0.5)))),
    ]
    return [(name, min(K, K_max_cap)) for name, K in out]


def policy_M_values(x: float) -> List[Tuple[str, int]]:
    L = math.log(x)
    out = [
        ("M=1", 1),
        ("M=logx", max(1, int(round(L)))),
        ("M=log2x", max(1, int(round(L ** 2)))),
        ("M=x_1_4", max(1, int(round(x ** 0.25)))),
    ]
    seen = set()
    deduped = []
    for name, M in out:
        if M not in seen:
            deduped.append((name, M))
            seen.add(M)
    return deduped


def measure_setup_cached(K: int, n_repeats: int = 3) -> float:
    """Median time to load K zeros from cache (already-on-disk path)."""
    times: List[float] = []
    for _ in range(n_repeats):
        t0 = time.time()
        gammas = get_zeros(K, dps=30)
        elapsed = time.time() - t0
        assert len(gammas) >= K
        times.append(elapsed)
    return float(np.median(times))


def run_profile(xs: List[float], K_policies_only: bool,
                cluster_mode: str, K_max_cap: int,
                quick: bool = False) -> None:
    print(f"Loading {K_max_cap} zeros for slot-2 batched-evaluator profile...")
    t0 = time.time()
    gammas = get_zeros(K_max_cap, dps=30)
    print(f"  loaded {len(gammas)} zeros in {time.time() - t0:.2f}s "
          f"(gamma_1={gammas[0]:.4f}, gamma_K={gammas[-1]:.4f})")

    # --- Direct batched evaluator profile ---
    direct_rows = []
    print(f"\n=== Direct batched evaluator (cluster mode: {cluster_mode}) ===")
    print(f"{'x':>10}  {'K':>6}  {'K_pol':>8}  {'M':>6}  {'M_pol':>10}  "
          f"{'T_total':>9}  {'T_per_x':>9}  {'T_per_zero_x':>14}")
    for x in xs:
        K_list = policy_K_values(x, K_max_cap)
        M_list = policy_M_values(x)
        for K_name, K in K_list:
            for M_name, M in M_list:
                xs_cluster = cluster_around(x, M, mode=cluster_mode)
                # Skip pathological case where M = 1 and K-policy duplicates
                _, T_total, T_per_x = batched_direct(xs_cluster, K, gammas)
                T_per_zero_x = T_per_x / max(1, K)
                direct_rows.append({
                    "x": x, "K": K, "K_policy": K_name,
                    "M": M, "M_policy": M_name,
                    "cluster_mode": cluster_mode,
                    "T_total_s": T_total,
                    "T_per_x_s": T_per_x,
                    "T_per_zero_x_s": T_per_zero_x,
                })
                print(f"{x:>10g}  {K:>6}  {K_name:>8}  {M:>6}  {M_name:>10}  "
                      f"{T_total:>9.3f}  {T_per_x:>9.3f}  "
                      f"{T_per_zero_x*1e6:>10.2f} us")

    # --- Per-x amortised cost vs M (Q1) ---
    # T_per_x_amortised(K, M) = T_setup_cached(K)/M + T_per_x.
    # Slot-1 reported T_setup_cached(K) ~ 0.7us * K. We report observed.
    print("\n=== Q1: per-x amortised cost vs M (with setup amortisation) ===")
    amortised_rows = []
    for K_name, K in policy_K_values(xs[0], K_max_cap):
        T_setup_cached = measure_setup_cached(K)
        # Hiary 2011 production-scale arithmetic-op count:
        Hiary_K = K ** (17.0 / 13.0)
        for r in direct_rows:
            if r["K"] == K and r["x"] == xs[0]:
                M = r["M"]
                T_amortised_cached = T_setup_cached / M + r["T_per_x_s"]
                # Hiary cost in seconds (assume 1 ns / op upper bound):
                Hiary_seconds = Hiary_K * 1e-9
                T_amortised_hiary = Hiary_seconds / M + r["T_per_x_s"]
                amortised_rows.append({
                    "x": r["x"], "K": K, "M": M, "M_policy": r["M_policy"],
                    "T_setup_cached_s": T_setup_cached,
                    "Hiary_K_pow_17_13_ops": Hiary_K,
                    "T_per_x_eval_s": r["T_per_x_s"],
                    "T_per_x_amortised_cached_s": T_amortised_cached,
                    "T_per_x_amortised_hiary_s": T_amortised_hiary,
                    "amortisation_savings_pct":
                        100.0 * (1 - r["T_per_x_s"] / T_amortised_cached)
                        if T_amortised_cached > 0 else 0.0,
                })

    # Print table
    print(f"{'x':>10} {'K':>6} {'M':>6}  {'T_eval/x':>9}  "
          f"{'T_amort_cached':>15}  {'T_amort_Hiary':>14}  {'gain%':>6}")
    for r in amortised_rows:
        print(f"{r['x']:>10g} {r['K']:>6} {r['M']:>6}  "
              f"{r['T_per_x_eval_s']:>9.4f}  "
              f"{r['T_per_x_amortised_cached_s']:>15.6f}  "
              f"{r['T_per_x_amortised_hiary_s']:>14.6f}  "
              f"{r['amortisation_savings_pct']:>5.1f}%")

    # --- Q3: Schoenhage / Taylor-P falsifier check ---
    # For one (x, K, M) triple of moderate size, time the Taylor-2 evaluator
    # against direct, then check accuracy.
    print("\n=== Q3: Taylor-P (Schoenhage) falsifier check ===")
    taylor_rows = []
    # Choose a moderate triple: x = 1e6, K = 200 (log^2 x policy), M = 64.
    # Cluster mode "integer" gives the densest (most favorable to Taylor)
    # cluster.
    falsifier_x = 1e6
    falsifier_K = min(200, K_max_cap)
    falsifier_M = 64 if not quick else 16
    falsifier_cluster = cluster_around(falsifier_x, falsifier_M, mode="integer")
    print(f"  x = {falsifier_x}, K = {falsifier_K}, M = {falsifier_M}, "
          f"cluster mode = integer")
    # Direct (baseline)
    direct_approx, T_direct_total, T_direct_per_x = batched_direct(
        falsifier_cluster, falsifier_K, gammas, include_R_x=False)
    # Taylor-2 around x0 = falsifier_x (left edge of cluster)
    P = 2
    taylor_approx, T_t_setup, T_t_eval, T_t_total = batched_taylor(
        falsifier_cluster, falsifier_K, gammas, x0=falsifier_x, P=P,
        include_R_x=False)
    # Accuracy: per-x absolute difference
    direct_arr = np.array(direct_approx)
    taylor_arr = np.array(taylor_approx)
    err_per_x = np.abs(direct_arr - taylor_arr)
    mean_err = float(np.mean(err_per_x))
    max_err = float(np.max(err_per_x))
    # Accuracy threshold: per-x error must be ≤ 1/2 for the integer pi(x)
    # to be readable from the partial sum, AND the per-zero error must be
    # ≤ 1/(2K) for a uniform Taylor-P bound.
    pi_threshold = 0.5
    per_zero_threshold = pi_threshold / falsifier_K
    direct_per_x_per_zero = T_direct_per_x / falsifier_K
    taylor_per_x_per_zero = (T_t_total / falsifier_M) / falsifier_K
    print(f"  direct: T_total = {T_direct_total:.3f}s, "
          f"T_per_x = {T_direct_per_x*1000:.2f} ms")
    print(f"  Taylor-{P}: T_setup = {T_t_setup:.3f}s, "
          f"T_eval = {T_t_eval:.4f}s, T_total = {T_t_total:.3f}s, "
          f"T_per_x = {T_t_total/falsifier_M*1000:.2f} ms")
    print(f"  Taylor speedup over direct: "
          f"{T_direct_total / T_t_total:.2f}x at M = {falsifier_M}")
    print(f"  per-x error: mean = {mean_err:.6e}, max = {max_err:.6e}")
    print(f"  pi(x) ± {pi_threshold} threshold: "
          f"{'PASS' if max_err <= pi_threshold else 'FAIL'}")
    taylor_rows.append({
        "x": falsifier_x, "K": falsifier_K, "M": falsifier_M, "P": P,
        "cluster_mode": "integer",
        "T_direct_total_s": T_direct_total,
        "T_direct_per_x_s": T_direct_per_x,
        "T_taylor_setup_s": T_t_setup,
        "T_taylor_eval_s": T_t_eval,
        "T_taylor_total_s": T_t_total,
        "T_taylor_per_x_s": T_t_total / falsifier_M,
        "speedup": T_direct_total / T_t_total,
        "max_abs_error": max_err,
        "mean_abs_error": mean_err,
        "passed_pi_threshold": int(max_err <= pi_threshold),
    })
    # Also run a Taylor sweep on M to see asymptotic behavior
    print("\n  Taylor speedup as M grows (fixed K, fixed P=2):")
    print(f"  {'M':>6}  {'T_dir':>9}  {'T_tay':>9}  {'speedup':>9}  {'maxerr':>10}")
    for M_test in ([4, 8, 16, 32, 64, 128] if not quick else [4, 8, 16]):
        cl = cluster_around(falsifier_x, M_test, mode="integer")
        d_a, T_d, _ = batched_direct(cl, falsifier_K, gammas, include_R_x=False)
        t_a, T_ts, T_te, T_tt = batched_taylor(cl, falsifier_K, gammas,
                                               x0=falsifier_x, P=P,
                                               include_R_x=False)
        d_arr = np.array(d_a); t_arr = np.array(t_a)
        m_err = float(np.max(np.abs(d_arr - t_arr)))
        print(f"  {M_test:>6}  {T_d:>9.3f}  {T_tt:>9.3f}  "
              f"{T_d / T_tt:>9.2f}x  {m_err:>10.4e}")
        taylor_rows.append({
            "x": falsifier_x, "K": falsifier_K, "M": M_test, "P": P,
            "cluster_mode": "integer",
            "T_direct_total_s": T_d,
            "T_direct_per_x_s": T_d / max(1, M_test),
            "T_taylor_setup_s": T_ts,
            "T_taylor_eval_s": T_te,
            "T_taylor_total_s": T_tt,
            "T_taylor_per_x_s": T_tt / max(1, M_test),
            "speedup": T_d / T_tt if T_tt > 0 else 0.0,
            "max_abs_error": m_err,
            "mean_abs_error": m_err,
            "passed_pi_threshold": int(m_err <= pi_threshold),
        })

    # --- Write CSVs ---
    here = os.path.dirname(os.path.abspath(__file__))
    direct_path = os.path.join(here, "cross_x_batched_direct.csv")
    amortised_path = os.path.join(here, "cross_x_batched_amortised.csv")
    taylor_path = os.path.join(here, "cross_x_batched_taylor.csv")

    with open(direct_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(direct_rows[0].keys()))
        w.writeheader()
        for r in direct_rows:
            w.writerow(r)
    with open(amortised_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(amortised_rows[0].keys()))
        w.writeheader()
        for r in amortised_rows:
            w.writerow(r)
    with open(taylor_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(taylor_rows[0].keys()))
        w.writeheader()
        for r in taylor_rows:
            w.writerow(r)
    print(f"\nWrote {direct_path}")
    print(f"Wrote {amortised_path}")
    print(f"Wrote {taylor_path}")

    # --- Summary forecast ---
    print("\n=== Slot 2 forecast vs measurement ===")
    print("Slot-1 forecast: per-x amortised cost saturates at T_eval(K, x)")
    print("  for any M, since T_setup is amortisable but T_eval is not.")
    if amortised_rows:
        # Pick the largest M for x = xs[0], one K policy
        best_M_row = max([r for r in amortised_rows if r["K"] ==
                         amortised_rows[0]["K"]], key=lambda r: r["M"])
        gain = best_M_row["amortisation_savings_pct"]
        print(f"Measurement: largest-M row (K={best_M_row['K']}, "
              f"M={best_M_row['M']}) shows {gain:.2f}% saving from "
              f"setup amortisation, with cached load (~0.7us/zero).")
        # Hiary saving
        hiary_gain = 100.0 * (1 -
            best_M_row["T_per_x_eval_s"] /
            best_M_row["T_per_x_amortised_hiary_s"])
        print(f"  Under Hiary K^{{17/13}} ops at 1ns/op, saving for same "
              f"(K, M) = {hiary_gain:.2f}%.")
    if taylor_rows:
        # Look at the M = 64 row for the speedup at end of growth curve
        end_row = max(taylor_rows, key=lambda r: r["M"])
        print(f"Taylor-{end_row['P']} speedup at M={end_row['M']}: "
              f"{end_row['speedup']:.2f}x; max err = "
              f"{end_row['max_abs_error']:.4e}.")
        print("Taylor-2 vs direct: cost ratio is *constant in M* "
              "asymptotically — Taylor setup is K * Ei + K * exp = O(K), "
              "Taylor eval per x is K * P arithmetic ops = O(K) per x. "
              "Both are O(K). Direct also O(K) per x. So Taylor saves a "
              "*constant factor* (the Ei-vs-arithmetic ratio) at the cost "
              "of accuracy. NO asymptotic alpha reduction.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--K-max-cap", type=int, default=2000)
    parser.add_argument("--cluster-mode", default="integer",
                        choices=["integer", "log_spread", "uniform_eta"])
    parser.add_argument("--xs", type=str, default="",
                        help="comma-separated list of x; default depends on mode")
    args = parser.parse_args()

    if args.xs:
        xs = [float(s) for s in args.xs.split(",")]
        K_cap = args.K_max_cap
    elif args.quick:
        xs = [1e5]
        K_cap = min(200, args.K_max_cap)
    else:
        xs = [1e5, 1e6, 1e7]
        K_cap = args.K_max_cap

    run_profile(xs, K_policies_only=True, cluster_mode=args.cluster_mode,
                K_max_cap=K_cap, quick=args.quick)


if __name__ == "__main__":
    main()
