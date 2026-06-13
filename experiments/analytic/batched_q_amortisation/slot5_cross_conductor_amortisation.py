"""
Thread 6, slot 5: cross-conductor (Q-batched) amortisation.

Mission (from .commit_state Thread 6 slot 5 PRIORITY (a) +
session230_commit_p1_composite_q_multi_axis_fft.md "next-action"):

  Slots 1-4 closed three axes of the Thread-6 amortisation problem:
    AXIS 1 (a-direction): trivial 8x at M=256 (zeros independent of a).
    AXIS 2 (chi-direction, fixed q): bounded constant ~1.7-2x via the
      cyclic-DFT primitive (slot 2 / slot 3).
    AXIS 3 (composite q): bounded constant ~1.75x via multi-axis FFT
      (slot 4).
  Slot 5 closes the FOURTH axis: cross-conductor batches.

  Question: given Q = {q_1, ..., q_M} of distinct prime conductors,
  can the per-(q, a) cost be made sub-linear in M? Equivalently, is
  there a structural primitive that shares zero-finding work across
  distinct conductors?

Theoretical structure:

  The slot-3 zero-finder pipeline at one prime conductor q decomposes:

    Stage A  cp_all = 1/n^{1/2 + i t}  for n = 1..N_q,  t in t-grid
    Stage B  W_q = aggregate cp_all by residue mod q in log-g order
             (per-conductor; q-dependent residue map)
    Stage C  L_main = phi(q) * ifft(W_q, axis=0)    (per-conductor;
             length-phi(q) FFT)
    Stage D  Reflected term: pointwise mult of L_main + W_chi factors
    Stage E  loggamma at t-grid (parity in {0,1}) -- conductor-independent!
    Stage F  Hardy Z = real(exp(i*theta) * L) -- depends on q via
             log(q/pi) only (cheap)
    Stage G  Sign-change bracket on Z; per-character lookup

  STAGES SHAREABLE ACROSS Q:
    Stage A (cp_all up to N_max = oversize * sqrt(q_max * t_max / (2pi)))
    Stage E (loggamma values: only depend on t-grid + parity)

  STAGES NON-SHAREABLE (per-conductor):
    Stage B, C (residue map and FFT length depend on q)
    Stage D (W_chi, q^{-it} factor depend on q)
    Stage F, G (per-conductor)

The slot-5 experiment:

  Build two end-to-end pipelines:

  (I) INDEPENDENT: run slot-3 zero finder once per q, sequentially.
      Cumulative cost = sum_i T_per_q(q_i).

  (II) Q-BATCHED-SHARED: share Stage A (cp_all up to N_max) and Stage E
      (loggamma table) across the family. Per-conductor work
      (Stages B, C, D, F, G) still runs per-q.

  Hypothesis: T_batched / T_independent = 1 - delta, where delta is a
  BOUNDED CONSTANT around 0.10-0.20 (saving Stage A + Stage E across
  M-1 conductors). NOT sub-linear in M.

What would falsify the bounded-constant prediction:

  (i)  T_batched / T_independent --> 0 as M --> oo for some cleverer
       sharing (would be a genuine partial-positive on the Q-axis).

  (ii) T_batched > T_independent (the shared Stage A is dominated by
       its larger N_max overhead and overflows the savings).

  (iii) Sharing Stage A actually breaks accuracy at the smaller q's
       (N_max is "too large" for q_min and breaks numerical
       conditioning of the smallest-q FFT).

Cross-domain ingredient:

  None new in slot 5. The slot-5 attack is rigorous-experimentation
  on the Q-axis using slot-3/4's existing primitives. The cross-domain
  technique register is unchanged. The slot-5 outcome IS the negative-
  shape closure of the Q-axis, which is the theoretical wrap of
  Thread 6.

References:
  - Slot 3 (S229): single-q FFT-shared pipeline. ~190 ms at q=1009, K=200.
  - Slot 4 (S230): composite-q multi-axis FFT. ~167 ms at q=1001, K=200.
  - S202 WRAP: the analogous cross-x amortisation analysis for Thread 3.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
import time
from typing import Dict, List, Tuple

import numpy as np
from scipy.special import loggamma

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from per_q_a_explicit_formula_profile import (  # noqa: E402
    euler_phi, primitive_root, get_chi0_zeros,
)
from slot2_afe_shared_l_eval import (  # noqa: E402
    coprimes_in_log_g_order, compute_complex_pow_all_n,
)
from slot3_zero_finder_dirichlet import (  # noqa: E402
    afe_truncation_oversized,
    compute_gauss_sums_fft,
    zero_brackets_vectorised,
)


# =====================================================================
# Stage A (shareable): cp_all = 1/n^{1/2+it}, n = 1..N, t in t_grid.
# Conductor-independent.
# =====================================================================


def build_shared_cp_all(t_grid: np.ndarray, N_max: int) -> np.ndarray:
    """Returns cp_all of shape (N_max, n_t)."""
    return compute_complex_pow_all_n(t_grid, N_max)


# =====================================================================
# Stage E (shareable): loggamma for parity 0/1, indexed by t-grid.
# =====================================================================


def build_shared_loggamma(t_grid: np.ndarray) -> Dict[str, np.ndarray]:
    """Pre-compute the loggamma values that all conductors will need.

    We need
      arg Gamma((1/2 + a + it)/2)  for a in {0, 1}
      Gamma((1/2 - it + a)/2) / Gamma((1/2 + it + a)/2)  for a in {0, 1}
        = exp(loggamma(z_-) - loggamma(z_+))

    Both depend only on t_grid + parity, NOT on q. So ONE pass per
    family.

    Returns dict with keys:
      arg_g_even  shape (n_t,)
      arg_g_odd   shape (n_t,)
      rg_even     shape (n_t,)
      rg_odd      shape (n_t,)
    """
    z_plus_even = (0.5 + 1j * t_grid) / 2
    z_plus_odd = (1.5 + 1j * t_grid) / 2
    z_minus_even = (0.5 - 1j * t_grid) / 2
    z_minus_odd = (1.5 - 1j * t_grid) / 2

    lg_p_even = loggamma(z_plus_even)
    lg_p_odd = loggamma(z_plus_odd)
    lg_m_even = loggamma(z_minus_even)
    lg_m_odd = loggamma(z_minus_odd)

    return {
        "arg_g_even": np.imag(lg_p_even),
        "arg_g_odd": np.imag(lg_p_odd),
        "rg_even": np.exp(lg_m_even - lg_p_even),
        "rg_odd": np.exp(lg_m_odd - lg_p_odd),
    }


# =====================================================================
# Per-conductor stages B, C, D, F, G with hooks for shared inputs.
# =====================================================================


def find_zeros_one_conductor_with_shared(
    q: int,
    t_grid: np.ndarray,
    N: int,
    cp_all_shared: np.ndarray,
    lg: Dict[str, np.ndarray],
    n_zeros_target: int | None = None,
    t_min: float = 2.0,
) -> Tuple[Dict[int, List[float]], Dict[str, float]]:
    """Run stages B, C, D, F, G for one conductor q using PRE-BUILT
    cp_all (Stage A) and PRE-BUILT loggamma table (Stage E).

    cp_all_shared has shape (N_max_shared, n_t) with N_max_shared >= N.
    We slice cp_all_shared[:N] for this conductor's AFE truncation.
    """
    t0_total = time.perf_counter()

    coprimes_logg, log_g = coprimes_in_log_g_order(q)
    phi = len(coprimes_logg)
    parities = (np.arange(phi) % 2).astype(np.int64)
    coprime_set = set(coprimes_logg)
    n_t = len(t_grid)

    # ----- Stage B: aggregate cp_all by residue mod q in log-g order -----
    t_B0 = time.perf_counter()
    W = np.zeros((phi, n_t), dtype=complex)
    for n_idx in range(N):
        r = (n_idx + 1) % q
        if r in coprime_set:
            W[log_g[r]] += cp_all_shared[n_idx]
    t_B = time.perf_counter() - t_B0

    # ----- Stage C: length-phi FFT over (Z/qZ)* -----
    t_C0 = time.perf_counter()
    M = phi * np.fft.ifft(W, axis=0)  # (phi, n_t)
    t_C = time.perf_counter() - t_C0

    # ----- Stage D: reflected term + Gauss sums + W_chi -----
    t_D0 = time.perf_counter()
    gauss_sums = compute_gauss_sums_fft(q)
    i_pow_a = np.where(parities == 0, 1.0 + 0j, 1j)
    W_chi = gauss_sums / (i_pow_a * math.sqrt(q))

    log_q_pi = math.log(q / math.pi)
    qp_factor = np.exp(-1j * log_q_pi * t_grid)

    parity_mask = parities[:, None].astype(bool)
    rg = np.where(parity_mask, lg["rg_odd"][None, :], lg["rg_even"][None, :])
    refl_factor = W_chi[:, None] * qp_factor[None, :] * rg
    L = M + refl_factor * np.conj(M)
    t_D = time.perf_counter() - t_D0

    # ----- Stage F: Hardy theta and Z -----
    t_F0 = time.perf_counter()
    arg_W = np.angle(W_chi)
    arg_gamma = np.where(
        parity_mask, lg["arg_g_odd"][None, :], lg["arg_g_even"][None, :]
    )
    theta = (t_grid[None, :] / 2.0) * log_q_pi + arg_gamma - arg_W[:, None] / 2.0
    Z = np.real(np.exp(1j * theta) * L)
    t_F = time.perf_counter() - t_F0

    # ----- Stage G: sign-change bracket -----
    t_G0 = time.perf_counter()
    zeros_db = zero_brackets_vectorised(Z, t_grid, t_min=t_min)
    K_principal = max(n_zeros_target if n_zeros_target else 1000, 100)
    zeta_zeros = get_chi0_zeros(K_principal)
    zeta_zeros_in_range = [g for g in zeta_zeros if t_min <= g <= t_grid[-1]]
    zeros_db[0] = zeta_zeros_in_range[:n_zeros_target] if n_zeros_target else zeta_zeros_in_range
    if n_zeros_target is not None:
        for j in range(1, phi):
            zeros_db[j] = zeros_db[j][:n_zeros_target]
    t_G = time.perf_counter() - t_G0

    t_total = time.perf_counter() - t0_total
    timing = {
        "t_total_s": t_total, "t_B_aggregate": t_B, "t_C_fft": t_C,
        "t_D_reflected": t_D, "t_F_hardy_z": t_F, "t_G_bracket": t_G,
        "phi": phi, "N_used": N,
    }
    return zeros_db, timing


# =====================================================================
# Two pipelines: INDEPENDENT vs Q-BATCHED-SHARED
# =====================================================================


def pipeline_independent(
    Q_list: List[int], t_max: float, n_t: int,
    n_zeros_target: int | None = None, t_min: float = 2.0,
    oversize: float = 12.0,
) -> Tuple[Dict[int, Dict[int, List[float]]], Dict[str, float]]:
    """Run slot-3 zero finder once per q, no cross-conductor sharing.
    Each q gets its own cp_all and its own loggamma table."""
    t0 = time.perf_counter()
    db_all: Dict[int, Dict[int, List[float]]] = {}
    cumulative = {
        "t_stage_A_setup_total": 0.0,
        "t_stage_E_setup_total": 0.0,
        "t_per_conductor_total": 0.0,
    }
    per_q_timings: List[Dict[str, float]] = []

    for q in Q_list:
        t_grid = np.linspace(max(t_min, 1.0), t_max, n_t)
        N = afe_truncation_oversized(q, t_max, oversize=oversize)

        tA0 = time.perf_counter()
        cp_all = build_shared_cp_all(t_grid, N)
        tA = time.perf_counter() - tA0
        cumulative["t_stage_A_setup_total"] += tA

        tE0 = time.perf_counter()
        lg = build_shared_loggamma(t_grid)
        tE = time.perf_counter() - tE0
        cumulative["t_stage_E_setup_total"] += tE

        tQ0 = time.perf_counter()
        zeros_db, timing = find_zeros_one_conductor_with_shared(
            q, t_grid, N, cp_all, lg,
            n_zeros_target=n_zeros_target, t_min=t_min,
        )
        tQ = time.perf_counter() - tQ0
        cumulative["t_per_conductor_total"] += tQ

        timing["q"] = q
        timing["N"] = N
        timing["t_stage_A"] = tA
        timing["t_stage_E"] = tE
        per_q_timings.append(timing)
        db_all[q] = zeros_db

    cumulative["t_total"] = time.perf_counter() - t0
    cumulative["per_q_timings"] = per_q_timings
    return db_all, cumulative


def pipeline_q_batched_shared(
    Q_list: List[int], t_max: float, n_t: int,
    n_zeros_target: int | None = None, t_min: float = 2.0,
    oversize: float = 12.0,
) -> Tuple[Dict[int, Dict[int, List[float]]], Dict[str, float]]:
    """Share Stage A (cp_all up to N_max for q_max) and Stage E
    (loggamma) across the family Q. Per-conductor stages still run."""
    t0 = time.perf_counter()
    q_max = max(Q_list)
    t_grid = np.linspace(max(t_min, 1.0), t_max, n_t)
    N_max = afe_truncation_oversized(q_max, t_max, oversize=oversize)

    tA0 = time.perf_counter()
    cp_all = build_shared_cp_all(t_grid, N_max)
    tA = time.perf_counter() - tA0

    tE0 = time.perf_counter()
    lg = build_shared_loggamma(t_grid)
    tE = time.perf_counter() - tE0

    db_all: Dict[int, Dict[int, List[float]]] = {}
    per_q_timings: List[Dict[str, float]] = []
    t_per_q_cumulative = 0.0

    for q in Q_list:
        N_q = afe_truncation_oversized(q, t_max, oversize=oversize)
        tQ0 = time.perf_counter()
        zeros_db, timing = find_zeros_one_conductor_with_shared(
            q, t_grid, N_q, cp_all, lg,
            n_zeros_target=n_zeros_target, t_min=t_min,
        )
        tQ = time.perf_counter() - tQ0
        t_per_q_cumulative += tQ
        timing["q"] = q
        timing["N"] = N_q
        per_q_timings.append(timing)
        db_all[q] = zeros_db

    cumulative = {
        "t_total": time.perf_counter() - t0,
        "t_stage_A_shared_once": tA,
        "t_stage_E_shared_once": tE,
        "t_per_conductor_total": t_per_q_cumulative,
        "per_q_timings": per_q_timings,
        "N_max": N_max,
    }
    return db_all, cumulative


# =====================================================================
# Cross-method correctness check
# =====================================================================


def verify_pipelines_agree(Q_list: List[int], t_max: float, n_t: int,
                           n_zeros_target: int = 10, t_min: float = 2.0) -> Dict[str, float]:
    """Run both pipelines on Q_list; check that they produce the SAME
    zero database. Stage-A sharing must not change accuracy."""
    db_indep, _ = pipeline_independent(
        Q_list, t_max, n_t,
        n_zeros_target=n_zeros_target, t_min=t_min,
    )
    db_batch, _ = pipeline_q_batched_shared(
        Q_list, t_max, n_t,
        n_zeros_target=n_zeros_target, t_min=t_min,
    )
    max_abs_diff = 0.0
    for q in Q_list:
        for j in db_indep[q]:
            if j == 0:
                continue
            zs_i = db_indep[q][j]
            zs_b = db_batch[q][j]
            n_min = min(len(zs_i), len(zs_b))
            for k in range(n_min):
                d = abs(zs_i[k] - zs_b[k])
                if d > max_abs_diff:
                    max_abs_diff = d
    return {"max_abs_diff": max_abs_diff}


# =====================================================================
# Main: run the M-sweep
# =====================================================================


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["sweep", "verify", "scaling"], default="sweep")
    parser.add_argument("--K", type=int, default=200)
    parser.add_argument("--n-t", type=int, default=823)
    parser.add_argument("--out", default=os.path.join(HERE, "slot5_cross_conductor.csv"))
    args = parser.parse_args()

    if args.mode == "verify":
        Q_small = [101, 251, 503, 1009]
        res = verify_pipelines_agree(Q_small, t_max=50.0, n_t=400, n_zeros_target=5)
        print(f"VERIFY: max abs diff between independent and Q-batched = {res['max_abs_diff']:.3e}")
        if res["max_abs_diff"] > 1e-10:
            print("  WARNING: pipelines disagree above 1e-10. Stage-A sharing breaks accuracy?")
        else:
            print("  OK. Q-batched-shared pipeline is numerically equivalent.")
        return

    # Default mode: M-sweep on Q = {1009, 2003, 5003, 10007}
    # Use first M elements for each M.
    Q_full = [1009, 2003, 5003, 10007]
    K = args.K
    # t_max chosen for K zeros at largest q: gamma_K ~ K * 2pi / log(q*K/(2pi e))
    # ~ 200 * 2pi / log(10007*200/(2pi*e)) ~ 130 + slack
    t_max_default = 200.0
    n_t = args.n_t

    rows = []
    headers = [
        "M", "Q", "K", "n_t", "t_max",
        "indep_total_s", "batch_total_s", "speedup",
        "indep_per_q_avg_ms", "batch_per_q_avg_ms",
        "indep_stageA_total_ms", "batch_stageA_once_ms",
        "indep_stageE_total_ms", "batch_stageE_once_ms",
    ]

    print("Slot 5: cross-conductor (Q-batched) amortisation sweep")
    print(f"Q_full = {Q_full}, K = {K}, n_t = {n_t}, t_max = {t_max_default}")
    print()

    # Warm-up the FFT planners by running one small instance
    _ = pipeline_independent([101], t_max=20.0, n_t=200, n_zeros_target=5)
    _ = pipeline_q_batched_shared([101], t_max=20.0, n_t=200, n_zeros_target=5)

    for M in [1, 2, 3, 4]:
        Q = Q_full[:M]
        # 3-run median for stability
        indep_times = []
        batch_times = []
        ind_meta = None
        batch_meta = None
        for run in range(3):
            db_i, ti = pipeline_independent(
                Q, t_max=t_max_default, n_t=n_t, n_zeros_target=K, t_min=2.0,
            )
            db_b, tb = pipeline_q_batched_shared(
                Q, t_max=t_max_default, n_t=n_t, n_zeros_target=K, t_min=2.0,
            )
            indep_times.append(ti["t_total"])
            batch_times.append(tb["t_total"])
            if run == 1:  # use median run for breakdown
                ind_meta = ti
                batch_meta = tb

        indep_med = float(np.median(indep_times))
        batch_med = float(np.median(batch_times))

        speedup = indep_med / batch_med
        row = {
            "M": M,
            "Q": ",".join(str(x) for x in Q),
            "K": K,
            "n_t": n_t,
            "t_max": t_max_default,
            "indep_total_s": indep_med,
            "batch_total_s": batch_med,
            "speedup": speedup,
            "indep_per_q_avg_ms": 1000 * indep_med / M,
            "batch_per_q_avg_ms": 1000 * batch_med / M,
            "indep_stageA_total_ms": 1000 * ind_meta["t_stage_A_setup_total"],
            "batch_stageA_once_ms": 1000 * batch_meta["t_stage_A_shared_once"],
            "indep_stageE_total_ms": 1000 * ind_meta["t_stage_E_setup_total"],
            "batch_stageE_once_ms": 1000 * batch_meta["t_stage_E_shared_once"],
        }
        rows.append(row)
        print(
            f"M={M:1d} | Q={row['Q']:25s} | "
            f"indep {1000*indep_med:7.1f} ms | "
            f"batch {1000*batch_med:7.1f} ms | "
            f"speedup {speedup:.3f}× | "
            f"per-q indep {row['indep_per_q_avg_ms']:6.1f} ms | "
            f"per-q batch {row['batch_per_q_avg_ms']:6.1f} ms"
        )

    # Save CSV
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=headers)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print()
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
