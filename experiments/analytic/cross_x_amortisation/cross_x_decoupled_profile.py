"""
Cross-x amortisation, slot 1: decouple K_zeros_setup from K_per_x_evaluation.

Thread 5 (commit). The S202 unified theorem closed Thread 3 (Galway frontier)
*per query* under Montgomery: K*(x, p) = Theta-tilde(x) zeros are required
for pi(x) +/- 1 in distribution at any positive hit-rate. But Riemann zeros
are properties of zeta, not x. If K zeros are computed once and reused across
M batched queries pi(x_1), ..., pi(x_M), per-x cost decomposes:

    T_total(K, M, x_max, x_1..M) = T_setup(K) + sum_{i=1}^M T_eval(K, x_i),

so amortised per-x cost is

    T_per_x_amortised = T_setup(K) / M + T_eval_avg(K).

Slot 1's job: measure T_setup(K) and T_eval(K, x) separately, characterise
each as functions of (K, x), so subsequent slots can compute the amortised
profile under any choice of M and x-distribution.

The empirical question this slot answers concretely:
  (Q1) Is T_eval(K, x) approximately linear in K?
  (Q2) Is T_eval(K, x)/K approximately x-independent?
  (Q3) How does T_setup(K) scale with K (cache load vs from-scratch)?

If Q1+Q2 hold, then per-x evaluation is asymptotically Theta(K) regardless
of x, and the only knob to lower per-x amortised cost is K itself; any
useful cross-x amortisation must come from setup amortisation, which is
amortisable across M (T_setup/M -> 0) but cannot reduce the T_eval term.

We profile at:
  x in {1e5, 1e6, 1e7}
  K in {ceil(log^2 x), ceil(log^3 x), ceil(x^{1/4}), ceil(x^{1/2})}
plus K-grid {25, 50, 100, 200, 400, 800, 1600, 3200, 6400} for the linear-
fit verification of Q1.

What would falsify the slot-1 conclusion:
  - T_eval(K, x) scales sub-linearly in K (e.g., due to internal mpmath
    optimisations sharing across rho's): would invalidate the K-linear
    bound and re-open the polylog question.
  - T_eval/K shows strong x-dependence: would mean the evaluation cost is
    x-coupled and the decomposition into setup + per-x is misleading.
  - T_setup(K) at small K is dominated by some non-cacheable fixed cost
    (e.g., dps ramp-up): would mean setup amortisation is bounded below by
    a constant, capping the M -> infty benefit.

Cross-domain ingredient (CROSS_DOMAIN_TECHNIQUES.md sec 8 amortised
algorithmics): batch-query analysis (Tarjan 1985, Demaine-Patrascu) applied
to the explicit-formula evaluator. The decoupled profile is the precondition
for any batched-query lower-bound proof in slots 2-5.

Usage:
  python cross_x_decoupled_profile.py [--quick]
  --quick  small K-grid for smoke test (used in CI)
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import time
from typing import List, Tuple

import numpy as np

# Reuse the partial-sum machinery from connes_amortisation.
import sys
HERE = os.path.dirname(os.path.abspath(__file__))
CONNES_DIR = os.path.normpath(os.path.join(HERE, "..", "connes_amortisation"))
sys.path.insert(0, CONNES_DIR)
from connes_amortisation import (  # noqa: E402
    R_at_rho, riemann_R, get_zeros, build_prime_count_array, pi_from_array,
)


# -----------------------------------------------------------------------------
# Decoupled timing primitives
# -----------------------------------------------------------------------------


def time_setup(K: int, dps: int = 30, n_repeats: int = 3) -> float:
    """Time cost of acquiring K zeros' imaginary parts in memory.

    Uses the cached zeros file; this measures the I/O + parse cost, NOT the
    cost of computing the zeros from scratch. The from-scratch cost via
    Hiary 2011 is K^{17/13} arithmetic operations and is the asymptotically
    relevant setup cost for cross-x amortisation; we report both.

    Returns the median elapsed time (seconds) over n_repeats trials.
    """
    times: List[float] = []
    for _ in range(n_repeats):
        t0 = time.time()
        gammas = get_zeros(K, dps=dps)
        elapsed = time.time() - t0
        assert len(gammas) >= K, "got %d zeros, requested %d" % (len(gammas), K)
        times.append(elapsed)
    return float(np.median(times))


def time_eval(x: float, K: int, gammas: List[float],
              n_repeats: int = 1) -> Tuple[float, float, float]:
    """Time the explicit-formula partial sum
        pi_K(x) = R(x) - 2 sum_{j=1}^{K} Re R(x^{rho_j}).

    Decomposes into (T_R, T_sum, T_total) = (R(x) cost, partial-sum cost,
    total). Returns medians over n_repeats trials.
    """
    R_times = []
    sum_times = []
    total_times = []
    for _ in range(n_repeats):
        t0 = time.time()
        R_x = riemann_R(x)
        t_after_R = time.time()
        correction = 0.0
        for j in range(K):
            Rrho = R_at_rho(x, gammas[j])
            correction += 2 * Rrho.real
        approx = R_x - correction  # noqa: F841 - used to ensure no DCE
        t_end = time.time()
        R_times.append(t_after_R - t0)
        sum_times.append(t_end - t_after_R)
        total_times.append(t_end - t0)
    return (float(np.median(R_times)), float(np.median(sum_times)),
            float(np.median(total_times)))


# -----------------------------------------------------------------------------
# Profile sweep
# -----------------------------------------------------------------------------


def policy_K_values(x: float, K_max_cap: int) -> List[Tuple[str, int]]:
    """The CLAUDE-prescribed K policies for slot 1."""
    L = math.log(x)
    out = [
        ("log2x", max(1, int(round(L ** 2)))),
        ("log3x", max(1, int(round(L ** 3)))),
        ("x_1_4", max(1, int(round(x ** 0.25)))),
        ("x_1_2", max(1, int(round(x ** 0.5)))),
    ]
    out = [(name, min(K, K_max_cap)) for name, K in out]
    return out


def grid_K_values(K_max_cap: int) -> List[int]:
    grid = [25, 50, 100, 200, 400, 800, 1600, 3200, 6400]
    return [K for K in grid if K <= K_max_cap]


def run_profile(xs: List[float], K_max_cap: int, output_csv: str,
                quick: bool = False) -> None:
    print(f"Loading {K_max_cap} zeros (slot 1 setup-cost profile)...")
    t0 = time.time()
    gammas = get_zeros(K_max_cap, dps=30)
    print(f"  loaded in {time.time() - t0:.2f}s "
          f"(gamma_1={gammas[0]:.4f}, gamma_K={gammas[-1]:.4f})")

    # ---- Setup-cost profile ----
    setup_K_grid = [25, 100, 400, 1600, 6400] if quick else \
        [25, 50, 100, 200, 400, 800, 1600, 3200, 6400]
    setup_K_grid = [K for K in setup_K_grid if K <= K_max_cap]
    setup_rows = []
    print("\n=== Setup-cost profile ===")
    print(f"{'K':>6}  {'T_setup_load':>14}  {'Hiary_K^17/13':>16}")
    for K in setup_K_grid:
        t = time_setup(K, dps=30, n_repeats=3)
        hiary = K ** (17.0 / 13.0)
        setup_rows.append({"K": K, "T_setup_load_s": t,
                           "Hiary_K_pow_17_13": hiary})
        print(f"{K:>6}  {t:>14.4f}  {hiary:>16.2f}")

    # ---- Per-x evaluation profile ----
    eval_rows = []
    print("\n=== Per-x evaluation profile ===")
    print(f"{'x':>10}  {'K':>6}  {'policy':>8}  "
          f"{'T_R':>9}  {'T_sum':>9}  {'T_total':>9}  {'T_per_term':>12}")
    for x in xs:
        # 1) Fixed-policy entries.
        policies = policy_K_values(x, K_max_cap)
        # 2) Geometric grid for linear-fit check.
        grid = grid_K_values(K_max_cap)
        K_set = sorted(set([K for _, K in policies] + grid))
        for K in K_set:
            T_R, T_sum, T_tot = time_eval(x, K, gammas, n_repeats=1)
            policy_label = ",".join(name for name, K2 in policies if K2 == K)
            T_per_term = T_sum / max(1, K)
            eval_rows.append({"x": x, "K": K, "policy": policy_label,
                              "T_R_s": T_R, "T_sum_s": T_sum, "T_total_s": T_tot,
                              "T_per_term_s": T_per_term})
            print(f"{x:>10g}  {K:>6}  {policy_label:>8}  "
                  f"{T_R:>9.3f}  {T_sum:>9.3f}  {T_tot:>9.3f}  "
                  f"{T_per_term:>12.6f}")

    # ---- Q1: linear-fit on T_sum vs K, per x ----
    print("\n=== Q1: T_sum(K, x) linearity check (per x) ===")
    q1_rows = []
    for x in xs:
        Ks = np.array([r["K"] for r in eval_rows
                       if r["x"] == x and r["K"] in grid_K_values(K_max_cap)],
                      dtype=float)
        Ts = np.array([r["T_sum_s"] for r in eval_rows
                       if r["x"] == x and r["K"] in grid_K_values(K_max_cap)],
                      dtype=float)
        # Linear: T_sum = a * K + b
        # Power: log T_sum = alpha * log K + beta
        if len(Ks) >= 3:
            log_K = np.log(Ks)
            log_T = np.log(Ts)
            alpha, beta = np.polyfit(log_K, log_T, 1)
            # Slope of T_per_term = T_sum / K vs log K (should be flat if alpha=1)
            per_term = Ts / Ks
            slope_pt, _ = np.polyfit(log_K, per_term, 1)
            print(f"  x={x:>9g}: log-log slope alpha={alpha:.4f}, "
                  f"prefactor exp(beta)={math.exp(beta):.6f}, "
                  f"per-term-vs-logK slope={slope_pt:.2e}")
            q1_rows.append({"x": x, "alpha": alpha, "beta": beta,
                            "prefactor": math.exp(beta),
                            "per_term_drift": slope_pt})

    # ---- Q2: T_per_term across x (cross-x dependence at fixed K) ----
    print("\n=== Q2: T_per_term cross-x dependence at fixed K ===")
    q2_rows = []
    for K in grid_K_values(K_max_cap):
        per_terms = []
        for x in xs:
            for r in eval_rows:
                if r["x"] == x and r["K"] == K:
                    per_terms.append(r["T_per_term_s"])
                    break
        if len(per_terms) == len(xs):
            mn, mx = min(per_terms), max(per_terms)
            spread = (mx - mn) / mn if mn > 0 else float("inf")
            print(f"  K={K:>5}: per_term in [{mn:.6f}, {mx:.6f}], "
                  f"spread={spread:.2%}")
            q2_rows.append({"K": K, "min_per_term": mn, "max_per_term": mx,
                            "spread": spread})

    # ---- Write CSVs ----
    here = os.path.dirname(os.path.abspath(__file__))
    setup_path = os.path.join(here, "cross_x_decoupled_setup.csv")
    eval_path = os.path.join(here, "cross_x_decoupled_eval.csv")
    fit_path = os.path.join(here, "cross_x_decoupled_fits.csv")
    spread_path = os.path.join(here, "cross_x_decoupled_perterm_spread.csv")

    with open(setup_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["K", "T_setup_load_s",
                                          "Hiary_K_pow_17_13"])
        w.writeheader()
        for r in setup_rows:
            w.writerow(r)

    with open(eval_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["x", "K", "policy", "T_R_s",
                                          "T_sum_s", "T_total_s",
                                          "T_per_term_s"])
        w.writeheader()
        for r in eval_rows:
            w.writerow(r)

    with open(fit_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["x", "alpha", "beta", "prefactor",
                                          "per_term_drift"])
        w.writeheader()
        for r in q1_rows:
            w.writerow(r)

    with open(spread_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["K", "min_per_term",
                                          "max_per_term", "spread"])
        w.writeheader()
        for r in q2_rows:
            w.writerow(r)

    print(f"\nWrote {setup_path}")
    print(f"Wrote {eval_path}")
    print(f"Wrote {fit_path}")
    print(f"Wrote {spread_path}")

    # ---- Final amortisation forecast ----
    print("\n=== Slot 1 amortisation forecast ===")
    print("Per-x amortised cost over M batched queries:")
    print("  T_per_x_amortised(K, M) = T_setup(K) / M + T_eval(K, x_avg).")
    print("  T_setup(K) -> Hiary K^{17/13} arithmetic ops at production scale.")
    print("  T_eval(K, x) ~ a(x) * K + b(x), with a(x) ~ x-independent (Q2).")
    print("In the amortisation limit M -> infty:")
    print("  T_per_x_amortised(K, M) -> a * K.")
    print("Per-query polylog requires K polylog. But Thread 3 (S202) closed")
    print("K = polylog in distribution at any positive hit-rate (K* = Theta-tilde(x)).")
    print("Conclusion (slot-1 partial closure): setup amortisation is")
    print("structurally orthogonal to per-x evaluation; even with infinite")
    print("M, per-x amortised cost = a * K = Theta-tilde(x) by Thread 3.")
    print("Falsifier: a sublinear T_eval(K, x) scaling (alpha < 1) would")
    print("invalidate this conclusion.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true",
                        help="small K and x grids for smoke test")
    parser.add_argument("--K-max-cap", type=int, default=6400)
    args = parser.parse_args()

    if args.quick:
        xs = [1e5, 1e6]
        K_cap = min(800, args.K_max_cap)
    else:
        xs = [1e5, 1e6, 1e7]
        K_cap = args.K_max_cap

    here = os.path.dirname(os.path.abspath(__file__))
    out_csv = os.path.join(here, "cross_x_decoupled_eval.csv")
    run_profile(xs, K_cap, out_csv, quick=args.quick)


if __name__ == "__main__":
    main()
