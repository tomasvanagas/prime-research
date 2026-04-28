"""
T_max sweep for the Delsarte LP on primes. The headline result from the main
script is that the LP optimum 1/f_hat^*(0) is asymptotically constant in N
(while prime density rho ~ 1/log N), so saturation S_N -> 0 like 1/log N.

This script tests whether 1/f_hat^*(0) varies with T_max. If 1/f_hat^*(0)
saturates at fixed T_max (ie further T does not tighten), the bound has
plateaued at a structural constant set by the singular-series shape on
[0, T_plateau]. Test T_max in {50, 100, 200, 400, 800, 1500}.
"""

import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

import sys
sys.path.insert(0, str(Path(__file__).parent))
from cohn_elkies_chi_p import (
    sieve_primes, autocorrelation_via_fft,
    build_lp_delsarte, hardy_littlewood_g,
)


def main():
    N = 1_000_000
    M = 4096
    T_values = [50, 100, 200, 400, 800, 1500]

    print(f"T_max sweep at N = {N}, M = {M}")
    print("Sieving primes...")
    chi = sieve_primes(N)
    pi_N = int(chi.sum())
    rho = pi_N / (N + 1)
    log_N = math.log(N)
    print(f"pi(N) = {pi_N:,}, rho = {rho:.6e}, 1/log N = {1.0/log_N:.6e}")

    # autocorrelation up to max T
    T_max = max(T_values)
    print(f"computing R_P(t) up to T_max = {T_max}...")
    R = autocorrelation_via_fft(chi, T_max)
    g_obs = R / pi_N
    g_obs[0] = 1.0

    g_HL = hardy_littlewood_g(N, T_max)

    results = []
    for T in T_values:
        print(f"\n--- T_max = {T} ---")

        prob_o = build_lp_delsarte(g_obs[: T + 1], T, M=M)
        t0 = time.time()
        res_o = linprog(
            prob_o["c"], A_ub=prob_o["A_ub"], b_ub=prob_o["b_ub"],
            A_eq=prob_o["A_eq"], b_eq=prob_o["b_eq"],
            bounds=prob_o["bounds"], method="highs",
        )
        dt_o = time.time() - t0
        f_hat_obs = -res_o.fun if res_o.success else None
        print(f"  g_obs : f_hat^*(0) = {f_hat_obs:.4f}, "
              f"bound = {1.0/f_hat_obs:.4e}, S_N = {rho*f_hat_obs:.4f} "
              f"(solve {dt_o:.1f}s)")

        prob_h = build_lp_delsarte(g_HL[: T + 1], T, M=M)
        t0 = time.time()
        res_h = linprog(
            prob_h["c"], A_ub=prob_h["A_ub"], b_ub=prob_h["b_ub"],
            A_eq=prob_h["A_eq"], b_eq=prob_h["b_eq"],
            bounds=prob_h["bounds"], method="highs",
        )
        dt_h = time.time() - t0
        f_hat_HL = -res_h.fun if res_h.success else None
        print(f"  g_HL  : f_hat^*(0) = {f_hat_HL:.4f}, "
              f"bound = {1.0/f_hat_HL:.4e}, S_N = {rho*f_hat_HL:.4f} "
              f"(solve {dt_h:.1f}s)")

        results.append({
            "T": T,
            "f_hat_obs": float(f_hat_obs) if f_hat_obs else None,
            "bound_obs": float(1.0 / f_hat_obs) if f_hat_obs else None,
            "S_obs": float(rho * f_hat_obs) if f_hat_obs else None,
            "f_hat_HL": float(f_hat_HL) if f_hat_HL else None,
            "bound_HL": float(1.0 / f_hat_HL) if f_hat_HL else None,
            "S_HL": float(rho * f_hat_HL) if f_hat_HL else None,
            "f_obs_first10": res_o.x[:10].tolist() if res_o.success else None,
        })

    print("\n\n=== T-sweep summary ===")
    print(f"{'T':>6} {'f_obs':>10} {'1/f_obs':>10} {'f_HL':>10} {'1/f_HL':>10}")
    for r in results:
        print(f"{r['T']:>6} {r['f_hat_obs']:>10.4f} {r['bound_obs']:>10.4e} "
              f"{r['f_hat_HL']:>10.4f} {r['bound_HL']:>10.4e}")

    out = Path(__file__).parent / "T_sweep_results.json"
    with open(out, "w") as f:
        json.dump({"N": N, "rho": rho, "log_N": log_N, "results": results}, f, indent=2)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
