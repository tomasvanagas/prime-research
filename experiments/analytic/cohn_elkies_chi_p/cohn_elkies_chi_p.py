"""
D29 — Cohn-Elkies / Delsarte LP bound applied to the prime indicator chi_P on Z.

Refined framing (per agent research): the experiment is a Delsarte-style LP on Z
with the prime pair-correlation as the distance distribution, not a literal
sphere-packing CEL. Cohn-Elkies sharing the Bochner-positivity + minimum-
separation skeleton; Delsarte-LP (Delsarte 1973 association schemes) is the
right citation for autocorrelation-aware density LPs.

LP formulation
--------------
Even f: Z -> R, supp f subset [-T, T]. Variables: f(0), f(1), ..., f(T).
Bochner: f_hat(xi) = f(0) + 2 sum_{t=1}^T f(t) cos(2 pi t xi) >= 0 for xi in [0, 1).

Vanilla CEL LP (no autocorrelation info, only minimum separation r):
  Constraints
    f_hat(xi) >= 0 at M points,
    f(0) = 1,
    f(t) <= 0 for t in [r, T].
  Objective: minimize f_hat(0) = f(0) + 2 sum_{t=1}^T f(t).
  ... but we want max(f_hat(0)) for tighter bound: |S|/N <= 1/f_hat(0).
  So: maximize f_hat(0) = f(0) + 2 sum f(t)
       subject to f(0) = 1, f(t) <= 0 for t >= r, f_hat >= 0.

Delsarte-aware LP (autocorrelation profile g(t) = R_P(t)/|P| known):
  drop f(t) <= 0 pointwise; replace with single aggregated constraint
    sum_{t=1}^T g(t) f(t) <= 0
  (so sum_{x,y in P} f(x-y) <= R_P(0) f(0) = |P|, recovering the bound.)
  Plancherel: sum_{x,y} f(x-y) = (1/N) sum_xi f_hat(xi) |1_S_hat(xi)|^2
            >= (1/N) f_hat(0) |P|^2.
  So |P|/N <= 1/f_hat(0).

Saturation ratio
----------------
S_N(P) := (|P|/N) * f_hat^*(0).  S_N -> 1 means the LP is tight on primes.
For the trivial CEL (only minimum separation), expect S_N ~ |P|/N * O(1) =
O(1/log N), i.e., loose by a log factor (vacuous bound at constant density).
Delsarte-aware version may tighten this if the prime autocorrelation profile
brings the LP optimum down toward |P|/N.

What this experiment actually measures
--------------------------------------
1. R_P(t) for t in [0, T_max] at N in {10^4, 10^5, 10^6} via FFT.
2. Vanilla LP optimum f_hat^*(0) using min-separation r=2 (twin-prime gap).
3. Delsarte-aware LP optimum using observed g_obs(t) = R_P(t)/|P|.
4. Compare to (a) actual |P|/N, (b) random Bernoulli matched-density control
   sampled from same support.
5. Inspect f^* coefficients for modular-form / Eisenstein-series structure.

Falsification: if both vanilla and Delsarte LPs give f_hat^*(0) within a
constant factor of log N (LP optimum scales as 1/log N matching prime density)
the LP ON PRIMES SATURATES (A-grade signal). If both LPs give f_hat^*(0) = O(1)
(LP optimum stays at constant density, factor log N looser than truth), LPs are
loose for primes (B-grade negative result, structural barrier).
"""

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


def sieve_primes(N):
    """Boolean array of length N+1: prime[k] = True iff k is prime."""
    s = np.ones(N + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(N**0.5) + 1):
        if s[i]:
            s[i * i :: i] = False
    return s


def autocorrelation_via_fft(chi, T_max):
    """Compute R(t) = sum_x chi(x) chi(x + t) for t in [0, T_max] via FFT."""
    n = len(chi)
    # zero-pad to avoid circular wrap when computing finite-Z autocorrelation
    L = 1
    while L < 2 * n:
        L *= 2
    F = np.fft.rfft(chi.astype(np.float64), n=L)
    P = (F * np.conj(F)).real
    full = np.fft.irfft(P, n=L)
    # full[t] = sum_x chi(x) chi((x+t) mod L); for t in [0, n-1] this equals R(t) on Z.
    return full[: T_max + 1]


def build_lp_vanilla(g_obs, R_P_0_over_P, T, M=2048, r_min=2, normalize="f0=1"):
    """
    Vanilla CEL LP: variables = f(0), f(1), ..., f(T).
      maximize f_hat(0) = f(0) + 2 sum_{t=1..T} f(t)
      subject to:
        f(0) = 1                              (normalization)
        f(t) <= 0 for t in [r_min, T]         (sphere packing min-sep r=r_min)
        f_hat(xi_k) >= 0 at M sample points   (Bochner)
        no upper bound on f(t) for t in [1, r_min-1]
    """
    nvar = T + 1
    xi = np.linspace(0.0, 0.5, M)
    # f_hat(xi) = f(0) + 2 sum_{t>=1} f(t) cos(2 pi t xi)
    # Bochner constraint: f_hat >= 0  =>  -f_hat <= 0
    A_bochner = np.zeros((M, nvar))
    A_bochner[:, 0] = 1.0
    for t in range(1, T + 1):
        A_bochner[:, t] = 2.0 * np.cos(2.0 * np.pi * t * xi)
    A_ub = -A_bochner  # constraint -f_hat <= 0
    b_ub = np.zeros(M)

    # f(t) <= 0 for t in [r_min, T]: identity rows.
    for t in range(r_min, T + 1):
        row = np.zeros(nvar)
        row[t] = 1.0
        A_ub = np.vstack([A_ub, row])
        b_ub = np.append(b_ub, 0.0)

    # Equality: f(0) = 1.
    A_eq = np.zeros((1, nvar))
    A_eq[0, 0] = 1.0
    b_eq = np.array([1.0])

    # Objective: maximize f_hat(0) = f(0) + 2 sum f(t).
    # linprog minimizes; minimize -f_hat(0).
    c = np.zeros(nvar)
    c[0] = -1.0
    for t in range(1, T + 1):
        c[t] = -2.0

    bounds = [(None, None)] * nvar
    return {
        "c": c,
        "A_ub": A_ub,
        "b_ub": b_ub,
        "A_eq": A_eq,
        "b_eq": b_eq,
        "bounds": bounds,
    }


def build_lp_delsarte(g_obs, T, M=2048):
    """
    Delsarte-aware LP using observed prime autocorrelation g_obs(t) = R_P(t)/|P|.

    Variables: f(0), f(1), ..., f(T). Even f, supp in [-T, T].
      maximize f_hat(0) = f(0) + 2 sum_{t=1..T} f(t)
      subject to:
        f(0) = 1                                          (normalization)
        sum_{t=1..T} g_obs(t) f(t) <= 0                   (aggregate Delsarte)
        f_hat(xi_k) >= 0 at M sample points               (Bochner)

    Bound: |P|/N <= 1/f_hat^*(0), where f_hat^* is the LP optimum.
    """
    nvar = T + 1
    xi = np.linspace(0.0, 0.5, M)
    A_bochner = np.zeros((M, nvar))
    A_bochner[:, 0] = 1.0
    for t in range(1, T + 1):
        A_bochner[:, t] = 2.0 * np.cos(2.0 * np.pi * t * xi)
    A_ub = -A_bochner
    b_ub = np.zeros(M)

    # Aggregate Delsarte: sum_{t>=1} g_obs(t) f(t) <= 0
    row = np.zeros(nvar)
    for t in range(1, T + 1):
        row[t] = float(g_obs[t])
    A_ub = np.vstack([A_ub, row])
    b_ub = np.append(b_ub, 0.0)

    A_eq = np.zeros((1, nvar))
    A_eq[0, 0] = 1.0
    b_eq = np.array([1.0])

    c = np.zeros(nvar)
    c[0] = -1.0
    for t in range(1, T + 1):
        c[t] = -2.0

    bounds = [(None, None)] * nvar
    return {
        "c": c,
        "A_ub": A_ub,
        "b_ub": b_ub,
        "A_eq": A_eq,
        "b_eq": b_eq,
        "bounds": bounds,
    }


def solve_lp(prob, label=""):
    t0 = time.time()
    res = linprog(
        prob["c"],
        A_ub=prob["A_ub"],
        b_ub=prob["b_ub"],
        A_eq=prob["A_eq"],
        b_eq=prob["b_eq"],
        bounds=prob["bounds"],
        method="highs",
    )
    dt = time.time() - t0
    if not res.success:
        print(f"  [{label}] LP FAILED: {res.message}")
        return None
    f = res.x
    # f_hat(0) = f(0) + 2 sum_{t>=1} f(t) = -res.fun (we minimised -f_hat(0)).
    f_hat_0 = -res.fun
    print(f"  [{label}] LP solved in {dt:.2f}s, f_hat(0) = {f_hat_0:.6e}, "
          f"density bound 1/f_hat(0) = {1.0/f_hat_0:.6e}")
    return {"f": f.tolist(), "f_hat_0": float(f_hat_0), "time_s": dt}


def hardy_littlewood_singular_series(t, t_cap=200):
    """HL singular series S(t) for prime gap t.

    S(t) = 0 if t is odd, else 2 C_2 prod_{p | t/2, p > 2} (p-1)/(p-2).
    Returned per the standard convention of E. Bombieri / HL conjecture.
    """
    if t == 0:
        return float("nan")  # not defined at 0 in usual sense
    if t % 2 == 1:
        return 0.0
    # twin prime constant 2 C_2 = 2 * prod_{p>2} (1 - 1/(p-1)^2) ~ 1.32032...
    C2 = 0.6601618158468695  # twin-prime constant
    val = 2.0 * C2
    m = t // 2
    # multiplicative factor over odd prime divisors of m (NOTE: of t, not t/2)
    # Actually the standard formula is over p | t, p > 2.
    p_list = []
    n = abs(t)
    p = 3
    while p * p <= n:
        if n % p == 0:
            p_list.append(p)
            while n % p == 0:
                n //= p
        p += 2
    if n > 1 and n != 1:
        if n > 2:
            p_list.append(n)
    for p in p_list:
        val *= (p - 1) / (p - 2)
    return val


def hardy_littlewood_g(N, T):
    """HL prediction for g(t) = R_P(t)/|P| at large N.

    For t even: g(t) ~ S(t) / log N. (Approximation; exact to leading order.)
    For t odd:  g(t) ~ 0 (only t=2 odd dominant, but t odd means one of p,q is 2.)
    """
    g = np.zeros(T + 1)
    g[0] = 1.0  # by definition
    log_N = math.log(N)
    for t in range(1, T + 1):
        S_t = hardy_littlewood_singular_series(t)
        # HL: R_P(t) ~ S(t) N / log^2 N, so g(t) = R_P(t) / |P| ~ S(t)/log N
        # (using |P| ~ N/log N).
        g[t] = S_t / log_N
    return g


def run_one_N(N, T_max=1000, M=2048, output_dir=None):
    print(f"\n=== N = {N:,} ===")
    print(f"  sieving primes up to {N}...")
    t0 = time.time()
    chi = sieve_primes(N)
    pi_N = int(chi.sum())
    print(f"  pi(N) = {pi_N:,}, |P|/N = {pi_N/(N+1):.6e} = 1/{(N+1)/pi_N:.3f}")
    print(f"  log N = {math.log(N):.4f}, expected 1/log N = {1.0/math.log(N):.6e}")

    print(f"  computing R_P(t) for t in [0, {T_max}] via FFT...")
    R = autocorrelation_via_fft(chi, T_max)
    print(f"  R_P(0) = {int(R[0])} (= pi(N) - 1 padding artifact, use |P|)")
    # use |P| from sieve directly
    Pcount = pi_N
    g_obs = R / Pcount
    g_obs[0] = 1.0  # exact
    print(f"  g_obs(2) = {g_obs[2]:.6e} (twin-prime ratio)")
    print(f"  g_obs(6) = {g_obs[6]:.6e}, g_obs(30) = {g_obs[30]:.6e}")

    g_HL = hardy_littlewood_g(N, T_max)

    # Vanilla CEL LP
    print(f"  building vanilla CEL LP (T={T_max}, M={M})...")
    prob_v = build_lp_vanilla(g_obs, 1.0, T_max, M=M, r_min=2)
    sol_v = solve_lp(prob_v, label="vanilla CEL r=2")

    # Delsarte LP using observed g_obs
    print(f"  building Delsarte LP (T={T_max}, M={M})...")
    prob_d = build_lp_delsarte(g_obs, T_max, M=M)
    sol_d = solve_lp(prob_d, label="Delsarte g_obs")

    # Delsarte LP using HL prediction g_HL
    print(f"  building Delsarte LP with HL g_HL...")
    prob_h = build_lp_delsarte(g_HL, T_max, M=M)
    sol_h = solve_lp(prob_h, label="Delsarte g_HL")

    # Random Bernoulli control: g_rand(t) = (|P|/N)^2 * N / |P| = |P|/N for t != 0.
    # Actually for Bernoulli with density rho: E[R(t)] = rho^2 (N - t), so
    # g_rand(t) = rho^2 (N-t) / (rho N) = rho * (1 - t/N) ~ rho.
    rho = pi_N / (N + 1)
    g_rand = np.full(T_max + 1, rho)
    g_rand[0] = 1.0
    prob_r = build_lp_delsarte(g_rand, T_max, M=M)
    sol_r = solve_lp(prob_r, label="Delsarte g_rand=rho")

    # Compute saturation ratios.
    rho_actual = pi_N / (N + 1)

    def sat(sol):
        if sol is None:
            return None
        return rho_actual * sol["f_hat_0"]

    sat_v = sat(sol_v)
    sat_d = sat(sol_d)
    sat_h = sat(sol_h)
    sat_r = sat(sol_r)

    print(f"\n  density rho = pi(N)/(N+1) = {rho_actual:.6e}")
    print(f"  saturation S_N(vanilla CEL)   = rho * f_hat^* = {sat_v}")
    print(f"  saturation S_N(Delsarte g_obs)= {sat_d}")
    print(f"  saturation S_N(Delsarte g_HL) = {sat_h}")
    print(f"  saturation S_N(Bernoulli rho) = {sat_r}")
    print(f"  bound 1/f_hat^* (vanilla)     = {1.0/sol_v['f_hat_0'] if sol_v else None}")
    print(f"  bound 1/f_hat^* (g_obs)       = {1.0/sol_d['f_hat_0'] if sol_d else None}")
    print(f"  bound 1/f_hat^* (g_HL)        = {1.0/sol_h['f_hat_0'] if sol_h else None}")
    print(f"  bound 1/f_hat^* (g_rand)      = {1.0/sol_r['f_hat_0'] if sol_r else None}")

    summary = {
        "N": N,
        "pi_N": pi_N,
        "rho": rho_actual,
        "log_N": math.log(N),
        "T_max": T_max,
        "M": M,
        "g_obs_first10": g_obs[:10].tolist(),
        "g_HL_first10": g_HL[:10].tolist(),
        "vanilla_CEL": {
            "f_hat_0": sol_v["f_hat_0"] if sol_v else None,
            "bound": 1.0 / sol_v["f_hat_0"] if sol_v else None,
            "saturation": sat_v,
            "f_first10": sol_v["f"][:10] if sol_v else None,
        },
        "delsarte_g_obs": {
            "f_hat_0": sol_d["f_hat_0"] if sol_d else None,
            "bound": 1.0 / sol_d["f_hat_0"] if sol_d else None,
            "saturation": sat_d,
            "f_first10": sol_d["f"][:10] if sol_d else None,
        },
        "delsarte_g_HL": {
            "f_hat_0": sol_h["f_hat_0"] if sol_h else None,
            "bound": 1.0 / sol_h["f_hat_0"] if sol_h else None,
            "saturation": sat_h,
        },
        "delsarte_g_rand": {
            "f_hat_0": sol_r["f_hat_0"] if sol_r else None,
            "bound": 1.0 / sol_r["f_hat_0"] if sol_r else None,
            "saturation": sat_r,
        },
    }

    if output_dir:
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        out = Path(output_dir) / f"results_N{N}.json"
        with open(out, "w") as f:
            json.dump(summary, f, indent=2, default=str)
        print(f"  -> wrote {out}")

        # also save full f vectors for the largest LP solution
        if sol_d is not None:
            out_f = Path(output_dir) / f"f_vector_N{N}_delsarte_g_obs.json"
            with open(out_f, "w") as f:
                json.dump(
                    {
                        "N": N,
                        "T_max": T_max,
                        "f": sol_d["f"],
                        "g_obs": g_obs.tolist(),
                    },
                    f,
                )
            print(f"  -> wrote {out_f}")

    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, nargs="+", default=[10_000, 100_000, 1_000_000])
    parser.add_argument("--T_max", type=int, default=500)
    parser.add_argument("--M", type=int, default=4096)
    parser.add_argument(
        "--output-dir",
        default="/apps/aplikacijos/prime-research/experiments/analytic/cohn_elkies_chi_p",
    )
    args = parser.parse_args()

    print(f"Cohn-Elkies / Delsarte LP on chi_P")
    print(f"N values: {args.N}")
    print(f"T_max = {args.T_max}, M = {args.M}")

    results = []
    for N in args.N:
        try:
            r = run_one_N(N, T_max=args.T_max, M=args.M, output_dir=args.output_dir)
            results.append(r)
        except Exception as e:
            print(f"  ERROR at N={N}: {e}")
            import traceback

            traceback.print_exc()

    # print summary table
    print(f"\n=== SUMMARY ===")
    print(f"{'N':>10} {'rho':>12} {'1/f_hat*_van':>15} {'1/f_hat*_obs':>15} "
          f"{'1/f_hat*_HL':>15} {'1/f_hat*_rand':>15} {'S_van':>8} {'S_obs':>8} "
          f"{'S_HL':>8} {'S_rand':>8}")
    for r in results:
        N = r["N"]
        rho = r["rho"]
        bv = r["vanilla_CEL"]["bound"]
        bo = r["delsarte_g_obs"]["bound"]
        bh = r["delsarte_g_HL"]["bound"]
        br = r["delsarte_g_rand"]["bound"]
        sv = r["vanilla_CEL"]["saturation"]
        so = r["delsarte_g_obs"]["saturation"]
        sh = r["delsarte_g_HL"]["saturation"]
        sr = r["delsarte_g_rand"]["saturation"]
        print(
            f"{N:>10} {rho:>12.6e} {bv:>15.6e} {bo:>15.6e} {bh:>15.6e} {br:>15.6e} "
            f"{sv:>8.4f} {so:>8.4f} {sh:>8.4f} {sr:>8.4f}"
        )

    if args.output_dir:
        out = Path(args.output_dir) / "summary.json"
        with open(out, "w") as f:
            json.dump(results, f, indent=2, default=str)
        print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
