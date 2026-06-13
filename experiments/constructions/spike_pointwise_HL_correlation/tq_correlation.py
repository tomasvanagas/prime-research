"""
tq_correlation.py
=================

C9.b: composition of the S191 pointwise spike approximator T_Q(n) with
the Hardy-Littlewood twin-prime singular series.

For positive integer N = 2^d, integer Q >= 1 and shift h >= 0, this
script computes:

  R_h(Q, N)        = (1/(N-h)) sum_{n=1}^{N-h} T_Q(n) T_Q(n+h),
  R_h^conn         = R_h - mean(T_Q)^2,
  S_Q(h)           = sum_{q sqf <= Q} mu(q)^2 / phi(q)^2 * c_q(h),
  pred_h(Q)        = (pi(N)/N)^2 * (S_Q(h) - 1),    # connected piece
  pi_h(N)          = (1/(N-h)) #{n <= N-h : n, n+h prime},
  pi_h^conn        = pi_h - (pi(N)/N)^2,
  HL_h_conn(N)     = (pi(N)/N)^2 * (S_HL(h) - 1)    # connected HL prediction.

The disconnected q=1 mode contributes the (pi/N)^2 squared-mean, which
is subtracted in both R_h^conn and the prediction by passing to S_Q - 1
(equiv. to summing only q >= 2 in the Ramanujan-Fourier expansion).
The diagonal q1 = q2 contribution to R_h^conn equals (pi/N)^2 * (S_Q(h) - 1)
exactly when N is divisible by every squarefree q <= Q (which holds for
N = 2^d only at q = 1, 2; for other q the finite-N drift contributes a
small correction).

Tests:

  F1 (identity): R_h^conn / pred_h(Q) in [0.85, 1.15] for h in {2,4,6,8,10,12,30,210}, Q in {30, 210, 2310, ~N^0.185}.
  F2 (HL recovery at large Q): R_h^conn / [(pi(N)/N)^2 * (S_HL(h) - 1)] -> 1 at Q ~ sqrt(N), even h.
  F3 (odd-h asymptote): for odd h, S_HL(h) - 1 = -1 and R_h^conn -> -(pi/N)^2 at large Q.
     Ratio R_h^conn / (-(pi/N)^2) -> 1.
  F4 (h=0 self-consistency): var(T_Q) = (pi(N)/N)^2 * (S_Q(0) - 1) where
     S_Q(0) - 1 = sum_{q sqf <= Q, q>=2} 1/phi(q).
  F5 (prime baseline): pi_h^conn / [(pi(N)/N)^2 * (S_HL(h) - 1)] close to 1 at d=20 for tested h.
"""

import json
import math
import time

import numpy as np


def sieve(n_max: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (chi_P, mu, phi) tables for indices 0..n_max."""
    is_prime = np.ones(n_max + 1, dtype=bool)
    is_prime[0] = is_prime[1] = False
    mu = np.ones(n_max + 1, dtype=np.int8)
    mu[0] = 0
    phi = np.arange(n_max + 1, dtype=np.int64)
    for p in range(2, n_max + 1):
        if is_prime[p]:
            for m in range(2 * p, n_max + 1, p):
                is_prime[m] = False
            # mu and phi via the "p-th prime updates"
            for m in range(p, n_max + 1, p):
                mu[m] = -mu[m]
            p2 = p * p
            for m in range(p2, n_max + 1, p2):
                mu[m] = 0
            for m in range(p, n_max + 1, p):
                phi[m] -= phi[m] // p
    chi_P = is_prime.astype(np.int8)
    return chi_P, mu, phi


def ramanujan_pattern(q: int, mu_arr: np.ndarray, phi_arr: np.ndarray) -> np.ndarray:
    """Return length-q array c_q(r) for r = 0..q-1.

    c_q(h) = mu(q/d) * phi(d) where d = gcd(q, h). c_q(0) = phi(q).
    """
    out = np.empty(q, dtype=np.int64)
    out[0] = int(phi_arr[q])
    for r in range(1, q):
        d = math.gcd(q, r)
        out[r] = int(mu_arr[q // d]) * int(phi_arr[d])
    return out


def build_T_Q(N: int, Q: int, mu: np.ndarray, phi: np.ndarray, prime_density: float) -> np.ndarray:
    """T_Q(n) = (pi(N)/N) sum_{q sqf <= Q} (mu(q)/phi(q)) c_q(n) for n in [0, N)."""
    T = np.zeros(N, dtype=np.float64)
    n_indices = np.arange(N, dtype=np.int64)
    for q in range(1, Q + 1):
        if mu[q] == 0:
            continue
        pattern = ramanujan_pattern(q, mu, phi).astype(np.float64)
        weight = float(mu[q]) / float(phi[q])
        T += weight * pattern[n_indices % q]
    T *= prime_density
    return T


def trunc_singular_series(Q: int, h: int, mu: np.ndarray, phi: np.ndarray) -> float:
    """S_Q(h) = sum_{q sqf <= Q} mu(q)^2 / phi(q)^2 * c_q(h)."""
    s = 0.0
    for q in range(1, Q + 1):
        if mu[q] == 0:
            continue
        h_mod = h % q if q > 1 else 0
        d = math.gcd(q, h_mod) if h_mod != 0 else q
        c_q_h = int(mu[q // d]) * int(phi[d])
        s += (1.0 / (float(phi[q]) ** 2)) * c_q_h
    return s


def textbook_HL(h: int, primes_for_C2: list[int]) -> float:
    """Full Hardy-Littlewood singular series S(h)."""
    if h == 0:
        return float('inf')
    if h % 2 == 1:
        return 0.0
    C2 = 1.0
    for p in primes_for_C2:
        if p == 2:
            continue
        C2 *= 1.0 - 1.0 / (p - 1) ** 2
    n = abs(h)
    while n % 2 == 0:
        n //= 2
    factors = []
    p = 3
    while p * p <= n:
        if n % p == 0:
            factors.append(p)
            while n % p == 0:
                n //= p
        p += 2
    if n > 1:
        factors.append(n)
    factor_prod = 1.0
    for p in factors:
        factor_prod *= (p - 1.0) / (p - 2.0)
    return 2.0 * C2 * factor_prod


def correlation_at_shift(T: np.ndarray, h: int) -> float:
    """(1/(N-h)) sum_{n=0}^{N-h-1} T[n] T[n+h]."""
    if h == 0:
        return float(np.mean(T * T))
    return float(np.mean(T[:-h] * T[h:]))


def prime_correlation(chi_P: np.ndarray, h: int) -> float:
    """(1/(N-h)) sum chi_P[n] chi_P[n+h]."""
    if h == 0:
        return float(np.mean(chi_P.astype(np.float64) ** 2))
    return float(np.mean(chi_P[:-h].astype(np.float64) * chi_P[h:].astype(np.float64)))


def main():
    D_LIST = [16, 18, 20]
    H_LIST = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 30, 210]
    # Q values: small (2, 6, 30 = primorial(3), 210 = primorial(4),
    # 2310 = primorial(5)), N^{0.185}, N^{0.21}, sqrt(N).

    out: dict = {}

    # Build a list of small primes for the C_2 calculation in textbook_HL.
    PMAX_C2 = 200_000
    is_prime_small = np.ones(PMAX_C2 + 1, dtype=bool)
    is_prime_small[:2] = False
    for p in range(2, int(PMAX_C2 ** 0.5) + 1):
        if is_prime_small[p]:
            is_prime_small[p * p :: p] = False
    primes_for_C2 = np.where(is_prime_small)[0].tolist()

    # Compute textbook HL once (independent of N).
    HL_textbook = {h: textbook_HL(h, primes_for_C2) for h in H_LIST}

    for d in D_LIST:
        N = 1 << d
        print(f"\n=== d = {d}, N = 2^{d} = {N:,} ===", flush=True)
        t0 = time.time()
        # Sieve up to N for chi_P, and up to enough for mu/phi tables.
        # For mu/phi we need up to max Q in this experiment; ensure
        # large enough.
        Q_VALS = [2, 6, 30, 210, 2310, max(2, round(N ** 0.185)), max(2, round(N ** 0.21)), max(2, round(N ** 0.5))]
        Q_VALS = sorted(set(Q_VALS))
        max_Q_here = max(Q_VALS)
        upper = max(N, max_Q_here)
        print(f"  sieving up to {upper:,} ...", flush=True)
        chi_P, mu, phi = sieve(upper)
        chi_P = chi_P[: N + 1]  # only need up to N
        pi_N = int(chi_P.sum())
        density = pi_N / N
        print(f"  pi(N) = {pi_N}, pi(N)/N = {density:.7f}", flush=True)

        out_d: dict = {
            "N": N,
            "pi_N": pi_N,
            "density": density,
            "Q_vals": Q_VALS,
            "H_vals": H_LIST,
            "results": {},  # results[Q][h] = {...}
            "prime_correlation": {},  # pi_h
            "HL_textbook": HL_textbook,
        }

        # Prime correlation (no Q dependence).
        print("  computing prime correlations ...", flush=True)
        for h in H_LIST:
            pi_h = prime_correlation(chi_P[:N], h)
            pi_h_conn = pi_h - density ** 2
            S_HL = HL_textbook[h]
            HL_full = density ** 2 * S_HL if S_HL != float('inf') else float('inf')
            HL_conn = density ** 2 * (S_HL - 1) if S_HL != float('inf') else float('inf')
            ratio_full = (
                pi_h / HL_full
                if HL_full not in (0, 0.0, float('inf'))
                else None
            )
            ratio_conn = (
                pi_h_conn / HL_conn
                if HL_conn not in (0, 0.0, float('inf'))
                else None
            )
            out_d["prime_correlation"][h] = {
                "pi_h": pi_h,
                "pi_h_conn": pi_h_conn,
                "HL_full": HL_full,
                "HL_conn": HL_conn,
                "ratio_full": ratio_full,
                "ratio_conn": ratio_conn,
            }

        # T_Q correlations.
        for Q in Q_VALS:
            print(f"  Q = {Q}: building T_Q ...", flush=True)
            tq0 = time.time()
            T = build_T_Q(N, Q, mu, phi, density)
            mean_T = float(np.mean(T))
            mean_T_sq = mean_T ** 2
            tq1 = time.time()
            results_Q = {}
            for h in H_LIST:
                R_h = correlation_at_shift(T, h)
                R_h_conn = R_h - mean_T_sq
                S_Q_h = trunc_singular_series(Q, h, mu, phi)
                # Connected prediction: only q >= 2 contributes, equiv. (S_Q(h) - 1).
                pred_conn = density ** 2 * (S_Q_h - 1.0)
                ratio_pred = R_h_conn / pred_conn if pred_conn != 0 else None
                S_HL = HL_textbook[h]
                pred_HL_conn = (
                    density ** 2 * (S_HL - 1.0)
                    if S_HL != float('inf')
                    else float('inf')
                )
                ratio_HL = (
                    R_h_conn / pred_HL_conn
                    if pred_HL_conn not in (0, 0.0, float('inf'))
                    else None
                )
                results_Q[h] = {
                    "R_h": R_h,
                    "R_h_conn": R_h_conn,
                    "S_Q_h": S_Q_h,
                    "pred_conn": pred_conn,
                    "pred_HL_conn": pred_HL_conn,
                    "ratio_pred": ratio_pred,
                    "ratio_HL": ratio_HL,
                }
            results_Q["_meta"] = {
                "mean_T": mean_T,
                "mean_T_sq": mean_T_sq,
                "build_time_s": tq1 - tq0,
            }
            out_d["results"][Q] = results_Q

        out[d] = out_d
        print(f"  d = {d} done in {time.time() - t0:.1f} s.", flush=True)

    # Save raw results.
    out_serial = {}
    for d, d_data in out.items():
        out_serial[str(d)] = {
            **{k: v for k, v in d_data.items() if k != "results"},
            "results": {str(Q): r for Q, r in d_data["results"].items()},
        }
        # Convert int keys in HL_textbook
        out_serial[str(d)]["HL_textbook"] = {str(h): v for h, v in d_data["HL_textbook"].items()}
        out_serial[str(d)]["prime_correlation"] = {
            str(h): v for h, v in d_data["prime_correlation"].items()
        }
        for Q_key, Q_data in out_serial[str(d)]["results"].items():
            new_Q_data = {"_meta": Q_data["_meta"]}
            for h in H_LIST:
                new_Q_data[str(h)] = Q_data[h]
            out_serial[str(d)]["results"][Q_key] = new_Q_data

    with open("tq_correlation_results.json", "w") as f:
        json.dump(out_serial, f, indent=2, default=lambda x: float(x) if isinstance(x, np.floating) else x)
    print("\nWrote tq_correlation_results.json")

    # Print summary tables to stdout.
    print("\n\n=== SUMMARY ===")
    for d in D_LIST:
        d_data = out[d]
        N = d_data["N"]
        density = d_data["density"]
        print(f"\nd = {d}, N = 2^{d}, pi(N)/N = {density:.6f}")
        print("  --- Prime correlation pi_h vs HL: full and connected ---")
        print(
            f"  {'h':>5}  {'pi_h':>13}  {'pi_h^conn':>13}  "
            f"{'HL_full':>13}  {'HL_conn':>13}  {'r_full':>8}  {'r_conn':>8}"
        )
        for h in H_LIST:
            pc = d_data["prime_correlation"][h]
            r_full = f"{pc['ratio_full']:>8.4f}" if pc['ratio_full'] is not None else "       -"
            r_conn = f"{pc['ratio_conn']:>8.4f}" if pc['ratio_conn'] is not None else "       -"
            HL_full_s = (
                f"{pc['HL_full']:>13.6e}" if pc['HL_full'] != float('inf') else "          inf"
            )
            HL_conn_s = (
                f"{pc['HL_conn']:>13.6e}" if pc['HL_conn'] != float('inf') else "          inf"
            )
            print(
                f"  {h:>5}  {pc['pi_h']:>13.6e}  {pc['pi_h_conn']:>+13.6e}  "
                f"{HL_full_s}  {HL_conn_s}  {r_full}  {r_conn}"
            )

        print("\n  --- T_Q autocorrelation: R_h^conn / [(pi/N)^2 (S_Q(h)-1)] (should be ~1) ---")
        Q_VALS = d_data["Q_vals"]
        header = "  h     " + "".join(f"{f'Q={Q}':>13}" for Q in Q_VALS)
        print(header)
        for h in H_LIST:
            row = f"  h={h:>4}"
            for Q in Q_VALS:
                r = d_data["results"][Q][h]
                if r["ratio_pred"] is None:
                    row += f"{'-':>13}"
                else:
                    row += f"{r['ratio_pred']:>13.5f}"
            print(row)

        print("\n  --- T_Q autocorrelation: R_h^conn / [(pi/N)^2 (S_HL(h)-1)] (should -> 1 at large Q) ---")
        print(header)
        for h in H_LIST:
            row = f"  h={h:>4}"
            for Q in Q_VALS:
                r = d_data["results"][Q][h]
                if r["ratio_HL"] is None:
                    row += f"{'-':>13}"
                else:
                    row += f"{r['ratio_HL']:>13.5f}"
            print(row)


if __name__ == "__main__":
    main()
