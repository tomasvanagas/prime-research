"""
tq_triple_correlation.py
========================

Three-point companion of the S205 (T_Q autocorrelation = truncated HL
twin-prime singular series) construction. We compute

    R_{h1, h2}(Q, N)  :=  (1/(N - max(h1, h2))) * sum_n T_Q(n) T_Q(n+h1) T_Q(n+h2)

for various (Q, h1, h2) and compare with the prime-by-prime closed form

    pred(Q, h1, h2)  =  (pi(N)/N)^3 * prod_{p sqf primes <= Q} G_p(h1, h2)

where  G_p(h1, h2) = (p - nu_p({0, h1, h2})) * p^2 / (p-1)^3,
       nu_p = number of distinct residues among {0, h1, h2} mod p.

For *primorial* W, T_W^{div}(n) = (pi(N)/N) * (W/phi(W)) * [gcd(n, W) = 1]
collapses pointwise (S208), and G_p truncated to p|W gives the
W-conductor truncation of the Hardy-Littlewood prime triple singular
series S_HL^{(W)}(0, h1, h2).

The script verifies five falsifications:

  F1  R_{h1,h2}(W, N) / [(pi/N)^3 prod_{p|W} G_p] -> 1 at primorial W.
  F2  R_{h1,h2}(Q, N) / [(pi/N)^3 prod_{p sqf primes <= Q} G_p] -> 1 at general Q.
  F3  HL recovery at large Q (Q ~ sqrt(N)): ratio to full
      S_HL(0, h1, h2) (over ALL primes).
  F4  h1 = h2 = 0 self-consistency: <T_W^3> = (pi/N)^3 (W/phi(W))^2.
  F5  Reduction to S205 2-point at h2 = 0: triple identity matches the
      pair identity multiplied by the constant (W/phi(W)).
  F6  prime triple density check at finite N.

Run time at d in {16, 18, 20} on a laptop: roughly 2-4 minutes total.
"""

import json
import math
import time
from itertools import product

import numpy as np


def sieve(n_max: int):
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
    """Return length-q array c_q(r) for r = 0..q-1.  c_q(h) = mu(q/d) phi(d), d = gcd(q, h)."""
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


def build_T_W_div(N: int, W: int, prime_density: float) -> np.ndarray:
    """T_W^{div}(n) = (pi(N)/N) * (W/phi(W)) * [gcd(n, W) = 1]. Squarefree W."""
    # Compute phi(W) directly.
    phi_W = W
    p = 2
    n = W
    seen = set()
    while p * p <= n:
        if n % p == 0:
            seen.add(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        seen.add(n)
    for q in seen:
        phi_W = phi_W * (q - 1) // q
    n_arr = np.arange(N, dtype=np.int64)
    coprime = np.ones(N, dtype=bool)
    for p in seen:
        coprime &= (n_arr % p != 0)
    out = np.where(coprime, 1.0, 0.0)
    out *= prime_density * (W / phi_W)
    return out


def nu_p(p: int, h_set: tuple) -> int:
    """Number of distinct residues mod p of the elements of h_set."""
    return len({h % p for h in h_set})


def G_p(p: int, h1: int, h2: int) -> float:
    """Prime-p factor of the HL prime triple singular series.

    G_p = (p - nu_p) * p^2 / (p - 1)^3
    """
    nu = nu_p(p, (0, h1, h2))
    return (p - nu) * (p * p) / ((p - 1) ** 3)


def f_p_direct(p: int, h1: int, h2: int) -> float:
    """Compute (1/p) sum_{r mod p} c_p(r) c_p(r+h1) c_p(r+h2) directly.

    For squarefree prime p, c_p(r) = -1 if p does not divide r, p-1 if p | r.
    """
    total = 0
    for r in range(p):
        cp_r = (p - 1) if r % p == 0 else -1
        cp_r1 = (p - 1) if (r + h1) % p == 0 else -1
        cp_r2 = (p - 1) if (r + h2) % p == 0 else -1
        total += cp_r * cp_r1 * cp_r2
    return total / p


def G_p_via_f(p: int, h1: int, h2: int) -> float:
    """Compute G_p via the Ramanujan-Fourier formula

    G_p = 1 + (1/(p-1)^2) [c_p(h1) + c_p(h2) + c_p(h2-h1)] - f_p / (p-1)^3
    """
    def cp(h):
        return (p - 1) if h % p == 0 else -1

    s = 1.0
    s += (cp(h1) + cp(h2) + cp(h2 - h1)) / ((p - 1) ** 2)
    s -= f_p_direct(p, h1, h2) / ((p - 1) ** 3)
    return s


def primes_up_to(N_max: int) -> list:
    sieve_arr = np.ones(N_max + 1, dtype=bool)
    sieve_arr[:2] = False
    for p in range(2, int(N_max ** 0.5) + 1):
        if sieve_arr[p]:
            sieve_arr[p * p :: p] = False
    return list(np.where(sieve_arr)[0])


def trunc_HL_triple_via_primes(prime_list: list, h1: int, h2: int) -> float:
    """Partial HL triple singular series over given primes p.

       prod_{p in prime_list} G_p(h1, h2) = prod (p - nu_p) p^2 / (p-1)^3.
    """
    s = 1.0
    for p in prime_list:
        s *= G_p(p, h1, h2)
    return s


def full_HL_triple(h1: int, h2: int, primes_for_C3: list) -> float:
    """Hardy-Littlewood prime triple singular series S(0, h1, h2)
    = prod_p (1 - nu_p/p) (1 - 1/p)^{-3}.

    The product converges; we truncate at primes_for_C3[-1] (~ 200 000)
    for ~ 1e-5 relative accuracy.
    """
    if 0 == h1 == h2:
        return float('inf')
    # Check admissibility: the tuple {0, h1, h2} is admissible iff for
    # every prime p, nu_p < p.
    for p in primes_for_C3:
        if nu_p(p, (0, h1, h2)) == p:
            # Inadmissible: triple covers all residues mod p, so HL = 0
            return 0.0
        if p > max(2, abs(h1), abs(h2), abs(h2 - h1), abs(h1 - h2)) + 5:
            # for p larger than max gap, nu_p = 3 (all distinct), and the
            # tail is bounded; we proceed.
            break
    s = 1.0
    for p in primes_for_C3:
        s *= G_p(p, h1, h2)
    return s


def correlation_triple(T: np.ndarray, h1: int, h2: int) -> float:
    """(1/L) sum T[n] T[n+h1] T[n+h2] over n with valid indices."""
    h_max = max(0, h1, h2)
    if h_max == 0:
        return float(np.mean(T * T * T))
    L = T.shape[0] - h_max
    a = T[:L]
    b = T[h1 : h1 + L]
    c = T[h2 : h2 + L]
    return float(np.mean(a * b * c))


def prime_triple_density(chi_P: np.ndarray, h1: int, h2: int) -> float:
    """Empirical (1/L) #{n: chi_P(n) = chi_P(n+h1) = chi_P(n+h2) = 1}."""
    h_max = max(0, h1, h2)
    if h_max == 0:
        return float(np.mean(chi_P.astype(np.float64) ** 3))
    L = chi_P.shape[0] - h_max
    a = chi_P[:L].astype(np.float64)
    b = chi_P[h1 : h1 + L].astype(np.float64)
    c = chi_P[h2 : h2 + L].astype(np.float64)
    return float(np.mean(a * b * c))


def squarefree_primes_up_to(Q: int, mu_arr: np.ndarray) -> list:
    """Primes p <= Q with mu(p) = -1."""
    return [p for p in range(2, Q + 1) if mu_arr[p] == -1 and p > 1]


def main():
    D_LIST = [16, 18, 20]
    # Triples (h1, h2) for the prime triple density. We focus on
    # *admissible* triples (those that don't cover all residues mod 2:
    # 0, h1, h2 should not all have distinct parities). Some inadmissible
    # ones included as F6 sanity.
    H_PAIRS = [
        (0, 0),     # F4: h1=h2=0 self-consistency
        (1, 0),     # h2=0 reduction to pair
        (2, 0),     # h2=0 reduction
        (6, 0),     # h2=0 reduction
        (2, 4),     # 0, 2, 4: even, ν_2=1; ν_3=2 (0,2,4 mod 3 = 0,2,1, distinct)
        (4, 2),     # symmetry
        (6, 12),    # 0, 6, 12: ν_2=1, ν_3=1, ν_5=2 (admissible iff at every p: ν_p<p)
        (2, 6),     # admissible: ν_2=1, ν_3=2 (0,2,0 mod 3 wait: 0,2,6 mod 3 = 0,2,0, ν_3=2)
        (2, 8),     # ν_2=1, ν_3=3 (0,2,8 mod 3 = 0,2,2, ν_3=2 actually)
        (4, 6),     # ν_2=1, ν_3=2 (0,4,6 mod 3 = 0,1,0, ν_3=2)
        (6, 30),    # 0,6,30: ν_2=1, ν_3=1, ν_5=1, ν_7=3, ...
        (30, 60),   # 0,30,60: many primes get ν_p=1
        (1, 3),     # F6: 0,1,3 mod 2 = {0,1,1} -> ν_2=2, mod 3 = {0,1,0} -> ν_3=2 admissible
        (1, 5),     # 0,1,5: mod 2 = {0,1,1}, mod 3 = {0,1,2} ν_3=3 INADMISSIBLE
        (1, 2),     # 0,1,2: mod 2={0,1,0}, ν_2=2; mod 3={0,1,2}, ν_3=3 INADMISSIBLE
    ]

    primes_for_C3 = primes_up_to(50000)

    out = {}

    for d in D_LIST:
        N = 1 << d
        print(f"\n=== d = {d}, N = 2^{d} = {N:,} ===", flush=True)
        t0 = time.time()
        # Q values to test for the GENERAL T_Q construction:
        Q_VALS = [2, 6, 30, 210, 2310, max(2, round(N ** 0.185)), max(2, round(N ** 0.5))]
        Q_VALS = sorted(set(Q_VALS))
        # Primorial conductors W to test for the DIVISOR-RESTRICTED T_W^div:
        W_VALS = [2, 6, 30, 210, 2310]
        upper = max(N, max(Q_VALS) + 1)
        print(f"  sieving up to {upper:,} ...", flush=True)
        chi_P, mu, phi = sieve(upper)
        chi_P = chi_P[: N + 1]
        pi_N = int(chi_P.sum())
        density = pi_N / N
        print(f"  pi(N) = {pi_N}, pi(N)/N = {density:.7f}", flush=True)

        out_d = {
            "N": N,
            "pi_N": pi_N,
            "density": density,
            "Q_vals": Q_VALS,
            "W_vals": W_VALS,
            "H_pairs": [list(p) for p in H_PAIRS],
            "G_p_consistency": {},  # one-off check that G_p formula matches Ramanujan-Fourier
            "primorial_W": {},      # F1: T_W^{div} triple correlations
            "general_Q":   {},      # F2: T_Q triple correlations
            "prime_triple_density": {},  # F6 baseline
            "HL_full":     {},      # full HL prediction
        }

        # G_p consistency check (algebraic, no N-dependence required):
        # for several primes and shifts, verify G_p(p, h1, h2) = G_p_via_f.
        print("  G_p formula consistency check (closed form vs Ramanujan-Fourier f_p)...", flush=True)
        gp_cells = []
        for p in [2, 3, 5, 7, 11]:
            for (h1, h2) in H_PAIRS:
                if h1 == 0 and h2 == 0:
                    continue
                a = G_p(p, h1, h2)
                b = G_p_via_f(p, h1, h2)
                gp_cells.append({
                    "p": p,
                    "h1": int(h1),
                    "h2": int(h2),
                    "G_closed_form": a,
                    "G_via_f": b,
                    "match": abs(a - b) < 1e-12,
                })
        out_d["G_p_consistency"] = gp_cells
        bad = [c for c in gp_cells if not c["match"]]
        if bad:
            print(f"  *** G_p consistency FAILED for {len(bad)} cells ***")
        else:
            print(f"  G_p consistency: {len(gp_cells)} cells OK", flush=True)

        # Full HL once per (h1, h2):
        print("  computing full HL triple singular series (truncated at p<=50000)...", flush=True)
        for (h1, h2) in H_PAIRS:
            if (h1, h2) == (0, 0):
                out_d["HL_full"][f"{h1},{h2}"] = None
                continue
            S_full = full_HL_triple(h1, h2, primes_for_C3)
            out_d["HL_full"][f"{h1},{h2}"] = S_full

        # Prime triple density baseline:
        print("  prime triple density baseline...", flush=True)
        for (h1, h2) in H_PAIRS:
            pi_3 = prime_triple_density(chi_P[:N], h1, h2)
            out_d["prime_triple_density"][f"{h1},{h2}"] = {
                "pi_3": pi_3,
                "rho3_HL": density ** 3 * out_d["HL_full"][f"{h1},{h2}"]
                            if out_d["HL_full"][f"{h1},{h2}"] is not None
                            and out_d["HL_full"][f"{h1},{h2}"] not in (float('inf'),)
                            else None,
            }

        # F1: T_W^{div} for primorial W -> exact prime-by-prime truncation
        print("  F1: T_W^{div} primorial-W triple correlations...", flush=True)
        for W in W_VALS:
            tw0 = time.time()
            T_W = build_T_W_div(N, W, density)
            mean_T = float(np.mean(T_W))
            tw1 = time.time()
            cells = {}
            # primes dividing W:
            primes_W = []
            n = W
            p = 2
            while p * p <= n:
                if n % p == 0:
                    primes_W.append(p)
                    while n % p == 0:
                        n //= p
                p += 1
            if n > 1:
                primes_W.append(n)
            for (h1, h2) in H_PAIRS:
                R = correlation_triple(T_W, h1, h2)
                pred = density ** 3 * trunc_HL_triple_via_primes(primes_W, h1, h2)
                ratio = R / pred if pred != 0 else None
                cells[f"{h1},{h2}"] = {
                    "R": R,
                    "pred": pred,
                    "ratio": ratio,
                }
            out_d["primorial_W"][W] = {
                "primes_W": primes_W,
                "mean_T": mean_T,
                "build_time_s": tw1 - tw0,
                "cells": cells,
            }
            print(f"    W = {W} (primes {primes_W}): built in {tw1-tw0:.1f}s", flush=True)

        # F2: T_Q for general Q
        print("  F2: T_Q general-Q triple correlations...", flush=True)
        for Q in Q_VALS:
            print(f"    Q = {Q}: building T_Q...", flush=True)
            tq0 = time.time()
            T_Q = build_T_Q(N, Q, mu, phi, density)
            mean_T_Q = float(np.mean(T_Q))
            tq1 = time.time()
            sqfree_primes_Q = squarefree_primes_up_to(Q, mu)
            cells = {}
            for (h1, h2) in H_PAIRS:
                R = correlation_triple(T_Q, h1, h2)
                pred = density ** 3 * trunc_HL_triple_via_primes(sqfree_primes_Q, h1, h2)
                ratio = R / pred if pred != 0 else None
                cells[f"{h1},{h2}"] = {
                    "R": R,
                    "pred": pred,
                    "ratio": ratio,
                }
            out_d["general_Q"][Q] = {
                "primes_Q": sqfree_primes_Q,
                "mean_T": mean_T_Q,
                "build_time_s": tq1 - tq0,
                "cells": cells,
            }

        out[d] = out_d
        print(f"  d = {d} done in {time.time() - t0:.1f}s", flush=True)

    # Save raw results.
    with open("tq_triple_correlation_results.json", "w") as f:
        json.dump(out, f, indent=2, default=lambda x: float(x) if isinstance(x, np.floating) else x)
    print("\nWrote tq_triple_correlation_results.json")

    # Print summary tables.
    print("\n\n=== SUMMARY ===")
    for d in D_LIST:
        d_data = out[d]
        print(f"\n# d = {d}, N = 2^{d}, pi(N)/N = {d_data['density']:.6f}\n")

        # G_p consistency
        cells = d_data["G_p_consistency"]
        bad = [c for c in cells if not c["match"]]
        print(f"  G_p formula consistency: {len(cells)} cells, {len(bad)} mismatches.")

        # F1: primorial W
        print(f"\n  F1 (primorial W, T_W^{{div}}): R / [(pi/N)^3 prod_{{p|W}} G_p]:")
        header = f"  {'(h1,h2)':>10}" + "".join(f"{f'W={W}':>11}" for W in d_data["W_vals"])
        print(header)
        for (h1, h2) in H_PAIRS:
            row = f"  ({h1:>3},{h2:>3})"
            for W in d_data["W_vals"]:
                cell = d_data["primorial_W"][W]["cells"][f"{h1},{h2}"]
                if cell["ratio"] is None:
                    row += f"{'-':>11}"
                else:
                    row += f"{cell['ratio']:>11.5f}"
            print(row)

        # F2: general Q
        print(f"\n  F2 (general Q, T_Q): R / [(pi/N)^3 prod_{{p sqf <= Q}} G_p]:")
        header = f"  {'(h1,h2)':>10}" + "".join(f"{f'Q={Q}':>11}" for Q in d_data["Q_vals"])
        print(header)
        for (h1, h2) in H_PAIRS:
            row = f"  ({h1:>3},{h2:>3})"
            for Q in d_data["Q_vals"]:
                cell = d_data["general_Q"][Q]["cells"][f"{h1},{h2}"]
                if cell["ratio"] is None:
                    row += f"{'-':>11}"
                else:
                    row += f"{cell['ratio']:>11.5f}"
            print(row)

        # F3: HL recovery at largest Q
        Q_max = max(d_data["Q_vals"])
        print(f"\n  F3 (HL recovery, Q = {Q_max} ~ N^0.5): R / [(pi/N)^3 S_HL(0,h1,h2)]:")
        for (h1, h2) in H_PAIRS:
            S_HL = d_data["HL_full"][f"{h1},{h2}"]
            if S_HL is None:
                continue
            R = d_data["general_Q"][Q_max]["cells"][f"{h1},{h2}"]["R"]
            pred = d_data["density"] ** 3 * S_HL
            ratio = R / pred if pred != 0 else None
            ratio_s = f"{ratio:.5f}" if ratio is not None else "-"
            print(f"    (h1,h2) = ({h1:>3},{h2:>3}): R = {R:.6e}, pred_HL = {pred:.6e}, ratio = {ratio_s}")

        # F6: prime triple density
        print(f"\n  F6 (prime triple density baseline):")
        for (h1, h2) in H_PAIRS:
            pi_3 = d_data["prime_triple_density"][f"{h1},{h2}"]["pi_3"]
            rho3_HL = d_data["prime_triple_density"][f"{h1},{h2}"]["rho3_HL"]
            ratio = pi_3 / rho3_HL if rho3_HL not in (None, 0, 0.0) else None
            ratio_s = f"{ratio:.4f}" if ratio is not None else "-"
            HL_s = f"{rho3_HL:.6e}" if rho3_HL is not None else "(degenerate)"
            print(f"    (h1,h2) = ({h1:>3},{h2:>3}): pi_3 = {pi_3:.6e}, rho^3 S_HL = {HL_s}, ratio = {ratio_s}")


if __name__ == "__main__":
    main()
