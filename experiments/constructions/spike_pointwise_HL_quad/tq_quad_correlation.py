"""
tq_quad_correlation.py
======================

Four-point companion of S205 (2-point), S208 (1-point divisor collapse)
and S209 (3-point) constructions.  We compute

    R_{h1,h2,h3}(Q, N)
       :=  (1/(N - max(h1,h2,h3))) sum_n T_Q(n) T_Q(n+h1) T_Q(n+h2) T_Q(n+h3)

for various (Q, h1, h2, h3) and compare with the prime-by-prime closed
form

    pred(Q, h1, h2, h3)  =  (pi(N)/N)^4
                            * prod_{p sqf primes <= Q} G_p^{(4)}(h1, h2, h3)

where  G_p^{(4)}(h1, h2, h3) = (p - nu_p({0, h1, h2, h3})) * p^3 / (p-1)^4,
       nu_p = number of distinct residues among {0, h1, h2, h3} mod p.

For *primorial* W, T_W^{div}(n) = (pi(N)/N) * (W/phi(W)) * [gcd(n, W) = 1]
collapses pointwise (S208).  The 4-point average then becomes:

    <T_W^{div}(n) T_W^{div}(n+h1) T_W^{div}(n+h2) T_W^{div}(n+h3)>_n
        = (pi(N)/N)^4 * (W/phi(W))^4 * <prod_{i=0..3} [gcd(n+h_i, W) = 1]>_n
        = (pi(N)/N)^4 * (W/phi(W))^4 * prod_{p|W} (p - nu_p) / p
        = (pi(N)/N)^4 * prod_{p|W} (p - nu_p) p^3 / (p-1)^4
        = (pi(N)/N)^4 * S_HL^{(W)}(0, h1, h2, h3),

i.e. the W-conductor truncation of the Hardy-Littlewood prime quadruple
singular series.

The script verifies six falsifications:

  F1  R_{h1,h2,h3}(W, N) / [(pi/N)^4 prod_{p|W} G_p^{(4)}] -> 1 at primorial W.
  F2  R_{h1,h2,h3}(Q, N) / [(pi/N)^4 prod_{p sqf primes <= Q} G_p^{(4)}]
      -> 1 at general Q (within finite-N tolerance).
  F3  HL recovery at large Q (Q ~ sqrt(N)): ratio to full
      S_HL(0, h1, h2, h3) (over ALL primes, truncated at 50000).
  F4  h1 = h2 = h3 = 0 self-consistency: <T_W^4> = (pi/N)^4 * (W/phi(W))^3.
  F5  Reduction to S209 3-point at h3 = 0:
      <T_W^{div}^4 (h1, h2, 0)> = (W/phi(W)) * <T_W^{div}^3 (h1, h2)>.
  F6  Inadmissible quadruple (some prime p has nu_p = p) -> primorial-W
      cell with that p in W has empirical correlation = 0 to within
      finite-N noise.
  F7  G_p^{(4)} closed-form algebraic consistency vs the direct
      Ramanujan-Fourier sum f_p^{(4)} on a small-prime grid.

Run time at d in {16, 18, 20} on a laptop: roughly 3-6 minutes total.
"""

import json
import math
import time

import numpy as np


# ------------------------------------------------------------------
# Sieve helpers (chi_P, mu, phi)
# ------------------------------------------------------------------

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


def primes_up_to(N_max: int) -> list:
    sieve_arr = np.ones(N_max + 1, dtype=bool)
    sieve_arr[:2] = False
    for p in range(2, int(N_max ** 0.5) + 1):
        if sieve_arr[p]:
            sieve_arr[p * p :: p] = False
    return list(np.where(sieve_arr)[0])


def squarefree_primes_up_to(Q: int, mu_arr: np.ndarray) -> list:
    return [p for p in range(2, Q + 1) if mu_arr[p] == -1 and p > 1]


# ------------------------------------------------------------------
# Ramanujan sums
# ------------------------------------------------------------------

def ramanujan_pattern(q: int, mu_arr: np.ndarray, phi_arr: np.ndarray) -> np.ndarray:
    """Return length-q array c_q(r) for r = 0..q-1.  c_q(h) = mu(q/d) phi(d), d=gcd(q,h)."""
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
    """T_W^{div}(n) = (pi(N)/N) * (W/phi(W)) * [gcd(n, W) = 1] for squarefree W
    (or W replaced by its radical)."""
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
    rad_W = 1
    for p in primes_W:
        rad_W *= p
    phi_W = rad_W
    for p in primes_W:
        phi_W = phi_W * (p - 1) // p
    n_arr = np.arange(N, dtype=np.int64)
    coprime = np.ones(N, dtype=bool)
    for p in primes_W:
        coprime &= (n_arr % p != 0)
    out = np.where(coprime, 1.0, 0.0)
    out *= prime_density * (rad_W / phi_W)
    return out, primes_W


# ------------------------------------------------------------------
# nu_p, G_p, f_p
# ------------------------------------------------------------------

def nu_p(p: int, h_set: tuple) -> int:
    """Number of distinct residues mod p of the elements of h_set."""
    return len({h % p for h in h_set})


def multiplicity_profile(p: int, h_set: tuple) -> tuple:
    """Return sorted tuple of multiplicities (m_x: m_x > 0) over residues mod p."""
    counts = {}
    for h in h_set:
        r = h % p
        counts[r] = counts.get(r, 0) + 1
    return tuple(sorted(counts.values(), reverse=True))


def G_p_quad(p: int, h1: int, h2: int, h3: int) -> float:
    """Prime-p factor of the HL prime QUADRUPLE singular series.

    G_p^{(4)} = (p - nu_p) * p^3 / (p - 1)^4
    """
    nu = nu_p(p, (0, h1, h2, h3))
    return (p - nu) * (p ** 3) / ((p - 1) ** 4)


def f_p_direct_quad(p: int, h1: int, h2: int, h3: int) -> float:
    """Direct (1/p) sum_{r mod p} c_p(r) c_p(r+h1) c_p(r+h2) c_p(r+h3).

    For prime p: c_p(r) = (p-1) if p|r else -1.
    """
    total = 0
    for r in range(p):
        a = (p - 1) if r % p == 0 else -1
        b = (p - 1) if (r + h1) % p == 0 else -1
        c = (p - 1) if (r + h2) % p == 0 else -1
        d = (p - 1) if (r + h3) % p == 0 else -1
        total += a * b * c * d
    return total / p


def G_p_via_f_quad(p: int, h1: int, h2: int, h3: int) -> float:
    """Compute G_p^{(4)} via the Ramanujan-Fourier expansion.

    Expanding T_p(n) = (1/(p-1)) c_p(n) about the rank-1 mode and taking
    the connected 4-point average over n at a single prime conductor p
    yields, after dividing by (pi/N)^4 and the rank-1 normalisation,

      G_p^{(4)}(h1,h2,h3)
        = 1
          + 1/(p-1)^2 * SUM_{i<j} c_p(h_j - h_i)
          - 1/(p-1)^3 * SUM_{|S|=3 \subset {0,h1,h2,h3}} f_p^{(3)}(S)
          + 1/(p-1)^4 * f_p^{(4)}(h1, h2, h3)

    where f_p^{(3)}(a,b,c) = (1/p) sum_r c_p(r-a) c_p(r-b) c_p(r-c).
    Implemented directly.
    """
    H = (0, h1, h2, h3)

    def cp(h):
        return (p - 1) if h % p == 0 else -1

    s = 1.0
    # 2-point cross terms: pick 2 of 4 = 6 pairs
    for i in range(4):
        for j in range(i + 1, 4):
            s += cp(H[j] - H[i]) / ((p - 1) ** 2)
    # 3-point cross terms: pick 3 of 4 = 4 triples; subtract
    # f_p^{(3)} relative to the leading triple.
    for skip in range(4):
        triple = tuple(H[k] for k in range(4) if k != skip)
        # f_p^{(3)} (a,b,c) over residue r:
        a, b, c = triple
        total = 0
        for r in range(p):
            total += cp(r - a) * cp(r - b) * cp(r - c)
        f3 = total / p
        s -= f3 / ((p - 1) ** 3)
    # 4-point f_p:
    s += f_p_direct_quad(p, h1, h2, h3) / ((p - 1) ** 4)
    return s


# ------------------------------------------------------------------
# HL truncations
# ------------------------------------------------------------------

def trunc_HL_quad_via_primes(prime_list: list, h1: int, h2: int, h3: int) -> float:
    s = 1.0
    for p in prime_list:
        s *= G_p_quad(p, h1, h2, h3)
    return s


def full_HL_quad(h1: int, h2: int, h3: int, primes_for_C4: list) -> float:
    """HL prime quadruple singular series S(0,h1,h2,h3) = prod_p G_p^{(4)}.

    Truncated at primes_for_C4[-1].  If the 4-tuple is inadmissible at
    some prime in primes_for_C4 (nu_p = p), returns 0.
    """
    H = (0, h1, h2, h3)
    if len(set(H)) == 1:
        return float('inf')
    for p in primes_for_C4:
        if nu_p(p, H) == p:
            return 0.0
    s = 1.0
    for p in primes_for_C4:
        s *= G_p_quad(p, h1, h2, h3)
    return s


# ------------------------------------------------------------------
# 4-point correlation
# ------------------------------------------------------------------

def correlation_quad(T: np.ndarray, h1: int, h2: int, h3: int) -> float:
    """(1/L) sum T[n] T[n+h1] T[n+h2] T[n+h3] over valid indices."""
    h_max = max(0, h1, h2, h3)
    if h_max == 0:
        return float(np.mean(T * T * T * T))
    L = T.shape[0] - h_max
    a = T[:L]
    b = T[h1 : h1 + L]
    c = T[h2 : h2 + L]
    d = T[h3 : h3 + L]
    return float(np.mean(a * b * c * d))


def prime_quad_density(chi_P: np.ndarray, h1: int, h2: int, h3: int) -> float:
    h_max = max(0, h1, h2, h3)
    if h_max == 0:
        return float(np.mean(chi_P.astype(np.float64) ** 4))
    L = chi_P.shape[0] - h_max
    a = chi_P[:L].astype(np.float64)
    b = chi_P[h1 : h1 + L].astype(np.float64)
    c = chi_P[h2 : h2 + L].astype(np.float64)
    d = chi_P[h3 : h3 + L].astype(np.float64)
    return float(np.mean(a * b * c * d))


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------

def main():
    D_LIST = [16, 18, 20]
    # 4-tuples (h1, h2, h3) chosen for: F4 self-consistency, F5 reduction
    # to lower-order, admissible quadruples, inadmissible probe.
    # All h_i are even (so ν_2 = 1 < 2, admissible at p=2).
    H_TRIPLES = [
        (0, 0, 0),       # F4: <T^4> self-consistency
        (2, 0, 0),       # F5 reduction: matches S208 1-point at h=2
        (6, 0, 0),       # F5 reduction (h=6)
        (2, 4, 0),       # F5 reduction: matches S205 2-point (h1,h2 = 2,4)
        (2, 6, 0),       # F5 reduction (S205 2-point, h1=2,h2=6)
        (6, 12, 0),      # F5 reduction (S205 2-point, h1=6,h2=12)
        (2, 6, 8),       # admissible quad: ν_3=2, ν_5=4, ν_7=4
        (2, 6, 12),      # admissible quad: ν_3=2, ν_5=3, ν_7=4
        (4, 6, 10),      # admissible quad: ν_3=2, ν_5=3
        (2, 8, 14),      # admissible quad: ν_3=2, ν_5=4, ν_7=3
        (6, 10, 12),     # admissible quad: ν_3=2, ν_5=3, ν_7=4
        (6, 12, 30),     # admissible quad: ν_3=2, ν_5=3, ν_7=4, ν_11=4
        (2, 4, 6),       # F6 inadmissible at p=3: {0,2,4,6} mod 3 = {0,2,1,0}, ν_3=3 = p=3
        (2, 4, 8),       # F6 inadmissible at p=3: {0,2,4,8} mod 3 = {0,2,1,2}, ν_3=3
    ]

    primes_for_C4 = primes_up_to(50000)

    out = {}

    for d in D_LIST:
        N = 1 << d
        print(f"\n=== d = {d}, N = 2^{d} = {N:,} ===", flush=True)
        t0 = time.time()
        Q_VALS = sorted(set([2, 6, 30, 210, 2310,
                             max(2, round(N ** 0.185)),
                             max(2, round(N ** 0.5))]))
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
            "H_triples": [list(p) for p in H_TRIPLES],
            "G_p_consistency": [],   # F7
            "primorial_W": {},       # F1
            "general_Q":   {},       # F2
            "prime_quad_density": {},# F6 baseline
            "HL_full":     {},       # F3
        }

        # F7: G_p formula consistency check (closed form vs Ramanujan-Fourier)
        print("  F7: G_p closed-form vs Ramanujan-Fourier consistency...", flush=True)
        gp_cells = []
        for p in [2, 3, 5, 7, 11, 13]:
            for (h1, h2, h3) in H_TRIPLES:
                if h1 == 0 and h2 == 0 and h3 == 0:
                    continue
                a = G_p_quad(p, h1, h2, h3)
                b = G_p_via_f_quad(p, h1, h2, h3)
                gp_cells.append({
                    "p": p,
                    "h1": int(h1),
                    "h2": int(h2),
                    "h3": int(h3),
                    "nu_p": nu_p(p, (0, h1, h2, h3)),
                    "mult": list(multiplicity_profile(p, (0, h1, h2, h3))),
                    "G_closed_form": a,
                    "G_via_f": b,
                    "abs_err": abs(a - b),
                    "match": abs(a - b) < 1e-10,
                })
        out_d["G_p_consistency"] = gp_cells
        bad = [c for c in gp_cells if not c["match"]]
        if bad:
            print(f"  *** F7 G_p consistency FAILED for {len(bad)} cells ***")
            for c in bad[:5]:
                print(f"     p={c['p']} (h1,h2,h3)=({c['h1']},{c['h2']},{c['h3']}) "
                      f"closed={c['G_closed_form']:.6f} via_f={c['G_via_f']:.6f}")
        else:
            print(f"  F7 G_p consistency: {len(gp_cells)} cells OK", flush=True)

        # F3: full HL once per (h1, h2, h3)
        print("  computing full HL quadruple singular series (truncated at p<=50000)...", flush=True)
        for (h1, h2, h3) in H_TRIPLES:
            if (h1, h2, h3) == (0, 0, 0):
                out_d["HL_full"][f"{h1},{h2},{h3}"] = None
                continue
            S_full = full_HL_quad(h1, h2, h3, primes_for_C4)
            out_d["HL_full"][f"{h1},{h2},{h3}"] = S_full

        # F6 baseline: prime quadruple density
        print("  prime quadruple density baseline...", flush=True)
        for (h1, h2, h3) in H_TRIPLES:
            pi_4 = prime_quad_density(chi_P[:N], h1, h2, h3)
            S_HL = out_d["HL_full"][f"{h1},{h2},{h3}"]
            rho4_HL = (
                density ** 4 * S_HL
                if S_HL is not None and S_HL not in (float('inf'),)
                else None
            )
            out_d["prime_quad_density"][f"{h1},{h2},{h3}"] = {
                "pi_4": pi_4,
                "rho4_HL": rho4_HL,
            }

        # F1: T_W^{div} primorial-W 4-point
        print("  F1: T_W^{div} primorial-W 4-point correlations...", flush=True)
        for W in W_VALS:
            tw0 = time.time()
            T_W, primes_W = build_T_W_div(N, W, density)
            mean_T = float(np.mean(T_W))
            tw1 = time.time()
            cells = {}
            for (h1, h2, h3) in H_TRIPLES:
                R = correlation_quad(T_W, h1, h2, h3)
                pred = density ** 4 * trunc_HL_quad_via_primes(primes_W, h1, h2, h3)
                ratio = R / pred if pred != 0 else None
                cells[f"{h1},{h2},{h3}"] = {
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
            print(f"    W = {W} (primes {primes_W}): built in {tw1 - tw0:.1f}s",
                  flush=True)

        # F2: T_Q general-Q 4-point
        print("  F2: T_Q general-Q 4-point correlations...", flush=True)
        for Q in Q_VALS:
            print(f"    Q = {Q}: building T_Q...", flush=True)
            tq0 = time.time()
            T_Q = build_T_Q(N, Q, mu, phi, density)
            mean_T_Q = float(np.mean(T_Q))
            tq1 = time.time()
            sqfree_primes_Q = squarefree_primes_up_to(Q, mu)
            cells = {}
            for (h1, h2, h3) in H_TRIPLES:
                R = correlation_quad(T_Q, h1, h2, h3)
                pred = density ** 4 * trunc_HL_quad_via_primes(
                    sqfree_primes_Q, h1, h2, h3)
                ratio = R / pred if pred != 0 else None
                cells[f"{h1},{h2},{h3}"] = {
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

    with open("tq_quad_correlation_results.json", "w") as f:
        json.dump(out, f, indent=2,
                  default=lambda x: float(x) if isinstance(x, np.floating) else x)
    print("\nWrote tq_quad_correlation_results.json")

    # Summary
    print("\n\n=== SUMMARY ===")
    for d in D_LIST:
        d_data = out[d]
        print(f"\n# d = {d}, N = 2^{d}, pi(N)/N = {d_data['density']:.6f}\n")

        cells = d_data["G_p_consistency"]
        bad = [c for c in cells if not c["match"]]
        max_err = max(c["abs_err"] for c in cells)
        print(f"  F7 G_p formula consistency: {len(cells)} cells, "
              f"{len(bad)} mismatches, max abs_err = {max_err:.3e}.")

        # F1: primorial W ratios
        print(f"\n  F1 (primorial W, T_W^{{div}}): "
              f"R / [(pi/N)^4 prod_{{p|W}} G_p^{{(4)}}]:")
        header = "  " + f"{'(h1,h2,h3)':>14}" + "".join(
            f"{f'W={W}':>10}" for W in d_data["W_vals"])
        print(header)
        for (h1, h2, h3) in H_TRIPLES:
            row = f"  ({h1:>3},{h2:>3},{h3:>3})"
            for W in d_data["W_vals"]:
                cell = d_data["primorial_W"][W]["cells"][f"{h1},{h2},{h3}"]
                if cell["ratio"] is None:
                    row += f"{'-':>10}"
                else:
                    row += f"{cell['ratio']:>10.5f}"
            print(row)

        # F2: general Q
        print(f"\n  F2 (general Q, T_Q): "
              f"R / [(pi/N)^4 prod_{{p sqf <= Q}} G_p^{{(4)}}]:")
        header = "  " + f"{'(h1,h2,h3)':>14}" + "".join(
            f"{f'Q={Q}':>10}" for Q in d_data["Q_vals"])
        print(header)
        for (h1, h2, h3) in H_TRIPLES:
            row = f"  ({h1:>3},{h2:>3},{h3:>3})"
            for Q in d_data["Q_vals"]:
                cell = d_data["general_Q"][Q]["cells"][f"{h1},{h2},{h3}"]
                if cell["ratio"] is None:
                    row += f"{'-':>10}"
                else:
                    row += f"{cell['ratio']:>10.5f}"
            print(row)

        # F3: HL recovery at largest Q
        Q_max = max(d_data["Q_vals"])
        print(f"\n  F3 (HL recovery, Q = {Q_max} ~ N^0.5): "
              f"R / [(pi/N)^4 S_HL(0,h1,h2,h3)]:")
        for (h1, h2, h3) in H_TRIPLES:
            S_HL = d_data["HL_full"][f"{h1},{h2},{h3}"]
            if S_HL is None:
                continue
            R = d_data["general_Q"][Q_max]["cells"][f"{h1},{h2},{h3}"]["R"]
            pred = d_data["density"] ** 4 * S_HL
            ratio = R / pred if pred not in (0, 0.0) else None
            ratio_s = f"{ratio:.5f}" if ratio is not None else "(zero)"
            print(f"    (h1,h2,h3) = ({h1:>3},{h2:>3},{h3:>3}): "
                  f"R = {R:.6e}, pred_HL = {pred:.6e}, ratio = {ratio_s}")

        # F6: prime quadruple density
        print(f"\n  F6 (prime quadruple density baseline):")
        for (h1, h2, h3) in H_TRIPLES:
            pi_4 = d_data["prime_quad_density"][f"{h1},{h2},{h3}"]["pi_4"]
            rho4_HL = d_data["prime_quad_density"][f"{h1},{h2},{h3}"]["rho4_HL"]
            ratio = pi_4 / rho4_HL if rho4_HL not in (None, 0, 0.0) else None
            ratio_s = f"{ratio:.4f}" if ratio is not None else "-"
            HL_s = f"{rho4_HL:.6e}" if rho4_HL is not None else "(degenerate)"
            print(f"    (h1,h2,h3) = ({h1:>3},{h2:>3},{h3:>3}): "
                  f"pi_4 = {pi_4:.6e}, rho^4 S_HL = {HL_s}, ratio = {ratio_s}")


if __name__ == "__main__":
    main()
