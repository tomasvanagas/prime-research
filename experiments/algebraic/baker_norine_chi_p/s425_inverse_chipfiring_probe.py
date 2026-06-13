"""
S425 — Re-verify-closure of E2.28 (Baker-Norine q-reduced form of D_P^N).

Adversarial frame: the S161 closure dismissed the chip-firing route as
"NOT algorithmically useful because building D_P^N requires knowing
primes." This script probes that claim by testing whether ANY
polylog-computable INPUT divisor on Γ_N or H_N has q-reduced form that
encodes π(N) (i.e., extracts π(N) without prime input).

Probes:
  (A) D = K · δ_1: trivially q-reduced -> K · δ_1.
  (B) D = Σ_n δ_n (uniform): q-reduces to N · δ_1 on Γ_N.
  (C) D = Σ_n τ(n) · δ_n (divisor function, polylog-computable per n).
  (D) D = Σ_n μ²(n) · δ_n (squarefree indicator, polylog).
  (E) D = Σ_n Λ(n) · δ_n via "is-prime-power-on-Γ_N-without-prime-input":
       fall back to D_P + D_pp (we DO need primes to construct this;
       included as a control to confirm the identity).

Verdict pattern:
  - Each non-prime divisor's q-reduced form has D'(1) equal to the
    SUM of input chips (conservation law) plus or minus the chip-firing
    redistribution, NOT π(N) in general.
  - π(N) emerges only from D_P (or arithmetically derived divisors that
    already encode primes).

Run: python3 s425_inverse_chipfiring_probe.py --N 32 64 128
"""

import argparse
import json
import time
from math import gcd

from baker_norine_chi_p import (
    divisibility_graph,
    q_reduce,
)
from sympy import isprime, divisor_count, mobius, factorint


def hasse_cover_graph(N):
    """V = [1, N], edges (a, p*a) for primes p with p*a <= N. Mirrors
    the S161 H_N construction so we can confirm the D_P -> pi(N) identity."""
    V = list(range(1, N + 1))
    adj = {v: [] for v in V}
    primes_le_N = [p for p in range(2, N + 1) if isprime(p)]
    for a in range(1, N + 1):
        for p in primes_le_N:
            b = p * a
            if b > N:
                break
            adj[a].append(b)
            adj[b].append(a)
    return V, adj


def divisor_DP(N):
    return {n: (1 if isprime(n) else 0) for n in range(1, N + 1)}


def divisor_uniform(N):
    return {n: 1 for n in range(1, N + 1)}


def divisor_sink_only(N, K):
    return {n: (K if n == 1 else 0) for n in range(1, N + 1)}


def divisor_tau(N):
    return {n: divisor_count(n) for n in range(1, N + 1)}


def divisor_mu_sq(N):
    return {n: int(mobius(n) != 0) for n in range(1, N + 1)}


def divisor_omega1(N):
    """1 if Omega(n) == 1 (prime powers), 0 else. Constructing this
    requires knowing prime powers, which is essentially knowing primes
    (a prime power is k = p^j; identifying p requires factoring, polylog
    via BPSW + perfect-power test). Included as a 'closer to prime
    indicator' control."""
    out = {}
    for n in range(1, N + 1):
        if n == 1:
            out[n] = 0
            continue
        f = factorint(n)
        out[n] = 1 if len(f) == 1 else 0
    return out


PROBES = [
    ("D_P (prime indicator)", divisor_DP),
    ("D_uniform (all-ones)", divisor_uniform),
    ("D_sink_K=N (chips at vertex 1 only)", lambda N: divisor_sink_only(N, N)),
    ("D_tau (divisor function)", divisor_tau),
    ("D_mu_sq (squarefree indicator)", divisor_mu_sq),
    ("D_omega1 (prime-power indicator)", divisor_omega1),
]


def run_probe(N, graph_name, V, adj):
    rows = []
    for label, fn in PROBES:
        D = fn(N)
        total = sum(D.values())
        t0 = time.time()
        try:
            D_red = q_reduce(D, q=1, V=V, adj=adj)
            wall = time.time() - t0
            sink = D_red[1]
            nonsink_support = [v for v, c in D_red.items() if v != 1 and c > 0]
            row = {
                "graph": graph_name,
                "divisor": label,
                "N": N,
                "input_total": total,
                "qreduced_sink_chips": sink,
                "qreduced_nonsink_support_size": len(nonsink_support),
                "qreduced_nonsink_total": total - sink,
                "wall_s": round(wall, 3),
            }
            print(
                f"  [{graph_name} N={N}] {label}: "
                f"input_total={total} -> sink={sink}, "
                f"non-sink chips on {len(nonsink_support)} vertices "
                f"(total {total - sink}), wall={wall:.2f}s"
            )
        except Exception as e:
            row = {
                "graph": graph_name,
                "divisor": label,
                "N": N,
                "error": str(e),
            }
            print(f"  [{graph_name} N={N}] {label}: ERROR {e}")
        rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", nargs="+", type=int, default=[16, 32, 64])
    ap.add_argument("--out", default="s425_inverse_chipfiring_results.json")
    args = ap.parse_args()

    all_rows = []
    for N in args.N:
        pi_N = sum(1 for k in range(2, N + 1) if isprime(k))
        pi_half = sum(1 for k in range(2, N // 2 + 1) if isprime(k))
        print(f"\n=== N = {N} (pi(N)={pi_N}, pi(N/2)={pi_half}) ===")

        # Gamma_N (full divisibility)
        V, adj = divisibility_graph(N)
        m = sum(len(adj[v]) for v in V) // 2
        print(f"Gamma_N: |V|={len(V)}, |E|={m}, genus={m - len(V) + 1}")
        all_rows.extend(run_probe(N, f"Gamma_{N}", V, adj))

        # H_N (Hasse cover) -- expensive to build at large N because we
        # need primes; this is meta-allowed since we are checking the
        # known D_P -> pi(N) identity, not seeking a polylog algorithm.
        Vh, adjh = hasse_cover_graph(N)
        mh = sum(len(adjh[v]) for v in Vh) // 2
        print(f"H_{N}: |V|={len(Vh)}, |E|={mh}, genus={mh - len(Vh) + 1}")
        all_rows.extend(run_probe(N, f"H_{N}", Vh, adjh))

    # Coerce sympy Integer to int for JSON
    def _coerce(o):
        if isinstance(o, dict):
            return {k: _coerce(v) for k, v in o.items()}
        if isinstance(o, list):
            return [_coerce(v) for v in o]
        try:
            return int(o)
        except (TypeError, ValueError):
            return o

    with open(args.out, "w") as f:
        json.dump([_coerce(r) for r in all_rows], f, indent=2)
    print(f"\nWrote {len(all_rows)} rows to {args.out}")


if __name__ == "__main__":
    main()
