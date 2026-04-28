"""
Richer signal extraction beyond rank: for each divisor D on graph G, compute:

1. q-reduced form D' — full chip distribution.
2. D'(q) margin — "winnability slack" at the sink.
3. The winnable-subtraction set W(D) := {v in V : D - delta_v is winnable}.
   For r=0, W(D) is a proper subset of V; its cardinality and structure
   provide a refined invariant.
4. A 2-step extension: W^(2)(D) := pairs (u,v) such that D - delta_u - delta_v
   is winnable; gives finer rank info.

Compares prime divisor signature to random matched-degree divisors.
"""

import sys
import json
import time
import argparse
import random
import statistics
from math import gcd
from collections import deque, Counter
from sympy import isprime

from baker_norine_chi_p import (
    divisibility_graph, coprimality_graph, graph_stats, is_connected,
    q_reduce, is_winnable, prime_divisor, random_matched_divisor
)


def winnable_subtraction_set(D, q, V, adj):
    """Return {v in V : D - delta_v is winnable}."""
    W = []
    for v in V:
        D_minus = dict(D)
        D_minus[v] -= 1
        if is_winnable(D_minus, q, V, adj):
            W.append(v)
    return W


def reduced_form_signature(D, q, V, adj):
    """Compute q-reduced form D' and extract signature features."""
    D_red = q_reduce(D, q, V, adj)
    return {
        "D_red": dict(D_red),
        "D_red_at_q": D_red[q],
        "support": [v for v in V if v != q and D_red[v] > 0],
        "max_chip": max(D_red[v] for v in V if v != q),
        "total_nonq_chips": sum(D_red[v] for v in V if v != q),
    }


def signature(D, q, V, adj):
    """Return full signature."""
    sig = reduced_form_signature(D, q, V, adj)
    sig["winnable_subtraction_set"] = winnable_subtraction_set(D, q, V, adj)
    sig["W_size"] = len(sig["winnable_subtraction_set"])
    return sig


def classify_vertices(V):
    """Tag each v with (is_prime, is_prime_power, omega(v))."""
    from sympy import factorint
    cls = {}
    for v in V:
        if v == 1:
            cls[v] = ("one", 1, 0)
            continue
        f = factorint(v)
        cls[v] = (
            "prime" if isprime(v) else "composite",
            "pp" if len(f) == 1 else "non_pp",
            sum(f.values()),  # Omega(v)
        )
    return cls


def run_signature_experiment(N, graph_type, n_random, q_choice, seed):
    print(f"\n=== N={N}, graph={graph_type}, q={q_choice} ===", flush=True)
    if graph_type == "div":
        V, adj = divisibility_graph(N)
    elif graph_type == "coprime":
        V, adj = coprimality_graph(N)
    else:
        raise ValueError(graph_type)
    n, m, g = graph_stats(V, adj)
    if not is_connected(V, adj):
        return {"N": N, "graph_type": graph_type, "skipped": "disconnected"}

    if q_choice == "low":
        q = V[0]
    elif q_choice == "max_deg":
        q = max(V, key=lambda v: len(adj[v]))
    else:
        q = q_choice

    primes = [v for v in V if isprime(v)]
    deg_DP = len(primes)
    print(f"  |V|={n}, |E|={m}, g={g}, deg(D_P)={deg_DP}, q={q}", flush=True)

    # Prime divisor signature
    D_P = prime_divisor(V)
    t0 = time.time()
    sig_P = signature(D_P, q, V, adj)
    t_P = time.time() - t0
    print(f"  D_P: D'(q)={sig_P['D_red_at_q']}, |supp|={len(sig_P['support'])}, "
          f"max_chip={sig_P['max_chip']}, |W|={sig_P['W_size']} ({t_P:.1f}s)", flush=True)

    # Random matched
    rng = random.Random(seed)
    randoms = []
    for i in range(n_random):
        D_R = random_matched_divisor(V, deg_DP, rng)
        t0 = time.time()
        sig_R = signature(D_R, q, V, adj)
        t_R = time.time() - t0
        randoms.append(sig_R)
        print(f"    rand {i+1}/{n_random}: D'(q)={sig_R['D_red_at_q']}, "
              f"|supp|={len(sig_R['support'])}, max_chip={sig_R['max_chip']}, "
              f"|W|={sig_R['W_size']} ({t_R:.1f}s)", flush=True)

    # Compute z-scores for each metric
    def z(metric_P, metrics_R):
        m = statistics.mean(metrics_R)
        s = statistics.stdev(metrics_R) if len(metrics_R) > 1 else 0
        return (metric_P - m) / s if s > 0 else (0.0 if metric_P == m else float('inf'))

    def stats(metrics):
        return {
            "mean": statistics.mean(metrics),
            "sd": statistics.stdev(metrics) if len(metrics) > 1 else 0,
            "min": min(metrics),
            "max": max(metrics),
            "values": metrics,
        }

    metrics = {}
    for key in ["D_red_at_q", "max_chip", "total_nonq_chips", "W_size"]:
        randoms_vals = [s[key] for s in randoms]
        metrics[key] = {
            "prime": sig_P[key],
            "random": stats(randoms_vals),
            "z": z(sig_P[key], randoms_vals),
        }

    # Compare W (the winnable subtraction set) by intersecting with primes / composites
    cls = classify_vertices(V)
    n_W_primes_P = sum(1 for v in sig_P["winnable_subtraction_set"] if cls[v][0] == "prime")
    n_W_comp_P = sum(1 for v in sig_P["winnable_subtraction_set"] if cls[v][0] == "composite")
    n_W_primes_R = [
        sum(1 for v in s["winnable_subtraction_set"] if cls[v][0] == "prime") for s in randoms
    ]
    n_W_comp_R = [
        sum(1 for v in s["winnable_subtraction_set"] if cls[v][0] == "composite") for s in randoms
    ]
    metrics["W_primes"] = {
        "prime": n_W_primes_P,
        "random": stats(n_W_primes_R),
        "z": z(n_W_primes_P, n_W_primes_R),
    }
    metrics["W_composites"] = {
        "prime": n_W_comp_P,
        "random": stats(n_W_comp_R),
        "z": z(n_W_comp_P, n_W_comp_R),
    }

    # Print summary z-scores
    print(f"  ---- Z-SCORES (prime vs {n_random} random) ----")
    for key in metrics:
        m = metrics[key]
        print(f"  {key:>20s}: D_P={m['prime']:5g}  rand mean={m['random']['mean']:6.2f}  "
              f"sd={m['random']['sd']:5.2f}  z={m['z']}")

    return {
        "N": N,
        "graph_type": graph_type,
        "q": q,
        "n_vertices": n,
        "n_edges": m,
        "genus": g,
        "deg_DP": deg_DP,
        "sig_DP": sig_P,
        "n_random": len(randoms),
        "metrics": metrics,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, nargs="+", default=[16, 32, 64])
    p.add_argument("--graph", choices=["div", "coprime", "both"], default="both")
    p.add_argument("--n_random", type=int, default=20)
    p.add_argument("--q", default="low")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out", default="signatures.json")
    args = p.parse_args()

    results = []
    graphs = ["div", "coprime"] if args.graph == "both" else [args.graph]
    for N in args.N:
        for gt in graphs:
            res = run_signature_experiment(N, gt, args.n_random, args.q, args.seed)
            results.append(res)
            with open(args.out, "w") as f:
                json.dump(results, f, indent=2, default=str)

    print("\n=== SUMMARY (z-scores) ===")
    for r in results:
        if "skipped" in r:
            continue
        print(f"\n  N={r['N']:3d} {r['graph_type']:>7s}: "
              f"deg={r['deg_DP']} g={r['genus']}")
        for key, m in r["metrics"].items():
            print(f"    {key:>20s}: D_P={m['prime']:6g}  z={m['z']}")


if __name__ == "__main__":
    main()
