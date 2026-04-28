"""
D31 — Adiprasito-Huh-Katz combinatorial Hodge theory of arithmetic prime-matroid.

Construction: arithmetic transversal matroid M_P^N
  - Ground set E = [2, N]
  - For each prime p <= N, block B_p = {n in E : p | n}
  - r(S) = max bipartite matching from S to {B_p : p in P_<=N}

Whitney expansion (sum over all subsets):
  chi_M(t) = sum_{S <= E} (-1)^|S| t^{r(M) - r(S)}

AHK 2018 (Annals 188, 381) prove that the absolute coefficients
|w_k| of chi_M(t) form a log-concave sequence:
  |w_k|^2 >= |w_{k-1}| |w_{k+1}|

Discriminating signal: slack delta_k := |w_k|^2 - |w_{k-1}| |w_{k+1}|.
Compare prime delta_k to Bernoulli-matched-density control.

Cross-domain ingredient: Adiprasito-Huh-Katz 2018 "Hodge theory for
combinatorial geometries" Annals 188, 381 = arXiv:1511.02888.

Falsification (E):  prime delta_k matches Bernoulli within +/- 2 sigma
                    (modulo W-trick mod-q residue correction)
Falsification (I):  prime delta_k deviates >= 3 sigma; new HL
                    singular-series identity
Falsification (A):  closed-form chi_M(t) in HL singular-series products
                    polylog-evaluable
"""

from __future__ import annotations

import argparse
import json
import math
import os
import random
import sys
import time
from typing import Dict, List, Sequence, Tuple

# -------------------- prime / divisor utilities --------------------

def primes_up_to(N: int) -> List[int]:
    if N < 2:
        return []
    sieve = [True] * (N + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(N ** 0.5) + 1):
        if sieve[i]:
            for j in range(i * i, N + 1, i):
                sieve[j] = False
    return [i for i in range(2, N + 1) if sieve[i]]


def prime_divisors(n: int, primes: Sequence[int]) -> List[int]:
    """Return prime divisors of n that are in `primes`."""
    return [p for p in primes if n % p == 0]


# -------------------- bipartite matching (augmenting-path) --------------------

def _augment(u: int, adj: Sequence[Sequence[int]],
             match_R: List[int], visited: List[bool]) -> bool:
    """Try to find augmenting path from left vertex u."""
    for v in adj[u]:
        if not visited[v]:
            visited[v] = True
            if match_R[v] == -1 or _augment(match_R[v], adj, match_R, visited):
                match_R[v] = u
                return True
    return False


def matroid_rank(mask: int, adj_all: Sequence[Sequence[int]],
                 n_primes: int) -> int:
    """Max matching size for the subset of ground encoded by `mask`."""
    match_R = [-1] * n_primes
    rank = 0
    bits = mask
    while bits:
        low = bits & -bits
        u = low.bit_length() - 1
        visited = [False] * n_primes
        if _augment(u, adj_all, match_R, visited):
            rank += 1
        bits ^= low
    return rank


# -------------------- Whitney expansion of characteristic polynomial --------------------

def whitney_char_poly_from_adj(adj: Sequence[Sequence[int]],
                               n_blocks: int,
                               verbose: bool = False) -> Tuple[List[int], int]:
    """
    Compute characteristic polynomial of the transversal matroid given
    a bipartite-graph adjacency adj[i] = list of block indices for ground
    element i.

    Returns (coeffs, rM) where chi_M(t) = sum_k coeffs[k] t^k and rM = r(M).
    """
    n = len(adj)
    if n == 0:
        return [1], 0
    rM = matroid_rank((1 << n) - 1, adj, n_blocks)

    coeffs = [0] * (rM + 1)
    total = 1 << n
    report_every = max(1, total // 20)

    t0 = time.time()
    for mask in range(total):
        rS = matroid_rank(mask, adj, n_blocks)
        sign = -1 if (bin(mask).count("1") & 1) else 1
        coeffs[rM - rS] += sign
        if verbose and (mask % report_every == 0) and mask > 0:
            elapsed = time.time() - t0
            frac = mask / total
            eta = elapsed * (1 - frac) / frac if frac > 0 else 0
            print(f"  ... {100 * frac:5.1f}%  elapsed {elapsed:6.1f}s  ETA {eta:6.1f}s",
                  file=sys.stderr, flush=True)

    return coeffs, rM


def reduced_char_poly_coeffs(coeffs: Sequence[int]) -> List[int]:
    """
    chi(t) is divisible by (t - 1) for any matroid of rank >= 1 (chi(1)=0).
    Return coefficients of chi(t) / (t - 1), with the same convention
    (index k = coefficient of t^k).
    """
    # chi(t) = sum coeffs[k] t^k.  Divide by (t-1) using synthetic division.
    # If chi(t) = (t - 1) q(t), then q(t) = sum b_k t^k where
    # b_{r-1} = coeffs[r], b_{k} = -coeffs[k+1] - b_{k+1}? Let's derive.
    # chi(t) / (t - 1):  using long division from highest degree.
    n = len(coeffs)
    q = [0] * (n - 1)
    # chi(t) = (t-1) q(t) ==> coeffs[k] = q[k-1] - q[k] (with q[n-1] = 0, q[-1] = 0)
    # So coeffs[n-1] = q[n-2], coeffs[k] = q[k-1] - q[k] for 1 <= k <= n-2,
    # coeffs[0] = -q[0].
    q[n - 2] = coeffs[n - 1]
    for k in range(n - 2, 0, -1):
        q[k - 1] = coeffs[k] + q[k]
    # Sanity: -q[0] should equal coeffs[0]
    assert -q[0] == coeffs[0], f"chi(1) != 0: {coeffs}, q[0]={q[0]}"
    return q


# -------------------- AHK log-concavity slack --------------------

def log_concavity_slack(abs_coeffs: Sequence[int]) -> List[float]:
    """
    Returns delta_k = |w_k|^2 - |w_{k-1}| |w_{k+1}| for k in 1..len-2.
    A non-negative sequence is log-concave (AHK).
    """
    deltas = []
    for k in range(1, len(abs_coeffs) - 1):
        d = abs_coeffs[k] ** 2 - abs_coeffs[k - 1] * abs_coeffs[k + 1]
        deltas.append(d)
    return deltas


# -------------------- experiment runner --------------------

def adj_from_blocks(N: int, block_indices: Sequence[int],
                    W: int = 1) -> Tuple[List[List[int]], int, List[int], List[int]]:
    """
    Build adjacency: ground = {n in [2,N] : gcd(n, W) = 1};
    effective blocks = {q in block_indices : gcd(q, W) = 1}.
    Returns (adj, n_blocks, ground, kept_blocks).
    """
    ground = [n for n in range(2, N + 1) if math.gcd(n, W) == 1]
    kept_blocks = sorted(q for q in block_indices if math.gcd(q, W) == 1)
    block_idx = {q: i for i, q in enumerate(kept_blocks)}
    adj = [[block_idx[q] for q in kept_blocks if x % q == 0] for x in ground]
    # Drop ground elements with no neighbours (would create loops in matroid)
    keep = [i for i, nbrs in enumerate(adj) if len(nbrs) >= 1]
    adj = [adj[i] for i in keep]
    ground = [ground[i] for i in keep]
    return adj, len(kept_blocks), ground, kept_blocks


def configuration_model_adj(adj_ref: Sequence[Sequence[int]], n_blocks: int,
                            rng: random.Random) -> List[List[int]]:
    """
    Configuration-model bipartite graph with same left and right degree
    sequences as adj_ref.  Multi-edges are dedup'd; isolated left nodes
    cannot occur because every left node has at least one stub.
    """
    n = len(adj_ref)
    deg_L = [len(nbrs) for nbrs in adj_ref]
    deg_R = [0] * n_blocks
    for nbrs in adj_ref:
        for r in nbrs:
            deg_R[r] += 1
    stubs_L = [i for i, d in enumerate(deg_L) for _ in range(d)]
    stubs_R = [j for j, d in enumerate(deg_R) for _ in range(d)]
    rng.shuffle(stubs_R)
    edges = set()
    for l, r in zip(stubs_L, stubs_R):
        edges.add((l, r))
    new_adj = [[] for _ in range(n)]
    for l, r in edges:
        new_adj[l].append(r)
    # Sanity: every left node has at least one neighbour (no loops)
    assert all(len(nbrs) >= 1 for nbrs in new_adj), \
        "config model produced an isolated left node (should not happen)"
    return new_adj


def run_one(N: int, adj: Sequence[Sequence[int]], n_blocks: int,
            label: str = "", verbose: bool = False) -> Dict:
    """Run characteristic-polynomial computation for one matroid and return summary."""
    coeffs, rM = whitney_char_poly_from_adj(adj, n_blocks, verbose=verbose)
    abs_coeffs = [abs(c) for c in coeffs]
    deltas = log_concavity_slack(abs_coeffs)
    n_neg = sum(1 for d in deltas if d < 0)
    # Reduced polynomial (chi(t) / (t-1)): only if r(M) >= 1 and ground non-empty
    try:
        red = reduced_char_poly_coeffs(coeffs)
        red_abs = [abs(c) for c in red]
        red_deltas = log_concavity_slack(red_abs)
        red_neg = sum(1 for d in red_deltas if d < 0)
    except (AssertionError, IndexError):
        red, red_abs, red_deltas, red_neg = [], [], [], -1
    return {
        "N": N,
        "n_blocks": n_blocks,
        "label": label,
        "ground_size": N - 1,
        "rM": rM,
        "coeffs": coeffs,
        "abs_coeffs": abs_coeffs,
        "deltas": deltas,
        "n_negative_deltas": n_neg,
        "reduced_coeffs": red,
        "reduced_abs_coeffs": red_abs,
        "reduced_deltas": red_deltas,
        "reduced_n_negative": red_neg,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--Ns", type=int, nargs="+", default=[8, 12, 16, 20])
    parser.add_argument("--n-controls", type=int, default=20,
                        help="number of degree-preserving config-model controls per N")
    parser.add_argument("--Ws", type=int, nargs="+", default=[1],
                        help="W-trick moduli to apply (1=no sieve, 2=odd, 6, 30, 210)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out", type=str,
                        default=os.path.join(os.path.dirname(__file__),
                                             "d31_ahk_matroid_hodge_data.json"))
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    rng = random.Random(args.seed)

    all_results = {"by_NW": {}}

    for N in args.Ns:
        for W in args.Ws:
            key = f"N={N}_W={W}"
            primes = primes_up_to(N)
            print(f"\n=== {key}  (|[2,N]|={N-1}, pi(N)={len(primes)}) ===",
                  file=sys.stderr, flush=True)

            # Prime matroid (W-tricked)
            adj_prime, n_blocks_prime, ground, kept_primes = \
                adj_from_blocks(N, primes, W=W)
            if len(adj_prime) == 0 or n_blocks_prime == 0:
                print(f"  W-trick produced empty matroid; skipping", file=sys.stderr)
                continue
            deg_L = [len(nbrs) for nbrs in adj_prime]
            deg_R = [0] * n_blocks_prime
            for nbrs in adj_prime:
                for r in nbrs:
                    deg_R[r] += 1
            print(f"  W={W}  ground={len(ground)}  blocks={kept_primes}  "
                  f"deg_R={deg_R}",
                  file=sys.stderr, flush=True)

            t0 = time.time()
            prime_res = run_one(N, adj_prime, n_blocks_prime,
                                label=f"prime_W{W}", verbose=args.verbose)
            t_prime = time.time() - t0
            prime_res["wall_time_s"] = t_prime
            prime_res["W"] = W
            prime_res["ground"] = ground
            prime_res["kept_primes"] = kept_primes
            prime_res["deg_L"] = deg_L
            prime_res["deg_R"] = deg_R
            print(f"  prime:  rM={prime_res['rM']}  |w|={prime_res['abs_coeffs']}  "
                  f"wall {t_prime:.2f} s", file=sys.stderr, flush=True)
            print(f"          deltas       ={prime_res['deltas']}  "
                  f"(neg: {prime_res['n_negative_deltas']})",
                  file=sys.stderr, flush=True)

            # Configuration-model controls
            controls: List[Dict] = []
            for c in range(args.n_controls):
                adj_ctrl = configuration_model_adj(adj_prime, n_blocks_prime, rng)
                t0 = time.time()
                ctrl_res = run_one(N, adj_ctrl, n_blocks_prime,
                                   label=f"config_model_{c}_W{W}", verbose=False)
                t_ctrl = time.time() - t0
                ctrl_res["wall_time_s"] = t_ctrl
                ctrl_res["W"] = W
                controls.append(ctrl_res)
            print(f"  controls: {len(controls)}  rMs={[c['rM'] for c in controls]}",
                  file=sys.stderr, flush=True)
            all_results["by_NW"][key] = {
                "prime": prime_res,
                "config_model": controls,
                "N": N, "W": W,
            }

    # Persist results
    with open(args.out, "w") as fh:
        json.dump(all_results, fh, indent=2)
    print(f"\nWritten {args.out}", file=sys.stderr, flush=True)

    # Summary statistics
    print("\n========== SUMMARY ==========", file=sys.stderr, flush=True)
    for key, bundle in all_results["by_NW"].items():
        pres = bundle["prime"]
        ctrls = bundle["config_model"]
        rM = pres["rM"]
        ctrl_rMs = [c["rM"] for c in ctrls]
        ctrl_coefs = [c["abs_coeffs"] for c in ctrls if c["rM"] == rM]
        n_match_rank = sum(1 for r in ctrl_rMs if r == rM)
        print(f"\n{key}  rM(prime)={rM}  rM(ctrl)={min(ctrl_rMs)}..{max(ctrl_rMs)}  "
              f"match_rank={n_match_rank}/{len(ctrl_rMs)}",
              file=sys.stderr, flush=True)
        if len(ctrl_coefs) >= 3:
            for k in range(rM + 1):
                vals = [c[k] for c in ctrl_coefs]
                mean = sum(vals) / len(vals)
                var = sum((v - mean) ** 2 for v in vals) / max(1, len(vals) - 1)
                std = math.sqrt(var)
                z = (pres["abs_coeffs"][k] - mean) / std if std > 0 else float('nan')
                print(f"  k={k}  prime={pres['abs_coeffs'][k]:>10}  "
                      f"ctrl_mean={mean:>12.2f}  ctrl_std={std:>10.2f}  z={z:>+7.2f}",
                      file=sys.stderr, flush=True)
            # Slack delta_k z-scores
            ctrl_slacks = [log_concavity_slack(c) for c in ctrl_coefs]
            for k in range(len(pres["deltas"])):
                vals = [s[k] for s in ctrl_slacks if k < len(s)]
                if not vals:
                    continue
                mean = sum(vals) / len(vals)
                var = sum((v - mean) ** 2 for v in vals) / max(1, len(vals) - 1)
                std = math.sqrt(var)
                z = (pres["deltas"][k] - mean) / std if std > 0 else float('nan')
                print(f"  delta_{k+1}  prime={pres['deltas'][k]:>12}  "
                      f"ctrl_mean={mean:>14.2f}  ctrl_std={std:>10.2f}  z={z:>+7.2f}",
                      file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
