"""
D31.a — Connected-part-only AHK deviation.

Tests the claim that the +3-5sigma deviation in D31 (S149) reduces
to the Bertrand `(t-1)^{nu(N)}` direct-summand inflation.

If true: when we restrict the matroid to its CONNECTED part (i.e.,
remove the isolated U_{1,1} factors contributed by primes p in (N/2, N]),
the prime-vs-config-model deviation should COLLAPSE to within +/- 2 sigma.

Connected-part matroid:
  - ground = {n in [2, N] : smallest prime factor of n is <= N/2}
            = [2, N] \\ {primes p with N/2 < p <= N}
  - blocks = {primes p <= N/2}
  - r(S) = max bipartite matching from S to blocks (same definition).

Cross-domain ingredient: same as D31 (AHK 2018 matroid Hodge); the
follow-up here is internal -- a different matroid on a subset of the
ground / block set.
"""

import argparse
import json
import math
import os
import random
import sys
import time

# Reuse machinery from main file
sys.path.insert(0, os.path.dirname(__file__))
from d31_ahk_matroid_hodge import (
    primes_up_to, matroid_rank, whitney_char_poly_from_adj,
    log_concavity_slack, configuration_model_adj, run_one,
    reduced_char_poly_coeffs,
)


def adj_for_connected_part(N: int):
    """Build the connected-part matroid: ground excludes isolated big primes."""
    primes = primes_up_to(N)
    big_primes = [p for p in primes if p > N // 2]
    small_primes = [p for p in primes if p <= N // 2]
    if not small_primes:
        return [], 0, [], []
    ground = [n for n in range(2, N + 1) if n not in big_primes]
    block_idx = {q: i for i, q in enumerate(small_primes)}
    adj = [[block_idx[q] for q in small_primes if x % q == 0] for x in ground]
    # Drop ground elements with no neighbours (shouldn't happen by construction)
    keep = [i for i, nb in enumerate(adj) if nb]
    adj = [adj[i] for i in keep]
    ground = [ground[i] for i in keep]
    return adj, len(small_primes), ground, small_primes


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--Ns", type=int, nargs="+", default=[12, 16, 20])
    p.add_argument("--n-controls", type=int, default=30)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out", type=str,
                   default=os.path.join(os.path.dirname(__file__),
                                        "d31a_connected_part_data.json"))
    args = p.parse_args()

    rng = random.Random(args.seed)
    out = {"by_N": {}}

    for N in args.Ns:
        adj_prime, n_blocks, ground, blocks = adj_for_connected_part(N)
        if n_blocks == 0:
            print(f"N={N}: skipping (no small primes)", file=sys.stderr)
            continue
        deg_R = [0] * n_blocks
        for nb in adj_prime:
            for r in nb:
                deg_R[r] += 1
        print(f"\n=== N={N}  conn matroid: |ground|={len(ground)} "
              f"blocks={blocks} deg_R={deg_R} ===", file=sys.stderr)
        # Prime conn matroid
        t0 = time.time()
        prime_res = run_one(N, adj_prime, n_blocks, label=f"prime_conn_{N}")
        t_prime = time.time() - t0
        print(f"  prime conn: rM={prime_res['rM']} |w|={prime_res['abs_coeffs']} "
              f"deltas={prime_res['deltas']} (neg: {prime_res['n_negative_deltas']}) "
              f"wall {t_prime:.2f}s", file=sys.stderr)
        # Config-model controls
        controls = []
        for c in range(args.n_controls):
            adj_c = configuration_model_adj(adj_prime, n_blocks, rng)
            ctrl = run_one(N, adj_c, n_blocks, label=f"conn_ctrl_{c}")
            controls.append(ctrl)
        print(f"  controls: rMs={[c['rM'] for c in controls]}", file=sys.stderr)

        out["by_N"][N] = {"prime_conn": prime_res, "controls": controls,
                          "ground": ground, "blocks": blocks, "deg_R": deg_R}

    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nWritten {args.out}", file=sys.stderr)

    # Stats summary
    print("\n========== SUMMARY (D31.a connected part) ==========", file=sys.stderr)
    for N in args.Ns:
        if N not in out["by_N"]:
            continue
        b = out["by_N"][N]
        pres = b["prime_conn"]
        ctrls = b["controls"]
        rM = pres["rM"]
        ctrl_match = [c["abs_coeffs"] for c in ctrls if c["rM"] == rM]
        print(f"\nN={N} conn  rM(prime)={rM}  rM(ctrl) min={min(c['rM'] for c in ctrls)} "
              f"max={max(c['rM'] for c in ctrls)}  match_rank={len(ctrl_match)}/{len(ctrls)}",
              file=sys.stderr)
        if len(ctrl_match) >= 3:
            for k in range(rM + 1):
                vals = [c[k] for c in ctrl_match]
                mean = sum(vals) / len(vals)
                var = sum((v - mean) ** 2 for v in vals) / max(1, len(vals) - 1)
                std = math.sqrt(var)
                z = (pres["abs_coeffs"][k] - mean) / std if std > 0 else float('nan')
                print(f"  k={k}  prime={pres['abs_coeffs'][k]:>8}  "
                      f"ctrl_mean={mean:>10.2f}  ctrl_std={std:>8.2f}  z={z:>+6.2f}",
                      file=sys.stderr)
            ctrl_slacks = [log_concavity_slack(c) for c in ctrl_match]
            for k in range(len(pres["deltas"])):
                vals = [s[k] for s in ctrl_slacks if k < len(s)]
                if not vals:
                    continue
                mean = sum(vals) / len(vals)
                var = sum((v - mean) ** 2 for v in vals) / max(1, len(vals) - 1)
                std = math.sqrt(var)
                z = (pres["deltas"][k] - mean) / std if std > 0 else float('nan')
                print(f"  delta_{k+1}  prime={pres['deltas'][k]:>10}  "
                      f"ctrl_mean={mean:>12.2f}  ctrl_std={std:>10.2f}  z={z:>+6.2f}",
                      file=sys.stderr)


if __name__ == "__main__":
    main()
