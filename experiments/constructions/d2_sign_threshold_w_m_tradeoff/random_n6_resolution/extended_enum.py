"""
C8.b — Random-control F4 resolution at N=6 via extended column enumeration.

Composition: E5.3 (PRIMES TC^0 open) + S84 column-enum framework
(`enum_d2_smart.py`) + C7-S89 oddness predictor + C8/S127 partial result.

S84 column-enum was W=1-only (K=1458 distinct sign-threshold truth tables
on N=6). This script extends to W >= 2 and runs `depth2_search` on:

  (a) PRIMES (sanity check: should reproduce S127 curve M*=4, 3, 3 at
      W=2, 3, 4),
  (b) matched-density random (seed parameter).

The S127 random ILP via direct (joint) encoding `sat_depth2_ilp.py`
returned UNKNOWN at (W=4, M=3) even at 600 s. Column enumeration
hands the bottom-layer threshold truth tables to the ILP as a fixed
catalog, which removes bilinear constraints (alpha * beta) and
typically resolves cleanly within seconds.

CLI:
  python extended_enum.py --N 6 --W 2 --target both --M-list 2,3,4 \\
                          --seed 42 --time-limit 120 --out out.json
"""
from __future__ import annotations
import argparse, json, os, sys, time
from itertools import product, combinations
from typing import Dict, List, Optional, Set, Tuple

# Reuse S84's depth2_search and primes/random tables.
HERE = os.path.dirname(os.path.abspath(__file__))
S84 = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(HERE))),
    "circuit_complexity", "sat_tc0_primes_n8",
)
sys.path.insert(0, S84)
from enum_d2_smart import (  # type: ignore
    primes_table, random_table, depth2_search,
)


def bits(x: int, N: int) -> List[int]:
    return [(x >> i) & 1 for i in range(N)]


def enumerate_thresholds(N: int, W: int, k_max: Optional[int] = None,
                         drop_complements: bool = True) -> List[Dict]:
    """
    Enumerate distinct sign-threshold truth tables on N inputs with weight
    bound W (each w_i in {-W, ..., W}) and at most k_max nonzero positions.

    Returns: list of {'w': tuple, 'T': int, 'tt': tuple-of-N=2^N bits}.
    Truth tables are deduplicated up to complementation (alpha = +-1
    covers complements at the top layer).
    """
    if k_max is None:
        k_max = N
    domain_bits = [bits(x, N) for x in range(2 ** N)]
    seen: Set[Tuple[int, ...]] = set()
    out: List[Dict] = []
    for k in range(1, k_max + 1):
        for support in combinations(range(N), k):
            # weights in {-W, ..., W}, nonzero on support
            for w_pat in product(range(-W, W + 1), repeat=k):
                # Skip patterns where any active position has w=0
                if any(v == 0 for v in w_pat):
                    continue
                w = [0] * N
                for ii, idx in enumerate(support):
                    w[idx] = w_pat[ii]
                w = tuple(w)
                sums = [sum(w[i] * b[i] for i in range(N))
                        for b in domain_bits]
                s_min, s_max = min(sums), max(sums)
                # T sweeps the meaningful thresholds; T = s_min gives
                # the all-1 TT, T = s_max + 1 gives all-0 TT.
                for T in range(s_min, s_max + 2):
                    tt = tuple(1 if s >= T else 0 for s in sums)
                    if all(v == 0 for v in tt) or all(v == 1 for v in tt):
                        continue
                    if drop_complements:
                        ctt = tuple(1 - v for v in tt)
                        canon = min(tt, ctt)
                        if canon in seen:
                            continue
                        seen.add(canon)
                    else:
                        if tt in seen:
                            continue
                        seen.add(tt)
                    out.append({"w": w, "T": T, "tt": tt})
    return out


def cleanup(o):
    if isinstance(o, dict):
        return {str(k): cleanup(v) for k, v in o.items()}
    if isinstance(o, list):
        return [cleanup(v) for v in o]
    if isinstance(o, tuple):
        return list(o)
    try:
        json.dumps(o); return o
    except (TypeError, ValueError):
        return str(o)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, default=6)
    p.add_argument("--W", type=int, default=2)
    p.add_argument("--k-max", type=int, default=None)
    p.add_argument("--M-list", type=str, default="2,3,4")
    p.add_argument("--target", choices=["primes", "random", "both"],
                   default="both")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--time-limit", type=int, default=120)
    p.add_argument("--out", type=str, required=True)
    args = p.parse_args()

    M_list = [int(x) for x in args.M_list.split(",")]
    print(f"N={args.N}, W={args.W}, k_max={args.k_max}, M_list={M_list}",
          flush=True)
    t0 = time.time()
    cands = enumerate_thresholds(args.N, args.W, args.k_max)
    print(f"Enumerated K={len(cands)} candidates in {time.time()-t0:.1f}s",
          flush=True)

    out: Dict = {"N": args.N, "W": args.W, "k_max": args.k_max,
                 "K": len(cands), "M_list": M_list,
                 "time_limit": args.time_limit, "seed": args.seed}

    primes = primes_table(args.N)
    weight = sum(primes)
    print(f"PRIMES weight: {weight}/{2**args.N}", flush=True)

    if args.target in ("primes", "both"):
        print("== PRIMES ==", flush=True)
        out["primes"] = depth2_search(primes, args.N, cands, M_list,
                                      args.time_limit)

    if args.target in ("random", "both"):
        rt = random_table(args.N, weight, seed=args.seed)
        print(f"== RANDOM (seed={args.seed}) ==", flush=True)
        out["random"] = depth2_search(rt, args.N, cands, M_list,
                                      args.time_limit)

    with open(args.out, "w") as f:
        json.dump(cleanup(out), f, indent=2)
    print(f"\nWrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
