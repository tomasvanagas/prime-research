"""
C8 — Depth-2 sign-threshold W-vs-M tradeoff for PRIMES.

Composes E5.3 (PRIMES TC^0 open) with the S84 ILP framework
(`experiments/circuit_complexity/sat_tc0_primes_n8/sat_depth2_ilp.py`),
extending the W=1-only column to a full (W, M) grid.

For each (W, M) cell we ask: does there exist a depth-2 sign-threshold
circuit with M bottom gates of input-weight bound W computing PRIMES on
{0, ..., 2^N - 1}? The smallest feasible M for each W gives the
tradeoff M*(W).

The encoder is identical to S84's (Big-M ILP, alpha = 2*beta - 1, AND
linearisation, T-symmetry breaking). Only the W parameter is swept.

CLI (see --help): scans (W, M) grid, writes JSON + log.

Usage:
  python d2_w_m_tradeoff.py --N 4 --Ws 1,2,4 --Ms 1,2,3,4 --time-limit 30
"""
from __future__ import annotations
import argparse, json, sys, time, os
from typing import List, Optional, Tuple, Dict

# Reuse the S84 encoder.
HERE = os.path.dirname(os.path.abspath(__file__))
S84 = os.path.join(
    os.path.dirname(os.path.dirname(HERE)),
    "circuit_complexity", "sat_tc0_primes_n8",
)
sys.path.insert(0, S84)
from sat_depth2_ilp import (  # type: ignore
    primes_table, random_table, encode_depth2, verify,
)
import pulp


def solve_one(target: List[int], N: int, M: int, W: int,
              time_limit: int) -> Dict:
    prob, vars_ = encode_depth2(target, N, M, W=W, symmetry_break=True)
    solver = pulp.PULP_CBC_CMD(msg=0, timeLimit=time_limit, threads=4)
    t0 = time.time()
    status_int = prob.solve(solver)
    elapsed = time.time() - t0
    status_name = pulp.LpStatus[status_int]
    out: Dict = {"elapsed_s": elapsed, "status_name": status_name,
                 "N": N, "M": M, "W": W}

    if status_name == "Optimal":
        try:
            T_top_v = int(round(vars_["T_top"].value()))
            T_v = [int(round(vars_["T"][j].value())) for j in range(M)]
            w_v = [[int(round(vars_["w"][j][i].value())) for i in range(N)]
                   for j in range(M)]
            beta_v = [int(round(vars_["beta"][j].value())) for j in range(M)]
            alpha_v = [2*b - 1 for b in beta_v]
            model = {"T_top": T_top_v, "T": T_v, "w": w_v, "alpha": alpha_v}
            ok, ncorr = verify(model, target, N)
            out["status"] = "sat" if ok else "verify_fail"
            out["model"] = model
            out["correct"] = ncorr
            out["total"] = 2 ** N
        except Exception as e:
            out["status"] = "unknown"
            out["err"] = str(e)
    elif status_name == "Infeasible":
        out["status"] = "unsat"
    elif status_name == "Not Solved":
        out["status"] = "unknown"
    else:
        out["status"] = "unknown"
    return out


def scan_grid(N: int, Ws: List[int], Ms: List[int],
              time_limit: int, target_kind: str = "primes",
              seed: int = 42) -> Dict:
    if target_kind == "primes":
        target = primes_table(N)
    elif target_kind == "random":
        target = random_table(N, sum(primes_table(N)), seed=seed)
    else:
        raise ValueError(target_kind)

    weight = sum(target)
    print(f"[{target_kind}] N={N}  weight={weight}/{2**N}", flush=True)
    print(f"  Ws={Ws}  Ms={Ms}  time_limit={time_limit}s", flush=True)

    cells: List[Dict] = []
    # Per-W early termination: once we hit a SAT M, we stop scanning
    # higher M for that W.
    M_star: Dict[int, Optional[int]] = {W: None for W in Ws}

    for W in Ws:
        for M in Ms:
            print(f"  N={N} W={W:2d} M={M:2d}: ", end="", flush=True)
            r = solve_one(target, N, M, W, time_limit)
            print(f"{r['status']:7s} ({r['elapsed_s']:6.1f}s)"
                  + (f"  corr={r.get('correct')}/{r.get('total')}"
                     if r['status'] == 'sat' else ""), flush=True)
            cells.append(r)
            if r["status"] == "sat" and M_star[W] is None:
                M_star[W] = M
                # Skip larger M for this W — already SAT.
                break
            if r["status"] == "unsat":
                continue
            # "unknown" — record and try larger M anyway.

    return {
        "N": N,
        "Ws": Ws,
        "Ms": Ms,
        "time_limit": time_limit,
        "target_kind": target_kind,
        "weight": weight,
        "cells": cells,
        "M_star": M_star,
    }


def cleanup(o):
    if isinstance(o, dict):
        return {str(k): cleanup(v) for k, v in o.items()}
    if isinstance(o, list):
        return [cleanup(v) for v in o]
    try:
        json.dumps(o); return o
    except (TypeError, ValueError):
        return str(o)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, default=4)
    p.add_argument("--Ws", type=str, default="1,2,4,8")
    p.add_argument("--Ms", type=str, default="1,2,3,4,5,6,7,8")
    p.add_argument("--time-limit", type=int, default=30)
    p.add_argument("--target", choices=["primes", "random", "both"],
                   default="primes")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out", type=str, required=True)
    args = p.parse_args()

    Ws = [int(x) for x in args.Ws.split(",")]
    Ms = [int(x) for x in args.Ms.split(",")]
    out: Dict = {"N": args.N, "Ws": Ws, "Ms": Ms,
                 "time_limit": args.time_limit}
    if args.target in ("primes", "both"):
        out["primes"] = scan_grid(args.N, Ws, Ms, args.time_limit,
                                  "primes", args.seed)
    if args.target in ("random", "both"):
        out["random"] = scan_grid(args.N, Ws, Ms, args.time_limit,
                                  "random", args.seed)

    with open(args.out, "w") as f:
        json.dump(cleanup(out), f, indent=2)
    print(f"\nWrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
