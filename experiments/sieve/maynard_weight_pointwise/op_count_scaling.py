"""
Operation count scaling: does Maynard weight evaluation fit in polylog(N)?

For TC^0 feasibility we need polylog(N) operations per single-n
evaluation. The divisor-tuple enumeration cost is bounded below by the
number of squarefree divisor triples (d_1, d_2, d_3) with d_i | n+h_i,
d_i ≤ R, gcd(d_i, d_j)=1, d_1 d_2 d_3 ≤ R.

Theoretical scaling:
  - per-coord squarefree divisors of n+h_i up to R: bounded by 2^omega(n+h_i),
    average d(n+h_i) ~ log n (Hardy-Ramanujan)
  - in product, mean tuple count ~ (log n)^k truncated by simplex
  - so mean ops per n = polylog(n) IS plausible
  - BUT max ops per n = max over n of d(n+h_i)^k can be n^o(1) (smooth n)
  - and the per-divisor cost is itself a search through ALL divisors
    of m up to R, which costs O(R) by trial division (or sqrt(m) by
    factorization, not polylog without Pratt certificate)

Question: in PRACTICE, does the median / 95th percentile op count
stay polylog as N grows?
"""
from __future__ import annotations
import argparse
import json
import math
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from maynard_weight_pointwise import (
    smallest_prime_factor, MaynardConfig, count_ops_for_evaluation,
)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--theta", type=float, default=0.20)
    p.add_argument("--Hs", type=str, default="0,2,6")
    p.add_argument("--out", type=str, default="op_count_scaling.json")
    args = p.parse_args()

    H = [int(x) for x in args.Hs.split(",")]

    # Test at increasing N. Sample 500 random odd n in [N, N + window].
    Ns = [10_000, 31_623, 100_000, 316_227, 1_000_000]
    sample_size = 500

    results = []
    for N in Ns:
        R = max(2.0, N ** args.theta)
        sieve_max = N + sample_size + max(H) + 100
        print(f"[setup] N={N}, R={R:.2f}, sieve...", flush=True)
        t0 = time.time()
        spf = smallest_prime_factor(sieve_max)
        print(f"  sieve in {time.time()-t0:.1f}s", flush=True)
        cfg = MaynardConfig(k=len(H), H=H, R=R)

        # Sample odd n in [N, N+sample_size*4]
        ns_sampled = []
        nn = N | 1  # ensure odd
        for _ in range(sample_size):
            ns_sampled.append(nn)
            nn += 2

        ops = []
        for nn in ns_sampled:
            r = count_ops_for_evaluation(nn, cfg, spf)
            ops.append(r["n_tuples_coprime"])

        ops_sorted = sorted(ops)
        results.append({
            "N": N,
            "R": R,
            "log_N": math.log(N),
            "log2_N": math.log2(N),
            "mean": sum(ops) / len(ops),
            "median": ops_sorted[len(ops_sorted) // 2],
            "p95": ops_sorted[int(0.95 * len(ops_sorted))],
            "p99": ops_sorted[int(0.99 * len(ops_sorted))],
            "max": max(ops),
        })
        last = results[-1]
        print(f"  N={N}: mean={last['mean']:.1f}, median={last['median']}, "
              f"p95={last['p95']}, p99={last['p99']}, max={last['max']}")

    Path(args.out).write_text(json.dumps(results, indent=2))

    # Fit log-log scaling
    print("\n=== Op count scaling (theta = " + f"{args.theta:.2f}" + ") ===")
    print(f"{'N':>10} {'R':>8} {'mean ops':>10} {'median':>8} "
          f"{'p95':>8} {'p99':>8} {'log_R':>7} {'mean/log_R^k':>14}")
    k = len(H)
    for r in results:
        log_R = math.log(r["R"]) if r["R"] > 1 else 1.0
        print(f"{r['N']:>10} {r['R']:>8.2f} {r['mean']:>10.2f} "
              f"{r['median']:>8} {r['p95']:>8} {r['p99']:>8} "
              f"{log_R:>7.2f} {r['mean']/(log_R**k):>14.4f}")

    # Attempt power-law fit: mean ~ N^alpha
    if len(results) >= 2:
        N0 = results[0]["N"]
        ops0 = results[0]["mean"]
        N1 = results[-1]["N"]
        ops1 = results[-1]["mean"]
        if ops0 > 0 and ops1 > 0:
            alpha = math.log(ops1 / ops0) / math.log(N1 / N0)
            print(f"\nPower-law fit: mean_ops ~ N^{alpha:.4f}")
            print(f"For polylog: alpha = 0; observed alpha = {alpha:.4f}")


if __name__ == "__main__":
    main()
