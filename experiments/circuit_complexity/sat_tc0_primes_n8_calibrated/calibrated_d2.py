"""
C7: calibrated 1-bit-bias random control for the S84 PRIMES-vs-random
depth-2 sign-threshold gap.

S84 found that at N=6 PRIMES needs M=6 depth-2 sign-threshold (W=1)
gates while ALL 10 unbiased matched-density random controls need
M ∈ {7, 8} (binomial null p < 0.001). The proposed mechanism: PRIMES
has a 75% bit_0 predictor (predict 1 if x odd, else 0) at N=6, while
matched-density unbiased random has bit_0 accuracy ~0.66.

This experiment constructs "calibrated random" boolean functions that
match PRIMES's class-conditional distribution exactly and asks
whether the M=6 result is fully explained by the oddness advantage.

Two calibration modes:
  - STRATIFIED: pick exactly k_odd random odd values + k_even random
    even values (k_odd, k_even = 17, 1 to match PRIMES at N=6).
    Total weight always 18; bit_0 predictor accuracy always 0.75.
  - BERNOULLI: independent Bernoulli per x with class-conditional
    probabilities (17/32 if odd, 1/32 if even). Variable weight.

We run depth2_search (S84's enum_d2_smart harness) on N_SAMPLES=20
samples in each mode and compare the min_M distribution to the
S84 baselines: PRIMES (M=6) and unbiased random (M ∈ {7,8}).
"""
from __future__ import annotations
import argparse, json, random, time, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                 "..", "sat_tc0_primes_n8"))
from enum_d2_smart import (enumerate_w1_thresholds_pruned,
                            depth2_search, primes_table)


def calibrated_stratified(N: int, k_odd: int, k_even: int, seed: int):
    """Pick k_odd random odd indices + k_even random even indices, set 1."""
    rng = random.Random(seed)
    odd_idxs  = [x for x in range(2**N) if  x & 1]
    even_idxs = [x for x in range(2**N) if not (x & 1)]
    chosen_odd  = rng.sample(odd_idxs,  k_odd)
    chosen_even = rng.sample(even_idxs, k_even)
    t = [0]*(2**N)
    for k in chosen_odd  + chosen_even:
        t[k] = 1
    return t


def calibrated_bernoulli(N: int, p_given_odd: float, p_given_even: float,
                         seed: int):
    """Independent Bernoulli(p_given_odd) on odd x, B(p_given_even) on even."""
    rng = random.Random(seed)
    t = [0]*(2**N)
    for x in range(2**N):
        p = p_given_odd if (x & 1) else p_given_even
        t[x] = 1 if rng.random() < p else 0
    return t


def bit0_accuracy(t):
    """Accuracy of the rule 'predict 1 iff x is odd'."""
    n = len(t)
    return sum(1 for x in range(n) if (1 if (x & 1) else 0) == t[x]) / n


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=6)
    parser.add_argument("--n-samples", type=int, default=20)
    parser.add_argument("--mode", choices=["stratified", "bernoulli", "both"],
                        default="both")
    parser.add_argument("--Mlist", type=str, default="3,4,5,6,7,8")
    parser.add_argument("--time-limit", type=int, default=120)
    parser.add_argument("--out", type=str,
                        default="calibrated_d2_n6.json")
    parser.add_argument("--seed-base", type=int, default=1000)
    args = parser.parse_args()

    M_list = [int(x) for x in args.Mlist.split(",")]
    N = args.N

    # Reference PRIMES stats.
    p = primes_table(N)
    n_p_odd  = sum(p[x] for x in range(2**N) if  x & 1)
    n_p_even = sum(p[x] for x in range(2**N) if not (x & 1))
    n_odd, n_even = 2**(N-1), 2**(N-1)
    p_odd, p_even = n_p_odd / n_odd, n_p_even / n_even
    print(f"N={N}; PRIMES weight={sum(p)}; "
          f"#prime∩odd={n_p_odd}, #prime∩even={n_p_even}")
    print(f"P(prime|odd) = {p_odd:.4f},  P(prime|even) = {p_even:.4f}")
    print(f"PRIMES bit_0 accuracy = {bit0_accuracy(p):.4f}")

    # Build candidate set once.
    t0 = time.time()
    cands = enumerate_w1_thresholds_pruned(N, k_max=N)
    print(f"K = {len(cands)} candidates "
          f"({time.time()-t0:.1f}s enumeration)")

    out = {"N": N, "M_list": M_list,
           "n_samples": args.n_samples,
           "ref": {"primes_weight": sum(p),
                   "p_odd": p_odd, "p_even": p_even,
                   "n_p_odd": n_p_odd, "n_p_even": n_p_even,
                   "primes_bit0_acc": bit0_accuracy(p)}}

    # === Stratified calibration ===
    if args.mode in ("stratified", "both"):
        out["stratified"] = {}
        print(f"\n=== Stratified (k_odd={n_p_odd}, k_even={n_p_even}) ===")
        for i in range(args.n_samples):
            seed = args.seed_base + i
            t = calibrated_stratified(N, n_p_odd, n_p_even, seed=seed)
            wt, b0 = sum(t), bit0_accuracy(t)
            print(f"\n-- stratified sample seed={seed}, "
                  f"weight={wt}, bit0_acc={b0:.4f} --")
            r = depth2_search(t, N, cands, M_list,
                              time_limit=args.time_limit)
            out["stratified"][seed] = {"weight": wt, "bit0_acc": b0,
                                        "search": r}

    # === Bernoulli calibration ===
    if args.mode in ("bernoulli", "both"):
        out["bernoulli"] = {}
        print(f"\n=== Bernoulli (p_odd={p_odd:.4f}, "
              f"p_even={p_even:.4f}) ===")
        for i in range(args.n_samples):
            seed = args.seed_base + 10000 + i
            t = calibrated_bernoulli(N, p_odd, p_even, seed=seed)
            wt, b0 = sum(t), bit0_accuracy(t)
            print(f"\n-- bernoulli sample seed={seed}, "
                  f"weight={wt}, bit0_acc={b0:.4f} --")
            r = depth2_search(t, N, cands, M_list,
                              time_limit=args.time_limit)
            out["bernoulli"][seed] = {"weight": wt, "bit0_acc": b0,
                                       "search": r}

    # Cleanup tuples & write.
    def cleanup(o):
        if isinstance(o, dict):  return {k: cleanup(v) for k, v in o.items()}
        if isinstance(o, list):  return [cleanup(v) for v in o]
        if isinstance(o, tuple): return list(o)
        try:
            json.dumps(o); return o
        except (TypeError, ValueError):
            return str(o)

    with open(args.out, "w") as f:
        json.dump(cleanup(out), f, indent=2)

    # === Summary ===
    print("\n\n=== Summary ===")
    print(f"PRIMES (S84): min_M = 6")
    print(f"Unbiased random (S84, 10 seeds): min_M ∈ {{7, 8}}, mean 7.6")

    for mode in ("stratified", "bernoulli"):
        if mode not in out: continue
        print(f"\nCalibrated {mode}:")
        ms = []
        for seed, entry in out[mode].items():
            mm = entry["search"].get("min_M", None)
            ms.append(mm)
            wt, b0 = entry["weight"], entry["bit0_acc"]
            print(f"  seed {seed}: weight={wt}, bit0_acc={b0:.4f}, "
                  f"min_M = {mm}")
        valid = [m for m in ms if isinstance(m, int)]
        if valid:
            print(f"  -> N={len(valid)}/{len(ms)} resolved, "
                  f"min_M distribution: mean={sum(valid)/len(valid):.2f}, "
                  f"min={min(valid)}, max={max(valid)}, "
                  f"hist={ {v: valid.count(v) for v in sorted(set(valid))} }")
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
