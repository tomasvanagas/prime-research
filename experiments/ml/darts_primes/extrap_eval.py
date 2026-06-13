"""
Re-instantiate the best PRIMES discrete circuit at N' > N and evaluate
generalisation. The discrete description names gates and prev-layer
indices, so it lifts to any input dimension N' >= max(selected) + 1
once we re-anchor inputs.

For input n in [0, 2^N'), we form the N'-bit representation, then take
the first N bits as the input to the trained N-bit circuit. Selected
indices in the bottom layer index those N bits. This tests whether the
trained circuit is a *projection* onto the low N bits or whether it
captures any structure that extends to higher bits.

We also evaluate on the *same* trained N-bit circuit applied directly
to n in [2^N, 2^N + 1000), where bits cycle (n mod 2^N).
"""
import json
import sys
from pathlib import Path

import numpy as np
from sympy import isprime

sys.path.insert(0, str(Path(__file__).parent))
from darts_primes import evaluate_discrete  # type: ignore


def chi_P_arr(n_arr):
    return np.array([1.0 if (n >= 2 and isprime(int(n))) else 0.0 for n in n_arr], dtype=np.float32)


def n_to_bits(n_arr, N):
    return ((n_arr.reshape(-1, 1) >> np.arange(N).reshape(1, N)) & 1).astype(np.float32)


def main(results_path: str):
    with open(results_path) as f:
        R = json.load(f)
    N = R["config"]["N"]
    primes = R["primes"]
    best = min(primes, key=lambda r: r["final_loss"])
    desc = best["discrete_desc"]

    # 1. Cyclic-bit extrapolation: n in [2^N, 2^(N+1)) using LOW N bits.
    print(f"=== Extrapolation: best PRIMES seed (final_loss={best['final_loss']:.4f}) ===")
    for label, (lo, hi) in [
        ("train", (0, 2 ** N)),
        ("[2^N, 2^N+1000)",  (2 ** N, 2 ** N + 1000)),
        ("[2^(N+1), 2^(N+1)+1000)", (2 ** (N + 1), 2 ** (N + 1) + 1000)),
        ("[10000, 11000)", (10000, 11000)),
        ("[10^5, 10^5+1000)", (100000, 101000)),
    ]:
        n_arr = np.arange(lo, hi, dtype=np.int64)
        y_arr = chi_P_arr(n_arr)
        # Cyclic: take low-N bits
        X = n_to_bits(n_arr & ((1 << N) - 1), N)
        res = evaluate_discrete(desc, X, y_arr)
        density = float(y_arr.mean())
        baseline = max(density, 1 - density)
        gap = res["accuracy"] - baseline
        print(f"  {label:30s}: discrete_acc={res['accuracy']:.4f}  density={density:.4f}  "
              f"majority_baseline={baseline:.4f}  gap={gap:+.4f}")


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "/apps/aplikacijos/prime-research/experiments/ml/darts_primes/run/results.json"
    main(path)
