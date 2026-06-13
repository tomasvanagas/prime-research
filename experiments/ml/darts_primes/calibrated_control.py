"""
Calibrated-1-bit control: a random Boolean function on [0, 2^N) with
matched parity-conditional distribution to PRIMES.

PRIMES at N=12: among 2048 odd n's, 563 are prime. Among 2048 even n's,
1 is prime (n=2). Calibrated control places the same 563 / 1 prime
indicators uniformly at random within the odd / even halves.

If DARTS achieves the SAME loss on PRIMES and on the calibrated-1-bit
control, then the entire PRIMES advantage over uncalibrated random is
explained by oddness — supporting S84's mechanism interpretation at
larger N. If DARTS shows a residual gap, it indicates structure beyond
oddness.
"""
import json
import time
from pathlib import Path

import numpy as np
import torch
from sympy import isprime

from darts_primes import DartsCircuit, n_to_bits, train_one


def primes_distribution(N: int):
    """Return arrays (n, chi_P) and the parity-conditional counts."""
    n_arr = np.arange(2 ** N, dtype=np.int64)
    y = np.array([1.0 if isprime(int(k)) else 0.0 for k in n_arr], dtype=np.float32)
    odd_count = int(y[n_arr % 2 == 1].sum())
    even_count = int(y[n_arr % 2 == 0].sum())
    return n_arr, y, odd_count, even_count


def calibrated_target(seed: int, N: int):
    """Random Boolean on [0, 2^N) with matched parity-conditional density to PRIMES."""
    n_arr, y, odd_count, even_count = primes_distribution(N)
    rng = np.random.RandomState(seed)
    table = np.zeros(2 ** N, dtype=np.float32)
    odd_idx = np.arange(2 ** N)[n_arr % 2 == 1]
    even_idx = np.arange(2 ** N)[n_arr % 2 == 0]
    odd_picks = rng.choice(odd_idx, size=odd_count, replace=False)
    even_picks = rng.choice(even_idx, size=even_count, replace=False)
    table[odd_picks] = 1.0
    table[even_picks] = 1.0

    def fn(n_array: np.ndarray) -> np.ndarray:
        out = np.zeros_like(n_array, dtype=np.float32)
        mask = (n_array >= 0) & (n_array < 2 ** N)
        out[mask] = table[n_array[mask]].astype(np.float32)
        return out

    return fn


def run(args):
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    results = {"config": vars(args), "calibrated": []}

    for seed in range(args.n_seeds):
        t0 = time.time()
        ctrl = calibrated_target(seed=2000 + seed, N=args.N)
        print(f"CALIBRATED seed {seed}")
        res = train_one(
            target_fn=ctrl,
            N=args.N,
            G1=args.G1,
            G2=args.G2,
            epochs=args.epochs,
            lr=args.lr,
            seed=seed + 200,
            verbose=args.verbose,
        )
        res["seed"] = seed
        res["wallclock_s"] = time.time() - t0
        results["calibrated"].append({
            "seed": seed,
            "final_loss": res["losses"][-1],
            "min_loss": min(res["losses"]),
            "soft_acc": res["soft_acc"],
            "discrete_acc": res["discrete_acc"],
            "discrete_desc": res["discrete_desc"],
            "wallclock_s": res["wallclock_s"],
            "losses": res["losses"],
        })
        print(f"  -> final_loss={res['losses'][-1]:.4f} soft_acc={res['soft_acc']:.4f} discrete_acc={res['discrete_acc']:.4f} ({res['wallclock_s']:.1f}s)")

    with open(out_dir / "calibrated_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"Calibrated results saved to {out_dir}/calibrated_results.json")


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, default=12)
    p.add_argument("--G1", type=int, default=16)
    p.add_argument("--G2", type=int, default=16)
    p.add_argument("--epochs", type=int, default=300)
    p.add_argument("--lr", type=float, default=1e-2)
    p.add_argument("--n_seeds", type=int, default=5)
    p.add_argument("--out-dir", type=str,
                   default="/apps/aplikacijos/prime-research/experiments/ml/darts_primes/run")
    p.add_argument("--verbose", action="store_true")
    args = p.parse_args()
    run(args)
