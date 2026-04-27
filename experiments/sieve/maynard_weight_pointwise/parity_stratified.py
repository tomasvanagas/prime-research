"""
Stratify by parity to isolate genuine Maynard sieve content from
trivial 'n odd' detection.

The base experiment shows AUC > 0.88 at theta=0.10 (R<3). At R<3 the
only divisors admitted are 1 and 2, so the weight reduces to a parity
detector: w(n) = 1 if n is odd, w(n) < 1 if n is even (and 2 | n+h_i
for h_i in H={0,2,6}, all of which are even, so all three are even).

To isolate Maynard content from parity content, restrict to odd n and
recompute AUC for prime-in-window. If AUC stays at chance there, the
'high AUC at low theta' is purely parity detection, not sieve content.
"""
from __future__ import annotations
import argparse
import json
import math
import sys
from pathlib import Path

# Reuse the main module
sys.path.insert(0, str(Path(__file__).parent))
from maynard_weight_pointwise import (
    sieve_of_eratosthenes, smallest_prime_factor,
    MaynardConfig, compute_maynard_weight,
    classify_window, mann_whitney_auc, roc_curve,
)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, default=100_000)
    p.add_argument("--window", type=int, default=20_000)
    p.add_argument("--theta", type=float, default=0.20)
    p.add_argument("--Hs", type=str, default="0,2,6")
    p.add_argument("--F", type=str, default="selberg_gpy")
    p.add_argument("--out", type=str, default="parity_results.json")
    args = p.parse_args()

    H = [int(x) for x in args.Hs.split(",")]
    N = args.N
    window = args.window
    theta = args.theta
    R = max(2.0, N ** theta)
    log_R = math.log(R) if R > 1 else 1.0
    sieve_max = N + window + max(H) + 100

    print(f"[setup] N={N}, window={window}, theta={theta}, R={R:.2f}, H={H}")
    is_prime = sieve_of_eratosthenes(sieve_max)
    spf = smallest_prime_factor(sieve_max)

    cfg = MaynardConfig(k=len(H), H=H, R=R, F_choice=args.F)

    # All n
    weights_all = []
    classes_all = []
    for n in range(N, N + window):
        w = compute_maynard_weight(n, cfg, spf, log_R)
        c = classify_window(n, H, is_prime)
        weights_all.append(w)
        classes_all.append((n, c))

    # Stratify by parity (n is odd vs even)
    weights_odd = [w for w, (n, _) in zip(weights_all, classes_all) if n % 2 == 1]
    classes_odd = [c for w, (n, c) in zip(weights_all, classes_all) if n % 2 == 1]
    weights_even = [w for w, (n, _) in zip(weights_all, classes_all) if n % 2 == 0]
    classes_even = [c for w, (n, c) in zip(weights_all, classes_all) if n % 2 == 0]

    def stratify_auc(weights, classes, label):
        pos = [w for w, c in zip(weights, classes) if c["any_prime"]]
        neg = [w for w, c in zip(weights, classes) if not c["any_prime"]]
        pos_strict = [w for w, c in zip(weights, classes) if c["primes_in_window"] >= 2]
        neg_strict = [w for w, c in zip(weights, classes) if c["primes_in_window"] <= 1]
        auc = mann_whitney_auc(pos, neg) if pos and neg else float("nan")
        auc_str = mann_whitney_auc(pos_strict, neg_strict) if pos_strict and neg_strict else float("nan")
        if not pos or not neg:
            roc = {"best_f1": float("nan"), "best_f1_tau": float("nan")}
        else:
            roc = roc_curve(pos, neg, n_points=200)
        return {
            "label": label,
            "n_total": len(weights),
            "n_pos_any": len(pos),
            "n_neg_any": len(neg),
            "n_pos_strict": len(pos_strict),
            "n_neg_strict": len(neg_strict),
            "auc_any": auc,
            "auc_strict": auc_str,
            "best_f1": roc["best_f1"],
            "best_f1_tau": roc["best_f1_tau"],
            "mean_w_pos": sum(pos) / len(pos) if pos else float("nan"),
            "mean_w_neg": sum(neg) / len(neg) if neg else float("nan"),
        }

    out = {
        "config": {"N": N, "window": window, "theta": theta, "R": R, "H": H, "F": args.F},
        "all": stratify_auc(weights_all, [c for _, c in classes_all], "all"),
        "odd_n": stratify_auc(weights_odd, classes_odd, "odd_n"),
        "even_n": stratify_auc(weights_even, classes_even, "even_n"),
    }
    Path(args.out).write_text(json.dumps(out, indent=2, default=str))

    print(f"\n=== Parity stratification at theta={theta:.2f} (R={R:.2f}) ===")
    for k in ["all", "odd_n", "even_n"]:
        d = out[k]
        print(f"  {k:>7}: n={d['n_total']:>6}, AUC_any={d['auc_any']:.4f}, "
              f"AUC_strict={d['auc_strict']:.4f}, F1={d['best_f1']:.3f}, "
              f"mean_w(prime)/mean_w(no-prime)="
              f"{d['mean_w_pos']:.3g}/{d['mean_w_neg']:.3g}")


if __name__ == "__main__":
    main()
