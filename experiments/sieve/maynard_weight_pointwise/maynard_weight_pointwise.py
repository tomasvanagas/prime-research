"""
Maynard multidimensional sieve weight w(n) as a single-n primality witness.

ATTACK_VECTORS.md §A5. Wild swing session.

The Maynard 2015 sieve weight is

    w(n) = (Σ_{d_1|n+h_1, ..., d_k|n+h_k, gcd(d_i,d_j)=1, d_1...d_k ≤ R}
              μ(d_1)...μ(d_k) F(log d_1/log R, ..., log d_k/log R))^2

with admissible tuple H = {h_1,...,h_k} and F a smooth function on the
simplex {x_i ≥ 0, Σ x_i ≤ 1}.

By Maynard's main theorem (Theorem 1.1, arxiv 1311.4600), for k large
and R = N^{1/4 - δ}, there is a constant c(k) > 0 such that

    Σ_{N ≤ n < 2N} w(n) χ_P(n+h_i) ≥ c(k) Σ_{N ≤ n < 2N} w(n)

for some i ∈ [k]. This is a POSITIVITY-AT-AGGREGATE result; it does NOT
state w(n) > τ* ⇒ n+h_i prime.

QUESTION (this experiment, A5 from ATTACK_VECTORS.md):
    Is w(n) > τ* a USEFUL POINTWISE WITNESS for "(n, n+2, n+6) contains
    a prime"?  If yes, Maynard weight is a candidate TC^0 primality test.
    If no, Maynard's most-refined-known sieve does NOT pointwise-separate
    primes — closing A5 as B-grade structural negative.

PROCEDURE:
    1. For k=3, H={0,2,6}, F(x) = (1-x_1-x_2-x_3)^2, R = N^θ
       (θ ∈ {0.1, 0.25, 0.4} tested).
    2. For each n in [N, 2N], compute w(n).
    3. Classify n by whether (n, n+2, n+6) contains ≥1 prime.
    4. Mann-Whitney AUC of w | n in admissible class vs not.
    5. ROC curve, find τ* maximizing precision*recall, F1 score, accuracy.

If AUC ≤ 0.55: pointwise no-information; B-grade structural negative.
If 0.55 < AUC < 0.85: weak pointwise signal; refine.
If AUC ≥ 0.85: A-grade — pointwise primality witness via the sieve.

Failure profile (E mode of A5): AUC ≤ 0.55 + structural reason
(weight reflects divisor tuple count, not prime location).
"""

from __future__ import annotations

import argparse
import json
import math
import random
import time
from collections import defaultdict
from dataclasses import dataclass
from itertools import product
from pathlib import Path
from typing import Iterable

# ---------------------------------------------------------------------------
# Number-theoretic primitives
# ---------------------------------------------------------------------------


def sieve_of_eratosthenes(n: int) -> list[bool]:
    is_prime = [True] * (n + 1)
    is_prime[0] = is_prime[1] = False
    for p in range(2, int(n**0.5) + 1):
        if is_prime[p]:
            for q in range(p * p, n + 1, p):
                is_prime[q] = False
    return is_prime


def smallest_prime_factor(n_max: int) -> list[int]:
    """spf[n] = smallest prime factor of n; spf[0]=spf[1]=0."""
    spf = [0] * (n_max + 1)
    for i in range(2, n_max + 1):
        if spf[i] == 0:  # i is prime
            for j in range(i, n_max + 1, i):
                if spf[j] == 0:
                    spf[j] = i
    return spf


def squarefree_divisors_up_to(n: int, limit: int, spf: list[int]) -> list[int]:
    """All squarefree divisors d of n with d ≤ limit. Includes 1."""
    factors = []
    m = n
    while m > 1 and m < len(spf):
        p = spf[m]
        factors.append(p)
        while m % p == 0:
            m //= p
    if m >= len(spf):
        # fallback (should not occur for n+h_i ≤ n_max)
        f = []
        x = m
        d = 2
        while d * d <= x:
            if x % d == 0:
                f.append(d)
                while x % d == 0:
                    x //= d
            d += 1
        if x > 1:
            f.append(x)
        factors.extend(f)
    # Squarefree divisors built from these primes
    divisors = [1]
    for p in factors:
        new = []
        for d in divisors:
            dp = d * p
            if dp <= limit:
                new.append(dp)
        divisors.extend(new)
    return [d for d in divisors if d <= limit]


def mobius_signed(divisors_factor_count: dict[int, int]) -> int:
    """μ(d) = (-1)^ω(d) for squarefree d."""
    # Not used directly; we track sign via length of factor list.
    raise NotImplementedError


def gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return a


# ---------------------------------------------------------------------------
# Maynard weight evaluator
# ---------------------------------------------------------------------------


@dataclass
class MaynardConfig:
    k: int
    H: list[int]
    R: float  # truncation threshold for d_1·...·d_k
    F_choice: str = "selberg_gpy"  # "selberg_gpy" | "linear" | "constant"

    def F(self, *xs: float) -> float:
        """Smooth function on the simplex Σ x_i ≤ 1.
        F=0 outside simplex; we enforce that via d_1...d_k ≤ R.
        """
        s = sum(xs)
        if s > 1.0:
            return 0.0
        if self.F_choice == "constant":
            return 1.0
        if self.F_choice == "selberg_gpy":
            # GPY weight (1 - Σ x_i)^2: standard Selberg majorant
            return (1.0 - s) ** 2
        if self.F_choice == "linear":
            # Linear weight (1 - Σ x_i): less standard, simpler
            return 1.0 - s
        if self.F_choice == "maynard_sym":
            # Maynard-style symmetric weight: (1 - Σx)^2 * (1 + a·Σx + b·Σ(x_i x_j))
            # Use a=0.5, b=2.0 (rough optimization for k=3)
            t = sum(xs)
            cross = 0.0
            for i in range(len(xs)):
                for j in range(i + 1, len(xs)):
                    cross += xs[i] * xs[j]
            return (1.0 - t) ** 2 * (1.0 + 0.5 * t + 2.0 * cross)
        raise ValueError(f"unknown F_choice {self.F_choice}")


def compute_maynard_weight(
    n: int,
    cfg: MaynardConfig,
    spf: list[int],
    log_R: float,
) -> float:
    """Evaluate w(n) = (Σ_{admissible (d_i)} μ(d_1)...μ(d_k) F(...))^2."""
    R = cfg.R
    # Per-coordinate: squarefree divisors of n+h_i, each with d_i ≤ R
    div_lists: list[list[tuple[int, int]]] = []
    for h in cfg.H:
        m = n + h
        ds = squarefree_divisors_up_to(m, int(R), spf)
        # decorate with omega(d) for sign
        decorated: list[tuple[int, int]] = []
        for d in ds:
            # count primes in factorization (squarefree → ω = number of prime factors)
            x = d
            omega = 0
            while x > 1 and x < len(spf):
                p = spf[x]
                while x % p == 0:
                    x //= p
                omega += 1
            decorated.append((d, omega))
        div_lists.append(decorated)

    S = 0.0
    # Enumerate triples with d_1·...·d_k ≤ R AND pairwise coprime
    for tup in product(*div_lists):
        ds = tuple(d for d, _ in tup)
        prod = 1
        ok = True
        for d in ds:
            prod *= d
            if prod > R:
                ok = False
                break
        if not ok:
            continue
        # Pairwise coprime
        coprime = True
        for i in range(len(ds)):
            for j in range(i + 1, len(ds)):
                if gcd(ds[i], ds[j]) != 1:
                    coprime = False
                    break
            if not coprime:
                break
        if not coprime:
            continue
        # μ-sign and F-value
        sign = 1
        for _, omega in tup:
            if omega & 1:
                sign = -sign
        if log_R > 0:
            xs = tuple(math.log(d) / log_R if d > 1 else 0.0 for d in ds)
        else:
            xs = tuple(0.0 for _ in ds)
        f_val = cfg.F(*xs)
        S += sign * f_val

    return S * S


def count_ops_for_evaluation(
    n: int,
    cfg: MaynardConfig,
    spf: list[int],
) -> dict:
    """Count operations to evaluate w(n) (admissible-tuple enumeration)."""
    R = cfg.R
    div_lists = []
    for h in cfg.H:
        m = n + h
        ds = squarefree_divisors_up_to(m, int(R), spf)
        div_lists.append(ds)
    # Total tuple count (after R-truncation)
    n_tuples_raw = 1
    for ds in div_lists:
        n_tuples_raw *= len(ds)
    n_tuples_simplex = 0
    n_tuples_coprime = 0
    for tup in product(*div_lists):
        prod = 1
        ok = True
        for d in tup:
            prod *= d
            if prod > R:
                ok = False
                break
        if not ok:
            continue
        n_tuples_simplex += 1
        # Coprime?
        coprime = True
        for i in range(len(tup)):
            for j in range(i + 1, len(tup)):
                if gcd(tup[i], tup[j]) != 1:
                    coprime = False
                    break
            if not coprime:
                break
        if coprime:
            n_tuples_coprime += 1
    return {
        "n": n,
        "div_counts": [len(ds) for ds in div_lists],
        "n_tuples_raw": n_tuples_raw,
        "n_tuples_simplex": n_tuples_simplex,
        "n_tuples_coprime": n_tuples_coprime,
    }


# ---------------------------------------------------------------------------
# Classification: does (n, n+2, n+6) contain a prime?
# ---------------------------------------------------------------------------


def classify_window(n: int, H: list[int], is_prime: list[bool]) -> dict:
    """How many of n+h_i (h_i in H) are prime?"""
    count = sum(1 for h in H if is_prime[n + h])
    return {
        "n": n,
        "primes_in_window": count,
        "any_prime": count >= 1,
        "primality_pattern": tuple(is_prime[n + h] for h in H),
    }


# ---------------------------------------------------------------------------
# AUC and ROC
# ---------------------------------------------------------------------------


def mann_whitney_auc(scores_pos: list[float], scores_neg: list[float]) -> float:
    """AUC = P(score_pos > score_neg)."""
    n_pos = len(scores_pos)
    n_neg = len(scores_neg)
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    # Combined ranking
    combined = [(s, 1) for s in scores_pos] + [(s, 0) for s in scores_neg]
    combined.sort()
    ranks = {}
    i = 0
    while i < len(combined):
        j = i
        while j < len(combined) and combined[j][0] == combined[i][0]:
            j += 1
        avg_rank = 0.5 * (i + 1 + j)  # 1-indexed mean rank
        for k in range(i, j):
            ranks[k] = avg_rank
        i = j
    sum_ranks_pos = 0.0
    for k, (_, lab) in enumerate(combined):
        if lab == 1:
            sum_ranks_pos += ranks[k]
    U = sum_ranks_pos - n_pos * (n_pos + 1) / 2.0
    return U / (n_pos * n_neg)


def roc_curve(scores_pos: list[float], scores_neg: list[float], n_points: int = 200) -> dict:
    all_scores = sorted(set(scores_pos + scores_neg), reverse=True)
    if len(all_scores) > n_points:
        idx = [int(i * (len(all_scores) - 1) / (n_points - 1)) for i in range(n_points)]
        thresholds = [all_scores[i] for i in idx]
    else:
        thresholds = all_scores

    pos_sorted = sorted(scores_pos, reverse=True)
    neg_sorted = sorted(scores_neg, reverse=True)
    n_pos = len(pos_sorted)
    n_neg = len(neg_sorted)

    points = []
    best_f1 = -1.0
    best_f1_tau = float("nan")
    best_acc = -1.0
    best_acc_tau = float("nan")
    best_pr = -1.0
    best_pr_tau = float("nan")

    for tau in thresholds:
        tp = sum(1 for s in pos_sorted if s >= tau)
        fp = sum(1 for s in neg_sorted if s >= tau)
        tn = n_neg - fp
        fn = n_pos - tp
        if tp + fp > 0:
            precision = tp / (tp + fp)
        else:
            precision = 0.0
        if tp + fn > 0:
            recall = tp / (tp + fn)
        else:
            recall = 0.0
        tpr = recall
        fpr = fp / n_neg if n_neg > 0 else 0.0
        accuracy = (tp + tn) / (n_pos + n_neg)
        if precision + recall > 0:
            f1 = 2 * precision * recall / (precision + recall)
        else:
            f1 = 0.0
        pr_product = precision * recall
        if f1 > best_f1:
            best_f1, best_f1_tau = f1, tau
        if accuracy > best_acc:
            best_acc, best_acc_tau = accuracy, tau
        if pr_product > best_pr:
            best_pr, best_pr_tau = pr_product, tau
        points.append({
            "tau": tau,
            "tp": tp, "fp": fp, "tn": tn, "fn": fn,
            "precision": precision, "recall": recall,
            "f1": f1, "accuracy": accuracy, "tpr": tpr, "fpr": fpr,
        })

    return {
        "points": points,
        "best_f1": best_f1, "best_f1_tau": best_f1_tau,
        "best_accuracy": best_acc, "best_accuracy_tau": best_acc_tau,
        "best_pr_product": best_pr, "best_pr_product_tau": best_pr_tau,
    }


# ---------------------------------------------------------------------------
# Main experiment driver
# ---------------------------------------------------------------------------


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, default=10_000, help="window start")
    p.add_argument("--window", type=int, default=5_000, help="window size")
    p.add_argument("--theta", type=float, default=0.25, help="R = N^theta")
    p.add_argument("--F", type=str, default="selberg_gpy",
                   choices=["selberg_gpy", "linear", "constant", "maynard_sym"])
    p.add_argument("--Hs", type=str, default="0,2,6", help="comma-separated H tuple")
    p.add_argument("--out", type=str, default="results.json")
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--quick", action="store_true", help="small mode for debugging")
    p.add_argument("--ops_sample", type=int, default=200,
                   help="sample size for op-count statistics")
    args = p.parse_args()

    random.seed(args.seed)
    H = [int(x) for x in args.Hs.split(",")]
    k = len(H)
    N = args.N
    window = args.window
    theta = args.theta
    if args.quick:
        window = min(window, 500)

    R = max(2.0, N**theta)
    log_R = math.log(R) if R > 1 else 1.0

    n_min = N
    n_max = N + window - 1
    sieve_max = n_max + max(H) + 100  # for primality
    spf_max = sieve_max  # for divisor enumeration of n+h_i

    print(f"[setup] N={N}, window={window}, theta={theta}, R={R:.2f}, "
          f"H={H}, k={k}, F_choice={args.F}")
    print(f"[setup] sieve up to {sieve_max}")
    t0 = time.time()
    is_prime = sieve_of_eratosthenes(sieve_max)
    spf = smallest_prime_factor(spf_max)
    print(f"[setup] sieve complete in {time.time()-t0:.2f}s")

    cfg = MaynardConfig(k=k, H=H, R=R, F_choice=args.F)

    # Compute w(n) for n in [n_min, n_max]
    weights = []
    classes = []
    counts_by_pattern = defaultdict(int)
    t0 = time.time()
    for n in range(n_min, n_max + 1):
        w = compute_maynard_weight(n, cfg, spf, log_R)
        cls = classify_window(n, H, is_prime)
        weights.append(w)
        classes.append(cls)
        counts_by_pattern[cls["primality_pattern"]] += 1
    elapsed = time.time() - t0
    print(f"[weights] computed {len(weights)} weights in {elapsed:.2f}s "
          f"({1000*elapsed/len(weights):.2f} ms/n)")

    # Op-count sampling
    print("[ops] sampling operation counts ...")
    op_samples = []
    sample_idxs = random.sample(range(n_min, n_max + 1),
                                min(args.ops_sample, n_max - n_min + 1))
    for n in sample_idxs:
        op_samples.append(count_ops_for_evaluation(n, cfg, spf))

    # Aggregate AUCs
    pos_any = [w for w, c in zip(weights, classes) if c["any_prime"]]
    neg_any = [w for w, c in zip(weights, classes) if not c["any_prime"]]
    pos_full = [w for w, c in zip(weights, classes) if c["primes_in_window"] >= 1]
    neg_full = [w for w, c in zip(weights, classes) if c["primes_in_window"] == 0]
    pos_strict = [w for w, c in zip(weights, classes) if c["primes_in_window"] >= 2]
    neg_strict = [w for w, c in zip(weights, classes) if c["primes_in_window"] <= 1]

    auc_any = mann_whitney_auc(pos_any, neg_any)
    auc_strict = mann_whitney_auc(pos_strict, neg_strict)

    # Per-position AUC (is n+h_i prime?)
    per_position_auc = {}
    for i, h in enumerate(H):
        pos_i = [w for w, c in zip(weights, classes) if c["primality_pattern"][i]]
        neg_i = [w for w, c in zip(weights, classes) if not c["primality_pattern"][i]]
        per_position_auc[h] = {
            "auc": mann_whitney_auc(pos_i, neg_i),
            "n_pos": len(pos_i),
            "n_neg": len(neg_i),
        }

    # Aggregate-positivity check (Maynard's actual claim):
    # Σ w(n)·χ_P(n+h_i) > c · Σ w(n) for some i?
    sum_w = sum(weights)
    aggregate = {}
    for i, h in enumerate(H):
        sum_wprime = sum(w for w, c in zip(weights, classes) if c["primality_pattern"][i])
        aggregate[h] = {
            "ratio": sum_wprime / sum_w if sum_w > 0 else float("nan"),
            "sum_w_chi_p": sum_wprime,
        }
    aggregate_max_ratio = max(d["ratio"] for d in aggregate.values())

    # ROC for "any prime in window"
    roc = roc_curve(pos_any, neg_any, n_points=400)

    # Statistics on weights
    def stats(xs):
        if not xs:
            return {"n": 0}
        s = sorted(xs)
        m = sum(xs) / len(xs)
        v = sum((x - m) ** 2 for x in xs) / len(xs)
        return {
            "n": len(xs),
            "mean": m,
            "std": math.sqrt(v),
            "min": s[0],
            "p10": s[len(s) // 10] if len(s) >= 10 else s[0],
            "median": s[len(s) // 2],
            "p90": s[(9 * len(s)) // 10] if len(s) >= 10 else s[-1],
            "max": s[-1],
            "frac_zero": sum(1 for x in xs if x == 0.0) / len(xs),
        }

    op_div_counts = [op["n_tuples_coprime"] for op in op_samples]
    op_raw = [op["n_tuples_raw"] for op in op_samples]

    out = {
        "config": {
            "N": N, "window": window, "theta": theta, "R": R,
            "log_R": log_R, "H": H, "k": k, "F_choice": args.F,
            "seed": args.seed,
        },
        "weights_stats": stats(weights),
        "weights_pos_any_stats": stats(pos_any),
        "weights_neg_any_stats": stats(neg_any),
        "auc_any_prime_in_window": auc_any,
        "auc_strict_2plus_primes": auc_strict,
        "per_position_auc": per_position_auc,
        "aggregate_positivity": aggregate,
        "aggregate_max_ratio": aggregate_max_ratio,
        "primality_pattern_counts": {str(k): v for k, v in counts_by_pattern.items()},
        "best_f1": roc["best_f1"],
        "best_f1_tau": roc["best_f1_tau"],
        "best_accuracy": roc["best_accuracy"],
        "best_accuracy_tau": roc["best_accuracy_tau"],
        "best_pr_product": roc["best_pr_product"],
        "best_pr_product_tau": roc["best_pr_product_tau"],
        "op_count_stats": {
            "raw_tuple_count": stats([float(x) for x in op_raw]),
            "coprime_simplex_tuple_count": stats([float(x) for x in op_div_counts]),
        },
        "elapsed_weights_sec": elapsed,
        "ms_per_n": 1000 * elapsed / len(weights),
    }

    out_path = Path(args.out)
    out_path.write_text(json.dumps(out, indent=2, default=str))
    print(f"[done] wrote {out_path}")
    print(f"[summary] AUC(any prime in window) = {auc_any:.4f}")
    print(f"[summary] AUC(strict, 2+ primes)   = {auc_strict:.4f}")
    print(f"[summary] aggregate max ratio      = {aggregate_max_ratio:.4f}")
    print(f"[summary] mean weight              = {out['weights_stats']['mean']:.4g}")
    print(f"[summary] mean coprime tuple count = "
          f"{out['op_count_stats']['coprime_simplex_tuple_count']['mean']:.2f}")
    print(f"[summary] best F1                  = {roc['best_f1']:.4f} at tau="
          f"{roc['best_f1_tau']:.4g}")


if __name__ == "__main__":
    main()
