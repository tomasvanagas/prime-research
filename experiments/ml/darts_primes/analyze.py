"""
Analysis of darts_primes results.

Computes:
  - Mean / std final BCE loss for PRIMES vs control across seeds
  - Welch t-test for difference of means
  - Mann-Whitney U (rank-based, non-parametric)
  - Discrete-circuit accuracy comparison
  - Best PRIMES circuit description (gates / connectivity)
  - Extrapolation accuracy on held-out [2^N, 2^N + 1000)
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats


def summarize(seeds_data, label):
    final_losses = np.array([s["final_loss"] for s in seeds_data])
    min_losses = np.array([s["min_loss"] for s in seeds_data])
    soft_acc = np.array([s["soft_acc"] for s in seeds_data])
    discrete_acc = np.array([s["discrete_acc"] for s in seeds_data])
    print(f"  {label}:")
    print(f"    final_loss: mean={final_losses.mean():.4f} std={final_losses.std(ddof=1):.4f} min={final_losses.min():.4f}")
    print(f"    min_loss:   mean={min_losses.mean():.4f} std={min_losses.std(ddof=1):.4f} min={min_losses.min():.4f}")
    print(f"    soft_acc:   mean={soft_acc.mean():.4f} std={soft_acc.std(ddof=1):.4f} max={soft_acc.max():.4f}")
    print(f"    discrete_acc: mean={discrete_acc.mean():.4f} std={discrete_acc.std(ddof=1):.4f} max={discrete_acc.max():.4f}")
    return {
        "final_losses": final_losses,
        "min_losses": min_losses,
        "soft_acc": soft_acc,
        "discrete_acc": discrete_acc,
    }


def main(results_path: str):
    with open(results_path) as f:
        R = json.load(f)
    cfg = R["config"]
    print(f"=== DARTS-PRIMES results ({results_path}) ===")
    print(f"Config: N={cfg['N']} G1={cfg['G1']} G2={cfg['G2']} epochs={cfg['epochs']} lr={cfg['lr']} n_seeds={cfg['n_seeds']}")
    print()

    P = summarize(R["primes"], "PRIMES")
    C = summarize(R["controls"], "CONTROL (matched-density random)")
    print()

    # Density-floor baseline.
    y_train = R["primes"][0].get("y_train")
    if y_train is None:
        # Could be excluded from JSON; recompute.
        from sympy import isprime
        n_all = np.arange(2 ** cfg["N"], dtype=np.int64)
        y_train = [(1 if isprime(int(n)) else 0) for n in n_all]
    y_train_np = np.asarray(y_train, dtype=np.float32)
    p_density = float(y_train_np.mean())
    h_density = -p_density * np.log(p_density + 1e-12) - (1 - p_density) * np.log(1 - p_density + 1e-12)
    baseline_acc = max(p_density, 1 - p_density)
    print(f"Baseline (density-only):")
    print(f"  prime density:   {p_density:.4f}")
    print(f"  entropy floor:   {h_density:.4f}  (constant predictor BCE)")
    print(f"  majority-class accuracy: {baseline_acc:.4f}")
    print()

    # Welch t-test on final loss.
    t_stat, t_p = stats.ttest_ind(P["final_losses"], C["final_losses"], equal_var=False)
    u_stat, u_p = stats.mannwhitneyu(P["final_losses"], C["final_losses"], alternative="two-sided")
    print(f"Welch t-test PRIMES vs CONTROL final_loss:  t={t_stat:.3f}  p={t_p:.4f}")
    print(f"Mann-Whitney U PRIMES vs CONTROL final_loss: U={u_stat:.1f} p={u_p:.4f}")
    print()

    t_stat2, t_p2 = stats.ttest_ind(P["soft_acc"], C["soft_acc"], equal_var=False)
    print(f"Welch t-test PRIMES vs CONTROL soft_acc:    t={t_stat2:.3f}  p={t_p2:.4f}")

    t_stat3, t_p3 = stats.ttest_ind(P["discrete_acc"], C["discrete_acc"], equal_var=False)
    print(f"Welch t-test PRIMES vs CONTROL discrete_acc:t={t_stat3:.3f}  p={t_p3:.4f}")
    print()

    # Best PRIMES architecture.
    best_idx = int(np.argmin([s["final_loss"] for s in R["primes"]]))
    best = R["primes"][best_idx]
    print(f"Best PRIMES seed: idx={best_idx}, final_loss={best['final_loss']:.4f}, soft_acc={best['soft_acc']:.4f}, discrete_acc={best['discrete_acc']:.4f}")
    desc = best["discrete_desc"]
    L1_ops = [n["op"] for n in desc["layer1"]]
    L2_ops = [n["op"] for n in desc["layer2"]]
    L3_op = desc["layer3"]["op"]
    from collections import Counter
    print(f"  Layer 1 ops (count): {dict(Counter(L1_ops))}")
    print(f"  Layer 2 ops (count): {dict(Counter(L2_ops))}")
    print(f"  Layer 3 op: {L3_op}")
    sel_sizes_l1 = [len(n["selected"]) for n in desc["layer1"]]
    sel_sizes_l2 = [len(n["selected"]) for n in desc["layer2"]]
    sel_size_l3 = len(desc["layer3"]["selected"])
    print(f"  Layer 1 |selected| stats: mean={np.mean(sel_sizes_l1):.1f} max={max(sel_sizes_l1)}")
    print(f"  Layer 2 |selected| stats: mean={np.mean(sel_sizes_l2):.1f} max={max(sel_sizes_l2)}")
    print(f"  Layer 3 |selected| = {sel_size_l3}")
    print()

    # Best CONTROL architecture for sanity.
    best_c_idx = int(np.argmin([s["final_loss"] for s in R["controls"]]))
    best_c = R["controls"][best_c_idx]
    desc_c = best_c["discrete_desc"]
    L1_ops_c = [n["op"] for n in desc_c["layer1"]]
    L2_ops_c = [n["op"] for n in desc_c["layer2"]]
    print(f"Best CONTROL seed: idx={best_c_idx}, final_loss={best_c['final_loss']:.4f}")
    print(f"  Layer 1 ops: {dict(Counter(L1_ops_c))}")
    print(f"  Layer 2 ops: {dict(Counter(L2_ops_c))}")
    print()

    # Extrapolation.
    extr = R["extrap"].get("primes_best", {})
    print(f"Extrapolation (best PRIMES discrete circuit on n in [{extr.get('n_lo')}, {extr.get('n_hi')})):")
    print(f"  discrete_acc on extrap window: {extr.get('discrete_acc'):.4f}")
    print(f"  prime density on extrap window: {extr.get('y_density'):.4f}")
    extr_baseline = max(extr.get('y_density', 0), 1 - extr.get('y_density', 1))
    print(f"  majority-class baseline on extrap: {extr_baseline:.4f}")
    print()

    # Verdict bookkeeping.
    print("--- Verdict markers ---")
    primes_better = bool(P['final_losses'].mean() < C['final_losses'].mean())
    primes_above_baseline = bool(P['discrete_acc'].mean() > baseline_acc)
    print(f"PRIMES outperformed CONTROL final_loss? {primes_better} ({P['final_losses'].mean():.4f} vs {C['final_losses'].mean():.4f})")
    print(f"PRIMES discrete_acc > majority-class baseline? {primes_above_baseline} ({P['discrete_acc'].mean():.4f} vs {baseline_acc:.4f})")
    if cfg.get("N") and extr:
        ext_above = bool(extr.get('discrete_acc', 0) > extr_baseline)
        print(f"Extrapolation discrete_acc > majority-class baseline? {ext_above} ({extr.get('discrete_acc', 0):.4f} vs {extr_baseline:.4f})")


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "/apps/aplikacijos/prime-research/experiments/ml/darts_primes/run/results.json"
    main(path)
