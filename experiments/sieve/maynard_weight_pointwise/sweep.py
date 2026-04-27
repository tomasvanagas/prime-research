"""Parameter sweep over (N, theta, F) for Maynard weight pointwise AUC."""
from __future__ import annotations
import json
import math
import subprocess
import sys
from pathlib import Path

CONFIGS = []
for N in [10_000, 31_623, 100_000]:
    for theta in [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40]:
        for F in ["selberg_gpy", "linear", "constant", "maynard_sym"]:
            CONFIGS.append({"N": N, "theta": theta, "F": F})

# Also test other admissible H tuples
H_OTHER = [
    "0,4,6",        # admissible 3-tuple, narrower
    "0,2,6,8,12",   # admissible 5-tuple
]
for N in [10_000, 100_000]:
    for theta in [0.20, 0.30]:
        for Hs in H_OTHER:
            CONFIGS.append({"N": N, "theta": theta, "F": "selberg_gpy", "Hs": Hs})

results = []
window_for_N = {10_000: 5_000, 31_623: 6_000, 100_000: 10_000}

for i, cfg in enumerate(CONFIGS):
    N = cfg["N"]
    theta = cfg["theta"]
    F = cfg["F"]
    Hs = cfg.get("Hs", "0,2,6")
    window = window_for_N[N]
    out_name = f"sweep_N{N}_t{theta:.2f}_F{F}_H{Hs.replace(',', '-')}.json"
    cmd = [
        sys.executable, "maynard_weight_pointwise.py",
        "--N", str(N), "--window", str(window),
        "--theta", str(theta), "--F", F, "--Hs", Hs,
        "--out", out_name, "--ops_sample", "100",
    ]
    print(f"[{i+1}/{len(CONFIGS)}] {' '.join(cmd[2:])}", flush=True)
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if p.returncode != 0:
        print(f"  FAILED: {p.stderr[:200]}", flush=True)
        continue
    data = json.loads(Path(out_name).read_text())
    results.append({
        "config": {"N": N, "theta": theta, "F": F, "Hs": Hs, "window": window},
        "R": data["config"]["R"],
        "auc_any": data["auc_any_prime_in_window"],
        "auc_strict": data["auc_strict_2plus_primes"],
        "agg_max_ratio": data["aggregate_max_ratio"],
        "mean_w": data["weights_stats"]["mean"],
        "best_f1": data["best_f1"],
        "best_acc": data["best_accuracy"],
        "ms_per_n": data["ms_per_n"],
        "mean_tuples": data["op_count_stats"]["coprime_simplex_tuple_count"]["mean"],
        "p90_tuples": data["op_count_stats"]["coprime_simplex_tuple_count"]["p90"],
    })

Path("sweep_summary.json").write_text(json.dumps(results, indent=2, default=str))
print("\n--- Summary ---")
print(f"{'N':>8} {'theta':>5} {'F':>15} {'H':>15} {'R':>7} {'AUC_any':>8} {'AUC_str':>8} "
      f"{'agg':>6} {'F1':>6} {'tup_p90':>7}")
for r in results:
    c = r["config"]
    print(f"{c['N']:>8} {c['theta']:>5.2f} {c['F']:>15} {c['Hs']:>15} "
          f"{r['R']:>7.2f} {r['auc_any']:>8.4f} {r['auc_strict']:>8.4f} "
          f"{r['agg_max_ratio']:>6.3f} {r['best_f1']:>6.3f} {r['p90_tuples']:>7.0f}")
