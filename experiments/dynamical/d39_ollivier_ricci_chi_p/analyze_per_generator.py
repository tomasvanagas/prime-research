"""Inspect per-generator κ for the prime-Cayley graph: does small-prime
generator (2, 3, ...) dominate, or is curvature structure-blind to which
prime is which?"""
import json
from pathlib import Path
import numpy as np

with open(Path(__file__).parent / "d39_results.json") as f:
    data = json.load(f)

print(
    f"{'N':>4} {'c':>4} {'g_canon':>8} {'kappa':>9} | sorted top 3, bot 3 κ across primes"
)
print("-" * 90)
for run in data["runs"]:
    N, c = run["N"], run["c"]
    p = run.get("prime")
    if not p:
        continue
    kpg = p["kappa_per_g"]  # list of [g, kappa]
    # Canonicalise g -> min(g, N-g)
    canon = [(min(g, N - g), k) for g, k in kpg]
    canon.sort()
    arr = np.array(canon)
    print(f"\n--- N={N}, c={c}, |S|={run['S_P_size']} ---")
    print(f"  κ̄={p['kappa_mean']:+.4f}, κ_min={p['kappa_min']:+.4f}")
    print(f"  per-g (sorted by g_canon):")
    for g, k in canon:
        print(f"    g={g:>4d}  κ={k:+.4f}")
    if len(canon) >= 6:
        sorted_by_k = sorted(canon, key=lambda x: x[1])
        print("  bottom 3:", [f"g={g}, κ={k:+.3f}" for g, k in sorted_by_k[:3]])
        print("  top 3:   ", [f"g={g}, κ={k:+.3f}" for g, k in sorted_by_k[-3:]])
