"""
At N=6, verify the PRIMES-vs-random gap is robust.
Run depth-2 sign-threshold (W=1, full candidate set) on:
  - PRIMES (M=6 expected)
  - 10 random matched (weight 18) controls
  - is_prime (different from pi(x) -- here pi just up to 64)
  - x is multiple-of-3 indicator (control)
"""
import sys, json, random, time
sys.path.insert(0, '.')
from enum_d2_smart import enumerate_w1_thresholds_pruned, depth2_search, primes_table, random_table

N = 6
cands = enumerate_w1_thresholds_pruned(N, k_max=N)
print(f"K = {len(cands)} candidates")

results = {}
# PRIMES
print("\n== PRIMES ==")
p = primes_table(N)
r = depth2_search(p, N, cands, [5, 6], time_limit=180)
results["primes"] = r

# Random with several seeds
results["random"] = {}
for seed in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
    print(f"\n== RANDOM seed {seed} ==")
    rt = random_table(N, sum(p), seed=seed)
    r = depth2_search(rt, N, cands, [5, 6, 7, 8], time_limit=180)
    results["random"][seed] = r

with open("n6_robust.json", "w") as f:
    def cleanup(o):
        if isinstance(o, dict): return {k: cleanup(v) for k, v in o.items()}
        if isinstance(o, list): return [cleanup(v) for v in o]
        if isinstance(o, tuple): return list(o)
        try: json.dumps(o); return o
        except (TypeError, ValueError): return str(o)
    json.dump(cleanup(results), f, indent=2)
print("\nWrote n6_robust.json")

# Summary
print("\n=== Summary ===")
print(f"  PRIMES: min_M = {results['primes'].get('min_M', '?')}")
for s, r in results["random"].items():
    print(f"  random seed {s}: min_M = {r.get('min_M', '?')}")
