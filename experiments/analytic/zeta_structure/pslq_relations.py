"""
PSLQ/LLL search for integer linear relations among Riemann zeta zeros.

Searches for:
1. Relations among small subsets (size 3,4,5) of the first 30 zeros,
   augmented with {1, pi, log(2*pi)}.
2. Pairwise relations: a*gamma_i + b*gamma_j + c*pi + d*log(2*pi) + e = 0
   for all pairs from the first 50 zeros.
3. Whether linear combinations of K zeros (K=5,10,20,50) approximate
   "nice" constants better than random.

Uses mpmath with 50+ decimal digits of precision.
"""

import sys
import itertools
import random
import time
import os

import mpmath

mpmath.mp.dps = 60  # 60 decimal digits

# ── Load zeros ──────────────────────────────────────────────────────────
DATA_PATH = "/apps/aplikacijos/prime-research/data/zeta_zeros_1000.txt"
OUT_DIR = "/apps/aplikacijos/prime-research/experiments/analytic/zeta_structure"

with open(DATA_PATH) as f:
    zeros = [mpmath.mpf(line.strip()) for line in f if line.strip()]

results = []

def log(msg):
    print(msg, flush=True)
    results.append(msg)

log(f"Loaded {len(zeros)} zeta zeros, precision = {mpmath.mp.dps} digits")
log(f"First 5: {[float(z) for z in zeros[:5]]}")

# Nice constants for augmenting PSLQ vectors
CONSTANTS = {
    "1": mpmath.mpf(1),
    "pi": mpmath.pi,
    "log(2*pi)": mpmath.log(2 * mpmath.pi),
}

MAX_COEFF = 1000

augment_names = ["1", "pi", "log(2*pi)"]
augment_vals = [CONSTANTS[n] for n in augment_names]

# ═══════════════════════════════════════════════════════════════════════
# PART 1: PSLQ on subsets of size 3, 4, 5 from first N zeros
#   Size 3: first 30 zeros -> C(30,3) = 4060
#   Size 4: first 20 zeros -> C(20,4) = 4845
#   Size 5: first 15 zeros -> C(15,5) = 3003
# ═══════════════════════════════════════════════════════════════════════
log("\n" + "=" * 72)
log("PART 1: PSLQ on subsets of zeros (sizes 3, 4, 5)")
log("  Each vector augmented with [1, pi, log(2*pi)]")
log("=" * 72)

relations_found_p1 = 0
total_tested_p1 = 0
t0 = time.time()

subset_configs = [
    (3, 30),  # C(30,3) = 4060
    (4, 20),  # C(20,4) = 4845
    (5, 15),  # C(15,5) = 3003
]

for subset_size, n_zeros in subset_configs:
    combos = list(itertools.combinations(range(n_zeros), subset_size))
    log(f"\n--- Subset size {subset_size} from first {n_zeros} zeros: {len(combos)} combinations ---")
    found_this_size = 0

    for idx, combo in enumerate(combos):
        if idx % 500 == 0 and idx > 0:
            print(f"  ... {idx}/{len(combos)} tested", flush=True)
        total_tested_p1 += 1
        vec = list(augment_vals) + [zeros[i] for i in combo]
        try:
            rel = mpmath.pslq(vec, maxcoeff=MAX_COEFF, maxsteps=300)
        except Exception:
            rel = None

        if rel is not None:
            names = augment_names + [f"gamma_{i+1}" for i in combo]
            residual = sum(c * v for c, v in zip(rel, vec))
            if abs(residual) < mpmath.mpf(10) ** (-40):
                terms = [f"({c})*{n}" for c, n in zip(rel, names) if c != 0]
                log(f"  RELATION FOUND: {' + '.join(terms)} = 0")
                log(f"    Residual: {float(residual):.2e}, Max|c|: {max(abs(c) for c in rel)}")
                relations_found_p1 += 1
                found_this_size += 1

    log(f"  Found {found_this_size} relations for size {subset_size}")

elapsed1 = time.time() - t0
log(f"\nPart 1 total: {relations_found_p1} relations from {total_tested_p1} tests ({elapsed1:.1f}s)")

# ═══════════════════════════════════════════════════════════════════════
# PART 2: Pairwise relations among first 50 zeros
#   a*gamma_i + b*gamma_j + c*pi + d*log(2*pi) + e = 0
# ═══════════════════════════════════════════════════════════════════════
log("\n" + "=" * 72)
log("PART 2: Pairwise PSLQ for first 50 zeros")
log("  Vector: [gamma_i, gamma_j, pi, log(2*pi), 1]")
log("=" * 72)

relations_found_p2 = 0
total_tested_p2 = 0
t0 = time.time()

pairs = list(itertools.combinations(range(50), 2))
log(f"Testing {len(pairs)} pairs...")

for idx, (i, j) in enumerate(pairs):
    if idx % 200 == 0 and idx > 0:
        print(f"  ... {idx}/{len(pairs)} pairs tested", flush=True)
    total_tested_p2 += 1
    vec = [zeros[i], zeros[j], mpmath.pi, mpmath.log(2 * mpmath.pi), mpmath.mpf(1)]
    try:
        rel = mpmath.pslq(vec, maxcoeff=MAX_COEFF, maxsteps=300)
    except Exception:
        rel = None

    if rel is not None:
        names = [f"gamma_{i+1}", f"gamma_{j+1}", "pi", "log(2*pi)", "1"]
        residual = sum(c * v for c, v in zip(rel, vec))
        if abs(residual) < mpmath.mpf(10) ** (-40):
            terms = [f"({c})*{n}" for c, n in zip(rel, names) if c != 0]
            log(f"  RELATION: {' + '.join(terms)} = 0")
            log(f"    Residual: {float(residual):.2e}")
            relations_found_p2 += 1

elapsed2 = time.time() - t0
log(f"\nPart 2 total: {relations_found_p2} relations from {total_tested_p2} pairs ({elapsed2:.1f}s)")

# ═══════════════════════════════════════════════════════════════════════
# PART 3: Can K-zero linear combinations approximate nice constants?
# ═══════════════════════════════════════════════════════════════════════
log("\n" + "=" * 72)
log("PART 3: Linear combinations of K zeros vs nice constants")
log("  Testing K = 5, 10, 20, 50 with PSLQ")
log("=" * 72)

target_constants = {
    "pi": mpmath.pi,
    "e": mpmath.e,
    "log(2*pi)": mpmath.log(2 * mpmath.pi),
    "euler_gamma": mpmath.euler,
    "pi^2": mpmath.pi ** 2,
    "sqrt(2*pi)": mpmath.sqrt(2 * mpmath.pi),
}

K_values = [5, 10, 20, 50]
NUM_RANDOM_TRIALS = 200
random.seed(42)

for K in K_values:
    log(f"\n--- K = {K} zeros ---")

    for const_name, const_val in target_constants.items():
        vec = [zeros[i] for i in range(K)] + [const_val]
        try:
            rel = mpmath.pslq(vec, maxcoeff=MAX_COEFF, maxsteps=2000)
        except Exception:
            rel = None

        if rel is not None:
            residual = sum(c * v for c, v in zip(rel, vec))
            max_c = max(abs(c) for c in rel)
            if abs(residual) < mpmath.mpf(10) ** (-40):
                log(f"  {const_name}: RELATION FOUND (max|c|={max_c})")
                nonzero = [(c, idx) for idx, c in enumerate(rel) if c != 0]
                if len(nonzero) <= 15:
                    terms = []
                    for c, idx in nonzero:
                        name = f"gamma_{idx+1}" if idx < K else const_name
                        terms.append(f"({c})*{name}")
                    log(f"    {' + '.join(terms)} = 0")
                else:
                    log(f"    ({len(nonzero)} nonzero coeffs, max|c|={max_c})")
                log(f"    Residual: {float(residual):.2e}")
            else:
                log(f"  {const_name}: PSLQ returned relation but residual too large ({float(residual):.2e})")
        else:
            log(f"  {const_name}: no relation found")

    # Random baseline
    log(f"\n  Random baseline (K={K}, {NUM_RANDOM_TRIALS} trials):")
    for const_name, const_val in target_constants.items():
        best_err = mpmath.mpf("inf")
        for _ in range(NUM_RANDOM_TRIALS):
            coeffs = [random.randint(-100, 100) for _ in range(K)]
            if all(c == 0 for c in coeffs):
                continue
            combo = sum(c * zeros[i] for i, c in enumerate(coeffs))
            if const_val != 0:
                ratio = combo / const_val
                nearest = round(float(ratio))
                if nearest != 0:
                    err = abs(combo - nearest * const_val)
                    if err < best_err:
                        best_err = err
        log(f"    {const_name}: best random approx error = {float(best_err):.6e}")

# ═══════════════════════════════════════════════════════════════════════
# PART 4: PSLQ on consecutive differences
# ═══════════════════════════════════════════════════════════════════════
log("\n" + "=" * 72)
log("PART 4: PSLQ on consecutive zero differences")
log("=" * 72)

diffs = [zeros[i+1] - zeros[i] for i in range(29)]
log(f"\nFirst 10 consecutive differences:")
for i in range(10):
    log(f"  gamma_{i+2} - gamma_{i+1} = {float(diffs[i]):.15f}")

log("\nPSLQ on groups of 5 consecutive differences + {1, pi, log(2*pi)}:")
for start in range(0, 25, 5):
    vec = list(augment_vals) + diffs[start:start+5]
    try:
        rel = mpmath.pslq(vec, maxcoeff=MAX_COEFF, maxsteps=500)
    except Exception:
        rel = None
    if rel is not None:
        residual = sum(c * v for c, v in zip(rel, vec))
        if abs(residual) < mpmath.mpf(10) ** (-40):
            names = augment_names + [f"d_{start+i+1}" for i in range(5)]
            terms = [f"({c})*{n}" for c, n in zip(rel, names) if c != 0]
            log(f"  RELATION at start={start}: {' + '.join(terms)} = 0")
        else:
            log(f"  Group start={start}: residual too large ({float(residual):.2e})")
    else:
        log(f"  Group start={start}: no relation found")

# ═══════════════════════════════════════════════════════════════════════
# SUMMARY & SAVE
# ═══════════════════════════════════════════════════════════════════════
log("\n" + "=" * 72)
log("SUMMARY")
log("=" * 72)

summary = f"""
Total relations found:
  Part 1 (subsets of 3,4,5): {relations_found_p1}
  Part 2 (all pairs from first 50 zeros): {relations_found_p2}

Interpretation:
  If NO relations found with |coeff| <= {MAX_COEFF}, this is strong numerical
  evidence that the zeta zeros are linearly independent over Q (as conjectured).
  This means no finite integer combination of zeros yields an algebraic number,
  reinforcing the information-theoretic barrier: each zero carries independent
  information that cannot be compressed via linear algebraic relations.
"""
log(summary)

report_path = os.path.join(OUT_DIR, "pslq_results.md")
with open(report_path, "w") as f:
    f.write("# PSLQ Search for Linear Relations Among Zeta Zeros\n\n")
    f.write(f"**Date:** 2026-04-05\n")
    f.write(f"**Precision:** {mpmath.mp.dps} decimal digits\n")
    f.write(f"**Max coefficient:** {MAX_COEFF}\n")
    f.write(f"**Zeros used:** {len(zeros)} (from zeta_zeros_1000.txt)\n\n")
    f.write("## Full Output\n\n```\n")
    f.write("\n".join(results))
    f.write("\n```\n")

log(f"\nResults saved to {report_path}")
