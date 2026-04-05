#!/usr/bin/env python3
"""
Kt(delta(n)|n) estimation for prime correction terms.

delta(n) = p(n) - round(R^{-1}(n))

Methods:
  (a) Arithmetic program search (BFS, n=1..500)
  (b) Description complexity via arithmetic expression fits (n=1..10000)
  (c) Conditional compression (block-level)
  (d) Incremental complexity via prediction residual entropy

Goal: determine if Kt grows as O(log n), O(sqrt(n)), or O(n^{1/3}).
"""

import numpy as np
import struct
import bz2
import math
import time
import sys
import os
from collections import deque
from io import BytesIO

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
DELTA_PATH = "/apps/aplikacijos/prime-research/experiments/information_theory/kt_delta_values.npy"
delta_all = np.load(DELTA_PATH)
N_TOTAL = len(delta_all)
print(f"Loaded {N_TOTAL} delta values. Range [{delta_all.min()}, {delta_all.max()}]")

OUT_DIR = "/apps/aplikacijos/prime-research/experiments/information_theory/kt_complexity"
os.makedirs(OUT_DIR, exist_ok=True)

results_lines = []
def log(msg):
    print(msg)
    results_lines.append(msg)

log("=" * 72)
log("Kt(delta(n)|n) ESTIMATION")
log("=" * 72)

# ===================================================================
# (a) Arithmetic program search -- BFS over small programs
# ===================================================================
log("\n--- (a) Arithmetic program search (n=1..500) ---")

def bfs_program_search(target, n, max_depth=6, timeout_s=0.05):
    """
    BFS over arithmetic expressions built from {n, 1, 2, 3, ...small_consts}
    and ops {+, -, *, //, %, neg}.  Returns (program_length, expr_str) or None.

    Program length in bits ~ log2(branching_factor) * depth.
    """
    small_consts = list(range(-10, 11))
    # Initial values on stack: each is (value, depth, description)
    initials = [(n, 0, "n")]
    for c in small_consts:
        initials.append((c, 0, str(c)))

    # BFS: state = (value, depth, description)
    # We want shortest depth to reach target
    seen = set()
    queue = deque()

    for val, d, desc in initials:
        if val == target:
            return (1, desc)
        if val not in seen:
            seen.add(val)
            queue.append((val, d, desc))

    deadline = time.time() + timeout_s

    while queue and time.time() < deadline:
        v1, d1, desc1 = queue.popleft()
        if d1 >= max_depth:
            continue

        for v2_init, _, desc2 in initials:
            new_depth = d1 + 1
            candidates = []
            # Binary ops
            candidates.append((v1 + v2_init, f"({desc1}+{desc2})"))
            candidates.append((v1 - v2_init, f"({desc1}-{desc2})"))
            candidates.append((v2_init - v1, f"({desc2}-{desc1})"))
            if abs(v2_init) <= 100 and abs(v1) <= 100:
                candidates.append((v1 * v2_init, f"({desc1}*{desc2})"))
            if v2_init != 0 and abs(v1) < 10**9:
                candidates.append((v1 // v2_init, f"({desc1}//{desc2})"))
                candidates.append((v1 % v2_init, f"({desc1}%{desc2})"))

            for val, desc in candidates:
                if val == target:
                    return (new_depth, desc)
                if abs(val) < 10**7 and val not in seen:
                    seen.add(val)
                    if new_depth < max_depth:
                        queue.append((val, new_depth, desc))

    return None


N_BFS = 500
bfs_results = np.full(N_BFS, np.nan)
bfs_times = []

t0 = time.time()
found_count = 0
for i in range(N_BFS):
    n = i + 1  # 1-indexed
    target = int(delta_all[i])
    t_start = time.time()
    result = bfs_program_search(target, n, max_depth=5, timeout_s=0.02)
    elapsed = time.time() - t_start
    bfs_times.append(elapsed)
    if result is not None:
        depth, expr = result
        # Kt estimate in bits: each step picks from ~25 ops/consts, so ~log2(25)*depth
        bits = depth * math.log2(25)
        bfs_results[i] = bits
        found_count += 1

t_total = time.time() - t0
log(f"  BFS search done in {t_total:.1f}s. Found programs for {found_count}/{N_BFS} values.")

# For values not found, use ceil(log2(|delta|+1)) + overhead as upper bound
for i in range(N_BFS):
    if np.isnan(bfs_results[i]):
        d = abs(int(delta_all[i]))
        bfs_results[i] = math.ceil(math.log2(d + 2)) + 10  # literal encoding + overhead

# Bin and average
bins_a = [1, 10, 50, 100, 200, 500]
log("  n_range -> avg Kt bits (BFS)")
for i_b in range(len(bins_a) - 1):
    lo, hi = bins_a[i_b], bins_a[i_b + 1]
    vals = bfs_results[lo-1:hi]
    log(f"    n=[{lo},{hi}]: mean={np.nanmean(vals):.2f} bits")


# ===================================================================
# (b) Description complexity via arithmetic expression fits
# ===================================================================
log("\n--- (b) Description complexity via expression fitting (n=1..10000) ---")

N_FIT = 10000
ns = np.arange(1, N_FIT + 1, dtype=float)
ds = delta_all[:N_FIT].astype(float)

# Fit various models and compute residuals -> bits needed
def bits_from_residual(residuals):
    """Given residuals (floats), estimate bits per value via entropy of rounded residuals."""
    r_int = np.round(residuals).astype(int)
    vals, counts = np.unique(r_int, return_counts=True)
    probs = counts / counts.sum()
    entropy = -np.sum(probs * np.log2(probs + 1e-30))
    return entropy

# Model 1: constant (0 params)
resid_const = ds - np.mean(ds)
bits_const = bits_from_residual(resid_const)
log(f"  Model: constant -> residual entropy = {bits_const:.2f} bits/value")

# Model 2: linear in n (2 params)
coeffs = np.polyfit(ns, ds, 1)
resid_lin = ds - np.polyval(coeffs, ns)
bits_lin = bits_from_residual(resid_lin)
log(f"  Model: linear(n) -> residual entropy = {bits_lin:.2f} bits/value")

# Model 3: polynomial deg 5 (6 params)
coeffs5 = np.polyfit(ns, ds, 5)
resid_p5 = ds - np.polyval(coeffs5, ns)
bits_p5 = bits_from_residual(resid_p5)
log(f"  Model: poly(5) -> residual entropy = {bits_p5:.2f} bits/value")

# Model 4: sqrt(n) * sin fit
# delta ~ a * sqrt(n) * sin(b * log(n) + c) + d
from scipy.optimize import curve_fit

def osc_model(n, a, b, c, d):
    return a * np.sqrt(n) * np.sin(b * np.log(n + 1) + c) + d

try:
    popt, _ = curve_fit(osc_model, ns, ds, p0=[1.0, 2.0, 0.0, -10.0], maxfev=5000)
    resid_osc = ds - osc_model(ns, *popt)
    bits_osc = bits_from_residual(resid_osc)
    log(f"  Model: sqrt(n)*sin(b*log(n)+c) -> residual entropy = {bits_osc:.2f} bits/value (4 params)")
except Exception as e:
    bits_osc = bits_const
    log(f"  Model: oscillatory fit failed: {e}")

# Now compute EFFECTIVE Kt per value at different n scales
# Use sliding window entropy estimation
window = 200
step = 50
kt_sliding = []
n_sliding = []

for start in range(0, N_FIT - window, step):
    block = ds[start:start + window]
    # Conditional: how many bits to describe this block given the index range?
    # Use empirical entropy of the block
    vals_int = np.round(block).astype(int)
    uvals, ucounts = np.unique(vals_int, return_counts=True)
    probs = ucounts / ucounts.sum()
    H = -np.sum(probs * np.log2(probs + 1e-30))
    kt_sliding.append(H)
    n_sliding.append(start + window // 2)

n_sliding = np.array(n_sliding)
kt_sliding = np.array(kt_sliding)

# Fit growth rate: Kt ~ alpha * log(n) + beta  vs  alpha * n^gamma
log_n = np.log(n_sliding + 1)
# Linear fit to log(n)
A_log = np.vstack([log_n, np.ones(len(log_n))]).T
sol_log, res_log, _, _ = np.linalg.lstsq(A_log, kt_sliding, rcond=None)
log(f"\n  Sliding-window entropy vs log(n): Kt ~ {sol_log[0]:.4f} * log(n) + {sol_log[1]:.4f}")
if len(res_log) > 0:
    log(f"    R^2 residual sum: {res_log[0]:.4f}")

# Fit Kt ~ a * n^gamma
valid = kt_sliding > 0
log_kt = np.log(kt_sliding[valid] + 1e-10)
log_ns = np.log(n_sliding[valid] + 1)
A_pow = np.vstack([log_ns, np.ones(len(log_ns))]).T
sol_pow, _, _, _ = np.linalg.lstsq(A_pow, log_kt, rcond=None)
gamma = sol_pow[0]
log(f"  Power-law fit: Kt ~ n^{gamma:.4f}")


# ===================================================================
# (c) Conditional compression
# ===================================================================
log("\n--- (c) Conditional compression (block-level) ---")

block_sizes = [100, 500, 1000, 5000, 10000, 50000, 100000]

for k in block_sizes:
    if k > N_TOTAL:
        continue
    block = delta_all[:k]

    # Raw delta bytes
    raw_bytes = block.astype(np.int32).tobytes()
    comp_raw = bz2.compress(raw_bytes, compresslevel=9)
    bits_raw = len(comp_raw) * 8
    bits_per_val_raw = bits_raw / k

    # Conditional: (n, delta(n)) pairs -- n is implicit from position, so same as raw
    # But let's also try: delta as sequence of differences (predictive coding)
    diffs = np.diff(block)
    diff_bytes = diffs.astype(np.int32).tobytes()
    comp_diff = bz2.compress(diff_bytes, compresslevel=9)
    bits_diff = len(comp_diff) * 8
    bits_per_val_diff = bits_diff / (k - 1)

    # Also: delta mod small primes pattern
    # Compress residuals after subtracting running mean
    running_mean = np.cumsum(block) / np.arange(1, k + 1)
    residuals = block - np.round(running_mean)
    res_bytes = residuals.astype(np.int32).tobytes()
    comp_res = bz2.compress(res_bytes, compresslevel=9)
    bits_res = len(comp_res) * 8
    bits_per_val_res = bits_res / k

    log(f"  k={k:>6d}: raw={bits_per_val_raw:.2f} b/v, diff={bits_per_val_diff:.2f} b/v, "
        f"resid={bits_per_val_res:.2f} b/v")

# Key comparison: does knowing n help?
log("\n  Conditional test: compress delta[i] vs compress (i, delta[i]):")
for k in [1000, 10000, 100000]:
    if k > N_TOTAL:
        continue
    block = delta_all[:k]

    # Just deltas
    raw = block.astype(np.int32).tobytes()
    c1 = len(bz2.compress(raw, 9)) * 8

    # Interleaved (i, delta[i])
    pairs = np.column_stack([np.arange(k, dtype=np.int32), block.astype(np.int32)])
    c2 = len(bz2.compress(pairs.tobytes(), 9)) * 8

    # The index sequence is highly compressible, so c2 overhead ~ log(k) per val
    # If delta is INDEPENDENT of n, c2 ~ c1 + k*log2(k) (indices)
    # If delta DEPENDS on n, c2 could be < c1 + k*log2(k)
    idx_only = np.arange(k, dtype=np.int32).tobytes()
    c_idx = len(bz2.compress(idx_only, 9)) * 8

    log(f"  k={k}: delta_only={c1/k:.2f} b/v, pairs={c2/k:.2f} b/v, "
        f"idx_only={c_idx/k:.2f} b/v, saving={(c1 + c_idx - c2)/k:.3f} b/v")


# ===================================================================
# (d) Incremental complexity via prediction residual
# ===================================================================
log("\n--- (d) Incremental complexity (prediction residual entropy) ---")

# Use sliding-window predictor: predict delta(n) from delta(n-W..n-1)
# Residual entropy = estimate of H(delta(n) | past, n)

def incremental_entropy(data, window_sizes=[10, 50, 200]):
    """Estimate bits(delta(n) | delta(n-W..n-1)) for various W."""
    results = {}
    for W in window_sizes:
        residuals = []
        for i in range(W, len(data)):
            # Simple predictor: mean of last W values
            pred = np.mean(data[i-W:i])
            residuals.append(int(round(data[i] - pred)))

        # Entropy of residuals in blocks
        residuals = np.array(residuals)
        block_size = min(1000, len(residuals) // 10)
        if block_size < 100:
            block_size = len(residuals)

        entropies = []
        positions = []
        for start in range(0, len(residuals) - block_size, block_size // 2):
            blk = residuals[start:start + block_size]
            vals, counts = np.unique(blk, return_counts=True)
            probs = counts / counts.sum()
            H = -np.sum(probs * np.log2(probs + 1e-30))
            entropies.append(H)
            positions.append(W + start + block_size // 2)

        results[W] = (np.array(positions), np.array(entropies))
    return results

# Run on full dataset
inc_results = incremental_entropy(delta_all, window_sizes=[10, 50, 200])

for W, (positions, entropies) in inc_results.items():
    log(f"\n  Window W={W}:")
    # Report at different n scales
    checkpoints = [500, 2000, 10000, 50000, 100000]
    for cp in checkpoints:
        if cp > N_TOTAL:
            continue
        mask = (positions > cp * 0.8) & (positions < cp * 1.2)
        if mask.sum() > 0:
            mean_H = np.mean(entropies[mask])
            log(f"    n~{cp}: H(resid) = {mean_H:.2f} bits/value")

# Also: compute the ratio of incremental complexity to log(n) and sqrt(n)
log("\n  Growth rate analysis of incremental entropy (W=50):")
positions, entropies = inc_results[50]
# Downsample for fitting
step = max(1, len(positions) // 100)
pos_fit = positions[::step].astype(float)
ent_fit = entropies[::step]

# Fit: H ~ a * log(n) + b
log_pos = np.log(pos_fit + 1)
A = np.vstack([log_pos, np.ones(len(log_pos))]).T
sol, res, _, _ = np.linalg.lstsq(A, ent_fit, rcond=None)
log(f"  Fit H ~ {sol[0]:.4f} * log(n) + {sol[1]:.4f}")

# Fit: H ~ a * n^gamma
log_ent = np.log(ent_fit + 1e-10)
A2 = np.vstack([np.log(pos_fit + 1), np.ones(len(pos_fit))]).T
sol2, _, _, _ = np.linalg.lstsq(A2, log_ent, rcond=None)
log(f"  Fit H ~ n^{sol2[0]:.6f}  (power law exponent)")

# Ratio tests
log("\n  Ratio H(n) / log(n) at selected n:")
for cp in [500, 2000, 10000, 50000, 90000]:
    mask = (positions > cp * 0.9) & (positions < cp * 1.1)
    if mask.sum() > 0:
        H_val = np.mean(entropies[mask])
        log(f"    n={cp}: H={H_val:.2f}, H/log2(n)={H_val/math.log2(cp):.3f}, "
            f"H/sqrt(n)={H_val/math.sqrt(cp):.5f}, H/n^(1/3)={H_val/cp**(1./3):.4f}")


# ===================================================================
# Summary and growth classification
# ===================================================================
log("\n" + "=" * 72)
log("SUMMARY: Kt(delta(n)|n) Growth Classification")
log("=" * 72)

# Collect key data points
log("\nMethod (a) BFS program search (n=1..500):")
for rng in [(1, 50), (50, 200), (200, 500)]:
    vals = bfs_results[rng[0]-1:rng[1]]
    log(f"  n=[{rng[0]},{rng[1]}]: Kt ~ {np.nanmean(vals):.1f} bits")

log(f"\nMethod (b) Expression fitting:")
log(f"  Sliding entropy growth: {sol_log[0]:.4f} * log(n) + {sol_log[1]:.4f}")
log(f"  Power-law exponent: {gamma:.4f}")

log(f"\nMethod (d) Incremental entropy growth:")
log(f"  Log fit: {sol[0]:.4f} * log(n)")
log(f"  Power-law exponent: {sol2[0]:.6f}")

# Classification
if abs(sol2[0]) < 0.05 and abs(gamma) < 0.05:
    classification = "APPROXIMATELY CONSTANT (O(1))"
elif abs(sol2[0]) < 0.15:
    classification = "SLOW GROWTH, possibly O(log log n)"
elif sol2[0] < 0.35:
    classification = f"CONSISTENT WITH O(log n) or O(n^{sol2[0]:.3f})"
elif sol2[0] < 0.45:
    classification = f"CONSISTENT WITH O(n^{sol2[0]:.3f}), between log(n) and sqrt(n)"
else:
    classification = f"GROWING AS O(n^{sol2[0]:.3f}), near sqrt(n) barrier"

log(f"\nCLASSIFICATION: {classification}")
log(f"\nInterpretation:")
if sol2[0] < 0.15:
    log("  delta(n) has nearly constant conditional entropy per value.")
    log("  This suggests the 'random' part is FIXED complexity, NOT growing with n.")
    log("  However, this measures STATISTICAL entropy, not Kt directly.")
    log("  The per-value entropy being ~constant means each delta(n) carries a fixed")
    log("  number of 'surprise bits', but the TOTAL information grows as O(n).")
elif sol2[0] < 0.35:
    log("  Kt grows slowly, consistent with O(log n).")
    log("  This would be GOOD news for polylog computation.")
else:
    log("  Kt grows polynomially in n, consistent with the sqrt(n) barrier.")
    log("  This is BAD news for polylog computation of individual p(n).")


# ===================================================================
# Plots
# ===================================================================
log("\n--- Generating plots ---")
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # (a) BFS results
    ax = axes[0, 0]
    ax.scatter(range(1, N_BFS + 1), bfs_results, s=1, alpha=0.3, c='blue')
    ax.set_xlabel('n')
    ax.set_ylabel('Kt estimate (bits)')
    ax.set_title('(a) BFS Program Search Kt(delta(n)|n)')
    # Overlay log(n) and sqrt(n) for reference
    nn = np.arange(1, N_BFS + 1)
    ax.plot(nn, 3 * np.log2(nn + 1), 'r-', label='3*log2(n)', alpha=0.7)
    ax.plot(nn, 0.5 * np.sqrt(nn), 'g-', label='0.5*sqrt(n)', alpha=0.7)
    ax.legend()

    # (b) Sliding window entropy
    ax = axes[0, 1]
    ax.plot(n_sliding, kt_sliding, 'b-', alpha=0.5, linewidth=0.5)
    ax.plot(n_sliding, sol_log[0] * np.log(n_sliding + 1) + sol_log[1], 'r--',
            label=f'{sol_log[0]:.3f}*log(n)+{sol_log[1]:.2f}')
    ax.set_xlabel('n')
    ax.set_ylabel('Entropy (bits/value)')
    ax.set_title('(b) Sliding Window Entropy')
    ax.legend()

    # (c) Compression bits per value
    ax = axes[1, 0]
    ks_plot = []
    bpv_raw = []
    bpv_diff = []
    for k in block_sizes:
        if k > N_TOTAL:
            continue
        block = delta_all[:k]
        raw = block.astype(np.int32).tobytes()
        c_r = len(bz2.compress(raw, 9)) * 8 / k
        diffs = np.diff(block).astype(np.int32).tobytes()
        c_d = len(bz2.compress(diffs, 9)) * 8 / (k - 1)
        ks_plot.append(k)
        bpv_raw.append(c_r)
        bpv_diff.append(c_d)
    ax.semilogx(ks_plot, bpv_raw, 'bo-', label='raw delta')
    ax.semilogx(ks_plot, bpv_diff, 'rs-', label='diff-coded')
    ax.set_xlabel('Block size k')
    ax.set_ylabel('Compressed bits/value')
    ax.set_title('(c) Conditional Compression')
    ax.legend()

    # (d) Incremental entropy at different scales
    ax = axes[1, 1]
    for W, (pos, ent) in inc_results.items():
        ax.plot(pos, ent, alpha=0.4, linewidth=0.5, label=f'W={W}')
    ax.set_xlabel('n')
    ax.set_ylabel('H(delta(n)|past) bits/value')
    ax.set_title('(d) Incremental Entropy')
    ax.legend()

    plt.tight_layout()
    plot_path = os.path.join(OUT_DIR, "kt_estimation_plots.png")
    plt.savefig(plot_path, dpi=150)
    log(f"Plots saved to {plot_path}")

    # Additional: growth-rate comparison plot
    fig2, ax2 = plt.subplots(1, 1, figsize=(10, 6))
    positions, entropies = inc_results[50]
    ax2.plot(positions, entropies, 'b-', alpha=0.3, linewidth=0.5, label='H(delta(n)|past, W=50)')
    nn = np.linspace(100, N_TOTAL, 500)
    ax2.plot(nn, sol[0] * np.log(nn) + sol[1], 'r-', linewidth=2,
             label=f'Fit: {sol[0]:.3f}*log(n)+{sol[1]:.2f}')
    ax2.plot(nn, np.full_like(nn, np.mean(entropies)), 'g--', linewidth=2,
             label=f'Constant: {np.mean(entropies):.2f}')
    ax2.set_xlabel('n')
    ax2.set_ylabel('bits per value')
    ax2.set_title('Kt growth rate: log(n) vs constant vs data')
    ax2.legend()
    growth_path = os.path.join(OUT_DIR, "kt_growth_rate.png")
    plt.savefig(growth_path, dpi=150)
    log(f"Growth rate plot saved to {growth_path}")

except ImportError:
    log("  matplotlib not available, skipping plots")
except Exception as e:
    log(f"  Plot error: {e}")

# Save results
results_path = os.path.join(OUT_DIR, "kt_estimation_results.txt")
with open(results_path, 'w') as f:
    f.write('\n'.join(results_lines))
print(f"\nResults saved to {results_path}")
