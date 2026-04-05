"""
mod_patterns.py -- Search for patterns in gamma_n mod various constants.

Loads 1000 zeta zeros and tests their distribution modulo several
mathematically significant constants. Checks equidistribution, circular
statistics, gap structure, and joint distributions.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ── Paths ──────────────────────────────────────────────────────────────
DATA   = Path("/apps/aplikacijos/prime-research/data/zeta_zeros_1000.txt")
OUTDIR = Path(__file__).parent
RESULTS = OUTDIR / "mod_patterns_results.md"

# ── Load zeros ─────────────────────────────────────────────────────────
gammas = np.loadtxt(str(DATA))
N = len(gammas)
print(f"Loaded {N} zeros, range [{gammas[0]:.4f}, {gammas[-1]:.4f}]")

# ── Moduli to test ─────────────────────────────────────────────────────
moduli = {
    "1":             1.0,
    "pi":            np.pi,
    "log(2pi)":      np.log(2 * np.pi),
    "2pi":           2 * np.pi,
    "sqrt(2pi)":     np.sqrt(2 * np.pi),
    "e":             np.e,
    "log(2)":        np.log(2),
    "log(3)":        np.log(3),
    "log(5)":        np.log(5),
    "log(7)":        np.log(7),
}

# ── Helper: circular stats ─────────────────────────────────────────────
def circular_stats(vals, m):
    """Circular mean and variance on [0, m)."""
    theta = 2 * np.pi * vals / m
    C = np.mean(np.cos(theta))
    S = np.mean(np.sin(theta))
    R_bar = np.sqrt(C**2 + S**2)
    circ_mean_angle = np.arctan2(S, C) % (2 * np.pi)
    circ_mean_val = circ_mean_angle * m / (2 * np.pi)
    circ_var = 1 - R_bar  # in [0, 1]
    return circ_mean_val, circ_var, R_bar

def rayleigh_test(vals, m):
    """Rayleigh test for uniformity on circle. Returns (Z, p-value)."""
    theta = 2 * np.pi * vals / m
    n = len(theta)
    C = np.sum(np.cos(theta))
    S = np.sum(np.sin(theta))
    R = np.sqrt(C**2 + S**2)
    Z = R**2 / n
    # p-value approximation
    p = np.exp(-Z) * (1 + (2*Z - Z**2) / (4*n) - (24*Z - 132*Z**2 + 76*Z**3 - 9*Z**4) / (288*n**2))
    p = max(0, min(1, p))
    return Z, p

# ── Part 1 & 2: Per-modulus analysis ──────────────────────────────────
results = {}

for name, m in moduli.items():
    vals = gammas % m

    # KS test against Uniform(0, m)
    ks_stat, ks_p = stats.kstest(vals / m, 'uniform')

    # Circular stats
    c_mean, c_var, R_bar = circular_stats(vals, m)

    # Rayleigh test
    ray_Z, ray_p = rayleigh_test(vals, m)

    results[name] = {
        'm': m,
        'ks_stat': ks_stat,
        'ks_p': ks_p,
        'circ_mean': c_mean,
        'circ_var': c_var,
        'R_bar': R_bar,
        'rayleigh_Z': ray_Z,
        'rayleigh_p': ray_p,
    }

    # Histogram
    fig, ax = plt.subplots(figsize=(8, 4))
    bins = 50
    ax.hist(vals, bins=bins, density=True, alpha=0.7, edgecolor='black', linewidth=0.3)
    ax.axhline(1.0/m, color='red', linestyle='--', label=f'Uniform 1/{m:.4f}')
    ax.set_title(f'gamma_n mod {name}  (N={N})')
    ax.set_xlabel(f'gamma_n mod {name}')
    ax.set_ylabel('Density')
    ax.legend()
    fig.tight_layout()
    safe = name.replace("(", "").replace(")", "").replace("*", "x")
    fig.savefig(str(OUTDIR / f"hist_mod_{safe}.png"), dpi=120)
    plt.close(fig)

    print(f"mod {name:12s}: KS={ks_stat:.4f} p={ks_p:.4f} | "
          f"circ_var={c_var:.4f} R_bar={R_bar:.4f} | "
          f"Rayleigh Z={ray_Z:.4f} p={ray_p:.4f}")

# ── Part 3: Detailed mod-1 analysis ───────────────────────────────────
print("\n=== MOD 1 DETAILED ANALYSIS ===")
frac = gammas % 1.0

# 3a. Weyl's criterion
weyl_results = []
for k in range(1, 51):
    S = np.mean(np.exp(2j * np.pi * k * gammas))
    weyl_results.append((k, abs(S)))
max_weyl = max(weyl_results, key=lambda x: x[1])
print(f"Weyl sums |S_k|: max at k={max_weyl[0]}, |S|={max_weyl[1]:.6f}")
print(f"  Mean |S_k| = {np.mean([w[1] for w in weyl_results]):.6f}")
print(f"  Expected for random: ~1/sqrt(N) = {1/np.sqrt(N):.6f}")

# Plot Weyl sums
fig, ax = plt.subplots(figsize=(10, 4))
ks_w = [w[0] for w in weyl_results]
vals_w = [w[1] for w in weyl_results]
ax.bar(ks_w, vals_w, color='steelblue', alpha=0.8)
ax.axhline(1/np.sqrt(N), color='red', linestyle='--', label=f'1/sqrt(N)={1/np.sqrt(N):.4f}')
ax.set_xlabel('k')
ax.set_ylabel('|S_k|')
ax.set_title('Weyl criterion: |(1/N) sum exp(2 pi i k gamma_n)| for k=1..50')
ax.legend()
fig.tight_layout()
fig.savefig(str(OUTDIR / "weyl_sums_mod1.png"), dpi=120)
plt.close(fig)

# 3b. Discrepancy D_N
frac_sorted = np.sort(frac)
# D_N = sup |F_N(x) - x|
F_N = np.arange(1, N+1) / N
D_plus = np.max(F_N - frac_sorted)
D_minus = np.max(frac_sorted - np.arange(0, N) / N)
D_N = max(D_plus, D_minus)
# Expected for truly random: ~sqrt(ln(ln(N))/(2N)) by law of iterated logarithm
expected_D = np.sqrt(np.log(np.log(N)) / (2 * N))
print(f"Discrepancy D_N = {D_N:.6f}  (expected random ~ {expected_D:.6f})")

# 3c. Gap analysis
gaps = np.diff(frac_sorted)
# Wrap-around gap
gaps_full = np.append(gaps, 1.0 - frac_sorted[-1] + frac_sorted[0])
mean_gap = 1.0 / N
print(f"Gaps: mean={np.mean(gaps_full):.6f} expected={mean_gap:.6f}")
print(f"  min gap={np.min(gaps_full):.6f}  max gap={np.max(gaps_full):.6f}")
print(f"  ratio max/mean={np.max(gaps_full)/mean_gap:.2f}  min/mean={np.min(gaps_full)/mean_gap:.4f}")

# For truly uniform, gaps follow Exp(N) -- check
# Normalize gaps by N: should be ~ Exp(1)
norm_gaps = gaps_full * N
gap_ks, gap_ks_p = stats.kstest(norm_gaps, 'expon')
print(f"  Normalized gaps vs Exp(1): KS={gap_ks:.4f} p={gap_ks_p:.4f}")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].hist(gaps_full, bins=40, density=True, alpha=0.7, edgecolor='black', linewidth=0.3)
axes[0].set_title('Gaps in gamma_n mod 1 (sorted)')
axes[0].set_xlabel('Gap size')
axes[0].set_ylabel('Density')

axes[1].hist(norm_gaps, bins=40, density=True, alpha=0.7, edgecolor='black', linewidth=0.3)
x_exp = np.linspace(0, np.max(norm_gaps), 200)
axes[1].plot(x_exp, np.exp(-x_exp), 'r--', label='Exp(1)')
axes[1].set_title('Normalized gaps (should be Exp(1) if uniform)')
axes[1].set_xlabel('Normalized gap')
axes[1].legend()
fig.tight_layout()
fig.savefig(str(OUTDIR / "gaps_mod1.png"), dpi=120)
plt.close(fig)

# ── Part 4: Joint distribution (mod 1, mod pi) ────────────────────────
print("\n=== JOINT DISTRIBUTION (mod 1, mod pi) ===")
x_joint = gammas % 1.0
y_joint = gammas % np.pi

fig, ax = plt.subplots(figsize=(7, 6))
h = ax.hist2d(x_joint, y_joint, bins=20, cmap='viridis')
plt.colorbar(h[3], ax=ax)
ax.set_xlabel('gamma_n mod 1')
ax.set_ylabel('gamma_n mod pi')
ax.set_title('Joint distribution: gamma_n mod 1 vs gamma_n mod pi')
fig.tight_layout()
fig.savefig(str(OUTDIR / "joint_mod1_modpi.png"), dpi=120)
plt.close(fig)

# Chi-squared test for independence
nbins = 15
H, xedges, yedges = np.histogram2d(x_joint, y_joint, bins=nbins)
# Expected under independence
row_marg = H.sum(axis=1)
col_marg = H.sum(axis=0)
expected = np.outer(row_marg, col_marg) / N
# Mask bins with small expected counts
mask = expected > 5
if mask.sum() > 0:
    chi2 = np.sum((H[mask] - expected[mask])**2 / expected[mask])
    # dof = (rows-1)*(cols-1) for full table, minus excluded cells
    dof = (nbins - 1) * (nbins - 1)
    chi2_p = 1 - stats.chi2.cdf(chi2, dof)
    print(f"Chi-squared independence test: chi2={chi2:.2f} dof={dof} p={chi2_p:.4f}")
else:
    chi2, dof, chi2_p = np.nan, 0, np.nan
    print("Chi-squared: insufficient data in bins")

# Also check correlation
corr, corr_p = stats.pearsonr(x_joint, y_joint)
print(f"Pearson correlation: r={corr:.6f} p={corr_p:.4f}")

# Spearman
scorr, scorr_p = stats.spearmanr(x_joint, y_joint)
print(f"Spearman correlation: rho={scorr:.6f} p={scorr_p:.4f}")

# ── Write results ──────────────────────────────────────────────────────
print(f"\nWriting results to {RESULTS}")

with open(str(RESULTS), 'w') as f:
    f.write("# Zeta Zeros mod Various Constants: Pattern Search\n\n")
    f.write(f"**Date**: 2026-04-05\n")
    f.write(f"**N**: {N} zeros, range [{gammas[0]:.4f}, {gammas[-1]:.4f}]\n\n")

    f.write("## 1. Per-Modulus Uniformity Tests\n\n")
    f.write("| Modulus | Value | KS stat | KS p-val | Circ Mean | Circ Var | R_bar | Rayleigh Z | Rayleigh p |\n")
    f.write("|---------|-------|---------|----------|-----------|----------|-------|------------|------------|\n")
    for name, r in results.items():
        f.write(f"| {name} | {r['m']:.6f} | {r['ks_stat']:.4f} | {r['ks_p']:.4f} | "
                f"{r['circ_mean']:.4f} | {r['circ_var']:.4f} | {r['R_bar']:.4f} | "
                f"{r['rayleigh_Z']:.4f} | {r['rayleigh_p']:.4f} |\n")

    # Flag any significant results
    sig_ks = [n for n, r in results.items() if r['ks_p'] < 0.05]
    sig_ray = [n for n, r in results.items() if r['rayleigh_p'] < 0.05]

    f.write("\n### Significance Summary\n\n")
    if sig_ks:
        f.write(f"- **KS test rejects uniformity (p<0.05)** for: {', '.join(sig_ks)}\n")
    else:
        f.write("- KS test: NO modulus shows significant deviation from uniformity\n")
    if sig_ray:
        f.write(f"- **Rayleigh test rejects circular uniformity (p<0.05)** for: {', '.join(sig_ray)}\n")
    else:
        f.write("- Rayleigh test: NO modulus shows significant circular non-uniformity\n")

    f.write("\n## 2. Mod 1 Detailed Analysis\n\n")
    f.write("### Weyl Criterion (Equidistribution)\n\n")
    f.write(f"- Checked |S_k| = |(1/N) sum exp(2 pi i k gamma_n)| for k=1..50\n")
    f.write(f"- Maximum |S_k| = {max_weyl[1]:.6f} at k={max_weyl[0]}\n")
    f.write(f"- Mean |S_k| = {np.mean([w[1] for w in weyl_results]):.6f}\n")
    f.write(f"- Expected for random: 1/sqrt(N) = {1/np.sqrt(N):.6f}\n")
    weyl_above = sum(1 for w in weyl_results if w[1] > 2/np.sqrt(N))
    f.write(f"- Weyl sums above 2/sqrt(N): {weyl_above}/50\n")
    if max_weyl[1] < 3/np.sqrt(N):
        f.write(f"- **Conclusion**: Sequence appears equidistributed mod 1\n")
    else:
        f.write(f"- **Conclusion**: Some Weyl sums are elevated -- possible non-equidistribution\n")

    f.write("\n### Discrepancy\n\n")
    f.write(f"- D_N = {D_N:.6f}\n")
    f.write(f"- Expected for random uniform: ~ {expected_D:.6f}\n")
    f.write(f"- Ratio D_N / expected: {D_N / expected_D:.2f}\n")

    f.write("\n### Gap Analysis (mod 1)\n\n")
    f.write(f"- Mean gap: {np.mean(gaps_full):.6f} (expected {mean_gap:.6f})\n")
    f.write(f"- Min gap: {np.min(gaps_full):.6f} (ratio to mean: {np.min(gaps_full)/mean_gap:.4f})\n")
    f.write(f"- Max gap: {np.max(gaps_full):.6f} (ratio to mean: {np.max(gaps_full)/mean_gap:.2f})\n")
    f.write(f"- Normalized gaps vs Exp(1): KS={gap_ks:.4f}, p={gap_ks_p:.4f}\n")
    if gap_ks_p < 0.05:
        f.write(f"- **Gap distribution deviates from exponential** (p={gap_ks_p:.4f})\n")
        f.write(f"  This suggests the mod-1 fractional parts are NOT like independent uniform samples.\n")
        f.write(f"  Consistent with GUE spacing statistics of zeta zeros.\n")
    else:
        f.write(f"- Gap distribution consistent with exponential (independent uniform)\n")

    f.write("\n## 3. Joint Distribution (mod 1 vs mod pi)\n\n")
    f.write(f"- Chi-squared independence: chi2={chi2:.2f}, dof={dof}, p={chi2_p:.4f}\n")
    f.write(f"- Pearson r = {corr:.6f}, p = {corr_p:.4f}\n")
    f.write(f"- Spearman rho = {scorr:.6f}, p = {scorr_p:.4f}\n")
    if chi2_p > 0.05 and corr_p > 0.05:
        f.write("- **Conclusion**: No evidence of dependence between mod 1 and mod pi\n")
    else:
        f.write("- **Conclusion**: Some evidence of dependence detected\n")

    f.write("\n## 4. Overall Conclusions\n\n")

    any_pattern = False
    notes = []
    if sig_ks:
        any_pattern = True
        notes.append(f"KS non-uniformity for mod {', '.join(sig_ks)}")
    if sig_ray:
        any_pattern = True
        notes.append(f"Rayleigh non-uniformity for mod {', '.join(sig_ray)}")
    if gap_ks_p < 0.05:
        any_pattern = True
        notes.append("Mod-1 gap distribution deviates from exponential (expected from GUE)")
    if max_weyl[1] > 3/np.sqrt(N):
        any_pattern = True
        notes.append(f"Elevated Weyl sum at k={max_weyl[0]}")
    if chi2_p < 0.05 or corr_p < 0.05:
        any_pattern = True
        notes.append("Joint (mod 1, mod pi) shows dependence")

    if any_pattern:
        f.write("Patterns detected:\n\n")
        for note in notes:
            f.write(f"- {note}\n")
        f.write("\nHowever, none of these patterns provide a route to computing p(n) in polylog time.\n")
        f.write("The gap deviation from exponential is expected: zeta zeros exhibit GUE repulsion,\n")
        f.write("which creates short-range correlations in their fractional parts.\n")
    else:
        f.write("**No significant patterns found.** The zeta zeros appear uniformly distributed\n")
        f.write("modulo all tested constants, consistent with theoretical expectations.\n")
        f.write("This is further evidence that the oscillatory part of pi(x) is pseudorandom.\n")

print("\nDone. All plots and results saved.")
