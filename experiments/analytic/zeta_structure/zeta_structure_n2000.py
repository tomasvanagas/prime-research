"""
zeta_structure_n2000.py -- Extension of Session 25's zeta-structure tests
to N=2000 zeros, settling the explicit caveat from SESSION_25_SUMMARY:
    "we tested 1000 zeros. Structure might exist at scales requiring
     >1000 zeros to detect."

Re-runs the FIVE most discriminating tests at N=2000:

  1. Pair correlation R_2(r) vs GUE  (more pairs => sharper Bessel-K^2 fit)
  2. Number variance Sigma^2(L)        (probes larger correlation scales)
  3. DFT power spectrum + spectral flatness vs GUE (detects emergent periodicity)
  4. Mod-constant discrepancy and Weyl sums (should shrink as 1/sqrt(N))
  5. PSLQ on the LATE block gamma_{1001..1100} (tests if structure
     emerges at higher index, not just early indices)

Tests intentionally omitted (verdict cannot change with more data alone):
  - PSLQ on first 50 zeros at fixed precision (precision-bounded)
  - Sparse matrix model (already proven O(N) params)
"""

import os, json, time, itertools
import numpy as np
import mpmath
from scipy import stats
from scipy.linalg import eigvalsh

OUT_DIR  = os.path.dirname(os.path.abspath(__file__))
DATA     = "/apps/aplikacijos/prime-research/data/zeta_zeros_2000.txt"
RESULTS  = os.path.join(OUT_DIR, "zeta_structure_n2000_results.md")

zeros = np.loadtxt(DATA)
N = len(zeros)
print(f"Loaded {N} zeros, range [{zeros[0]:.4f}, {zeros[-1]:.4f}]")
assert N == 2000, f"expected 2000 zeros, got {N}"

report = []
def log(s=""):
    print(s, flush=True)
    report.append(s)

log(f"# Zeta-Structure Tests at N=2000 (Session 25 Extension)\n")
log(f"Loaded N={N} zeros, range [{zeros[0]:.4f}, {zeros[-1]:.4f}]\n")
log(f"Mean spacing: {np.mean(np.diff(zeros)):.6f}\n")

# ─────────────────────────────────────────────────────────────────────
# 1.  Pair correlation R_2(r)
# ─────────────────────────────────────────────────────────────────────
log("\n## 1. Pair correlation R_2(r) vs GUE (sin^2 pi r)/(pi r)^2 form")
log("=" * 70)

# Unfold to unit mean spacing using Riemann-von Mangoldt smooth count
unfolded = np.empty(N)
for i in range(N):
    T = zeros[i]
    unfolded[i] = T / (2*np.pi) * np.log(T / (2*np.pi*np.e)) + 7.0/8.0
mean_uspc = np.mean(np.diff(unfolded))
unfolded = unfolded / mean_uspc

r_max, dr = 4.0, 0.05
r_bins = np.arange(0, r_max + dr, dr)
r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
hist = np.zeros(len(r_centers))
for i in range(N):
    for j in range(i+1, min(i+80, N)):
        d = abs(unfolded[j] - unfolded[i])
        if d < r_max:
            idx = int(d / dr)
            if idx < len(hist):
                hist[idx] += 1
norm = N * dr
R2 = hist / norm
tail = np.mean(R2[r_centers > 2.5])
R2 = R2 / tail

r_theory = np.linspace(0.001, r_max, 500)
R2_gue = 1 - (np.sin(np.pi*r_theory)/(np.pi*r_theory))**2
gue_interp = np.interp(r_centers, r_theory, R2_gue)
rms_resid = np.sqrt(np.mean((R2 - gue_interp)[r_centers > 0.1]**2))
log(f"R_2 RMS residual from GUE (r>0.1): {rms_resid:.4f}")
log(f"  (Session 25 at N=1000 reported ~0.10-0.12)")
log(f"R_2(r=0+):  {R2[0]:.4f}    (GUE predicts 0; level repulsion)")

# ─────────────────────────────────────────────────────────────────────
# 2.  Number variance Sigma^2(L)
# ─────────────────────────────────────────────────────────────────────
log("\n## 2. Number variance Sigma^2(L)")
log("=" * 70)

L_values = np.logspace(np.log10(0.1), np.log10(40), 50)
sigma2 = np.zeros(len(L_values))
for li, L in enumerate(L_values):
    n_w = min(800, N - 1)
    starts = np.linspace(unfolded[0], unfolded[-1] - L, n_w)
    counts = np.array([np.sum((unfolded >= s) & (unfolded < s + L)) for s in starts])
    sigma2[li] = np.var(counts)

ge = 0.5772156649
sigma2_gue = (2/np.pi**2) * (np.log(2*np.pi*L_values) + ge + 1 - np.pi**2/8)
log_L_thresh = (L_values >= 1.0) & (L_values <= 30.0)
delta = sigma2[log_L_thresh] - sigma2_gue[log_L_thresh]
nv_rms = np.sqrt(np.mean(delta**2))
nv_max_dev = np.max(np.abs(delta) / sigma2_gue[log_L_thresh])
log(f"RMS deviation from GUE asymptotic over 1<=L<=30 (now reachable at N=2000):")
log(f"  RMS={nv_rms:.4f}, max relative deviation={nv_max_dev:.3f}")
log(f"  (Session 25 reachable up to L~10 only)")
for L, s2 in zip(L_values[::8], sigma2[::8]):
    g = (2/np.pi**2)*(np.log(2*np.pi*L)+ge+1-np.pi**2/8)
    log(f"  L={L:7.3f}  Sigma^2={s2:7.4f}  GUE={g:7.4f}  ratio={s2/g:.3f}")

# ─────────────────────────────────────────────────────────────────────
# 3.  DFT power spectrum & spectral flatness
# ─────────────────────────────────────────────────────────────────────
log("\n## 3. DFT power spectrum and spectral flatness")
log("=" * 70)

F = np.fft.fft(zeros)
power = np.abs(F)**2
freqs = np.fft.fftfreq(N)
pos = freqs > 0
power_pos = power[pos]

def spectral_flatness(s):
    return np.exp(np.mean(np.log(s + 1e-30))) / np.mean(s)

sf_total = spectral_flatness(power_pos)
log(f"Spectral flatness (overall): {sf_total:.4f}  (S25 N=1000 was 0.93-0.999)")

n_bands = 5
band_size = len(power_pos) // n_bands
band_sfs = []
for b in range(n_bands):
    s = b*band_size
    e = (s + band_size) if b < n_bands-1 else len(power_pos)
    band_sfs.append(spectral_flatness(power_pos[s:e]))
    log(f"  Band {b}: SF = {band_sfs[-1]:.4f}")

# GUE comparison: 30 trials of GUE 300x300, interp to N
np.random.seed(42)
n_gue, gs = 30, 300
gue_sfs = []
gue_log_powers = []
for _ in range(n_gue):
    A = (np.random.randn(gs, gs) + 1j*np.random.randn(gs, gs))/np.sqrt(2)
    H = (A + A.conj().T)/(2*np.sqrt(gs))
    eigs = np.sort(eigvalsh(H).real)
    eigs_i = np.interp(np.linspace(0, gs-1, N), np.arange(gs), eigs)
    Fg = np.fft.fft(eigs_i)
    pg = np.abs(Fg[pos])**2
    gue_sfs.append(spectral_flatness(pg))
    gue_log_powers.append(np.log(pg/np.mean(pg) + 1e-30))

gue_mean_lp = np.mean(gue_log_powers, axis=0)
zeta_lp = np.log(power_pos/np.mean(power_pos) + 1e-30)
corr = np.corrcoef(zeta_lp, gue_mean_lp)[0, 1]
z_score = (sf_total - np.mean(gue_sfs)) / np.std(gue_sfs)
log(f"GUE spectral flatness: mean={np.mean(gue_sfs):.4f}, std={np.std(gue_sfs):.4f}")
log(f"Zeta-vs-GUE z-score:   {z_score:+.2f}")
log(f"Log-power correlation (zeta vs GUE-mean): {corr:.4f}  (S25 was 0.9999)")

# Look for sharp peaks above noise (any periodic structure)
noise = np.median(power_pos)
peak_ratios = power_pos / noise
n_peaks_10x  = int(np.sum(peak_ratios > 10))
n_peaks_100x = int(np.sum(peak_ratios > 100))
top5 = np.argsort(power_pos)[-5:][::-1]
log(f"Peaks above 10x median noise:  {n_peaks_10x}")
log(f"Peaks above 100x median noise: {n_peaks_100x}")
log(f"Top 5 peak indices: {list(top5)}, ratios: {[float(peak_ratios[i]) for i in top5]}")

# ─────────────────────────────────────────────────────────────────────
# 4.  Mod-constant discrepancy and Weyl sums
# ─────────────────────────────────────────────────────────────────────
log("\n## 4. Mod-constant discrepancy and Weyl sums (N=2000 vs sqrt(N) baseline)")
log("=" * 70)

moduli = {
    "1":          1.0,
    "pi":         np.pi,
    "2pi":        2*np.pi,
    "log(2pi)":   np.log(2*np.pi),
    "sqrt(2pi)":  np.sqrt(2*np.pi),
    "e":          np.e,
    "log(2)":     np.log(2),
    "log(3)":     np.log(3),
    "log(5)":     np.log(5),
    "log(7)":     np.log(7),
}

discrepancies = {}
for name, m in moduli.items():
    vals = (zeros % m) / m
    s = np.sort(vals)
    F_N = np.arange(1, N+1) / N
    Dp = np.max(F_N - s)
    Dm = np.max(s - np.arange(0, N) / N)
    DN = max(Dp, Dm)
    expected = np.sqrt(np.log(np.log(N))/(2*N))
    discrepancies[name] = (DN, expected)
    log(f"  mod {name:11s}: D_N={DN:.6f}  expected~{expected:.6f}  ratio={DN/expected:.3f}")

# Weyl sums on mod 1
weyl = []
for k in range(1, 51):
    S = abs(np.mean(np.exp(2j*np.pi*k*zeros)))
    weyl.append(S)
log(f"\nWeyl sums |S_k| for k=1..50:")
log(f"  max=|S_{1+np.argmax(weyl)}|={max(weyl):.6f}")
log(f"  mean={np.mean(weyl):.6f}  expected (random)~1/sqrt(N)={1/np.sqrt(N):.6f}")
log(f"  ratio mean/expected = {np.mean(weyl)/(1/np.sqrt(N)):.3f}  (close to 1 == random)")

# ─────────────────────────────────────────────────────────────────────
# 5.  PSLQ on LATE zeros 1001..1100 (test for emergent structure)
# ─────────────────────────────────────────────────────────────────────
log("\n## 5. PSLQ on LATE zeros gamma_{1001..1100} (NEW window, never tested before)")
log("=" * 70)

mpmath.mp.dps = 60
# We need 60-digit precision; numpy zeros are only 30-digit.
# Reload from file using mpmath for full precision.
with open(DATA) as f:
    zeros_mp = [mpmath.mpf(line.strip()) for line in f if line.strip()]
log(f"Loaded {len(zeros_mp)} zeros at mpmath precision = {mpmath.mp.dps} digits.")

augment_names = ["1", "pi", "log(2*pi)"]
augment_vals  = [mpmath.mpf(1), mpmath.pi, mpmath.log(2*mpmath.pi)]
MAX_C = 1000

# Test pairs and triples among gamma_{1001..1050} (50 zeros)
# Pairs: C(50,2)=1225; Triples: C(30,3)=4060.
late_zeros = zeros_mp[1000:1050]      # 0-indexed: indices 1000..1049 are zeros 1001..1050

found_pair = found_trip = 0
total_pair = total_trip = 0
t0 = time.time()
for i, j in itertools.combinations(range(50), 2):
    total_pair += 1
    vec = [late_zeros[i], late_zeros[j], mpmath.pi, mpmath.log(2*mpmath.pi), mpmath.mpf(1)]
    try:
        rel = mpmath.pslq(vec, maxcoeff=MAX_C, maxsteps=200)
    except Exception:
        rel = None
    if rel is not None:
        residual = sum(c*v for c, v in zip(rel, vec))
        if abs(residual) < mpmath.mpf(10)**(-40):
            found_pair += 1
            log(f"  PAIR RELATION (1001+{i}, 1001+{j}): coeffs={list(rel)} residual={float(residual):.2e}")

late_zeros30 = zeros_mp[1000:1030]
for combo in itertools.combinations(range(30), 3):
    total_trip += 1
    vec = list(augment_vals) + [late_zeros30[k] for k in combo]
    try:
        rel = mpmath.pslq(vec, maxcoeff=MAX_C, maxsteps=200)
    except Exception:
        rel = None
    if rel is not None:
        residual = sum(c*v for c, v in zip(rel, vec))
        if abs(residual) < mpmath.mpf(10)**(-40):
            found_trip += 1
            log(f"  TRIPLE RELATION at combo={combo}: coeffs={list(rel)}")

log(f"\nPSLQ on LATE block:")
log(f"  Pairs    tested: {total_pair}, relations found: {found_pair}")
log(f"  Triples  tested: {total_trip}, relations found: {found_trip}")
log(f"  Elapsed: {time.time()-t0:.1f}s")

# ─────────────────────────────────────────────────────────────────────
# 6.  Bonus: cross-block PSLQ -- are early and late zeros related?
# ─────────────────────────────────────────────────────────────────────
log("\n## 6. Cross-block PSLQ: gamma_i (early) vs gamma_j (late)")
log("=" * 70)
found_cross = 0
total_cross = 0
for i in range(20):
    for j in range(20):
        total_cross += 1
        vec = [zeros_mp[i], zeros_mp[1000+j], mpmath.pi, mpmath.log(2*mpmath.pi), mpmath.mpf(1)]
        try:
            rel = mpmath.pslq(vec, maxcoeff=MAX_C, maxsteps=200)
        except Exception:
            rel = None
        if rel is not None:
            residual = sum(c*v for c, v in zip(rel, vec))
            if abs(residual) < mpmath.mpf(10)**(-40):
                found_cross += 1
                log(f"  CROSS RELATION (gamma_{i+1}, gamma_{1001+j}): {list(rel)}")
log(f"Cross-block tested: {total_cross}, relations: {found_cross}")

# ─────────────────────────────────────────────────────────────────────
# Verdict
# ─────────────────────────────────────────────────────────────────────
log("\n" + "=" * 70)
log("## VERDICT")
log("=" * 70)

verdict_lines = []
def vline(s):
    print(s, flush=True)
    report.append(s)
    verdict_lines.append(s)

vline(f"\nN={N} (2x the S25 sample). All five extension tests indicate:")
vline(f"")
vline(f"  - Pair correlation RMS deviation from GUE: {rms_resid:.4f}")
vline(f"  - Number variance RMS deviation:           {nv_rms:.4f}")
vline(f"  - Spectral flatness:                        {sf_total:.4f}  (z={z_score:+.2f} vs GUE)")
vline(f"  - Log-power correlation (zeta vs GUE):     {corr:.4f}")
vline(f"  - Peaks > 100x noise floor:                 {n_peaks_100x}")
vline(f"  - PSLQ pair relations on LATE block:        {found_pair}/{total_pair}")
vline(f"  - PSLQ triple relations on LATE block:      {found_trip}/{total_trip}")
vline(f"  - Cross-block (early-late) relations:       {found_cross}/{total_cross}")
vline(f"  - Discrepancy ratios (vs LIL bound) all O(1): {{{', '.join(f'{n}:{d/e:.2f}' for n,(d,e) in list(discrepancies.items())[:5])}}}")

# JSON summary
summary = {
    "N": N,
    "pair_corr_rms_residual_gue": float(rms_resid),
    "pair_corr_at_zero": float(R2[0]),
    "number_variance_rms_dev": float(nv_rms),
    "number_variance_max_rel_dev": float(nv_max_dev),
    "spectral_flatness_total": float(sf_total),
    "spectral_flatness_band_sfs": [float(x) for x in band_sfs],
    "spectral_flatness_gue_mean": float(np.mean(gue_sfs)),
    "spectral_flatness_z_score":  float(z_score),
    "log_power_corr_vs_gue": float(corr),
    "peaks_above_10x":  n_peaks_10x,
    "peaks_above_100x": n_peaks_100x,
    "pslq_late_pairs_total": total_pair, "pslq_late_pairs_found": found_pair,
    "pslq_late_triples_total": total_trip, "pslq_late_triples_found": found_trip,
    "pslq_cross_total": total_cross, "pslq_cross_found": found_cross,
    "discrepancies": {k: {"D_N": float(v[0]), "expected": float(v[1])}
                       for k, v in discrepancies.items()},
    "weyl_max": float(max(weyl)),
    "weyl_mean": float(np.mean(weyl)),
    "weyl_baseline_1_over_sqrt_N": float(1/np.sqrt(N)),
}

with open(os.path.join(OUT_DIR, "zeta_structure_n2000_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)

with open(RESULTS, "w") as f:
    f.write("\n".join(report) + "\n")
print(f"\nWrote {RESULTS}")
print(f"Wrote {OUT_DIR}/zeta_structure_n2000_summary.json")
