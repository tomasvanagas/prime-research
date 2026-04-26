"""
FOCUS-2 (TODO.md) E2.11 pre-test: iterated finite-difference signature of
the six concrete candidate intermediate quantities for a fourth encoding
of pi(x).

Methodology
-----------
For pi(x) - R(x), iterated finite differences Delta^k grow with RMS ratio
RMS(Delta^{k+1})/RMS(Delta^k) -> 2.0 (the canonical white-noise signature;
see E2.11 in EDGES.md / experiments/algebraic/identity_search/wz_definite_sum.py).
A candidate T_i(x) whose detrended residual exhibits the same signature is
just another encoding of the SAME GUE-random information as pi(x), and we
close it as mode I (info loss / equivalence) without running PSLQ.

Candidates (FOCUS-2):
  1. T_1(x) = sum_{n<=x} {log Gamma(n)}     (fractional part)
  2. T_2(x) = sum_{n<=x} H_n
     T_3(x) = sum_{n<=x} H_n^2
  3. Psi(x, B) for B = log^c(x), c in {2, 3, 4}
  4. sigma_k summatories for k = 2, 3 (k=1 covered in S64)
  5. Q(x) squarefree count
  6. phi(n) summatory       (covered S64; re-checked with E2.11)

Three controls bracket the verdict scale:
  C1 - white-noise random: ratio -> 2.0
  C2 - smooth polynomial 0.001 x^2 + x: annihilated by Delta^3
  C3 - pi(x) - R(x): expected ratio -> 2.0 (literal E2.11 reference)

Output
------
Per-candidate row: residual std, RMS(Delta^k), RMS ratios at k = 1..7.
Verdict { WHITE | SMOOTH | INTERMEDIATE } from the ratio asymptote.
"""
import numpy as np
import scipy.special as sp
import time
from mpmath import mp, mpf, mpc, exp as mp_exp, log as mp_log, ei as mp_ei, pi as mp_pi
mp.dps = 40

# ---------- Parameters ----------
N = 200_000          # sieve up to here
PROBE_LO = 100_000   # contiguous probe window
PROBE_HI = 150_000   # 50000 consecutive integers -> robust LSQ + Delta^7
MAX_K = 7            # iterate Delta up to this order
DETREND_DEG = 6      # polynomial detrend degree (absorbs x, x log x, x^2, ...)
print(f"# E2.11 pre-test on FOCUS-2 candidates")
print(f"# sieve N={N}, probe [{PROBE_LO}, {PROBE_HI}), max Delta^{MAX_K}, detrend deg {DETREND_DEG}")

x_probe = np.arange(PROBE_LO, PROBE_HI, dtype=np.int64)
W = len(x_probe)

# ---------- Sieves ----------
print(f"\n## Sieving multiplicative functions to N={N} ...")
t0 = time.time()

# smallest prime factor
spf = np.zeros(N+1, dtype=np.int64)
for p in range(2, N+1):
    if spf[p] == 0:
        for k in range(p, N+1, p):
            if spf[k] == 0:
                spf[k] = p

is_prime = (spf == np.arange(N+1))
is_prime[0] = is_prime[1] = False

# Largest prime factor (for smoothness indicator)
lpf = np.zeros(N+1, dtype=np.int64)
for p in range(2, N+1):
    if is_prime[p]:
        for k in range(p, N+1, p):
            lpf[k] = p

# Multiplicative functions via spf factorisation pass
mu = np.zeros(N+1, dtype=np.int64); mu[1] = 1
phi = np.zeros(N+1, dtype=np.int64); phi[1] = 1
Q = np.zeros(N+1, dtype=np.int64); Q[1] = 1   # squarefree indicator (1 if mu^2 = 1)
for n in range(2, N+1):
    p = spf[n]
    m = n // p
    if m % p == 0:
        # p^2 | n: squarefull factor
        mu[n] = 0
        Q[n] = 0
        phi[n] = p * phi[m]
    else:
        mu[n] = -mu[m]
        Q[n] = Q[m]                       # squarefree iff m squarefree
        phi[n] = (p - 1) * phi[m]

# sigma_k via direct double-loop (k=2,3 — k=1 done in S64)
sigma_2 = np.zeros(N+1, dtype=np.float64)
sigma_3 = np.zeros(N+1, dtype=np.float64)
for d in range(1, N+1):
    d2 = d*d
    d3 = d2*d
    for n in range(d, N+1, d):
        sigma_2[n] += d2
        sigma_3[n] += d3

print(f"  sieves done in {time.time()-t0:.2f}s")

# ---------- Helpers ----------
def difference_signature(y, max_k=MAX_K):
    """Return (rms_list, ratio_list) for Delta^1..Delta^max_k."""
    delta = np.asarray(y, dtype=np.float64).copy()
    rms = []
    for _ in range(max_k):
        delta = np.diff(delta)
        rms.append(float(np.sqrt(np.mean(delta**2))))
    ratios = [rms[k]/rms[k-1] if rms[k-1] > 0 else float('inf') for k in range(1, len(rms))]
    return rms, ratios

def detrend(y, deg=DETREND_DEG):
    """Polynomial detrend on normalised x. Robust to leading x, x log x, x^2 terms."""
    x = np.arange(len(y), dtype=np.float64)
    xn = (x - x.mean()) / x.std()
    coefs = np.polyfit(xn, y, deg)
    return y - np.polyval(coefs, xn)

def verdict(ratios):
    """Classify based on Delta^k RMS ratio asymptote.

    WHITE       : last ratios ~ 2.0 (within 0.15)
    SMOOTH      : ratios diverge or collapse (annihilated)
    INTERMEDIATE: bounded but not 2.0 (e.g. ~ sqrt(2), 1.5)
    """
    last = ratios[-3:]
    avg = sum(last) / len(last)
    spread = max(last) - min(last)
    if any(r > 5 or r < 0.5 for r in last):
        return f"SMOOTH (avg={avg:.3f}, spread={spread:.3f})"
    if abs(avg - 2.0) < 0.15:
        return f"WHITE (avg={avg:.3f})"
    return f"INTERMEDIATE (avg={avg:.3f})"

def report(name, y_full, leading_form_note=""):
    """Take y on x_probe, detrend, run E2.11 signature, print row."""
    y_det = detrend(y_full)
    rms, ratios = difference_signature(y_det)
    print(f"\n### {name}")
    if leading_form_note:
        print(f"  leading form: {leading_form_note}")
    print(f"  std(residual)   = {y_det.std():.6e}")
    print(f"  RMS(Delta^k):   " + ", ".join(f"k={k+1}:{r:.3e}" for k, r in enumerate(rms)))
    print(f"  Ratios:         " + ", ".join(f"{r:.4f}" for r in ratios))
    print(f"  VERDICT: {verdict(ratios)}")
    return rms, ratios

# ---------- Controls ----------
print("\n\n# === CONTROLS ===")

rng = np.random.default_rng(2026)
report("C1 white noise (random N(0,1) on probe range)",
       rng.standard_normal(W),
       "i.i.d. N(0,1) — expected ratio ~ 2.0 asymptotically")

# smooth polynomial: detrend kills it -> Delta^k annihilated
xp = x_probe.astype(np.float64)
report("C2 smooth polynomial 0.001 x^2 + x",
       0.001 * xp**2 + xp,
       "detrend deg 6 absorbs entire signal; Delta^k of residual ~ 0")

# pi(x) - R(x): the canonical E2.11 reference
print("\n## C3 pi(x) - R(x) ...")
t1 = time.time()
pi_arr = np.cumsum(is_prime[1:N+1].astype(np.int64))   # pi[i] = pi(i+1) (offset by 1)
# pi_arr[k] = pi(k+1) for k=0..N-1
# we want pi(x) for x in x_probe
pi_probe = pi_arr[x_probe - 1].astype(np.float64)

# R(x) via mpmath: R(x) = sum_k mu(k)/k * li(x^{1/k})
def R_mp(x, K=20):
    """Riemann's R function. K=20 terms suffices for x up to 10^9 to high precision."""
    xv = mpf(x)
    s = mpf(0)
    for k in range(1, K+1):
        mk = int(mu[k]) if k < len(mu) else 0
        if mk == 0:
            continue
        # li(x^{1/k}) = ei(log(x)/k)
        s += mpf(mk) / k * mp_ei(mp_log(xv) / k)
    return s

# Sample R only at endpoints + a coarse grid; smooth so polyfit ok
print(f"  computing R at 200 grid points then interpolating ...")
sample_idx = np.linspace(0, W-1, 200).astype(int)
R_samples = np.array([float(R_mp(int(x_probe[i]))) for i in sample_idx])
# linear-spline interpolation onto full grid; R is C^infinity so this is fine
from numpy import interp
R_probe = interp(np.arange(W).astype(np.float64), sample_idx.astype(np.float64), R_samples)

f_x = pi_probe - R_probe
report("C3 pi(x) - R(x) (E2.11 reference)",
       f_x,
       "should produce ratio -> 2.0 per E2.11")
print(f"  (R computation: {time.time()-t1:.1f}s)")

# ---------- Candidates ----------
print("\n\n# === FOCUS-2 CANDIDATES ===")

# 1. T_1(x) = sum_{n<=x} {log Gamma(n)}
print("\n## Candidate 1: T_1(x) = sum_{n<=x} {log Gamma(n)} (fractional part)")
t2 = time.time()
ln_gamma = sp.gammaln(np.arange(1, N+1)).astype(np.float64)
frac_ln_gamma = ln_gamma - np.floor(ln_gamma)        # fractional part of log Gamma
T1_full = np.cumsum(frac_ln_gamma)
T1_probe = T1_full[x_probe - 1]
report("T_1 = sum {log Gamma(n)}", T1_probe,
       "leading ~ x/2 (mean of {log Gamma} -> 1/2 by equidistribution mod 1)")

# 2a. T_2(x) = sum_{n<=x} H_n
print("\n## Candidate 2a: T_2(x) = sum_{n<=x} H_n (harmonic-number summatory)")
inv_n = 1.0 / np.arange(1, N+1, dtype=np.float64)
H = np.cumsum(inv_n)            # H[i] = H_{i+1}
T2_full = np.cumsum(H)          # T2[i] = sum_{n=1..i+1} H_n
T2_probe = T2_full[x_probe - 1]
report("T_2 = sum H_n", T2_probe,
       "closed form: (x+1) H_x - x; leading x log x")

# 2b. T_3(x) = sum_{n<=x} H_n^2
print("\n## Candidate 2b: T_3(x) = sum_{n<=x} H_n^2")
T3_full = np.cumsum(H**2)
T3_probe = T3_full[x_probe - 1]
report("T_3 = sum H_n^2", T3_probe,
       "leading x (log x)^2; bounded oscillation only from gamma constant")

# 3. Psi(x, B) for B = log^c(x), c in {2, 3, 4}
print("\n## Candidate 3: Psi(x, log^c x) for c in {2, 3, 4} (smooth-number count, varying B)")
# Smoothness: lpf(n) <= B
# We'll evaluate for *multiple* B values across the probe range (B varies with x).
# Concretely: run with B = log^c(midpoint) -- a fixed B for the whole window, then
# run several windows. For E2.11 it suffices to evaluate one window with a
# fixed B large enough to admit non-trivial dependence on x.
def psi_xB_summatory(B, n_full=N):
    """Sum_{n<=N} 1[lpf(n) <= B]; n=1 contributes 1."""
    smooth_indic = np.zeros(n_full + 1, dtype=np.int64)
    smooth_indic[1] = 1
    smooth_indic[2:n_full+1] = (lpf[2:n_full+1] <= B).astype(np.int64)
    return np.cumsum(smooth_indic[1:n_full+1])

x_mid = (PROBE_LO + PROBE_HI) // 2
log_x_mid = np.log(x_mid)
for c in (2, 3, 4):
    B = max(2, int(log_x_mid ** c))
    if B > N:
        B = N
    print(f"\n  c = {c}: B = floor(log(x_mid)^{c}) = {B}")
    Psi_full = psi_xB_summatory(B)
    Psi_probe = Psi_full[x_probe - 1].astype(np.float64)
    report(f"Psi(x, B={B})", Psi_probe,
           f"smooth-number count with B = log^{c}(x); main term Dickman-rho")

# 4. sigma_k summatory for k = 2, 3
print("\n## Candidate 4: sigma_k summatory for k in {2, 3}")
S2_full = np.cumsum(sigma_2)        # offset: S2[i] = sum_{n=1..i} sigma_2(n) (since sigma_2[0]=0)
S2_full = S2_full[1:]               # drop the 0 entry; now S2[i] = sum_{n=1..i+1} but sigma_2[0] is 0 so it's just cumsum
# simpler:
S2_correct = np.cumsum(sigma_2[1:])  # sum_{n=1..k+1} sigma_2(n) for k=0..N-1, length N
S3_correct = np.cumsum(sigma_3[1:])
S2_probe = S2_correct[x_probe - 1].astype(np.float64)
S3_probe = S3_correct[x_probe - 1].astype(np.float64)
report("sigma_2 summatory", S2_probe, "leading zeta(3)/3 * x^3")
report("sigma_3 summatory", S3_probe, "leading zeta(4)/4 * x^4")

# 5. Q(x) squarefree count
print("\n## Candidate 5: Q(x) = #{n<=x : mu(n)^2 = 1}")
Q_full = np.cumsum(Q[1:N+1])
Q_probe = Q_full[x_probe - 1].astype(np.float64)
report("Q(x) squarefree count", Q_probe,
       "leading 6/pi^2 * x; residual is Mobius/Mertens-type")

# 6. phi summatory (also tested S64 by ρ; here by E2.11)
print("\n## Candidate 6: T_6(x) = sum_{n<=x} phi(n)")
T6_full = np.cumsum(phi[1:N+1])
T6_probe = T6_full[x_probe - 1].astype(np.float64)
report("T_6 = sum phi(n)", T6_probe,
       "leading 3/pi^2 * x^2; residual Mertens-type")

print(f"\n\n# Total experiment time: {time.time() - t0:.1f}s")
print("# Summary: all WHITE => 4th-encoding pre-test confirms mode I closure")
print("#          (residuals are GUE-random like pi(x) - R(x), no structural distinction).")
