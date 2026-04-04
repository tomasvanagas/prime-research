#!/usr/bin/env python3
"""
SESSION 4: HIERARCHICAL DECOMPOSITION OF THE PRIME CORRECTION
==============================================================

Goal: decompose delta(n) = p(n) - R^{-1}(n) into computable layers.

Four approaches:
  A) Wavelet decomposition of delta(n) - multi-scale structure
  B) Cumulative correction Sigma(n) = sum_{k=1}^{n} delta(k) - smoother integral
  C) Chebyshev bias exploitation - residue class constraints
  D) Multi-scale R^{-1} with zeta zeros - oscillatory correction from rho
"""

import math
import time
import numpy as np
from scipy import signal, stats, optimize
from scipy.fft import fft, fftfreq
import pywt
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# CORE: Riemann R function and inverse (float precision)
# ============================================================

_EULER_GAMMA = 0.5772156649015328606065120900824024
_MU = [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1, 0,
       -1, 0, -1, 0, 1, 1, -1, 0, 0, 1, 0, 0, -1, -1, -1, 0, 1,
       1, 1, 0, -1, 1, 1, 0, -1, -1, -1, 0, 0, -1, -1, 0, 0, 0]

def _li(x):
    """Logarithmic integral via Ramanujan series."""
    if x <= 1.0: return 0.0
    lnx = math.log(x)
    r = _EULER_GAMMA + math.log(abs(lnx))
    t = 1.0
    for k in range(1, 200):
        t *= lnx / k
        c = t / k
        r += c
        if abs(c) < 1e-15 * max(1.0, abs(r)):
            break
    return r

def R_func(x):
    """Riemann R(x) = sum mu(k)/k * li(x^{1/k})."""
    if x <= 1: return 0.0
    r = 0.0
    for k in range(1, len(_MU)):
        if _MU[k] == 0: continue
        xk = x ** (1.0 / k)
        if xk <= 1.0001: break
        r += _MU[k] / k * _li(xk)
    return r

def R_prime(x):
    """R'(x) = 1/log(x) * sum mu(k)/k^2 * x^{1/k - 1}... approx 1/log(x)."""
    if x <= 2: return 1.0
    return 1.0 / math.log(x)

def inv_R(n):
    """R^{-1}(n): Newton iteration."""
    if n <= 5: return [0, 2, 3, 5, 7, 11][n]
    x = float(n) * math.log(n) + float(n) * math.log(math.log(n))
    for _ in range(100):
        rx = R_func(x)
        dx = (n - rx) * math.log(x)
        x += dx
        if abs(dx) < 1e-10: break
    return x


# ============================================================
# PRECOMPUTE: primes, delta, cumulative
# ============================================================

def sieve_primes(limit):
    """Simple sieve of Eratosthenes."""
    is_p = bytearray(b'\x01') * (limit + 1)
    is_p[0] = is_p[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_p[i]:
            is_p[i*i::i] = bytearray(len(is_p[i*i::i]))
    return [i for i in range(2, limit + 1) if is_p[i]]

print("=" * 72)
print("SESSION 4: HIERARCHICAL DECOMPOSITION OF PRIME CORRECTION")
print("=" * 72)

t0 = time.time()
N_MAX = 100000
# Sieve enough primes
print(f"\nSieving primes up to p({N_MAX})...")
# p(100000) ~ 1.3M
primes = sieve_primes(1500000)
if len(primes) < N_MAX:
    primes = sieve_primes(2000000)
print(f"  Found {len(primes)} primes in {time.time()-t0:.2f}s")
print(f"  p(1) = {primes[0]}, p({N_MAX}) = {primes[N_MAX-1]}")

# Compute delta(n) = p(n) - R^{-1}(n)
# Use vectorized approach: R^{-1}(n) ~ n*(ln n + ln ln n) then Newton
print(f"\nComputing delta(n) for n=1..{N_MAX}...")
t0 = time.time()
delta = np.zeros(N_MAX)
rinv = np.zeros(N_MAX)

# Batch compute: start from good initial guess, iterate
ns_arr = np.arange(1, N_MAX + 1, dtype=float)
# Initial guess for R^{-1}(n)
x_arr = np.where(ns_arr > 5,
                 ns_arr * np.log(ns_arr) + ns_arr * np.log(np.maximum(np.log(ns_arr), 1.0)),
                 np.array([2.0, 3.0, 5.0, 7.0, 11.0] + [11.0] * (N_MAX - 5))[:N_MAX])

# Newton iteration in batch (using scalar R_func but with good init)
for i in range(N_MAX):
    n = i + 1
    x = x_arr[i]
    for _ in range(30):  # reduced from 100
        rx = R_func(x)
        dx = (n - rx) * math.log(max(x, 2.0))
        x += dx
        if abs(dx) < 1e-8: break
    rinv[i] = x
    delta[i] = primes[i] - x
    if (i + 1) % 20000 == 0:
        print(f"    {i+1}/{N_MAX} done ({time.time()-t0:.1f}s)")
print(f"  Done in {time.time()-t0:.2f}s")
print(f"  delta range: [{delta.min():.2f}, {delta.max():.2f}]")
print(f"  delta mean: {delta.mean():.4f}, std: {delta.std():.4f}")

ns = np.arange(1, N_MAX + 1, dtype=float)


# ============================================================
# APPROACH A: WAVELET DECOMPOSITION
# ============================================================

def approach_A_wavelets():
    """Decompose delta(n) into wavelet scales."""
    print("\n" + "=" * 72)
    print("APPROACH A: WAVELET DECOMPOSITION OF delta(n)")
    print("=" * 72)

    # Use a power-of-2 length for clean decomposition
    N = 2**16  # 65536
    d = delta[:N]

    # Multi-level wavelet decomposition
    wavelet = 'db8'  # Daubechies-8: good frequency localization
    max_level = pywt.dwt_max_level(N, pywt.Wavelet(wavelet).dec_len)
    levels = min(max_level, 12)
    print(f"\n  Wavelet: {wavelet}, levels: {levels}")

    coeffs = pywt.wavedec(d, wavelet, level=levels)
    # coeffs[0] = approximation (lowest freq), coeffs[1..] = details (high to low freq)

    print(f"\n  Scale decomposition:")
    print(f"  {'Level':<8} {'# Coeffs':<12} {'Energy':<15} {'% of Total':<12} {'Max|c|':<12}")
    print(f"  {'-'*60}")

    total_energy = np.sum(d**2)
    energies = []
    for i, c in enumerate(coeffs):
        e = np.sum(c**2)
        energies.append(e)
        label = f"A{levels}" if i == 0 else f"D{levels - i + 1}"
        print(f"  {label:<8} {len(c):<12} {e:<15.2f} {100*e/total_energy:<12.2f} {np.max(np.abs(c)):<12.4f}")

    # Reconstruct each scale separately
    print(f"\n  Reconstruction by scale:")
    reconstructions = []
    for i in range(len(coeffs)):
        c_zero = [np.zeros_like(c) for c in coeffs]
        c_zero[i] = coeffs[i].copy()
        rec = pywt.waverec(c_zero, wavelet)[:N]
        reconstructions.append(rec)
        rmse = np.sqrt(np.mean((rec)**2))
        label = f"A{levels}" if i == 0 else f"D{levels - i + 1}"
        print(f"  {label}: RMSE={rmse:.4f}, range=[{rec.min():.2f}, {rec.max():.2f}]")

    # KEY TEST: is the approximation (smoothest part) predictable?
    approx = reconstructions[0]  # The smoothest component
    nn = np.arange(1, N + 1, dtype=float)

    # Try fitting the smooth approximation to known functions
    print(f"\n  Fitting smooth approximation (A{levels}) to analytic forms:")

    # Form: a * sqrt(n) * log(n) + b * sqrt(n) + c * n^{1/3} + d
    def model_smooth(n, a, b, c, d, e):
        return a * np.sqrt(n) * np.log(n) + b * np.sqrt(n) + c * n**(1./3) + d * np.log(n) + e

    try:
        popt, pcov = optimize.curve_fit(model_smooth, nn, approx, p0=[0.01, 0.1, 0.1, 1, 0], maxfev=10000)
        fitted = model_smooth(nn, *popt)
        residual = approx - fitted
        r2 = 1 - np.sum(residual**2) / np.sum((approx - approx.mean())**2)
        print(f"    sqrt(n)*log(n) + sqrt(n) + n^(1/3) + log(n) + const:")
        print(f"    Coeffs: {popt}")
        print(f"    R^2 = {r2:.8f}")
        print(f"    Residual RMSE = {np.sqrt(np.mean(residual**2)):.6f}")
    except Exception as e:
        print(f"    Fit failed: {e}")

    # Check autocorrelation of detail coefficients at each scale
    print(f"\n  Autocorrelation (lag-1) of wavelet detail coefficients:")
    for i in range(1, min(len(coeffs), 8)):
        c = coeffs[i]
        if len(c) > 10:
            ac = np.corrcoef(c[:-1], c[1:])[0, 1]
            label = f"D{levels - i + 1}"
            print(f"    {label} (n={len(c)}): autocorr = {ac:.6f}")

    # Check if detail coefficients at each scale follow a pattern
    print(f"\n  Distribution of detail coefficients:")
    for i in range(1, min(len(coeffs), 8)):
        c = coeffs[i]
        if len(c) > 30:
            stat, p = stats.normaltest(c)
            label = f"D{levels - i + 1}"
            print(f"    {label}: mean={np.mean(c):.4f}, std={np.std(c):.4f}, "
                  f"normal_test p={p:.4e} ({'normal' if p > 0.05 else 'NOT normal'})")

    return reconstructions, coeffs


# ============================================================
# APPROACH B: CUMULATIVE CORRECTION
# ============================================================

def approach_B_cumulative():
    """Analyze Sigma(n) = sum_{k=1}^{n} delta(k) - should be smoother."""
    print("\n" + "=" * 72)
    print("APPROACH B: CUMULATIVE CORRECTION Sigma(n) = sum delta(k)")
    print("=" * 72)

    Sigma = np.cumsum(delta)
    nn = ns

    print(f"\n  Sigma(n) range: [{Sigma.min():.2f}, {Sigma.max():.2f}]")
    print(f"  Sigma(100) = {Sigma[99]:.4f}")
    print(f"  Sigma(1000) = {Sigma[999]:.4f}")
    print(f"  Sigma(10000) = {Sigma[9999]:.4f}")
    print(f"  Sigma(100000) = {Sigma[-1]:.4f}")

    # Sigma(n) / n = mean delta up to n
    mean_delta = Sigma / nn
    print(f"\n  Mean delta = Sigma(n)/n:")
    for k in [100, 1000, 10000, 50000, 100000]:
        print(f"    n={k}: {mean_delta[k-1]:.6f}")

    # Key: can we model Sigma(n)?
    # Theoretical: sum of (p(k) - R^{-1}(k)) should relate to theta(x) - x type sums
    # where theta(x) = sum_{p <= x} log p

    # Try fitting: Sigma(n) ~ a * n^alpha + b * n^beta + c * n * log(n) + ...
    print(f"\n  Fitting Sigma(n) to power-law + logarithmic models:")

    # Simple power law: Sigma(n) ~ a * n^alpha
    # Use log-log regression on |Sigma| to find alpha
    valid = np.abs(Sigma) > 0.01
    if np.sum(valid) > 100:
        # Sign of Sigma can flip, so look at the growth envelope
        abs_sig = np.maximum(np.abs(Sigma), 0.01)
        log_n = np.log(nn[1000:])
        log_sig = np.log(abs_sig[1000:])
        slope, intercept, r, p, se = stats.linregress(log_n, log_sig)
        print(f"    Log-log regression (n > 1000): Sigma ~ e^{intercept:.2f} * n^{slope:.4f}")
        print(f"    R^2 = {r**2:.6f}")

    # Model: Sigma(n) = a * n * sqrt(log(n)) + b * n + c * sqrt(n) + d
    def model_sigma(n, a, b, c, d):
        return a * n * np.sqrt(np.log(n)) + b * n + c * np.sqrt(n) + d

    try:
        popt, _ = optimize.curve_fit(model_sigma, nn[10:], Sigma[10:],
                                     p0=[0.1, 0.0, 0.1, 0], maxfev=10000)
        fitted = model_sigma(nn[10:], *popt)
        residual = Sigma[10:] - fitted
        r2 = 1 - np.sum(residual**2) / np.sum((Sigma[10:] - Sigma[10:].mean())**2)
        print(f"\n    Model: a*n*sqrt(log n) + b*n + c*sqrt(n) + d")
        print(f"    Coeffs: a={popt[0]:.6f}, b={popt[1]:.6f}, c={popt[2]:.6f}, d={popt[3]:.6f}")
        print(f"    R^2 = {r2:.10f}")
        print(f"    Residual RMSE = {np.sqrt(np.mean(residual**2)):.4f}")
        print(f"    Residual max = {np.max(np.abs(residual)):.4f}")
    except Exception as e:
        print(f"    Fit failed: {e}")

    # Model: Sigma(n) = a * n^alpha + b * n^beta
    def model_sigma2(n, a, alpha, b, beta, c):
        return a * n**alpha + b * n**beta + c

    try:
        popt2, _ = optimize.curve_fit(model_sigma2, nn[10:], Sigma[10:],
                                      p0=[1.0, 1.0, 0.1, 0.5, 0], maxfev=50000)
        fitted2 = model_sigma2(nn[10:], *popt2)
        residual2 = Sigma[10:] - fitted2
        r2_2 = 1 - np.sum(residual2**2) / np.sum((Sigma[10:] - Sigma[10:].mean())**2)
        print(f"\n    Model: a*n^alpha + b*n^beta + c")
        print(f"    Coeffs: a={popt2[0]:.6f}, alpha={popt2[1]:.6f}, b={popt2[2]:.6f}, beta={popt2[3]:.6f}, c={popt2[4]:.6f}")
        print(f"    R^2 = {r2_2:.10f}")
        print(f"    Residual RMSE = {np.sqrt(np.mean(residual2**2)):.4f}")
    except Exception as e:
        print(f"    Power-law fit failed: {e}")

    # Check: is the RESIDUAL of Sigma(n) (after removing smooth trend) random?
    # This would tell us if there's more structure to exploit
    try:
        resid = residual if 'residual' in dir() else Sigma[10:] - np.polyval(np.polyfit(nn[10:], Sigma[10:], 3), nn[10:])
        # Autocorrelation of residual
        ac1 = np.corrcoef(resid[:-1], resid[1:])[0, 1]
        ac5 = np.corrcoef(resid[:-5], resid[5:])[0, 1]
        ac10 = np.corrcoef(resid[:-10], resid[10:])[0, 1]
        print(f"\n    Residual autocorrelation: lag1={ac1:.6f}, lag5={ac5:.6f}, lag10={ac10:.6f}")
    except:
        pass

    # Compute second difference: delta(n) = Sigma(n) - Sigma(n-1)
    # So delta(n) is recoverable. But can we get delta(n) accurately?
    # If Sigma(n) has relative error eps, then delta(n) has absolute error ~ 2*eps*Sigma(n)/n
    print(f"\n  Error propagation analysis:")
    for n_test in [1000, 10000, 100000]:
        sig_val = abs(Sigma[n_test-1])
        delta_val = abs(delta[n_test-1])
        ratio = sig_val / (n_test * max(delta_val, 0.01))
        print(f"    n={n_test}: |Sigma|={sig_val:.1f}, |delta|={delta_val:.2f}, amplification ratio={ratio:.2f}")

    return Sigma


# ============================================================
# APPROACH C: CHEBYSHEV BIAS EXPLOITATION
# ============================================================

def approach_C_chebyshev():
    """Use residue class biases to constrain p(n)."""
    print("\n" + "=" * 72)
    print("APPROACH C: CHEBYSHEV BIAS EXPLOITATION")
    print("=" * 72)

    # For small moduli, count primes in each residue class
    moduli = [3, 4, 5, 6, 7, 8, 10, 12, 16, 20, 24, 30]
    N_test = N_MAX

    print(f"\n  Prime distribution in residue classes (up to p({N_test})):")
    for m in moduli:
        counts = {}
        for p in primes[:N_test]:
            r = p % m
            if math.gcd(r, m) == 1:  # Only coprime residues
                counts[r] = counts.get(r, 0) + 1
        if len(counts) > 1:
            total = sum(counts.values())
            expected = total / len(counts)
            bias = {r: c - expected for r, c in sorted(counts.items())}
            print(f"    mod {m}: {dict(sorted(counts.items()))}")
            max_bias = max(abs(v) for v in bias.values())
            print(f"      bias: {dict(sorted(bias.items()))}, max_bias={max_bias:.0f} ({100*max_bias/expected:.2f}%)")

    # KEY: Can we predict which residue class p(n) falls in?
    # For mod 4: primes are either 1 or 3 mod 4
    print(f"\n  Residue prediction accuracy (mod 4):")
    mod4 = np.array([p % 4 for p in primes[:N_test]])
    # Running ratio of (3 mod 4) / total
    cum_3 = np.cumsum(mod4 == 3)
    cum_1 = np.cumsum(mod4 == 1)
    for k in [100, 1000, 10000, 50000, N_test]:
        ratio = cum_3[k-1] / k
        print(f"    n={k}: P(p=3 mod 4) = {ratio:.6f}  (bias = {ratio - 0.5:.6f})")

    # For each prime, what fraction of the gap can residue constraints explain?
    print(f"\n  Gap constraint analysis:")
    print(f"  How much of the search space can be eliminated by residue constraints?")

    for m in [6, 30, 210]:
        # Primes > m must be coprime to m
        coprime_to_m = [r for r in range(m) if math.gcd(r, m) == 1]
        density = len(coprime_to_m) / m
        eliminated = 1 - density
        print(f"    mod {m}: {len(coprime_to_m)} coprime residues out of {m}, "
              f"density={density:.4f}, eliminated={100*eliminated:.1f}%")

    # Can we use CRT + biases to narrow candidate window?
    # Test: knowing p(n-1), how many candidates in [p(n-1)+1, p(n-1)+gap]
    # are consistent with all residue constraints?
    print(f"\n  CRT narrowing test (sample of 1000 primes around n=50000):")
    test_range = range(50000, 51000)
    total_gap = 0
    total_candidates = 0
    coprime30 = set(r for r in range(30) if math.gcd(r, 30) == 1)
    for i in test_range:
        if i >= len(primes) - 1:
            break
        gap = primes[i] - primes[i-1]
        # Count integers in gap that are coprime to 30 (=2*3*5)
        candidates_30 = sum(1 for x in range(primes[i-1]+1, primes[i]+1)
                           if (x % 30) in coprime30)
        total_gap += gap
        total_candidates += candidates_30
    avg_gap = total_gap / len(test_range)
    avg_cand = total_candidates / len(test_range)
    print(f"    Average gap: {avg_gap:.2f}")
    print(f"    Average candidates (coprime to 30): {avg_cand:.2f}")
    print(f"    Reduction factor: {avg_cand/avg_gap:.4f}")
    print(f"    NOTE: This is just 8/30 = {8/30:.4f} - no new information beyond coprimality")

    # More sophisticated: use the RUNNING bias to predict residue class
    print(f"\n  Sequential residue prediction (mod 6):")
    mod6 = np.array([p % 6 for p in primes[:N_test]])  # always 1 or 5
    # Predict: if running count of (1 mod 6) > (5 mod 6), predict 5, else 1
    correct = 0
    cum1 = 0
    cum5 = 0
    for i in range(1000, N_test):
        if mod6[i-1] == 1: cum1 += 1
        else: cum5 += 1
        pred = 5 if cum1 >= cum5 else 1  # Predict the less common one (Chebyshev)
        if mod6[i] == pred:
            correct += 1
    accuracy = correct / (N_test - 1000)
    print(f"    Naive Chebyshev prediction accuracy: {accuracy:.6f}")
    print(f"    (Baseline random = 0.5)")


# ============================================================
# APPROACH D: MULTI-SCALE R^{-1} WITH ZETA ZEROS
# ============================================================

# First 30 nontrivial zeros of the Riemann zeta function (imaginary parts)
# rho = 1/2 + i*gamma
ZETA_ZEROS = [
    14.134725141734693, 21.022039638771555, 25.010857580145688,
    30.424876125859513, 32.935061587739189, 37.586178158825671,
    40.918719012147495, 43.327073280914999, 48.005150881167159,
    49.773832477672302, 52.970321477714460, 56.446247697063394,
    59.347044002602353, 60.831778524609809, 65.112544048081606,
    67.079810529494173, 69.546401711173979, 72.067157674481907,
    75.704690699083933, 77.144840068874805, 79.337375020249367,
    82.910380854086030, 84.735492980517050, 87.425274613125229,
    88.809111207634465, 92.491899270558484, 94.651344040519838,
    95.870634228245309, 98.831194218193692, 101.31785100573139,
]

def li_complex(x, rho_im):
    """
    Approximate li(x^rho) where rho = 1/2 + i*gamma.
    li(x^rho) = Ei(rho * log(x))
    For large x, li(x^rho) ~ x^rho / (rho * log(x))
    """
    if x <= 1:
        return 0.0
    lnx = math.log(x)
    # rho = 0.5 + i * gamma
    # x^rho = x^{1/2} * (cos(gamma*lnx) + i*sin(gamma*lnx))
    sqrtx = math.sqrt(x)
    phase = rho_im * lnx
    # li(x^rho) ~ x^rho / (rho * lnx)
    # Real part: Re[x^{1/2} * e^{i*phase} / ((0.5 + i*gamma) * lnx)]
    # = x^{1/2} / lnx * Re[e^{i*phase} / (0.5 + i*gamma)]
    # = x^{1/2} / lnx * Re[(cos(phase) + i*sin(phase)) * (0.5 - i*gamma) / (0.25 + gamma^2)]
    denom = 0.25 + rho_im**2
    real_part = sqrtx / lnx * (0.5 * math.cos(phase) + rho_im * math.sin(phase)) / denom
    imag_part = sqrtx / lnx * (0.5 * math.sin(phase) - rho_im * math.cos(phase)) / denom
    return real_part, imag_part

def R_with_zeros(x, K):
    """
    R(x) with correction from first K zeta zeros.
    pi(x) ~ R(x) - sum_{rho} R(x^rho)
    where sum is over first K pairs of conjugate zeros.
    The contribution of rho and conj(rho) gives -2 * Re[li(x^rho)].
    """
    base = R_func(x)
    correction = 0.0
    for j in range(min(K, len(ZETA_ZEROS))):
        gamma = ZETA_ZEROS[j]
        re_li, _ = li_complex(x, gamma)
        correction += 2 * re_li  # pair of conjugate zeros
    return base - correction

def approach_D_multiscale():
    """Test multi-scale R^{-1} with varying numbers of zeta zeros."""
    print("\n" + "=" * 72)
    print("APPROACH D: MULTI-SCALE R^{-1} WITH ZETA ZEROS")
    print("=" * 72)

    # Test how adding zeta zeros improves accuracy
    test_ns = [100, 500, 1000, 5000, 10000, 50000, 100000]

    print(f"\n  Accuracy of R_K(x) = R(x) - sum_{{j=1}}^{{K}} 2*Re[li(x^rho_j)]")
    print(f"  Measured at x = p(n) for various n and K:")
    print(f"\n  {'n':<10}", end="")
    K_values = [0, 1, 2, 5, 10, 20, 30]
    for K in K_values:
        print(f"{'K='+str(K):<12}", end="")
    print("  actual pi")
    print(f"  {'-'*100}")

    for n in test_ns:
        x = primes[n-1]  # x = p(n), so pi(x) = n
        print(f"  {n:<10}", end="")
        for K in K_values:
            est = R_with_zeros(x, K)
            err = n - est
            print(f"{err:>+10.2f}  ", end="")
        print(f"  {n}")

    # Now the KEY question: can we invert R_K?
    # For each K, compute R_K^{-1}(n) and measure |p(n) - R_K^{-1}(n)|
    print(f"\n  Inversion: |p(n) - R_K^{{-1}}(n)| for various K:")
    print(f"\n  {'n':<10}", end="")
    for K in K_values:
        print(f"{'K='+str(K):<12}", end="")
    print()
    print(f"  {'-'*90}")

    errors_by_K = {K: [] for K in K_values}
    for n in test_ns:
        print(f"  {n:<10}", end="")
        for K in K_values:
            # Newton iteration to find x such that R_K(x) = n
            x = inv_R(n)  # start from plain R^{-1}
            for _ in range(50):
                rk = R_with_zeros(x, K)
                dx = (n - rk) * math.log(x)
                x += dx
                if abs(dx) < 0.001:
                    break
            err = primes[n-1] - x
            errors_by_K[K].append(abs(err))
            print(f"{err:>+10.2f}  ", end="")
        print()

    # Summarize: RMSE for each K
    print(f"\n  Summary RMSE across test points:")
    for K in K_values:
        rmse = np.sqrt(np.mean(np.array(errors_by_K[K])**2))
        print(f"    K={K:>3}: RMSE = {rmse:.4f}")

    # Contribution of each zero to the correction
    print(f"\n  Contribution of individual zeros at n=50000 (x={primes[49999]}):")
    x = primes[49999]
    for j in range(min(30, len(ZETA_ZEROS))):
        gamma = ZETA_ZEROS[j]
        re_li, _ = li_complex(x, gamma)
        print(f"    rho_{j+1} (gamma={gamma:.4f}): 2*Re[li(x^rho)] = {2*re_li:>+10.4f}")

    # Optimal K analysis
    print(f"\n  Finding optimal K for each n:")
    for n in [1000, 10000, 50000, 100000]:
        x = primes[n-1]
        best_K = 0
        best_err = abs(n - R_func(x))
        for K in range(1, 31):
            est = R_with_zeros(x, K)
            err = abs(n - est)
            if err < best_err:
                best_err = err
                best_K = K
        print(f"    n={n:>6}: optimal K={best_K:>3}, residual error={best_err:.4f}")


# ============================================================
# APPROACH E: FOURIER ANALYSIS OF DELTA
# ============================================================

def approach_E_fourier():
    """Look for periodic structure in delta(n)."""
    print("\n" + "=" * 72)
    print("APPROACH E (BONUS): FOURIER ANALYSIS OF delta(n)")
    print("=" * 72)

    N = 2**16
    d = delta[:N]

    # FFT
    D = fft(d)
    freqs = fftfreq(N)
    power = np.abs(D)**2

    # Find dominant frequencies
    half = N // 2
    power_half = power[1:half]
    freqs_half = freqs[1:half]

    # Top 20 frequencies by power
    top_idx = np.argsort(power_half)[-20:][::-1]
    print(f"\n  Top 20 frequencies in delta(n):")
    print(f"  {'Rank':<6} {'Frequency':<15} {'Period':<15} {'Power':<15} {'% Total':<10}")
    total_power = np.sum(power_half)
    for rank, idx in enumerate(top_idx):
        f = abs(freqs_half[idx])
        period = 1.0 / f if f > 0 else float('inf')
        pw = power_half[idx]
        print(f"  {rank+1:<6} {f:<15.8f} {period:<15.2f} {pw:<15.2f} {100*pw/total_power:<10.4f}")

    # Check if any frequency corresponds to known mathematical constants
    print(f"\n  Checking against zeta zeros / (2*pi):")
    for j, gamma in enumerate(ZETA_ZEROS[:10]):
        expected_freq = gamma / (2 * math.pi)
        # Find closest frequency in FFT
        closest_idx = np.argmin(np.abs(freqs_half - expected_freq / N))
        # We need freq in terms of index: period in n-space
        # gamma contributes oscillation with "period" 2*pi/gamma in log(x) space
        # In n-space, this maps to... complicated
        print(f"    rho_{j+1}: gamma={gamma:.4f}, 2*pi/gamma = {2*math.pi/gamma:.4f} (log-space period)")

    # Power spectrum decay rate
    sorted_power = np.sort(power_half)[::-1]
    cumulative = np.cumsum(sorted_power) / total_power
    for frac in [0.5, 0.9, 0.95, 0.99]:
        n_needed = np.searchsorted(cumulative, frac) + 1
        print(f"\n  Frequencies needed for {100*frac:.0f}% of power: {n_needed} out of {half}")


# ============================================================
# SYNTHESIS: MULTI-LAYER RECONSTRUCTION TEST
# ============================================================

def synthesis():
    """Try combining the best elements of all approaches."""
    print("\n" + "=" * 72)
    print("SYNTHESIS: COMBINED MULTI-LAYER RECONSTRUCTION")
    print("=" * 72)

    # Layer 1: R^{-1}(n) with optimal zeros
    # Layer 2: Smooth wavelet approximation correction
    # Layer 3: Residual

    # Test at sample points
    test_ns = list(range(1000, 65001, 5000))
    errors_plain = []
    errors_zeros = []
    errors_combined = []

    # Precompute wavelet smooth correction at full resolution
    N = 2**16
    d = delta[:N]
    coeffs = pywt.wavedec(d, 'db8', level=12)
    # Keep only the smooth approximation
    c_smooth = [coeffs[0]] + [np.zeros_like(c) for c in coeffs[1:]]
    smooth_part = pywt.waverec(c_smooth, 'db8')[:N]

    # Fit smooth part to an analytic formula
    nn_fit = np.arange(1, N + 1, dtype=float)
    def model(n, a, b, c, d, e, f):
        return a * np.sqrt(n) * np.log(n) + b * np.sqrt(n) + c * n**(1./3) + d * np.log(n)**2 + e * np.log(n) + f

    try:
        popt, _ = optimize.curve_fit(model, nn_fit, smooth_part,
                                     p0=[0.01, 0.1, 0.1, 0.01, 1, 0], maxfev=20000)
    except:
        popt = None

    print(f"\n  Testing multi-layer reconstruction at 100 points:")
    print(f"  {'n':<10} {'plain err':<14} {'zeros err':<14} {'combined err':<14} {'actual delta':<14}")
    print(f"  {'-'*66}")

    for n in test_ns:
        idx = n - 1
        # Plain R^{-1}
        ri = rinv[idx]
        err_plain = primes[idx] - ri

        # R^{-1} with 10 zeros
        x = ri
        for _ in range(30):
            rk = R_with_zeros(x, 10)
            dx = (n - rk) * math.log(x)
            x += dx
            if abs(dx) < 0.001: break
        err_zeros = primes[idx] - x

        # Combined: zeros + smooth wavelet correction
        if popt is not None and n <= N:
            smooth_correction = model(float(n), *popt)
            err_combined = err_zeros - smooth_correction
        else:
            err_combined = err_zeros

        errors_plain.append(abs(err_plain))
        errors_zeros.append(abs(err_zeros))
        errors_combined.append(abs(err_combined))

        if n <= 10000 or n % 10000 == 0:
            print(f"  {n:<10} {err_plain:>+12.2f}  {err_zeros:>+12.2f}  {err_combined:>+12.2f}  {delta[idx]:>+12.2f}")

    print(f"\n  RMSE comparison:")
    print(f"    Plain R^{{-1}}:          {np.sqrt(np.mean(np.array(errors_plain)**2)):.4f}")
    print(f"    R^{{-1}} + 10 zeros:     {np.sqrt(np.mean(np.array(errors_zeros)**2)):.4f}")
    print(f"    Zeros + wavelet smooth: {np.sqrt(np.mean(np.array(errors_combined)**2)):.4f}")

    # The fundamental question: what fraction of delta is in each layer?
    print(f"\n  VARIANCE DECOMPOSITION of delta(n):")
    var_total = np.var(delta[:N])

    # Sample zero corrections at 200 points instead of all N
    sample_idx = np.linspace(100, N-1, 200, dtype=int)
    zero_corr_samples = np.zeros(len(sample_idx))
    for si, i in enumerate(sample_idx):
        x = rinv[i]
        base = R_func(x)
        with_zeros = R_with_zeros(x, 10)
        zero_corr_samples[si] = with_zeros - base
    print(f"    Zero correction (10 zeros) variance at 200 samples: {np.var(zero_corr_samples):.4f}")

    var_smooth = np.var(smooth_part)
    var_residual = np.var(d - smooth_part)

    print(f"    Total variance:              {var_total:.4f}")
    print(f"    Wavelet smooth component:    {var_smooth:.4f} ({100*var_smooth/var_total:.1f}%)")
    print(f"    Wavelet detail (residual):   {var_residual:.4f} ({100*var_residual/var_total:.1f}%)")

    # Autocorrelation of residual after removing smooth
    residual = d - smooth_part
    ac1 = np.corrcoef(residual[:-1], residual[1:])[0, 1]
    ac2 = np.corrcoef(residual[:-2], residual[2:])[0, 1]
    ac5 = np.corrcoef(residual[:-5], residual[5:])[0, 1]
    print(f"\n    Residual autocorrelation: lag1={ac1:.6f}, lag2={ac2:.6f}, lag5={ac5:.6f}")
    print(f"    (If near 0 -> residual is unpredictable noise)")
    print(f"    (If near 1 -> more structure to exploit)")


# ============================================================
# FINAL ASSESSMENT
# ============================================================

def final_assessment():
    """Theoretical assessment of each approach for p(10^100)."""
    print("\n" + "=" * 72)
    print("FINAL ASSESSMENT: VIABILITY FOR p(10^100)")
    print("=" * 72)

    print("""
  For p(10^100), we need absolute error < 1 (to identify the exact prime).
  p(10^100) ~ 10^100 * ln(10^100) ~ 2.3 * 10^102

  Current best: R^{-1}(n) gives ~47% of digits -> error ~ 10^55

  What each approach contributes:

  A) WAVELET DECOMPOSITION:
     - The smooth component captures a meaningful trend in delta(n)
     - BUT: the smooth component grows with n and cannot be extrapolated
       reliably to n = 10^100 from data up to n = 10^5
     - VERDICT: Useful for understanding, NOT for computation at 10^100

  B) CUMULATIVE CORRECTION:
     - Sigma(n) is indeed smoother than delta(n)
     - BUT: recovering delta(n) from Sigma requires numerical differentiation,
       which amplifies errors by factor ~n
     - For n=10^100, need relative error < 10^{-102} in Sigma -> IMPOSSIBLE
     - VERDICT: Mathematically interesting, computationally useless

  C) CHEBYSHEV BIAS:
     - Reduces candidate space by constant factor (e.g., 8/30 for mod 30)
     - Does NOT help identify the exact prime - only eliminates ~73% of integers
     - VERDICT: Trivial improvement, does not change complexity class

  D) MULTI-SCALE R^{-1} WITH ZETA ZEROS:
     - Each zero reduces error by O(sqrt(x)/log(x)) contribution
     - With K zeros: error ~ O(sqrt(x) * K^{-alpha}) for some alpha
     - To reach error < 1 at x ~ 10^102: need sqrt(10^102) = 10^51 reduction
     - Would need BILLIONS of zeros (and we only know ~10^13 zeros)
     - Computing zeta zeros at that precision is itself O(T * log T)
     - VERDICT: Best approach theoretically, but need O(10^51) zeros -> impossible

  FUNDAMENTAL BARRIER:
     The error |p(n) - R^{-1}(n)| ~ sqrt(p(n)) * (something involving zeros)
     This is NOT a deficiency of the method - it reflects the ACTUAL
     irregularity of primes. The Riemann Hypothesis says this is O(x^{1/2+eps}).

     To get the exact prime, one needs to account for ALL oscillatory
     contributions from ALL zeta zeros - this is equivalent to computing
     pi(x) exactly, which requires O(x^{2/3}) work.

     CONCLUSION: No hierarchical decomposition can bypass the fundamental
     barrier that the prime correction delta(n) has irreducible complexity
     requiring O(p(n)^{2/3}) computation for exact results.
""")


# ============================================================
# RUN ALL
# ============================================================

if __name__ == '__main__':
    t_start = time.time()

    reconstructions, wavelet_coeffs = approach_A_wavelets()
    Sigma = approach_B_cumulative()
    approach_C_chebyshev()
    approach_D_multiscale()
    approach_E_fourier()
    synthesis()
    final_assessment()

    print(f"\nTotal runtime: {time.time() - t_start:.2f}s")
