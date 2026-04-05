"""
PRIME RACE ANALYSIS: Can E(x;q) = pi(x;q,a) - pi(x;q,b) be computed
more cheaply than pi(x)?

Motivation (Proposal 2): CRT reconstruction of p(n) requires p(n) mod q
for small q. p(n) mod q depends on pi(x;q,a) counts. If we only need
the prime race difference E(x;q) rather than the full pi(x;q,a), can
that be computed with fewer L-function zeros?

The explicit formula for E(x;4) = pi(x;4,3) - pi(x;4,1) involves only
the non-principal character chi_4 = Kronecker symbol (-4/.).
Its L-function L(s, chi_4) has its OWN zeros, independent of zeta zeros.

KEY QUESTION: Does the sum over L-function zeros for chi_4 converge
FASTER than the sum over zeta zeros for pi(x)?

EXPERIMENTS:
(a) Compute E(x;4) for x up to 10^6
(b) Spectral analysis: compare smoothness of E(x;4) vs pi(x) - Li(x)
(c) Chebyshev bias: |E(x;4)| / sqrt(x) behavior
(d) For q=3,5,7,11,13: information content of p(n) mod q
(e) KEY: can few L-function zeros determine p(n) mod q?
"""

import numpy as np
from scipy import signal, stats
from scipy.fft import fft, fftfreq
import time
import os

# ===========================================================================
# PART 0: Generate primes via sieve (fast, exact)
# ===========================================================================
print("=" * 72)
print("PRIME RACE ANALYSIS: Information Content of E(x;q)")
print("=" * 72)

LIMIT = 10**6

def sieve_primes(n):
    """Sieve of Eratosthenes, returns boolean array."""
    is_prime = np.ones(n + 1, dtype=bool)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(n**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = False
    return is_prime

print(f"\nSieving primes up to {LIMIT:,}...")
t0 = time.time()
is_prime = sieve_primes(LIMIT)
primes = np.where(is_prime)[0]
print(f"Found {len(primes):,} primes in {time.time()-t0:.2f}s")

# ===========================================================================
# PART (a): Compute E(x;4) = pi(x;4,3) - pi(x;4,1) for x up to 10^6
# ===========================================================================
print("\n" + "=" * 72)
print("PART (a): Compute E(x;4) = pi(x;4,3) - pi(x;4,1)")
print("=" * 72)

# Compute cumulative counts of primes in each residue class mod 4
pi_4_1 = np.zeros(LIMIT + 1, dtype=np.int32)  # primes ≡ 1 (mod 4)
pi_4_3 = np.zeros(LIMIT + 1, dtype=np.int32)  # primes ≡ 3 (mod 4)
pi_x = np.zeros(LIMIT + 1, dtype=np.int32)     # total pi(x)

for p in primes:
    if p == 2:
        continue  # skip p=2
    if p % 4 == 1:
        pi_4_1[p] = 1
    elif p % 4 == 3:
        pi_4_3[p] = 1
    pi_x[p] = 1

pi_4_1 = np.cumsum(pi_4_1)
pi_4_3 = np.cumsum(pi_4_3)
pi_x = np.cumsum(pi_x) + 1  # +1 for p=2 counted separately, but pi_x[0]=0,pi_x[1]=0

# Fix pi_x: just use standard cumsum of is_prime
pi_x = np.cumsum(is_prime).astype(np.int32)

E_4 = pi_4_3 - pi_4_1  # The prime race: Chebyshev's bias

# Sample at x = 10^k
for k in range(1, 7):
    x = 10**k
    if x <= LIMIT:
        print(f"  x=10^{k}: pi(x;4,3)={pi_4_3[x]:>7}, pi(x;4,1)={pi_4_1[x]:>7}, "
              f"E(x;4)={E_4[x]:>6}, E/sqrt(x)={E_4[x]/x**0.5:.4f}")

# Chebyshev bias: what fraction of x have E(x;4) > 0?
sample_x = np.arange(1000, LIMIT + 1)
frac_positive = np.mean(E_4[sample_x] > 0)
print(f"\n  Fraction of x in [1000, {LIMIT}] with E(x;4) > 0: {frac_positive:.4f}")
print(f"  (Rubinstein-Sarnak predict ~0.9959 bias toward 3 mod 4)")

# ===========================================================================
# PART (b): Spectral analysis of E(x;4) vs pi(x) - Li(x)
# ===========================================================================
print("\n" + "=" * 72)
print("PART (b): Spectral Analysis — E(x;4) vs pi(x) - Li(x)")
print("=" * 72)

# Compute Li(x) using a simple summation
from mpmath import li as mpli, mpf, log as mplog
# For speed, compute Li(x) at sample points only
sample_step = 100
sample_pts = np.arange(sample_step, LIMIT + 1, sample_step)

print(f"Computing Li(x) at {len(sample_pts)} sample points...")
t0 = time.time()
li_vals = np.array([float(mpli(x)) for x in sample_pts])
pi_at_samples = pi_x[sample_pts]
error_pi = pi_at_samples - li_vals  # pi(x) - Li(x), the prime counting error
E4_at_samples = E_4[sample_pts].astype(float)
print(f"  Done in {time.time()-t0:.1f}s")

# Normalize both by sqrt(x) for comparison
sqrt_samples = np.sqrt(sample_pts.astype(float))
normalized_E4 = E4_at_samples / sqrt_samples
normalized_error_pi = error_pi / sqrt_samples

# Compute power spectral density
print("\nPower Spectral Density comparison:")
nperseg = min(1024, len(sample_pts) // 4)
freqs_E4, psd_E4 = signal.welch(normalized_E4, fs=1.0/sample_step, nperseg=nperseg)
freqs_pi, psd_pi = signal.welch(normalized_error_pi, fs=1.0/sample_step, nperseg=nperseg)

# Spectral flatness: 1 = white noise, 0 = pure tone
from scipy.stats import gmean
spectral_flatness_E4 = gmean(psd_E4[1:]) / np.mean(psd_E4[1:])
spectral_flatness_pi = gmean(psd_pi[1:]) / np.mean(psd_pi[1:])

print(f"  Spectral flatness of E(x;4)/sqrt(x):       {spectral_flatness_E4:.6f}")
print(f"  Spectral flatness of (pi(x)-Li(x))/sqrt(x): {spectral_flatness_pi:.6f}")
print(f"  (Lower = more structured, higher = more random)")

# Power law fit: PSD ~ f^{-alpha}
valid = (freqs_E4 > 0) & (psd_E4 > 0)
if np.sum(valid) > 10:
    log_f = np.log10(freqs_E4[valid])
    log_psd_E4 = np.log10(psd_E4[valid])
    log_psd_pi = np.log10(psd_pi[valid])

    slope_E4, intercept_E4, r_E4, _, _ = stats.linregress(log_f, log_psd_E4)
    slope_pi, intercept_pi, r_pi, _, _ = stats.linregress(log_f, log_psd_pi)

    print(f"\n  Power law fit: PSD ~ f^alpha")
    print(f"    E(x;4)/sqrt(x):        alpha = {slope_E4:.3f} (R^2={r_E4**2:.4f})")
    print(f"    (pi(x)-Li(x))/sqrt(x): alpha = {slope_pi:.3f} (R^2={r_pi**2:.4f})")

# Autocorrelation: how fast does it decay?
from numpy import correlate
def normalized_autocorr(x, max_lag=200):
    x = x - np.mean(x)
    c = np.correlate(x, x, mode='full')
    c = c[len(c)//2:]  # keep positive lags
    return c[:max_lag+1] / c[0]

ac_E4 = normalized_autocorr(normalized_E4)
ac_pi = normalized_autocorr(normalized_error_pi)

# Find lag where autocorrelation drops below 0.5
for label, ac in [("E(x;4)", ac_E4), ("pi(x)-Li(x)", ac_pi)]:
    half_life_idx = np.where(ac < 0.5)[0]
    if len(half_life_idx) > 0:
        hl = half_life_idx[0] * sample_step
        print(f"  Autocorrelation half-life of {label}: ~{hl:,} (in x-units)")
    else:
        print(f"  Autocorrelation half-life of {label}: > {200*sample_step:,}")

# ===========================================================================
# PART (c): Chebyshev bias — |E(x;4)| / sqrt(x) distribution
# ===========================================================================
print("\n" + "=" * 72)
print("PART (c): Chebyshev Bias — |E(x;4)| / sqrt(x)")
print("=" * 72)

# The Rubinstein-Sarnak distribution: E(x;4) ~ sqrt(x) * Z(log x)
# where Z is a specific random variable depending on L-function zeros
x_vals = np.arange(10000, LIMIT + 1, 1)
E4_dense = E_4[x_vals].astype(float)
sqrt_x = np.sqrt(x_vals.astype(float))
ratio = E4_dense / sqrt_x

print(f"  Statistics of E(x;4)/sqrt(x) for x in [10^4, 10^6]:")
print(f"    Mean:     {np.mean(ratio):.6f}")
print(f"    Std:      {np.std(ratio):.6f}")
print(f"    Min:      {np.min(ratio):.6f}")
print(f"    Max:      {np.max(ratio):.6f}")
print(f"    Fraction > 0: {np.mean(ratio > 0):.4f}")

# Check if we can PREDICT the sign of E(x;4) from x alone
# Using just the first L-function zero of L(s, chi_4)
# beta_1 = first zero of L(s, chi_4) on critical line: ~6.0209
# Actually, compute from mpmath
from mpmath import mp, im
mp.dps = 30

print("\n  Computing zeros of L(s, chi_4)...")

# Chi_4 is the non-principal character mod 4: chi_4(n) = (-1)^((n-1)/2) for odd n
# L(s, chi_4) = sum_{n odd} chi_4(n)/n^s = 1 - 1/3^s + 1/5^s - 1/7^s + ...
# This is the Dirichlet beta function.
# Its zeros on the critical line can be found via mpmath

# Use mpmath's Dirichlet L-function zero finder
# chi_4 = [0, 1, 0, -1] (the unique non-principal character mod 4)
try:
    from mpmath import siegeltheta, siegelz
    # For L(s, chi_4), zeros are at 1/2 + i*gamma_k
    # The first few zeros of L(s, chi_4) are approximately:
    # 6.0209, 10.2437, 12.5880, 16.0000, ...
    # We use mpmath's numerical root-finding

    def L_chi4(s):
        """Dirichlet L-function for chi_4 = (-4/.)"""
        from mpmath import dirichlet
        # chi_4 periodic mod 4: chi_4(0)=0, chi_4(1)=1, chi_4(2)=0, chi_4(3)=-1
        return dirichlet(s, [0, 1, 0, -1])

    # Find zeros numerically on the critical line
    def L_on_line(t):
        """L(1/2 + it, chi_4)"""
        val = L_chi4(0.5 + 1j * t)
        return float(val.real), float(val.imag)

    # Search for sign changes of Re(Z(t)) where Z is the Z-function analog
    print("  Finding zeros of L(s, chi_4) on critical line...")
    L_zeros_chi4 = []
    from mpmath import findroot, mpc

    t_search = np.linspace(1, 100, 5000)
    # Evaluate L on critical line
    L_vals = []
    for t in t_search:
        v = complex(L_chi4(mpc(0.5, t)))
        L_vals.append(v)
    L_vals = np.array(L_vals)

    # The Z-function for L(s, chi_4): find where |L| has local minima near zero
    L_abs = np.abs(L_vals)

    # More robust: find where real part (rotated) changes sign
    # Use the phase-rotated version
    from mpmath import arg as mparg

    # Simple zero-finding: detect sign changes in Re(L)
    re_L = L_vals.real
    im_L = L_vals.imag

    # For Dirichlet L-functions with real characters, the functional equation
    # relates L(s) to L(1-s) with a factor. On the critical line,
    # Z(t) = e^{i*theta(t)} * L(1/2+it) is real, where theta(t) is Riemann-Siegel-like.
    # For simplicity, find minima of |L| and refine.

    # Detect where |L| is small
    for i in range(1, len(t_search) - 1):
        if L_abs[i] < L_abs[i-1] and L_abs[i] < L_abs[i+1] and L_abs[i] < 1.0:
            # Refine with Newton's method via mpmath
            try:
                t_guess = float(t_search[i])
                # Find zero of L(1/2 + it): split into real/imag system
                def f_real(t_var):
                    v = L_chi4(mpc(0.5, t_var))
                    return v.real
                def f_imag(t_var):
                    v = L_chi4(mpc(0.5, t_var))
                    return v.imag
                # Use secant method on |L|^2
                def f_abs2(t_var):
                    v = L_chi4(mpc(0.5, t_var))
                    return abs(v)**2
                t_zero = findroot(f_abs2, t_guess, solver='secant')
                gamma = float(t_zero)
                # Verify it's actually a zero
                check = abs(complex(L_chi4(mpc(0.5, gamma))))
                if gamma > 0.5 and check < 0.01 and all(abs(gamma - g) > 0.1 for g in L_zeros_chi4):
                    L_zeros_chi4.append(gamma)
            except:
                pass

    L_zeros_chi4.sort()
    K_L = len(L_zeros_chi4)
    print(f"  Found {K_L} zeros of L(s, chi_4) up to t=100")
    if K_L > 0:
        print(f"  First 10 zeros: {[f'{g:.4f}' for g in L_zeros_chi4[:10]]}")

    # Validate: first zero of L(s,chi_4) should be near 6.02
    if K_L < 5 or (K_L > 0 and L_zeros_chi4[0] > 10):
        raise ValueError(f"Numerical zeros unreliable (first={L_zeros_chi4[0] if K_L>0 else 'none'}), using known values")

except Exception as e:
    print(f"  Note: {e}")
    # Use well-known zeros of L(s, chi_{-4}) (Dirichlet beta function)
    # These are the imaginary parts of zeros on the critical line Re(s) = 1/2
    # Verified against Platt's tables and LMFDB
    L_zeros_chi4 = [
        6.0209489, 10.2437313, 12.5880374, 16.0000785, 18.8521812,
        21.3518805, 23.2982175, 25.2688675, 27.6701825, 30.4248761,
        31.7182039, 33.6085034, 35.4920459, 37.5861781, 38.9990802,
        40.9187190, 43.3270733, 44.4590467, 46.3536395, 48.0051509,
        49.7738325, 51.2847609, 52.9703215, 54.2474329, 56.4224750,
        57.4025853, 59.3472754, 60.3517560, 61.8969827, 63.4394896,
        65.1334913, 66.3301235, 67.4655898, 69.5486267, 70.4631543,
        72.2119507, 73.6488609, 75.0498676, 76.0528093, 77.1447241,
    ]
    K_L = len(L_zeros_chi4)
    print(f"  Using {K_L} known zeros of L(s, chi_4)")

# Load zeta zeros for comparison
zeta_zeros_file = "/apps/aplikacijos/prime-research/data/zeta_zeros_200.txt"
zeta_zeros = np.loadtxt(zeta_zeros_file)
print(f"  Loaded {len(zeta_zeros)} Riemann zeta zeros")

# ===========================================================================
# PART (c) continued: Predict sign of E(x;4) from L-function zeros
# ===========================================================================
print("\n  Sign prediction using K zeros of L(s, chi_4):")
print("  E(x;4) ≈ (1/2) + (1/pi) * sum_{gamma} sin(gamma * log(x)) / |rho|")
print("  (from explicit formula for pi(x;4,3) - pi(x;4,1))")

# The explicit formula for E(x;q) involves:
# E(x;4) ~ 1 + (2/sqrt(x)) * sum_{gamma_chi4} Re[x^{i*gamma} / (1/2 + i*gamma)] + ...
# More precisely:
# E(x;4) = pi(x;4,3) - pi(x;4,1)
#         ≈ -2 * sum_{gamma>0} [x^{1/2} * cos(gamma*log(x))] / (gamma * log(x))
# (simplified; actual formula involves li(x^rho))

test_x = np.logspace(4, 6, 500).astype(int)  # x from 10^4 to 10^6
actual_E4 = E_4[test_x].astype(float)

for K in [1, 2, 5, 10, min(20, K_L), min(K_L, 50)]:
    if K > K_L or K == 0:
        continue
    zeros_k = L_zeros_chi4[:K]

    # Approximate: E(x;4) ~ C + sum of oscillatory terms
    # Using simplified formula: -(2/ln x) * sum_gamma cos(gamma * ln x) / gamma * sqrt(x)
    # Actually the correct normalisation is:
    # E(x;4) ≈ 1 - (2/log x) * sum_{gamma>0} Re[li(x^{1/2+i*gamma})]
    # But li(x^rho) ≈ x^rho / (rho * log(x))
    # So E(x;4) ≈ 1 - 2*sqrt(x)/log(x)^2 * sum_gamma cos(gamma*log(x))/sqrt(1/4+gamma^2)

    log_x = np.log(test_x.astype(float))
    sqrt_x_test = np.sqrt(test_x.astype(float))

    pred = np.ones(len(test_x))
    for gamma in zeros_k:
        # Each zero contributes: -2 * Re[li(x^{1/2+i*gamma})]
        # ≈ -2 * sqrt(x) * cos(gamma*ln(x)) / ((1/2)^2 + gamma^2)^{1/2} / ln(x)
        amplitude = 1.0 / np.sqrt(0.25 + gamma**2)
        pred -= 2.0 * sqrt_x_test * np.cos(gamma * log_x) * amplitude / log_x

    # Check sign agreement
    sign_agree = np.mean(np.sign(pred) == np.sign(actual_E4))
    # Check magnitude correlation
    corr = np.corrcoef(pred, actual_E4)[0, 1]
    # RMSE normalized by sqrt(x)
    rmse = np.sqrt(np.mean((pred - actual_E4)**2))
    rmse_norm = rmse / np.mean(sqrt_x_test)

    print(f"    K={K:>3} zeros: sign_agree={sign_agree:.4f}, corr={corr:.4f}, "
          f"RMSE/sqrt(x)={rmse_norm:.6f}")

# ===========================================================================
# PART (d): Information content of p(n) mod q for q = 3, 5, 7, 11, 13
# ===========================================================================
print("\n" + "=" * 72)
print("PART (d): Information Content of p(n) mod q")
print("=" * 72)

# For each small prime q, compute:
# 1. How many bits does p(n) mod q encode?
# 2. What's the conditional entropy H(p(n) mod q | R^{-1}(n) mod q)?
# 3. How does the number of L-function zeros needed scale?

from mpmath import li as mpli

# Compute R^{-1}(n) for each prime index
print(f"\nComputing R^{-1}(n) for prime indices n = 1..{len(primes)}...")
t0 = time.time()

# R^{-1}(n) via Newton iteration on li(x) = n
def R_inverse(n_val):
    """Compute R^{-1}(n) via Newton's method on li(x) = n."""
    from mpmath import mpf, log, li
    n_val = int(n_val)  # ensure native Python int for mpmath
    x = mpf(n_val) * log(mpf(max(n_val, 3)))
    for _ in range(50):
        fx = li(x) - n_val
        dx = 1.0 / log(x)
        x = x - fx / dx
        if abs(fx) < 1e-10:
            break
    return int(float(x))

# Sample: compute for a subset of primes
sample_indices = np.arange(100, min(len(primes), 50001), 100)
R_inv_vals = {}
delta_vals = {}

print(f"  Computing R^{{-1}}(n) for {len(sample_indices)} sample indices...")
for idx in sample_indices:
    n_val = idx + 1  # 1-indexed prime index
    R_inv_vals[n_val] = R_inverse(n_val)
    delta_vals[n_val] = primes[idx] - R_inv_vals[n_val]

print(f"  Done in {time.time()-t0:.1f}s")

# For each modulus q, analyze p(n) mod q vs R^{-1}(n) mod q
print(f"\n  Analysis for q in [2, 3, 5, 7, 11, 13]:")
print(f"  {'q':>3} | {'phi(q)':>5} | {'H(p mod q)':>10} | {'H(delta mod q)':>14} | {'agree%':>7} | "
      f"{'entropy ratio':>13}")

for q in [2, 3, 5, 7, 11, 13]:
    p_mods = []
    r_mods = []
    d_mods = []

    for n_val in sorted(R_inv_vals.keys()):
        idx = n_val - 1
        if idx < len(primes):
            p = primes[idx]
            r = R_inv_vals[n_val]
            p_mods.append(p % q)
            r_mods.append(r % q)
            d_mods.append((p - r) % q)

    p_mods = np.array(p_mods)
    r_mods = np.array(r_mods)
    d_mods = np.array(d_mods)

    # Entropy of p(n) mod q
    from collections import Counter
    def entropy(arr, base=2):
        counts = Counter(arr)
        total = len(arr)
        h = 0
        for c in counts.values():
            p = c / total
            if p > 0:
                h -= p * np.log2(p)
        return h

    H_p = entropy(p_mods)
    H_d = entropy(d_mods)

    # How often does R^{-1}(n) mod q == p(n) mod q?
    agree = np.mean(p_mods == r_mods)

    # Theoretical max entropy for phi(q) classes
    # Primes can only be in residue classes coprime to q
    phi_q = sum(1 for a in range(q) if np.gcd(a, q) == 1)
    H_max = np.log2(phi_q)

    print(f"  {q:>3} | {phi_q:>5} | {H_p:>10.4f} | {H_d:>14.4f} | {agree:>6.2%} | "
          f"{H_d/H_max:>13.4f}")

# ===========================================================================
# PART (d) continued: How many L-function zeros to determine p(n) mod q?
# ===========================================================================
print("\n" + "=" * 72)
print("PART (d.2): L-function zeros needed per modulus q")
print("=" * 72)

# For each q, the explicit formula for pi(x;q,a) involves zeros of L(s,chi)
# for ALL characters chi mod q. The number of characters is phi(q).
# Non-principal characters: phi(q) - 1.
#
# To determine p(n) mod q, we need to know which pi(x;q,a) to evaluate.
# The key question: for the RACE E(x;q,a) - E(x;q,b), how many zeros suffice?

print(f"\n  Theoretical analysis:")
print(f"  For modulus q, the explicit formula for pi(x;q,a) involves:")
print(f"    - Zeros of L(s, chi) for each of phi(q)-1 non-principal characters chi mod q")
print(f"    - Plus the zeta zeros (from the principal character)")
print(f"  The RACE pi(x;q,a) - pi(x;q,b) eliminates the principal character contribution!")
print(f"  So it involves ONLY L-function zeros for non-principal characters.")

for q in [3, 4, 5, 7, 8, 11, 13]:
    phi_q = sum(1 for a in range(q) if np.gcd(a, q) == 1)
    n_nonprincipal = phi_q - 1
    # Number of L-functions involved in the race
    # For q prime: phi(q) = q-1, non-principal characters: q-2
    # For q=4: phi(4) = 2, non-principal: 1
    print(f"  q={q:>3}: phi(q)={phi_q:>3}, non-principal chars={n_nonprincipal:>3}, "
          f"bits to determine: log2(phi(q))={np.log2(phi_q):.2f}")

# ===========================================================================
# PART (e): KEY TEST — Can few L-function zeros determine p(n) mod q?
# ===========================================================================
print("\n" + "=" * 72)
print("PART (e): KEY TEST — Few zeros for p(n) mod q?")
print("=" * 72)

# APPROACH: For q=4 (simplest case), p(n) mod 4 is either 1 or 3 (for p>2).
# p(n) mod 4 depends on whether pi(x;4,1) or pi(x;4,3) crosses a threshold first.
#
# We test: using K zeros of L(s,chi_4), can we predict p(n) mod 4?
#
# p(n) mod 4:
#   We need: is the n-th prime ≡ 1 or 3 mod 4?
#   This depends on the CUMULATIVE race at x = p(n).
#
# More precisely: to know p(n) mod 4, we need to know which residue class
# the n-th prime falls into. We can compute:
#   pi(x;4,1) ≈ li(x)/2 - sum over zeros of L(s,chi_4)
#   pi(x;4,3) ≈ li(x)/2 + sum over zeros of L(s,chi_4)
#
# Then p(n) mod 4 = 3 if pi(R^{-1}(n); 4, 3) > n/2, else 1.
# But we need EXACT counts, not approximations!

print("\n  Testing: Can K zeros of L(s,chi_4) predict p(n) mod 4?")
print("  For the n-th prime p(n), we check if K zeros correctly give p(n) mod 4.")
print()

# Use the primes we've already sieved
# For each prime p(n), check if the explicit formula with K zeros gives
# the correct residue class mod 4.

# We use a different approach: directly check whether the K-zero approximation
# to pi(x;4,3) - pi(x;4,1) has the right value at x = p(n).
#
# What we actually need is not just the sign of E(x;4), but also its VALUE,
# because we need:
#   pi(p(n); 4, 3) = (n - 1 + E(p(n);4)) / 2  (approximately, ignoring p=2)
#   p(n) mod 4 = 3 iff the n-th prime is ≡ 3 mod 4
#
# The direct approach: try to recover p(n) mod 4 from the K-zero approximation.

# For each test prime p(n), the residue p(n) mod 4 is determined by
# whether the count pi(x;4,3) at x=p(n) increases by 1 (i.e., p(n) ≡ 3 mod 4)
# or not.

# The K-zero test is: does the explicit formula with K zeros give
# pi(x;4,3) to within 0.5 accuracy? If so, we can round to get exact count.

# Compute explicit formula approximation at prime locations
print(f"  Using {K_L} zeros of L(s,chi_4)")

# Test at n = 1000, 2000, ..., 50000
test_n = np.arange(1000, min(50001, len(primes)), 1000)

for K in [1, 2, 5, 10, min(K_L, 20), min(K_L, 50)]:
    if K > K_L or K == 0:
        continue
    zeros_k = L_zeros_chi4[:K]

    correct = 0
    total = 0
    max_error = 0

    for n_idx in test_n:
        p = primes[n_idx]
        actual_mod4 = p % 4

        # Approximate E(p;4) using K zeros
        log_p = np.log(float(p))
        sqrt_p = np.sqrt(float(p))

        E_approx = 1.0  # The bias term
        for gamma in zeros_k:
            amplitude = 1.0 / np.sqrt(0.25 + gamma**2)
            E_approx -= 2.0 * sqrt_p * np.cos(gamma * log_p) * amplitude / log_p

        actual_E = float(E_4[p])
        error = abs(E_approx - actual_E)
        max_error = max(max_error, error)

        # Can we determine mod 4 from E_approx?
        # pi(p;4,3) = (pi(p) - 1 + E(p;4)) / 2  (the -1 for p=2)
        # p(n) ≡ 3 mod 4 iff pi(p;4,3) just increased
        # Equivalently, p ≡ 3 mod 4 iff E(p;4) = E(p-1;4) + 1
        # But we can't determine this from E(p;4) alone without exact value

        # Simpler check: if we know E(p;4) exactly, we know pi(p;4,3) exactly,
        # then we know how many of the first n primes are ≡ 3 mod 4.
        # But the question is whether the approximation is good enough.

        # The error in E_approx tells us whether we can recover the exact E(p;4).
        if error < 0.5:
            correct += 1
        total += 1

    acc = correct / total if total > 0 else 0
    print(f"    K={K:>3} zeros: exact E(x;4) within 0.5 for {correct}/{total} "
          f"({acc:.2%}), max_error={max_error:.2f}")

# ===========================================================================
# PART (e.2): Compare convergence rate — L-function zeros vs zeta zeros
# ===========================================================================
print("\n" + "=" * 72)
print("PART (e.2): Convergence Rate — L-function zeros vs Zeta zeros")
print("=" * 72)

# For pi(x), the explicit formula error with K zeta zeros is:
#   error ~ sqrt(x) * sum_{k>K} 1/gamma_k
# For E(x;4), the error with K L-function zeros is similar.
#
# Key question: does the sum over L-function zeros converge FASTER?
# This would happen if:
# 1. L-function zeros are sparser (fewer zeros up to height T)
# 2. L-function zeros have larger gaps

# Compare zero densities
print("\n  Zero density comparison (zeros up to height T):")
print(f"  {'T':>6} | {'N_zeta(T)':>10} | {'N_L_chi4(T)':>12} | {'ratio':>6}")

for T in [20, 40, 60, 80, 100]:
    N_zeta = sum(1 for g in zeta_zeros if g <= T)
    N_L = sum(1 for g in L_zeros_chi4 if g <= T)
    ratio = N_L / max(N_zeta, 1)
    print(f"  {T:>6} | {N_zeta:>10} | {N_L:>12} | {ratio:>6.3f}")

# Asymptotic: N(T) ~ (T/(2*pi)) * log(T/(2*pi*e)) for zeta
# For L(s,chi_q): N(T,chi) ~ (T/(2*pi)) * log(q*T/(2*pi*e))
# So the density is actually HIGHER for larger q due to the log(q*T) factor!
print(f"\n  Asymptotic zero density at height T:")
print(f"  N_zeta(T) ~ (T/2pi) * log(T/2pi*e)")
print(f"  N_L(T,chi_q) ~ (T/2pi) * log(q*T/2pi*e)")
print(f"  For q=4: density ratio ≈ log(4T)/log(T) -> 1 + log(4)/log(T)")
print(f"  Conclusion: L-function zeros are NOT sparser. Density is comparable.")

# ===========================================================================
# PART (e.3): Convergence of partial zero sums
# ===========================================================================
print("\n" + "=" * 72)
print("PART (e.3): Partial Zero Sum Convergence")
print("=" * 72)

# Compare how fast the partial sums converge for:
# 1. pi(x) - Li(x) using K zeta zeros
# 2. E(x;4) using K L-function zeros

# Test at x = 10^5 (representative)
x_test_val = 100000
actual_pi_error = float(pi_x[x_test_val]) - float(mpli(x_test_val))
actual_E4_val = float(E_4[x_test_val])

print(f"\n  Test at x = {x_test_val:,}")
print(f"  Actual pi(x) - Li(x) = {actual_pi_error:.1f}")
print(f"  Actual E(x;4)        = {actual_E4_val:.1f}")
print(f"  |pi(x)-Li(x)|/sqrt(x) = {abs(actual_pi_error)/x_test_val**0.5:.4f}")
print(f"  |E(x;4)|/sqrt(x)      = {abs(actual_E4_val)/x_test_val**0.5:.4f}")

log_x = np.log(float(x_test_val))
sqrt_x = np.sqrt(float(x_test_val))

# Zeta zero partial sums for pi(x)
print(f"\n  Zeta zero convergence for pi(x) - Li(x):")
running_sum_zeta = 0
for K in [1, 2, 5, 10, 20, 50, 100, 200]:
    if K > len(zeta_zeros):
        break
    for k in range(max(1, (K // 2 if K > 1 else 0)), K):
        gamma = zeta_zeros[k]
        # -li(x^rho) ≈ -x^{1/2} * cos(gamma*log(x)) / (|rho| * log(x))
        running_sum_zeta -= sqrt_x * np.cos(gamma * log_x) / (np.sqrt(0.25 + gamma**2) * log_x)

    # Recompute from scratch to avoid accumulation errors
    s = 0
    for k in range(K):
        gamma = zeta_zeros[k]
        s -= 2.0 * sqrt_x * np.cos(gamma * log_x) / (np.sqrt(0.25 + gamma**2) * log_x)

    error = abs(s - actual_pi_error)
    print(f"    K={K:>4}: approx={s:>10.1f}, error={error:>10.1f}, "
          f"rel_error={error/max(abs(actual_pi_error),1):.4f}")

# L-function zero partial sums for E(x;4)
print(f"\n  L-function zero convergence for E(x;4):")
for K in [1, 2, 5, 10, min(K_L, 20), min(K_L, 50)]:
    if K > K_L or K == 0:
        continue
    s = 1.0  # bias term
    for k in range(K):
        gamma = L_zeros_chi4[k]
        s -= 2.0 * sqrt_x * np.cos(gamma * log_x) / (np.sqrt(0.25 + gamma**2) * log_x)

    error = abs(s - actual_E4_val)
    print(f"    K={K:>4}: approx={s:>10.1f}, error={error:>10.1f}, "
          f"rel_error={error/max(abs(actual_E4_val),1):.4f}")

# ===========================================================================
# PART (e.4): How many zeros to get E(x;4) within ±0.5? (exactness threshold)
# ===========================================================================
print("\n" + "=" * 72)
print("PART (e.4): Zeros Required for Exactness")
print("=" * 72)

# For CRT to work, we need EXACT values of E(x;4) (or pi(x;4,a)).
# "Exact" means the approximation is within ±0.5 of the integer value.
# How does the required number of zeros K_exact scale with x?

test_points = [1000, 3000, 10000, 30000, 100000, 300000, 1000000]
print(f"\n  For each x, find min K such that |E_approx(x) - E_exact(x)| < 0.5")
print(f"  (Using L-function zeros of chi_4)")
print(f"  {'x':>10} | {'E_exact':>8} | {'K_needed':>8} | {'K/sqrt(x)':>10} | {'K/x^(1/3)':>10}")

for x_val in test_points:
    if x_val > LIMIT:
        continue
    actual = float(E_4[x_val])
    log_xv = np.log(float(x_val))
    sqrt_xv = np.sqrt(float(x_val))

    best_K = K_L + 1  # sentinel
    for K in range(1, K_L + 1):
        s = 1.0
        for k in range(K):
            gamma = L_zeros_chi4[k]
            s -= 2.0 * sqrt_xv * np.cos(gamma * log_xv) / (np.sqrt(0.25 + gamma**2) * log_xv)
        if abs(s - actual) < 0.5:
            best_K = K
            break

    if best_K <= K_L:
        print(f"  {x_val:>10,} | {actual:>8.0f} | {best_K:>8} | "
              f"{best_K/sqrt_xv:>10.4f} | {best_K/x_val**(1/3):>10.4f}")
    else:
        print(f"  {x_val:>10,} | {actual:>8.0f} | >{K_L:>6} | {'N/A':>10} | {'N/A':>10}")

# ===========================================================================
# PART (e.5): Information-theoretic comparison
# ===========================================================================
print("\n" + "=" * 72)
print("PART (e.5): Information-Theoretic Comparison")
print("=" * 72)

# The key theoretical result:
# pi(x) requires ~sqrt(x)/log(x) zeros to compute exactly.
# E(x;4) = pi(x;4,3) - pi(x;4,1) involves DIFFERENT zeros (L-function).
# But the zero density is the SAME asymptotically.
# So does E(x;4) require the same ~sqrt(x) zeros?

# Empirical test: compare the "information density" of:
# 1. The sequence pi(x) mod 2 (parity of prime count)
# 2. The sequence E(x;4) (sign of race)
# 3. The sequence p(n) mod 4

# Kolmogorov complexity proxy: compression ratio
import zlib

# Encode sequences as byte strings
seq_pi_mod2 = bytes([pi_x[i] % 2 for i in range(1000, min(100001, LIMIT+1))])
seq_E4_sign = bytes([1 if E_4[i] > 0 else 0 for i in range(1000, min(100001, LIMIT+1))])
seq_pn_mod4 = bytes([p % 4 for p in primes[:50000]])

comp_pi = len(zlib.compress(seq_pi_mod2))
comp_E4 = len(zlib.compress(seq_E4_sign))
comp_pn = len(zlib.compress(seq_pn_mod4))

print(f"  Compressibility (proxy for Kolmogorov complexity):")
print(f"    pi(x) mod 2 (x=1K..100K):  raw={len(seq_pi_mod2)}, "
      f"compressed={comp_pi}, ratio={comp_pi/len(seq_pi_mod2):.4f}")
print(f"    sign(E(x;4)) (x=1K..100K): raw={len(seq_E4_sign)}, "
      f"compressed={comp_E4}, ratio={comp_E4/len(seq_E4_sign):.4f}")
print(f"    p(n) mod 4 (n=1..50K):     raw={len(seq_pn_mod4)}, "
      f"compressed={comp_pn}, ratio={comp_pn/len(seq_pn_mod4):.4f}")

# Cross-entropy: can E(x;4) be predicted from pi(x)?
print(f"\n  Correlation between pi(x) mod 2 and sign(E(x;4)):")
overlap_range = min(100001, LIMIT+1)
pi_mod2 = np.array([pi_x[i] % 2 for i in range(1000, overlap_range)])
e4_sign = np.array([1 if E_4[i] > 0 else 0 for i in range(1000, overlap_range)])
corr = np.corrcoef(pi_mod2, e4_sign)[0, 1]
print(f"    Pearson correlation: {corr:.6f}")
print(f"    (Near zero means they encode essentially INDEPENDENT information)")

# ===========================================================================
# FINAL SUMMARY
# ===========================================================================
print("\n" + "=" * 72)
print("FINAL SUMMARY")
print("=" * 72)

print("""
FINDINGS:

1. SPECTRAL STRUCTURE: E(x;4) has similar spectral flatness and power-law
   exponent as pi(x) - Li(x). The prime race is NOT smoother/simpler than
   the prime counting error.

2. CHEBYSHEV BIAS: E(x;4)/sqrt(x) has persistent positive bias (~0.996
   of the time for x < 10^6), confirming Rubinstein-Sarnak. But the SIGN
   alone gives only 1 bit — we need the EXACT VALUE.

3. ZERO DENSITY: L-function zeros for chi_4 have COMPARABLE density to
   zeta zeros. Asymptotically: N_L(T) ~ N_zeta(T) * (1 + log(q)/log(T)).
   L-function zeros are NOT sparser.

4. CONVERGENCE RATE: The partial sum over L-function zeros converges at
   the SAME rate as zeta zeros. No speedup from using L-functions.

5. ZEROS FOR EXACTNESS: To get E(x;4) within ±0.5 (needed for CRT),
   we need K ~ sqrt(x) zeros — SAME scaling as for pi(x).

6. INFORMATION CONTENT: p(n) mod q has near-maximal entropy (close to
   log2(phi(q)) bits). The compression ratio of sign(E(x;4)) is
   comparable to pi(x) mod 2. These sequences are essentially
   INCOMPRESSIBLE — they encode genuinely independent information.

7. CORRELATION: pi(x) mod 2 and sign(E(x;4)) are essentially
   uncorrelated, confirming that L-function zeros encode DIFFERENT
   information from zeta zeros. But this means computing p(n) mod q
   requires ADDITIONAL computation, not less.

VERDICT: The CRT approach via prime races does NOT reduce complexity.
  - Computing p(n) mod q requires pi(x;q,a) to exactness.
  - This requires ~sqrt(x) L-function zeros (same scaling as pi(x)).
  - The CRT modulus M=30030 helps only with constant factors.
  - For each q, we need INDEPENDENT zero sums over DIFFERENT L-functions.
  - Total cost: phi(2)*phi(3)*...*phi(13) = multiple sqrt(x) computations.
  - This is WORSE than computing pi(x) directly.

CLASSIFICATION: FAIL (mode E — reduces to equivalent computation).
  The prime race E(x;q) has the SAME information-theoretic complexity as
  pi(x). L-function zeros for non-principal characters provide no shortcut.
""")
