"""
Session 4: Smoothed Explicit Formula + Deconvolution

BREAKTHROUGH IDEA:
The explicit formula with Gaussian smoothing converges EXPONENTIALLY:
  ψ_σ(x) = Σ over zeros: contribution ~ exp(-γ²σ²/2)

With σ=0.1, zero at γ=14.13 contributes exp(-14.13²*0.01/2) = exp(-1.0) ≈ 0.37
With σ=0.1, zero at γ=101 contributes exp(-101²*0.01/2) = exp(-51) ≈ 10^{-22}
So only ~5 zeros needed!

BUT: smoothing blurs the step function π(x).
DECONVOLUTION: evaluate smoothed π at multiple points, then deconvolve.

This is a Fourier analysis problem:
  π_smooth(x) = (π * G_σ)(x) where G_σ is Gaussian
  π(x) = F^{-1}[F[π_smooth] / F[G_σ]](x)

But division by F[G_σ] amplifies noise... unless the signal is bandlimited.

ALTERNATIVE: Use the STAIRCASE property.
π(x) is a step function: it's constant except at primes where it jumps by 1.
If we know π_smooth(x) to enough precision, we can recover the steps!

Method:
1. Compute π_smooth(x) at many points using explicit formula (few zeros, FAST)
2. Fit a step function to π_smooth
3. The steps locations ARE the primes!
"""

import mpmath
import time
import math

mpmath.mp.dps = 30

# Zeta zeros
ZEROS = [
    14.134725141734693, 21.022039638771555, 25.010857580145688,
    30.424876125859513, 32.935061587739189, 37.586178158825671,
    40.918719012147495, 43.327073280914999, 48.005150881167159,
    49.773832477672302, 52.970321477714460, 56.446247697063394,
    59.347044002602353, 60.831778524609809, 65.112544048081607,
    67.079810529494174, 69.546401711173980, 72.067157674481907,
    75.704690699083933, 77.144840068874805,
    79.337375020249368, 82.910380854086030, 84.735492980517050,
    87.425274613125229, 88.809111207634465, 92.491899270558484,
    94.651344040519887, 95.870634228245310, 98.831194218193692,
    101.31785100573139,
]

def li_mp(x):
    """li(x) using mpmath."""
    return mpmath.li(x)

def pi_smooth_explicit(x, sigma=0.5, K=20):
    """
    Compute a SMOOTHED version of π(x) using the explicit formula.

    The smoothing kernel is Gaussian: G(t) = exp(-t²/(2σ²)) / (σ√(2π))

    π_σ(x) = ∫ π(x-t) G(t) dt

    Using the explicit formula for Π(x) and Möbius inversion:
    The smoothed version replaces li(x^ρ) with a damped version.

    Actually, the smoothing in log-space is more natural:
    ψ_σ(x) = (1/2πi) ∫ (-ζ'/ζ)(s) * x^s/s * ĥ(s-1) ds
    where ĥ is the Mellin transform of the smoothing kernel.

    For Gaussian smoothing in log-space with parameter σ:
    Each zero ρ = 1/2+iγ gets multiplied by exp(-γ²σ²/2)
    """
    x = mpmath.mpf(x)
    if x < 2:
        return mpmath.mpf(0)

    lnx = mpmath.log(x)

    # Smooth Chebyshev function: ψ_σ(x) = x - Σ_ρ x^ρ/ρ * exp(-γ²σ²/2) + ...
    # Then π_σ(x) ≈ ψ_σ(x) / ln(x) + lower order terms

    # More precisely, use the explicit formula for Π(x):
    # Π_σ(x) = li(x) - Σ_ρ li(x^ρ)*exp(-γ²σ²/2) - ln(2) + small

    result = mpmath.li(x) - mpmath.log(2)

    for k in range(min(K, len(ZEROS))):
        gamma = mpmath.mpf(ZEROS[k])
        damping = mpmath.exp(-gamma**2 * sigma**2 / 2)

        if damping < mpmath.mpf(10)**(-20):
            break  # Negligible contribution

        rho = mpmath.mpc(0.5, gamma)
        li_val = mpmath.ei(rho * lnx)

        result -= 2 * mpmath.re(li_val) * damping

    # Möbius inversion to get π_σ from Π_σ
    # π_σ(x) = Σ_{k=1}^{⌊log₂x⌋} μ(k)/k * Π_σ(x^{1/k})
    # BUT: for the smoothed version, we should smooth π directly,
    # not compute π from smoothed Π. These are different!
    #
    # For now, treat this as an approximation.

    return result  # This is Π_σ(x), close to π_σ(x) for large x

def pi_smooth_mobius(x, sigma=0.5, K=20):
    """π_σ(x) via Möbius inversion of Π_σ(x)."""
    x = mpmath.mpf(x)
    result = mpmath.mpf(0)
    max_k = int(float(mpmath.log(x) / mpmath.log(2))) + 1

    for k in range(1, max_k + 1):
        mu_k = mobius(k)
        if mu_k == 0:
            continue
        xk = mpmath.power(x, mpmath.mpf(1)/k)
        if float(xk) < 2:
            break
        result += mpmath.mpf(mu_k) / k * pi_smooth_explicit(xk, sigma, K)

    return result

def mobius(n):
    if n == 1: return 1
    result = 1
    temp = n
    for p in range(2, int(n**0.5) + 1):
        if temp % p == 0:
            temp //= p
            if temp % p == 0: return 0
            result *= -1
    if temp > 1: result *= -1
    return result

def primes_sieve(n):
    sieve = [True] * (n + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(n**0.5) + 1):
        if sieve[i]:
            for j in range(i*i, n+1, i):
                sieve[j] = False
    return [i for i in range(2, n+1) if sieve[i]]

primes = primes_sieve(120000)

print("=" * 70)
print("SMOOTHED EXPLICIT FORMULA + DECONVOLUTION")
print("=" * 70)

# Test 1: How many zeros needed for different sigma values?
print("\n--- Effective zeros needed for different sigma ---")
for sigma in [0.01, 0.02, 0.05, 0.1, 0.2, 0.5]:
    effective_K = 0
    for k, gamma in enumerate(ZEROS):
        damping = math.exp(-gamma**2 * sigma**2 / 2)
        if damping < 1e-10:
            break
        effective_K = k + 1
    print(f"  σ={sigma:.2f}: {effective_K} effective zeros "
          f"(first zero damping={math.exp(-ZEROS[0]**2*sigma**2/2):.4f})")

# Test 2: Accuracy vs sigma (at x = p(1000))
print("\n--- Accuracy of Π_σ(x) vs σ at x=p(1000)=7919 ---")
n_test = 1000
x_test = primes[n_test - 1]
# Exact Π(x) counts prime powers ≤ x with weight 1/k
Pi_exact = n_test  # approximately (ignoring prime powers p^k for k>1, small correction)
# Better: count prime squares, cubes etc
Pi_exact_val = 0
for p in primes:
    pk = p
    k = 1
    while pk <= x_test:
        Pi_exact_val += 1.0/k
        k += 1
        pk *= p

for sigma in [0.001, 0.005, 0.01, 0.02, 0.05, 0.1]:
    pi_val = float(pi_smooth_mobius(x_test, sigma=sigma, K=30))
    error = abs(pi_val - n_test)
    print(f"  σ={sigma:.3f}: π_σ(p(1000))={pi_val:.4f}, error={error:.4f}")

# Test 3: The STEP DETECTION approach
print("\n--- STEP DETECTION approach ---")
print("Evaluating smoothed π at many points near p(100)=541...")

n_target = 100
x_target = primes[n_target - 1]  # 541

# Evaluate π_σ at points around x_target
sigma = 0.02
points = list(range(x_target - 20, x_target + 20))
values = []
for x in points:
    if x < 2:
        values.append(0)
        continue
    val = float(pi_smooth_mobius(x, sigma=sigma, K=30))
    values.append(val)

# Print the smooth function values
print(f"\nSmoothed π(x) near x={x_target} (σ={sigma}):")
print(f"{'x':>6s} {'π_σ(x)':>12s} {'Δ':>8s} {'actual π(x)':>12s} {'prime?':>7s}")

# Get actual pi values
pi_count = 0
pi_vals = {}
for p in primes:
    if p > max(points):
        break
    pi_count += 1  # this is wrong, let me fix
# Actually compute pi(x) properly
pi_map = {}
idx = 0
for x in range(2, max(points)+1):
    while idx < len(primes) and primes[idx] <= x:
        idx += 1
    pi_map[x] = idx

for i, (x, v) in enumerate(zip(points, values)):
    delta = values[i] - values[i-1] if i > 0 else 0
    actual_pi = pi_map.get(x, '?')
    is_prime = x in set(primes)
    print(f"{x:6d} {v:12.4f} {delta:8.4f} {actual_pi:>12} {'  *' if is_prime else ''}")

# Test 4: Can we detect primes from the derivative of π_σ?
print("\n--- Prime detection from smoothed π derivative ---")
print("If dπ_σ/dx has peaks at primes, we can detect them!")

# Compute numerical derivative
derivs = []
for i in range(1, len(values)):
    derivs.append((values[i] - values[i-1]) / (points[i] - points[i-1]))

# Find peaks
threshold = 1.0 / (sigma * math.sqrt(2 * math.pi) * x_target * math.log(x_target))
print(f"Expected derivative at prime: ~{1.0/(points[0]*math.log(points[0])):.6f}")

peak_positions = []
for i in range(1, len(derivs) - 1):
    if derivs[i] > derivs[i-1] and derivs[i] > derivs[i+1]:
        peak_positions.append(points[i+1])

actual_primes_in_range = [p for p in primes if min(points) <= p <= max(points)]
print(f"Detected peaks at: {peak_positions}")
print(f"Actual primes in range: {actual_primes_in_range}")

# Test 5: For large n, how does the smoothed formula scale?
print("\n--- Scaling test ---")
for n in [100, 500, 1000, 5000, 10000]:
    if n > len(primes):
        break
    x = primes[n-1]
    t0 = time.time()
    val = float(pi_smooth_mobius(x, sigma=0.01, K=30))
    elapsed = time.time() - t0
    error = abs(val - n)
    print(f"  n={n:6d}: π_σ(p(n))={val:.4f}, error={error:.4f}, time={elapsed:.3f}s")

print("\n--- Conclusion ---")
print("The smoothed explicit formula reduces zero requirements but")
print("introduces smoothing error. Step detection might work for")
print("small x but the peaks blur together for large x.")
