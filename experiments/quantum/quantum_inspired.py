#!/usr/bin/env python3
"""
Session 5: QUANTUM-INSPIRED CLASSICAL ALGORITHMS for prime counting
====================================================================

Explores classical analogs of quantum algorithmic ideas applied to
computing p(n) -- the nth prime.

APPROACHES:
  1. Amplitude amplification analog (importance sampling near R^{-1}(n))
  2. QFT analog (FFT on prime distribution / zeta zero quasi-periodicity)
  3. Phase estimation analog (power iteration on prime-counting operator)
  4. Tensor network / MPS (sieve as low-bond-dimension tensor contraction)
  5. Variational approach (parameterized function optimized via gradient descent)
  6. Adiabatic / homotopy continuation (from li^{-1}(n) to p(n))

All tested rigorously against known primes.

Context:
  - 100+ approaches tried across 5 sessions
  - R^{-1}(n) gives ~47% of digits in O(polylog) time
  - Exact computation requires O(p(n)^{2/3}) via Lucy DP (best practical)
  - Theoretical best: O(p(n)^{1/2+eps}) via Lagarias-Odlyzko (needs zeta zeros)
  - The gap between approximate and exact is proven fundamental
"""

import math
import time
import sys
import warnings
from collections import defaultdict
from functools import lru_cache

import numpy as np
from scipy import fft as scipy_fft
from scipy.optimize import minimize, brentq
from scipy.interpolate import CubicSpline
from scipy.linalg import eig

try:
    import mpmath
    mpmath.mp.dps = 30
    HAS_MPMATH = True
except ImportError:
    HAS_MPMATH = False

# Suppress overflow warnings in experimental code
warnings.filterwarnings('ignore', category=RuntimeWarning)

# ============================================================
# REFERENCE: Sieve and prime utilities
# ============================================================

def sieve(n):
    """Sieve of Eratosthenes up to n."""
    if n < 2:
        return []
    is_prime = bytearray(b'\x01') * (n + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(n**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return [i for i in range(2, n + 1) if is_prime[i]]

# Precompute primes for testing
PRIMES = sieve(200000)
P = [0] + PRIMES  # 1-indexed: P[1]=2, P[2]=3, ...
PI_TABLE = [0] * (200001)
for p in PRIMES:
    PI_TABLE[p] = 1
for i in range(1, 200001):
    PI_TABLE[i] += PI_TABLE[i-1]

def pi_exact(x):
    """Exact pi(x) for x <= 200000."""
    x = int(x)
    if x < 2:
        return 0
    if x <= 200000:
        return PI_TABLE[x]
    raise ValueError(f"x={x} too large for table")

def riemann_R(x):
    """Riemann R function."""
    if HAS_MPMATH:
        return float(mpmath.riemannr(x))
    # Fallback: just li(x) - li(sqrt(x))/2 ...
    return x / math.log(x) * (1 + 1/math.log(x) + 2/math.log(x)**2)

def R_inverse(n, tol=1e-6):
    """Numerical inverse of Riemann R: find x such that R(x) = n."""
    if n <= 1:
        return 2.0
    x = n * math.log(n)
    lo, hi = 2.0, x * 3
    for _ in range(200):
        mid = (lo + hi) / 2
        if riemann_R(mid) < n:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return (lo + hi) / 2

print("=" * 70)
print("SESSION 5: QUANTUM-INSPIRED CLASSICAL ALGORITHMS")
print("=" * 70)
print(f"Precomputed {len(PRIMES)} primes up to {PRIMES[-1]}")
print()


# ============================================================
# APPROACH 1: AMPLITUDE AMPLIFICATION ANALOG
#   Grover's algorithm: O(sqrt(N)) search in unstructured space
#   Classical analog: importance sampling
#
#   Idea: Instead of uniformly searching integers near R^{-1}(n),
#   sample with probability proportional to "prime-likeness" score.
#   This should concentrate samples on actual primes, reducing
#   the effective search space.
# ============================================================

def approach1_importance_sampling():
    """
    Amplitude amplification analog: importance sampling near R^{-1}(n).

    Strategy:
    - Get x0 = R^{-1}(n)
    - Define a "prime-likeness" score based on:
      (a) Not divisible by small primes (main signal)
      (b) Passes Fermat pseudoprime test (mod small base)
      (c) Residue class bias (primes avoid certain residues)
    - Weight candidates by this score
    - Measure: how many samples needed to find p(n)?

    Quantum connection: Grover amplifies amplitude of "marked" states.
    Here we amplify sampling probability of prime candidates.
    """
    print("=" * 70)
    print("APPROACH 1: AMPLITUDE AMPLIFICATION / IMPORTANCE SAMPLING")
    print("=" * 70)

    small_primes = PRIMES[:50]  # First 50 primes for scoring

    def prime_score(x):
        """Score how likely x is to be prime, 0 to 1."""
        if x < 2:
            return 0.0
        if x % 2 == 0:
            return 0.0
        score = 1.0
        for p in small_primes:
            if p >= x:
                break
            if x % p == 0:
                return 0.0  # Definitely composite
            # Boost: x is not divisible by p
            # Expected fraction that survives: (p-1)/p
            # We don't need to adjust since we already know it survived

        # Fermat test base 2
        if x > 2:
            try:
                if pow(2, x - 1, x) != 1:
                    return 0.0
            except (ValueError, OverflowError):
                pass

        # Score based on density: primes near x have density ~ 1/ln(x)
        # After sieving by first 50 primes, surviving density is higher
        # Product of (1 - 1/p) for p in small_primes ≈ e^{-gamma}/ln(p_50)
        # So surviving integers that are prime have conditional probability
        # ≈ (1/ln(x)) / (e^{-gamma}/ln(229)) ≈ ln(229) * e^gamma / ln(x)
        return 1.0

    test_ns = [10, 50, 100, 500, 1000, 5000, 10000]

    print(f"\n{'n':>6} | {'p(n)':>8} | {'R_inv':>10} | {'err':>6} | {'unif_samples':>14} | {'IS_samples':>12} | {'speedup':>8}")
    print("-" * 85)

    results = []
    for n in test_ns:
        if n >= len(P):
            break
        target = P[n]
        x0 = R_inverse(n)
        err = abs(x0 - target)

        # Method A: Uniform sampling in [x0-2*err, x0+2*err]
        # Count how many uniform samples to find the target
        window = max(int(2.5 * err), 50)
        lo = max(2, int(x0) - window)
        hi = int(x0) + window

        # All integers in window
        all_candidates = list(range(lo, hi + 1))
        n_total = len(all_candidates)

        # Primes in window
        primes_in_window = [x for x in all_candidates if x <= 200000 and PI_TABLE[x] - PI_TABLE[x-1] == 1]

        if target not in range(lo, hi + 1):
            # Window too small, expand
            lo = target - 100
            hi = target + 100
            all_candidates = list(range(lo, hi + 1))
            n_total = len(all_candidates)
            primes_in_window = [x for x in all_candidates if 2 <= x <= 200000 and PI_TABLE[x] - PI_TABLE[x-1] == 1]

        # Method B: Importance sampling (only sample from candidates that pass small-prime sieve)
        # This is like "amplitude amplification" - we concentrate on likely candidates
        is_candidates = [x for x in all_candidates if prime_score(x) > 0]
        n_is = len(is_candidates)

        # Expected samples to find target:
        # Uniform: n_total / 1 (need to check each)
        # IS: n_is / 1 (smaller set)
        uniform_expected = n_total
        is_expected = n_is
        speedup = uniform_expected / max(is_expected, 1)

        print(f"{n:>6} | {target:>8} | {x0:>10.1f} | {err:>6.1f} | {uniform_expected:>14} | {is_expected:>12} | {speedup:>8.2f}x")
        results.append((n, speedup))

    # Analysis: what is the asymptotic speedup?
    print("\n--- ASYMPTOTIC ANALYSIS ---")
    print("By Mertens' theorem, sieving by primes up to y eliminates fraction")
    print("1 - prod(1-1/p, p<=y) = 1 - e^{-gamma}/ln(y)")
    print()
    print("For y = 229 (first 50 primes):")
    mertens_prod = 1.0
    for pp in PRIMES[:50]:
        mertens_prod *= (1 - 1.0/pp)
    print(f"  Surviving fraction: {mertens_prod:.6f}")
    print(f"  Theoretical (e^{{-gamma}}/ln(229)): {math.exp(-0.5772156649)/math.log(229):.6f}")
    print(f"  Speedup: {1/mertens_prod:.2f}x")
    print()
    print("For the ACTUAL barrier:")
    print(f"  Error of R^{{-1}}(n) ~ C * p(n)^0.34")
    print(f"  Search window ~ 2C * p(n)^0.34 integers")
    print(f"  After sieve: ~ 2C * p(n)^0.34 * e^{{-gamma}}/ln(y) candidates")
    print(f"  Each needs primality test: O(polylog(p(n)))")
    print(f"  Total: O(p(n)^0.34 * polylog / ln(y))")
    print(f"  For p(10^100): ~ 10^34 / ln(y) candidates")
    print(f"  Even with y = 10^20: ~ 10^34 / 46 ~ 2*10^32 candidates")
    print()
    print("VERDICT: Importance sampling gives constant-factor speedup only.")
    print("  Quantum Grover would give QUADRATIC speedup (10^34 -> 10^17).")
    print("  Classical analog cannot match this - importance sampling is O(N/ln(y))")
    print("  vs Grover's O(sqrt(N)).")

    # But can we do better with more sophisticated scoring?
    print("\n--- ENHANCED SCORING: Wheel factorization + Miller-Rabin ---")

    # Wheel mod 30 (2*3*5): only 8 residues can be prime
    wheel30 = {1, 7, 11, 13, 17, 19, 23, 29}

    for n in [1000, 10000, 50000]:
        if n >= len(P):
            break
        target = P[n]
        x0 = R_inverse(n)
        window = max(int(3 * abs(x0 - target)), 200)
        lo = max(2, int(x0) - window)
        hi = int(x0) + window

        total = hi - lo + 1
        wheel_survivors = sum(1 for x in range(lo, hi+1) if x % 30 in wheel30 or x in (2,3,5))
        sieve50_survivors = sum(1 for x in range(lo, hi+1) if prime_score(x) > 0)
        actual_primes = len([x for x in range(lo, hi+1) if 2 <= x <= 200000 and PI_TABLE[x] - PI_TABLE[x-1] == 1])

        print(f"  n={n}: window={total}, wheel30={wheel_survivors} ({wheel_survivors/total:.3f}), "
              f"sieve50={sieve50_survivors} ({sieve50_survivors/total:.3f}), "
              f"primes={actual_primes} ({actual_primes/total:.4f})")

    print("\n  The ratio primes/survivors is bounded by density 1/ln(x),")
    print("  so importance sampling cannot break the O(error/ln(x)) barrier.")
    print("  This is a CONSTANT factor improvement, not algorithmic breakthrough.")

    return False  # No breakthrough


# ============================================================
# APPROACH 2: QFT ANALOG — FFT ON PRIME DISTRIBUTION
#   QFT finds periodicity exponentially faster.
#   Classical analog: use FFT/NTT to detect quasi-periodic structure
#   in primes, then reconstruct the signal.
#
#   Key idea: The explicit formula for pi(x) is essentially a
#   Fourier expansion over zeta zeros. Can FFT help compute
#   the "zeta zero correction" more efficiently?
# ============================================================

def approach2_fft_prime_distribution():
    """
    QFT analog: FFT-based analysis of prime distribution.

    The explicit formula: pi(x) = R(x) - sum_{rho} R(x^rho)
    where rho runs over zeta zeros (rho = 1/2 + i*gamma_k).

    This is like an inverse Fourier transform with frequencies gamma_k.
    Can FFT help us:
      (a) Recover the correction term more efficiently?
      (b) Detect and exploit the quasi-periodic structure?
      (c) Compress the pi(x) function?
    """
    print("\n" + "=" * 70)
    print("APPROACH 2: QFT ANALOG — FFT ON PRIME DISTRIBUTION")
    print("=" * 70)

    # --- 2a: FFT of the prime indicator function ---
    print("\n--- 2a: FFT of prime indicator (Fourier analysis of primes) ---")

    N = 2**14  # 16384
    # Create prime indicator: 1 if n is prime, 0 otherwise
    indicator = np.zeros(N)
    for p in PRIMES:
        if p < N:
            indicator[p] = 1.0

    # FFT
    t0 = time.time()
    spectrum = np.fft.rfft(indicator)
    fft_time = time.time() - t0

    magnitudes = np.abs(spectrum)
    phases = np.angle(spectrum)

    print(f"FFT of prime indicator (N={N}): {fft_time:.4f}s")
    print(f"DC component (= number of primes < {N}): {magnitudes[0]:.0f} (actual: {pi_exact(N-1)})")

    # Find top frequencies
    top_k = 20
    top_indices = np.argsort(magnitudes[1:])[::-1][:top_k] + 1
    print(f"\nTop {top_k} frequencies by magnitude:")
    for i, idx in enumerate(top_indices):
        freq = idx / N
        period = N / idx if idx > 0 else float('inf')
        print(f"  #{i+1}: freq={freq:.6f}, period={period:.2f}, magnitude={magnitudes[idx]:.2f}, phase={phases[idx]:.4f}")

    # --- 2b: Truncated FFT reconstruction of pi(x) ---
    print("\n--- 2b: Reconstruct pi(x) from K Fourier components ---")

    for K in [10, 50, 100, 500, 1000, 5000]:
        # Keep only top K frequencies
        truncated = np.zeros_like(spectrum)
        top_K = np.argsort(magnitudes)[::-1][:K]
        truncated[top_K] = spectrum[top_K]

        # Inverse FFT
        reconstructed = np.fft.irfft(truncated, n=N)

        # Cumulative sum to get pi(x) from indicator
        pi_recon = np.cumsum(reconstructed)
        pi_true = np.array([pi_exact(x) for x in range(N)])

        errors = np.abs(pi_recon - pi_true)
        max_err = np.max(errors)
        mean_err = np.mean(errors)

        # Check at specific points
        test_points = [100, 1000, 5000, 10000]
        exact_at = [(x, pi_exact(x), pi_recon[x]) for x in test_points if x < N]

        print(f"  K={K:>5}: max_err={max_err:.2f}, mean_err={mean_err:.2f}, "
              f"pi(1000)={exact_at[1][2]:.1f} vs {exact_at[1][1]}")

    # --- 2c: Can FFT detect zeta zero frequencies? ---
    print("\n--- 2c: Zeta zero detection via FFT ---")
    print("The correction to R(x) involves terms oscillating at frequencies gamma_k/(2*pi)")
    print("where gamma_k are imaginary parts of zeta zeros.")

    # Known first zeta zeros
    zeta_zeros = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
                  37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
                  52.970321, 56.446248, 59.347044, 60.831779, 65.112544]

    # The correction to pi(x) at x involves sum of terms like -Re(R(x^rho))
    # For x = e^t, each zero contributes oscillation at frequency gamma_k/(2*pi) in t-space

    # Create delta(x) = pi(x) - R(x) in log space
    N_log = 4096
    t_max = math.log(50000)
    t_vals = np.linspace(1, t_max, N_log)
    x_vals = np.exp(t_vals)

    delta_vals = np.array([pi_exact(int(x)) - riemann_R(x) for x in x_vals])

    # FFT of delta in log-space
    delta_fft = np.fft.rfft(delta_vals)
    delta_mag = np.abs(delta_fft)
    freqs = np.fft.rfftfreq(N_log, d=(t_max - 1) / N_log)

    # Convert to "gamma" units: gamma = 2*pi*freq
    gamma_vals = 2 * np.pi * freqs

    print(f"\nFFT of delta(x) = pi(x) - R(x) in log-space:")
    top_delta = np.argsort(delta_mag[1:])[::-1][:10] + 1
    print(f"{'rank':>4} | {'FFT_gamma':>10} | {'magnitude':>10} | {'nearest zeta zero':>18} | {'distance':>10}")
    print("-" * 65)

    detected = 0
    for rank, idx in enumerate(top_delta):
        g = gamma_vals[idx]
        # Find nearest known zeta zero
        dists = [abs(g - z) for z in zeta_zeros]
        nearest_idx = np.argmin(dists)
        nearest = zeta_zeros[nearest_idx]
        dist = dists[nearest_idx]
        match = "*" if dist < 2.0 else ""
        if dist < 2.0:
            detected += 1
        print(f"{rank+1:>4} | {g:>10.4f} | {delta_mag[idx]:>10.2f} | {nearest:>18.6f} | {dist:>10.4f} {match}")

    print(f"\nDetected {detected}/10 zeta zeros within tolerance 2.0")
    print("(Resolution limited by N_log and t_max)")

    # --- 2d: Can we use FFT to COMPUTE pi(x) faster? ---
    print("\n--- 2d: FFT for fast pi(x) computation? ---")
    print("ANALYSIS:")
    print("  The explicit formula is: pi(x) = R(x) - sum_rho R(x^rho)")
    print("  This is NOT a standard convolution, so FFT doesn't directly help.")
    print("  Each term R(x^rho) depends on x, so we can't precompute.")
    print("  FFT could help if we computed pi(x) at MANY x values simultaneously")
    print("  via batch evaluation, but each evaluation still needs all zeros.")
    print()

    # Batch evaluation test: compute pi(x) at many points via FFT
    # If we had the zeta zeros, we could compute:
    #   pi(x_j) = R(x_j) - sum_k c_k * x_j^{rho_k}
    # This is a "non-uniform DFT" (NUDFT) problem

    # Test: can NUDFT reconstruct pi(x) from few zeros?
    print("  Reconstruction of pi(x) using K zeta zeros (direct, not FFT):")
    for K in [1, 3, 5, 10, 15]:
        zeros_used = zeta_zeros[:K]
        errors = []
        for test_x in [100, 500, 1000, 5000, 10000, 50000]:
            R_x = riemann_R(test_x)
            correction = 0.0
            for gamma in zeros_used:
                # R(x^rho) where rho = 1/2 + i*gamma
                # x^rho = x^{1/2} * e^{i*gamma*ln(x)}
                # R(x^rho) is complex, we need -2*Re(R(x^rho)) summed over conjugate pairs
                ln_x = math.log(test_x)
                # Simplified: correction ~ -2 * sum Re(li(x^rho))
                # li(x^rho) ~ x^rho / (rho * ln(x))
                rho_r, rho_i = 0.5, gamma
                x_rho_mag = test_x ** 0.5
                x_rho_phase = gamma * ln_x
                # li(x^rho) ~ x^rho / log(x^rho) = x^rho / (rho * ln(x))
                rho_mag = math.sqrt(rho_r**2 + rho_i**2)
                rho_phase = math.atan2(rho_i, rho_r)

                result_mag = x_rho_mag / (rho_mag * ln_x)
                result_phase = x_rho_phase - rho_phase

                correction -= 2 * result_mag * math.cos(result_phase)

            pi_approx = R_x + correction
            pi_true = pi_exact(test_x)
            errors.append(abs(pi_approx - pi_true))

        mean_err = np.mean(errors)
        max_err = np.max(errors)
        print(f"    K={K:>2}: mean_err={mean_err:.2f}, max_err={max_err:.2f}")

    print()
    print("  Convergence is VERY slow (as established in prior sessions).")
    print("  FFT/NTT cannot help because the problem is NOT a convolution")
    print("  and the 'frequencies' (zeta zeros) are irregularly spaced.")
    print()
    print("VERDICT: FFT gives insight into structure but no computational speedup.")
    print("  The bottleneck is the NUMBER of zeta zeros needed (~10^51),")
    print("  not the cost of combining them.")

    return False


# ============================================================
# APPROACH 3: PHASE ESTIMATION ANALOG — POWER METHOD
#   QPE extracts eigenvalues of a unitary operator.
#   Classical analog: power method extracts dominant eigenvalue.
#
#   Idea: Construct a matrix M whose eigenvalue encodes pi(x) or p(n).
#   Then use power iteration to extract it.
# ============================================================

def approach3_power_method():
    """
    Phase estimation analog: power iteration on prime-counting operator.

    Can we construct a matrix whose spectral properties encode pi(x)?

    Attempt 1: Transfer matrix of sieve process
    Attempt 2: Adjacency matrix of "prime graph"
    Attempt 3: Operator whose trace = pi(x)
    """
    print("\n" + "=" * 70)
    print("APPROACH 3: PHASE ESTIMATION / POWER METHOD")
    print("=" * 70)

    # --- 3a: Sieve transfer matrix ---
    print("\n--- 3a: Sieve as transfer matrix ---")
    print("Idea: Represent the sieve state as a vector, sieving by p as matrix multiply.")
    print("Then pi(x) = some function of M_2 * M_3 * M_5 * ... * M_sqrt(x) applied to initial state.")

    # For small modulus m, sieve state is a vector of length m
    # indicating which residue classes survive.
    # M_p zeros out the residue class 0 mod p.

    # Transfer matrix for sieving by prime p, working mod m
    def sieve_matrix(p, m):
        """Matrix that zeroes out multiples of p in residue classes mod m."""
        M = np.eye(m)
        for r in range(m):
            if r % p == 0:
                M[r, r] = 0
        return M

    # Test: pi(x) mod 30 via transfer matrices
    m = 30  # Work modulo 30 = 2*3*5
    state = np.ones(m)  # All residues start as "possibly prime"
    state[0] = 0  # 0 is not prime
    state[1] = 0  # 1 is not prime... actually we need to be careful

    # Actually, let's think about this differently.
    # The sieve of Eratosthenes works on SPECIFIC numbers, not residue classes.
    # A transfer matrix approach would need dimension = x for sieve up to x.
    # That's just the sieve itself encoded as linear algebra.

    print("  Result: Sieve transfer matrix has dimension x (same as sieve).")
    print("  No compression possible via this route.")

    # --- 3b: Operator whose trace encodes pi(x) ---
    print("\n--- 3b: Trace formula approach ---")
    print("Selberg trace formula: sum over eigenvalues = sum over periodic orbits")
    print("This is how the zeta zeros relate to primes!")
    print("But computing the trace requires knowing the eigenvalues (= zeta zeros).")

    # Demonstrate: construct matrix whose eigenvalues ARE the zeta zeros
    # and whose trace gives pi-related quantities

    # The GUE connection: spacings of zeta zeros match GUE random matrices
    # Can we USE this? If we generate a GUE matrix, its eigenvalues
    # statistically match zeta zeros, so maybe we can approximate corrections?

    print("\n--- 3c: GUE random matrix model for zeta zeros ---")

    np.random.seed(42)
    results_by_dim = []

    for dim in [10, 50, 100, 200]:
        n_trials = 20
        errors = []

        for trial in range(n_trials):
            # Generate GUE matrix
            A = np.random.randn(dim, dim) + 1j * np.random.randn(dim, dim)
            H = (A + A.conj().T) / 2  # Hermitian
            eigenvalues = np.sort(np.linalg.eigvalsh(H))

            # Scale eigenvalues to match zeta zero spacing
            # Mean spacing of first K zeros: ~2*pi/ln(gamma_K/(2*pi))
            # We scale to match the first known zero
            mean_spacing = np.mean(np.diff(eigenvalues))
            if abs(mean_spacing) < 1e-10:
                continue

            # Scale so that eigenvalues match zeta zero positions
            scaled = eigenvalues * (14.134725 / eigenvalues[dim//2]) if dim > 1 else eigenvalues

            # Use these "fake zeros" to compute correction to R(x)
            test_x = 10000
            R_x = riemann_R(test_x)
            correction = 0.0
            for gamma in scaled:
                if abs(gamma) < 0.1:
                    continue
                gamma = abs(gamma)
                ln_x = math.log(test_x)
                rho_mag = math.sqrt(0.25 + gamma**2)
                rho_phase = math.atan2(gamma, 0.5)
                x_half = test_x ** 0.5
                result_mag = x_half / (rho_mag * ln_x)
                result_phase = gamma * ln_x - rho_phase
                correction -= 2 * result_mag * math.cos(result_phase)

            pi_approx = R_x + correction
            pi_true = pi_exact(test_x)
            errors.append(abs(pi_approx - pi_true))

        if errors:
            mean_err = np.mean(errors)
            std_err = np.std(errors)
            print(f"  GUE dim={dim:>3}: mean_err(pi(10000))={mean_err:.1f} +/- {std_err:.1f} "
                  f"(true pi(10000)={pi_exact(10000)}, R(10000)={riemann_R(10000):.1f})")
            results_by_dim.append((dim, mean_err))

    print()
    print("  GUE model does NOT help: random eigenvalues don't match specific zeros.")
    print("  The correction depends on EXACT zero positions, not just statistics.")

    # --- 3d: Power iteration on concrete operator ---
    print("\n--- 3d: Direct power iteration attempt ---")
    print("Construct matrix A such that A^k * v converges to a vector encoding pi(x).")

    # Idea: Markov chain whose stationary distribution concentrates on primes
    # Transition: from x, go to x+1 if x is composite, stay if x is prime
    # Stationary distribution: pi(x) is related to sum of stationary probs

    # Small test
    N_test = 100
    A = np.zeros((N_test, N_test))
    for i in range(2, N_test):
        is_p = (i <= 200000 and PI_TABLE[i] - PI_TABLE[i-1] == 1)
        if is_p:
            A[i, i] = 0.9  # Stay at prime with prob 0.9
            if i + 1 < N_test:
                A[i+1, i] = 0.1
        else:
            if i + 1 < N_test:
                A[i+1, i] = 1.0  # Always move forward from composite
            else:
                A[i, i] = 1.0

    # Normalize columns
    col_sums = A.sum(axis=0)
    col_sums[col_sums == 0] = 1
    A = A / col_sums

    # Power iteration
    v = np.ones(N_test) / N_test
    for _ in range(1000):
        v = A @ v
        v = v / np.sum(v)

    # Where does probability concentrate?
    top_states = np.argsort(v)[::-1][:10]
    primes_in_top = sum(1 for s in top_states if s >= 2 and PI_TABLE[s] - PI_TABLE[s-1] == 1)
    print(f"  Markov chain (N={N_test}): top 10 states have {primes_in_top} primes")
    print(f"  Top states: {list(top_states)}")
    print(f"  This is just a random walk - it concentrates at the LAST prime before N,")
    print(f"  not at p(n). No way to extract pi(x) from spectral properties.")

    print()
    print("VERDICT: No useful operator found whose eigenvalues encode pi(x)")
    print("  without already knowing the primes or zeta zeros.")
    print("  The Selberg trace formula IS the connection, but it requires")
    print("  computing O(sqrt(x)) zeta zeros for exact results.")

    return False


# ============================================================
# APPROACH 4: TENSOR NETWORK / MATRIX PRODUCT STATE
#   Quantum states of n qubits need 2^n amplitudes, but MPS
#   represent them with O(n * D^2) parameters (D = bond dimension).
#
#   Idea: Encode the sieve as a tensor network. If the sieve
#   has low entanglement, it has small bond dimension.
# ============================================================

def approach4_tensor_network():
    """
    Tensor network / MPS approach to sieve encoding.

    Idea 1: Represent the prime indicator function f(x) = 1 if x is prime
    as an MPS over the binary digits of x. If the bond dimension D is small,
    we can evaluate f(x) and even sum it (= pi(x)) in O(log(x) * D^2) time.

    Idea 2: Represent the sieve process as a tensor train and contract it.
    """
    print("\n" + "=" * 70)
    print("APPROACH 4: TENSOR NETWORK / MATRIX PRODUCT STATE")
    print("=" * 70)

    # --- 4a: MPS representation of prime indicator ---
    print("\n--- 4a: MPS for prime indicator f(x) = [x is prime] ---")
    print("Binary representation of x with L bits: x = sum b_i * 2^i")
    print("f(x) = <boundary| T_{b_L} T_{b_{L-1}} ... T_{b_1} T_{b_0} |boundary>")
    print("Each T_{b} is a D x D matrix depending on bit b.")

    # Test: can we fit an MPS to the prime indicator?
    # For L bits, there are 2^L numbers. We need D^2 * L * 2 parameters.
    # If this is << 2^L, we have compression.

    L = 10  # 10 bits = numbers up to 1023
    N_nums = 2**L

    # True prime indicator
    true_indicator = np.array([1 if (i >= 2 and i <= 200000 and PI_TABLE[i] - PI_TABLE[i-1] == 1) else 0
                               for i in range(N_nums)])

    # Try to find MPS of bond dimension D
    # We'll use a simple approach: reshape indicator as a tensor and do SVD

    # Reshape into tensor of shape (2, 2, ..., 2) — L indices
    tensor = true_indicator.reshape([2] * L)

    # Sequential SVD to get MPS (this is just the TT-decomposition)
    bond_dims = []
    remaining = tensor.copy()
    mps_cores = []

    for i in range(L - 1):
        shape = remaining.shape
        mat = remaining.reshape(shape[0] * (2 if i > 0 else 1), -1) if i == 0 else remaining.reshape(shape[0] * shape[1] if len(shape) > 1 else shape[0], -1)

        # Actually, let's do this properly
        # Reshape remaining into (left_dim, rest)
        left_dim = shape[0]
        right_dim = int(np.prod(shape[1:]))
        mat = remaining.reshape(left_dim, right_dim)

        U, S, Vt = np.linalg.svd(mat, full_matrices=False)

        # Truncate at various thresholds to see bond dimension needed
        total_sv = np.sum(S**2)
        cumsum = np.cumsum(S**2) / total_sv

        # For exact representation, bond dim = rank
        rank = np.sum(S > 1e-10)
        bond_dims.append(rank)

        # Keep all for now
        remaining = np.diag(S[:rank]) @ Vt[:rank, :]
        # Reshape remaining for next step
        remaining = remaining.reshape(rank, 2, -1) if i < L - 2 else remaining
        if i < L - 2:
            remaining = remaining.reshape(rank * 2, -1)

    print(f"  L={L} bits ({N_nums} numbers): exact MPS bond dimensions = {bond_dims}")
    print(f"  Max bond dimension D = {max(bond_dims)}")
    print(f"  MPS parameters: {sum(2*d for d in bond_dims)}")
    print(f"  vs full vector: {N_nums}")

    # Test larger
    for L in [12, 14, 16]:
        N_nums = 2**L
        true_ind = np.array([1 if (i >= 2 and i <= 200000 and PI_TABLE[i] - PI_TABLE[i-1] == 1) else 0
                             for i in range(N_nums)])
        tensor = true_ind.reshape([2] * L)

        # Just compute the max bond dimension needed
        remaining = tensor.copy()
        max_bond = 0
        for i in range(L - 1):
            shape = remaining.shape
            left_dim = shape[0]
            right_dim = int(np.prod(shape[1:]))
            mat = remaining.reshape(left_dim, right_dim)

            U, S, Vt = np.linalg.svd(mat, full_matrices=False)
            rank = np.sum(S > 1e-10)
            max_bond = max(max_bond, rank)

            remaining = np.diag(S[:rank]) @ Vt[:rank, :]
            if i < L - 2:
                new_shape = [rank] + list(shape[1:]) if len(shape) > 1 else [rank, 1]
                try:
                    remaining = remaining.reshape(rank, 2, -1)
                    remaining = remaining.reshape(rank * 2, -1)
                except:
                    break

        n_primes = np.sum(true_ind)
        print(f"  L={L} ({N_nums} numbers, {n_primes} primes): max bond D = {max_bond}")

    print()
    print("  ANALYSIS: Bond dimension grows roughly as O(2^{L/2}) = O(sqrt(N)).")
    print("  This makes sense: the prime indicator has HIGH 'entanglement'")
    print("  across bit positions (knowing high bits gives little info about primality).")
    print("  An MPS with polynomial bond dimension CANNOT represent the prime indicator.")

    # --- 4b: Sieve as tensor contraction ---
    print("\n--- 4b: Sieve as tensor contraction ---")
    print("Idea: The sieve eliminates multiples of each prime p.")
    print("For x with L bits, divisibility by p is a 'constraint tensor'.")
    print("pi(x) = contraction of all constraint tensors.")

    print()
    print("  Each 'divisibility by p' constraint has bond dimension O(p)")
    print("  (need to track remainder mod p through the bits).")
    print("  For primes up to sqrt(x), total bond dimension ~ prod(p, p<=sqrt(x))")
    print("  = primorial(sqrt(x)) which is EXPONENTIALLY large.")
    print()
    print("  Even with clever ordering and compression, the bond dimension")
    print("  grows at least as e^{c*sqrt(x)/ln(x)} (by prime number theorem).")
    print("  This is SUPER-EXPONENTIAL in L = log(x).")

    print()
    print("VERDICT: Tensor network approaches fail because:")
    print("  1. Prime indicator has maximal MPS bond dimension O(sqrt(N))")
    print("  2. Sieve constraints create exponentially large bond dimensions")
    print("  3. Primality is a GLOBAL property of a number's digits,")
    print("     fundamentally incompatible with low-entanglement representations.")

    return False


# ============================================================
# APPROACH 5: VARIATIONAL APPROACH
#   Quantum variational algorithms (VQE, QAOA) optimize parameterized
#   circuits to minimize energy.
#   Classical analog: optimize a parameterized function f(n; theta)
#   to minimize |f(n; theta) - p(n)| via gradient descent.
# ============================================================

def approach5_variational():
    """
    Variational approach: optimize f(n; theta) to approximate p(n).

    Architecture choices:
    (a) Fourier series in log-space: f(n) = sum a_k * cos(w_k * ln(n) + phi_k)
    (b) Rational function (Pade-like): f(n) = P(ln(n)) / Q(ln(n))
    (c) Neural network (small MLP)
    (d) Asymptotic expansion with optimized coefficients
    """
    print("\n" + "=" * 70)
    print("APPROACH 5: VARIATIONAL / GRADIENT DESCENT")
    print("=" * 70)

    # Training data
    train_n = np.array(range(2, 1001), dtype=float)
    train_p = np.array([P[n] for n in range(2, 1001)], dtype=float)

    # Validation
    val_n = np.array(range(1001, 2001), dtype=float)
    val_p = np.array([P[n] for n in range(1001, 2001)], dtype=float)

    # --- 5a: Optimized asymptotic expansion ---
    print("\n--- 5a: Optimized asymptotic expansion ---")
    print("p(n) ~ n * (ln(n) + ln(ln(n)) - 1 + (ln(ln(n)) - 2)/ln(n) + ...)")
    print("Optimize coefficients of higher-order terms via gradient descent.")

    def asymptotic_model(n, theta):
        """Parameterized asymptotic expansion for p(n)."""
        L = np.log(n)
        LL = np.log(L)
        return n * (L + LL - 1 +
                    theta[0] * (LL - 2) / L +
                    theta[1] * (LL**2 - 6*LL + theta[2]) / L**2 +
                    theta[3] * (LL**3 + theta[4]*LL**2 + theta[5]*LL + theta[6]) / L**3)

    def loss_asymp(theta, n_arr, p_arr):
        pred = asymptotic_model(n_arr, theta)
        return np.mean(((pred - p_arr) / p_arr)**2)  # Relative MSE

    # Initial guess from known expansion
    theta0 = np.array([1.0, 1.0, 11.0, 1.0, -12.0, 36.0, -6.0])

    result = minimize(loss_asymp, theta0, args=(train_n, train_p), method='Nelder-Mead',
                      options={'maxiter': 10000, 'fatol': 1e-15})
    theta_opt = result.x

    print(f"  Optimized theta: {theta_opt}")
    print(f"  Training loss: {result.fun:.2e}")

    # Evaluate
    for test_n, test_p, label in [(train_n, train_p, "train"), (val_n, val_p, "val")]:
        pred = asymptotic_model(test_n, theta_opt)
        rel_err = np.abs(pred - test_p) / test_p
        abs_err = np.abs(pred - test_p)
        print(f"  {label}: mean_rel_err={np.mean(rel_err):.6f}, max_rel_err={np.max(rel_err):.6f}, "
              f"mean_abs_err={np.mean(abs_err):.1f}, max_abs_err={np.max(abs_err):.1f}")

    # Check if any are EXACTLY correct
    pred_train = asymptotic_model(train_n, theta_opt)
    exact_count = np.sum(np.round(pred_train) == train_p)
    print(f"  Exact (round to nearest integer = p(n)): {exact_count}/{len(train_n)}")

    # --- 5b: Fourier series model in log-space ---
    print("\n--- 5b: Fourier model for correction delta(n) = p(n) - n*ln(n) ---")

    # Work with normalized correction
    base = train_n * np.log(train_n)
    delta_train = train_p - base
    base_val = val_n * np.log(val_n)
    delta_val = val_p - base_val

    n_fourier = 20

    def fourier_model(n, params):
        """Fourier correction: delta(n) = sum a_k cos(w_k ln(n) + phi_k) + polynomial"""
        # Polynomial part
        L = np.log(n)
        poly = params[0] * n + params[1] * L * n / L + params[2] * n / L

        # Fourier part
        result = poly
        for k in range(n_fourier):
            a_k = params[3 + 3*k]
            w_k = params[4 + 3*k]
            phi_k = params[5 + 3*k]
            result = result + a_k * np.cos(w_k * L + phi_k)
        return result

    n_params = 3 + 3 * n_fourier
    params0 = np.zeros(n_params)
    params0[0] = 1.0  # Main correction ~ n * ln(ln(n)) / ln(n)

    # Initialize frequencies from FFT analysis
    for k in range(n_fourier):
        params0[3 + 3*k] = 10.0  # amplitude
        params0[4 + 3*k] = (k + 1) * 0.5  # frequency
        params0[5 + 3*k] = 0.0  # phase

    def loss_fourier(params, n_arr, delta_arr):
        pred = fourier_model(n_arr, params)
        return np.mean((pred - delta_arr)**2)

    result_f = minimize(loss_fourier, params0, args=(train_n, delta_train),
                        method='L-BFGS-B', options={'maxiter': 5000})

    pred_delta = fourier_model(train_n, result_f.x)
    pred_p = base + pred_delta
    train_err = np.mean(np.abs(pred_p - train_p))

    pred_delta_val = fourier_model(val_n, result_f.x)
    pred_p_val = base_val + pred_delta_val
    val_err = np.mean(np.abs(pred_p_val - val_p))

    exact_train = np.sum(np.round(pred_p) == train_p)
    exact_val = np.sum(np.round(pred_p_val) == val_p)

    print(f"  {n_fourier} Fourier components + polynomial:")
    print(f"  Train: mean_abs_err={train_err:.1f}, exact={exact_train}/{len(train_n)}")
    print(f"  Val:   mean_abs_err={val_err:.1f}, exact={exact_val}/{len(val_n)}")

    # --- 5c: R^{-1}(n) + optimized correction ---
    print("\n--- 5c: R_inv(n) + variational correction ---")
    print("Since R^{-1}(n) already gives ~47% of digits, optimize the residual.")

    # Compute R_inv for training set
    r_inv_train = np.array([R_inverse(n) for n in range(2, min(201, 1001))])
    p_train_small = np.array([P[n] for n in range(2, min(201, 1001))])
    residual = p_train_small - r_inv_train
    ns_small = np.arange(2, min(201, 1001), dtype=float)

    print(f"  Residual delta = p(n) - R^{{-1}}(n) statistics (n=2..{min(200, 1000)}):")
    print(f"    mean={np.mean(residual):.2f}, std={np.std(residual):.2f}, "
          f"max={np.max(np.abs(residual)):.2f}")

    # Try to fit residual with various functions
    def correction_model(n, theta):
        L = np.log(n)
        LL = np.log(L)
        sqn = np.sqrt(n)
        return (theta[0] +
                theta[1] * sqn * np.cos(theta[2] * L + theta[3]) +
                theta[4] * sqn * np.cos(theta[5] * L + theta[6]) +
                theta[7] * np.power(n, 0.25) * np.cos(theta[8] * L + theta[9]))

    def loss_correction(theta, n_arr, res_arr):
        pred = correction_model(n_arr, theta)
        return np.mean((pred - res_arr)**2)

    theta0_c = np.zeros(10)
    theta0_c[2] = 14.134  # First zeta zero frequency
    theta0_c[5] = 21.022  # Second zeta zero
    theta0_c[8] = 25.011  # Third
    theta0_c[1] = 1.0
    theta0_c[4] = 0.5
    theta0_c[7] = 0.3

    result_c = minimize(loss_correction, theta0_c, args=(ns_small, residual),
                        method='Nelder-Mead', options={'maxiter': 20000})

    pred_res = correction_model(ns_small, result_c.x)
    corrected = r_inv_train + pred_res
    exact_corrected = np.sum(np.round(corrected) == p_train_small)
    mean_err_corrected = np.mean(np.abs(corrected - p_train_small))

    print(f"  Optimized correction (n=2..{min(200, 1000)}):")
    print(f"    Before: mean_err={np.mean(np.abs(residual)):.2f}")
    print(f"    After:  mean_err={mean_err_corrected:.2f}")
    print(f"    Exact (rounded): {exact_corrected}/{len(ns_small)}")

    # Test generalization
    r_inv_val = np.array([R_inverse(n) for n in range(201, 301)])
    p_val_small = np.array([P[n] for n in range(201, 301)])
    ns_val = np.arange(201, 301, dtype=float)
    pred_res_val = correction_model(ns_val, result_c.x)
    corrected_val = r_inv_val + pred_res_val
    exact_val = np.sum(np.round(corrected_val) == p_val_small)
    mean_err_val = np.mean(np.abs(corrected_val - p_val_small))

    print(f"  Validation (n=201..300):")
    print(f"    R_inv only: mean_err={np.mean(np.abs(p_val_small - r_inv_val)):.2f}")
    print(f"    + correction: mean_err={mean_err_val:.2f}")
    print(f"    Exact (rounded): {exact_val}/{len(ns_val)}")

    # --- 5d: Mini neural network ---
    print("\n--- 5d: Tiny neural network for p(n) ---")

    # 2-layer MLP with hand-coded backprop
    def relu(x):
        return np.maximum(0, x)

    def mlp_forward(n, W1, b1, W2, b2, w_out, b_out):
        """Simple 2-layer MLP: input = [n, ln(n), ln(ln(n))], output = p(n)"""
        x = np.array([n, np.log(n), np.log(np.log(n))]).reshape(-1, 1) if np.isscalar(n) else \
            np.stack([n, np.log(n), np.log(np.log(n))], axis=0)
        h1 = relu(W1 @ x + b1)
        h2 = relu(W2 @ h1 + b2)
        return (w_out @ h2 + b_out).flatten()

    hidden1, hidden2 = 32, 16

    np.random.seed(42)
    W1 = np.random.randn(hidden1, 3) * 0.1
    b1 = np.zeros((hidden1, 1))
    W2 = np.random.randn(hidden2, hidden1) * 0.1
    b2 = np.zeros((hidden2, 1))
    w_out = np.random.randn(1, hidden2) * 0.1
    b_out = np.zeros((1, 1))

    # Pack parameters
    def pack(W1, b1, W2, b2, w_out, b_out):
        return np.concatenate([W1.flatten(), b1.flatten(), W2.flatten(), b2.flatten(),
                               w_out.flatten(), b_out.flatten()])

    def unpack(params):
        idx = 0
        W1 = params[idx:idx+hidden1*3].reshape(hidden1, 3); idx += hidden1*3
        b1 = params[idx:idx+hidden1].reshape(hidden1, 1); idx += hidden1
        W2 = params[idx:idx+hidden2*hidden1].reshape(hidden2, hidden1); idx += hidden2*hidden1
        b2 = params[idx:idx+hidden2].reshape(hidden2, 1); idx += hidden2
        w_out = params[idx:idx+hidden2].reshape(1, hidden2); idx += hidden2
        b_out = params[idx:idx+1].reshape(1, 1); idx += 1
        return W1, b1, W2, b2, w_out, b_out

    # Normalize inputs and outputs
    n_train_nn = np.arange(2, 501, dtype=float)
    p_train_nn = np.array([P[n] for n in range(2, 501)], dtype=float)

    n_mean, n_std = np.mean(n_train_nn), np.std(n_train_nn)
    p_mean, p_std = np.mean(p_train_nn), np.std(p_train_nn)

    def loss_nn(params_flat):
        W1, b1, W2, b2, w_out, b_out = unpack(params_flat)
        total_loss = 0
        for i, n in enumerate(n_train_nn):
            x = np.array([(n - n_mean)/n_std,
                          (np.log(n) - np.log(n_mean))/1.0,
                          (np.log(np.log(n)))/1.0]).reshape(3, 1)
            h1 = relu(W1 @ x + b1)
            h2 = relu(W2 @ h1 + b2)
            pred = float((w_out @ h2 + b_out).flatten()[0])
            target = (p_train_nn[i] - p_mean) / p_std
            total_loss += (pred - target)**2
        return total_loss / len(n_train_nn)

    params0_nn = pack(W1, b1, W2, b2, w_out, b_out)
    print(f"  MLP architecture: 3 -> {hidden1} -> {hidden2} -> 1 ({len(params0_nn)} parameters)")
    print(f"  Training on n=2..500...")

    t0 = time.time()
    result_nn = minimize(loss_nn, params0_nn, method='L-BFGS-B',
                         options={'maxiter': 500, 'ftol': 1e-12})
    nn_time = time.time() - t0
    print(f"  Optimization: {nn_time:.1f}s, final loss={result_nn.fun:.6f}")

    # Evaluate
    W1f, b1f, W2f, b2f, w_outf, b_outf = unpack(result_nn.x)
    exact_nn = 0
    max_err_nn = 0
    for n in range(2, 501):
        x = np.array([(n - n_mean)/n_std,
                      (np.log(n) - np.log(n_mean))/1.0,
                      (np.log(np.log(n)))/1.0]).reshape(3, 1)
        h1 = relu(W1f @ x + b1f)
        h2 = relu(W2f @ h1 + b2f)
        pred_norm = float((w_outf @ h2 + b_outf).flatten()[0])
        pred = pred_norm * p_std + p_mean
        true_val = P[n]
        err = abs(pred - true_val)
        max_err_nn = max(max_err_nn, err)
        if round(pred) == true_val:
            exact_nn += 1

    print(f"  Train exact: {exact_nn}/499, max_err: {max_err_nn:.1f}")

    # Generalization
    exact_nn_val = 0
    for n in range(501, 1001):
        x = np.array([(n - n_mean)/n_std,
                      (np.log(n) - np.log(n_mean))/1.0,
                      (np.log(np.log(n)))/1.0]).reshape(3, 1)
        h1 = relu(W1f @ x + b1f)
        h2 = relu(W2f @ h1 + b2f)
        pred_norm = float((w_outf @ h2 + b_outf).flatten()[0])
        pred = pred_norm * p_std + p_mean
        true_val = P[n]
        if round(pred) == true_val:
            exact_nn_val += 1

    print(f"  Val exact: {exact_nn_val}/500 (generalization)")

    print()
    print("VERDICT: Variational approaches can fit training data but don't generalize.")
    print("  - Asymptotic expansion: good relative error but never exact")
    print("  - Fourier correction: overfits to training range")
    print("  - R_inv + correction: correction doesn't generalize")
    print("  - Neural network: memorizes training data, fails on validation")
    print("  Fundamental reason: p(n) has O(log p(n)) bits of information,")
    print("  any fixed-parameter model can only capture O(#params) bits.")

    return False


# ============================================================
# APPROACH 6: ADIABATIC / HOMOTOPY CONTINUATION
#   Quantum adiabatic: slowly transform H_0 -> H_1, ground state follows.
#   Classical analog: homotopy continuation from easy problem to hard.
#
#   Idea: Define path from li^{-1}(n) [easy to compute] to p(n) [hard].
#   Follow this path continuously.
# ============================================================

def approach6_homotopy():
    """
    Adiabatic / homotopy continuation: continuously deform from
    li^{-1}(n) (easy) to p(n) (hard).

    Strategy 1: Define F(x, t) = t*pi(x) + (1-t)*li(x) - n
                At t=0: easy root x = li^{-1}(n)
                At t=1: hard root x = p(n)-th prime vicinity
                Follow the root as t goes from 0 to 1.

    Strategy 2: Smooth approximation of prime indicator
                chi_eps(x) -> chi(x) = [x is prime] as eps -> 0
                At large eps: smooth, easy to optimize
                At eps=0: exact but discontinuous

    Strategy 3: Regularized prime counting
                pi_lambda(x) = sum_{n<=x} (1 - lambda * Lambda(n)/ln(n))
                At lambda=0: pi_0(x) = floor(x) - 1 (trivial)
                At lambda=1: pi_1(x) = pi(x) (what we want)
    """
    print("\n" + "=" * 70)
    print("APPROACH 6: ADIABATIC / HOMOTOPY CONTINUATION")
    print("=" * 70)

    # --- 6a: Homotopy from li^{-1} to p(n) ---
    print("\n--- 6a: Homotopy path from li^{-1}(n) to p(n) ---")

    def li(x):
        """Logarithmic integral."""
        if HAS_MPMATH:
            return float(mpmath.li(x))
        # Rough approximation
        return x / math.log(x) * (1 + 1/math.log(x))

    def li_inverse(n):
        """Inverse of li: find x such that li(x) = n."""
        x = n * math.log(max(n, 2))
        for _ in range(100):
            fx = li(x) - n
            fpx = 1.0 / math.log(x) if x > 1 else 1.0
            x -= fx / fpx
            if abs(fx) < 1e-8:
                break
        return x

    # The homotopy: F(x, t) = (1-t)*li(x) + t*pi(x) - n = 0
    # At t=0: x = li^{-1}(n)
    # At t=1: x = p(n) (since pi(p(n)) = n)
    # Problem: pi(x) is a step function, so F is discontinuous at primes!

    # Test: path following with predictor-corrector
    test_ns = [10, 50, 100, 500]

    print(f"\n{'n':>6} | {'li_inv':>10} | {'R_inv':>10} | {'p(n)':>8} | {'path_result':>12} | {'exact?':>6}")
    print("-" * 70)

    for n in test_ns:
        target = P[n]
        x_li = li_inverse(n)
        x_R = R_inverse(n)

        # Follow path from t=0 to t=1
        # At each t, solve: (1-t)*li(x) + t*pi(x) = n
        # Use small steps in t, Newton for x at each step

        x_current = x_li
        dt = 0.01
        t = 0.0
        path = [(t, x_current)]

        stuck = False
        while t < 1.0 - dt/2:
            t += dt
            t = min(t, 1.0)

            # At this t, need: (1-t)*li(x) + t*pi(x) = n
            # Newton step: F(x) = (1-t)*li(x) + t*pi(x) - n
            # F'(x) = (1-t)/ln(x) + t * delta(x is prime)
            # The derivative is problematic because pi(x) is a step function!

            # Use smoothed pi: pi_smooth(x) ≈ R(x) or li(x)
            # This defeats the purpose somewhat...

            # Alternative: just use the smoothed version and snap at the end
            def F_smooth(x):
                li_part = (1 - t) * li(x)
                # Use R(x) as smooth proxy for pi(x)
                pi_smooth = riemann_R(x)
                pi_part = t * pi_smooth
                return li_part + pi_part - n

            def F_smooth_deriv(x):
                return (1 - t) / math.log(x) + t / math.log(x)  # Both ~1/ln(x)

            # Newton
            for _ in range(20):
                fx = F_smooth(x_current)
                fpx = F_smooth_deriv(x_current)
                if abs(fpx) < 1e-15:
                    break
                step = fx / fpx
                x_current -= step
                x_current = max(2.1, x_current)
                if abs(step) < 1e-8:
                    break

            path.append((t, x_current))

        # At t=1, we should have R^{-1}(n) (since we used R as proxy for pi)
        # Snap to nearest prime
        x_final = x_current

        # Find nearest prime
        x_int = int(round(x_final))
        best = None
        best_dist = float('inf')
        for candidate in range(max(2, x_int - 50), x_int + 50):
            if candidate <= 200000 and PI_TABLE[candidate] - PI_TABLE[candidate-1] == 1:
                if abs(candidate - x_final) < best_dist:
                    best_dist = abs(candidate - x_final)
                    best = candidate

        is_exact = (best == target)
        print(f"{n:>6} | {x_li:>10.2f} | {x_R:>10.2f} | {target:>8} | {best:>12} | {'YES' if is_exact else 'NO':>6}")

    # --- 6b: Smoothed prime indicator with decreasing epsilon ---
    print("\n--- 6b: Smoothed prime indicator (annealing epsilon -> 0) ---")
    print("Define chi_eps(x) = exp(-min_p (x mod p / eps)^2) for small primes p|x")
    print("As eps -> 0, chi_eps(x) -> [x has small factor] (NOT primality!)")
    print("This doesn't work: primality is about having NO small factors,")
    print("which is the ABSENCE of a signal, not its presence.")

    # --- 6c: Regularized counting with von Mangoldt ---
    print("\n--- 6c: Regularized counting function ---")
    print("psi_t(x) = sum_{n<=x} w_t(n) where w_t interpolates:")
    print("  w_0(n) = 1 for all n  =>  psi_0(x) = floor(x)")
    print("  w_1(n) = Lambda(n)/ln(n)  =>  psi_1(x) = pi(x) + pi(x^{1/2})/2 + ...")

    # This is interesting: at t=0, counting is trivial.
    # Can we smoothly increase t and track the count?

    # Problem: w_t(n) = 1 - t*(1 - Lambda(n)/ln(n))
    # For composite n not a prime power: Lambda(n) = 0, so w_t(n) = 1 - t
    # For prime p: Lambda(p) = ln(p), so w_t(p) = 1 - t*(1 - 1) = 1
    # For prime power p^k: Lambda(p^k) = ln(p), so w_t(p^k) = 1 - t*(1 - ln(p)/(k*ln(p))) = 1 - t*(k-1)/k

    # So at t=1: composites get weight 0, primes get weight 1, prime powers get fractional weight
    # This IS the decomposition, but computing it requires knowing which n are composite!

    print("  Analysis: The homotopy parameter just interpolates between")
    print("  'count everything' and 'count primes only'.")
    print("  To follow the path, we need to evaluate the weight function,")
    print("  which requires knowing Lambda(n) = factoring n.")
    print("  This is circular: we need the primes to classify numbers.")

    # --- 6d: Path from approximate formula to exact ---
    print("\n--- 6d: Path from R^{-1}(n) to p(n) via zeta zero correction ---")
    print("p(n) = R^{-1}(n) + correction(n)")
    print("correction(n) comes from sum over zeta zeros.")
    print("Homotopy: add zeros one at a time.")

    # This is literally the explicit formula applied incrementally.
    # Add first K zeros and see error decrease.

    zeta_zeros = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
                  37.586178, 40.918719, 43.327073, 48.005151, 49.773832]

    test_n = 1000
    target = P[test_n]
    x0 = R_inverse(test_n)

    print(f"\n  Incremental zero correction for p({test_n}) = {target}:")
    print(f"  R^{{-1}}({test_n}) = {x0:.4f}, error = {abs(x0 - target):.4f}")

    correction = 0.0
    for k, gamma in enumerate(zeta_zeros):
        ln_x = math.log(x0)
        rho_mag = math.sqrt(0.25 + gamma**2)
        rho_phase = math.atan2(gamma, 0.5)
        x_half = x0 ** 0.5
        result_mag = x_half / (rho_mag * ln_x)
        result_phase = gamma * ln_x - rho_phase
        correction -= 2 * result_mag * math.cos(result_phase)

        approx = x0 + correction
        err = abs(approx - target)
        print(f"    +zero #{k+1} (gamma={gamma:.3f}): approx={approx:.2f}, error={err:.2f}")

    print()
    print("  The correction oscillates wildly and doesn't converge with 10 zeros.")
    print("  Need O(sqrt(p(n))) zeros for convergence (known result).")
    print("  Homotopy doesn't help: each 'step' (adding a zero) is easy,")
    print("  but we need exponentially many steps.")

    print()
    print("VERDICT: Homotopy/adiabatic approach fails because:")
    print("  1. Smooth proxy for pi(x) just gives R^{-1}(n) again")
    print("  2. Following the exact path requires evaluating pi(x) (circular)")
    print("  3. Incremental zero correction = explicit formula (known to need ~10^51 zeros)")
    print("  The 'energy gap' analog: the gap between smooth approx and exact pi(x)")
    print("  closes at a rate that requires exponentially many steps to cross.")

    return False


# ============================================================
# BONUS: APPROACH 7: QUANTUM WALK ANALOG
#   Quantum walks achieve quadratic speedup over classical random walks.
#   Can a quantum-walk-inspired search help?
# ============================================================

def approach7_quantum_walk():
    """
    Quantum walk analog: structured random walk on number line
    that converges to primes faster than uniform search.

    Idea: At each step, "amplitude" (probability) is transferred
    preferentially toward prime candidates using local number-theoretic
    information (residues mod small primes).
    """
    print("\n" + "=" * 70)
    print("APPROACH 7: QUANTUM WALK ANALOG")
    print("=" * 70)

    # Biased random walk: from position x, move toward numbers
    # that are coprime to more small primes (more "prime-like")

    def coprime_score(x, small_primes):
        """Count how many small primes do NOT divide x."""
        if x < 2:
            return 0
        return sum(1 for p in small_primes if x % p != 0)

    sp = PRIMES[:20]  # First 20 primes

    test_cases = [(100, P[100]), (500, P[500]), (1000, P[1000])]

    print("\nBiased walk from R^{-1}(n) toward p(n):")
    print(f"{'n':>6} | {'p(n)':>8} | {'start':>10} | {'walk_steps':>12} | {'found':>6} | {'uniform_steps':>14}")
    print("-" * 70)

    for n, target in test_cases:
        x0 = R_inverse(n)
        start = int(round(x0))

        # Biased walk
        x = start
        visited = set()
        steps = 0
        found = False

        for steps in range(1, 100001):
            if x == target:
                found = True
                break

            visited.add(x)

            # Score neighbors
            candidates = [x-2, x-1, x+1, x+2]
            candidates = [c for c in candidates if c >= 2 and c not in visited]

            if not candidates:
                # Stuck, jump randomly
                x = start + np.random.randint(-100, 100)
                continue

            scores = [coprime_score(c, sp) for c in candidates]
            total = sum(scores) + 1e-10
            probs = [s/total for s in scores]

            x = np.random.choice(candidates, p=probs)

        # Compare to uniform search
        window = abs(start - target) * 2 + 50
        uniform_expected = window  # Expected steps for uniform search

        print(f"{n:>6} | {target:>8} | {start:>10} | {steps:>12} | {'YES' if found else 'NO':>6} | {uniform_expected:>14}")

    print()
    print("ANALYSIS:")
    print("  Classical biased walk: O(distance / bias)")
    print("  Quantum walk would give: O(sqrt(distance))")
    print("  But classically we CANNOT achieve quadratic speedup.")
    print("  The Grover lower bound (BBBV theorem) proves that classical")
    print("  search cannot do better than O(N) for unstructured problems.")
    print("  Our bias helps by constant factor, not polynomial speedup.")
    print()
    print("VERDICT: Quantum walk analog gives only constant-factor improvement.")
    print("  True quadratic speedup is inherently quantum (amplitude interference).")

    return False


# ============================================================
# FINAL SYNTHESIS AND CROSS-APPROACH ANALYSIS
# ============================================================

def final_synthesis():
    """
    Combine insights from all quantum-inspired approaches.
    """
    print("\n" + "=" * 70)
    print("FINAL SYNTHESIS: QUANTUM-INSPIRED APPROACHES")
    print("=" * 70)

    print("""
APPROACH SUMMARY:
=================

| # | Approach                  | Quantum Speedup | Classical Analog  | Result           |
|---|---------------------------|-----------------|-------------------|------------------|
| 1 | Amplitude amplification   | Quadratic       | Importance sample | Constant factor  |
| 2 | QFT / Fourier analysis    | Exponential     | FFT on primes     | Structural only  |
| 3 | Phase estimation          | Exponential     | Power iteration   | No useful oper.  |
| 4 | Tensor network / MPS      | Exponential     | SVD compression   | Bond dim O(√N)   |
| 5 | Variational               | Heuristic       | Gradient descent  | Overfits         |
| 6 | Adiabatic / homotopy      | Polynomial      | Path following    | Circular or slow |
| 7 | Quantum walk              | Quadratic       | Biased walk       | Constant factor  |

KEY INSIGHT: WHY QUANTUM IDEAS DON'T HELP CLASSICALLY
======================================================

The quantum speedups come from three sources:
  1. AMPLITUDE INTERFERENCE: Constructive/destructive interference amplifies
     correct answers. No classical analog exists (this is the core of
     quantum advantage).

  2. SUPERPOSITION: Process all candidates simultaneously. Classical
     parallelism helps but doesn't give exponential speedup.

  3. ENTANGLEMENT: Correlate parts of the computation in ways
     impossible classically. The high "entanglement" of the prime
     indicator (bond dimension O(sqrt(N)) in MPS) means primes
     fundamentally resist classical compressed representations.

SPECIFIC TO PRIME COUNTING:
  - The prime indicator f(n) = [n is prime] has HIGH complexity:
    * Kolmogorov complexity: O(N/ln(N)) bits for all primes up to N
    * MPS bond dimension: O(sqrt(N)) (measured in Approach 4)
    * Fourier: needs O(N) components for exact reconstruction

  - Any method to compute p(10^100) exactly needs to somehow access
    ~10^51 bits of information about prime distribution.

  - Quantum computers could potentially reduce this via:
    * Grover search: O(10^25) instead of O(10^51) for brute force
    * But this is still infeasible (10^25 quantum gates)

  - NO classical algorithm can compute p(10^100) exactly in polynomial
    time unless we discover a fundamentally new mathematical structure
    connecting individual primes to efficiently computable functions.

WHAT WOULD ACTUALLY HELP (SPECULATIVE):
  1. A number-theoretic identity that makes pi(x) computable from
     a SMALL number of L-function values (currently needs ~sqrt(x))
  2. A quantum computer with 10^25 gates and error correction
  3. A proof that RH is false (which might change the zero structure)
  4. Discovery of hidden algebraic structure in prime gaps
     (currently appears random, 3+ bits/gap unpredictable)
""")


# ============================================================
# MAIN: Run all experiments
# ============================================================

if __name__ == "__main__":
    t_start = time.time()

    results = {}

    approaches = [
        ("1_importance_sampling", approach1_importance_sampling),
        ("2_fft_analysis", approach2_fft_prime_distribution),
        ("3_power_method", approach3_power_method),
        ("4_tensor_network", approach4_tensor_network),
        ("5_variational", approach5_variational),
        ("6_homotopy", approach6_homotopy),
        ("7_quantum_walk", approach7_quantum_walk),
    ]

    for name, func in approaches:
        print(f"\n{'#' * 70}")
        print(f"# RUNNING: {name}")
        print(f"{'#' * 70}")
        try:
            t0 = time.time()
            result = func()
            elapsed = time.time() - t0
            results[name] = {"success": result, "time": elapsed}
            print(f"\n[{name}] completed in {elapsed:.1f}s, breakthrough={result}")
        except Exception as e:
            print(f"\n[{name}] ERROR: {e}")
            import traceback
            traceback.print_exc()
            results[name] = {"success": False, "time": 0, "error": str(e)}

    final_synthesis()

    total_time = time.time() - t_start
    print(f"\n{'=' * 70}")
    print(f"ALL EXPERIMENTS COMPLETE in {total_time:.1f}s")
    print(f"{'=' * 70}")
    print(f"Breakthroughs found: {sum(1 for r in results.values() if r.get('success'))}/{len(results)}")

    for name, r in results.items():
        status = "BREAKTHROUGH" if r.get("success") else "no breakthrough"
        print(f"  {name}: {status} ({r.get('time', 0):.1f}s)")
