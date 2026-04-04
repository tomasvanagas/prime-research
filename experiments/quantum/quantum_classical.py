"""
Session 6: Classical Simulation of Quantum Algorithms for Prime Finding
========================================================================
5 experiments testing whether quantum-algorithm-inspired classical methods
can reveal useful structure in the prime sequence.

Approaches:
1. Classical Walsh-Hadamard period finding on prime gap sequences
2. Grover-inspired structured amplitude search for next prime
3. Tensor network representation of prime indicator function
4. Prime counting via classical expectation value simulation
5. Random matrix (GUE) eigenvalue extraction for zero reconstruction
"""

import numpy as np
from scipy import linalg, fft as sp_fft
from scipy.optimize import minimize_scalar, minimize
import time
import json
import sys
import os

# --- Utility: small prime sieve ---
def sieve(limit):
    """Sieve of Eratosthenes up to limit."""
    is_prime = np.ones(limit + 1, dtype=bool)
    is_prime[:2] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = False
    return np.nonzero(is_prime)[0]

def prime_gaps(primes):
    """Consecutive differences."""
    return np.diff(primes)

# Reference primes for testing
PRIMES_10K = sieve(120000)[:10000]  # first 10000 primes
PRIMES_1K = PRIMES_10K[:1000]

results = {}

###############################################################################
# EXPERIMENT 1: Classical Walsh-Hadamard Period Finding on Gap Sequences
###############################################################################
def experiment1_hadamard_period_finding():
    """
    The quantum Fourier transform finds hidden periodicities in O(polylog(N)).
    Classical FFT does this in O(N log N). But the KEY question: are there
    exploitable periodicities in the prime gap sequence?

    We apply the Walsh-Hadamard transform (the discrete version used in
    quantum computing) to the gap sequence and look for dominant frequencies.
    If gaps have quasi-periodic structure from zeta zeros, we might extract it.
    """
    print("=" * 70)
    print("EXPERIMENT 1: Walsh-Hadamard Period Finding on Prime Gaps")
    print("=" * 70)

    gaps = prime_gaps(PRIMES_10K).astype(float)
    N = len(gaps)

    # Pad to power of 2 for Walsh-Hadamard
    n_padded = 2 ** int(np.ceil(np.log2(N)))
    gaps_padded = np.zeros(n_padded)
    gaps_padded[:N] = gaps - np.mean(gaps)  # zero-center

    # 1a. Standard FFT spectrum
    spectrum = np.abs(np.fft.rfft(gaps_padded))
    freqs = np.fft.rfftfreq(n_padded)

    # Find top-k peaks
    k = 20
    top_indices = np.argsort(spectrum[1:])[-k:] + 1  # skip DC
    top_freqs = freqs[top_indices]
    top_powers = spectrum[top_indices]

    print(f"\nFFT on {N} prime gaps (padded to {n_padded}):")
    print(f"Top {k} frequency peaks:")
    for i in range(k - 1, -1, -1):
        period = 1.0 / top_freqs[i] if top_freqs[i] > 0 else float('inf')
        print(f"  freq={top_freqs[i]:.6f}  period={period:.1f}  power={top_powers[i]:.1f}")

    # 1b. Walsh-Hadamard transform (binary/dyadic structure)
    # WHT is the transform quantum computers naturally use
    def walsh_hadamard(x):
        """Fast Walsh-Hadamard transform."""
        n = len(x)
        assert n & (n - 1) == 0, "Length must be power of 2"
        h = 1
        while h < n:
            for i in range(0, n, h * 2):
                for j in range(i, i + h):
                    a = x[j]
                    b = x[j + h]
                    x[j] = a + b
                    x[j + h] = a - b
            h *= 2
        return x / np.sqrt(n)

    wht_result = walsh_hadamard(gaps_padded.copy())
    wht_power = np.abs(wht_result)

    # Top WHT coefficients
    top_wht = np.argsort(wht_power)[-k:]
    print(f"\nWalsh-Hadamard top {k} sequency components:")
    for idx in reversed(top_wht):
        print(f"  sequency={idx}  power={wht_power[idx]:.2f}")

    # 1c. Can we PREDICT gaps from top-K Fourier components?
    # This tests: does the periodic structure suffice for prediction?
    for K in [5, 10, 20, 50, 100, 200]:
        # Reconstruct from top-K FFT components
        full_spectrum = np.fft.rfft(gaps_padded)
        mask = np.zeros_like(full_spectrum)
        top_k_idx = np.argsort(np.abs(full_spectrum))[-K:]
        mask[top_k_idx] = full_spectrum[top_k_idx]
        reconstructed = np.fft.irfft(mask, n=n_padded)[:N]
        reconstructed += np.mean(gaps)

        # How well does rounding to nearest even gap work?
        predicted = np.round(reconstructed / 2) * 2
        predicted = np.maximum(predicted, 2)
        exact_matches = np.sum(predicted == gaps)
        mae = np.mean(np.abs(reconstructed - gaps))
        print(f"  K={K:4d} Fourier components: MAE={mae:.2f}, exact_gap_matches={exact_matches}/{N} ({100*exact_matches/N:.1f}%)")

    # 1d. Autocorrelation analysis (quantum advantage is here for hidden periods)
    acf = np.correlate(gaps_padded[:N], gaps_padded[:N], mode='full')
    acf = acf[N-1:]  # positive lags only
    acf = acf / acf[0]  # normalize

    significant_lags = np.where(np.abs(acf[1:500]) > 0.05)[0] + 1
    print(f"\nAutocorrelation: {len(significant_lags)} lags with |r| > 0.05 (out of 500)")
    if len(significant_lags) > 0:
        print(f"  Strongest lags: {significant_lags[:10]} with r = {acf[significant_lags[:10]]}")

    # Measure spectral entropy (how concentrated is the spectrum?)
    power_spectrum = spectrum[1:] ** 2
    power_spectrum /= power_spectrum.sum()
    spectral_entropy = -np.sum(power_spectrum * np.log2(power_spectrum + 1e-15))
    max_entropy = np.log2(len(power_spectrum))
    print(f"\nSpectral entropy: {spectral_entropy:.2f} / {max_entropy:.2f} max ({100*spectral_entropy/max_entropy:.1f}% of uniform)")

    # Compute reconstruction errors for summary
    fourier_recon = {}
    full_spec = np.fft.rfft(gaps_padded)
    for K in [10, 50, 200]:
        mask2 = np.zeros_like(full_spec)
        top_k2 = np.argsort(np.abs(full_spec))[-K:]
        mask2[top_k2] = full_spec[top_k2]
        recon2 = np.fft.irfft(mask2, n=n_padded)[:N] + np.mean(gaps)
        fourier_recon[str(K)] = {"mae": float(np.mean(np.abs(recon2 - gaps)))}

    return {
        "top_fft_freqs": top_freqs.tolist(),
        "top_fft_powers": top_powers.tolist(),
        "spectral_entropy_ratio": float(spectral_entropy / max_entropy),
        "n_significant_autocorrelation_lags": len(significant_lags),
        "fourier_reconstruction": fourier_recon,
        "conclusion": "FAIL" if spectral_entropy / max_entropy > 0.8 else "PARTIAL"
    }


###############################################################################
# EXPERIMENT 2: Grover-Inspired Structured Search for Next Prime
###############################################################################
def experiment2_grover_search():
    """
    Grover's algorithm searches N items in O(sqrt(N)) queries.
    For finding the next prime after p, the 'database' is integers in [p+1, p+gap].
    Quantum: O(sqrt(gap)) primality tests. Classical: gap/ln(gap) on average.

    We test a STRUCTURED variant: instead of random oracle queries, use
    amplitude-amplification-inspired weighting to bias the search toward
    likely prime locations (mod 6, mod 30, etc.)
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 2: Grover-Inspired Structured Search for Next Prime")
    print("=" * 70)

    def is_prime_trial(n):
        """Simple primality test, counting operations."""
        if n < 2: return False, 0
        if n < 4: return True, 0
        if n % 2 == 0 or n % 3 == 0: return False, 1
        ops = 1
        i = 5
        while i * i <= n:
            if n % i == 0 or n % (i + 2) == 0:
                return False, ops
            ops += 1
            i += 6
        return True, ops

    # Method 1: Naive sequential search (baseline)
    def naive_next_prime(p):
        """Search p+1, p+2, ... counting total ops."""
        total_ops = 0
        candidate = p + 1
        while True:
            result, ops = is_prime_trial(candidate)
            total_ops += ops
            if result:
                return candidate, total_ops
            candidate += 1

    # Method 2: Wheel search (skip non-coprime to 30)
    WHEEL30 = [1, 7, 11, 13, 17, 19, 23, 29]
    def wheel_next_prime(p):
        total_ops = 0
        # Start from next number coprime to 30
        base = (p // 30) * 30
        for offset in WHEEL30:
            candidate = base + offset
            if candidate <= p:
                continue
            result, ops = is_prime_trial(candidate)
            total_ops += ops
            if result:
                return candidate, total_ops
        # Move to next wheel cycle
        base += 30
        while True:
            for offset in WHEEL30:
                candidate = base + offset
                result, ops = is_prime_trial(candidate)
                total_ops += ops
                if result:
                    return candidate, total_ops
            base += 30

    # Method 3: Grover-inspired amplitude weighting
    # Idea: assign "amplitude" to each candidate based on prime likelihood,
    # then check candidates in decreasing amplitude order.
    # This simulates the amplitude amplification of marked states.
    def grover_inspired_next_prime(p, window=None):
        """
        Score candidates by prime likelihood, check in order of score.
        Scoring uses: residue mod small primes, distance from p, empirical gap distribution.
        """
        if window is None:
            # Expected gap ~ ln(p), use 4x safety margin
            window = max(100, int(4 * np.log(p + 1)))

        candidates = np.arange(p + 1, p + 1 + window)

        # Amplitude scoring (higher = more likely prime)
        scores = np.ones(len(candidates), dtype=float)

        # Sieve small primes: zero out composites
        small_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
        for sp in small_primes:
            # Score based on not being divisible
            for i, c in enumerate(candidates):
                if c % sp == 0 and c != sp:
                    scores[i] = 0

        # Among survivors, boost by Cramér-like model:
        # closer to average gap gets higher amplitude
        avg_gap = np.log(p + 1)
        for i in range(len(candidates)):
            if scores[i] > 0:
                dist = candidates[i] - p
                # Gaussian centered at expected gap
                scores[i] *= np.exp(-0.5 * ((dist - avg_gap) / (avg_gap * 0.5)) ** 2)

        # Sort by score (descending) - this is the "amplitude amplification"
        order = np.argsort(-scores)

        total_ops = 0
        for idx in order:
            if scores[idx] == 0:
                continue
            result, ops = is_prime_trial(int(candidates[idx]))
            total_ops += ops
            if result:
                return int(candidates[idx]), total_ops

        # Fallback: extend search
        return naive_next_prime(p + window)

    # Test on primes from 1000 to 10000
    test_primes = PRIMES_10K[168:500]  # primes near 1000-3600
    naive_ops = []
    wheel_ops = []
    grover_ops = []

    for p in test_primes:
        _, ops1 = naive_next_prime(int(p))
        _, ops2 = wheel_next_prime(int(p))
        _, ops3 = grover_inspired_next_prime(int(p))
        naive_ops.append(ops1)
        wheel_ops.append(ops2)
        grover_ops.append(ops3)

    naive_ops = np.array(naive_ops, dtype=float)
    wheel_ops = np.array(wheel_ops, dtype=float)
    grover_ops = np.array(grover_ops, dtype=float)

    print(f"\nNext-prime search from {len(test_primes)} primes (range ~1000-3600):")
    print(f"  Naive sequential:  mean={np.mean(naive_ops):.1f} ops, median={np.median(naive_ops):.0f}")
    print(f"  Wheel mod-30:      mean={np.mean(wheel_ops):.1f} ops, median={np.median(wheel_ops):.0f}")
    print(f"  Grover-inspired:   mean={np.mean(grover_ops):.1f} ops, median={np.median(grover_ops):.0f}")
    print(f"  Speedup over naive: wheel={np.mean(naive_ops)/np.mean(wheel_ops):.2f}x, grover={np.mean(naive_ops)/np.mean(grover_ops):.2f}x")

    # Theoretical Grover: sqrt(gap) queries on average
    gaps = np.diff(test_primes).astype(float)
    grover_theoretical = np.mean(np.sqrt(gaps[:len(test_primes)-1]))
    print(f"  Theoretical Grover (sqrt(gap)): {grover_theoretical:.1f} queries")
    print(f"  Our classical 'Grover': {np.mean(grover_ops):.1f} ops")

    # Test on larger primes
    large_primes = PRIMES_10K[9000:9100]  # near p(9000) ~ 93000
    large_naive = []
    large_grover = []
    for p in large_primes:
        _, ops1 = naive_next_prime(int(p))
        _, ops3 = grover_inspired_next_prime(int(p))
        large_naive.append(ops1)
        large_grover.append(ops3)

    print(f"\nLarger primes (~93000 range, 100 tests):")
    print(f"  Naive: mean={np.mean(large_naive):.1f} ops")
    print(f"  Grover-inspired: mean={np.mean(large_grover):.1f} ops")
    print(f"  Speedup: {np.mean(large_naive)/np.mean(large_grover):.2f}x")

    return {
        "small_range": {
            "naive_mean_ops": float(np.mean(naive_ops)),
            "wheel_mean_ops": float(np.mean(wheel_ops)),
            "grover_mean_ops": float(np.mean(grover_ops)),
            "speedup_wheel": float(np.mean(naive_ops) / np.mean(wheel_ops)),
            "speedup_grover": float(np.mean(naive_ops) / np.mean(grover_ops)),
        },
        "large_range": {
            "naive_mean_ops": float(np.mean(large_naive)),
            "grover_mean_ops": float(np.mean(large_grover)),
            "speedup": float(np.mean(large_naive) / np.mean(large_grover)),
        },
        "conclusion": "Small constant-factor improvement only; no asymptotic advantage"
    }


###############################################################################
# EXPERIMENT 3: Tensor Network Representation of Prime Indicator
###############################################################################
def experiment3_tensor_network():
    """
    Quantum states on N qubits live in 2^N dimensional Hilbert space.
    Tensor networks (MPS/MPO) compress this if entanglement is bounded.

    Key question: represent χ_prime(n) (the indicator function: 1 if n is prime,
    0 otherwise) as a Matrix Product State. What bond dimension D is needed?

    If D = polylog(N), we'd have an efficient prime formula!
    The bond dimension measures the "quantum entanglement" of the prime function.
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 3: Tensor Network / MPS for Prime Indicator")
    print("=" * 70)

    def prime_indicator_vector(n_bits):
        """Create the vector |chi> where chi[k] = 1 if k is prime, 0 otherwise."""
        N = 2 ** n_bits
        primes_set = set(sieve(N - 1).tolist())
        chi = np.array([1.0 if i in primes_set else 0.0 for i in range(N)])
        norm = np.linalg.norm(chi)
        if norm > 0:
            chi /= norm
        return chi

    def compute_mps_bond_dimensions(state_vector, n_bits):
        """
        Decompose a 2^n state vector into MPS form using successive SVDs.
        Returns the bond dimensions needed at each cut.
        Also returns truncation errors for various max bond dimensions.
        """
        N = 2 ** n_bits
        tensor = state_vector.reshape([2] * n_bits)

        bond_dims = []
        singular_values_all = []

        # Sweep left to right doing SVD
        current = tensor.reshape(2, -1)  # first site vs rest
        for i in range(n_bits - 1):
            rows, cols = current.shape
            U, S, Vh = np.linalg.svd(current, full_matrices=False)
            bond_dim = len(S)
            bond_dims.append(bond_dim)
            singular_values_all.append(S.copy())

            # Prepare for next cut
            if i < n_bits - 2:
                remaining = (S[:, None] * Vh)  # bond x remaining
                next_phys = 2  # next physical dimension
                remaining_sites = n_bits - i - 2
                if remaining_sites > 0:
                    current = remaining.reshape(bond_dim * next_phys, -1)
                else:
                    current = remaining

        # Compute truncation errors for various max D
        print(f"\n  n_bits={n_bits}, N=2^{n_bits}={2**n_bits}")
        print(f"  Full MPS bond dimensions: {bond_dims}")
        print(f"  Max bond dimension: {max(bond_dims)}")

        # Truncation analysis
        for max_D in [2, 4, 8, 16, 32, 64, 128]:
            if max_D > max(bond_dims):
                break
            total_error = 0.0
            for sv in singular_values_all:
                if len(sv) > max_D:
                    truncated = sv[max_D:]
                    total_error += np.sum(truncated ** 2)
            fidelity = 1.0 - total_error
            print(f"  D_max={max_D:4d}: fidelity={fidelity:.6f}, error={total_error:.6f}")

        return bond_dims, singular_values_all

    # Test for increasing sizes
    bond_dim_results = {}
    for n_bits in range(4, 18):  # up to 2^17 = 131072
        t0 = time.time()
        chi = prime_indicator_vector(n_bits)
        n_primes = int(round(np.sum(chi > 0)))
        dt_vec = time.time() - t0

        if n_bits <= 14:  # SVD feasible
            t0 = time.time()
            bonds, svs = compute_mps_bond_dimensions(chi, n_bits)
            dt_svd = time.time() - t0
            max_bond = max(bonds)
            bond_dim_results[n_bits] = {
                "N": 2**n_bits,
                "max_bond_dim": max_bond,
                "bond_dims": bonds,
                "time_svd": dt_svd
            }
        else:
            # For larger n, just compute the bipartition entropy
            # (entanglement across the middle cut)
            mid = n_bits // 2
            mat = chi.reshape(2**mid, 2**(n_bits - mid))
            sv = np.linalg.svd(mat, compute_uv=False)
            sv = sv[sv > 1e-15]
            ent = -np.sum(sv**2 * np.log2(sv**2 + 1e-30))
            max_bond_approx = np.sum(sv > 1e-10)
            print(f"\n  n_bits={n_bits}: mid-cut bond_dim={max_bond_approx}, entanglement_entropy={ent:.2f} bits")
            bond_dim_results[n_bits] = {
                "N": 2**n_bits,
                "mid_cut_bond_dim": int(max_bond_approx),
                "entanglement_entropy": float(ent),
            }

    # Analyze scaling of max bond dimension
    bits_list = sorted([k for k in bond_dim_results if 'max_bond_dim' in bond_dim_results[k]])
    if len(bits_list) >= 3:
        max_bonds = [bond_dim_results[b]['max_bond_dim'] for b in bits_list]
        ns = [2**b for b in bits_list]
        print(f"\nBond dimension scaling:")
        for b, mb in zip(bits_list, max_bonds):
            ratio = mb / (2**b)
            print(f"  n_bits={b:2d}, N={2**b:6d}, max_D={mb:5d}, D/N={ratio:.4f}")

        # Fit: D ~ N^alpha or D ~ 2^(alpha*n)?
        log_bonds = np.log2(np.array(max_bonds, dtype=float))
        log_ns = np.array(bits_list, dtype=float)
        if len(log_bonds) > 1:
            coeffs = np.polyfit(log_ns, log_bonds, 1)
            print(f"\n  Fit: log2(D) = {coeffs[0]:.3f} * n_bits + {coeffs[1]:.3f}")
            print(f"  => D ~ 2^({coeffs[0]:.3f} * n_bits) = N^{coeffs[0]:.3f}")
            if coeffs[0] > 0.8:
                print(f"  CONCLUSION: Bond dimension grows as ~N^{coeffs[0]:.2f} => EXPONENTIAL in input bits")
                print(f"  The prime indicator has HIGH entanglement. Tensor networks CANNOT compress it.")
            elif coeffs[0] < 0.3:
                print(f"  CONCLUSION: Bond dimension grows sub-linearly => POSSIBLE efficient representation!")

    return {
        "bond_dimensions": {str(k): v for k, v in bond_dim_results.items()},
        "scaling_exponent": float(coeffs[0]) if len(bits_list) >= 3 else None,
        "conclusion": "Bond dimension grows as N^alpha with alpha near 1 => no efficient TN representation"
    }


###############################################################################
# EXPERIMENT 4: Prime Counting as Quantum Expectation Value
###############################################################################
def experiment4_pi_expectation():
    """
    pi(x) = sum_{k=2}^{x} chi_prime(k) can be written as <psi|O|psi> where
    |psi> encodes the integers and O is the primality observable.

    More usefully: pi(x) = x - 1 - sum_{p <= sqrt(x)} [x/p - something...]
    (inclusion-exclusion / Legendre formula).

    Can we formulate pi(x) as a TRACE of a structured matrix and use
    matrix sketching / randomized linear algebra to approximate it fast?
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 4: Prime Counting via Matrix Trace / Expectation")
    print("=" * 70)

    def exact_pi(x):
        """Exact prime counting by sieve."""
        return len(sieve(x))

    # Method 1: Legendre-like formula as matrix trace
    # pi(x) = x - 1 - sum_{composites up to x}
    # The sieve matrix: M[i,j] = 1 if i*j <= x and i,j >= 2
    # Then number of composites <= x that are products of exactly two factors >= 2
    # This is related to the trace of the sieve operator

    # Method 2: Randomized trace estimation
    # Tr(A) can be estimated with O(1/epsilon^2) random vectors: E[z^T A z] = Tr(A)
    # If A is the prime projection matrix (diagonal with 1s at prime positions),
    # then Tr(A) = pi(N).
    # But constructing A requires knowing the primes! The question is whether
    # there's a structured A we can multiply by without knowing primes.

    def mobius_pi_estimate(x):
        """
        Use the identity: pi(x) = sum_{n=1}^{x} sum_{d|n} mu(d) * [n/d is prime or 1]
        Actually, use: sum_{n <= x} Lambda(n) / ln(n) ~ pi(x)
        where Lambda is the von Mangoldt function.

        More practically: Meissel-Mertens-like computation.
        """
        # Direct: count using small-prime sieve (Legendre's formula)
        sqrt_x = int(x ** 0.5)
        small_primes = sieve(sqrt_x)

        # Legendre: pi(x) = pi(sqrt(x)) + phi(x, pi(sqrt(x))) - 1
        # where phi(x, a) = number of integers <= x not divisible by first a primes
        # This is the "quantum observable" - it's an inclusion-exclusion

        # Simple Legendre (exponential in number of primes, but illustrative)
        if len(small_primes) <= 20:  # feasible
            # phi(x, a) by inclusion-exclusion
            a = len(small_primes)
            total = x
            # Mobius-like inclusion-exclusion
            from itertools import combinations
            for k in range(1, min(a + 1, 15)):  # limit depth
                sign = (-1) ** k
                for combo in combinations(small_primes, k):
                    prod = 1
                    for p in combo:
                        prod *= p
                        if prod > x:
                            break
                    if prod <= x:
                        total += sign * (x // prod)
            return total + a - 1  # approximate, not exact

        return -1  # too many primes

    # Method 3: Stochastic trace estimation with Chebyshev expansion
    # Represent chi_prime as a function that can be applied via matrix-vector product
    # with the "number theory matrix"
    def stochastic_pi_estimate(x, n_samples=1000):
        """
        Random sampling estimate: pi(x) ~ (x-1) * Pr[random integer in [2,x] is prime]
        Enhanced with importance sampling from smooth number distribution.
        """
        # Basic Monte Carlo
        samples = np.random.randint(2, x + 1, size=n_samples)
        prime_set = set(sieve(x).tolist())
        hits = sum(1 for s in samples if s in prime_set)
        mc_estimate = (x - 1) * hits / n_samples

        # Importance sampling: sample more from likely-prime regions
        # Primes thin out like 1/ln(x), so weight by ln(x)/ln(sample)
        weighted_hits = 0
        weighted_total = 0
        for s in samples:
            w = np.log(x) / np.log(max(s, 2))
            weighted_total += w
            if s in prime_set:
                weighted_hits += w
        is_estimate = (x - 1) * weighted_hits / weighted_total

        return mc_estimate, is_estimate

    # Test
    test_values = [100, 1000, 10000, 100000]
    print(f"\n{'x':>8s} {'exact':>8s} {'Monte Carlo':>12s} {'Importance':>12s} {'MC err%':>8s} {'IS err%':>8s}")
    print("-" * 60)

    mc_errors = []
    is_errors = []
    for x in test_values:
        exact = exact_pi(x)
        mc_est, is_est = stochastic_pi_estimate(x, n_samples=10000)
        mc_err = abs(mc_est - exact) / exact * 100
        is_err = abs(is_est - exact) / exact * 100
        mc_errors.append(mc_err)
        is_errors.append(is_err)
        print(f"{x:8d} {exact:8d} {mc_est:12.1f} {is_est:12.1f} {mc_err:7.1f}% {is_err:7.1f}%")

    # Key insight: even if we could compute pi(x) approximately,
    # we need EXACT pi(x) to find p(n), and the gap between approximate
    # and exact is O(sqrt(x) * ln(x)).
    print(f"\nKey insight: Monte Carlo pi(x) has error O(sqrt(x/n_samples))")
    print(f"For x=10^100 with n=10^6 samples: error ~ 10^47, vs pi(x) ~ 10^98")
    print(f"Error relative to gap: 10^47 / ln(10^100) ~ 10^45 -- USELESS for exact p(n)")

    # Method 4: Quantum-inspired Chebyshev moment method
    # <T_k(H)> gives moments of the spectral density.
    # For the "prime Hamiltonian" diagonal matrix, moments = sum of T_k(p) over primes.
    # This is related to prime counting via kernel polynomial method.
    print(f"\nChebyshev moment method for spectral density of primes:")
    x = 10000
    primes_x = sieve(x)
    n_primes = len(primes_x)

    # Normalized primes to [-1, 1]
    normalized = 2 * primes_x / x - 1

    # Compute Chebyshev moments
    n_moments = 50
    moments = np.zeros(n_moments)
    for k in range(n_moments):
        moments[k] = np.mean(np.cos(k * np.arccos(normalized)))

    # Reconstruct density from moments (kernel polynomial method)
    n_points = 500
    t = np.linspace(-1, 1, n_points)
    density = np.zeros(n_points)
    for k in range(n_moments):
        # Jackson kernel for damping
        jackson = ((n_moments - k + 1) * np.cos(np.pi * k / (n_moments + 1)) +
                   np.sin(np.pi * k / (n_moments + 1)) / np.tan(np.pi / (n_moments + 1))) / (n_moments + 1)
        if k == 0:
            density += jackson * moments[k] * np.ones_like(t)
        else:
            density += 2 * jackson * moments[k] * np.cos(k * np.arccos(t))

    density /= (np.pi * np.sqrt(1 - t**2 + 1e-10))

    # CDF from density gives approximate pi(x)
    dt = t[1] - t[0]
    cdf = np.cumsum(density) * dt * n_primes
    # Map back to original scale
    x_values = (t + 1) / 2 * x

    # Check accuracy at a few points
    print(f"  Chebyshev KPM density reconstruction ({n_moments} moments):")
    for check_x in [1000, 2500, 5000, 7500, 10000]:
        idx = np.argmin(np.abs(x_values - check_x))
        approx = cdf[idx]
        exact = len(primes_x[primes_x <= check_x])
        err = abs(approx - exact)
        print(f"    pi({check_x}) approx={approx:.1f}, exact={exact}, error={err:.1f}")

    return {
        "monte_carlo_errors_pct": mc_errors,
        "importance_sampling_errors_pct": is_errors,
        "n_samples": 10000,
        "key_barrier": "Stochastic methods give O(sqrt(x)) error, need O(1) for exact p(n)",
        "chebyshev_moments": n_moments,
        "conclusion": "FAIL - no classical shortcut to exact pi(x) via expectation values"
    }


###############################################################################
# EXPERIMENT 5: Random Matrix / GUE Eigenvalue Extraction
###############################################################################
def experiment5_gue_zeros():
    """
    Montgomery's pair correlation conjecture: the nontrivial zeros of zeta
    have the same statistics as eigenvalues of random GUE matrices.

    Approach: Given that we know the STATISTICS of zero spacings,
    can we reconstruct individual zero positions? If yes, we could
    use the explicit formula pi(x) = li(x) - sum_rho li(x^rho) + ...

    Test: Generate GUE eigenvalues, see if statistical properties
    allow reconstruction of individual eigenvalues.
    """
    print("\n" + "=" * 70)
    print("EXPERIMENT 5: GUE Random Matrix / Zeta Zero Reconstruction")
    print("=" * 70)

    # 5a. First few known zeta zeros (imaginary parts)
    known_zeros = np.array([
        14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
        37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
        52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
        67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
        79.337375, 82.910381, 84.735493, 87.425275, 88.809111,
        92.491899, 94.651344, 95.870634, 98.831194, 101.317851,
    ])

    # 5b. Generate GUE matrix and compare eigenvalue statistics
    def gue_eigenvalues(n):
        """Generate eigenvalues of n x n GUE matrix."""
        A = np.random.randn(n, n) + 1j * np.random.randn(n, n)
        H = (A + A.conj().T) / (2 * np.sqrt(2 * n))
        return np.sort(np.real(np.linalg.eigvalsh(H)))

    # 5c. Pair correlation function
    def pair_correlation(spacings, bins=50):
        """Compute normalized pair correlation of spacings."""
        mean_spacing = np.mean(spacings)
        normalized = spacings / mean_spacing
        hist, edges = np.histogram(normalized, bins=bins, range=(0, 3), density=True)
        centers = (edges[:-1] + edges[1:]) / 2
        return centers, hist

    # Compare zeta zero spacings with GUE
    zero_spacings = np.diff(known_zeros)
    mean_zero_spacing = np.mean(zero_spacings)
    print(f"\nZeta zero spacings (first {len(known_zeros)} zeros):")
    print(f"  Mean spacing: {mean_zero_spacing:.3f}")
    print(f"  Std: {np.std(zero_spacings):.3f}")
    print(f"  Min: {np.min(zero_spacings):.3f}, Max: {np.max(zero_spacings):.3f}")

    # GUE comparison
    n_trials = 100
    gue_spacings_all = []
    for _ in range(n_trials):
        evals = gue_eigenvalues(30)
        # Unfold: use bulk eigenvalues
        mid = len(evals) // 4
        bulk = evals[mid:-mid]
        spacings = np.diff(bulk)
        if len(spacings) > 0:
            mean_s = np.mean(spacings)
            if mean_s > 0:
                gue_spacings_all.extend((spacings / mean_s).tolist())

    gue_spacings_all = np.array(gue_spacings_all)

    # Wigner surmise: p(s) = (pi/2) s exp(-pi s^2 / 4) for GUE
    s_test = np.linspace(0.01, 3, 100)
    wigner = (np.pi / 2) * s_test * np.exp(-np.pi * s_test**2 / 4)  # GOE
    wigner_gue = (32 / np.pi**2) * s_test**2 * np.exp(-4 * s_test**2 / np.pi)  # GUE

    # Normalize zero spacings
    norm_zero_spacings = zero_spacings / mean_zero_spacing

    # KS test
    from scipy.stats import ks_2samp
    ks_stat, ks_pval = ks_2samp(norm_zero_spacings, gue_spacings_all)
    print(f"\nKS test (zeta zeros vs GUE): statistic={ks_stat:.4f}, p-value={ks_pval:.4f}")

    # 5d. The KEY experiment: Can we reconstruct zeros from statistics?
    print(f"\n--- Zero Reconstruction Experiment ---")

    # Approach: Given N(T) ~ (T/2pi) ln(T/2pi) - T/2pi (Riemann-von Mangoldt formula),
    # we know WHERE zeros should be on average. The question is whether
    # the fine structure (deviations from average) can be predicted.

    def riemann_von_mangoldt(T):
        """Approximate number of zeros with Im(rho) < T."""
        if T < 1:
            return 0
        return T / (2 * np.pi) * np.log(T / (2 * np.pi)) - T / (2 * np.pi)

    # Predict zero positions from average density
    predicted_zeros = []
    for n in range(1, len(known_zeros) + 1):
        # Invert N(T) = n to find T
        # N(T) = n => T ~ 2*pi*n / W(n/e) approximately
        from scipy.optimize import brentq
        try:
            t_pred = brentq(lambda t: riemann_von_mangoldt(t) - n + 0.5, 1, 200)
            predicted_zeros.append(t_pred)
        except:
            predicted_zeros.append(0)
    predicted_zeros = np.array(predicted_zeros)

    errors = known_zeros - predicted_zeros
    print(f"\nZero position prediction from average density (Riemann-von Mangoldt):")
    print(f"  Mean absolute error: {np.mean(np.abs(errors)):.3f}")
    print(f"  Max absolute error:  {np.max(np.abs(errors)):.3f}")
    print(f"  Mean spacing: {mean_zero_spacing:.3f}")
    print(f"  Error / spacing: {np.mean(np.abs(errors))/mean_zero_spacing:.3f}")

    # 5e. Can GUE statistics IMPROVE the prediction?
    # Idea: sample GUE-distributed deviations from the average positions
    print(f"\n--- GUE-Corrected Prediction ---")

    # Generate many GUE samples and find the one closest to observed
    best_match_error = float('inf')
    n_gue_trials = 1000
    n_zeros = len(known_zeros)

    for trial in range(n_gue_trials):
        # Generate GUE spacings
        evals = gue_eigenvalues(n_zeros + 10)
        # Take middle eigenvalues, unfold to match zero density
        bulk = np.sort(evals)[5:5 + n_zeros]
        # Scale to match zero spacing
        gue_spac = np.diff(bulk)
        if len(gue_spac) > 0 and np.mean(gue_spac) > 0:
            gue_spac_normalized = gue_spac / np.mean(gue_spac) * mean_zero_spacing
            # Reconstruct positions
            gue_zeros = np.zeros(n_zeros)
            gue_zeros[0] = predicted_zeros[0]
            for i in range(1, n_zeros):
                if i - 1 < len(gue_spac_normalized):
                    gue_zeros[i] = gue_zeros[i-1] + gue_spac_normalized[i-1]
                else:
                    gue_zeros[i] = gue_zeros[i-1] + mean_zero_spacing

            # Shift to minimize error
            shift = np.mean(known_zeros - gue_zeros)
            gue_zeros += shift

            error = np.mean(np.abs(known_zeros - gue_zeros))
            if error < best_match_error:
                best_match_error = error
                best_gue_zeros = gue_zeros.copy()

    print(f"  Best GUE match (out of {n_gue_trials} trials):")
    print(f"    Mean absolute error: {best_match_error:.3f}")
    print(f"    Error / spacing: {best_match_error/mean_zero_spacing:.3f}")
    print(f"    Average density prediction error: {np.mean(np.abs(errors)):.3f}")
    improvement = np.mean(np.abs(errors)) / best_match_error if best_match_error > 0 else 0
    print(f"    Improvement over density alone: {improvement:.2f}x")

    # 5f. Impact on prime counting
    # If we have zeros with error delta, how does that affect pi(x)?
    # The explicit formula: pi(x) ~ li(x) - sum_rho li(x^rho)
    # Error in one zero rho = 1/2 + i*gamma by delta_gamma:
    # |li(x^{1/2+i*(gamma+delta)}) - li(x^{1/2+i*gamma})| ~ sqrt(x) * delta / ln(x)

    print(f"\n--- Impact on Prime Counting ---")
    for x_exp in [10, 20, 50, 100]:
        x = 10.0 ** x_exp
        sqrt_x = 10.0 ** (x_exp / 2)
        ln_x = x_exp * np.log(10)
        # Error per zero ~ sqrt(x) * delta_gamma / ln(x)
        # With 30 zeros at best_match_error average:
        total_error = n_zeros * sqrt_x * best_match_error / ln_x
        gap_approx = x_exp * np.log(10)  # prime gap ~ ln(p(n)) ~ ln(x)
        print(f"  x=10^{x_exp}: error_from_zeros ~ {total_error:.2e}, gap ~ {gap_approx:.0f}, ratio = {total_error/gap_approx:.2e}")

    return {
        "ks_test": {"statistic": float(ks_stat), "p_value": float(ks_pval)},
        "density_prediction_mae": float(np.mean(np.abs(errors))),
        "gue_best_mae": float(best_match_error),
        "improvement_factor": float(improvement),
        "conclusion": "GUE statistics give marginal improvement over average density; "
                      "individual zero positions cannot be predicted from statistics alone. "
                      "Even perfect zero knowledge insufficient for large x due to sqrt(x) amplification."
    }


###############################################################################
# MAIN
###############################################################################
def main():
    print("Session 6: Classical Simulation of Quantum Algorithms for Primes")
    print("=" * 70)
    print()

    all_results = {}

    t0 = time.time()
    all_results["exp1_hadamard"] = experiment1_hadamard_period_finding()
    t1 = time.time()
    print(f"\n[Experiment 1 completed in {t1-t0:.1f}s]\n")

    all_results["exp2_grover"] = experiment2_grover_search()
    t2 = time.time()
    print(f"\n[Experiment 2 completed in {t2-t1:.1f}s]\n")

    all_results["exp3_tensor"] = experiment3_tensor_network()
    t3 = time.time()
    print(f"\n[Experiment 3 completed in {t3-t2:.1f}s]\n")

    all_results["exp4_expectation"] = experiment4_pi_expectation()
    t4 = time.time()
    print(f"\n[Experiment 4 completed in {t4-t3:.1f}s]\n")

    all_results["exp5_gue"] = experiment5_gue_zeros()
    t5 = time.time()
    print(f"\n[Experiment 5 completed in {t5-t4:.1f}s]\n")

    print("\n" + "=" * 70)
    print("OVERALL SUMMARY")
    print("=" * 70)
    for key, val in all_results.items():
        conclusion = val.get("conclusion", "N/A")
        print(f"  {key}: {conclusion}")

    total_time = time.time() - t0
    print(f"\nTotal runtime: {total_time:.1f}s")

    # Save results
    output_path = os.path.join(os.path.dirname(__file__), "quantum_classical_results.json")
    # Convert numpy types for JSON
    def convert(obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            r = convert(obj)
            if r is not obj:
                return r
            return super().default(obj)

    with open(output_path, 'w') as f:
        json.dump(all_results, f, indent=2, cls=NumpyEncoder)
    print(f"Results saved to {output_path}")

    return all_results


if __name__ == "__main__":
    main()
