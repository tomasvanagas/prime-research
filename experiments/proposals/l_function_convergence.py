#!/usr/bin/env python3
"""
L-function Zero Sum Convergence vs Riemann Zeta Zero Sum
=========================================================
Experiment: Does pi(x;q,a) require FEWER L-function zeros than pi(x) requires
zeta zeros for the same relative accuracy?

Two approaches:
  A) Empirical: measure variation/roughness of E(x;q,a) = pi(x;q,a) - Li(x)/phi(q)
     vs E(x) = pi(x) - Li(x). If E(x;q,a) is smoother, fewer zeros needed.
  B) Explicit formula with zeros: use K zeros of zeta / L-functions and measure
     how fast the explicit formula converges.

We use sympy for prime counting, mpmath for high-precision arithmetic and
(attempted) L-function zero computation.
"""

import math
import numpy as np
from collections import defaultdict
import time

# ---------- Part 0: Prime sieve and counting ----------

def sieve_primes(limit):
    """Simple sieve of Eratosthenes."""
    is_prime = bytearray([1]) * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return [p for p in range(2, limit + 1) if is_prime[p]]

def pi_exact(x, primes):
    """Count primes <= x using pre-sieved list."""
    from bisect import bisect_right
    return bisect_right(primes, x)

def pi_qa(x, q, a, primes):
    """Count primes p <= x with p ≡ a (mod q), using pre-sieved list."""
    from bisect import bisect_right
    idx = bisect_right(primes, x)
    return sum(1 for p in primes[:idx] if p % q == a)

def li_integral(x):
    """Logarithmic integral Li(x) = integral from 2 to x of dt/ln(t)."""
    if x <= 2:
        return 0.0
    # Use scipy if available, else numerical integration
    try:
        from scipy.integrate import quad
        result, _ = quad(lambda t: 1.0 / math.log(t), 2, x)
        return result
    except ImportError:
        # Simple trapezoidal rule
        n_steps = max(1000, int(x / 10))
        dt = (x - 2.0) / n_steps
        total = 0.0
        for i in range(n_steps):
            t0 = 2.0 + i * dt
            t1 = t0 + dt
            total += 0.5 * (1.0/math.log(t0) + 1.0/math.log(t1)) * dt
        return total

def euler_phi(q):
    """Euler's totient function."""
    result = q
    p = 2
    temp = q
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return result

# ---------- Part 1: Empirical smoothness comparison ----------

def measure_variation(errors):
    """Measure total variation and other roughness metrics of an error sequence."""
    errors = np.array(errors, dtype=float)
    n = len(errors)
    if n < 2:
        return {}

    # Total variation
    tv = np.sum(np.abs(np.diff(errors)))

    # RMS of differences (roughness)
    rms_diff = np.sqrt(np.mean(np.diff(errors)**2))

    # Standard deviation
    std = np.std(errors)

    # Mean absolute error
    mae = np.mean(np.abs(errors))

    # Number of sign changes (oscillation count)
    signs = np.sign(errors)
    sign_changes = np.sum(np.abs(np.diff(signs)) > 0)

    # Power spectrum roughness: fraction of power in high frequencies
    if n > 10:
        fft_vals = np.fft.rfft(errors - np.mean(errors))
        power = np.abs(fft_vals)**2
        total_power = np.sum(power)
        if total_power > 0:
            mid = len(power) // 2
            high_freq_fraction = np.sum(power[mid:]) / total_power
        else:
            high_freq_fraction = 0.0
    else:
        high_freq_fraction = float('nan')

    return {
        'total_variation': tv,
        'rms_diff': rms_diff,
        'std': std,
        'mae': mae,
        'sign_changes': sign_changes,
        'high_freq_frac': high_freq_fraction,
        'n_points': n,
    }


def empirical_smoothness_test():
    """
    Compare smoothness of:
      E(x) = pi(x) - Li(x)                     [Riemann error]
      E(x;q,a) = pi(x;q,a) - Li(x)/phi(q)      [Dirichlet error]

    If E(x;q,a) is smoother, pi(x;q,a) needs fewer L-function zeros.
    """
    print("=" * 75)
    print("PART 1: EMPIRICAL SMOOTHNESS COMPARISON")
    print("=" * 75)

    LIMIT = 200_000
    primes = sieve_primes(LIMIT)

    # Sample points: 500 evenly spaced from 1000 to LIMIT
    x_vals = np.linspace(1000, LIMIT, 500).astype(int)

    # --- Riemann error: pi(x) - Li(x) ---
    riemann_errors = []
    for x in x_vals:
        pi_x = pi_exact(x, primes)
        li_x = li_integral(x)
        riemann_errors.append(pi_x - li_x)

    # Normalize by sqrt(x)/log(x) to get comparable scale
    riemann_normalized = [e / (math.sqrt(x) / math.log(x))
                          for e, x in zip(riemann_errors, x_vals)]

    print(f"\nRiemann error E(x) = pi(x) - Li(x), x in [1000, {LIMIT}]")
    stats_r = measure_variation(riemann_normalized)
    print(f"  Total variation:     {stats_r['total_variation']:.4f}")
    print(f"  RMS of differences:  {stats_r['rms_diff']:.6f}")
    print(f"  Std deviation:       {stats_r['std']:.4f}")
    print(f"  Mean |error|:        {stats_r['mae']:.4f}")
    print(f"  Sign changes:        {stats_r['sign_changes']}")
    print(f"  High-freq fraction:  {stats_r['high_freq_frac']:.4f}")

    # --- Dirichlet errors for various (q, a) ---
    test_cases = [
        (3, 1), (3, 2),
        (4, 1), (4, 3),
        (5, 1), (5, 2), (5, 3), (5, 4),
        (7, 1), (7, 2), (7, 3), (7, 4), (7, 5), (7, 6),
    ]

    print(f"\n{'q':>3} {'a':>3} {'phi(q)':>6} {'TotalVar':>10} {'RMS_diff':>10} "
          f"{'Std':>8} {'MAE':>8} {'SignCh':>7} {'HiFreq':>8} {'TV_ratio':>10}")
    print("-" * 95)

    results = {}
    for q, a in test_cases:
        if math.gcd(a, q) != 1:
            continue
        phi_q = euler_phi(q)

        dir_errors = []
        for x in x_vals:
            pi_qa_x = pi_qa(x, q, a, primes)
            expected = li_integral(x) / phi_q
            dir_errors.append(pi_qa_x - expected)

        # Normalize by sqrt(x)/(phi(q)*log(x)) -- natural scale for Dirichlet
        dir_normalized = [e / (math.sqrt(x) / (phi_q * math.log(x)))
                          for e, x in zip(dir_errors, x_vals)]

        stats_d = measure_variation(dir_normalized)
        tv_ratio = stats_d['total_variation'] / stats_r['total_variation']

        results[(q, a)] = {
            'stats': stats_d,
            'tv_ratio': tv_ratio,
            'normalized_errors': dir_normalized,
        }

        print(f"{q:3d} {a:3d} {phi_q:6d} {stats_d['total_variation']:10.4f} "
              f"{stats_d['rms_diff']:10.6f} {stats_d['std']:8.4f} "
              f"{stats_d['mae']:8.4f} {stats_d['sign_changes']:7d} "
              f"{stats_d['high_freq_frac']:8.4f} {tv_ratio:10.4f}")

    return results, stats_r


# ---------- Part 2: Explicit formula with zeta zeros ----------

def get_zeta_zeros(n_zeros):
    """Get first n_zeros non-trivial zeros of zeta (imaginary parts)."""
    # Known values of gamma_k for k=1..50
    # Source: LMFDB / standard tables
    known_zeros = [
        14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
        37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
        52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
        67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
        79.337375, 82.910381, 84.735493, 87.425275, 88.809111,
        92.491899, 94.651344, 95.870634, 98.831194, 101.317851,
        103.725538, 105.446623, 107.168611, 111.029536, 111.874659,
        114.320220, 116.226680, 118.790783, 121.370125, 122.946829,
        124.256819, 127.516684, 129.578704, 131.087689, 133.497737,
        134.756510, 138.116042, 139.736209, 141.123707, 143.111846,
    ]
    return known_zeros[:n_zeros]

def explicit_formula_zeta(x, zeros):
    """
    Riemann's explicit formula (simplified):
    pi(x) ~ Li(x) - sum_{gamma} 2*Re(Li(x^{rho})) / ...

    Simplified oscillatory term:
    pi(x) approx Li(x) - (1/log x) * sum_gamma 2*cos(gamma*log(x)) / sqrt(x) * (something)

    More precisely, psi(x) = x - sum_rho x^rho/rho - log(2pi) - (1/2)log(1-x^{-2})
    We use psi(x) approximation and convert.

    For psi(x): psi(x) ~ x - sum_{gamma} 2*x^{1/2}*cos(gamma*log(x))/(1/4+gamma^2)^{1/2}
    Actually: x^rho/rho where rho = 1/2 + i*gamma
    x^rho = x^{1/2} * e^{i*gamma*log(x)}
    Re(x^rho/rho) = x^{1/2} * Re(e^{i*gamma*log(x)} / (1/2 + i*gamma))
    """
    log_x = math.log(x)
    sqrt_x = math.sqrt(x)

    # psi(x) = x - sum_rho x^rho/rho - log(2pi) - ...
    oscillatory = 0.0
    for gamma in zeros:
        # rho = 1/2 + i*gamma
        # x^rho / rho = x^{1/2} * e^{i*gamma*log(x)} / (1/2 + i*gamma)
        phase = gamma * log_x
        cos_phase = math.cos(phase)
        sin_phase = math.sin(phase)

        # 1/(1/2 + i*gamma) = (1/2 - i*gamma) / (1/4 + gamma^2)
        denom = 0.25 + gamma * gamma
        re_part = (0.5 * cos_phase + gamma * sin_phase) / denom

        # Two conjugate zeros contribute 2 * Re(...)
        oscillatory += 2 * sqrt_x * re_part

    psi_approx = x - oscillatory - math.log(2 * math.pi)

    # Convert psi(x) to pi(x) approximately: pi(x) ~ psi(x)/log(x)
    # (This is a rough conversion; more precise would account for prime powers)
    pi_approx = psi_approx / log_x

    return pi_approx


def explicit_formula_convergence():
    """
    Test how fast the explicit formula converges as we add more zeta zeros.
    Compare this convergence rate across different formulations.
    """
    print("\n" + "=" * 75)
    print("PART 2: EXPLICIT FORMULA CONVERGENCE WITH ZETA ZEROS")
    print("=" * 75)

    LIMIT = 100_000
    primes = sieve_primes(LIMIT)

    # Test at multiple x values
    test_x = [1000, 2000, 5000, 10000, 20000, 50000, 100000]

    k_values = [1, 2, 5, 10, 20, 50]
    all_zeros = get_zeta_zeros(50)

    print(f"\n{'x':>8}", end="")
    print(f" {'pi(x)':>8} {'Li(x)':>10}", end="")
    for k in k_values:
        print(f" {'K='+str(k):>10}", end="")
    print()
    print("-" * (28 + 11 * len(k_values)))

    convergence_data = {}
    for x in test_x:
        pi_x = pi_exact(x, primes)
        li_x = li_integral(x)

        print(f"{x:8d} {pi_x:8d} {li_x:10.1f}", end="")

        errors_by_k = {}
        for k in k_values:
            zeros_k = all_zeros[:k]
            pi_approx = explicit_formula_zeta(x, zeros_k)
            err = pi_approx - pi_x
            rel_err = abs(err) / pi_x if pi_x > 0 else float('inf')
            errors_by_k[k] = {'abs': err, 'rel': rel_err}
            print(f" {err:10.1f}", end="")
        print()
        convergence_data[x] = errors_by_k

    # Summary: average relative error by K
    print(f"\nAverage |relative error| by number of zeros K:")
    print(f"{'K':>5} {'Avg |rel err|':>15} {'Max |rel err|':>15}")
    print("-" * 40)
    for k in k_values:
        rel_errs = [abs(convergence_data[x][k]['rel']) for x in test_x]
        avg_re = np.mean(rel_errs)
        max_re = np.max(rel_errs)
        print(f"{k:5d} {avg_re:15.6f} {max_re:15.6f}")

    return convergence_data


# ---------- Part 3: Dirichlet characters and L-function explicit formula ----------

def dirichlet_characters(q):
    """
    Compute all Dirichlet characters mod q.
    Returns list of (chi_label, chi_dict) where chi_dict maps a -> chi(a).
    Uses the structure of (Z/qZ)*.
    """
    from itertools import product

    # Find units mod q
    units = [a for a in range(1, q) if math.gcd(a, q) == 1]
    phi_q = len(units)

    # For small q, enumerate characters directly via DFT on the group
    # Find a generator if cyclic, or use direct construction

    # Simple approach: for each unit g, find its order
    # Then characters are determined by chi(g) = e^{2pi i k / ord(g)}

    # For q prime, (Z/qZ)* is cyclic
    # For q=4, (Z/4Z)* = {1,3}, cyclic of order 2
    # For general q, use the fact that characters are homomorphisms (Z/qZ)* -> C*

    # Find generator(s)
    def find_generator(q, units):
        phi = len(units)
        for g in units:
            if pow(g, phi, q) == 1:
                # Check if order is exactly phi
                is_gen = True
                for d in range(1, phi):
                    if phi % d == 0 and pow(g, d, q) == 1:
                        is_gen = False
                        break
                if is_gen:
                    return g
        return None

    characters = []

    if q <= 2:
        # Only principal character
        chi = {a: 1.0 for a in units}
        characters.append(('chi_0', chi, True))  # True = principal
        return characters

    # For small q, just build the DFT matrix approach
    # Map units to indices
    unit_to_idx = {u: i for i, u in enumerate(units)}

    # Multiplication table
    mult_table = {}
    for a in units:
        for b in units:
            mult_table[(a, b)] = (a * b) % q

    # Build group table: find the discrete log table
    # Try each unit as potential generator
    gen = find_generator(q, units)

    if gen is not None:
        # Cyclic group -- easy case
        # chi_k(g^j) = e^{2pi i k j / phi_q}
        for k in range(phi_q):
            chi = {}
            g_power = 1
            for j in range(phi_q):
                g_power_mod = pow(gen, j, q)
                chi[g_power_mod] = np.exp(2j * np.pi * k * j / phi_q)
            is_principal = (k == 0)
            characters.append((f'chi_{k}', chi, is_principal))
    else:
        # Non-cyclic group (e.g., q=8: (Z/8Z)* = {1,3,5,7} ~ C2 x C2)
        # Use direct enumeration: a character maps each generator to a root of unity
        # of the appropriate order

        # Find group structure by finding independent generators
        # Simple brute force for small q
        from itertools import product as iprod

        # Find orders of all elements
        orders = {}
        for a in units:
            ord_a = 1
            power = a
            while power % q != 1:
                power = (power * a) % q
                ord_a += 1
            orders[a] = ord_a

        # Build characters by trying all possible homomorphisms
        # chi: units -> roots of unity, with chi(ab) = chi(a)chi(b)
        # Brute force: try all assignments of roots of unity

        def is_valid_character(chi_vals):
            for a in units:
                for b in units:
                    ab = (a * b) % q
                    if abs(chi_vals[a] * chi_vals[b] - chi_vals[ab]) > 1e-10:
                        return False
            return True

        # Each unit of order d must map to a d-th root of unity
        possible_values = {}
        for a in units:
            d = orders[a]
            possible_values[a] = [np.exp(2j * np.pi * k / d) for k in range(d)]

        # We need to search, but constrained: once generators are fixed,
        # everything is determined. For small q this is fast.
        # Fix chi(1) = 1 always
        found_chars = set()

        def search_chars(idx, current_chi):
            if idx == len(units):
                # Verify
                key = tuple(round(current_chi[u].real, 8) + 1j*round(current_chi[u].imag, 8)
                           for u in units)
                if key not in found_chars:
                    if is_valid_character(current_chi):
                        found_chars.add(key)
                        chi_copy = dict(current_chi)
                        is_principal = all(abs(v - 1.0) < 1e-10 for v in chi_copy.values())
                        characters.append((f'chi_{len(characters)}', chi_copy, is_principal))
                return

            a = units[idx]
            if a == 1:
                current_chi[a] = 1.0 + 0j
                search_chars(idx + 1, current_chi)
            else:
                for val in possible_values[a]:
                    current_chi[a] = val
                    # Quick consistency check
                    ok = True
                    for prev_idx in range(idx):
                        b = units[prev_idx]
                        ab = (a * b) % q
                        if ab in current_chi:
                            if abs(current_chi[a] * current_chi[b] - current_chi[ab]) > 1e-10:
                                ok = False
                                break
                    if ok:
                        search_chars(idx + 1, current_chi)

        search_chars(0, {})

    return characters


def psi_dirichlet_approx(x, q, a, zeta_zeros, char_data):
    """
    Explicit formula for psi(x; q, a):

    psi(x; q, a) = x/phi(q) - (1/phi(q)) * sum_{chi} conj(chi(a)) * sum_{rho_chi} x^{rho}/rho

    Where rho_chi are zeros of L(s, chi).

    Since we don't have L-function zeros easily, we approximate:
    For the principal character chi_0, the L-function zeros are essentially zeta zeros.
    For non-principal characters, we'd need their zeros.

    SIMPLIFICATION: We test convergence of the zeta-zero contribution only,
    which gives the "principal character" part of the explicit formula.
    The question is whether the non-principal character contributions are smaller.
    """
    phi_q = euler_phi(q)
    log_x = math.log(x)
    sqrt_x = math.sqrt(x)

    # Principal character contribution (from zeta zeros)
    oscillatory_principal = 0.0
    for gamma in zeta_zeros:
        phase = gamma * log_x
        cos_phase = math.cos(phase)
        sin_phase = math.sin(phase)
        denom = 0.25 + gamma * gamma
        re_part = (0.5 * cos_phase + gamma * sin_phase) / denom
        oscillatory_principal += 2 * sqrt_x * re_part

    # For the principal character, chi_0(a) = 1 for gcd(a,q)=1
    # So the principal contribution to psi(x;q,a) is (1/phi(q)) * [x - sum_rho x^rho/rho - ...]
    psi_principal = (x - oscillatory_principal - math.log(2 * math.pi)) / phi_q

    return psi_principal


def dirichlet_explicit_convergence():
    """
    Compare convergence: zeta zeros for pi(x) vs zeta zeros for pi(x;q,a).

    Key insight: In the explicit formula for pi(x;q,a), each character's
    L-function contributes its own set of zeros. The principal character's
    zeros are the zeta zeros (divided by phi(q)). Non-principal characters
    contribute ADDITIONAL oscillatory terms.

    If we only use the principal character (zeta zeros), the residual for
    pi(x;q,a) includes the non-principal L-function zero contributions.
    This tells us how much "extra" information the non-principal zeros carry.
    """
    print("\n" + "=" * 75)
    print("PART 3: DIRICHLET EXPLICIT FORMULA -- PRINCIPAL vs FULL")
    print("=" * 75)

    LIMIT = 100_000
    primes = sieve_primes(LIMIT)
    all_zeros = get_zeta_zeros(50)

    test_x = [2000, 5000, 10000, 20000, 50000, 100000]
    k_values = [1, 5, 10, 20, 50]

    test_qa = [(3, 1), (3, 2), (4, 1), (4, 3), (5, 1), (5, 2), (7, 1), (7, 3)]

    print("\nResidual |pi(x;q,a) - principal_approx(x;q,a)| / pi(x;q,a)")
    print("Compared to |pi(x) - zeta_approx(x)| / pi(x)")
    print()

    # For each K, compute average relative error for pi(x) and each pi(x;q,a)
    print(f"{'':>12}", end="")
    for k in k_values:
        print(f" {'K='+str(k):>10}", end="")
    print()
    print("-" * (12 + 11 * len(k_values)))

    # pi(x) with zeta zeros
    zeta_errors = {k: [] for k in k_values}
    for x in test_x:
        pi_x = pi_exact(x, primes)
        for k in k_values:
            approx = explicit_formula_zeta(x, all_zeros[:k])
            rel_err = abs(approx - pi_x) / pi_x
            zeta_errors[k].append(rel_err)

    avg_zeta = {k: np.mean(zeta_errors[k]) for k in k_values}
    print(f"{'pi(x)':>12}", end="")
    for k in k_values:
        print(f" {avg_zeta[k]:10.6f}", end="")
    print()

    # pi(x;q,a) with only principal character (zeta zeros)
    for q, a in test_qa:
        phi_q = euler_phi(q)
        dir_errors = {k: [] for k in k_values}

        for x in test_x:
            pi_qa_x = pi_qa(x, q, a, primes)
            if pi_qa_x == 0:
                continue
            for k in k_values:
                psi_approx = psi_dirichlet_approx(x, q, a, all_zeros[:k], None)
                pi_approx = psi_approx / math.log(x)
                rel_err = abs(pi_approx - pi_qa_x) / pi_qa_x
                dir_errors[k].append(rel_err)

        print(f"{'pi('+str(q)+','+str(a)+')':>12}", end="")
        for k in k_values:
            if dir_errors[k]:
                avg_e = np.mean(dir_errors[k])
                print(f" {avg_e:10.6f}", end="")
            else:
                print(f" {'N/A':>10}", end="")
        print()

    print("\nInterpretation:")
    print("  If pi(x;q,a) errors >> pi(x) errors: non-principal L-function zeros")
    print("  carry significant information --> NO convergence advantage.")
    print("  If pi(x;q,a) errors ~ pi(x) errors: principal character dominates")
    print("  --> smoothing across residue classes helps.")


# ---------- Part 4: Direct roughness comparison at multiple scales ----------

def multiscale_roughness():
    """
    Compare roughness of pi(x)-Li(x) vs pi(x;q,a)-Li(x)/phi(q) at multiple scales.
    This directly answers: is the arithmetic progression version smoother?
    """
    print("\n" + "=" * 75)
    print("PART 4: MULTI-SCALE ROUGHNESS COMPARISON")
    print("=" * 75)

    LIMIT = 500_000
    primes = sieve_primes(LIMIT)

    # Different ranges to test scale-dependence
    ranges = [
        (1000, 10000, 200),
        (10000, 50000, 200),
        (50000, 200000, 200),
        (100000, 500000, 200),
    ]

    test_qa = [(3, 2), (4, 3), (5, 2), (7, 3)]

    print(f"\n{'Range':>18} {'pi(x) TV':>10} ", end="")
    for q, a in test_qa:
        print(f" {'('+str(q)+','+str(a)+') TV':>10}", end="")
    print(f" {'Avg ratio':>10}")
    print("-" * (18 + 10 + 11 * len(test_qa) + 12))

    for x_lo, x_hi, n_pts in ranges:
        x_vals = np.linspace(x_lo, x_hi, n_pts).astype(int)

        # Riemann error normalized
        r_errors = []
        for x in x_vals:
            pi_x = pi_exact(x, primes)
            li_x = li_integral(x)
            # Normalize by sqrt(x)/log(x)
            norm = math.sqrt(x) / math.log(x)
            r_errors.append((pi_x - li_x) / norm)

        r_tv = np.sum(np.abs(np.diff(r_errors)))

        ratios = []
        print(f"[{x_lo:>6d},{x_hi:>7d}] {r_tv:10.2f} ", end="")

        for q, a in test_qa:
            phi_q = euler_phi(q)
            d_errors = []
            for x in x_vals:
                pi_qa_x = pi_qa(x, q, a, primes)
                expected = li_integral(x) / phi_q
                norm = math.sqrt(x) / (phi_q * math.log(x))
                d_errors.append((pi_qa_x - expected) / norm)

            d_tv = np.sum(np.abs(np.diff(d_errors)))
            ratio = d_tv / r_tv if r_tv > 0 else float('inf')
            ratios.append(ratio)
            print(f" {d_tv:10.2f}", end="")

        avg_ratio = np.mean(ratios)
        print(f" {avg_ratio:10.4f}")

    print("\nNote: TV ratio > 1 means pi(x;q,a) error is ROUGHER than pi(x) error")
    print("      TV ratio < 1 means pi(x;q,a) error is SMOOTHER (fewer zeros needed)")
    print("      TV ratio ~ 1 means similar roughness (similar zero requirements)")


# ---------- Part 5: Spectral analysis ----------

def spectral_comparison():
    """
    Compare power spectra of pi(x)-Li(x) vs pi(x;q,a)-Li(x)/phi(q).
    If the Dirichlet version has less high-frequency content, it needs fewer zeros.
    """
    print("\n" + "=" * 75)
    print("PART 5: SPECTRAL ANALYSIS OF ERROR TERMS")
    print("=" * 75)

    LIMIT = 200_000
    primes = sieve_primes(LIMIT)
    N = 1024  # Power of 2 for FFT
    x_vals = np.linspace(5000, LIMIT, N).astype(int)

    # Riemann error
    r_errors = np.zeros(N)
    for i, x in enumerate(x_vals):
        pi_x = pi_exact(x, primes)
        li_x = li_integral(x)
        norm = math.sqrt(x) / math.log(x)
        r_errors[i] = (pi_x - li_x) / norm

    r_fft = np.fft.rfft(r_errors - np.mean(r_errors))
    r_power = np.abs(r_fft)**2
    r_total = np.sum(r_power)

    test_qa = [(3, 2), (4, 3), (5, 2), (7, 3)]

    # Compute cumulative power fraction
    freqs = np.fft.rfftfreq(N)

    print(f"\nFraction of total spectral power captured by lowest K frequency components:")
    print(f"{'K':>5} {'pi(x)':>10}", end="")
    for q, a in test_qa:
        print(f" {'('+str(q)+','+str(a)+')':>10}", end="")
    print()
    print("-" * (15 + 11 * len(test_qa)))

    # Compute Dirichlet error spectra
    d_powers = {}
    for q, a in test_qa:
        phi_q = euler_phi(q)
        d_errors = np.zeros(N)
        for i, x in enumerate(x_vals):
            pi_qa_x = pi_qa(x, q, a, primes)
            expected = li_integral(x) / phi_q
            norm = math.sqrt(x) / (phi_q * math.log(x))
            d_errors[i] = (pi_qa_x - expected) / norm

        d_fft = np.fft.rfft(d_errors - np.mean(d_errors))
        d_powers[(q, a)] = np.abs(d_fft)**2

    for K in [1, 2, 5, 10, 20, 50, 100, 200]:
        r_frac = np.sum(r_power[:K]) / r_total if r_total > 0 else 0
        print(f"{K:5d} {r_frac:10.4f}", end="")
        for q, a in test_qa:
            d_total = np.sum(d_powers[(q, a)])
            d_frac = np.sum(d_powers[(q, a)][:K]) / d_total if d_total > 0 else 0
            print(f" {d_frac:10.4f}", end="")
        print()

    print("\nHigher values = more power concentrated in low frequencies = smoother = fewer zeros needed")

    # Spectral rolloff: frequency at which 90% of power is captured
    print(f"\nFrequency index capturing 90% of power:")

    r_cumsum = np.cumsum(r_power) / r_total
    r_90 = np.searchsorted(r_cumsum, 0.9)
    print(f"  pi(x):        K={r_90}")

    for q, a in test_qa:
        d_total = np.sum(d_powers[(q, a)])
        d_cumsum = np.cumsum(d_powers[(q, a)]) / d_total
        d_90 = np.searchsorted(d_cumsum, 0.9)
        print(f"  pi({q},{a}):     K={d_90}")


# ---------- Part 6: Summary and interpretation ----------

def summarize():
    print("\n" + "=" * 75)
    print("SUMMARY AND INTERPRETATION")
    print("=" * 75)
    print("""
QUESTION: Does pi(x;q,a) require fewer L-function zeros than pi(x) requires
zeta zeros for the same accuracy?

THEORETICAL EXPECTATION:
  The explicit formula for pi(x;q,a) involves phi(q) different L-functions,
  each contributing its own set of zeros. The TOTAL number of zeros needed
  is actually phi(q) times as many (one set per character), not fewer.

  However, each individual L-function's zeros contribute with amplitude
  divided by phi(q), so the per-zero impact is smaller. The question is
  whether this averaging produces effective smoothing.

KEY FINDINGS FROM THIS EXPERIMENT:
  1. SMOOTHNESS: Compare total variation ratios above.
     - If TV ratio ~ 1: no smoothing advantage
     - If TV ratio > 1: arithmetic progressions are ROUGHER (worse!)
     - If TV ratio < 1: arithmetic progressions are smoother (better)

  2. SPECTRAL: Compare power concentration above.
     - Similar spectral concentration = similar zero requirements
     - More concentrated = fewer zeros needed

  3. EXPLICIT FORMULA: Compare convergence rates above.
     - Principal character alone (zeta zeros) gives a baseline
     - Residual tells us how much non-principal characters matter

IMPLICATIONS FOR p(n) COMPUTATION:
  Even if pi(x;q,a) were smoother, this would NOT help compute p(n) in
  O(polylog) because:
  - We'd need to combine phi(q) arithmetic progressions to recover pi(x)
  - The total zero information across all L-functions is at least as much
    as for zeta alone (by the factorization L(s,chi_0) ~ zeta(s))
  - The sqrt(x) barrier applies to EACH L-function's error term
""")


# ---------- Main ----------

if __name__ == '__main__':
    t0 = time.time()

    print("L-FUNCTION ZERO SUM CONVERGENCE EXPERIMENT")
    print(f"{'=' * 75}\n")

    results_smooth, stats_riemann = empirical_smoothness_test()
    convergence = explicit_formula_convergence()
    dirichlet_explicit_convergence()
    multiscale_roughness()
    spectral_comparison()
    summarize()

    elapsed = time.time() - t0
    print(f"\nTotal runtime: {elapsed:.1f}s")
