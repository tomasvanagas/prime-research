#!/usr/bin/env python3
"""
Zeta Zero Truncation Analysis
==============================
Key question: How many zeta zeros T are needed so that
  |pi(x) - R(x) + sum_{k=1}^{T} 2*Re(R(x^{rho_k}))| < 0.5
allowing exact rounding to pi(x)?

The explicit formula (Riemann):
  pi(x) = R(x) - sum_{rho} R(x^rho) + small corrections

We use the Gram series for R with mpmath at sufficient precision.
For R(x^rho) where rho = 1/2 + i*gamma, we need enough precision to
handle the oscillating terms (gamma * log(x) can be large).

Also: For what fraction of x in [N,2N] is round(R(x)) = pi(x)?
"""

import os
import sys
import time
import numpy as np
from collections import defaultdict

try:
    import mpmath
except ImportError:
    print("ERROR: mpmath required. pip install mpmath")
    sys.exit(1)

DATA_DIR = "/apps/aplikacijos/prime-research/data"


def load_zeros(max_zeros=1000):
    """Load imaginary parts of zeta zeros from data files."""
    for n in [1000, 500, 300, 200]:
        path = os.path.join(DATA_DIR, f"zeta_zeros_{n}.txt")
        if os.path.exists(path):
            zeros = []
            with open(path) as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split()
                    if len(parts) >= 2:
                        zeros.append(parts[1])  # Keep as string for precision
                    else:
                        zeros.append(parts[0])
            zeros = zeros[:max_zeros]
            print(f"Loaded {len(zeros)} zero strings from {path}")
            return zeros
    raise FileNotFoundError("No zeta zero files found in data/")


def pi_exact(x):
    """Count primes <= x using sieve of Eratosthenes."""
    x = int(x)
    if x < 2:
        return 0
    sieve = bytearray(b'\x01') * (x + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(x**0.5) + 1):
        if sieve[i]:
            sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
    return sum(sieve)


def R_gram_complex(z, n_terms=300):
    """
    Compute R(e^z) via Gram series:
      R(e^z) = 1 + sum_{k=1}^{n_terms} z^k / (k * k! * zeta(k+1))

    z can be real or complex. All arithmetic done in current mpmath precision.
    """
    one = mpmath.mpf(1)
    result = mpmath.mpc(one, 0) if isinstance(z, mpmath.mpc) else one
    z_power = mpmath.mpc(one, 0) if isinstance(z, mpmath.mpc) else one
    factorial = one

    for k in range(1, n_terms + 1):
        z_power *= z
        factorial *= k
        zeta_val = mpmath.zeta(k + 1)
        term = z_power / (k * factorial * zeta_val)
        result += term
        # Convergence check: relative to overall size
        if abs(term) < mpmath.power(10, -mpmath.mp.dps + 10):
            break

    return result


def compute_zero_contrib(log_x, gamma_str):
    """
    Compute 2*Re(R(x^rho)) where rho = 1/2 + i*gamma.
    log(x^rho) = rho * log(x) = (1/2)*log(x) + i*gamma*log(x)

    We set precision based on the magnitude of gamma*log(x) to ensure
    the oscillating Gram series converges accurately.
    """
    gamma = mpmath.mpf(gamma_str)
    # The argument to Gram series is z = rho * log_x
    z = mpmath.mpc(log_x / 2, gamma * log_x)

    # Need enough terms for convergence. |z|^k / (k * k!) must become small.
    # |z| ~ gamma * log_x for large gamma.
    # k! grows faster, but we need k ~ |z| terms.
    mag = float(abs(z))
    n_terms = max(200, int(mag * 2.5) + 50)

    R_val = R_gram_complex(z, n_terms=n_terms)
    return 2 * mpmath.re(R_val)


def set_precision_for(x, max_gamma):
    """
    Set mpmath precision high enough for the computation.
    The Gram series terms involve z^k / k! where |z| ~ max_gamma * log(x).
    Peak term is around k ~ |z|, with magnitude ~ |z|^|z| / |z|! ~ e^|z| / sqrt(2*pi*|z|).
    We need enough decimal digits to handle this intermediate growth.
    """
    log_x = float(np.log(x))
    z_mag = max_gamma * log_x
    # Stirling: peak ~ exp(z_mag) roughly, so we need z_mag / ln(10) extra digits
    extra_digits = int(z_mag / 2.302585) + 20
    dps = max(50, extra_digits + 30)
    mpmath.mp.dps = dps
    return dps


def analyze_truncation(x_values, T_values, gamma_strs):
    """For each x, compute error as function of number of zeros used."""

    print(f"\n{'='*80}")
    print(f"ZETA ZERO TRUNCATION ANALYSIS")
    print(f"{'='*80}")
    print(f"Available zeros: {len(gamma_strs)}")
    print()

    results = {}

    for x in x_values:
        print(f"\n--- x = {x:,} ---")
        t0 = time.time()

        pi_val = pi_exact(x)

        # For R(x) computation, moderate precision suffices
        mpmath.mp.dps = 50
        log_x_mp = mpmath.log(mpmath.mpf(x))
        R_val = R_gram_complex(log_x_mp, n_terms=300)
        R_float = float(R_val)
        err_R_only = R_float - pi_val

        print(f"  pi({x:,}) = {pi_val}")
        print(f"  R({x:,}) = {R_float:.6f}")
        print(f"  R(x) - pi(x) = {err_R_only:.6f}")
        print(f"  round(R(x)) = {round(R_float)}, correct = {round(R_float) == pi_val}")

        errors = {}
        T_min = None

        # Compute contributions one-by-one, accumulating
        # For each T threshold, report the accumulated error
        cumulative = mpmath.mpf(0)
        prev_n = 0

        print(f"\n  {'T':>6} | {'Approx':>14} | {'Error':>14} | {'< 0.5?':>7} | {'dps':>5} | {'Time':>7}")
        print(f"  {'-'*6}-+-{'-'*14}-+-{'-'*14}-+-{'-'*7}-+-{'-'*5}-+-{'-'*7}")

        for T in T_values:
            if T > len(gamma_strs):
                break

            t1 = time.time()

            for k in range(prev_n, T):
                gamma_val = float(gamma_strs[k])
                # Set precision for THIS zero
                dps = set_precision_for(x, gamma_val)
                log_x_hp = mpmath.log(mpmath.mpf(x))
                contrib = compute_zero_contrib(log_x_hp, gamma_strs[k])
                cumulative += contrib

            prev_n = T
            dps_used = mpmath.mp.dps

            # Reset to moderate precision for final arithmetic
            mpmath.mp.dps = 50
            approx_val = float(R_val - cumulative)
            err = abs(approx_val - pi_val)
            dt = time.time() - t1

            exact = err < 0.5
            marker = "YES" if exact else "no"
            if exact and T_min is None:
                T_min = T
                marker = "YES *"

            errors[T] = err
            print(f"  {T:>6} | {approx_val:>14.6f} | {err:>14.6f} | {marker:>7} | {dps_used:>5} | {dt:>6.2f}s")

        elapsed = time.time() - t0
        T_max_tested = max(T for T in T_values if T <= len(gamma_strs))
        print(f"\n  T_min for exact rounding: {T_min if T_min else '>' + str(T_max_tested)}")
        print(f"  Total time: {elapsed:.1f}s")

        results[x] = {
            'pi': pi_val,
            'R': R_float,
            'errors': errors,
            'T_min': T_min,
            'R_error': abs(err_R_only),
        }

    return results


def analyze_lucky_fraction(N_values):
    """
    For what fraction of x in [N, 2N] is round(R(x)) = pi(x)?
    """
    print(f"\n{'='*80}")
    print(f"'LUCKY' FRACTION: How often is round(R(x)) = pi(x)?")
    print(f"{'='*80}")

    mpmath.mp.dps = 30  # Enough for R(x) with small x

    results = {}
    for N in N_values:
        lo, hi = N, 2 * N
        t0 = time.time()

        sieve = bytearray(b'\x01') * (hi + 1)
        sieve[0] = sieve[1] = 0
        for i in range(2, int(hi**0.5) + 1):
            if sieve[i]:
                sieve[i*i::i] = bytearray(len(sieve[i*i::i]))

        pi_cum = [0] * (hi + 1)
        for i in range(1, hi + 1):
            pi_cum[i] = pi_cum[i - 1] + sieve[i]

        correct = 0
        total = hi - lo + 1
        max_err = 0.0
        errors_list = []

        for x in range(lo, hi + 1):
            pi_val = pi_cum[x]
            log_x = mpmath.log(mpmath.mpf(x))
            R_val = float(R_gram_complex(log_x, n_terms=150))
            err = abs(R_val - pi_val)
            if round(R_val) == pi_val:
                correct += 1
            max_err = max(max_err, err)
            errors_list.append(err)

        frac = correct / total
        elapsed = time.time() - t0

        print(f"\n  N = {N:,}, range [{lo:,}, {hi:,}], time: {elapsed:.1f}s")
        print(f"  Correct: {correct}/{total} = {frac:.4f} ({frac*100:.1f}%)")
        print(f"  Max |R(x) - pi(x)|: {max_err:.4f}")
        print(f"  Mean |R(x) - pi(x)|: {np.mean(errors_list):.4f}")
        print(f"  Median |R(x) - pi(x)|: {np.median(errors_list):.4f}")

        print(f"  Error distribution:")
        for threshold in [0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]:
            cnt = sum(1 for e in errors_list if e < threshold)
            print(f"    |error| < {threshold:>4.1f}: {cnt:>6}/{total} ({cnt/total*100:5.1f}%)")

        results[N] = {'frac': frac, 'max_err': max_err, 'mean_err': np.mean(errors_list)}

    return results


def print_scaling_summary(trunc_results, lucky_results):
    """Print final scaling analysis."""
    print(f"\n{'='*80}")
    print(f"SCALING SUMMARY")
    print(f"{'='*80}")

    print(f"\n  Truncation (how many zeros for exact rounding):")
    print(f"  {'x':>12} | {'pi(x)':>10} | {'|R-pi|':>10} | {'T_min':>8} | {'sqrt(x)':>10}")
    print(f"  {'-'*12}-+-{'-'*10}-+-{'-'*10}-+-{'-'*8}-+-{'-'*10}")
    for x in sorted(trunc_results.keys()):
        r = trunc_results[x]
        T_min = r['T_min']
        T_str = str(T_min) if T_min else ">1000"
        print(f"  {x:>12,} | {r['pi']:>10,} | {r['R_error']:>10.3f} | {T_str:>8} | {x**0.5:>10.1f}")

    # Error at fixed T
    print(f"\n  Error at fixed T:")
    T_show = [10, 50, 100, 200, 500, 1000]
    hdr = f"  {'x':>12} |"
    for T in T_show:
        hdr += f" T={T:>4} |"
    print(hdr)
    for x in sorted(trunc_results.keys()):
        r = trunc_results[x]
        line = f"  {x:>12,} |"
        for T in T_show:
            if T in r['errors']:
                e = r['errors'][T]
                if e < 0.5:
                    line += f" {e:>5.2f}*|"
                else:
                    line += f" {e:>5.2f} |"
            else:
                line += f"   N/A |"
        print(line)

    print(f"\n  Lucky fraction (round(R(x))=pi(x)):")
    print(f"  {'N':>10} | {'Lucky %':>10} | {'Max Err':>10} | {'Mean Err':>10} | {'~sqrt(N)/ln(N)':>14}")
    print(f"  {'-'*10}-+-{'-'*10}-+-{'-'*10}-+-{'-'*10}-+-{'-'*14}")
    for N in sorted(lucky_results.keys()):
        r = lucky_results[N]
        expected_err = N**0.5 / np.log(N)
        print(f"  {N:>10,} | {r['frac']*100:>9.1f}% | {r['max_err']:>10.3f} | {r['mean_err']:>10.3f} | {expected_err:>14.2f}")


def main():
    print("=" * 80)
    print("ZETA ZERO TRUNCATION ANALYSIS FOR EXPLICIT FORMULA")
    print("How many zeros are needed for pi(x) = round(R(x) - sum_rho R(x^rho))?")
    print("=" * 80)

    gamma_strs = load_zeros(1000)

    # Use fewer, smaller T values since high-gamma zeros require enormous precision
    T_values = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]

    # Start with small x where 1000 zeros might suffice
    # For x=100, sqrt(x)~10, so maybe 10-20 zeros enough
    # For x=1000, sqrt(x)~31.6
    x_values = [100, 500, 1000, 5000, 10000]

    # Part 1: Truncation analysis
    trunc_results = analyze_truncation(x_values, T_values, gamma_strs)

    # Part 2: Lucky fraction
    print("\n\nComputing lucky fractions...")
    lucky_N = [50, 100, 500, 1000, 5000, 10000, 50000]
    lucky_results = analyze_lucky_fraction(lucky_N)

    # Part 3: Summary
    print_scaling_summary(trunc_results, lucky_results)

    print(f"\n{'='*80}")
    print("CONCLUSIONS")
    print(f"{'='*80}")
    print("1. The lucky fraction (zero zeros needed) decreases as x grows.")
    print("   This is because |R(x)-pi(x)| ~ sqrt(x)/log(x) grows unboundedly.")
    print("2. The number of zeros needed for exact rounding reveals the true scaling.")
    print("3. If T_min ~ C*sqrt(x), the explicit formula inherently requires O(sqrt(x)) work.")
    print("4. Even adaptive approaches cannot avoid this if the error is spread across")
    print("   many independent oscillatory contributions.")


if __name__ == "__main__":
    main()
