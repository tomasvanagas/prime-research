"""
Zero Scaling Analysis: How does the number of zeros needed scale with x?
=========================================================================
Session 11 experiment.

For exact pi(x) via explicit formula, we need the truncation error < 0.5.
The error at N zeros is approximately |E(x,N)| ~ C * sqrt(x) * ln(N) / N.

This experiment measures the EXACT scaling empirically:
- For each x, find minimum N where round(S(N)) = pi(x) (for plain and accelerated)
- Fit N_min(x) to determine the scaling law
- Compare: N ~ sqrt(x)? N ~ x^{1/3}? N ~ log(x)^k?

This determines whether ANY constant-factor improvement from acceleration
translates to an asymptotic improvement.
"""

import mpmath
import math
import os
import sys
import time
import numpy as np

mpmath.mp.dps = 50

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
from riemann_explicit import R_function, R_at_rho, MOBIUS


def load_zeros(count):
    data_dir = os.path.join(SCRIPT_DIR, '..', '..', 'data')
    for try_count in [1000, 500, 300, 200]:
        fname = os.path.join(data_dir, f"zeta_zeros_{try_count}.txt")
        if os.path.exists(fname):
            with open(fname) as f:
                zeros = [line.strip() for line in f if line.strip()]
            if len(zeros) >= count:
                return zeros[:count]
    raise FileNotFoundError(f"Need {count} zeros")


def richardson_2pt(S1, S2, N1, N2, alpha=1):
    """2-point Richardson extrapolation."""
    n1a = N1 ** alpha
    n2a = N2 ** alpha
    return (n2a * S2 - n1a * S1) / (n2a - n1a)


def richardson_3pt(S1, S2, S3, N1, N2, N3, alpha=1):
    """3-point Richardson."""
    R12 = richardson_2pt(S1, S2, N1, N2, alpha)
    R23 = richardson_2pt(S2, S3, N2, N3, alpha)
    n2a = N2 ** (alpha + 1)
    n3a = N3 ** (alpha + 1)
    return (n3a * R23 - n2a * R12) / (n3a - n2a)


def find_min_N(exact, Rx, terms, method='plain', max_N=None):
    """Find minimum N where method gives exact rounding."""
    if max_N is None:
        max_N = len(terms)

    # Binary search for minimum N
    partial = 0.0
    partial_sums = []
    for t in terms[:max_N]:
        partial += t
        partial_sums.append(Rx - partial)

    if method == 'plain':
        for N in range(1, max_N + 1):
            if round(partial_sums[N-1]) == exact:
                # Check stability: verify next 5 values also round correctly
                stable = True
                for k in range(1, min(6, max_N - N + 1)):
                    if round(partial_sums[N-1+k]) != exact:
                        stable = False
                        break
                if stable:
                    return N
        return None

    elif method == 'richardson1':
        # Richardson order 1 with N//2 and N
        for N in range(10, max_N + 1, 2):
            N1 = N // 2
            N2 = N
            val = richardson_2pt(partial_sums[N1-1], partial_sums[N2-1], N1, N2)
            if round(val) == exact:
                # Check stability
                stable = True
                for k in range(2, min(12, max_N - N + 1), 2):
                    N1b = (N + k) // 2
                    N2b = N + k
                    if N2b <= max_N:
                        val2 = richardson_2pt(partial_sums[N1b-1], partial_sums[N2b-1], N1b, N2b)
                        if round(val2) != exact:
                            stable = False
                            break
                if stable:
                    return N
        return None

    elif method == 'richardson2':
        # Richardson order 2 with N//3, 2N//3, N
        for N in range(15, max_N + 1, 3):
            N1 = max(5, N // 3)
            N2 = max(N1 + 1, 2 * N // 3)
            N3 = N
            if N1 < N2 < N3 <= max_N:
                try:
                    val = richardson_3pt(
                        partial_sums[N1-1], partial_sums[N2-1], partial_sums[N3-1],
                        N1, N2, N3)
                    if round(val) == exact:
                        return N
                except:
                    pass
        return None

    return None


# Exact values of pi(x) for testing
EXACT_PI = {
    100: 25,
    200: 46,
    500: 95,
    1000: 168,
    2000: 303,
    5000: 669,
    10000: 1229,
    20000: 2262,
    50000: 5133,
    100000: 9592,
    200000: 17984,
    500000: 41538,
    1000000: 78498,
}


def main():
    print("Zero Scaling Analysis")
    print("=" * 80)
    print("Finding minimum N (# of zeta zeros) for exact pi(x) as x grows")
    print()

    max_zeros = 1000
    zeros_str = load_zeros(max_zeros)

    methods = ['plain', 'richardson1', 'richardson2']
    results = {m: [] for m in methods}

    for x_val in sorted(EXACT_PI.keys()):
        exact = EXACT_PI[x_val]
        sqrt_x = math.sqrt(x_val)

        print(f"\nx = {x_val:>10,}  pi(x) = {exact:>6d}  sqrt(x) = {sqrt_x:.1f}")

        # Compute terms
        t0 = time.time()
        x_mpf = mpmath.mpf(x_val)
        Rx = float(R_function(x_mpf))
        terms = []
        for gs in zeros_str[:max_zeros]:
            gamma = mpmath.mpf(gs)
            rho = mpmath.mpc(0.5, gamma)
            val = 2 * float(mpmath.re(R_at_rho(x_mpf, rho)))
            terms.append(val)
        elapsed = time.time() - t0

        for method in methods:
            min_N = find_min_N(exact, Rx, terms, method=method, max_N=max_zeros)
            if min_N is not None:
                ratio = min_N / sqrt_x
                results[method].append((x_val, min_N))
                print(f"  {method:15s}: N_min = {min_N:4d}  "
                      f"N/sqrt(x) = {ratio:.4f}  "
                      f"N/x^(1/3) = {min_N / x_val**(1/3):.4f}  "
                      f"N/ln(x)^2 = {min_N / math.log(x_val)**2:.4f}")
            else:
                print(f"  {method:15s}: NOT EXACT with {max_zeros} zeros")

    # Scaling analysis
    print("\n" + "=" * 80)
    print("SCALING ANALYSIS")
    print("=" * 80)

    for method in methods:
        data = results[method]
        if len(data) < 3:
            print(f"\n{method}: insufficient data for scaling analysis")
            continue

        xs = np.array([d[0] for d in data], dtype=float)
        Ns = np.array([d[1] for d in data], dtype=float)

        # Fit: log(N) = alpha * log(x) + C
        log_x = np.log(xs)
        log_N = np.log(Ns)

        # Use last 2/3 of data points for fit (more reliable at larger x)
        start = len(log_x) // 3
        coeffs = np.polyfit(log_x[start:], log_N[start:], 1)
        alpha = coeffs[0]
        C = coeffs[1]

        print(f"\n{method}:")
        print(f"  Fit: N_min ~ {math.exp(C):.3f} * x^{alpha:.4f}")
        print(f"  For x^{1/2} scaling, alpha should be 0.500")
        print(f"  For x^{1/3} scaling, alpha should be 0.333")
        print(f"  For polylog scaling, alpha should be ~0")
        print(f"  MEASURED alpha = {alpha:.4f}")

        # Extrapolate: at x = 10^100, how many zeros needed?
        N_10_100 = math.exp(C) * (1e100 ** alpha)
        print(f"  Extrapolated N for x=10^100: {N_10_100:.2e}")

        # Also fit N_min / sqrt(x) to see if ratio is constant
        ratios = Ns / np.sqrt(xs)
        print(f"  N/sqrt(x) ratios: {', '.join(f'{r:.4f}' for r in ratios[-5:])}")

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print("""
If alpha ≈ 0.5: Zeros needed scale as sqrt(x) — NO improvement possible via acceleration
If alpha < 0.5: Some acceleration is effective — worth investigating
If alpha ≈ 0: BREAKTHROUGH — polylog zeros suffice

For p(10^100), x ~ 10^102:
  sqrt(x) = 10^51 zeros → 10^51 zero computations → infeasible
  x^{1/3} = 10^34 zeros → still infeasible
  polylog(x) = 100^k zeros → FEASIBLE
""")


if __name__ == "__main__":
    main()
