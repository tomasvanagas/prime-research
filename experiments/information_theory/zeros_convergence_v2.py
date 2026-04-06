"""
Zero Convergence Test v2
=========================
Use the standard explicit formula with proper li(x^rho) computation.

pi(x) = R(x) - sum_{gamma>0} 2*Re(li(x^{1/2+i*gamma})) - 1/ln(2) + small terms

Session 41.
"""

import math
import mpmath

mpmath.mp.dps = 30

def sieve_primes(limit):
    is_prime = bytearray([1]) * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return [i for i in range(2, limit + 1) if is_prime[i]]

def R_real(x):
    """Gram series for R(x)."""
    x = mpmath.mpf(x)
    ln_x = mpmath.log(x)
    result = mpmath.mpf(1)
    term = mpmath.mpf(1)
    for n in range(1, 200):
        term *= ln_x / n
        zeta_contrib = mpmath.mpf(1) / (n * mpmath.zeta(n + 1))
        result += term * zeta_contrib
        if abs(term * zeta_contrib) < mpmath.mpf(10)**(-20):
            break
    return float(result)

def li_complex(rho, ln_x):
    """Compute li(x^rho) = Ei(rho * ln(x)) using mpmath.
    rho is complex, ln_x is real."""
    z = rho * ln_x
    # li(x^rho) = Ei(rho * ln(x))
    return mpmath.ei(z)

def load_zeros(filename):
    zeros = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line:
                try:
                    zeros.append(float(line))
                except ValueError:
                    pass
    return zeros

def main():
    print("=" * 70)
    print("ZERO CONVERGENCE TEST v2")
    print("=" * 70)

    zeros = load_zeros("data/zeta_zeros_1000.txt")
    print(f"Loaded {len(zeros)} zeros (γ_1={zeros[0]:.4f} to γ_{len(zeros)}={zeros[-1]:.4f})")

    primes = sieve_primes(200000)

    test_cases = [
        (100, primes[99]),
        (500, primes[499]),
        (1000, primes[999]),
        (5000, primes[4999]),
        (10000, primes[9999]),
    ]

    summary = []

    for n, actual_p in test_cases:
        print(f"\n{'='*70}")
        print(f"p({n}) = {actual_p}, testing at x = {actual_p} + 0.5")
        print(f"{'='*70}")

        x = mpmath.mpf(actual_p) + mpmath.mpf(0.5)
        ln_x = mpmath.log(x)

        # Smooth part: R(x) via Gram series
        Rx = R_real(float(x))

        # Constant correction
        const = 1.0 / math.log(2)

        pi_base = Rx - const
        print(f"  R(x) = {Rx:.6f}")
        print(f"  R(x) - 1/ln2 = {pi_base:.6f} (target π(x) = {n})")

        # Add zero corrections: -2*Re(li(x^rho)) per zero pair
        checkpoints = sorted(set(
            list(range(0, 11)) + [15, 20, 25, 30, 40, 50, 75, 100,
            150, 200, 300, 500, 750, 1000]
        ))
        checkpoints = [c for c in checkpoints if c <= len(zeros)]

        total_zero_correction = mpmath.mpf(0)
        prev_k = 0
        found_K = None

        print(f"\n  {'K':>6} {'Zero corr':>12} {'π(x) est':>12} {'Error':>10} {'Round':>6} {'OK?':>4}")

        for K in checkpoints:
            for k in range(prev_k, K):
                gamma = mpmath.mpf(zeros[k])
                rho = mpmath.mpc(0.5, gamma)
                # li(x^rho) = Ei(rho * ln_x)
                li_val = li_complex(rho, ln_x)
                # Conjugate pair: 2*Re(li(x^rho))
                total_zero_correction += 2 * mpmath.re(li_val)
            prev_k = K

            pi_est = float(pi_base - total_zero_correction)
            error = pi_est - n
            rounded = round(pi_est)
            exact = "YES" if rounded == n else ""

            show = K <= 10 or K in [15,20,25,30,40,50,75,100,150,200,300,500,750,1000]
            if show or (found_K is None and exact == "YES"):
                print(f"  {K:6d} {float(total_zero_correction):+12.4f} {pi_est:12.4f} {error:+10.4f} {rounded:6d} {exact:>4}")

            if found_K is None and exact == "YES":
                found_K = K

        if found_K is not None:
            summary.append((n, actual_p, found_K))
            print(f"\n  ★ EXACT after {found_K} zeros!")
        else:
            summary.append((n, actual_p, f">{len(zeros)}"))
            print(f"\n  ✗ Not exact with {len(zeros)} zeros (error: {error:+.4f})")

    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"\n  {'n':>8} {'p(n)':>10} {'√p(n)':>8} {'Zeros needed':>13}")
    for n, p, k in summary:
        print(f"  {n:8d} {p:10d} {math.sqrt(p):8.1f} {str(k):>13}")

if __name__ == "__main__":
    main()
