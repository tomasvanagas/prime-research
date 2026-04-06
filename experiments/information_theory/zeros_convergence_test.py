"""
Zero Convergence Test
======================
For p(1000) = 7919, add zeta zeros one at a time and watch
pi(x) converge to the correct value. How many zeros until exact?

Uses mpmath for reliable complex arithmetic.

Session 41.
"""

import math
import sys

try:
    import mpmath
    mpmath.mp.dps = 25
except ImportError:
    print("ERROR: mpmath required. Install with: pip install mpmath")
    sys.exit(1)

def sieve_primes(limit):
    is_prime = bytearray([1]) * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return [i for i in range(2, limit + 1) if is_prime[i]]

def R_real(x):
    """Riemann R(x) for real x using mpmath."""
    x = mpmath.mpf(x)
    result = mpmath.mpf(0)
    for n in range(1, 100):
        mu_n = mobius(n)
        if mu_n == 0:
            continue
        xn = x ** (mpmath.mpf(1) / n)
        if xn <= 1:
            break
        result += mpmath.mpf(mu_n) / n * mpmath.li(xn)
    return float(result)

def R_zero_contrib(x, gamma):
    """Contribution of zero pair rho, conj(rho) to the explicit formula.
    Returns 2*Re(R(x^rho)) where rho = 1/2 + i*gamma."""
    x = mpmath.mpf(x)
    gamma = mpmath.mpf(gamma)
    rho = mpmath.mpc(0.5, gamma)

    result = mpmath.mpc(0)
    for n in range(1, 50):
        mu_n = mobius(n)
        if mu_n == 0:
            continue
        # x^{rho/n}
        rho_n = rho / n
        xrn = mpmath.power(x, rho_n)
        if abs(xrn) < 1.001:
            break
        li_val = mpmath.li(xrn)
        result += mpmath.mpf(mu_n) / n * li_val

    # Return 2*Re(result) for the conjugate pair
    return 2 * float(result.real)

_mobius_cache = {}
def mobius(n):
    if n in _mobius_cache:
        return _mobius_cache[n]
    if n == 1:
        _mobius_cache[1] = 1
        return 1
    # Factor n
    result = 1
    temp = n
    for p in range(2, int(n**0.5) + 2):
        if temp % p == 0:
            temp //= p
            if temp % p == 0:
                _mobius_cache[n] = 0
                return 0
            result *= -1
    if temp > 1:
        result *= -1
    _mobius_cache[n] = result
    return result

def load_zeros(filename):
    zeros = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                try:
                    zeros.append(float(line))
                except ValueError:
                    pass
    return zeros

def main():
    print("=" * 70)
    print("ZERO CONVERGENCE TEST: How many zeros for exact p(n)?")
    print("=" * 70)

    zeros = load_zeros("data/zeta_zeros_1000.txt")
    print(f"Loaded {len(zeros)} zeta zeros")

    primes = sieve_primes(200000)

    test_cases = [
        (100, primes[99]),
        (500, primes[499]),
        (1000, primes[999]),
        (2000, primes[1999]),
        (5000, primes[4999]),
        (10000, primes[9999]),
    ]

    summary = []

    for n, actual_p in test_cases:
        print(f"\n{'='*70}")
        print(f"Finding p({n}) = {actual_p}")
        print(f"{'='*70}")

        x = actual_p + 0.5

        # Smooth part
        R_x = R_real(x)
        const_term = 1.0 / math.log(2)  # ≈ 1.4427
        pi_base = R_x - const_term

        print(f"  R({x:.1f}) = {R_x:.6f}")
        print(f"  R(x) - 1/ln2 = {pi_base:.6f} (target: {n})")
        print(f"  Error without zeros: {pi_base - n:+.6f}")

        # Add zeros incrementally
        checkpoints = sorted(set(
            list(range(0, 11)) + [15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 300, 500, 750, 1000]
        ))
        checkpoints = [c for c in checkpoints if c <= len(zeros)]

        corrections = []  # store individual corrections for incremental sum
        total_correction = 0.0
        found_K = None
        prev_k = 0

        print(f"\n  {'K':>6} {'Total corr':>12} {'π(x) est':>12} {'Error':>10} {'Round':>6} {'OK?':>4}")

        for K in checkpoints:
            for k in range(prev_k, K):
                contrib = R_zero_contrib(x, zeros[k])
                total_correction += contrib
            prev_k = K

            pi_est = pi_base - total_correction
            error = pi_est - n
            rounded = round(pi_est)
            exact = "YES" if rounded == n else ""

            if K <= 10 or K in [15,20,25,30,40,50,75,100,150,200,300,500,750,1000] or (found_K is None and exact == "YES"):
                print(f"  {K:6d} {total_correction:+12.6f} {pi_est:12.6f} {error:+10.6f} {rounded:6d} {exact:>4}")

            if found_K is None and exact == "YES":
                found_K = K

        if found_K is not None:
            summary.append((n, actual_p, found_K))
            print(f"\n  ★ EXACT after {found_K} zeros")
        else:
            summary.append((n, actual_p, ">1000"))
            print(f"\n  ✗ NOT exact with {len(zeros)} zeros (error: {error:+.4f})")

    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY: Zeros needed for exact p(n)")
    print(f"{'='*70}")
    print(f"\n  {'n':>8} {'p(n)':>10} {'√p(n)':>8} {'Zeros needed':>13}")
    for n, p, k in summary:
        sqrtp = math.sqrt(p)
        print(f"  {n:8d} {p:10d} {sqrtp:8.1f} {str(k):>13}")

if __name__ == "__main__":
    main()
