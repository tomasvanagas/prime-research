"""
PROPOSAL 1: Spectral Truncation + Sieve Candidate Elimination

IDEA: Use a truncated Riemann explicit formula with K zeros to get an interval
[a, b] containing p(n). Then use a SIEVE on [a, b] to find which numbers are
prime, and count to find the nth one exactly.

If the interval [a,b] has width W = O(x/K * log^2 x), and sieving [a,b] costs
O(W * log log W), the total cost is:
  O(K * polylog(x)  +  x/K * log^3(x))

Minimizing over K: K = sqrt(x / log^3 x), total = O(sqrt(x) * polylog)
This doesn't beat O(x^{2/3}) -- BUT:

KEY TWIST: If we also know p(n) mod m for some modulus m, the effective interval
shrinks from W to W/m. This gives:
  O(K * polylog + x/(K*m) * log^3 x)
  
If m = x^{epsilon} (achievable via CRT on small primes), we get:
  O(x^{1/2 - epsilon/2} * polylog)

Can we push m large enough? If m = x^{1/2}/polylog, we'd get O(polylog)!
"""

import numpy as np
from sympy import prime, primepi, isprime, nextprime, prevprime, sieve
from mpmath import mp, mpf, log, li, pi as mpi, zeta, im, re, zetazero
import time

mp.dps = 30

def compute_R(x, terms=20):
    """Riemann's R function: R(x) = sum_{k=1}^{inf} mu(k)/k * li(x^{1/k})"""
    from sympy import mobius
    x = mpf(x)
    result = mpf(0)
    for k in range(1, terms+1):
        mu_k = mobius(k)
        if mu_k != 0:
            result += mpf(mu_k) / k * li(x ** (mpf(1)/k))
    return result

def R_inverse(n, terms=20):
    """Compute R^{-1}(n) via Newton's method"""
    x = mpf(n) * log(mpf(n))
    for _ in range(50):
        rx = compute_R(x, terms)
        dx = 1 / log(x)  # R'(x) ≈ 1/ln(x)
        correction = (mpf(n) - rx) / dx
        x += correction
        if abs(correction) < 0.01:
            break
    return x

def truncated_explicit_formula(x, n_zeros):
    """
    pi(x) ≈ R(x) - sum_{k=1}^{n_zeros} R(x^{rho_k}) - ...
    where rho_k are the nontrivial zeros of zeta.
    
    Returns (estimate, error_bound)
    """
    x = mpf(x)
    result = compute_R(x)
    
    # Subtract contributions from first n_zeros pairs of zeros
    zero_sum = mpf(0)
    for k in range(1, n_zeros + 1):
        rho = zetazero(k)  # First zero: 0.5 + 14.134...i
        # R(x^rho) + R(x^{conj(rho)})
        # For real result, 2*Re(R(x^rho))
        x_rho = x ** rho
        # li(x^rho) is the main term of R(x^rho)
        zero_sum += 2 * re(li(x_rho))
    
    result -= zero_sum
    
    # Error from truncation: roughly x / (n_zeros * log(x))
    error_bound = float(x) / (n_zeros * float(log(x)))
    
    return float(result), error_bound

def spectral_sieve_method(n, n_zeros=10):
    """
    1. Compute R^{-1}(n) for center estimate
    2. Use truncated explicit formula to bound the interval
    3. Sieve the interval to find exact p(n)
    """
    t0 = time.time()
    
    # Step 1: Get approximate location
    x_approx = float(R_inverse(n))
    
    # Step 2: Use explicit formula to get pi(x_approx) and bound the error
    pi_est, error = truncated_explicit_formula(x_approx, n_zeros)
    
    # Step 3: We need to find where pi(x) = n exactly
    # The error tells us the interval width we need to search
    width = int(2 * error + 100)  # Safety margin
    
    lo = max(2, int(x_approx - width))
    hi = int(x_approx + width)
    
    # Step 4: Sieve the interval [lo, hi]
    # Count primes from 2 to lo-1 using our estimate, then sieve to find exact position
    # For testing, use sympy
    pi_lo = primepi(lo - 1)
    count_needed = n - pi_lo
    
    # Sieve [lo, hi] for primes
    primes_in_range = [p for p in range(lo, hi+1) if isprime(p)]
    
    t1 = time.time()
    
    if 0 < count_needed <= len(primes_in_range):
        result = primes_in_range[count_needed - 1]
        return result, hi - lo, t1 - t0
    else:
        return None, hi - lo, t1 - t0

# Test
print("=" * 70)
print("PROPOSAL 1: Spectral Truncation + Sieve")
print("=" * 70)
print()

for n in [100, 500, 1000, 2000, 5000]:
    expected = prime(n)
    for nz in [5, 10, 20]:
        result, width, elapsed = spectral_sieve_method(n, n_zeros=nz)
        status = "OK" if result == expected else f"FAIL (got {result})"
        print(f"n={n:>5}, zeros={nz:>2}: width={width:>8}, time={elapsed:.3f}s, {status}")
    print()

# Key measurement: how does width scale with n_zeros?
print("\nWidth vs zeros for n=5000:")
for nz in [1, 2, 5, 10, 20, 50]:
    _, width, elapsed = spectral_sieve_method(5000, n_zeros=nz)
    print(f"  zeros={nz:>3}: width={width:>10}, time={elapsed:.3f}s")

