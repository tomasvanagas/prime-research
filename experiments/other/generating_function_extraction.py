"""
Experiment A: Extracting pi(x) from the Prime Zeta Function P(s)

The prime zeta function: P(s) = sum_{p prime} p^{-s}
Related to Riemann zeta via: ln(zeta(s)) = sum_{k=1}^inf P(ks)/k
Mobius inversion: P(s) = sum_{k=1}^inf mu(k)/k * ln(zeta(ks))

Question: Can we extract pi(x) from P(s) values at specific points
WITHOUT reducing to the standard explicit formula (which requires O(sqrt(x)) zeros)?

Approach:
1. Compute P(s) via Mobius inversion using known zeta values
2. Try to extract pi(x) via:
   a) Perron's formula applied to sum_{p} p^{-s} (Dirichlet series for pi)
   b) Evaluating P(s) at multiple s values and inverting
   c) Finite differences / Newton series approach

Date: 2026-04-04
"""

import numpy as np
from mpmath import mp, mpf, zeta, log, pi, fac, gamma, exp, inf
import sympy
from sympy import primepi, isprime, nextprime, mobius as sym_mobius
import time

mp.dps = 50  # high precision

def true_pi(x):
    """Ground truth pi(x) using sympy."""
    return int(primepi(x))

def compute_prime_zeta_mobius(s, K_max=20):
    """
    P(s) = sum_{k=1}^{K_max} mu(k)/k * ln(zeta(ks))

    Uses Mobius inversion of ln(zeta(s)) = sum_k P(ks)/k.
    Only needs zeta evaluations at integer multiples of s.
    """
    result = mpf(0)
    for k in range(1, K_max + 1):
        mu_k = int(sym_mobius(k))
        if mu_k == 0:
            continue
        z = zeta(k * s)
        result += mpf(mu_k) / k * log(z)
    return result

def compute_prime_zeta_direct(s, x_max=10000):
    """Direct computation of P(s) = sum_{p <= x_max} p^{-s} for verification."""
    result = mpf(0)
    p = 2
    while p <= x_max:
        result += mpf(p) ** (-s)
        p = int(nextprime(p))
    return result

# ============================================================
# Test 1: Verify Mobius inversion gives correct P(s) values
# ============================================================
print("=" * 70)
print("TEST 1: Verifying P(s) via Mobius inversion")
print("=" * 70)

for s_val in [2, 3, 4, 5, 6, 8, 10]:
    P_mobius = compute_prime_zeta_mobius(s_val, K_max=30)
    P_direct = compute_prime_zeta_direct(s_val, x_max=100000)
    diff = abs(P_mobius - P_direct)
    print(f"s={s_val:2d}: P_mob={float(P_mobius):.12f}  P_dir={float(P_direct):.12f}  diff={float(diff):.2e}")

# ============================================================
# Test 2: Can we extract pi(x) from P(s) at multiple points?
# ============================================================
print("\n" + "=" * 70)
print("TEST 2: Extracting pi(x) from P(s) values")
print("=" * 70)
print("""
Idea: P(s) = sum_{p prime} p^{-s}
If we know P(s) for s = s_1, ..., s_K, we have K equations in the unknowns {p}.
But pi(x) = count of primes <= x, so we need:
  pi(x) = #{p : p <= x} = sum_{p prime} 1_{p <= x}

Perron's formula: pi(x) = (1/2pi i) integral P(s) * x^s / s ds

This IS equivalent to the standard approach -- let's verify numerically.
""")

def perron_integral_prime_zeta(x, T=50, sigma=2.0, N_points=1000):
    """
    Attempt: pi(x) via Perron's formula on the prime zeta function.
    pi(x) = (1/2*pi*i) * integral_{sigma-iT}^{sigma+iT} P(s) * x^s / s ds

    But P(s) has a natural boundary at Re(s)=0, so this integral
    is MUCH harder than the standard one with -zeta'/zeta.
    """
    dt = 2 * T / N_points
    result = mpf(0)

    for j in range(N_points):
        t = -T + (j + 0.5) * dt
        s = sigma + 1j * t
        # P(s) via Mobius inversion
        Ps = compute_prime_zeta_mobius(complex(s.real, s.imag), K_max=15)
        integrand = Ps * mpf(x) ** s / s
        result += integrand * dt

    result = result / (2 * pi * 1j)
    return float(result.real)

# This is expensive, test small values
print("\nPerron integral on P(s) -- testing small x:")
for x in [10, 20, 50, 100]:
    pi_true = true_pi(x)
    try:
        pi_perron = perron_integral_prime_zeta(x, T=30, sigma=2.0, N_points=200)
        print(f"x={x:5d}: pi(x)={pi_true}, Perron estimate={pi_perron:.2f}, error={abs(pi_perron - pi_true):.2f}")
    except Exception as e:
        print(f"x={x:5d}: pi(x)={pi_true}, Perron FAILED: {e}")

# ============================================================
# Test 3: Newton series / finite differences on P(s)
# ============================================================
print("\n" + "=" * 70)
print("TEST 3: Finite difference approach on P(s)")
print("=" * 70)
print("""
Idea: The forward differences Delta^k P(s_0) at s_0 encode information
about the distribution of primes. Can we extract pi(x)?

Delta^0 P(s) = P(s) = sum_p p^{-s}
Delta^1 P(s) = P(s+1) - P(s) = sum_p p^{-s}(p^{-1} - 1)
Delta^k P(s) = sum_p p^{-s} * (p^{-1} - 1)^k

So Delta^k P(s) / (-1)^k = sum_p p^{-s} * (1 - p^{-1})^k

As k -> inf, (1 - p^{-1})^k -> 0 for all p, but largest primes decay slowest.
This does NOT give a way to count primes up to x.
""")

# Compute forward differences of P at integer points
s0 = 2
P_values = []
for j in range(20):
    P_values.append(float(compute_prime_zeta_mobius(s0 + j, K_max=30)))

# Forward differences
diffs = [P_values[:]]
for order in range(1, 15):
    new_diffs = []
    for j in range(len(diffs[-1]) - 1):
        new_diffs.append(diffs[-1][j+1] - diffs[-1][j])
    diffs.append(new_diffs)

print(f"\nP(s) at s = {s0}, {s0+1}, ..., {s0+9}:")
for j in range(10):
    print(f"  P({s0+j}) = {P_values[j]:.12f}")

print(f"\nForward differences Delta^k P({s0}):")
for k in range(10):
    if diffs[k]:
        print(f"  Delta^{k} P({s0}) = {diffs[k][0]:.12e}")

# ============================================================
# Test 4: P(s) partial sums and their relationship to pi(x)
# ============================================================
print("\n" + "=" * 70)
print("TEST 4: Relationship between P(s) partial sums and pi(x)")
print("=" * 70)
print("""
Key test: How many P(s) evaluations do we need to determine pi(x)?

If P_N(s) = sum_{p <= N} p^{-s}, then:
  P_N(s) - P_{N-1}(s) = N^{-s} if N is prime, 0 otherwise.

So recovering pi(x) from P(s) values is equivalent to recovering
individual terms N^{-s}, which requires solving an exponential sum problem.

For K evaluation points s_1,...,s_K, the system is:
  P(s_j) = sum_{p <= x} p^{-s_j} for j=1,...,K

This is a Vandermonde-like system with nodes p^{-1}.
Condition number grows exponentially with pi(x).
""")

# Numerical test: try to recover primes from P(s) at multiple points
from numpy.linalg import lstsq, cond

x_test = 50
primes_up_to_x = [p for p in range(2, x_test+1) if isprime(p)]
pi_x = len(primes_up_to_x)
print(f"\nx={x_test}, pi(x)={pi_x}, primes={primes_up_to_x}")

# Build the system: evaluate P(s_j) for s_j = 2, 3, ..., 2+K
K_values = [5, 10, 15, 20, 30, 50]
candidates = list(range(2, x_test + 1))  # all integers 2..x as potential primes

for K in K_values:
    s_points = [2 + j * 0.5 for j in range(K)]

    # Matrix: A[j,i] = candidates[i]^{-s_j}
    A = np.zeros((K, len(candidates)))
    for j, s in enumerate(s_points):
        for i, n in enumerate(candidates):
            A[j, i] = n ** (-s)

    # RHS: P(s_j) computed from actual primes
    b = np.zeros(K)
    for j, s in enumerate(s_points):
        b[j] = sum(p ** (-s) for p in primes_up_to_x)

    # Condition number
    cn = cond(A)

    # Try to solve for indicator vector (1 if prime, 0 if not)
    x_sol, residuals, rank, sv = lstsq(A, b, rcond=None)

    # Round to 0/1 and check
    x_rounded = (x_sol > 0.5).astype(int)
    recovered_primes = [candidates[i] for i in range(len(candidates)) if x_rounded[i] == 1]
    correct = (recovered_primes == primes_up_to_x)

    print(f"K={K:3d} evals: cond={cn:.2e}, rank={rank}, recovered={'CORRECT' if correct else 'WRONG'} "
          f"(got {len(recovered_primes)} primes)")
    if not correct and K <= 15:
        print(f"  Expected: {primes_up_to_x}")
        print(f"  Got:      {recovered_primes}")

# ============================================================
# Test 5: Does this path differ from explicit formula?
# ============================================================
print("\n" + "=" * 70)
print("TEST 5: Computational path comparison")
print("=" * 70)
print("""
ANALYSIS: The prime zeta approach P(s) = sum mu(k)/k * ln(zeta(ks))

To extract pi(x) from P(s), we need Perron's formula:
  pi(x) = (1/2pi*i) * int P(s) * x^s / s ds

This integral runs over a vertical line in the complex plane.
The poles of P(s) come from the Mobius inversion:
  P(s) = ln(zeta(s)) - P(2s)/2 - ...

The poles of ln(zeta(s)) are at s=1 and at the zeta zeros rho.
After Mobius inversion, the poles of P(s) are at:
  s = 1, s = rho, s = 1/2, s = rho/2, s = 1/3, s = rho/3, ...

So the residues involve ALL zeta zeros (at rho, rho/2, rho/3, ...).
This is STRICTLY WORSE than the standard explicit formula, which only
has zeros at rho.

VERDICT: The P(s) path gives a MORE COMPLICATED sum over zeta zeros,
not a simpler one. It is EQUIVALENT (failure mode E) but with worse constants.
""")

# Verify: count the effective number of "terms" needed
print("Number of terms in Perron residue expansion:")
print("  Standard explicit formula: sum over rho (zeta zeros)")
print("  Prime zeta Perron: sum over rho, rho/2, rho/3, ... (MORE terms)")
print()
print("Specifically, for T zeros in standard formula,")
print("prime zeta approach needs ~T * sum_{k=1}^{K} 1 = T*K terms.")
print()
print("CONCLUSION: Generating function P(s) approach is EQUIVALENT to")
print("explicit formula, with WORSE computational properties (more poles).")
print("Failure mode: E (Equivalence)")
