"""
Session 4: HIGH-PRECISION Möbius inversion of explicit formula

KEY FINDING: Möbius inversion of Π(x) dramatically reduces error from zeta zeros.
With 20 zeros, direct formula gives error ~40, but Möbius inversion gives error ~1.

Hypothesis: With enough zeros and proper computation of li(x^ρ),
can we get error < 0.5 and thus EXACT π(x)?

If so: we have a formula for π(x) using only O(K) zeta zeros,
where K might be polylog(x) instead of sqrt(x).
"""

import mpmath
import math
import time

mpmath.mp.dps = 50  # 50 decimal digits

# Known zeta zeros (imaginary parts), high precision
# We'll compute more with mpmath if needed
ZETA_ZEROS_STR = [
    "14.134725141734693790457251983562",
    "21.022039638771554992628479593897",
    "25.010857580145688763213790992563",
    "30.424876125859513210311897530584",
    "32.935061587739189690662368964075",
    "37.586178158825671257217763480021",
    "40.918719012147495187398126914634",
    "43.327073280914999519496122165406",
    "48.005150881167159727942472749427",
    "49.773832477672302181916784678564",
    "52.970321477714460644147296608881",
    "56.446247697063394804367759476706",
    "59.347044002602353381775021659007",
    "60.831778524609809844259901824524",
    "65.112544048081607691227685668403",
    "67.079810529494174714996758556791",
    "69.546401711173979252926857526554",
    "72.067157674481907581936107079767",
    "75.704690699083933168326916762030",
    "77.144840068874805372682664856305",
    "79.337375020249367922763592877116",
    "82.910380854086030183164837494770",
    "84.735492980517050105735311206827",
    "87.425274613125229406531667850919",
    "88.809111207634465423682348079509",
    "92.491899270558484296259725241810",
    "94.651344040519886966597925815200",
    "95.870634228245309758741029219246",
    "98.831194218193692233324420138622",
    "101.31785100573139122878544794833",
]

ZETA_ZEROS = [mpmath.mpf(z) for z in ZETA_ZEROS_STR]

def li_complex(z):
    """Logarithmic integral li(z) = Ei(ln z) for complex z."""
    if z == 0:
        return mpmath.mpf(0)
    return mpmath.ei(mpmath.log(z))

def Pi_explicit_hp(x, K=30):
    """
    Riemann's Π(x) using the EXACT explicit formula with K zeta zeros.
    Π(x) = li(x) - Σ_ρ li(x^ρ) - ln(2) + ∫_x^∞ dt/(t(t²-1)ln t)

    The integral term is negligible for large x.
    """
    x = mpmath.mpf(x)
    if x < 2:
        return mpmath.mpf(0)

    result = mpmath.li(x) - mpmath.log(2)

    # Subtract contribution from each pair of zeros ρ, ρ̄
    for k in range(min(K, len(ZETA_ZEROS))):
        gamma = ZETA_ZEROS[k]
        rho = mpmath.mpc(mpmath.mpf('0.5'), gamma)
        rho_conj = mpmath.mpc(mpmath.mpf('0.5'), -gamma)

        # li(x^ρ) + li(x^ρ̄) = 2 Re[li(x^ρ)]
        x_rho = mpmath.power(x, rho)
        li_val = mpmath.ei(mpmath.log(x) * rho)  # = li(x^ρ) = Ei(ρ ln x)

        result -= 2 * mpmath.re(li_val)

    # Integral term (very small for x > 2)
    # ∫_x^∞ dt/(t(t²-1)ln t) ≈ 0 for large x
    # For x ≈ 2: about 0.1...
    # We can compute it for accuracy
    if x < 100:
        integral = mpmath.quad(lambda t: 1/(t*(t**2-1)*mpmath.log(t)), [x, mpmath.inf])
        result += integral

    return result

def pi_from_Pi_hp(x, K=30):
    """π(x) via Möbius inversion of Π(x)."""
    x = mpmath.mpf(x)
    result = mpmath.mpf(0)
    max_k = int(mpmath.log(x) / mpmath.log(2)) + 1

    for k in range(1, max_k + 1):
        mu_k = mobius_simple(k)
        if mu_k == 0:
            continue
        xk = mpmath.power(x, mpmath.mpf(1)/k)
        if xk < 2:
            break
        result += mpmath.mpf(mu_k) / k * Pi_explicit_hp(xk, K)

    return result

def mobius_simple(n):
    """Möbius function."""
    if n == 1:
        return 1
    result = 1
    temp = n
    for p in range(2, int(n**0.5) + 1):
        if temp % p == 0:
            temp //= p
            if temp % p == 0:
                return 0
            result *= -1
    if temp > 1:
        result *= -1
    return result

def primes_sieve(n):
    """Simple sieve for verification."""
    sieve = [True] * (n + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(n**0.5) + 1):
        if sieve[i]:
            for j in range(i*i, n+1, i):
                sieve[j] = False
    return [i for i in range(2, n+1) if sieve[i]]

# Generate reference
primes = primes_sieve(200000)
pi_exact = {}  # x -> π(x) lookup
count = 0
for i, p in enumerate(primes):
    pi_exact[p] = i + 1

print("=" * 70)
print("HIGH-PRECISION Möbius Inversion of Explicit Formula")
print(f"Using {len(ZETA_ZEROS)} zeta zeros, {mpmath.mp.dps} digit precision")
print("=" * 70)

# Test with increasing zeros
print("\n--- Accuracy vs number of zeros ---")
test_points = [primes[n-1] for n in [100, 500, 1000, 2000, 5000, 10000]]
test_ns = [100, 500, 1000, 2000, 5000, 10000]

for K in [5, 10, 15, 20, 25, 30]:
    errors = []
    for x, n in zip(test_points, test_ns):
        t0 = time.time()
        pi_val = pi_from_Pi_hp(x, K)
        elapsed = time.time() - t0
        error = float(abs(pi_val - n))
        errors.append(error)

    max_err = max(errors)
    avg_err = sum(errors) / len(errors)
    exact_count = sum(1 for e in errors if e < 0.5)
    print(f"K={K:2d}: max_err={max_err:.4f}, avg_err={avg_err:.4f}, "
          f"exact={exact_count}/{len(errors)}")

# Detailed test for K=30
print(f"\n--- Detailed results with K=30 zeros ---")
correct = 0
total = 0
max_error = 0
for n in [10, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 15000]:
    if n > len(primes):
        break
    x = primes[n-1]
    t0 = time.time()
    pi_val = pi_from_Pi_hp(x, K=30)
    elapsed = time.time() - t0
    error = float(abs(pi_val - n))
    rounded = int(round(float(pi_val)))
    is_exact = (rounded == n)
    if is_exact:
        correct += 1
    total += 1
    max_error = max(max_error, error)
    print(f"  n={n:6d}: π(p(n))={float(pi_val):.6f}, error={error:.6f}, "
          f"exact={'YES' if is_exact else 'NO'}, time={elapsed:.3f}s")

print(f"\nExact: {correct}/{total}")
print(f"Max error: {max_error:.6f}")

# KEY QUESTION: Does error grow with x? Or stay bounded?
print("\n--- Error scaling analysis ---")
print("Does error grow as sqrt(x)/T or stay bounded?")
for n in [100, 1000, 10000]:
    x = primes[n-1]
    pi_val = pi_from_Pi_hp(x, K=30)
    error = float(abs(pi_val - n))
    sqrtx = float(mpmath.sqrt(x))
    T = float(ZETA_ZEROS[-1])  # height of highest zero used
    theoretical_error = sqrtx / T
    print(f"  n={n:6d}, x={x:10d}: error={error:.4f}, "
          f"sqrt(x)/T={theoretical_error:.4f}, ratio={error/theoretical_error:.4f}")
