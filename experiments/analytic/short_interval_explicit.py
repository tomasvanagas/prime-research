"""
Short-interval explicit formula analysis.

Question: For counting primes in [x-W, x+W], how many zeta zeros are needed
for error < 1? Does the scaling match the theoretical K ~ x/W prediction,
or is there unexpected structure?

This is a COST ANALYSIS experiment — testing whether the zero requirement
for short-interval counting differs from the full pi(x) computation.

Key theoretical prediction:
- Full pi(x): needs K ~ sqrt(x) zeros for error < 1
- Short interval width W: needs K ~ x/W zeros (from sinc decay)
- Hybrid optimum: W = sqrt(x) → K = sqrt(x), total cost O(sqrt(x)*polylog)
- Iteration SHOULD NOT help: each round needs ~same zeros

We test this empirically for x up to 10^6.
"""

import numpy as np
from sympy import primepi, prime, isprime, primerange
import warnings
warnings.filterwarnings('ignore')

try:
    import mpmath
    mpmath.mp.dps = 50
    HAS_MPMATH = True
except ImportError:
    HAS_MPMATH = False

def get_zeta_zeros(num_zeros):
    """Get first num_zeros nontrivial zeta zeros (imaginary parts)."""
    if not HAS_MPMATH:
        # Hardcode first 30 zeros for fallback
        known = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
                 37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
                 52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
                 67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
                 79.337375, 82.910381, 84.735493, 87.425275, 88.809111,
                 92.491899, 94.651344, 95.870634, 98.831194, 101.317851]
        return known[:num_zeros]

    zeros = []
    for k in range(1, num_zeros + 1):
        try:
            z = float(mpmath.zetazero(k).imag)
            zeros.append(z)
        except:
            break
    return zeros

def li_rho(x, gamma):
    """Compute the real part of li(x^rho) for rho = 1/2 + i*gamma."""
    if x <= 1:
        return 0.0
    log_x = np.log(x)
    # x^rho = x^{1/2} * exp(i*gamma*log(x))
    sqrt_x = np.sqrt(x)
    phase = gamma * log_x
    # li(x^rho) ≈ x^rho / (rho * log(x)) for large x
    # More precisely: li(x^rho) = Ei(rho*log(x))
    # For numerical stability, use the approximation:
    # Re[li(x^rho)] ≈ (sqrt(x)/log(x)) * (cos(phase)/(0.25+gamma^2))*(0.5*cos(phase)+gamma*sin(phase))
    # Actually let's use the integral representation more carefully

    # Simplified: Re[li(x^{1/2+ig})] ≈ sqrt(x) * [0.5*cos(g*lnx) + g*sin(g*lnx)] / ((0.25+g^2)*ln(x))
    denom = (0.25 + gamma**2) * log_x
    numerator = sqrt_x * (0.5 * np.cos(phase) + gamma * np.sin(phase))
    return numerator / denom

def pi_explicit_formula(x, zeros, use_all=True):
    """Compute pi(x) using the truncated explicit formula."""
    if x < 2:
        return 0

    log_x = np.log(x)

    # Main term: R(x) = sum_{k=1}^{inf} mu(k)/k * li(x^{1/k})
    # Approximate as li(x) - li(x^{1/2})/2 - li(x^{1/3})/3 + ...
    from sympy import mobius, li as sympy_li
    R_x = 0
    for k in range(1, int(log_x / np.log(2)) + 2):
        mu_k = int(mobius(k))
        if mu_k == 0:
            continue
        xk = x ** (1.0/k)
        if xk < 2:
            break
        # li(y) ≈ y/ln(y) for large y
        li_val = float(sympy_li(xk))
        R_x += mu_k / k * li_val

    # Subtract zero contributions
    zero_sum = 0
    for gamma in zeros:
        zero_sum += 2 * li_rho(x, gamma)  # pair: rho and conj(rho)

    return R_x - zero_sum

def pi_short_interval_explicit(a, b, zeros):
    """Count primes in (a, b] using explicit formula difference."""
    return pi_explicit_formula(b, zeros) - pi_explicit_formula(a, zeros)

def count_primes_sieve(a, b):
    """Exact count of primes in (a, b] using sympy."""
    return int(primepi(b)) - int(primepi(a))

def test_short_interval_zero_requirement():
    """Main experiment: for various x and W, how many zeros give error < 1?"""

    print("=" * 80)
    print("SHORT-INTERVAL EXPLICIT FORMULA: ZERO REQUIREMENT ANALYSIS")
    print("=" * 80)

    # Get zeros
    max_zeros = 30
    zeros = get_zeta_zeros(max_zeros)
    print(f"\nLoaded {len(zeros)} zeta zeros")

    test_points = [1000, 5000, 10000, 50000, 100000, 500000]

    print("\n--- PART 1: Full pi(x) error vs number of zeros ---")
    print(f"{'x':>10} {'K zeros':>8} {'pi_approx':>12} {'pi_exact':>10} {'error':>10} {'|err|/sqrt(x)':>14}")
    print("-" * 75)

    for x in test_points:
        exact = int(primepi(x))
        for K in [1, 3, 5, 10, 15, 20, 25, 30]:
            if K > len(zeros):
                break
            approx = pi_explicit_formula(x, zeros[:K])
            err = approx - exact
            ratio = abs(err) / np.sqrt(x) if x > 0 else 0
            print(f"{x:>10} {K:>8} {approx:>12.2f} {exact:>10} {err:>10.2f} {ratio:>14.4f}")
        print()

    print("\n--- PART 2: Short-interval error for various widths ---")
    print("For x=100000, interval [x-W, x+W]")
    print(f"{'W':>10} {'W/sqrt(x)':>10} {'K zeros':>8} {'pred_pi':>10} {'exact_pi':>10} {'error':>8} {'theory_K':>10}")
    print("-" * 80)

    x = 100000
    exact_pi_x = int(primepi(x))
    sqrt_x = np.sqrt(x)

    for W in [10, 50, 100, 500, 1000, 5000, 10000, 50000]:
        a = max(2, x - W)
        b = x + W
        exact_count = count_primes_sieve(a, b)

        # Theoretical number of zeros needed: ~ x / (W * log(x))
        theory_K = x / (W * np.log(x))

        for K in [1, 3, 5, 10, 15, 20, 25, 30]:
            if K > len(zeros):
                break
            pred = pi_short_interval_explicit(a, b, zeros[:K])
            err = pred - exact_count
            if abs(err) < 1.5:  # Found sufficient K
                print(f"{W:>10} {W/sqrt_x:>10.2f} {K:>8} {pred:>10.1f} {exact_count:>10} {err:>8.2f} {theory_K:>10.1f} ← SUFFICIENT")
                break
        else:
            # Report best result with max zeros
            pred = pi_short_interval_explicit(a, b, zeros[:len(zeros)])
            err = pred - exact_count
            print(f"{W:>10} {W/sqrt_x:>10.2f} {len(zeros):>8} {pred:>10.1f} {exact_count:>10} {err:>8.2f} {theory_K:>10.1f} (max zeros)")

    print("\n--- PART 3: Sinc decay verification ---")
    print("For each zero, measure the contribution to interval [x-W, x+W] vs [0, x]")
    print(f"{'gamma':>10} {'full_contrib':>14} {'W=1000':>12} {'W=100':>12} {'W=10':>12} {'theory_cutoff':>14}")
    print("-" * 80)

    x = 100000
    for gamma in zeros[:15]:
        full = 2 * li_rho(x, gamma)
        # For short interval, contribution scales as sinc(gamma*W/x)
        contribs = []
        for W in [1000, 100, 10]:
            # Short interval contribution
            c_upper = 2 * li_rho(x + W, gamma)
            c_lower = 2 * li_rho(x - W, gamma)
            contribs.append(c_upper - c_lower)

        cutoff = np.pi * x / 1000  # W=1000 cutoff
        print(f"{gamma:>10.3f} {full:>14.6f} {contribs[0]:>12.6f} {contribs[1]:>12.6f} {contribs[2]:>12.6f} {cutoff:>14.1f}")

    print("\n--- PART 4: Iterative refinement test ---")
    print("Round 1: Use K1 zeros to narrow interval. Round 2: Use K2 zeros for smaller interval.")
    print("Does round 2 need FEWER zeros than round 1?")

    x = 100000
    exact_pi = int(primepi(x))

    for K1 in [5, 10, 15, 20]:
        if K1 > len(zeros):
            break
        # Round 1: compute pi(x) with K1 zeros
        approx1 = pi_explicit_formula(x, zeros[:K1])
        err1 = approx1 - exact_pi
        W1 = abs(err1) + 5  # interval width after round 1

        # Round 2: short interval [x - W1, x + W1]
        a2, b2 = max(2, x - int(W1)), x + int(W1)
        exact_short = count_primes_sieve(a2, b2)

        # How many zeros needed for round 2?
        for K2 in range(1, len(zeros) + 1):
            pred2 = pi_short_interval_explicit(a2, b2, zeros[:K2])
            if abs(pred2 - exact_short) < 1.5:
                break

        theory_K2 = x / (W1 * np.log(x)) if W1 > 0 else float('inf')
        print(f"  K1={K1:>3}: err1={err1:>8.2f}, W1={W1:>8.1f}, K2_needed={K2:>3}, theory_K2={theory_K2:>8.1f}, K2/K1={K2/K1:>6.2f}")

    print("\n--- PART 5: Cost comparison ---")
    print("Comparing: (A) full pi(x) with K zeros vs (B) hybrid: K1 zeros + sieve W")
    print(f"{'x':>10} {'K_full':>8} {'K1_hybrid':>10} {'W_sieve':>10} {'cost_A':>12} {'cost_B':>12} {'ratio':>8}")
    print("-" * 75)

    for x in [1000, 10000, 100000]:
        exact = int(primepi(x))
        sqrt_x = np.sqrt(x)

        # Method A: Find min K for error < 1
        K_full = len(zeros)
        for K in range(1, len(zeros) + 1):
            approx = pi_explicit_formula(x, zeros[:K])
            if abs(approx - exact) < 1.0:
                K_full = K
                break

        # Method B: Hybrid — use K1 zeros to get width W, then sieve
        best_cost_B = float('inf')
        best_K1 = 0
        best_W = 0
        for K1 in range(1, len(zeros) + 1):
            approx = pi_explicit_formula(x, zeros[:K1])
            W = abs(approx - exact) + 2
            # Sieve cost ~ W * log(log(x)) + sqrt(x)
            cost_sieve = W * max(1, np.log(np.log(x))) + sqrt_x
            cost_total = K1 * np.log(x) + cost_sieve  # K1 * polylog for computing zeros
            if cost_total < best_cost_B:
                best_cost_B = cost_total
                best_K1 = K1
                best_W = W

        cost_A = K_full * np.log(x)  # simplified cost
        print(f"{x:>10} {K_full:>8} {best_K1:>10} {best_W:>10.1f} {cost_A:>12.1f} {best_cost_B:>12.1f} {cost_A/best_cost_B:>8.2f}")

    print("\n" + "=" * 80)
    print("CONCLUSIONS")
    print("=" * 80)

if __name__ == "__main__":
    test_short_interval_zero_requirement()
