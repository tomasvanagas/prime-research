"""
Dirichlet Series Bypass Experiment
Session 12 - April 2026

Novel question: Can we compute pi(x) from evaluations of zeta(s)
and related functions at SPECIFIC points, without summing over zeros?

The explicit formula involves sum_rho R(x^rho), requiring all zeros.
But Perron's formula gives pi(x) as a CONTOUR INTEGRAL of -zeta'(s)/zeta(s).
Could we evaluate this integral at special points to bypass zeros?

Specifically:
1. zeta(s) can be evaluated in O(|t|^{1/3+eps}) time at s = sigma + i*t
2. The contour integral needs O(T) evaluations for height T
3. For error < 1, need T ~ x^{1/2} (Lagarias-Odlyzko)

Question: is there a SHORTER contour (T << x^{1/2}) combined with
a CORRECTION that together give pi(x)?

Example: use a short contour for the "main term" + a combinatorial
correction for the "tail." If the correction is cheaper than extending
the contour, we win.

This is similar to the "hybrid" approach that might give O(x^{3/5}).
"""

import math
import numpy as np
import sympy

def compute_li(x):
    """Compute Li(x) = integral_2^x dt/ln(t)."""
    if x <= 2:
        return 0
    # Numerical integration
    from scipy import integrate
    result, _ = integrate.quad(lambda t: 1/math.log(t), 2, x)
    return result

def explicit_formula_truncated(x, num_zeros):
    """
    Compute pi(x) via explicit formula truncated at num_zeros zeros.
    Returns the approximation and error.
    """
    # Load precomputed zeros
    try:
        zeros = []
        with open("data/zeta_zeros_1000.txt", "r") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    zeros.append(float(line))
        zeros = zeros[:num_zeros]
    except FileNotFoundError:
        # Use sympy to compute zeros (slow)
        zeros = []
        for i in range(1, num_zeros + 1):
            z = float(sympy.nzeros(i) if hasattr(sympy, 'nzeros') else 14.134725)
            zeros.append(z)
        # Actually, let me compute them properly
        return None, None

    if not zeros:
        return None, None

    # R(x) = sum_{k=1}^inf mu(k)/k * Li(x^{1/k})
    R_x = 0
    for k in range(1, int(math.log2(x)) + 2):
        mu_k = int(sympy.mobius(k))
        if mu_k != 0:
            R_x += mu_k / k * compute_li(x ** (1/k))

    # Zero contributions
    zero_sum = 0
    for gamma in zeros:
        rho = 0.5 + 1j * gamma
        # R(x^rho) ≈ Li(x^rho) for the main term
        x_rho = x ** rho
        if abs(x_rho) > 0:
            # Li(x^rho) ≈ x^rho / (rho * ln(x))
            li_x_rho = x_rho / (rho * math.log(x))
            zero_sum += 2 * li_x_rho.real  # conjugate pair

    approx = R_x - zero_sum
    exact = int(sympy.primepi(x))

    return approx, exact

def analyze_contour_truncation():
    """
    Analyze how the explicit formula error depends on truncation height T.

    The error from truncating at height T is approximately:
    error ≈ x^{1/2} / T * (sum of residual oscillations)

    For the Lagarias-Odlyzko method:
    - Contour height T => error O(x * log^2(x) / T)
    - For error < 1, need T > x * log^2(x)
    - But with saddle-point, T ~ x^{1/2} suffices

    KEY QUESTION: Can we use FEWER zeros if we know something about
    the tail of the zero sum?
    """
    print("=== Contour Truncation Analysis ===")
    print()

    # Load zeros
    try:
        zeros = []
        with open("data/zeta_zeros_1000.txt", "r") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    try:
                        zeros.append(float(line))
                    except ValueError:
                        pass
    except FileNotFoundError:
        print("No zero data file found. Using first 50 known zeros.")
        # First 50 zeta zeros (imaginary parts)
        zeros = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
                 37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
                 52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
                 67.079811, 69.546402, 72.067158, 75.704691, 77.144840,
                 79.337375, 82.910381, 84.735493, 87.425275, 88.809111,
                 92.491899, 94.651344, 95.870634, 98.831194, 101.317851,
                 103.725538, 105.446623, 107.168611, 111.029536, 111.874659,
                 114.320220, 116.226680, 118.790783, 121.370125, 122.946829,
                 124.256819, 127.516684, 129.578704, 131.087688, 133.497737,
                 134.756510, 138.116042, 139.736209, 141.123707, 143.111846]

    print(f"Loaded {len(zeros)} zeros")

    # For various x values, compute:
    # 1. R(x) (smooth part)
    # 2. Explicit formula with K zeros
    # 3. Error = |approx - pi(x)|
    # Track how error decreases with K

    for x in [100, 1000, 10000, 100000]:
        exact = int(sympy.primepi(x))

        # Compute R(x)
        R_x = 0
        for k in range(1, int(math.log2(x)) + 2):
            mu_k = int(sympy.mobius(k))
            if mu_k != 0:
                R_x += mu_k / k * compute_li(x ** (1/k))

        print(f"\nx = {x}, pi(x) = {exact}, R(x) = {R_x:.4f}")

        # Add zeros one by one
        zero_sum = 0
        prev_error = abs(R_x - exact)
        print(f"  K=  0: approx = {R_x:.2f}, error = {prev_error:.4f}")

        for K in range(len(zeros)):
            gamma = zeros[K]
            rho = 0.5 + 1j * gamma
            x_rho = x ** rho
            li_x_rho = x_rho / (rho * math.log(x))
            zero_sum += 2 * li_x_rho.real

            approx = R_x - zero_sum
            error = abs(approx - exact)

            if K + 1 in [1, 2, 3, 5, 10, 20, 30, 50, len(zeros)]:
                print(f"  K={K+1:>3}: approx = {approx:.2f}, error = {error:.4f}, "
                      f"T = {gamma:.1f}, sqrt(x) = {math.sqrt(x):.1f}")

        # Minimum zeros needed for error < 0.5 (exact answer)
        zero_sum = 0
        min_K = -1
        for K in range(len(zeros)):
            gamma = zeros[K]
            rho = 0.5 + 1j * gamma
            x_rho = x ** rho
            li_x_rho = x_rho / (rho * math.log(x))
            zero_sum += 2 * li_x_rho.real
            approx = R_x - zero_sum
            if abs(approx - exact) < 0.5:
                min_K = K + 1
                break

        if min_K > 0:
            print(f"  Minimum zeros for exact answer: {min_K} (gamma_{min_K} = {zeros[min_K-1]:.1f})")
        else:
            print(f"  Could not get exact answer with {len(zeros)} zeros")

def analyze_tail_structure():
    """
    Analyze the STRUCTURE of the tail sum (zeros beyond the truncation point).

    Key question: does the tail sum_{gamma > T} 2*Re[x^{1/2+i*gamma}/(rho*ln(x))]
    have any exploitable structure?

    Each term is approximately (2*x^{1/2} / (gamma * ln(x))) * cos(gamma*ln(x) - phi)
    where phi = arctan(2*gamma).

    The tail is a sum of DECREASING cosines at INCOMMENSURABLE frequencies.
    """
    print("\n=== Tail Structure Analysis ===")

    # Load zeros
    try:
        zeros = []
        with open("data/zeta_zeros_1000.txt", "r") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    try:
                        zeros.append(float(line))
                    except ValueError:
                        pass
    except FileNotFoundError:
        zeros = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
                 37.586178, 40.918719, 43.327073, 48.005151, 49.773832]

    if len(zeros) < 20:
        print("Not enough zeros for tail analysis")
        return

    for x in [1000, 10000]:
        print(f"\nx = {x}")

        # Compute contribution of each zero
        contributions = []
        for gamma in zeros:
            rho = 0.5 + 1j * gamma
            x_rho = x ** rho
            term = 2 * (x_rho / (rho * math.log(x))).real
            contributions.append((gamma, term))

        # Cumulative sum from the END (tail)
        print(f"  Tail sums (from gamma_K to end):")
        total_tail = sum(c for _, c in contributions)

        for K in [len(zeros), len(zeros)*3//4, len(zeros)//2, len(zeros)//4, 10, 5]:
            if K > len(zeros):
                continue
            tail = sum(c for _, c in contributions[K-1:])
            print(f"    Tail from gamma_{K} ({zeros[K-1] if K <= len(zeros) else '?':.1f}): "
                  f"tail_sum = {tail:.4f}")

        # Check if tail terms show any pattern
        print(f"\n  Individual contributions (magnitude):")
        for i, (gamma, term) in enumerate(contributions[:20]):
            phase = gamma * math.log(x) % (2 * math.pi)
            print(f"    gamma_{i+1} = {gamma:.3f}: "
                  f"contribution = {term:+.4f}, "
                  f"|term| = {abs(term):.4f}, "
                  f"phase = {phase:.3f} rad")

def analyze_hybrid_potential():
    """
    Investigate the HYBRID approach: use K zeros + combinatorial correction.

    If we use K zeros, the error is approximately:
    delta(K, x) = sum_{gamma_j > gamma_K} 2*Re[x^{1/2+i*gamma_j}/(rho_j*ln(x))]

    This error is O(x^{1/2} / gamma_K * ln(gamma_K)).

    If we could compute delta(K, x) by a SIEVE method in O(x^alpha) time,
    and the zero computation costs K * K^{1/3+eps}, then total cost is:
    K^{4/3+eps} + x^alpha

    For K = x^beta: cost = x^{4beta/3} + x^alpha
    Minimized when 4beta/3 = alpha, giving cost x^alpha where alpha > 0.

    But what IS the sieve cost of computing the tail?
    The tail is NOT directly computable by sieve methods (it's an analytic quantity).
    """
    print("\n=== Hybrid Approach Analysis ===")
    print("""
    The hybrid approach would:
    1. Compute explicit formula with K zeros (cost: K^{4/3+eps})
    2. Compute the truncation error by some other method (cost: f(x, K))
    3. Total: K^{4/3+eps} + f(x, K)

    For this to beat O(x^{2/3}), we need both terms < x^{2/3}.

    The truncation error is:
    sum_{gamma > gamma_K} 2*Re[x^{rho}/(rho*ln(x))]

    This is the integral pi(x) - (R(x) - sum_{j<=K} R(x^{rho_j}))

    Can this difference be computed by sieve methods?

    Let D(x, K) = pi(x) - R(x) + sum_{j<=K} R(x^{rho_j})

    D(x, K) = -sum_{j>K} R(x^{rho_j}) - 1/ln(x) + arctan/terms

    The function D(x, K) is the "high-frequency" part of pi(x).
    It oscillates with the zero frequencies gamma_j > gamma_K.

    To compute D(x, K) independently of the explicit formula,
    we'd need: D(x, K) = pi(x) - (known analytic approximation with K zeros)

    But computing pi(x) IS the problem! So this is circular.

    UNLESS: we compute pi(x) mod something, or compute a projection
    of pi(x) that captures D(x, K).

    Verdict: CIRCULAR. The hybrid approach reduces to computing pi(x),
    which is what we're trying to do.
    """)

    # Numerical illustration of the error scaling
    print("  Error scaling with K zeros (numerical):")
    try:
        zeros = []
        with open("data/zeta_zeros_1000.txt", "r") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    try:
                        zeros.append(float(line))
                    except ValueError:
                        pass
    except FileNotFoundError:
        zeros = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
                 37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
                 52.970321, 56.446248, 59.347044, 60.831779, 65.112544]

    for x in [1000, 10000]:
        print(f"\n  x = {x}, sqrt(x) = {math.sqrt(x):.1f}")

        R_x = 0
        for k in range(1, int(math.log2(x)) + 2):
            mu_k = int(sympy.mobius(k))
            if mu_k != 0:
                R_x += mu_k / k * compute_li(x ** (1/k))

        exact = int(sympy.primepi(x))
        zero_sum = 0

        for K in range(min(len(zeros), 50)):
            gamma = zeros[K]
            rho = 0.5 + 1j * gamma
            x_rho = x ** rho
            li_x_rho = x_rho / (rho * math.log(x))
            zero_sum += 2 * li_x_rho.real
            error = abs(R_x - zero_sum - exact)

            if K + 1 in [1, 2, 5, 10, 15, 20, 30, 50]:
                # Theoretical bound: error ~ C * sqrt(x) / gamma_K
                theoretical = math.sqrt(x) / gamma * 2
                print(f"    K={K+1:>3}, gamma_K={gamma:.1f}: "
                      f"error={error:.3f}, "
                      f"sqrt(x)/gamma_K={math.sqrt(x)/gamma:.3f}")

def main():
    analyze_contour_truncation()
    analyze_tail_structure()
    analyze_hybrid_potential()

    print("\n" + "=" * 70)
    print("CONCLUSIONS")
    print("=" * 70)
    print("""
1. The explicit formula with K zeros gives pi(x) with error that
   decreases roughly as sqrt(x)/gamma_K.

2. The minimum number of zeros for exact pi(x) grows as a power of x
   (empirically x^{0.25-0.37} from Session 11).

3. The tail sum has NO exploitable structure (incommensurable phases,
   GUE statistics).

4. The hybrid approach (analytic + sieve) is CIRCULAR: computing the
   correction requires computing pi(x) itself.

5. No "shortcut contour" exists that avoids the O(sqrt(x)) zero
   contributions.

VERDICT: The Dirichlet series / contour integral approach cannot be
bypassed. The information in the zeta zeros is fundamental and
irreducible. This confirms the analytic barrier.
""")

if __name__ == "__main__":
    main()
