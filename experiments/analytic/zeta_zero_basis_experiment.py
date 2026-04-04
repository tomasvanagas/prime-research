#!/usr/bin/env python3
"""
Session 13 Experiment: Zeta zero basis functions for pi(x)

Instead of li(x^{1/k}), use the EXACT explicit formula terms:
  R(x^rho) for each zeta zero rho = 1/2 + i*gamma

The question: if we use K zeta zeros as basis functions, how does the
approximation error scale with K and x?

More importantly: is there a REORDERING or WEIGHTING of the zeros that
converges faster than the standard ordering by |gamma|?

If we could find an ordering where O(polylog(x)) terms suffice for
error < 0.5, that would be a breakthrough.
"""

import numpy as np
from scipy.special import expi
from sympy import primepi
import warnings
warnings.filterwarnings('ignore')

def li(x):
    """Logarithmic integral"""
    if x <= 1.01:
        return 0.0
    return float(expi(np.log(x)))

def R_term(x, gamma):
    """
    Compute the real part of R(x^rho) where rho = 1/2 + i*gamma.

    R(x^rho) ≈ li(x^rho) = Ei(rho * ln(x))

    We compute 2*Re[li(x^rho)] since zeros come in conjugate pairs.
    """
    if x <= 1:
        return 0.0
    log_x = np.log(x)
    # rho = 1/2 + i*gamma
    # rho * log(x) = log(x)/2 + i*gamma*log(x)
    z = complex(log_x / 2, gamma * log_x)

    # Ei(z) for complex z - use series expansion for small |z|
    # and asymptotic expansion for large |z|
    result = complex_Ei(z)
    return 2 * result.real  # Factor of 2 for conjugate pair

def complex_Ei(z, terms=200):
    """Compute Ei(z) for complex z using series expansion.
    Ei(z) = gamma_euler + ln(z) + sum_{k=1}^{inf} z^k / (k * k!)
    """
    gamma_euler = 0.5772156649015329

    # For large |z|, use asymptotic
    if abs(z) > 50:
        # Asymptotic: Ei(z) ~ e^z/z * (1 + 1/z + 2/z^2 + ...)
        result = np.exp(z) / z
        term = 1.0 / z
        for k in range(1, 20):
            term *= (k) / z
            result += np.exp(z) / z * term
            if abs(term) < 1e-15:
                break
        return result

    result = gamma_euler + np.log(z)
    term = z
    for k in range(1, terms):
        result += term / (k * np.math.factorial(k) if k < 170 else float('inf'))
        term *= z
        if k >= 170 or abs(term) / (k * np.math.factorial(min(k, 169))) < 1e-15:
            break
        term = z ** (k + 1)

    # Better: use recursion for the series
    result2 = gamma_euler + np.log(z)
    term = 1.0
    for k in range(1, terms):
        term *= z / k
        result2 += term / k
        if abs(term / k) < 1e-15:
            break

    return result2

def load_zeta_zeros(n=200):
    """Load precomputed zeta zeros."""
    import os
    zero_file = "/apps/aplikacijos/prime-research/data/zeta_zeros_200.txt"
    if os.path.exists(zero_file):
        zeros = []
        with open(zero_file) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    try:
                        zeros.append(float(line))
                    except ValueError:
                        continue
        return zeros[:n]

    # Fallback: first 50 zeros hardcoded
    return [
        14.134725141734693, 21.022039638771555, 25.010857580145688,
        30.424876125859513, 32.935061587739189, 37.586178158825671,
        40.918719012147495, 43.327073280914999, 48.005150881167159,
        49.773832477672302, 52.970321477714460, 56.446247697063394,
        59.347044002602353, 60.831778524609809, 65.112544048081606,
        67.079810529494174, 69.546401711173979, 72.067157674481907,
        75.704690699083933, 77.144840068874805, 79.337375020249367,
        82.910380854086030, 84.735492980517050, 87.425274613125229,
        88.809111207634465, 92.491899270558484, 94.651344040519838,
        95.870634228245309, 98.831194218193692, 101.31785100573139,
        103.72553804532511, 105.44662305232542, 107.16861118427640,
        111.02953554316967, 111.87465917699263, 114.32022091545271,
        116.22668032085755, 118.79078286597621, 121.37012500242066,
        122.94682929355258, 124.25681855434576, 127.51668387959649,
        129.57870420047855, 131.08768853093265, 133.49773720299758,
        134.75650975337387, 138.11604205453344, 139.73620895212138,
        141.12370740402112, 143.11184580762063
    ][:n]

def riemann_R(x, K=50):
    """Standard Riemann R(x)"""
    result = 0.0
    log_x = np.log(max(x, 2))
    term = 1.0
    for k in range(1, K+1):
        from sympy import mobius
        mu_k = int(mobius(k))
        if mu_k != 0:
            result += mu_k / k * li(x ** (1.0/k))
    return result

def experiment_zero_convergence():
    """How does the explicit formula converge as we add zeros?"""
    print("="*70)
    print("EXPERIMENT: Explicit formula convergence with zeta zeros")
    print("="*70)

    zeros = load_zeta_zeros(50)
    print(f"Loaded {len(zeros)} zeta zeros")

    # Test points
    test_x = [100, 500, 1000, 5000, 10000]

    for x in test_x:
        pi_exact = int(primepi(x))
        R_x = riemann_R(x)

        # Correction from 1/ln(x) and arctan term
        correction = -1.0/np.log(x) + (1/np.pi)*np.arctan(np.pi/np.log(x))

        # Cumulative zero contribution
        cumulative = R_x + correction

        print(f"\nx = {x}, pi(x) = {pi_exact}")
        print(f"  R(x) = {R_x:.4f} (error = {R_x - pi_exact:.4f})")

        errors = []
        for K in [1, 2, 5, 10, 20, 50]:
            K_use = min(K, len(zeros))
            zero_sum = sum(R_term(x, zeros[k]) for k in range(K_use))
            approx = R_x - zero_sum + correction
            err = approx - pi_exact
            errors.append((K, err))
            print(f"  K={K_use:3d} zeros: approx = {approx:.4f}, error = {err:.4f}")

        # Check if error decays monotonically
        abs_errors = [abs(e) for _, e in errors]
        monotone = all(abs_errors[i] >= abs_errors[i+1] for i in range(len(abs_errors)-1))
        print(f"  Monotone convergence: {monotone}")

def experiment_zero_reordering():
    """What if we use zeros in a DIFFERENT order?"""
    print("\n" + "="*70)
    print("EXPERIMENT: Zero reordering for faster convergence")
    print("="*70)

    zeros = load_zeta_zeros(50)

    x_values = [1000, 5000, 10000]

    for x in x_values:
        pi_exact = int(primepi(x))
        R_x = riemann_R(x)
        correction = -1.0/np.log(x) + (1/np.pi)*np.arctan(np.pi/np.log(x))

        # Compute all individual contributions
        contributions = []
        for i, gamma in enumerate(zeros):
            contrib = R_term(x, gamma)
            contributions.append((i, gamma, contrib))

        print(f"\nx = {x}, pi(x) = {pi_exact}, R(x) = {R_x:.4f}")

        # Strategy 1: Standard order (ascending gamma)
        cum_standard = R_x + correction
        standard_errors = []
        for i, gamma, contrib in contributions:
            cum_standard -= contrib
            standard_errors.append(cum_standard - pi_exact)

        # Strategy 2: Greedy order (largest |contribution| first)
        sorted_by_size = sorted(contributions, key=lambda t: abs(t[2]), reverse=True)
        cum_greedy = R_x + correction
        greedy_errors = []
        for i, gamma, contrib in sorted_by_size:
            cum_greedy -= contrib
            greedy_errors.append(cum_greedy - pi_exact)

        # Strategy 3: Greedy minimize error at each step
        remaining = list(contributions)
        cum_optimal = R_x + correction
        optimal_errors = []
        used = []
        for step in range(len(zeros)):
            best_err = float('inf')
            best_idx = 0
            for j, (i, gamma, contrib) in enumerate(remaining):
                trial = cum_optimal - contrib - pi_exact
                if abs(trial) < abs(best_err):
                    best_err = trial
                    best_idx = j
            i, gamma, contrib = remaining.pop(best_idx)
            cum_optimal -= contrib
            optimal_errors.append(cum_optimal - pi_exact)
            used.append(gamma)

        # Compare at K = 5, 10, 20, 50
        print(f"  {'K':>5} | {'Standard':>10} | {'|Contrib| order':>15} | {'Greedy optimal':>15}")
        for K in [5, 10, 20, min(50, len(zeros))]:
            k_idx = K - 1
            print(f"  {K:5d} | {standard_errors[k_idx]:10.4f} | {greedy_errors[k_idx]:15.4f} | {optimal_errors[k_idx]:15.4f}")

        # When does each strategy first achieve |error| < 0.5?
        for name, errors in [("Standard", standard_errors),
                              ("|Contrib|", greedy_errors),
                              ("Greedy", optimal_errors)]:
            for k, err in enumerate(errors):
                if abs(err) < 0.5:
                    print(f"  {name:12s} first |error| < 0.5 at K = {k+1}")
                    break
            else:
                print(f"  {name:12s} never achieves |error| < 0.5 with {len(zeros)} zeros")

        # Show the greedy-optimal ordering
        print(f"  Greedy-optimal first 10 zeros: gamma = {[f'{g:.2f}' for g in used[:10]]}")

def experiment_minimum_zeros_scaling():
    """How does the minimum number of zeros for exactness scale with x?"""
    print("\n" + "="*70)
    print("EXPERIMENT: Minimum zeros for exact pi(x) vs x")
    print("="*70)

    zeros = load_zeta_zeros(50)

    results = []
    for x in range(100, 10001, 100):
        pi_exact = int(primepi(x))
        R_x = riemann_R(x)
        correction = -1.0/np.log(x) + (1/np.pi)*np.arctan(np.pi/np.log(x))

        # Standard order
        cum = R_x + correction
        K_needed = None
        for k, gamma in enumerate(zeros):
            contrib = R_term(x, gamma)
            cum -= contrib
            if abs(cum - pi_exact) < 0.5:
                K_needed = k + 1
                break

        results.append((x, K_needed))

    print(f"\n{'x':>8} | {'K_min':>6} | {'sqrt(x)':>8} | {'K/sqrt(x)':>10} | {'ln(x)^2':>8} | {'K/ln(x)^2':>10}")
    print("-" * 70)
    for x, K in results:
        if K is not None:
            sq = np.sqrt(x)
            ln2 = np.log(x)**2
            print(f"{x:8d} | {K:6d} | {sq:8.1f} | {K/sq:10.4f} | {ln2:8.2f} | {K/ln2:10.4f}")
        else:
            print(f"{x:8d} | {'> '+str(len(zeros)):>6} | {np.sqrt(x):8.1f} | {'---':>10}")

    # Fit K_min = C * x^alpha
    valid = [(x, K) for x, K in results if K is not None]
    if len(valid) > 10:
        log_x = np.log([x for x, K in valid])
        log_K = np.log([K for x, K in valid])
        alpha, log_C = np.polyfit(log_x, log_K, 1)
        C = np.exp(log_C)
        print(f"\nFit: K_min ≈ {C:.4f} * x^{alpha:.4f}")
        print(f"(For polylog, need alpha → 0; for known barrier, alpha ≈ 0.25-0.5)")

if __name__ == "__main__":
    experiment_zero_convergence()
    experiment_zero_reordering()
    experiment_minimum_zeros_scaling()

    print("\n" + "="*70)
    print("KEY QUESTION: Can clever reordering/weighting of zeros reduce")
    print("the number needed from O(sqrt(x)) to O(polylog(x))?")
    print("If K_min ~ x^alpha with alpha > 0, the answer is NO.")
    print("="*70)
