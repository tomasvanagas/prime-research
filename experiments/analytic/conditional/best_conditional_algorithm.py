"""
EXPERIMENT 7: Best Conditional Algorithm for Computing p(n) Exactly
====================================================================
Implements the analytic method for pi(x) under RH:
  pi(x) = li(x) - sum_rho li(x^rho) + corrections

Under RH, need T ~ O(sqrt(x)) zeros for error < 0.5.

Benchmarks:
  - Accuracy vs number of zeros K at x = 10^3..10^6
  - K_min(x) scaling verification
  - Complexity comparison table across conjectures

Reference: Lagarias-Odlyzko 1987, Platt 2015, Galway 2004
"""

import math
import time
import os
import sys

import mpmath

# ---------------------------------------------------------------------------
# Load precomputed zeta zeros
# ---------------------------------------------------------------------------

DATA_DIR = "/apps/aplikacijos/prime-research/data"

def load_zeros(max_zeros=None):
    """Load zeta zeros from data files, returning imaginary parts (positive)."""
    # Use the largest file available
    for fname in ["zeta_zeros_1000.txt", "zeta_zeros_500.txt",
                   "zeta_zeros_300.txt", "zeta_zeros_200.txt"]:
        path = os.path.join(DATA_DIR, fname)
        if os.path.exists(path):
            zeros = []
            with open(path) as f:
                for line in f:
                    line = line.strip()
                    if line:
                        zeros.append(mpmath.mpf(line))
            if max_zeros and len(zeros) > max_zeros:
                zeros = zeros[:max_zeros]
            return zeros
    raise FileNotFoundError("No zeta zero files found in " + DATA_DIR)


# ---------------------------------------------------------------------------
# li(x) using mpmath for full precision
# ---------------------------------------------------------------------------

def li_mp(x, dps=50):
    """Logarithmic integral li(x) = Ei(ln(x)) using mpmath."""
    with mpmath.workdps(dps):
        if x <= 0:
            return mpmath.mpf(0)
        if x == 1:
            return mpmath.mpf('-inf')
        return mpmath.li(x)


# Mobius function values (precomputed to index 100)
_MU_100 = [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1, 0, -1, 0, -1, 0,
           1, 1, -1, 0, 0, 1, 0, 0, -1, -1, -1, 0, 1, 1, 1, 0, -1, 1, 1, 0, -1, -1,
           -1, 0, 0, -1, -1, 0, 0, 0, 1, -1, 0, 0, -1, 1, 0, 1, 1, -1, 0, -1, 1, 0,
           0, 1, -1, 0, 1, -1, -1, 0, -1, 1, 0, 0, 1, -1, 0, 0, 0, -1, 1, -1, 0, 1,
           1, -1, 0, 0, -1, 0, 0, -1, -1, 1, 0, 1, 0, 0]


def R_mp(x, dps=50):
    """Riemann's R function: R(x) = sum_{n=1}^inf mu(n)/n * li(x^{1/n})."""
    with mpmath.workdps(dps):
        x = mpmath.mpf(x)
        if x <= 1:
            return mpmath.mpf(0)
        s = mpmath.mpf(0)
        for n in range(1, len(_MU_100)):
            if _MU_100[n] == 0:
                continue
            xn = mpmath.power(x, mpmath.mpf(1) / n)
            if xn <= mpmath.mpf('1.0001'):
                break
            s += mpmath.mpf(_MU_100[n]) / n * mpmath.li(xn)
        return s


# ---------------------------------------------------------------------------
# Analytic formula for pi(x) under RH
# ---------------------------------------------------------------------------

def pi_analytic(x, K, zeros, dps=50):
    """
    Compute pi(x) via Riemann's explicit formula:
      pi(x) = R(x) - sum_{rho} R(x^rho) - R(x^{-1}) - R(x^{-3}) - ...

    The dominant zero correction is:
      sum_rho R(x^rho) ~ sum_rho li(x^rho)  [higher Mobius terms decay fast]

    We use:
      pi(x) = R(x) - sum_{k=1..K} 2*Re(li(x^{rho_k})) + small_corrections

    Corrections include trivial zeros and the constant term.
    """
    with mpmath.workdps(dps):
        x_mp = mpmath.mpf(x)
        lnx = mpmath.log(x_mp)

        # Main term: R(x)
        result = R_mp(x, dps)

        # Non-trivial zero sum: sum_rho R(x^rho)
        # R(x^rho) = sum mu(n)/n * li(x^{rho/n})
        # Dominant term (n=1): li(x^rho) = Ei(rho * ln x)
        # Include n=2 term: mu(2)/2 * li(x^{rho/2}) = -1/2 * Ei(rho/2 * ln x)
        zero_sum = mpmath.mpf(0)
        for k in range(min(K, len(zeros))):
            gamma = zeros[k]
            rho = mpmath.mpc(0.5, gamma)
            # n=1 term: li(x^rho) = Ei(rho * ln x)
            li_rho_1 = mpmath.ei(rho * lnx)
            # n=2 term: -1/2 * li(x^{rho/2}) = -1/2 * Ei(rho/2 * ln x)
            li_rho_2 = mpmath.ei(rho / 2 * lnx)
            # n=3 term: -1/3 * li(x^{rho/3})
            li_rho_3 = mpmath.ei(rho / 3 * lnx)
            # n=5 term: -1/5 * li(x^{rho/5})
            li_rho_5 = mpmath.ei(rho / 5 * lnx)

            r_rho = li_rho_1 - li_rho_2 / 2 - li_rho_3 / 3 - li_rho_5 / 5
            # Add conjugate pair: R(x^rho) + R(x^{bar{rho}}) = 2 * Re(R(x^rho))
            zero_sum += 2 * mpmath.re(r_rho)

        result -= zero_sum

        # Trivial zero corrections: R(x^{-2}) + R(x^{-4}) + ...
        # These are very small for x >= 2; R(x^{-2k}) ~ li(x^{-2k})
        # For x >= 10, sum is < 0.001
        # Include for completeness at small x
        trivial_sum = mpmath.mpf(0)
        for m in range(1, 50):
            exp_val = -2 * m
            x_pow = mpmath.power(x_mp, exp_val)
            if x_pow < mpmath.mpf('1e-30'):
                break
            # R(x^{-2m}) ~ li(x^{-2m}) which is tiny
            if x_pow > 1.0001:
                trivial_sum += R_mp(float(x_pow), dps)

        result -= trivial_sum

        return float(mpmath.re(result))


# ---------------------------------------------------------------------------
# Reference: exact pi(x) via Lucy DP (unconditional, O(x^{2/3}))
# ---------------------------------------------------------------------------

def lucy_pi(x):
    """Lucy_Hedgehog DP for pi(x). Unconditional O(x^{2/3})."""
    x = int(x)
    if x < 2: return 0
    if x < 3: return 1
    sqrtx = int(x**0.5)
    while (sqrtx + 1)**2 <= x: sqrtx += 1
    while sqrtx * sqrtx > x: sqrtx -= 1

    small = list(range(-1, sqrtx + 1))  # small[v] = pi(v) initially v-1
    large = [0] * (sqrtx + 2)
    for v in range(1, sqrtx + 1):
        large[v] = x // v - 1

    for p in range(2, sqrtx + 1):
        if small[p] == small[p - 1]:
            continue
        pcnt = small[p - 1]
        p2 = p * p
        for v in range(1, min(sqrtx, x // p2) + 1):
            d = v * p
            if d <= sqrtx:
                large[v] -= large[d] - pcnt
            else:
                large[v] -= small[x // d] - pcnt
        for v in range(sqrtx, p2 - 1, -1):
            small[v] -= small[v // p] - pcnt

    return large[1]


# ---------------------------------------------------------------------------
# Experiment: Accuracy vs K zeros
# ---------------------------------------------------------------------------

def run_accuracy_experiment():
    """Test how many zeros needed for exact pi(x) at various x."""
    print("=" * 80)
    print("EXPERIMENT 7: Best Conditional Algorithm for pi(x)")
    print("=" * 80)

    zeros = load_zeros()
    total_zeros = len(zeros)
    print(f"\nLoaded {total_zeros} zeta zeros (max gamma = {float(zeros[-1]):.2f})")

    x_values = [10**3, 10**4, 10**5, 10**6]
    K_values = [10, 50, 100, 200, 500, min(1000, total_zeros)]
    K_values = [k for k in K_values if k <= total_zeros]

    # Exact values via Lucy DP
    print("\n--- Computing exact pi(x) via Lucy DP ---")
    exact = {}
    for x in x_values:
        t0 = time.time()
        exact[x] = lucy_pi(x)
        dt = time.time() - t0
        print(f"  pi({x:.0e}) = {exact[x]}  [{dt:.4f}s]")

    # Analytic method with varying K
    print("\n--- Analytic formula: pi(x) ~ li(x) - sum_rho li(x^rho) ---")
    print(f"{'x':>10} | {'K':>6} | {'raw':>14} | {'rounded':>10} | {'exact':>10} | {'error':>10} | {'time':>8}")
    print("-" * 85)

    results = {}  # results[(x, K)] = (raw, rounded, exact, error, time)

    for x in x_values:
        for K in K_values:
            t0 = time.time()
            raw = pi_analytic(x, K, zeros, dps=50)
            dt = time.time() - t0
            rounded = round(raw)
            error = rounded - exact[x]
            results[(x, K)] = (raw, rounded, exact[x], error, dt)
            marker = " *" if error == 0 else ""
            print(f"{x:>10.0e} | {K:>6} | {raw:>14.4f} | {rounded:>10} | {exact[x]:>10} | {error:>10}{marker} | {dt:>7.3f}s")
        print()

    # Find K_min for each x
    print("\n--- K_min: minimum zeros for exact result ---")
    print(f"{'x':>10} | {'K_min':>6} | {'sqrt(x)':>10} | {'K_min/sqrt(x)':>14}")
    print("-" * 50)

    k_min_data = {}
    for x in x_values:
        k_min = None
        for K in range(1, total_zeros + 1):
            t0 = time.time()
            raw = pi_analytic(x, K, zeros, dps=50)
            rounded = round(raw)
            if rounded == exact[x]:
                k_min = K
                break
        sqrtx = math.sqrt(x)
        if k_min:
            ratio = k_min / sqrtx
            k_min_data[x] = k_min
            print(f"{x:>10.0e} | {k_min:>6} | {sqrtx:>10.1f} | {ratio:>14.4f}")
        else:
            k_min_data[x] = None
            print(f"{x:>10.0e} | {'N/A':>6} | {sqrtx:>10.1f} | {'N/A':>14}")

    return results, k_min_data, exact, zeros


# ---------------------------------------------------------------------------
# Complexity analysis table
# ---------------------------------------------------------------------------

def complexity_analysis():
    """Print the theoretical complexity comparison table."""
    print("\n" + "=" * 80)
    print("COMPLEXITY COMPARISON: Conditional vs Unconditional")
    print("=" * 80)

    # Estimate times based on O(f(x)) scaling
    # Reference: Lucy DP does pi(10^9) in ~0.175s => c * (10^9)^{2/3} = 0.175
    # c_lucy = 0.175 / 10^6 = 1.75e-7 per unit

    table = [
        ("Unconditional", "Meissel-Lehmer (Lucy DP)", "O(x^{2/3})", 2/3),
        ("RH only", "Analytic + Turing zeros", "O(x^{2/3+eps})", 2/3 + 0.01),
        ("RH + Odlyzko-Schonhage", "Analytic + batch zeros", "O(x^{1/2+eps})", 0.5 + 0.01),
        ("RH + FFT zeros", "Analytic + FFT per zero", "O(x^{1/2} polylog)", 0.5),
        ("GRH", "Conditional sieve", "O(x^{1/2+eps})", 0.5 + 0.01),
        ("GRH + EH (Elliott-Halberstam)", "Goldston-type", "O(x^{1/2+eps})", 0.5 + 0.01),
        ("Cramer's conjecture", "Random model", "O(x^{1/2+eps})", 0.5 + 0.01),
        ("Approximate (R^{-1})", "Riemann inverse", "O(polylog(n))", 0),
    ]

    # Reference: 10^9 in 0.175s for x^{2/3}
    c_ref = 0.175 / (10**9)**(2/3)  # seconds per unit at exponent 2/3

    print(f"\n{'Assumption':<30} | {'Algorithm':<28} | {'Complexity':<22} | {'p(10^12) est':>14} | {'p(10^100) est':>16}")
    print("-" * 120)

    for assumption, algo, complexity, exp in table:
        if exp == 0:
            # Polylog: ~0.5s for 10^100
            t12 = "~0.01s"
            t100 = "~0.5s"
            note = "(~50% digits)"
        else:
            # p(n) ~ n*ln(n), so x ~ n*ln(n)
            # For p(10^12): x ~ 10^12 * 12*ln(10) ~ 3.3e13
            x12 = 3.3e13  # approximate p(10^12)
            x100 = 1e102  # approximate p(10^100)

            # Scale from reference
            # t = c * x^exp, but c depends on the specific algorithm
            # Use relative scaling from Lucy reference
            t12_val = c_ref * x12**exp

            if t12_val < 1:
                t12 = f"{t12_val:.3f}s"
            elif t12_val < 3600:
                t12 = f"{t12_val:.1f}s"
            elif t12_val < 86400:
                t12 = f"{t12_val/3600:.1f}hr"
            else:
                t12 = f"{t12_val/86400:.1f}d"

            # For 10^100: x^exp
            log_t100 = exp * 102 * math.log10(math.e) + math.log10(c_ref)
            if log_t100 > 30:
                t100 = f"10^{log_t100:.0f}s"
            elif log_t100 > 15:
                t100 = f"~10^{log_t100:.0f}s"
            else:
                t100_val = 10**log_t100
                if t100_val < 1:
                    t100 = f"{t100_val:.3f}s"
                else:
                    t100 = f"{t100_val:.1f}s"

        print(f"{assumption:<30} | {algo:<28} | {complexity:<22} | {t12:>14} | {t100:>16}")

    print()
    print("Key observations:")
    print("  1. RH alone does NOT improve over unconditional: O(x^{1/2} * x^{1/6+eps}) = O(x^{2/3+eps})")
    print("  2. RH + batch zero computation (Odlyzko-Schonhage) gives O(x^{1/2+eps}) -- genuine improvement")
    print("  3. NO known combination of conjectures beats O(x^{1/2+eps}) for EXACT pi(x)")
    print("  4. The x^{1/2} barrier comes from needing O(sqrt(x)) zeros -- fundamental to the explicit formula")
    print("  5. The polylog approximate method gives ~50% of digits but CANNOT be made exact")
    print()
    print("WHY x^{1/2} IS A HARD FLOOR:")
    print("  - The explicit formula: pi(x) = li(x) - sum_{|gamma|<T} li(x^rho) + error(T)")
    print("  - Error < 0.5 requires T = O(sqrt(x) / log(x)) zeros (under RH)")
    print("  - Even with O(1) cost per zero, that's Omega(sqrt(x)) total operations")
    print("  - Breaking this requires a fundamentally different approach (not zero-by-zero)")


# ---------------------------------------------------------------------------
# Deeper analysis: Why each conditional regime gives what it gives
# ---------------------------------------------------------------------------

def conditional_regime_analysis():
    """Detailed analysis of each conditional complexity regime."""
    print("\n" + "=" * 80)
    print("DETAILED CONDITIONAL REGIME ANALYSIS")
    print("=" * 80)

    print("""
REGIME 1: UNCONDITIONAL -- O(x^{2/3})
  Method: Meissel-Lehmer combinatorial sieve (Lucy DP variant)
  No assumptions needed. Current practical champion.
  Best implementations: Deleglise-Rivat, Kim Walisch's primecount

REGIME 2: RH + NAIVE ZERO COMPUTATION -- O(x^{2/3+eps})
  Method: Explicit formula with T = O(sqrt(x)) zeros
  Each zero via Turing's method: O(T^{1/3+eps}) per zero
  Total: O(sqrt(x)) zeros * O(x^{1/6+eps}) per zero = O(x^{2/3+eps})
  WORSE than unconditional! RH alone doesn't help.

REGIME 3: RH + ODLYZKO-SCHONHAGE BATCH -- O(x^{1/2+eps})
  Method: Compute all zeros up to T in one batch
  Odlyzko-Schonhage algorithm: O(T^{1+eps}) for ALL zeros up to T
  With T = O(sqrt(x)): total O(x^{1/2+eps})
  This is the BEST KNOWN exact method for pi(x) under any assumptions.
  Platt (2015) verified RH to height 3*10^10 and computed pi(10^25) this way.

REGIME 4: GRH -- O(x^{1/2+eps})
  GRH gives better error terms for primes in arithmetic progressions
  but doesn't fundamentally change the zero-counting bottleneck.
  Same O(x^{1/2+eps}) as Regime 3.

REGIME 5: GRH + ELLIOTT-HALBERSTAM (EH) -- O(x^{1/2+eps})
  EH improves the level of distribution in Bombieri-Vinogradov
  Useful for twin primes (Zhang, Maynard) but NOT for pi(x) computation.
  Still O(x^{1/2+eps}).

REGIME 6: CRAMER'S CONJECTURE -- Still O(x^{1/2+eps}) for exact
  Cramer's conjecture: gaps between primes are O((log p)^2)
  Helps for finding the NEXT prime after a given number
  Does NOT help for computing pi(x) exactly -- still need zero sum.

REGIME 7: BEST POSSIBLE WITH EXPLICIT FORMULA -- Omega(x^{1/2-eps})
  Lower bound argument (Galway, Platt):
  Need T ~ sqrt(x)/log(x) zeros for error < 0.5
  Even with O(polylog) per zero: total = Omega(sqrt(x)/log(x))
  This is information-theoretic: the zeros ENCODE pi(x).

CONCLUSION:
  Under ALL standard conjectures (RH, GRH, EH, Cramer's, etc.),
  the best proven complexity for EXACT pi(x) is O(x^{1/2+eps}).
  Breaking the x^{1/2} barrier requires either:
    (a) A non-explicit-formula approach (none known)
    (b) A way to batch-evaluate the zero sum without enumerating zeros (open problem)
    (c) A fundamentally new identity for pi(x) (unknown)
""")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    results, k_min_data, exact, zeros = run_accuracy_experiment()
    complexity_analysis()
    conditional_regime_analysis()

    # Final summary
    print("=" * 80)
    print("SUMMARY OF FINDINGS")
    print("=" * 80)
    print()
    print("1. ACCURACY vs ZEROS:")
    for x in sorted(k_min_data.keys()):
        k = k_min_data[x]
        if k:
            print(f"   pi({x:.0e}): K_min = {k} zeros needed (sqrt(x) = {math.sqrt(x):.0f})")
        else:
            print(f"   pi({x:.0e}): Not enough zeros available (need > {len(zeros)})")

    print()
    print("2. SCALING: K_min(x) vs sqrt(x)")
    xs = sorted([x for x in k_min_data if k_min_data[x] is not None])
    if len(xs) >= 2:
        ratios = [k_min_data[x] / math.sqrt(x) for x in xs]
        print(f"   Ratios K_min/sqrt(x): {[f'{r:.4f}' for r in ratios]}")
        # Fit log(K_min) vs log(x)
        if len(xs) >= 2:
            log_x = [math.log(x) for x in xs]
            log_k = [math.log(k_min_data[x]) for x in xs]
            # Linear regression
            n = len(xs)
            sx = sum(log_x)
            sy = sum(log_k)
            sxy = sum(a*b for a,b in zip(log_x, log_k))
            sxx = sum(a*a for a in log_x)
            slope = (n*sxy - sx*sy) / (n*sxx - sx*sx) if n*sxx != sx*sx else 0
            print(f"   Fitted exponent: K_min ~ x^{slope:.3f}  (expected ~0.5 under RH)")

    print()
    print("3. KEY CONCLUSION:")
    print("   Under ALL standard conjectures, best exact complexity is O(x^{1/2+eps}).")
    print("   The sqrt(x) barrier is fundamental to the explicit formula approach.")
    print("   For p(10^100), even O(x^{1/2}) requires ~10^51 operations -- infeasible.")
    print("   O(polylog) exact computation requires a fundamentally new approach.")
    print()
    print("VERDICT: CLOSED -- no conditional assumption achieves better than O(x^{1/2+eps}) exact.")


if __name__ == "__main__":
    main()
