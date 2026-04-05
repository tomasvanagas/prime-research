"""
Approximate polynomial degree of the prime indicator chi_P and counting function pi(x).

REAL-VALUED approximate degree (not GF(p)):
  adeg_eps(f) = min degree d such that exists real polynomial p of degree d
                with |p(x) - f(x)| <= eps for all x in {0,1}^N.

This is the standard notion from:
  - Nisan-Szegedy (1994): degree vs approximate degree
  - Beals et al. (2001): relation to quantum query complexity
  - Bun-Thaler (2017): approximate degree and communication complexity

Key questions:
  1. Does adeg(chi_P) scale as Theta(N), sqrt(N), or polylog(N)?
  2. How does adeg(pi(x)) compare to adeg(chi_P)?
  3. What does this imply for quantum query complexity of primality?

Distinct from approx_degree.py (Session 23) which studied GF(p) agreement fractions.
This experiment studies the real-valued L_infinity approximation -- the standard notion.
"""

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csc_matrix
from sympy import isprime, primepi
from itertools import combinations
from math import comb
import time
import sys

FLUSH = lambda: sys.stdout.flush()
P = lambda *a, **kw: (print(*a, **kw), FLUSH())


# ============================================================
# CORE: Build monomial evaluation matrix
# ============================================================

def monomial_matrix(N, max_degree):
    """
    Build matrix M where M[x, j] = prod_{i in S_j} x_i
    for all x in {0,1}^N and all subsets S_j with |S_j| <= max_degree.

    Returns M (2^N x num_monomials) and list of subset bitmasks.
    """
    n_points = 2 ** N
    # Collect bitmasks for subsets of size <= max_degree
    masks = []
    for d in range(max_degree + 1):
        for S in combinations(range(N), d):
            mask = 0
            for bit in S:
                mask |= (1 << bit)
            masks.append(mask)

    n_mono = len(masks)
    M = np.zeros((n_points, n_mono), dtype=np.float64)

    for j, mask in enumerate(masks):
        for x in range(n_points):
            if (x & mask) == mask:
                M[x, j] = 1.0

    return M, masks


def truth_table_indicator(N):
    """chi_P(x) = 1 if x is prime, 0 otherwise, for x in 0..2^N-1."""
    return np.array([1.0 if isprime(x) else 0.0 for x in range(2 ** N)])


def truth_table_pi_normalized(N):
    """pi(x) / 2^N for x in 0..2^N-1, normalized to [0,1]."""
    vals = [float(primepi(x)) for x in range(2 ** N)]
    arr = np.array(vals)
    return arr / (2 ** N)


def truth_table_pi_raw(N):
    """pi(x) for x in 0..2^N-1."""
    return np.array([float(primepi(x)) for x in range(2 ** N)])


# ============================================================
# LP solver for approximate degree
# ============================================================

def min_approx_error_lp(M, f, timeout=120):
    """
    Given monomial matrix M (n_points x n_monomials) and target f (n_points,),
    find coefficients c and minimum epsilon such that:
        |M @ c - f| <= epsilon  componentwise

    LP: minimize eps
        subject to:  M @ c - f <=  eps
                    -M @ c + f <=  eps
    Variables: [c_0, ..., c_{m-1}, eps]
    """
    n_points, n_mono = M.shape
    n_vars = n_mono + 1

    # Objective: minimize eps (last variable)
    c_obj = np.zeros(n_vars)
    c_obj[-1] = 1.0

    # Build constraint matrix [M, -1; -M, -1] x <= [f; -f]
    ones_col = -np.ones((n_points, 1))
    A_upper = np.hstack([M, ones_col])
    A_lower = np.hstack([-M, ones_col])
    A_ub = np.vstack([A_upper, A_lower])
    b_ub = np.concatenate([f, -f])

    bounds = [(None, None)] * n_mono + [(0, None)]

    result = linprog(c_obj, A_ub=A_ub, b_ub=b_ub, bounds=bounds,
                     method='highs', options={'presolve': True, 'time_limit': timeout})

    if result.success:
        return result.x[-1], result.x[:n_mono]
    else:
        return None, None


# ============================================================
# Monomial count helper
# ============================================================

def n_monomials(N, d):
    return sum(comb(N, k) for k in range(d + 1))


# ============================================================
# EXPERIMENT 1: Approximate degree of chi_P (prime indicator)
# ============================================================

def experiment1_indicator():
    P("=" * 72)
    P("EXPERIMENT 1: APPROXIMATE DEGREE OF chi_P (PRIME INDICATOR)")
    P("=" * 72)
    P()
    P("For each N, find min degree d such that a real polynomial of degree d")
    P("epsilon-approximates chi_P on {0,1}^N. Threshold: eps < 0.49 (for rounding).")
    P()

    results = {}
    MAX_MONO = 4500  # limit monomial count for tractability

    for N in range(4, 13):
        n_pts = 2 ** N
        f = truth_table_indicator(N)
        n_primes = int(np.sum(f))
        density = n_primes / n_pts

        P(f"--- N = {N}  (2^N = {n_pts}, primes = {n_primes}, density = {density:.4f}) ---")

        errors_by_degree = {}
        t0 = time.time()

        for d in range(0, N + 1):
            nm = n_monomials(N, d)
            if nm > MAX_MONO:
                P(f"  deg {d:>3}: {nm} monomials -- skipping (limit {MAX_MONO})")
                break

            M, _ = monomial_matrix(N, d)
            eps, coeffs = min_approx_error_lp(M, f)

            if eps is not None:
                errors_by_degree[d] = eps
                sufficient = eps < 0.49
                marker = " *** SUFFICIENT (eps < 0.49)" if sufficient else ""
                P(f"  deg {d:>3}: eps = {eps:.6f}  ({nm} monomials){marker}")
                if sufficient:
                    break
            else:
                P(f"  deg {d:>3}: LP FAILED ({nm} monomials)")
                break

        elapsed = time.time() - t0

        adeg = None
        for d in sorted(errors_by_degree.keys()):
            if errors_by_degree[d] < 0.49:
                adeg = d
                break
        if adeg is None:
            adeg = max(errors_by_degree.keys()) if errors_by_degree else N
            # Mark as lower bound if we didn't reach sufficiency
            P(f"  => adeg_0.49(chi_P, N={N}) >= {adeg} (did not converge)  [{elapsed:.2f}s]")
        else:
            P(f"  => adeg_0.49(chi_P, N={N}) = {adeg}  [{elapsed:.2f}s]")

        results[N] = {
            'adeg': adeg,
            'errors': errors_by_degree,
            'n_primes': n_primes,
            'time': elapsed,
            'converged': any(e < 0.49 for e in errors_by_degree.values())
        }
        P()

    # Scaling summary
    P("\n" + "=" * 72)
    P("SCALING SUMMARY: adeg(chi_P) vs N")
    P("=" * 72)
    P(f"{'N':>4} {'adeg':>5} {'conv':>5} {'adeg/N':>8} {'adeg/sqrtN':>11} {'n_primes':>10}")
    P("-" * 55)
    for N in sorted(results.keys()):
        r = results[N]
        ratio_n = r['adeg'] / N
        ratio_sqrt = r['adeg'] / np.sqrt(N)
        conv = "yes" if r['converged'] else "no"
        P(f"{N:>4} {r['adeg']:>5} {conv:>5} {ratio_n:>8.3f} {ratio_sqrt:>11.3f} {r['n_primes']:>10}")

    return results


# ============================================================
# EXPERIMENT 2: Approximate degree of pi(x) (counting function)
# ============================================================

def experiment2_counting():
    P("\n" + "=" * 72)
    P("EXPERIMENT 2: APPROXIMATE DEGREE OF pi(x) (COUNTING FUNCTION)")
    P("=" * 72)
    P()
    P("f(x) = pi(x) / 2^N  (normalized to [0,1]).")
    P("Find min degree d with max |p(x) - pi(x)| < 0.5 (sufficient for rounding).")
    P()

    results = {}
    MAX_MONO = 4500

    for N in range(4, 13):
        n_pts = 2 ** N
        f_norm = truth_table_pi_normalized(N)
        f_raw = truth_table_pi_raw(N)
        max_pi = f_raw[-1]

        P(f"--- N = {N}  (2^N = {n_pts}, pi(2^N-1) = {int(max_pi)}) ---")

        errors_by_degree = {}
        t0 = time.time()

        for d in range(0, N + 1):
            nm = n_monomials(N, d)
            if nm > MAX_MONO:
                P(f"  deg {d:>3}: {nm} monomials -- skipping")
                break

            M, _ = monomial_matrix(N, d)
            eps, coeffs = min_approx_error_lp(M, f_norm)

            if eps is not None:
                raw_err = eps * (2 ** N)
                errors_by_degree[d] = {'eps_norm': eps, 'raw_err': raw_err}
                sufficient = raw_err < 0.5
                marker = " *** EXACT" if sufficient else ""
                P(f"  deg {d:>3}: norm_eps = {eps:.6f}, raw_err = {raw_err:.3f}{marker}")
                if sufficient:
                    break
            else:
                P(f"  deg {d:>3}: LP FAILED")
                break

        elapsed = time.time() - t0

        adeg_exact = None
        for d in sorted(errors_by_degree.keys()):
            if errors_by_degree[d]['raw_err'] < 0.5:
                adeg_exact = d
                break
        if adeg_exact is None:
            adeg_exact = N

        results[N] = {
            'adeg_exact': adeg_exact,
            'errors': errors_by_degree,
            'max_pi': int(max_pi),
            'time': elapsed
        }
        P(f"  => adeg_exact(pi, N={N}) = {adeg_exact}  [{elapsed:.2f}s]")
        P()

    # Scaling summary
    P("\n" + "=" * 72)
    P("SCALING SUMMARY: adeg(pi) vs N")
    P("=" * 72)
    P(f"{'N':>4} {'adeg_exact':>11} {'ratio_N':>8} {'ratio_sqrtN':>12} {'pi(2^N)':>8}")
    P("-" * 50)
    for N in sorted(results.keys()):
        r = results[N]
        ratio_n = r['adeg_exact'] / N
        ratio_sqrt = r['adeg_exact'] / np.sqrt(N)
        P(f"{N:>4} {r['adeg_exact']:>11} {ratio_n:>8.3f} {ratio_sqrt:>12.3f} {r['max_pi']:>8}")

    return results


# ============================================================
# EXPERIMENT 3: Quantum query complexity bounds
# ============================================================

def experiment3_quantum(indicator_results):
    P("\n" + "=" * 72)
    P("EXPERIMENT 3: QUANTUM QUERY COMPLEXITY BOUNDS")
    P("=" * 72)
    P()
    P("Beals et al. (2001): Q(f) >= adeg(f)/2")
    P("where Q(f) = bounded-error quantum query complexity of f.")
    P()
    P("Also: D(f) <= O(adeg(f)^6)  [Aaronson-Shi refinements]")
    P("      D(f) >= adeg(f)")
    P()

    P("Reference approximate degrees (known results):")
    P("  OR_N:       adeg = Theta(sqrt(N))")
    P("  AND_N:      adeg = Theta(sqrt(N))")
    P("  MAJORITY_N: adeg = Theta(sqrt(N))")
    P("  PARITY_N:   adeg = N (exact)")
    P("  ADDRESS_N:  adeg = Theta(sqrt(N))")
    P()

    P("Prime indicator results:")
    P(f"{'N':>4} {'adeg':>6} {'sqrtN':>6} {'N':>4} {'adeg/N':>7} {'adeg/sqrtN':>11} {'conv':>5}")
    P("-" * 50)
    for N in sorted(indicator_results.keys()):
        r = indicator_results[N]
        ad = r['adeg']
        sqrtN = np.sqrt(N)
        conv = "yes" if r['converged'] else "no"
        P(f"{N:>4} {ad:>6} {sqrtN:>6.2f} {N:>4} {ad/N:>7.3f} {ad/sqrtN:>11.3f} {conv:>5}")

    P()

    # Power-law fit on converged data only
    Ns = [n for n in sorted(indicator_results.keys()) if indicator_results[n]['converged']]
    if len(Ns) >= 3:
        log_N = np.log([float(n) for n in Ns])
        log_adeg = np.log([float(indicator_results[n]['adeg']) for n in Ns])

        if np.std(log_adeg) > 0.01:
            coeffs = np.polyfit(log_N, log_adeg, 1)
            alpha = coeffs[0]
            c = np.exp(coeffs[1])
            P(f"Power-law fit (converged points): adeg(chi_P) ~ {c:.3f} * N^{alpha:.3f}")
            P()

            if alpha > 0.85:
                P("CONCLUSION: adeg(chi_P) ~ Theta(N)  [linear, like PARITY]")
                P("  => Quantum query complexity Q(chi_P) = Omega(N)")
                P("  => No quantum speedup for primality via bit-queries")
                P("  => Classical deterministic query complexity D(chi_P) = Theta(N)")
            elif alpha > 0.4 and alpha < 0.6:
                P("CONCLUSION: adeg(chi_P) ~ Theta(sqrt(N))  [like MAJORITY/OR]")
                P("  => Quantum query complexity Q(chi_P) = Omega(sqrt(N))")
                P("  => Grover-like quadratic speedup possible")
            elif alpha < 0.3:
                P("CONCLUSION: adeg(chi_P) = o(sqrt(N))  [potentially polylog!]")
                P("  => Efficient quantum algorithm may exist")
                P("  => Classical shortcuts via lifting theorems possible")
            else:
                P(f"CONCLUSION: adeg(chi_P) ~ N^{alpha:.2f}  [intermediate]")


# ============================================================
# EXPERIMENT 4: Partial function (promise problem)
# ============================================================

def experiment4_partial():
    P("\n" + "=" * 72)
    P("EXPERIMENT 4: PARTIAL FUNCTION (PROMISE PROBLEM)")
    P("=" * 72)
    P()
    P("Approximate chi_P only on inputs coprime to small primes.")
    P("Exclude trivially composite numbers (even, div by 3, etc).")
    P()

    small_primes_list = [2, 3, 5, 7, 11, 13]
    MAX_MONO = 4500
    results = {}

    for N in range(4, 12):
        n_pts = 2 ** N
        f_full = truth_table_indicator(N)

        for k in [1, 2, 3]:
            primes_k = small_primes_list[:k]

            # Find inputs coprime to first k primes (keep the primes themselves)
            valid_indices = []
            for x in range(n_pts):
                if x <= 1:
                    valid_indices.append(x)
                    continue
                coprime = True
                for p in primes_k:
                    if x % p == 0 and x != p:
                        coprime = False
                        break
                if coprime:
                    valid_indices.append(x)

            n_valid = len(valid_indices)
            if n_valid < 5:
                continue

            f_partial = f_full[valid_indices]

            best_d = None
            for d in range(0, N + 1):
                nm = n_monomials(N, d)
                if nm > MAX_MONO:
                    break

                M_full, _ = monomial_matrix(N, d)
                M_partial = M_full[valid_indices, :]
                eps, _ = min_approx_error_lp(M_partial, f_partial, timeout=30)
                if eps is not None and eps < 0.49:
                    best_d = d
                    break

            if N not in results:
                results[N] = {}
            adeg_val = best_d if best_d is not None else N
            results[N][k] = {
                'n_valid': n_valid,
                'frac_valid': n_valid / n_pts,
                'adeg': adeg_val
            }

            primes_str = ",".join(str(p) for p in primes_k)
            P(f"  N={N}, coprime to {{{primes_str}}}: "
              f"{n_valid}/{n_pts} ({n_valid/n_pts:.1%}), adeg = {adeg_val}")

    # Summary table
    P()
    P("Summary: adeg on promise vs full")
    P(f"{'N':>4} {'cop(2)':>7} {'cop(2,3)':>9} {'cop(2,3,5)':>11}")
    P("-" * 35)
    for N in sorted(results.keys()):
        c1 = results[N].get(1, {}).get('adeg', '?')
        c2 = results[N].get(2, {}).get('adeg', '?')
        c3 = results[N].get(3, {}).get('adeg', '?')
        P(f"{N:>4} {c1!s:>7} {c2!s:>9} {c3!s:>11}")

    return results


# ============================================================
# EXPERIMENT 5: SOS / compositeness certificates
# ============================================================

def experiment5_sos():
    P("\n" + "=" * 72)
    P("EXPERIMENT 5: SOS CERTIFICATES FOR COMPOSITENESS")
    P("=" * 72)
    P()
    P("adeg(1 - chi_P) = adeg(chi_P) since both are {0,1}-valued.")
    P("More interesting: for each degree d, what fraction of composites")
    P("can be CERTIFIED composite by a degree-d polynomial certificate?")
    P()
    P("For composite n: find degree-d poly g with g(n)^2 > 0 and g(p) = 0 for all primes p.")
    P("This is a 'refutation' certificate. We check if the LP dual gives certificates.")
    P()

    MAX_MONO = 4500

    for N in range(4, 11):
        n_pts = 2 ** N
        f = truth_table_indicator(N)
        g = 1.0 - f

        n_composite = int(np.sum(g))
        P(f"--- N = {N}  (composites+{{0,1}}: {n_composite}/{n_pts}) ---")

        for d in range(1, N + 1):
            nm = n_monomials(N, d)
            if nm > MAX_MONO:
                P(f"  deg {d:>3}: skipped ({nm} monomials)")
                break

            M, _ = monomial_matrix(N, d)
            eps, _ = min_approx_error_lp(M, g)
            if eps is not None:
                sufficient = eps < 0.49
                marker = " *** SUFFICIENT" if sufficient else ""
                P(f"  deg {d:>3}: eps(1-chi_P) = {eps:.6f}{marker}")
                if sufficient:
                    break
            else:
                P(f"  deg {d:>3}: LP FAILED")
                break
        P()

    P("NOTE: adeg(1-chi_P) = adeg(chi_P) as expected (both {0,1}-valued,")
    P("complementary). The Lasserre hierarchy at level d certifies compositeness")
    P("when adeg <= d.")


# ============================================================
# EXPERIMENT 6: Detailed error decay for select N values
# ============================================================

def experiment6_error_curves():
    P("\n" + "=" * 72)
    P("EXPERIMENT 6: DETAILED ERROR DECAY CURVES")
    P("=" * 72)
    P()
    P("Full error vs degree for N=8,10: characterize the decay shape.")
    P()

    MAX_MONO = 4500

    for N in [8, 10]:
        n_pts = 2 ** N
        f = truth_table_indicator(N)
        P(f"--- N = {N} ---")
        P(f"{'deg':>5} {'eps':>12} {'log2(eps)':>10} {'monomials':>10}")
        P("-" * 42)

        for d in range(0, N + 1):
            nm = n_monomials(N, d)
            if nm > MAX_MONO:
                P(f"  {d:>3}    (skipped, {nm} monomials)")
                continue

            M, _ = monomial_matrix(N, d)
            eps, _ = min_approx_error_lp(M, f)
            if eps is not None and eps > 1e-15:
                P(f"  {d:>3} {eps:>12.8f} {np.log2(eps):>10.3f} {nm:>10}")
            elif eps is not None:
                P(f"  {d:>3} {eps:>12.2e} {'< -50':>10} {nm:>10}")
            else:
                P(f"  {d:>3}    LP FAILED")

        P()


# ============================================================
# MAIN
# ============================================================

def main():
    P("*" * 72)
    P("* APPROXIMATE POLYNOMIAL DEGREE OF PRIME INDICATOR & pi(x)")
    P("* Real-valued L_inf (Nisan-Szegedy / Beals et al. framework)")
    P("*" * 72)
    P()

    t_start = time.time()

    # Experiment 1: prime indicator
    ind_results = experiment1_indicator()

    # Experiment 2: counting function pi(x)
    pi_results = experiment2_counting()

    # Experiment 3: quantum query complexity
    experiment3_quantum(ind_results)

    # Experiment 4: promise / partial function
    experiment4_partial()

    # Experiment 5: SOS certificates
    experiment5_sos()

    # Experiment 6: error decay curves
    experiment6_error_curves()

    total_time = time.time() - t_start

    # ============================================================
    # FINAL SYNTHESIS
    # ============================================================
    P("\n" + "=" * 72)
    P("FINAL SYNTHESIS")
    P("=" * 72)

    P()
    P("1. APPROXIMATE DEGREE OF chi_P (INDICATOR):")
    Ns = sorted(ind_results.keys())
    adegs = [ind_results[n]['adeg'] for n in Ns]
    convs = [ind_results[n]['converged'] for n in Ns]
    P(f"   N:    {Ns}")
    P(f"   adeg: {adegs}")
    P(f"   conv: {convs}")

    conv_Ns = [n for n in Ns if ind_results[n]['converged']]
    if len(conv_Ns) >= 3:
        log_N = np.log([float(n) for n in conv_Ns])
        log_adeg = np.log([float(ind_results[n]['adeg']) for n in conv_Ns])
        if np.std(log_adeg) > 0.01:
            alpha = np.polyfit(log_N, log_adeg, 1)[0]
            P(f"   Scaling fit: adeg ~ N^{alpha:.3f}")
            if alpha > 0.85:
                P("   VERDICT: LINEAR (like PARITY). No shortcut via polynomial method.")
            elif alpha > 0.4:
                P(f"   VERDICT: Sublinear N^{alpha:.2f}. Interesting intermediate regime.")
            else:
                P("   VERDICT: Sub-sqrt. Potential breakthrough!")

    P()
    P("2. APPROXIMATE DEGREE OF pi(x) (COUNTING):")
    Ns_pi = sorted(pi_results.keys())
    adegs_pi = [pi_results[n]['adeg_exact'] for n in Ns_pi]
    P(f"   N:    {Ns_pi}")
    P(f"   adeg: {adegs_pi}")

    P()
    P("3. IMPLICATIONS FOR p(n) COMPUTATION:")
    P("   - adeg(chi_P) gives lower bound on quantum query complexity: Q >= adeg/2")
    P("   - If adeg = Theta(N): need to read all N bits, no shortcut exists")
    P("   - If adeg = o(N): some bits are redundant, shortcuts may exist")
    P("   - The approximate degree of pi(x) is the more relevant quantity")
    P("     for computing p(n), since p(n) = pi^{-1}(n).")

    P()
    P(f"Total runtime: {total_time:.1f}s")


if __name__ == "__main__":
    main()
