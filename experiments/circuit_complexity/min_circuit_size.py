#!/usr/bin/env python3
"""
Minimum Boolean Circuit Size for chi_P(x) = [x is prime]

For N-bit inputs (x in {0, ..., 2^N - 1}), investigate:
1. Truth table of primality indicator
2. Complexity measures (sensitivity, block sensitivity, certificate complexity)
3. BDD size upper bounds (trying multiple variable orderings)
4. Comparison with random functions of same density
5. Growth rate analysis: polynomial vs exponential
6. Decision tree depth (deterministic query complexity)
"""

import itertools
import math
import random
import sys
import time
from collections import defaultdict
from functools import lru_cache

import numpy as np
from scipy.optimize import curve_fit

# ---------------------------------------------------------------------------
# Primality
# ---------------------------------------------------------------------------

def sieve(limit):
    """Return set of primes up to limit."""
    if limit < 2:
        return set()
    is_prime = [False, False] + [True] * (limit - 1)
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, limit + 1, i):
                is_prime[j] = False
    return set(i for i, v in enumerate(is_prime) if v)


def truth_table(N):
    """Return list of length 2^N: tt[x] = 1 if x is prime, else 0."""
    primes = sieve(2**N - 1)
    return [1 if x in primes else 0 for x in range(2**N)]


# ---------------------------------------------------------------------------
# Complexity measures
# ---------------------------------------------------------------------------

def essential_variables(tt, N):
    """Count variables that actually affect the output."""
    n = 2**N
    count = 0
    for bit in range(N):
        affects = False
        for x in range(n):
            neighbour = x ^ (1 << bit)
            if tt[x] != tt[neighbour]:
                affects = True
                break
        if affects:
            count += 1
    return count


def sensitivity(tt, N):
    """Max over all x of number of single-bit flips that change output."""
    n = 2**N
    max_s = 0
    for x in range(n):
        s = 0
        for bit in range(N):
            if tt[x] != tt[x ^ (1 << bit)]:
                s += 1
        max_s = max(max_s, s)
    return max_s


def sensitivity_on_0(tt, N):
    """Sensitivity restricted to inputs where f=0."""
    n = 2**N
    max_s = 0
    for x in range(n):
        if tt[x] == 1:
            continue
        s = 0
        for bit in range(N):
            if tt[x ^ (1 << bit)] != tt[x]:
                s += 1
        max_s = max(max_s, s)
    return max_s


def sensitivity_on_1(tt, N):
    """Sensitivity restricted to inputs where f=1."""
    n = 2**N
    max_s = 0
    for x in range(n):
        if tt[x] == 0:
            continue
        s = 0
        for bit in range(N):
            if tt[x ^ (1 << bit)] != tt[x]:
                s += 1
        max_s = max(max_s, s)
    return max_s


def block_sensitivity(tt, N):
    """
    Block sensitivity: max over x of the max number of DISJOINT blocks B_1,...,B_k
    of variables such that flipping each block changes f(x).
    We use a greedy approach for each x.
    """
    n = 2**N
    max_bs = 0
    for x in range(n):
        # Find all subsets of bits whose flip changes f(x)
        # Greedy: try single bits first, then pairs, etc.
        used = set()
        count = 0
        # Try single bits
        singles = []
        for bit in range(N):
            if bit not in used and tt[x ^ (1 << bit)] != tt[x]:
                singles.append(bit)
        # Try to find disjoint sensitive blocks greedily
        # Start with larger blocks for potentially better coverage
        # But for efficiency, just use single-bit blocks (which gives sensitivity)
        # then try to find additional multi-bit blocks

        # Actually, for correctness: enumerate all possible sensitive blocks
        # and find maximum disjoint packing. For small N this is feasible.

        sensitive_blocks = []
        # Enumerate subsets up to size N (all non-empty subsets of {0,...,N-1})
        # For N <= 12 this is 2^12 = 4096 subsets max per x -- feasible
        for mask in range(1, 2**N):
            bits_in_block = []
            flip_mask = 0
            for bit in range(N):
                if mask & (1 << bit):
                    bits_in_block.append(bit)
                    flip_mask |= (1 << bit)
            if tt[x ^ flip_mask] != tt[x]:
                sensitive_blocks.append(mask)

        # Greedy max disjoint packing (by smallest block first for better packing)
        sensitive_blocks.sort(key=lambda m: bin(m).count('1'))
        used_mask = 0
        count = 0
        for block in sensitive_blocks:
            if block & used_mask == 0:  # disjoint
                used_mask |= block
                count += 1
        max_bs = max(max_bs, count)
    return max_bs


def certificate_complexity(tt, N):
    """
    Certificate complexity: for each x, the minimum number of bits that must be
    fixed to certify f(x). We find the min certificate for each x and return max.

    A certificate for x with f(x)=v is a subset S of variables such that
    for all y agreeing with x on S, f(y)=v.

    We find minimum certificate size by trying subsets of increasing size.
    For N <= 12 we use a smarter approach: for each x, find which bits are
    "necessary" by checking if removing a bit from the full certificate still works.
    """
    n = 2**N
    max_cert = 0

    for x in range(n):
        v = tt[x]
        # Start with all N bits as certificate, then try to remove bits
        cert_bits = list(range(N))

        # Greedy removal: try removing each bit
        changed = True
        while changed:
            changed = False
            for bit in list(cert_bits):
                # Try removing this bit
                trial = [b for b in cert_bits if b != bit]
                # Check if trial is still a valid certificate
                # We need: for all y agreeing with x on trial bits, f(y) = v
                trial_mask = sum(1 << b for b in trial)
                x_pattern = x & trial_mask
                valid = True
                # Check all y that agree with x on trial bits
                free_bits = [b for b in range(N) if b not in trial]
                for combo in range(2**len(free_bits)):
                    y = x_pattern
                    for i, fb in enumerate(free_bits):
                        if combo & (1 << i):
                            y |= (1 << fb)
                    if tt[y] != v:
                        valid = False
                        break
                if valid:
                    cert_bits = trial
                    changed = True
                    break  # restart the removal loop

        max_cert = max(max_cert, len(cert_bits))

    return max_cert


# ---------------------------------------------------------------------------
# BDD construction
# ---------------------------------------------------------------------------

class BDDNode:
    """A node in a Binary Decision Diagram."""
    __slots__ = ['var', 'low', 'high', 'id']
    _counter = 0

    def __init__(self, var, low, high):
        self.var = var
        self.low = low
        self.high = high
        BDDNode._counter += 1
        self.id = BDDNode._counter


def build_bdd(tt, N, var_order):
    """
    Build a reduced ordered BDD for the truth table with given variable ordering.
    Returns the number of non-terminal nodes.

    var_order: list of variable indices, from root to leaves.
    """
    # Terminal nodes
    ZERO = 0
    ONE = 1

    # Unique table: (var, low_id, high_id) -> node_id
    unique_table = {}
    node_count = [0]  # mutable counter

    # Cache for node creation
    def mk(var, low, high):
        """Create or retrieve a BDD node."""
        if low == high:
            return low
        key = (var, low, high)
        if key in unique_table:
            return unique_table[key]
        node_count[0] += 1
        node_id = node_count[0] + 1  # +1 to avoid clash with 0,1 terminals
        unique_table[key] = node_id
        return node_id

    # Build bottom-up using the variable ordering
    # For each assignment to variables in var_order[level:], we know the sub-function

    # Map from (partial assignment as tuple) -> terminal value
    # Start from full truth table, reduce level by level from bottom

    # Actually, build top-down recursively with memoization

    memo = {}

    def build(level, assignment_mask, assignment_val):
        """
        Build BDD for sub-function at given level.
        assignment_mask: which variables are already assigned
        assignment_val: their values (as a single int, bit positions = variable indices)
        """
        if level == N:
            # All variables assigned, look up truth table
            return ONE if tt[assignment_val] else ZERO

        key = (level, assignment_val & assignment_mask)
        if key in memo:
            return memo[key]

        var = var_order[level]
        bit = 1 << var
        new_mask = assignment_mask | bit

        low = build(level + 1, new_mask, assignment_val & ~bit)
        high = build(level + 1, new_mask, assignment_val | bit)

        result = mk(var, low, high)
        memo[key] = result
        return result

    build(0, 0, 0)
    return node_count[0]


def min_bdd_size(tt, N, max_orderings=None):
    """
    Find minimum BDD size over variable orderings.
    For N <= 8: try all N! orderings.
    For N > 8: try natural, reverse, and random orderings.
    Returns (min_size, best_order).
    """
    if max_orderings is None:
        if math.factorial(N) <= 5040:  # N <= 7
            max_orderings = math.factorial(N)
        else:
            max_orderings = 102  # natural + reverse + 100 random

    best_size = float('inf')
    best_order = None

    if math.factorial(N) <= max_orderings:
        # Try all orderings
        for order in itertools.permutations(range(N)):
            order = list(order)
            size = build_bdd(tt, N, order)
            if size < best_size:
                best_size = size
                best_order = order
    else:
        # Try natural order
        order = list(range(N))
        size = build_bdd(tt, N, order)
        if size < best_size:
            best_size = size
            best_order = list(order)

        # Try reverse order
        order = list(range(N-1, -1, -1))
        size = build_bdd(tt, N, order)
        if size < best_size:
            best_size = size
            best_order = list(order)

        # Try random orderings
        for _ in range(max_orderings - 2):
            order = list(range(N))
            random.shuffle(order)
            size = build_bdd(tt, N, order)
            if size < best_size:
                best_size = size
                best_order = list(order)

    return best_size, best_order


# ---------------------------------------------------------------------------
# Decision tree depth (deterministic query complexity)
# ---------------------------------------------------------------------------

def decision_tree_depth(tt, N):
    """
    Compute the minimum-depth decision tree for the Boolean function.
    Uses dynamic programming over sub-functions (identified by partial assignments).

    A sub-function is defined by which variables are fixed and to what values.
    We represent it as: (free_vars_frozenset, fixed_pattern) where fixed_pattern
    gives the values of non-free variables.

    For N <= 12 this is feasible because the number of distinct sub-functions
    is bounded.
    """
    # For efficiency, represent sub-functions by their truth tables over free variables
    # Use memoization on the sub-function truth table

    memo = {}

    def min_depth(fixed_mask, fixed_val):
        """
        Return minimum query depth needed to determine f, given that
        variables in fixed_mask are already set to values in fixed_val.
        """
        # Determine the set of outputs for all completions
        free_bits = [b for b in range(N) if not (fixed_mask & (1 << b))]
        n_free = len(free_bits)

        # Compute truth table of sub-function
        vals = set()
        sub_tt = []
        for combo in range(2**n_free):
            x = fixed_val
            for i, fb in enumerate(free_bits):
                if combo & (1 << i):
                    x |= (1 << fb)
            sub_tt.append(tt[x])
            vals.add(tt[x])

        if len(vals) == 1:
            return 0  # constant function, no queries needed

        # Memoize on sub-function truth table
        key = tuple(sub_tt)
        if key in memo:
            return memo[key]

        # Try querying each free variable, take the best
        best = n_free  # worst case: query all
        for bit in free_bits:
            new_mask = fixed_mask | (1 << bit)
            # Query bit, branch on 0 and 1
            d0 = min_depth(new_mask, fixed_val)  # bit = 0
            d1 = min_depth(new_mask, fixed_val | (1 << bit))  # bit = 1
            depth = 1 + max(d0, d1)
            best = min(best, depth)
            if best <= 1:
                break  # can't do better

        memo[key] = best
        return best

    return min_depth(0, 0)


# ---------------------------------------------------------------------------
# Random function comparison
# ---------------------------------------------------------------------------

def random_tt_same_density(N, density, rng):
    """Generate a random truth table with approximately the same density of 1s."""
    n = 2**N
    num_ones = round(density * n)
    tt = [0] * n
    ones_pos = rng.sample(range(n), num_ones)
    for p in ones_pos:
        tt[p] = 1
    return tt


# ---------------------------------------------------------------------------
# Growth rate fitting
# ---------------------------------------------------------------------------

def fit_polynomial(Ns, sizes):
    """Fit size = a * N^c. Return (a, c, R^2)."""
    Ns = np.array(Ns, dtype=float)
    sizes = np.array(sizes, dtype=float)

    # log(size) = log(a) + c * log(N)
    log_N = np.log(Ns)
    log_s = np.log(sizes)

    try:
        coeffs = np.polyfit(log_N, log_s, 1)
        c = coeffs[0]
        a = np.exp(coeffs[1])
        predicted = a * Ns**c
        ss_res = np.sum((sizes - predicted)**2)
        ss_tot = np.sum((sizes - np.mean(sizes))**2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        return a, c, r2
    except Exception:
        return None, None, 0


def fit_exponential(Ns, sizes):
    """Fit size = a * b^N. Return (a, b, R^2)."""
    Ns = np.array(Ns, dtype=float)
    sizes = np.array(sizes, dtype=float)

    # log(size) = log(a) + N * log(b)
    log_s = np.log(sizes)

    try:
        coeffs = np.polyfit(Ns, log_s, 1)
        b = np.exp(coeffs[0])
        a = np.exp(coeffs[1])
        predicted = a * b**Ns
        ss_res = np.sum((sizes - predicted)**2)
        ss_tot = np.sum((sizes - np.mean(sizes))**2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        return a, b, r2
    except Exception:
        return None, None, 0


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    random.seed(42)
    rng = random.Random(42)

    Ns = [4, 5, 6, 7, 8, 9, 10, 11, 12]

    results = []

    print("=" * 100)
    print("MINIMUM BOOLEAN CIRCUIT SIZE FOR chi_P(x) = [x is prime]")
    print("=" * 100)
    print()

    for N in Ns:
        t0 = time.time()
        n = 2**N
        tt = truth_table(N)
        num_primes = sum(tt)
        density = num_primes / n

        print(f"--- N = {N} (x in [0, {n-1}]), {num_primes} primes, density = {density:.4f} ---")

        # Complexity measures
        ess = essential_variables(tt, N)
        sens = sensitivity(tt, N)
        s0 = sensitivity_on_0(tt, N)
        s1 = sensitivity_on_1(tt, N)

        # Block sensitivity and certificate complexity -- skip for large N
        if N <= 10:
            bs = block_sensitivity(tt, N)
        else:
            bs = None

        if N <= 10:
            cert = certificate_complexity(tt, N)
        else:
            cert = None

        # BDD size
        if N <= 7:
            max_ord = None  # all orderings
        elif N == 8:
            max_ord = 5000
        else:
            max_ord = 102

        bdd_size, bdd_order = min_bdd_size(tt, N, max_orderings=max_ord)

        # Decision tree depth -- expensive for large N
        if N <= 10:
            dt_depth = decision_tree_depth(tt, N)
        else:
            dt_depth = None

        # Random comparison
        rand_bdd_sizes = []
        rand_sens = []
        rand_bs = []
        n_random = 20

        for _ in range(n_random):
            rtt = random_tt_same_density(N, density, rng)

            # BDD with just natural + reverse + 10 random orderings for speed
            rs, _ = min_bdd_size(rtt, N, max_orderings=12 if N > 7 else max_ord)
            rand_bdd_sizes.append(rs)

            rs_val = sensitivity(rtt, N)
            rand_sens.append(rs_val)

            if N <= 10:
                rbs = block_sensitivity(rtt, N)
                rand_bs.append(rbs)

        elapsed = time.time() - t0

        # Store results
        result = {
            'N': N, 'n': n, 'num_primes': num_primes, 'density': density,
            'essential': ess, 'sensitivity': sens, 's0': s0, 's1': s1,
            'block_sensitivity': bs, 'certificate': cert,
            'bdd_size': bdd_size, 'bdd_order': bdd_order,
            'dt_depth': dt_depth,
            'rand_bdd_mean': np.mean(rand_bdd_sizes),
            'rand_bdd_std': np.std(rand_bdd_sizes),
            'rand_sens_mean': np.mean(rand_sens),
            'rand_sens_std': np.std(rand_sens),
            'rand_bs_mean': np.mean(rand_bs) if rand_bs else None,
            'rand_bs_std': np.std(rand_bs) if rand_bs else None,
            'elapsed': elapsed,
        }
        results.append(result)

        print(f"  Essential variables:    {ess}/{N}")
        print(f"  Sensitivity:           {sens}  (s0={s0}, s1={s1})")
        if bs is not None:
            print(f"  Block sensitivity:     {bs}")
        if cert is not None:
            print(f"  Certificate complexity:{cert}")
        print(f"  BDD size (min found):  {bdd_size}  (order={bdd_order})")
        if dt_depth is not None:
            print(f"  Decision tree depth:   {dt_depth}")
        print(f"  Random BDD (same density): {np.mean(rand_bdd_sizes):.1f} +/- {np.std(rand_bdd_sizes):.1f}")
        print(f"  Random sensitivity:        {np.mean(rand_sens):.1f} +/- {np.std(rand_sens):.1f}")
        if rand_bs:
            print(f"  Random block sensitivity:  {np.mean(rand_bs):.1f} +/- {np.std(rand_bs):.1f}")
        print(f"  Time: {elapsed:.2f}s")
        print()

    # -----------------------------------------------------------------------
    # Summary table
    # -----------------------------------------------------------------------
    print()
    print("=" * 100)
    print("SUMMARY TABLE")
    print("=" * 100)
    hdr = f"{'N':>3} {'2^N':>6} {'#prm':>5} {'dens':>6} {'ess':>4} {'sens':>5} {'bs':>4} {'cert':>5} {'BDD':>6} {'DT':>4} {'rBDD':>10} {'rSens':>10}"
    print(hdr)
    print("-" * len(hdr))
    for r in results:
        bs_str = str(r['block_sensitivity']) if r['block_sensitivity'] is not None else '-'
        cert_str = str(r['certificate']) if r['certificate'] is not None else '-'
        dt_str = str(r['dt_depth']) if r['dt_depth'] is not None else '-'
        rbdd_str = f"{r['rand_bdd_mean']:.1f}+/-{r['rand_bdd_std']:.1f}"
        rsens_str = f"{r['rand_sens_mean']:.1f}+/-{r['rand_sens_std']:.1f}"
        print(f"{r['N']:>3} {r['n']:>6} {r['num_primes']:>5} {r['density']:>6.4f} {r['essential']:>4} {r['sensitivity']:>5} {bs_str:>4} {cert_str:>5} {r['bdd_size']:>6} {dt_str:>4} {rbdd_str:>10} {rsens_str:>10}")

    # -----------------------------------------------------------------------
    # Growth rate analysis
    # -----------------------------------------------------------------------
    print()
    print("=" * 100)
    print("GROWTH RATE ANALYSIS: BDD size of chi_P vs N")
    print("=" * 100)

    bdd_Ns = [r['N'] for r in results]
    bdd_sizes = [r['bdd_size'] for r in results]

    a_p, c_p, r2_p = fit_polynomial(bdd_Ns, bdd_sizes)
    a_e, b_e, r2_e = fit_exponential(bdd_Ns, bdd_sizes)

    print(f"\nPolynomial fit:   BDD_size = {a_p:.4f} * N^{c_p:.4f}")
    print(f"  R^2 = {r2_p:.6f}")
    print(f"\nExponential fit:  BDD_size = {a_e:.4f} * {b_e:.4f}^N")
    print(f"  R^2 = {r2_e:.6f}")

    if r2_e > r2_p:
        print(f"\n>> Exponential fit is BETTER (R^2 = {r2_e:.6f} vs {r2_p:.6f})")
        print(f">> BDD size appears to grow EXPONENTIALLY: ~{b_e:.3f}^N")
    else:
        print(f"\n>> Polynomial fit is BETTER (R^2 = {r2_p:.6f} vs {r2_e:.6f})")
        print(f">> BDD size appears to grow POLYNOMIALLY: ~N^{c_p:.3f}")

    # Also fit random BDD sizes
    rand_bdd_means = [r['rand_bdd_mean'] for r in results]
    a_rp, c_rp, r2_rp = fit_polynomial(bdd_Ns, rand_bdd_means)
    a_re, b_re, r2_re = fit_exponential(bdd_Ns, rand_bdd_means)

    print(f"\nRandom functions (same density):")
    print(f"  Polynomial fit:  R^2 = {r2_rp:.6f},  ~N^{c_rp:.3f}" if c_rp else "  Polynomial fit: failed")
    print(f"  Exponential fit: R^2 = {r2_re:.6f},  ~{b_re:.3f}^N" if b_re else "  Exponential fit: failed")

    # -----------------------------------------------------------------------
    # BDD size ratio: chi_P vs random
    # -----------------------------------------------------------------------
    print()
    print("=" * 100)
    print("BDD SIZE RATIO: chi_P / random(same density)")
    print("=" * 100)
    for r in results:
        ratio = r['bdd_size'] / r['rand_bdd_mean'] if r['rand_bdd_mean'] > 0 else float('inf')
        print(f"  N={r['N']:>2}:  chi_P={r['bdd_size']:>6}  random={r['rand_bdd_mean']:>8.1f}  ratio={ratio:.3f}")

    # -----------------------------------------------------------------------
    # Key findings
    # -----------------------------------------------------------------------
    print()
    print("=" * 100)
    print("KEY FINDINGS")
    print("=" * 100)

    # Check if all variables are essential
    all_essential = all(r['essential'] == r['N'] for r in results)
    print(f"\n1. All variables essential for all N? {all_essential}")

    # Sensitivity trend
    print(f"\n2. Sensitivity values: {[r['sensitivity'] for r in results]}")
    print(f"   (For N bits, max possible sensitivity = N)")

    # BDD size trend
    print(f"\n3. BDD sizes: {bdd_sizes}")
    print(f"   Growth: {'exponential' if r2_e > r2_p else 'polynomial'} (better R^2)")

    # Decision tree depths
    dt_vals = [(r['N'], r['dt_depth']) for r in results if r['dt_depth'] is not None]
    print(f"\n4. Decision tree depths: {dt_vals}")
    print(f"   (= minimum worst-case queries to determine primality)")

    # Comparison with random
    ratios = [r['bdd_size'] / r['rand_bdd_mean'] for r in results if r['rand_bdd_mean'] > 0]
    print(f"\n5. BDD size ratios (chi_P / random): {[f'{r:.2f}' for r in ratios]}")
    if all(r < 1.5 for r in ratios):
        print("   chi_P has SIMILAR complexity to random functions of same density.")
    elif all(r < 1 for r in ratios):
        print("   chi_P is SIMPLER than random functions (lower BDD size).")
    else:
        trend = "increasing" if ratios[-1] > ratios[0] else "decreasing"
        print(f"   Ratio trend: {trend}")


if __name__ == '__main__':
    main()
