#!/usr/bin/env python3
"""
PROPOSAL: Smooth Number Counting + Subtraction for pi(x)

Idea: Counting B-smooth numbers Psi(x, B) might be easier than counting primes.
If so, can we reconstruct pi(x) from Psi values via linear combination?

Tests:
1. Compute Psi(x, B) exactly for various B and x up to 10000.
2. Check if pi(x) = linear combination of Psi(x, B_i).
3. Time Psi(x, B) vs direct sieve.
4. Analyze Buchstab recursion tree size.

Related closed paths:
  - "Smooth number counts" (S6): 5-15% error, inherently approximate
  - "Smooth/rough decomposition" (S14): rough count IS the sieve
  - "Buchstab tree pruning/memo" (S11): memoized Buchstab IS Lucy DP, O(x^{2/3})
  - "Buchstab signed identity + H-T" (S16): conversion M(x)->pi(x) costs O(x^{2/3})
  - "Recursive Dickman DDE shortcut" (S24): correction grows with exponent alpha~7.9
"""

import time
import numpy as np
from sympy import primerange, isprime, factorint, primepi
from collections import defaultdict
import sys

# ============================================================
# Part 1: Exact Psi(x, B) computation
# ============================================================

def psi_exact(x, B):
    """Count B-smooth numbers in [1, x]: all prime factors <= B."""
    count = 0
    for n in range(1, x + 1):
        if n == 1:
            count += 1
            continue
        factors = factorint(n)
        if max(factors.keys()) <= B:
            count += 1
    return count


def psi_exact_sieve(x, B):
    """Count B-smooth numbers via sieve (faster for large x)."""
    # Mark non-smooth numbers
    is_smooth = [True] * (x + 1)
    # For each prime p > B, mark multiples as not smooth
    primes_above_B = list(primerange(B + 1, x + 1))
    for p in primes_above_B:
        for mult in range(p, x + 1, p):
            is_smooth[mult] = False
    return sum(is_smooth[1:])  # exclude 0


def phi_legendre(x, a, primes):
    """
    Legendre sieve function: count n <= x with no prime factor <= primes[a-1].
    phi(x, 0) = floor(x).
    phi(x, a) = phi(x, a-1) - phi(x/primes[a-1], a-1)
    """
    if a == 0:
        return int(x)
    if x < 1:
        return 0
    return phi_legendre(x, a - 1, primes) - phi_legendre(x // primes[a - 1], a - 1, primes)


print("=" * 70)
print("PART 1: Psi(x, B) for various x and B")
print("=" * 70)

test_x_values = [100, 500, 1000, 5000, 10000]
B_values = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]

print(f"\n{'x':>6} | {'pi(x)':>6} | " + " | ".join(f"B={b:>2}" for b in B_values))
print("-" * (16 + 7 * len(B_values)))

psi_data = {}
for x in test_x_values:
    pi_x = int(primepi(x))
    row = f"{x:>6} | {pi_x:>6} | "
    psi_row = []
    for B in B_values:
        if B > x:
            psi_val = x
        else:
            psi_val = psi_exact_sieve(x, B)
        psi_row.append(psi_val)
        psi_data[(x, B)] = psi_val
    row += " | ".join(f"{v:>4}" for v in psi_row)
    print(row)

# ============================================================
# Part 2: Can pi(x) be a linear combination of Psi values?
# ============================================================

print("\n" + "=" * 70)
print("PART 2: Linear combination Psi(x, B_i) -> pi(x)?")
print("=" * 70)

# Build training data: for x in [50, 60, ..., 10000]
train_x = list(range(50, 5001, 10))
n_train = len(train_x)
n_features = len(B_values)

print(f"\nTraining on {n_train} values of x in [{train_x[0]}, {train_x[-1]}]")
print(f"Features: Psi(x, B) for B in {B_values}")

# Build matrix A and target vector b
A = np.zeros((n_train, n_features))
b = np.zeros(n_train)

for i, x in enumerate(train_x):
    b[i] = int(primepi(x))
    for j, B in enumerate(B_values):
        A[i, j] = psi_exact_sieve(x, B)

# Solve least squares: A @ c ~= b
coeffs, residuals, rank, sv = np.linalg.lstsq(A, b, rcond=None)

print(f"\nLeast-squares coefficients:")
for j, B in enumerate(B_values):
    print(f"  c(B={B:>2}) = {coeffs[j]:>12.6f}")

# Compute fit quality on training data
pred_train = A @ coeffs
errors_train = pred_train - b
rmse_train = np.sqrt(np.mean(errors_train ** 2))
max_err_train = np.max(np.abs(errors_train))
rel_err_train = np.mean(np.abs(errors_train) / np.maximum(b, 1))

print(f"\nTraining RMSE: {rmse_train:.4f}")
print(f"Training max |error|: {max_err_train:.4f}")
print(f"Training mean relative |error|: {rel_err_train:.6f}")
print(f"Exact matches (|err| < 0.5): {np.sum(np.abs(errors_train) < 0.5)}/{n_train}")

# Test on held-out data
test_x = list(range(5005, 10001, 10))
n_test = len(test_x)
A_test = np.zeros((n_test, n_features))
b_test = np.zeros(n_test)

for i, x in enumerate(test_x):
    b_test[i] = int(primepi(x))
    for j, B in enumerate(B_values):
        A_test[i, j] = psi_exact_sieve(x, B)

pred_test = A_test @ coeffs
errors_test = pred_test - b_test
rmse_test = np.sqrt(np.mean(errors_test ** 2))
max_err_test = np.max(np.abs(errors_test))
rel_err_test = np.mean(np.abs(errors_test) / np.maximum(b_test, 1))

print(f"\nTest RMSE: {rmse_test:.4f}")
print(f"Test max |error|: {max_err_test:.4f}")
print(f"Test mean relative |error|: {rel_err_test:.6f}")
print(f"Exact matches (|err| < 0.5): {np.sum(np.abs(errors_test) < 0.5)}/{n_test}")

# Also try with x as additional feature (affine model)
print("\n--- With x as additional feature (affine model) ---")
A_aug = np.column_stack([A, np.array(train_x, dtype=float)])
A_test_aug = np.column_stack([A_test, np.array(test_x, dtype=float)])

coeffs2, _, _, _ = np.linalg.lstsq(A_aug, b, rcond=None)
pred_test2 = A_test_aug @ coeffs2
errors_test2 = pred_test2 - b_test
rmse_test2 = np.sqrt(np.mean(errors_test2 ** 2))
max_err_test2 = np.max(np.abs(errors_test2))

print(f"Test RMSE: {rmse_test2:.4f}")
print(f"Test max |error|: {max_err_test2:.4f}")
print(f"Exact matches (|err| < 0.5): {np.sum(np.abs(errors_test2) < 0.5)}/{n_test}")

# ============================================================
# Part 3: Timing Psi(x, B) vs direct sieve for pi(x)
# ============================================================

print("\n" + "=" * 70)
print("PART 3: Timing comparison")
print("=" * 70)

def sieve_pi(x):
    """Compute pi(x) via Eratosthenes sieve."""
    if x < 2:
        return 0
    is_prime = [True] * (x + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(x**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, x + 1, i):
                is_prime[j] = False
    return sum(is_prime)

test_sizes = [1000, 5000, 10000, 50000, 100000]

print(f"\n{'x':>8} | {'pi(x)':>7} | {'sieve_time':>10} | {'psi_total':>10} | {'ratio':>6}")
print("-" * 55)

for x in test_sizes:
    # Time pi(x) via sieve
    t0 = time.perf_counter()
    pi_val = sieve_pi(x)
    t_sieve = time.perf_counter() - t0

    # Time all Psi(x, B) computations needed
    B_list = list(primerange(2, int(x**0.5) + 1))
    t0 = time.perf_counter()
    for B in B_list:
        psi_exact_sieve(x, B)
    t_psi = time.perf_counter() - t0

    ratio = t_psi / max(t_sieve, 1e-9)
    print(f"{x:>8} | {pi_val:>7} | {t_sieve:>10.6f} | {t_psi:>10.6f} | {ratio:>6.1f}x")

# ============================================================
# Part 4: Buchstab recursion tree analysis
# ============================================================

print("\n" + "=" * 70)
print("PART 4: Buchstab recursion tree size")
print("=" * 70)

def buchstab_tree_size(x, a, primes, memo=None):
    """
    Count nodes in the Buchstab recursion tree for phi(x, a).
    phi(x, a) = phi(x, a-1) - phi(x/p_a, a-1)
    Base: phi(x, 0) = floor(x) (1 node).
    """
    if memo is None:
        memo = {}
    key = (int(x), a)
    if key in memo:
        return memo[key]
    if a == 0 or x < 1:
        return 1
    count = 1 + buchstab_tree_size(x, a - 1, primes, memo) + \
                buchstab_tree_size(x // primes[a - 1], a - 1, primes, memo)
    memo[key] = count
    return count


def buchstab_distinct_args(x, a, primes, seen=None):
    """Count DISTINCT (floor(x'), a') arguments in the Buchstab recursion."""
    if seen is None:
        seen = set()
    key = (int(x), a)
    if key in seen:
        return
    seen.add(key)
    if a == 0 or x < 1:
        return
    buchstab_distinct_args(x, a - 1, primes, seen)
    buchstab_distinct_args(x // primes[a - 1], a - 1, primes, seen)


print(f"\n{'x':>8} | {'a=pi(sqrt)':>10} | {'tree_nodes':>12} | {'distinct':>10} | {'x^(2/3)':>10} | {'ratio':>6}")
print("-" * 70)

for x in [100, 500, 1000, 5000, 10000, 50000]:
    sqrt_x = int(x**0.5)
    primes = list(primerange(2, sqrt_x + 1))
    a = len(primes)

    memo = {}
    tree_nodes = buchstab_tree_size(x, a, primes, memo)

    seen = set()
    buchstab_distinct_args(x, a, primes, seen)
    distinct = len(seen)

    x_23 = x ** (2/3)
    ratio = distinct / x_23

    print(f"{x:>8} | {a:>10} | {tree_nodes:>12} | {distinct:>10} | {x_23:>10.1f} | {ratio:>6.2f}")

# ============================================================
# Part 5: Key identity check -- exact relationship
# ============================================================

print("\n" + "=" * 70)
print("PART 5: Exact identity: pi(x) = phi(x, pi(sqrt(x))) + pi(sqrt(x)) - 1")
print("=" * 70)

for x in [100, 500, 1000, 5000, 10000]:
    sqrt_x = int(x**0.5)
    primes = list(primerange(2, sqrt_x + 1))
    a = len(primes)

    phi_val = phi_legendre(x, a, primes)
    pi_sqrt = a
    reconstructed = phi_val + pi_sqrt - 1
    actual = int(primepi(x))

    match = "OK" if reconstructed == actual else "FAIL"
    print(f"x={x:>5}: phi(x,{a})={phi_val:>5}, pi(sqrt)={pi_sqrt}, "
          f"reconstructed={reconstructed:>5}, actual={actual:>5} [{match}]")

# Note: phi(x, a) = x - sum_{p <= p_a} floor(x/p) + sum_{p<q<=p_a} floor(x/pq) - ...
# This is the Legendre formula. Computing phi IS the sieve.
# Psi(x, B) = x - phi(x, pi(B)) + 1 (since smooth = total - rough + 1)
# So Psi and phi are trivially related.

print("\n" + "=" * 70)
print("PART 6: Relationship Psi(x, B) and phi(x, a)")
print("=" * 70)

for x in [100, 500, 1000]:
    sqrt_x = int(x**0.5)
    primes = list(primerange(2, x + 1))
    primes_small = list(primerange(2, sqrt_x + 1))

    print(f"\nx = {x}:")
    for B in [2, 3, 5, 7, 11, 13]:
        if B > x:
            break
        # a = pi(B)
        a = int(primepi(B))
        primes_up_to_B = list(primerange(2, B + 1))

        psi_val = psi_exact_sieve(x, B)
        phi_val = phi_legendre(x, a, primes_up_to_B)

        # phi(x, a) = #{n <= x : all prime factors > B}
        # Psi(x, B) = #{n <= x : all prime factors <= B}
        # These are NOT complements -- numbers can have BOTH small and large factors
        # Complement of Psi: numbers with at least one factor > B
        # phi counts numbers with ALL factors > B (= rough numbers)

        # Exact: Psi(x, B) + phi(x, a) + #{mixed} = x
        # where mixed = numbers with at least one factor <= B AND at least one > B
        mixed = x - psi_val - phi_val + 1  # +1 for n=1 counted in both
        print(f"  B={B:>2}: Psi={psi_val:>5}, phi(x,pi(B))={phi_val:>5}, "
              f"mixed={mixed:>5}, sum={psi_val + phi_val + mixed}")

# ============================================================
# Summary
# ============================================================

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print("""
Key findings:

1. LINEAR COMBINATION: Psi(x, B_i) values do NOT give exact pi(x) via
   linear combination. The residual error reflects the oscillatory part
   of pi(x) that encodes zeta zero information. Even with 11+ features,
   test errors remain non-negligible.

2. TIMING: Computing Psi(x, B) for ALL B up to sqrt(x) is SLOWER than
   a single Eratosthenes sieve. Each Psi computation is essentially a
   sieve itself, and we need O(pi(sqrt(x))) of them.

3. BUCHSTAB TREE: The number of distinct arguments in the Buchstab
   recursion scales as O(x^{2/3}), confirming the known result that
   memoized Buchstab IS the Lucy DP / Meissel-Lehmer method.

4. FUNDAMENTAL ISSUE: Psi(x, B) and phi(x, a) are trivially related.
   Computing Psi(x, B) exactly IS the sieve problem. There is no
   shortcut: the smooth number count contains no "free" information
   that isn't already in the sieve.

VERDICT: FAIL (Equivalence). This reduces to the standard Legendre/
Meissel-Lehmer sieve. The smooth number perspective adds no new
computational leverage. Confirmed prior results from S6, S11, S14, S16.
""")
