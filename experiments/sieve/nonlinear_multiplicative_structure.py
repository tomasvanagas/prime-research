"""
Experiment: Multiplicative structure in floor values.

Key identity: floor(floor(x/a)/b) = floor(x/(ab)) always.
So there's no "correction" from iterated floor division.

But: floor(x/a) * floor(x/b) != floor(x*floor(x/b)/a) in general.
And: floor(x/a) mod floor(x/b) carries information.

We investigate:
1. The lattice structure of {floor(x/k)} and its relation to primes
2. Nonlinear combinations that might extract pi(x) with fewer operations
3. Whether products/GCDs of floor values create useful "sieve" effects
"""

import numpy as np
from sympy import primepi, isprime, primerange, factorint, gcd
import math
from collections import Counter

# ============================================================
# Experiment 1: GCD structure of floor values
# ============================================================
print("=" * 70)
print("EXPERIMENT 1: GCD of pairs of floor values")
print("=" * 70)

for x in [100, 500, 1000]:
    pi_x = primepi(x)
    sqrt_x = int(math.sqrt(x))

    # Compute floor(x/k) for k=1..sqrt(x)
    fv = [x // k for k in range(1, sqrt_x + 1)]

    # GCD of consecutive floor values
    gcds = [math.gcd(fv[i], fv[i+1]) for i in range(len(fv)-1)]

    # How often does GCD = 1?
    coprime_count = sum(1 for g in gcds if g == 1)

    print(f"\nx={x}, sqrt(x)={sqrt_x}, pi(x)={pi_x}")
    print(f"  Consecutive GCDs: {gcds[:20]}...")
    print(f"  Coprime consecutive pairs: {coprime_count}/{len(gcds)}")
    print(f"  Sum of GCDs: {sum(gcds)}")

    # GCD of floor(x/k) and k
    gcd_with_k = [math.gcd(x // k, k) for k in range(1, sqrt_x + 1)]
    prime_gcd = [math.gcd(x // p, p) for p in primerange(2, sqrt_x + 1)]
    comp_gcd = [math.gcd(x // k, k) for k in range(4, sqrt_x + 1) if not isprime(k)]

    print(f"  GCD(floor(x/k), k) at primes: {prime_gcd[:15]}")
    print(f"  GCD(floor(x/k), k) at composites: {comp_gcd[:15]}")

    # For prime p: gcd(floor(x/p), p) = p if p|x, else gcd(floor(x/p), p)
    # Since p is prime, this is either 1 or p.
    prime_divides = sum(1 for g in prime_gcd if g > 1)
    comp_divides = sum(1 for g in comp_gcd if g > 1)
    print(f"  GCD > 1 at primes: {prime_divides}/{len(prime_gcd)}")
    print(f"  GCD > 1 at composites: {comp_divides}/{len(comp_gcd)}")

# ============================================================
# Experiment 2: Floor value collisions and prime counting
# ============================================================
print("\n" + "=" * 70)
print("EXPERIMENT 2: Floor value collisions")
print("=" * 70)

# Two values k, k' have a "collision" if floor(x/k) = floor(x/k')
# The collision structure is intimately related to the Dirichlet hyperbola

for x in [100, 500, 1000]:
    pi_x = primepi(x)

    # Group k values by their floor(x/k)
    groups = {}
    for k in range(1, x + 1):
        v = x // k
        if v not in groups:
            groups[v] = []
        groups[v].append(k)

    # Statistics on group sizes
    sizes = [len(g) for g in groups.values()]

    print(f"\nx={x}: pi(x)={pi_x}")
    print(f"  Distinct floor values: {len(groups)}")
    print(f"  Group size distribution: min={min(sizes)}, max={max(sizes)}, mean={np.mean(sizes):.1f}")

    # Groups containing at least one prime
    prime_groups = sum(1 for g in groups.values() if any(isprime(k) for k in g))
    # Groups that are entirely prime
    all_prime_groups = sum(1 for g in groups.values() if all(isprime(k) for k in g) and min(g) >= 2)

    print(f"  Groups with >= 1 prime: {prime_groups}")
    print(f"  Groups that are all-prime: {all_prime_groups}")

    # For singleton groups (size 1): how many are prime?
    singletons = [g[0] for g in groups.values() if len(g) == 1]
    singleton_primes = sum(1 for k in singletons if isprime(k) and k >= 2)
    print(f"  Singleton groups: {len(singletons)}, of which prime: {singleton_primes}")

    # Key insight check: each "group" corresponds to a range k in [x/(v+1)+1, x/v]
    # The number of primes in this range is exactly pi(x/v) - pi(x/(v+1))
    # Summing these over v gives pi(x).
    # This is just the standard decomposition, nothing new.

# ============================================================
# Experiment 3: Quadratic sieve analog
# ============================================================
print("\n" + "=" * 70)
print("EXPERIMENT 3: Quadratic sieve analog for counting primes")
print("=" * 70)

# In the quadratic sieve for factoring N:
# Find y such that y^2 mod N is smooth
# Combine smooth relations to get x^2 = y^2 (mod N)

# Analog for counting: can we find "relations" among floor values
# that reveal prime structure?

# Idea: floor(x/a) * floor(x/b) mod floor(x/c) for various a,b,c
# If this detects something about divisibility by primes in [a,b,c]...

for x in [100, 200, 500]:
    pi_x = primepi(x)
    sqrt_x = int(math.sqrt(x))

    # Build "relation matrix": R[i,j] = floor(x/i) mod floor(x/j) for small i,j
    K = min(20, sqrt_x)
    fv = {k: x // k for k in range(1, K + 1)}

    print(f"\nx={x}, pi(x)={pi_x}")

    # Check: floor(x/p) mod p for prime p
    for p in primerange(2, K + 1):
        val = (x // p) % p
        remainder = x % p
        print(f"  floor({x}/{p}) = {x // p}, mod {p} = {val}, x mod {p} = {remainder}")
        # Note: floor(x/p) = (x - (x mod p)) / p
        # So floor(x/p) mod p = ((x - r) / p) mod p where r = x mod p

    # Quadratic relation: floor(x/a)^2 - x*floor(x/a)/a
    # = floor(x/a) * (floor(x/a) - x/a) = floor(x/a) * (-{x/a})
    # where {x/a} is the fractional part. This is always <= 0.

    # Try: sum_{p prime} floor(x/p)^2 vs sum_{k} floor(x/k)^2
    sum_sq_prime = sum((x // p) ** 2 for p in primerange(2, x + 1))
    sum_sq_all = sum((x // k) ** 2 for k in range(2, x + 1))

    # Dirichlet series analog: sum_{k} floor(x/k)^2 / k^s
    print(f"  Sum floor(x/p)^2 over primes: {sum_sq_prime}")
    print(f"  Sum floor(x/k)^2 over all k: {sum_sq_all}")
    print(f"  Ratio: {sum_sq_prime / sum_sq_all:.6f}")

# ============================================================
# Experiment 4: Nonlinear Mobius-like inversion
# ============================================================
print("\n" + "=" * 70)
print("EXPERIMENT 4: Nonlinear Mobius-like formulas")
print("=" * 70)

# Standard: pi(x) = sum_{n<=x} Lambda(n)/log(n) (approximately)
# where Lambda(n) = log(p) if n=p^k, 0 otherwise
# And: psi(x) = sum_{n<=x} Lambda(n) = x - sum_rho x^rho/rho - ...

# Nonlinear variant: Can we express pi(x) using products of Mobius-type sums?
# E.g., pi(x) = f(M(x), M(x/2), M(x/3), ...) for some nonlinear f?

# Where M(x) = sum_{n<=x} mu(n) = Mertens function

# First, let's see what information M(x/k) values carry
from sympy import mobius

for x in [50, 100, 200]:
    pi_x = primepi(x)
    sqrt_x = int(math.sqrt(x))

    # Compute M(y) = sum_{n<=y} mu(n) for y = x/k
    def mertens(y):
        return sum(mobius(n) for n in range(1, int(y) + 1))

    K_max = min(30, sqrt_x + 1)
    M_vals = {k: mertens(x / k) for k in range(1, K_max)}

    print(f"\nx={x}, pi(x)={pi_x}")
    print(f"  M(x/k) for k=1..{min(10,K_max-1)}: {[M_vals[k] for k in range(1, min(11, K_max))]}")

    # Try: pi(x) = sum_k c_k * M(x/k) (this is standard Mobius inversion)
    # pi(x) = sum_{k<=x} M(x/k) * (something about log)
    # Actually: sum_{k=1}^{x} M(x/k) = 1 (always)

    check = sum(mertens(x / k) for k in range(1, x + 1))
    print(f"  sum_k M(x/k) = {check} (should be 1)")

    # Now nonlinear: M(x/a) * M(x/b)
    # Can we fit pi(x) = sum_{a,b} w_{a,b} * M(x/a) * M(x/b)?
    K = 8
    features = []
    for a in range(1, K + 1):
        for b in range(a, K + 1):
            features.append(M_vals.get(a, 0) * M_vals.get(b, 0))

    print(f"  M(x/a)*M(x/b) features (a,b<=8): {features[:10]}...")

# Fit across multiple x values
print("\n  Fitting pi(x) = sum w_{a,b} * M(x/a)*M(x/b) ...")
from sympy import mobius

K = 6
X_range = list(range(10, 150))

def mertens_fast(limit):
    """Compute Mertens function for all values up to limit"""
    mu = [0] * (limit + 1)
    mu[1] = 1
    for i in range(1, limit + 1):
        for j in range(2 * i, limit + 1, i):
            mu[j] -= mu[i]
    # Prefix sum
    M = [0] * (limit + 1)
    for i in range(1, limit + 1):
        M[i] = M[i - 1] + mu[i]
    return M

M_table = mertens_fast(max(X_range) + 1)

feature_matrix = []
target = []
for x in X_range:
    M_vals = [M_table[x // k] for k in range(1, K + 1)]
    row = []
    row.append(1)  # constant
    for a in range(K):
        row.append(M_vals[a])  # linear
    for a in range(K):
        for b in range(a, K):
            row.append(M_vals[a] * M_vals[b])  # quadratic
    feature_matrix.append(row)
    target.append(primepi(x))

A = np.array(feature_matrix, dtype=float)
y = np.array(target, dtype=float)
coeffs, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
pred = A @ coeffs
max_err = np.max(np.abs(pred - y))
exact = sum(1 for i in range(len(X_range)) if abs(round(pred[i]) - y[i]) < 0.01)
print(f"  K={K}, features={A.shape[1]}, max_err={max_err:.4f}, exact={exact}/{len(X_range)}")

# Test generalization
X_test = list(range(150, 300))
M_table2 = mertens_fast(max(X_test) + 1)
feat_test = []
target_test = []
for x in X_test:
    M_vals = [M_table2[x // k] for k in range(1, K + 1)]
    row = [1]
    for a in range(K):
        row.append(M_vals[a])
    for a in range(K):
        for b in range(a, K):
            row.append(M_vals[a] * M_vals[b])
    feat_test.append(row)
    target_test.append(primepi(x))

A_test = np.array(feat_test, dtype=float)
y_test = np.array(target_test, dtype=float)
pred_test = A_test @ coeffs
max_err_test = np.max(np.abs(pred_test - y_test))
exact_test = sum(1 for i in range(len(X_test)) if abs(round(pred_test[i]) - y_test[i]) < 0.01)
print(f"  Generalization (150-299): max_err={max_err_test:.4f}, exact={exact_test}/{len(X_test)}")

print("\nDone.")
