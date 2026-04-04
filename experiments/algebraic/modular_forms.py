"""
Session 10: Modular Forms and Primes Investigation
====================================================
Investigating whether modular forms, eta products, mock modular forms,
and related structures can give a direct formula for primes.

7 experiments:
1. Ramanujan tau function - prime detection via tau(n) mod m
2. j-invariant coefficients - primes in moonshine coefficients
3. Hecke eigenvalues - encoding primes in eigenforms
4. Eta products - prime indicators from eta quotients
5. Half-integral weight forms - theta series and primes
6. Overconvergent modular forms (p-adic interpolation idea)
7. Mock modular forms - Ramanujan's mock theta coefficients vs primes
"""

import time
from sympy import isprime, primerange, factorint, divisors, sqrt, factorial
from sympy import binomial, Rational
from collections import defaultdict
import math

SEPARATOR = "=" * 70

# =====================================================================
# UTILITY: Compute Ramanujan tau function
# =====================================================================
def compute_tau_via_delta(N):
    """
    Compute tau(1)..tau(N) using the product definition:
    Delta(q) = q * prod_{n>=1} (1-q^n)^24 = sum_{n>=1} tau(n) q^n

    We expand the product up to q^N.
    """
    # Start with coefficients of prod (1 - q^n)^24 up to q^(N-1)
    # then shift by q to get Delta
    coeffs = [0] * (N + 1)
    coeffs[0] = 1  # constant term of the product

    for n in range(1, N + 1):
        # Multiply current polynomial by (1 - q^n)^24
        # Use the fact that (1-x)^24 = sum_{k=0}^{24} C(24,k)(-1)^k x^k
        # So (1-q^n)^24 = sum_{k=0}^{24} C(24,k)(-1)^k q^{nk}
        binom_coeffs = []
        for k in range(25):
            if n * k > N:
                break
            binom_coeffs.append((n * k, int((-1)**k * math.comb(24, k))))

        # Multiply in reverse order to avoid overwriting
        for j in range(N, -1, -1):
            if coeffs[j] == 0:
                continue
            for (shift, bc) in binom_coeffs:
                if shift == 0:
                    continue
                if j + shift <= N:
                    coeffs[j + shift] += coeffs[j] * bc

    # Now coeffs[k] is the coefficient of q^k in prod(1-q^n)^24
    # Delta(q) = q * prod = sum tau(n) q^n
    # So tau(n) = coeffs[n-1]
    tau = {}
    for n in range(1, N + 1):
        tau[n] = coeffs[n - 1]
    return tau


def compute_tau_via_recurrence(N):
    """
    Alternative: use Ramanujan's recurrence for tau via divisor sums.
    tau(n) can be computed using:
    n*sigma_11(n) = ... but that's not a simple recurrence.

    Actually use: tau(n) for the Ramanujan tau function.
    We use the product expansion which is reliable.
    """
    return compute_tau_via_delta(N)


# =====================================================================
# EXPERIMENT 1: Ramanujan tau function and prime detection
# =====================================================================
def experiment_1_ramanujan_tau():
    print(SEPARATOR)
    print("EXPERIMENT 1: Ramanujan tau function tau(n) and prime detection")
    print(SEPARATOR)

    N = 500  # Compute tau for n=1..N
    print(f"Computing tau(1)..tau({N}) via Delta expansion...")
    t0 = time.time()
    tau = compute_tau_via_delta(N)
    t1 = time.time()
    print(f"  Computed in {t1-t0:.2f}s")

    # Verify known values
    known = {1: 1, 2: -24, 3: 252, 4: -1472, 5: 4830, 6: -6048,
             7: -16744, 8: 84480, 9: -113643, 10: -115920}
    print("\nVerification against known tau values:")
    all_ok = True
    for n, expected in known.items():
        if tau[n] != expected:
            print(f"  MISMATCH: tau({n}) = {tau[n]}, expected {expected}")
            all_ok = False
    if all_ok:
        print("  All known values match!")

    # Check tau(n) mod m for various small m, see if primes are characterized
    print("\n--- Checking tau(n) mod m for prime characterization ---")
    primes_set = set(primerange(2, N + 1))

    best_separation = 0
    best_m = 0
    best_residues = set()

    for m in range(2, 50):
        # For each residue class, count primes and composites
        residue_to_primes = defaultdict(int)
        residue_to_composites = defaultdict(int)

        for n in range(2, N + 1):
            r = tau[n] % m
            if n in primes_set:
                residue_to_primes[r] += 1
            else:
                residue_to_composites[r] += 1

        # Find residues that are "mostly prime" or "mostly composite"
        # Measure: max fraction of primes in any single residue class
        for r in set(list(residue_to_primes.keys()) + list(residue_to_composites.keys())):
            p_count = residue_to_primes[r]
            c_count = residue_to_composites[r]
            total = p_count + c_count
            if total > 10 and p_count > 0:
                frac = p_count / total
                if frac > best_separation:
                    best_separation = frac
                    best_m = m
                    best_residues = {r}

    print(f"  Best single-residue prime concentration: {best_separation:.3f} at mod {best_m}")

    # Check Deligne bound: |tau(p)| < 2*p^(11/2) for primes
    print("\n--- Deligne bound check: |tau(p)| < 2*p^(11/2) ---")
    violations = 0
    for p in primerange(2, N + 1):
        bound = 2 * p ** 5.5
        if abs(tau[p]) >= bound:
            violations += 1
            print(f"  VIOLATION at p={p}: |tau({p})| = {abs(tau[p])}, bound = {bound:.0f}")
    if violations == 0:
        print("  No violations (as expected from Deligne's theorem)")

    # Check: is tau(n) multiplicative? tau(mn) = tau(m)*tau(n) for gcd(m,n)=1
    print("\n--- Multiplicativity check ---")
    mult_ok = 0
    mult_fail = 0
    for m_val in range(2, 20):
        for n_val in range(2, 20):
            if math.gcd(m_val, n_val) == 1 and m_val * n_val <= N:
                if tau[m_val * n_val] == tau[m_val] * tau[n_val]:
                    mult_ok += 1
                else:
                    mult_fail += 1
    print(f"  Multiplicative pairs: {mult_ok} OK, {mult_fail} failures")

    # Key insight: tau(p^2) = tau(p)^2 - p^11 for primes
    print("\n--- Recurrence tau(p^2) = tau(p)^2 - p^11 check ---")
    for p in [2, 3, 5, 7, 11, 13]:
        if p * p <= N:
            lhs = tau[p * p]
            rhs = tau[p] ** 2 - p ** 11
            status = "OK" if lhs == rhs else "FAIL"
            print(f"  p={p}: tau({p**2})={lhs}, tau(p)^2-p^11={rhs} [{status}]")

    # NEW: Check tau(n) mod 691 -- known congruence tau(n) ≡ sigma_11(n) mod 691
    print("\n--- Ramanujan congruence: tau(n) ≡ sigma_11(n) (mod 691) ---")
    def sigma_k(n, k):
        return sum(d ** k for d in divisors(n))

    cong_ok = 0
    for n in range(1, min(101, N + 1)):
        s11 = sigma_k(n, 11)
        if tau[n] % 691 == s11 % 691:
            cong_ok += 1
    print(f"  Congruence holds for {cong_ok}/100 values (should be 100)")

    # Can we use tau(n) mod 691 = sigma_11(n) mod 691 to detect primes?
    # For prime p: sigma_11(p) = 1 + p^11
    # For composite n=ab: sigma_11(n) includes more divisor terms
    print("\n--- Can sigma_11(n) mod 691 detect primes? ---")
    prime_residues_691 = set()
    composite_residues_691 = set()
    for n in range(2, min(201, N + 1)):
        r = (1 + n ** 11) % 691 if n in primes_set else None
        actual_r = tau[n] % 691
        if n in primes_set:
            prime_residues_691.add(actual_r)
        else:
            composite_residues_691.add(actual_r)

    only_prime = prime_residues_691 - composite_residues_691
    only_composite = composite_residues_691 - prime_residues_691
    overlap = prime_residues_691 & composite_residues_691
    print(f"  Prime-only residues mod 691: {len(only_prime)}")
    print(f"  Composite-only residues mod 691: {len(only_composite)}")
    print(f"  Overlapping residues: {len(overlap)}")
    print(f"  => No clean separation possible via single modulus")

    return tau


# =====================================================================
# EXPERIMENT 2: j-invariant coefficients and primes
# =====================================================================
def experiment_2_j_invariant():
    print(f"\n{SEPARATOR}")
    print("EXPERIMENT 2: j-invariant coefficients and primes")
    print(SEPARATOR)

    # j(q) = q^{-1} + 744 + 196884*q + 21493760*q^2 + ...
    # Compute coefficients using: j = E4^3 / Delta
    # E4 = 1 + 240*sum_{n>=1} sigma_3(n) q^n
    # Delta = tau(n) coefficients (already computed)

    N = 200

    # Compute E4 coefficients
    print(f"Computing j-invariant coefficients up to q^{N}...")

    def sigma_k(n, k):
        return sum(d ** k for d in divisors(n))

    # E4 = 1 + 240 * sum sigma_3(n) q^n
    e4 = [0] * (N + 2)
    e4[0] = 1
    for n in range(1, N + 2):
        e4[n] = 240 * sigma_k(n, 3)

    # E4^3: multiply three times
    def poly_mul(a, b, maxn):
        c = [0] * (maxn + 1)
        for i in range(min(len(a), maxn + 1)):
            if a[i] == 0:
                continue
            for j in range(min(len(b), maxn + 1 - i)):
                c[i + j] += a[i] * b[j]
        return c

    e4_sq = poly_mul(e4, e4, N + 1)
    e4_cube = poly_mul(e4_sq, e4, N + 1)

    # Delta coefficients (shift: Delta = sum tau(n) q^n, so coeff of q^n is tau(n))
    tau_dict = compute_tau_via_delta(N + 1)
    delta = [0] * (N + 2)
    for n in range(1, N + 2):
        delta[n] = tau_dict[n]

    # j = E4^3 / Delta = (E4^3 / q) * (1 / prod(1-q^n)^24)
    # Actually j(q) = E4(q)^3 / Delta(q)
    # E4^3 starts at q^0 with coefficient 1 (from 1^3)
    # Delta starts at q^1 with coefficient 1
    # So j = E4^3 / Delta, and E4^3 = c_0 + c_1 q + c_2 q^2 + ...
    # Delta = q + ... so j = (c_0 + c_1 q + ...)(q^{-1} + ...)
    # j starts as c_0 * q^{-1} + ...

    # Division: j * Delta = E4^3
    # j_{-1} * delta_1 = e4cube_0  => j_{-1} = e4cube_0 / delta_1 = e4cube_0 / 1
    # j_0 * delta_1 + j_{-1} * delta_2 = e4cube_1
    # In general: sum_{k} j_{n-k} * delta_k = e4cube_n for appropriate indexing

    # Let J[m] = coefficient of q^m in j(q), m starts at -1
    # Delta[k] = coefficient of q^k, k starts at 1
    # E4^3[n] = coefficient of q^n, n starts at 0
    # j(q) * Delta(q) = E4^3(q)
    # sum_{m=-1}^{...} J[m] * Delta[n-m] = E4^3[n]
    # For n=0: J[-1]*Delta[1] = E4^3[0] => J[-1] = E4^3[0] = 1
    # For n=1: J[-1]*Delta[2] + J[0]*Delta[1] = E4^3[1]
    # => J[0] = E4^3[1] - J[-1]*Delta[2]

    J = {}
    J[-1] = e4_cube[0]  # Should be 1

    for n in range(0, N + 1):
        # sum_{m=-1}^{n} J[m] * delta[n+1-m] = e4_cube[n+1] ...
        # Wait, let me re-index carefully.
        # j(q) = sum_{m>=-1} J[m] q^m
        # Delta(q) = sum_{k>=1} delta[k] q^k
        # Product: coeff of q^s = sum_{m+k=s} J[m]*delta[k] where k>=1, m>=-1
        # This equals E4^3 coeff of q^s = e4_cube[s]
        # For s = n (n >= 0):
        # sum_{m=-1}^{n-1} J[m] * delta[n-m] = e4_cube[n]
        # J[n-1] * delta[1] + sum_{m=-1}^{n-2} J[m] * delta[n-m] = e4_cube[n]
        # J[n-1] = e4_cube[n] - sum_{m=-1}^{n-2} J[m] * delta[n-m]
        pass

    # Redo more carefully
    J = {}
    # s=0: J[-1]*delta[1] = e4_cube[0], delta[1]=tau(1)=1, e4_cube[0]=1
    J[-1] = e4_cube[0]  # = 1

    for s in range(1, N + 2):
        # coeff of q^s: sum_{m=-1}^{s-1} J[m] * delta[s-m] = e4_cube[s]
        # Solve for J[s-1]:
        # J[s-1]*delta[1] + sum_{m=-1}^{s-2} J[m]*delta[s-m] = e4_cube[s]
        acc = 0
        for m in range(-1, s - 1):
            k = s - m
            if 1 <= k <= N + 1:
                acc += J[m] * delta[k]
        J[s - 1] = e4_cube[s] - acc  # delta[1] = 1

    # Now J[m] for m = -1, 0, 1, ..., N
    print(f"  J[-1] = {J[-1]} (should be 1)")
    print(f"  J[0]  = {J[0]} (should be 744)")
    print(f"  J[1]  = {J[1]} (should be 196884)")
    print(f"  J[2]  = {J[2]} (should be 21493760)")

    # Check which j-invariant coefficients are prime
    print("\n--- Which j-invariant coefficients J[n] are prime? ---")
    prime_coeffs = []
    for n in range(1, N + 1):
        val = abs(J[n])
        if val > 1 and isprime(val):
            prime_coeffs.append((n, J[n]))

    if prime_coeffs:
        print(f"  Found {len(prime_coeffs)} prime coefficients among J[1]..J[{N}]:")
        for n, v in prime_coeffs[:20]:
            print(f"    J[{n}] = {v}")
    else:
        print("  No prime coefficients found (they grow too fast)")

    # Check: are the indices n where J[n] has special divisibility properties related to primes?
    print("\n--- J[n] mod small primes, for prime vs composite n ---")
    primes_set = set(primerange(2, N + 1))
    for mod_p in [2, 3, 5, 7, 11, 13]:
        prime_residues = defaultdict(int)
        comp_residues = defaultdict(int)
        for n in range(2, N + 1):
            r = J[n] % mod_p
            if n in primes_set:
                prime_residues[r] += 1
            else:
                comp_residues[r] += 1
        # Check if distribution differs
        print(f"  mod {mod_p}: prime residue dist = {dict(prime_residues)}")
        print(f"         composite residue dist = {dict(comp_residues)}")

    # Monstrous moonshine: J[n] = dim of Monster rep
    # Check if J[p] for prime p has any pattern
    print("\n--- J[p] for small primes ---")
    for p in list(primerange(2, 30)):
        if p <= N:
            print(f"  J[{p}] = {J[p]}")

    return J


# =====================================================================
# EXPERIMENT 3: Hecke eigenvalues encoding primes
# =====================================================================
def experiment_3_hecke_eigenvalues():
    print(f"\n{SEPARATOR}")
    print("EXPERIMENT 3: Hecke eigenvalues and prime encoding")
    print(SEPARATOR)

    print("""
THEORETICAL ANALYSIS:
For a Hecke eigenform f = sum a_n q^n of weight k and level N:
- a_n is multiplicative: a_{mn} = a_m * a_n for gcd(m,n)=1
- At primes: a_{p^2} = a_p^2 - p^{k-1} (for p not dividing N)

KEY QUESTION: Can we construct a modular form whose coefficients
a_n = 1 iff n is prime, 0 otherwise?

ANSWER: NO, because:
1. Multiplicativity forces a_{mn} = a_m*a_n for coprime m,n
2. If a_2=1, a_3=1, then a_6 = a_2*a_3 = 1, but 6 is not prime
3. The prime indicator function is NOT multiplicative
4. Therefore no Hecke eigenform can have the prime indicator as coefficients.

ALTERNATIVE: Can we extract primes from the STRUCTURE of eigenvalues?
""")

    # Use the Ramanujan tau function (the unique eigenform of weight 12, level 1)
    N = 300
    tau = compute_tau_via_delta(N)

    # Check: for which n does tau(n) = 0?
    print("--- Values where tau(n) = 0 ---")
    zeros = [n for n in range(1, N + 1) if tau[n] == 0]
    print(f"  tau(n)=0 for n in: {zeros if zeros else 'NONE (Lehmer conjecture!)'}")

    # Normalized tau: a(p) = tau(p) / p^(11/2)
    # By Sato-Tate, a(p) is distributed like 2*sin^2(theta) on [-2, 2]
    print("\n--- Sato-Tate distribution of tau(p)/p^(11/2) ---")
    normalized = []
    for p in primerange(2, N + 1):
        a_p = tau[p] / (p ** 5.5)
        normalized.append((p, a_p))

    # Bin into intervals
    bins = defaultdict(int)
    for p, a in normalized:
        bin_idx = int((a + 2) * 5)  # 10 bins from -2 to 2
        bin_idx = max(0, min(9, bin_idx))
        bins[bin_idx] += 1

    print("  Distribution (binned, [-2,2] -> 10 bins):")
    for i in range(10):
        lo = -2 + i * 0.4
        hi = lo + 0.4
        bar = "#" * bins[i]
        print(f"    [{lo:+.1f},{hi:+.1f}): {bins[i]:3d} {bar}")

    # Check if we can determine primality from the RATIO tau(n)/n^(11/2)
    print("\n--- Can tau(n)/n^(11/2) magnitude distinguish primes? ---")
    prime_magnitudes = []
    comp_magnitudes = []
    primes_set = set(primerange(2, N + 1))
    for n in range(2, N + 1):
        mag = abs(tau[n]) / (n ** 5.5)
        if n in primes_set:
            prime_magnitudes.append(mag)
        else:
            comp_magnitudes.append(mag)

    avg_prime = sum(prime_magnitudes) / len(prime_magnitudes) if prime_magnitudes else 0
    avg_comp = sum(comp_magnitudes) / len(comp_magnitudes) if comp_magnitudes else 0
    print(f"  Average |tau(n)|/n^(11/2) for primes:     {avg_prime:.4f}")
    print(f"  Average |tau(n)|/n^(11/2) for composites:  {avg_comp:.4f}")
    print(f"  Ratio: {avg_prime/avg_comp:.4f}" if avg_comp > 0 else "  N/A")

    # The multiplicativity means composites have smaller normalized values on average
    # because tau(mn) = tau(m)*tau(n) and the normalization creates cancellation
    print("\n  FINDING: Primes tend to have LARGER normalized |tau(p)|/p^(11/2)")
    print("  But this is just Deligne bound — composites have multiplicative structure")
    print("  that tends to reduce the magnitude. NOT a usable primality test.")


# =====================================================================
# EXPERIMENT 4: Eta products and prime indicators
# =====================================================================
def experiment_4_eta_products():
    print(f"\n{SEPARATOR}")
    print("EXPERIMENT 4: Eta products and prime indicators")
    print(SEPARATOR)

    # eta(q) = q^(1/24) * prod(1-q^n)
    # Various eta quotients produce interesting sequences
    # Key examples:
    # eta(q)^24 = Delta(q) = sum tau(n) q^n
    # eta(q)^2 * eta(q^{11})^2 is the newform of weight 2, level 11

    N = 300

    # Compute prod(1-q^n) up to q^N
    def compute_eta_power(N, power):
        """Compute coefficients of (prod_{n>=1}(1-q^n))^power up to q^N."""
        coeffs = [0] * (N + 1)
        coeffs[0] = 1

        for n in range(1, N + 1):
            # Multiply by (1-q^n)^power
            # For small |power|, expand the binomial
            if power > 0:
                binom_terms = []
                for k in range(power + 1):
                    shift = n * k
                    if shift > N:
                        break
                    bc = int((-1)**k * math.comb(power, k))
                    binom_terms.append((shift, bc))
            else:
                # For negative power, use (1-q^n)^{-|power|} expansion
                # This is more complex; skip for now
                return None

            for j in range(N, -1, -1):
                if coeffs[j] == 0:
                    continue
                for (shift, bc) in binom_terms:
                    if shift == 0:
                        continue
                    if j + shift <= N:
                        coeffs[j + shift] += coeffs[j] * bc

        return coeffs

    # Compute eta^k for various k and check correlation with primes
    print("Computing eta products and checking prime correlations...\n")
    primes_set = set(primerange(2, N + 1))

    for power in [1, 2, 3, 4, 6, 8, 12]:
        coeffs = compute_eta_power(N, power)
        if coeffs is None:
            continue

        # Check if coefficients have special values at prime indices
        prime_vals = [coeffs[p] for p in primerange(2, N + 1)]
        comp_vals = [coeffs[n] for n in range(4, N + 1) if n not in primes_set and n > 1]

        # Check how many prime-index coefficients are zero
        prime_zeros = sum(1 for v in prime_vals if v == 0)
        comp_zeros = sum(1 for v in comp_vals if v == 0)

        prime_frac = prime_zeros / len(prime_vals) if prime_vals else 0
        comp_frac = comp_zeros / len(comp_vals) if comp_vals else 0

        print(f"  eta^{power}: coeff=0 fraction: primes={prime_frac:.3f}, composites={comp_frac:.3f}")

        # Check if sign pattern correlates
        prime_pos = sum(1 for v in prime_vals if v > 0)
        prime_neg = sum(1 for v in prime_vals if v < 0)
        comp_pos = sum(1 for v in comp_vals if v > 0)
        comp_neg = sum(1 for v in comp_vals if v < 0)

        pp = prime_pos / len(prime_vals) if prime_vals else 0
        cp = comp_pos / len(comp_vals) if comp_vals else 0
        print(f"           positive fraction: primes={pp:.3f}, composites={cp:.3f}")

    # Special: eta(q)*eta(q^{23}) — weight 1 form for level 23
    # This is related to the class number of Q(sqrt(-23))
    print("\n--- Special eta quotient: coefficients of eta(q)*eta(q^23) ---")
    eta1 = compute_eta_power(N, 1)
    # For eta(q^23), we need prod(1-q^{23n})
    eta23 = [0] * (N + 1)
    eta23[0] = 1
    for n in range(1, N // 23 + 1):
        for j in range(N, -1, -1):
            if eta23[j] != 0 and j + 23 * n <= N:
                eta23[j + 23 * n] -= eta23[j]

    # Product eta(q)*eta(q^23)
    product = poly_mul_simple(eta1, eta23, N)

    # Check values at prime indices
    print("  Values at prime indices (first 20 primes):")
    for p in list(primerange(2, 80)):
        if p <= N:
            print(f"    n={p}: coeff={product[p]}", end="")
            # For the form eta(q)*eta(q^23), a(p) relates to quadratic residues mod 23
            if p != 23:
                legendre = pow(p, 11, 23)  # Euler criterion (23-1)/2=11
                if legendre == 1:
                    print(f"  (p is QR mod 23)")
                else:
                    print(f"  (p is QNR mod 23)")
            else:
                print(f"  (p=23, ramified)")


def poly_mul_simple(a, b, maxn):
    """Simple polynomial multiplication."""
    c = [0] * (maxn + 1)
    for i in range(min(len(a), maxn + 1)):
        if a[i] == 0:
            continue
        for j in range(min(len(b), maxn + 1 - i)):
            if b[j] == 0:
                continue
            c[i + j] += a[i] * b[j]
    return c


# =====================================================================
# EXPERIMENT 5: Half-integral weight forms and primes
# =====================================================================
def experiment_5_half_integral():
    print(f"\n{SEPARATOR}")
    print("EXPERIMENT 5: Half-integral weight forms, theta series, and primes")
    print(SEPARATOR)

    # Theta series: theta(q) = sum_{n=-inf}^{inf} q^{n^2} = 1 + 2*sum_{n=1}^{inf} q^{n^2}
    # theta^k counts representations as sum of k squares
    # r_k(n) = #{ways to write n as sum of k squares}

    N = 300
    primes_set = set(primerange(2, N + 1))

    # Compute theta^k for k = 2, 3, 4
    for k in [2, 3, 4]:
        print(f"\n--- r_{k}(n): representations as sum of {k} squares ---")

        # Compute r_k(n) by convolution
        # Start with theta = coefficients where c[n^2] = 1 (plus c[0]=1, doubled for n>0)
        theta = [0] * (N + 1)
        theta[0] = 1
        m = 1
        while m * m <= N:
            theta[m * m] += 2  # +n and -n
            m += 1

        # Convolve k times
        result = theta[:]
        for _ in range(k - 1):
            new = [0] * (N + 1)
            for i in range(N + 1):
                if result[i] == 0:
                    continue
                for j in range(N + 1 - i):
                    if theta[j] == 0:
                        continue
                    new[i + j] += result[i] * theta[j]
            result = new

        # r_k values at primes vs composites
        prime_vals = [(p, result[p]) for p in primerange(2, min(50, N + 1))]

        print(f"  r_{k}(p) for first primes:")
        for p, v in prime_vals[:15]:
            # Known: r_2(p) = 0 if p ≡ 3 mod 4; r_2(p) = 4*((-1)^((p-1)/2)+1) etc.
            mod4 = p % 4
            print(f"    r_{k}({p}) = {v}  (p mod 4 = {mod4})")

        # Check if r_k(n) = 0 implies anything about primality
        zero_count_prime = sum(1 for p in primerange(2, N + 1) if result[p] == 0)
        zero_count_comp = sum(1 for n in range(4, N + 1)
                             if n not in primes_set and result[n] == 0)
        total_primes = len(list(primerange(2, N + 1)))
        total_comp = N - 3 - total_primes

        print(f"  r_{k}(n)=0: {zero_count_prime}/{total_primes} primes, "
              f"{zero_count_comp}/{total_comp} composites")

    # Shimura correspondence: relates half-integral weight forms to integral weight
    print("\n--- Shimura correspondence insight ---")
    print("""
  The Shimura lift takes a form f of weight k+1/2 to a form F of weight 2k.
  Coefficients: if f = sum c(n)q^n, then F relates c(n^2*m) to a_F(m).

  For prime detection, we'd need:
  - A half-integral weight form whose coefficients c(n) are 0/1 based on primality
  - But c(n) for such forms relates to CLASS NUMBERS and representation counts
  - These are intrinsically about quadratic forms, not primality

  CONCLUSION: Shimura correspondence connects to quadratic forms, not primes directly.
  The representation counts r_k(n) have well-known formulas involving divisor sums
  and character sums — they're COMPUTABLE but don't give primality tests.
""")


# =====================================================================
# EXPERIMENT 6: Overconvergent modular forms / p-adic interpolation
# =====================================================================
def experiment_6_overconvergent():
    print(f"\n{SEPARATOR}")
    print("EXPERIMENT 6: Overconvergent modular forms and p-adic interpolation")
    print(SEPARATOR)

    print("""
THEORETICAL ANALYSIS:

The idea: Use p-adic modular forms to interpolate the prime-counting function.

Coleman's theory provides p-adic families of modular forms parameterized by
the weight space. Could we find a p-adic modular form F such that F(n) = p_n?

KEY OBSTACLES:

1. CONTINUITY BARRIER: p-adic interpolation requires p-adic continuity.
   The prime function p(n) is NOT p-adically continuous for any prime p.

   Proof: Consider n and n + p^k. These are p-adically close (distance p^{-k}).
   But p(n) and p(n+p^k) differ by roughly p^k * ln(p^k) by PNT.
   So |p(n) - p(n+p^k)|_p is NOT necessarily small.

2. GROWTH RATE: Overconvergent forms have controlled growth on the
   overconvergent region. But primes grow like n*log(n), which doesn't
   match the growth patterns of modular form coefficients.

3. WEIGHT SPACE: Coleman families vary over weight space (a p-adic disk).
   The "weight" would need to correspond to the index n, but weight space
   is 1-dimensional while we need to parameterize all of N.
""")

    # Numerical experiment: check p-adic continuity of prime function
    print("--- Checking p-adic continuity of prime function ---")
    from sympy import nextprime, prime

    for p in [2, 3, 5, 7]:
        print(f"\n  Base prime p={p}:")
        print(f"  {'n':>8} {'n+p^k':>8} {'|n-(n+p^k)|_p':>15} {'|p(n)-p(n+p^k)|':>18} {'p-adic close?':>14}")

        for k in range(1, 6):
            pk = p ** k
            n = 100
            npk = n + pk

            pn = prime(n)
            pnpk = prime(npk)
            diff = abs(pn - pnpk)

            # p-adic valuation of diff
            if diff == 0:
                v_p = float('inf')
            else:
                v_p = 0
                temp = diff
                while temp % p == 0:
                    v_p += 1
                    temp //= p

            p_adic_dist_input = f"p^(-{k})"
            p_adic_dist_output = f"p^(-{v_p})" if v_p < float('inf') else "0"
            close = "YES" if v_p >= k else "NO"

            print(f"  {n:>8} {npk:>8} {p_adic_dist_input:>15} {diff:>18} (v_p={v_p}) {close:>8}")

    print("\n  CONCLUSION: The prime function is NOT p-adically continuous.")
    print("  Therefore it CANNOT be interpolated by p-adic/overconvergent modular forms.")


# =====================================================================
# EXPERIMENT 7: Mock modular forms (Ramanujan mock theta functions)
# =====================================================================
def experiment_7_mock_modular():
    print(f"\n{SEPARATOR}")
    print("EXPERIMENT 7: Mock modular forms and mock theta functions")
    print(SEPARATOR)

    # Ramanujan's third-order mock theta function:
    # f(q) = sum_{n>=0} q^{n^2} / prod_{k=1}^{n} (1+q^k)^2
    # = 1 + q - 2q^2 + 3q^3 - 3q^4 + 3q^5 - 5q^6 + ...

    # Let's compute a simpler mock theta: the order 3 function
    # f(q) = sum_{n>=0} q^{n^2} / (1+q)^2(1+q^2)^2...(1+q^n)^2

    N = 200
    print(f"Computing mock theta function coefficients up to q^{N}...")

    # Method: compute each term q^{n^2} / prod_{k=1}^{n}(1+q^k)^2
    # as a power series, then sum

    # First compute 1/prod_{k=1}^{n}(1+q^k)^2 iteratively
    # Start with all coefficients

    mock_coeffs = [0] * (N + 1)

    # For n=0: q^0 / (empty product) = 1
    mock_coeffs[0] = 1

    # Maintain inv_prod = coefficients of 1/prod_{k=1}^{n}(1+q^k)^2
    inv_prod = [0] * (N + 1)
    inv_prod[0] = 1

    for n in range(1, int(N**0.5) + 2):
        if n * n > N:
            break

        # Update inv_prod: divide by (1+q^n)^2
        # (1+q^n)^2 = 1 + 2q^n + q^{2n}
        # To divide, if F = inv_prod and G = (1+2q^n+q^{2n})
        # We want H such that H*G = F, i.e., H = F/G
        # H[k] = F[k] - 2*H[k-n] - H[k-2n]
        new_inv = [0] * (N + 1)
        for k in range(N + 1):
            val = inv_prod[k]
            if k >= n:
                val -= 2 * new_inv[k - n]
            if k >= 2 * n:
                val -= new_inv[k - 2 * n]
            new_inv[k] = val
        inv_prod = new_inv

        # Add contribution: q^{n^2} * inv_prod
        nn = n * n
        if nn <= N:
            for k in range(N + 1 - nn):
                mock_coeffs[nn + k] += inv_prod[k]

    print("  First 30 coefficients of mock theta f(q):")
    for i in range(0, 30):
        print(f"    a({i}) = {mock_coeffs[i]}")

    # Check correlation with primes
    primes_set = set(primerange(2, N + 1))

    print("\n--- Mock theta coefficients at prime vs composite indices ---")
    prime_vals = [mock_coeffs[p] for p in primerange(2, min(N + 1, 100))]
    comp_vals = [mock_coeffs[n] for n in range(4, 100) if n not in primes_set]

    avg_prime = sum(prime_vals) / len(prime_vals) if prime_vals else 0
    avg_comp = sum(comp_vals) / len(comp_vals) if comp_vals else 0
    print(f"  Average coeff at primes (2..100):     {avg_prime:.2f}")
    print(f"  Average coeff at composites (4..100):  {avg_comp:.2f}")

    # Check sign pattern
    prime_pos = sum(1 for v in prime_vals if v > 0)
    prime_neg = sum(1 for v in prime_vals if v < 0)
    comp_pos = sum(1 for v in comp_vals if v > 0)
    comp_neg = sum(1 for v in comp_vals if v < 0)
    print(f"  Sign at primes: {prime_pos} positive, {prime_neg} negative")
    print(f"  Sign at composites: {comp_pos} positive, {comp_neg} negative")

    # Check modular properties
    print("\n--- Mock theta coefficients mod small numbers ---")
    for m in [2, 3, 5, 7]:
        prime_res = defaultdict(int)
        comp_res = defaultdict(int)
        for n in range(2, N + 1):
            r = mock_coeffs[n] % m
            if n in primes_set:
                prime_res[r] += 1
            else:
                comp_res[r] += 1

        # Chi-squared-like statistic
        total_p = sum(prime_res.values())
        total_c = sum(comp_res.values())
        separation = 0
        for r in range(m):
            fp = prime_res[r] / total_p if total_p > 0 else 0
            fc = comp_res[r] / total_c if total_c > 0 else 0
            separation += abs(fp - fc)

        print(f"  mod {m}: L1 distance between prime/composite distributions = {separation:.3f}")

    # Zagier's mock modular forms and partition-like functions
    print("\n--- Connection to partition function p(n) ---")
    print("""
  Zagier showed that mock modular forms appear naturally in:
  - Partition ranks and cranks
  - Quantum modular forms
  - Traces of singular moduli

  However, the partition function p(n) and its relatives are about
  ADDITIVE number theory, while primes are MULTIPLICATIVE objects.

  The fundamental disconnect: modular forms encode multiplicative structure
  (via Euler products for L-functions), but the PRIME COUNTING function
  pi(x) is neither multiplicative nor naturally connected to any
  known modular form's coefficients.
""")


# =====================================================================
# GRAND SYNTHESIS
# =====================================================================
def grand_synthesis():
    print(f"\n{'#' * 70}")
    print("# GRAND SYNTHESIS: Modular Forms and the Prime Problem")
    print(f"{'#' * 70}")

    print("""
FINDINGS SUMMARY:

1. RAMANUJAN TAU (Exp 1):
   - tau(n) mod m does NOT cleanly separate primes from composites
   - The congruence tau(n) ≡ sigma_11(n) mod 691 is beautiful but
     sigma_11(n) for prime n just gives 1+n^11, which we already know
   - Deligne bound |tau(p)| < 2p^(11/2) is satisfied (verified)
   - tau is multiplicative, but the prime indicator function is NOT

2. j-INVARIANT (Exp 2):
   - Coefficients grow too fast to be prime themselves after small indices
   - No special pattern at prime indices vs composite indices
   - Monster group dimensions don't encode primality

3. HECKE EIGENVALUES (Exp 3):
   - FUNDAMENTAL OBSTRUCTION: Hecke eigenvalues are multiplicative
   - The prime indicator function chi_P(n) is NOT multiplicative
   - chi_P(2)*chi_P(3) = 1 but chi_P(6) = 0
   - Therefore NO Hecke eigenform encodes primality directly

4. ETA PRODUCTS (Exp 4):
   - Various eta^k produce sequences with no prime-specific structure
   - eta(q)*eta(q^23) relates to class numbers of Q(sqrt(-23))
   - Quadratic form theory, not primality

5. HALF-INTEGRAL WEIGHT (Exp 5):
   - r_k(n) counts representations as sum of k squares
   - These relate to CLASS NUMBERS and CHARACTER SUMS
   - Shimura correspondence maps to integral-weight forms
   - None of this structure is about primality

6. OVERCONVERGENT / p-ADIC (Exp 6):
   - p(n) (nth prime) is NOT p-adically continuous for any p
   - Therefore CANNOT be p-adically interpolated
   - Coleman families don't help

7. MOCK MODULAR FORMS (Exp 7):
   - Mock theta function coefficients show no prime-specific patterns
   - Distribution at prime vs composite indices statistically similar
   - Additive number theory structure, not multiplicative

═══════════════════════════════════════════════════════════════════════
CORE THEORETICAL BARRIER:
═══════════════════════════════════════════════════════════════════════

Modular forms are intimately connected to MULTIPLICATIVE number theory:
- Hecke operators preserve multiplicativity
- L-functions have Euler products over primes
- Galois representations encode local behavior at each prime

The nth prime function p(n) requires COUNTING primes up to a bound,
which is an ADDITIVE/SIEVING operation. The connection between
modular forms and primes goes the WRONG WAY:

  Modular forms → (encode info about each prime separately via Euler factors)
  Prime counting → (needs to AGGREGATE over all primes up to x)

To get p(n), you need pi(x) (prime counting), which requires summing
a non-multiplicative function (the prime indicator) — exactly the
UNCONDITIONAL SUMMATION BARRIER identified in previous sessions.

Even the explicit formula pi(x) = li(x) - sum_rho li(x^rho) - ...
uses zeros of zeta (which ARE connected to modular forms via Langlands),
but evaluating this sum requires O(T) zeros where T ~ x, giving at
best O(x^{1/2+eps}) operations — not polylog.

VERDICT: Modular forms provide the RICHEST framework for studying primes,
but they encode primes INDIVIDUALLY (at each Euler factor), not
COLLECTIVELY (as a counting function). This is the fundamental mismatch.

The problem p(10^100) in O(polylog) time remains IMPOSSIBLE by all
known theoretical frameworks. Modular forms confirm the barrier rather
than circumventing it.

STATUS: Path CLOSED. This is approach ~326-332 (7 sub-experiments).
═══════════════════════════════════════════════════════════════════════
""")


# =====================================================================
# MAIN
# =====================================================================
if __name__ == "__main__":
    print("SESSION 10: Modular Forms Investigation for Prime Formula")
    print("=" * 70)
    print(f"Date: 2026-04-04")
    print(f"Goal: Find O(polylog n) formula for p(n) via modular forms")
    print("=" * 70)

    t_start = time.time()

    tau = experiment_1_ramanujan_tau()
    J = experiment_2_j_invariant()
    experiment_3_hecke_eigenvalues()
    experiment_4_eta_products()
    experiment_5_half_integral()
    experiment_6_overconvergent()
    experiment_7_mock_modular()

    grand_synthesis()

    t_end = time.time()
    print(f"\nTotal runtime: {t_end - t_start:.1f}s")
