"""
FOCUS-1 Sub-attack 2 (TODO.md): Non-cyclotomic ring AKS variant.

Standard AKS works in Z_n[x]/(x^r - 1) and the proof exploits the cyclotomic
factorization to embed an "introspective monoid" into F_p[x]/h(x) where h is
an irreducible factor of Phi_r(x) mod p of degree > log^2(n).

S47 (CLOSED_PATHS line 690) closed the cyclotomic-CRT shortcut: when r is the
AKS-prescribed prime, x^r - 1 factors only as (x-1)*Phi_r(x), so the working
sub-ring still has dimension r-1. To break out of cyclotomic-trivial-factor
bind, we test whether AKS congruence

    (x + b)^n  ≡  x^n + b   (mod f(x), n)

works in a NON-CYCLOTOMIC ring Z_n[x]/f(x), with f(x) = x^d + a (Eisenstein-
flavoured) of degree d ~ log^2(n) — much smaller than r.

This script measures EMPIRICALLY:
  Q1  Fraction of composites in [3, M] that pass the congruence for given d,a,b.
  Q2  How does the false-positive rate scale with d (the polylog parameter)?
  Q3  Is there a (d, a) choice where O(1) values of b distinguish primes
      from composites with vanishing false-positive rate?
  Q4  Cost of certifying f(x) = x^d + a is irreducible mod p (p|n).
      (Rabin's irreducibility test is in NC^2 — too deep for TC^0.
      For x^d + a the criterion is simpler: irreducible over F_p iff
      every prime divisor q of d satisfies (-a) is not a q-th power in F_p,
      and (if 4 | d) -4a is not a 4th power. — Capelli's theorem.)

Construction discipline (CLAUDE.md S60 rules): produce code that runs;
report findings under construction-incoherent / FAIL / PARTIAL category.
"""

import math
import sys
from sympy import isprime, primerange, factorint

# ------------------------------------------------------------
# Polynomial arithmetic in Z_n[x] / (x^d + a)
# Polynomials represented as length-d lists of ints (coeff of x^i at index i).
# Reduction: x^d  ≡  -a  (mod f).
# ------------------------------------------------------------

def poly_mul_mod(p, q, d, a, n):
    """Multiply two polynomials of degree < d, reduce mod (x^d + a, n)."""
    out = [0] * (2 * d - 1)
    for i, pi in enumerate(p):
        if pi == 0:
            continue
        for j, qj in enumerate(q):
            if qj == 0:
                continue
            out[i + j] = (out[i + j] + pi * qj) % n
    # reduce x^d, x^{d+1}, ..., x^{2d-2} using x^d ≡ -a
    for k in range(2 * d - 2, d - 1, -1):
        c = out[k]
        if c == 0:
            continue
        out[k - d] = (out[k - d] - a * c) % n
        out[k] = 0
    return out[:d]


def poly_pow_mod(base, exp, d, a, n):
    """Compute base^exp mod (x^d + a, n) by square-and-multiply."""
    result = [0] * d
    result[0] = 1
    cur = list(base)
    while exp > 0:
        if exp & 1:
            result = poly_mul_mod(result, cur, d, a, n)
        exp >>= 1
        if exp > 0:
            cur = poly_mul_mod(cur, cur, d, a, n)
    return result


def aks_congruence_holds(n, d, a, b):
    """Test (x + b)^n ≡ x^n + b (mod x^d + a, n)."""
    base = [0] * d
    base[1] = 1            # x
    base[0] = b % n        # + b
    lhs = poly_pow_mod(base, n, d, a, n)

    x_poly = [0] * d
    x_poly[1] = 1
    rhs = poly_pow_mod(x_poly, n, d, a, n)
    rhs[0] = (rhs[0] + b) % n

    return lhs == rhs


# ------------------------------------------------------------
# Capelli's irreducibility test for x^d + a over F_p
# ------------------------------------------------------------

def is_qth_power_modp(c, q, p):
    """Check whether c is a q-th power residue in F_p^*. Uses Euler-style test."""
    c = c % p
    if c == 0:
        return True
    # c is a q-th power in F_p^* iff c^((p-1)/gcd(q, p-1)) == 1
    g = math.gcd(q, p - 1)
    return pow(c, (p - 1) // g, p) == 1


def x_d_plus_a_irreducible_modp(d, a, p):
    """Capelli: x^d + a irreducible over F_p iff for every prime q | d,
    -a is NOT a q-th power in F_p; and if 4 | d, -4a is not a 4th power."""
    minus_a = (-a) % p
    for q in factorint(d).keys():
        if is_qth_power_modp(minus_a, q, p):
            return False
    if d % 4 == 0:
        if is_qth_power_modp((-4 * a) % p, 4, p):
            return False
    return True


# ------------------------------------------------------------
# Experiments
# ------------------------------------------------------------

def experiment_q1_q2(M=200, max_d=None, b_max=4):
    """For composites n in [4, M], test the AKS congruence in Z_n[x]/(x^d+a)
    for several (d, a, b) and report false-positive rates."""
    print(f"\n=== Q1/Q2: false-positive rate of non-cyclotomic AKS, n in [4, {M}] ===")
    composites = [n for n in range(4, M + 1) if not isprime(n)]
    primes = [n for n in primerange(4, M + 1)]
    print(f"primes={len(primes)}, composites={len(composites)}")

    rows = []
    for n in [49, 51, 91, 121, 169, 187, 221, 289, 341, 561]:  # spotcheck
        L = max(2, int(math.ceil(math.log2(n))))
        for d_choice in (3, 4, 5, max(3, int(math.ceil(L * L / 4)))):
            for a in (1, 2, 3):
                # If n < some prime factor, skip
                results = []
                for b in range(1, b_max + 1):
                    try:
                        ok = aks_congruence_holds(n, d_choice, a, b)
                        results.append((b, ok))
                    except Exception as e:
                        results.append((b, f"err:{e}"))
                rows.append((n, d_choice, a, results))

    print("\nSpot-check (composite n, congruence outcome per b):")
    print(f"{'n':>5} {'d':>3} {'a':>3}  results (b: pass?)")
    for n, d, a, res in rows:
        s = " ".join(f"b={b}:{'P' if v is True else ('F' if v is False else 'E')}"
                     for b, v in res)
        print(f"{n:>5} {d:>3} {a:>3}  {s}")

    # Aggregate FP rate across all composites in [4, M]
    print("\nAggregate sweep (b=1, a=1):")
    print(f"{'d':>3} {'pass_prime':>11} {'pass_comp':>10} {'FP_rate':>9}")
    for d in (2, 3, 4, 5, 6, 8):
        prime_pass = sum(1 for n in primes if aks_congruence_holds(n, d, 1, 1))
        comp_pass = sum(1 for n in composites if aks_congruence_holds(n, d, 1, 1))
        fp = comp_pass / max(1, len(composites))
        print(f"{d:>3}   {prime_pass}/{len(primes):<6}  {comp_pass}/{len(composites):<5}  {fp:>9.4f}")

    print("\nAggregate sweep (b=1..4 ALL must pass, a=1):")
    print(f"{'d':>3} {'pass_prime':>11} {'pass_comp':>10} {'FP_rate':>9}")
    for d in (2, 3, 4, 5, 6, 8):
        def all_pass(n_):
            return all(aks_congruence_holds(n_, d, 1, b) for b in (1, 2, 3, 4))
        prime_pass = sum(1 for n in primes if all_pass(n))
        comp_pass = sum(1 for n in composites if all_pass(n))
        fp = comp_pass / max(1, len(composites))
        print(f"{d:>3}   {prime_pass}/{len(primes):<6}  {comp_pass}/{len(composites):<5}  {fp:>9.4f}")


def experiment_q3_compose_a(M=300):
    """Does varying a help separate primes from composites at fixed d?"""
    print(f"\n=== Q3: variation in a (fixed small d), n in [4, {M}] ===")
    primes = list(primerange(4, M + 1))
    composites = [n for n in range(4, M + 1) if not isprime(n)]
    for d in (3, 4):
        print(f"\n d={d}, b=1, sweeping a in [1..6]")
        print(f"{'a':>3} {'prime_pass':>11} {'comp_pass':>10} {'FP_rate':>9}")
        for a in range(1, 7):
            prime_pass = sum(1 for n in primes if aks_congruence_holds(n, d, a, 1))
            comp_pass = sum(1 for n in composites if aks_congruence_holds(n, d, a, 1))
            fp = comp_pass / max(1, len(composites))
            print(f"{a:>3}   {prime_pass}/{len(primes):<6}  {comp_pass}/{len(composites):<5}  {fp:>9.4f}")


def experiment_q4_irreducibility():
    """Probe Capelli irreducibility cost for x^d + a over small primes p."""
    print("\n=== Q4: irreducibility of x^d + a over F_p (Capelli) ===")
    print("Frequency of (d, a) pairs giving irreducible x^d + a over random small primes:")
    for d in (2, 3, 4, 5, 6, 7, 8, 12, 16):
        irr_count = 0
        total = 0
        for a in range(1, 11):
            for p in primerange(7, 200):
                total += 1
                if x_d_plus_a_irreducible_modp(d, a, p):
                    irr_count += 1
        print(f" d={d:>3}  irreducible fraction = {irr_count}/{total} = {irr_count/total:.3f}")

    print("\nCircuit-depth estimate: Capelli reduces to {(p-1)/gcd(q,p-1)}-th-power")
    print("residue tests for each prime q | d. Each test is one modular exponentiation,")
    print("which is in TC^0 (Hesse-Allender-Barrington 2002). Number of tests = omega(d)")
    print("= O(log d / log log d) = O(log log n) for d = polylog(n). Total depth O(1)*O(1).")
    print("=> Capelli IS in TC^0 PROVIDED p is given. Issue: we need p | n, requiring")
    print("at least partial factorisation of n — circular.")


def experiment_q5_correctness_failure(d=4, a=1, M=500):
    """Stress-test: do composites n EVER fail the test?
    If yes, what is the smallest b that catches each composite?
    If b=1 alone suffices for all composites in [4, M], the test is interesting."""
    print(f"\n=== Q5: smallest b that catches each composite, d={d}, a={a}, M={M} ===")
    catches = {}
    survived = []
    for n in range(4, M + 1):
        if isprime(n):
            continue
        smallest = None
        for b in range(1, 12):
            if not aks_congruence_holds(n, d, a, b):
                smallest = b
                break
        if smallest is None:
            survived.append(n)
        else:
            catches.setdefault(smallest, []).append(n)
    for k in sorted(catches.keys()):
        v = catches[k]
        print(f" caught at b={k}: {len(v)} composites  e.g. {v[:6]}")
    print(f" survived all b in [1..11]: {len(survived)} composites")
    if survived:
        print(f"   (these are AKS-style false positives in this ring, sample:")
        print(f"    {survived[:20]}")
        # factor a few survivors
        for n in survived[:8]:
            print(f"     n={n}, factors={factorint(n)}")


def experiment_q6_adversarial(d_list=(3, 4, 5, 6, 8, 12)):
    """Probe Carmichael numbers and strong pseudoprimes.
    These are the canonical adversarial cases — they fool naive Fermat tests."""
    print("\n=== Q6: ADVERSARIAL test on Carmichael / pseudoprimes ===")
    # First 30 Carmichael numbers (OEIS A002997)
    carmichael = [561, 1105, 1729, 2465, 2821, 6601, 8911, 10585, 15841, 29341,
                  41041, 46657, 52633, 62745, 63973, 75361, 101101, 115921, 126217, 162401,
                  172081, 188461, 252601, 278545, 294409, 314821, 334153, 340561, 399001, 410041]
    # A few Fermat pseudoprimes base 2
    psp_base2 = [341, 561, 645, 1105, 1387, 1729, 1905, 2047, 2465, 2701, 2821, 3277, 4033]
    # Some powers of small primes (square-trap)
    prime_powers = [25, 49, 121, 169, 289, 361, 529, 841, 961, 1369, 1681]

    cases = {"Carmichael": carmichael,
             "Fermat-PSP (base 2)": psp_base2,
             "Prime powers p^2": prime_powers}

    for d in d_list:
        for a in (1, 2):
            print(f"\n--- d={d}, a={a} ---")
            for cname, nums in cases.items():
                passes_b1 = []
                passes_all_b14 = []
                for n in nums:
                    if aks_congruence_holds(n, d, a, 1):
                        passes_b1.append(n)
                    if all(aks_congruence_holds(n, d, a, b) for b in (1, 2, 3, 4)):
                        passes_all_b14.append(n)
                print(f"  {cname:24s} (n={len(nums)}):  pass(b=1) = {len(passes_b1)},"
                      f"  pass(b=1..4) = {len(passes_all_b14)}")
                if passes_b1:
                    print(f"      passes_b1 sample: {passes_b1[:6]}")
                if passes_all_b14:
                    print(f"      passes_b14 sample: {passes_all_b14[:6]}")


def experiment_q7_scale(d_fixed=4, a=1, M_cap=5000, sample=300):
    """How does FP rate scale up to n = M_cap?"""
    print(f"\n=== Q7: scaling FP rate to large n (d={d_fixed}, a={a}, b=1) ===")
    composites = [n for n in range(4, M_cap + 1) if not isprime(n)]
    if len(composites) > sample:
        # Take a stratified sample (every k-th)
        step = len(composites) // sample
        composites = composites[::step][:sample]
    primes = list(primerange(4, M_cap + 1))
    if len(primes) > sample:
        step = len(primes) // sample
        primes = primes[::step][:sample]

    cp = sum(1 for n in primes if aks_congruence_holds(n, d_fixed, a, 1))
    cc = sum(1 for n in composites if aks_congruence_holds(n, d_fixed, a, 1))
    print(f"  primes pass: {cp}/{len(primes)}")
    print(f"  composites pass (FP): {cc}/{len(composites)} = {cc/len(composites):.4f}")
    if cc > 0:
        survivors = [n for n in composites if aks_congruence_holds(n, d_fixed, a, 1)]
        print(f"  survivor sample: {survivors[:20]}")
        for n in survivors[:8]:
            print(f"     n={n}, factors={factorint(n)}")


if __name__ == "__main__":
    print("Non-cyclotomic ring AKS variant: empirical probe")
    print("=" * 60)
    experiment_q1_q2(M=200, b_max=4)
    experiment_q3_compose_a(M=300)
    experiment_q4_irreducibility()
    experiment_q5_correctness_failure(d=3, a=1, M=500)
    experiment_q5_correctness_failure(d=4, a=1, M=500)
    experiment_q5_correctness_failure(d=4, a=2, M=500)
    experiment_q5_correctness_failure(d=6, a=1, M=500)
    experiment_q6_adversarial()
    experiment_q7_scale(d_fixed=4, a=1, M_cap=5000, sample=300)
    experiment_q7_scale(d_fixed=4, a=2, M_cap=5000, sample=300)
    experiment_q7_scale(d_fixed=6, a=1, M_cap=5000, sample=300)
    print("\nDone.")
