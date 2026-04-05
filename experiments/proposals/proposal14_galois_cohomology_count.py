"""
PROPOSAL 6: Galois Cohomology Point Counting Transfer
======================================================

IDEA: In algebraic geometry, the Weil conjectures (proven by Deligne) give
EXACT formulas for counting F_q-rational points on varieties:

  |V(F_q)| = sum_{i} (-1)^i Tr(Frob | H^i_et(V, Q_l))

For the "variety" whose F_p-points are primes, this is trivial (primes ARE
the structure). But there's a TRANSFER principle:

If we can encode pi(x) as a point count on a suitable algebraic variety,
then the Lefschetz trace formula gives it as a sum over eigenvalues of
Frobenius, which might be computable.

CONCRETE APPROACH: The Chebyshev psi function
  psi(x) = sum_{p^k <= x} ln(p)
is related to the trace of a "prime-counting operator."

The EXPLICIT FORMULA is:
  psi(x) = x - sum_rho x^rho / rho - ln(2*pi) - (1/2)*ln(1 - x^{-2})

This IS a trace formula! The analogy:
  - psi(x) <-> |V(F_q)| (point count)
  - x <-> q (the "base field")
  - rho <-> eigenvalues of Frobenius
  - The "variety" is the "arithmetic site" (Connes-Consani)

CAN WE EXPLOIT THIS?

In algebraic geometry, there exist fast algorithms for computing
|V(F_q)| for specific varieties:
  - Elliptic curves: O(log(q)^4) via Schoof's algorithm
  - Hyperelliptic curves: O(g^3 * log(q)^3) via Kedlaya's algorithm

These work because H^1 is FINITE-DIMENSIONAL (dimension = 2g).

For pi(x), the "H^1" is INFINITE-dimensional (all zeta zeros).
But what if we could TRUNCATE or APPROXIMATE it?

SCHOOF-TYPE APPROACH FOR PI(X):
1. Compute pi(x) mod l for small primes l, using the structure of
   primes in arithmetic progressions
2. For each l: pi(x) mod l relates to the number of primes ≡ a (mod l)
3. Use Dirichlet L-functions to compute pi(x; l, a) mod l
4. CRT reconstruction

THE KEY QUESTION: Can we compute pi(x) mod l in O(polylog(x)) time,
analogous to how Schoof computes |E(F_p)| mod l?

SCHOOF'S KEY TRICK: For elliptic curves, |E(F_p)| mod l is determined
by the action of Frobenius on E[l] (the l-torsion points), which lives
in F_{p^{O(l^2)}}. The computation is polynomial in log(p) and l.

ANALOGUE: For primes, "Frobenius on the l-torsion" would be the action
of the Galois group Gal(Q(zeta_l)/Q) on primes. This IS Dirichlet's
theorem! The character chi mod l plays the role of the l-torsion
representation.

So: pi(x) mod l = f(characters mod l, L-function values)

The Dirichlet L-function at s=1 is computable in O(polylog(l)) time.
But we need L-function VALUES related to counting primes up to x,
which requires the zeros of L(s, chi) up to height ~x^{1/2}.

UNLESS: We only need the TRACE (sum over all characters), which might
telescope or simplify.

TEST: Compute sum_{chi mod q} chi(a) * pi(x, q, a) for various q, x
and check for simplification.
"""

import numpy as np
from sympy import (prime, primepi, isprime, nextprime, primerange,
                   totient, factorint, primitive_root, Mod)
from mpmath import mp, mpf, li, log, exp, pi as mpi, cos, sin, zeta
import math
from collections import Counter

def dirichlet_characters(q):
    """
    Generate all Dirichlet characters mod q.
    Returns list of dicts mapping residues -> character values.
    """
    if q == 1:
        return [{0: 1}]

    phi_q = int(totient(q))

    try:
        g = int(primitive_root(q))
    except:
        # q is not prime, use multiplicative group structure
        # Simplified: only handle prime q
        return None

    # Characters are determined by chi(g) = e^{2*pi*i*k/phi(q)}
    characters = []
    for k in range(phi_q):
        omega = np.exp(2j * np.pi * k / phi_q)
        char = {}
        power = 1
        for j in range(phi_q):
            char[power] = omega ** (k * j)
            power = (power * g) % q
        # Set chi(a) = 0 for gcd(a,q) > 1
        for a in range(q):
            if a not in char:
                char[a] = 0
        characters.append(char)

    return characters

def pi_in_class(x, q, a):
    """Count primes p <= x with p ≡ a (mod q)."""
    count = 0
    for p in primerange(2, x + 1):
        if int(p) % q == a:
            count += 1
    return count

def character_sum_test(x, q):
    """
    Compute S(chi) = sum_{p <= x} chi(p) for each character chi mod q.

    By the explicit formula for Dirichlet L-functions:
    sum_{p <= x} chi(p) * ln(p) = delta_{chi=chi_0} * x - sum_rho x^rho/rho + ...

    The principal character gives ~ x, all others are small (GRH: O(sqrt(x)*log(x))).

    KEY: If we can compute S(chi) for all chi mod q, we can recover
    pi(x; q, a) = (1/phi(q)) * sum_chi chi_bar(a) * S(chi)
    """
    chars = dirichlet_characters(q)
    if chars is None:
        return None

    phi_q = int(totient(q))

    # Compute S(chi) for each character
    char_sums = []
    for chi in chars:
        s = 0
        for p in primerange(2, x + 1):
            s += chi[int(p) % q]
        char_sums.append(s)

    # Recover pi(x; q, a) from character sums
    recovered = {}
    for a in range(q):
        if math.gcd(a, q) == 1:
            val = 0
            for k, chi in enumerate(chars):
                chi_bar_a = np.conj(chi.get(a, 0))
                val += chi_bar_a * char_sums[k]
            val = val / phi_q
            recovered[a] = round(val.real)

    # Direct computation for verification
    direct = {}
    for a in range(q):
        if math.gcd(a, q) == 1:
            direct[a] = pi_in_class(x, q, a)

    return {
        'q': q,
        'x': x,
        'char_sums': char_sums,
        'recovered': recovered,
        'direct': direct,
        'match': all(recovered.get(a) == direct.get(a) for a in direct),
        'principal_sum': char_sums[0],
        'max_non_principal': max(abs(s) for s in char_sums[1:]) if len(char_sums) > 1 else 0,
    }

def schoof_analogue_test(x, max_q=20):
    """
    Test the Schoof analogue:
    For each prime q, compute pi(x) mod q using the character-theoretic formula.

    pi(x) = sum_a pi(x; q, a) for a in (Z/qZ)*
    plus pi(x; q, 0) for primes dividing q.

    We need pi(x) mod q. Can the character sums help?

    pi(x) mod q = [sum_a pi(x; q, a)] mod q + [count of p|q, p<=x] mod q
                = [(1/phi(q)) * sum_chi sum_{a coprime to q} chi_bar(a) * S(chi)] mod q + small
                = [(1/phi(q)) * S(chi_0) * phi(q)] mod q + correction
                = [S(chi_0)] mod q + correction from non-principal characters

    S(chi_0) = number of primes <= x coprime to q = pi(x) - #{p <= x, p | q}
    So pi(x) mod q = S(chi_0) mod q + #{p | q, p <= x} mod q

    This is CIRCULAR: S(chi_0) = pi(x) - small correction.

    BUT: the non-principal character sums give the DISTRIBUTION across classes.
    If we know the distribution, we know each class count, and their sum is pi(x).

    The non-principal sums are bounded by O(sqrt(x)*log(x)) (GRH).
    So S(chi) for chi != chi_0 has at most ~log(x)/2 bits of information.
    These bits encode the distribution, and by summing we recover pi(x).

    This is just another way of saying: pi(x) = smooth + oscillatory.
    """
    pi_x = int(primepi(x))
    results = {}

    for q in primerange(2, max_q + 1):
        q = int(q)
        test = character_sum_test(x, q)
        if test is None:
            continue

        # pi(x) mod q
        pi_mod_q = pi_x % q

        # Can we get this from character sums alone?
        # S(chi_0) = pi(x) - #{p|q, p<=x}
        primes_dividing_q = sum(1 for p in range(2, q+1) if isprime(p) and q % p == 0 and p <= x)
        s_chi0 = pi_x - primes_dividing_q

        results[q] = {
            'pi_mod_q': pi_mod_q,
            's_chi0': int(test['char_sums'][0].real),
            'max_non_principal': float(test['max_non_principal']),
            'sqrt_x': math.sqrt(x),
            'ratio': float(test['max_non_principal']) / math.sqrt(x) if x > 0 else 0,
            'match': test['match'],
        }

    return results

def trace_formula_connection(x_values):
    """
    Explore the connection between the Selberg trace formula and prime counting.

    The Selberg trace formula for SL(2,Z)\H gives:
    sum_n h(r_n) = (Area/4pi) * integral h(r) r*tanh(pi*r) dr
                   + sum_{T} (log N(T0)) / (N(T)^{1/2} - N(T)^{-1/2}) * g(log N(T))
                   + ...

    where r_n are spectral parameters (related to Maass form eigenvalues)
    and T ranges over conjugacy classes (related to primes via N(T) = p^k).

    The ANALOGY: spectral parameters <-> zeta zeros, conjugacy classes <-> primes.

    In the Selberg trace formula, the SUM OVER PRIMES appears on one side.
    If we could evaluate the spectral side efficiently, we'd get the prime side.

    The spectral side involves Maass form eigenvalues, which are discrete
    but known to high precision for small eigenvalues.
    """
    results = []
    for x in x_values:
        pi_x = int(primepi(x))
        # The "test function" h(r) = indicator of |r| <= T
        # gives a count related to pi(x)

        # Spectral side: sum over Laplacian eigenvalues lambda_n = 1/4 + r_n^2
        # For SL(2,Z), first few eigenvalues are known
        # lambda_1 ~ 91.14... (Maass form)
        # r_1 ~ 9.534...

        # The geometric side relates to class numbers and prime geodesics
        # Prime geodesics of length ln(p) for primes p

        # Count: sum_{p^k, p prime} ln(p) * f(ln(p^k))
        # = psi(x) for appropriate f

        results.append({
            'x': x,
            'pi_x': pi_x,
            'psi_x': sum(math.log(p) * (math.floor(math.log(x) / math.log(p)))
                        for p in primerange(2, x + 1)),
            'x_minus_psi': x - sum(math.log(p) * (math.floor(math.log(x) / math.log(p)))
                                   for p in primerange(2, x + 1)),
        })

    return results


if __name__ == '__main__':
    print("=" * 80)
    print("PROPOSAL 6: Galois Cohomology / Schoof Analogue")
    print("=" * 80)

    print("\n--- Part A: Character sum analysis ---")
    for q in [3, 5, 7, 11]:
        for x in [100, 1000]:
            result = character_sum_test(x, q)
            if result:
                print(f"  q={q}, x={x}: match={result['match']}, "
                      f"|S_principal|={abs(result['principal_sum']):.0f}, "
                      f"|S_max_nonprinc|={result['max_non_principal']:.1f}, "
                      f"sqrt(x)={math.sqrt(x):.1f}")

    print("\n--- Part B: Schoof analogue - pi(x) mod q from characters ---")
    for x in [1000, 5000]:
        results = schoof_analogue_test(x, max_q=20)
        print(f"\n  x={x}, pi(x)={int(primepi(x))}:")
        for q, info in sorted(results.items()):
            print(f"    q={q:2d}: pi mod q = {info['pi_mod_q']}, "
                  f"|non-principal|/sqrt(x) = {info['ratio']:.3f}, "
                  f"match={info['match']}")

    print("\n--- Part C: Trace formula connection ---")
    trace_results = trace_formula_connection([100, 1000, 10000])
    for r in trace_results:
        print(f"  x={r['x']:6d}: pi(x)={r['pi_x']:5d}, "
              f"psi(x)={r['psi_x']:.2f}, x-psi(x)={r['x_minus_psi']:.2f}")

    print("\n--- Part D: Non-principal character sum scaling ---")
    print("  GRH predicts |S(chi)| = O(sqrt(x) * log(x))")
    for x in [100, 500, 1000, 5000, 10000]:
        for q in [3, 5, 7]:
            result = character_sum_test(x, q)
            if result:
                ratio = result['max_non_principal'] / (math.sqrt(x) * math.log(x))
                print(f"    x={x:6d}, q={q}: |S_max|={result['max_non_principal']:.1f}, "
                      f"sqrt(x)*ln(x)={math.sqrt(x)*math.log(x):.1f}, "
                      f"ratio={ratio:.4f}")
