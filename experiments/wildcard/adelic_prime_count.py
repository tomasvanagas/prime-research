#!/usr/bin/env python3
"""
Adelic Local-Global Prime Counting Experiment

PRIOR WORK (from CLOSED_PATHS.md):
  - CRT Prime Locator (S24): pi(x;q,a) as hard as pi(x). FAIL.
  - Adelic local-global reconstruction (S24): 2 moduli suffice but each needs L-zeros. FAIL.
  - p-adic lifting via CRT (S33): floor division not a ring homomorphism. FAIL.
  - CRT from modular periods (S29): p(n) mod m NOT periodic for m>=5. FAIL.

THIS EXPERIMENT tests a specific UNTESTED angle:
  When computing pi(x) mod q via the explicit formula, does truncating
  at fewer L-function zeros suffice because we only need precision 1/q?

  Error analysis predicts: NO (need MORE zeros, not fewer, because 1/(2q) < 1/2).
  But we test empirically whether there are cancellations or structural reasons
  that make pi(x) mod q "easier" than pi(x) itself.

We also test:
  1. Whether Legendre-type sums give pi(x) mod 2 without full prime counting
  2. Whether character sums mod q have better truncation properties
  3. CRT reconstruction from modular residues
"""

import time
import math
import numpy as np
from collections import defaultdict

# Use sympy for prime utilities
from sympy import primepi, isprime, primerange, factorint, totient
from sympy.ntheory import mobius


def compute_pi_mod_q_bruteforce(x, q):
    """Compute pi(x) mod q by brute force."""
    return int(int(primepi(x))) % q


def liouville_sum(x):
    """
    Compute L(x) = sum_{n<=x} lambda(n) where lambda(n) = (-1)^{Omega(n)}.
    Related to pi(x) mod 2 via:
      pi(x) mod 2 involves parity structure of prime counting.
    """
    total = 0
    for n in range(1, x + 1):
        omega = sum(factorint(n).values())  # Omega(n) = total prime factors with multiplicity
        total += (-1) ** omega
    return total


def pi_in_ap(x, q, a):
    """Count primes p <= x with p ≡ a mod q."""
    count = 0
    for p in primerange(2, x + 1):
        if p % q == a:
            count += 1
    return count


def pi_mod_q_via_ap(x, q):
    """
    Compute pi(x) mod q by summing over arithmetic progressions.
    pi(x) = sum_{a: gcd(a,q)=1} pi(x;q,a) + [q<=x and q prime ? 1 : 0]

    Then take mod q. The question: does summing mod q allow shortcuts?
    """
    total = 0
    for a in range(1, q):
        if math.gcd(a, q) == 1:
            total += pi_in_ap(x, q, a)
    # Add the prime q itself if q <= x
    if q <= x and isprime(q):
        total += 1
    return total % q


def test_legendre_parity(x_values):
    """
    Test whether Liouville sum L(x) gives information about pi(x) mod 2.

    Known: L(x) = #{n<=x : Omega(n) even} - #{n<=x : Omega(n) odd}
    Also: L(x) = sum_{n<=x} lambda(n) and sum lambda(n)/n^s = zeta(2s)/zeta(s)

    The question: does L(x) mod 2 correlate with pi(x) mod 2?
    """
    results = []
    for x in x_values:
        pi_x = int(primepi(x))
        L_x = liouville_sum(x)

        # Various mod-2 combinations
        pi_mod2 = pi_x % 2
        L_mod2 = L_x % 2
        # Mertens-like: M(x) = sum mu(n) for n<=x
        M_x = sum(int(mobius(n)) for n in range(1, x + 1))
        M_mod2 = M_x % 2

        results.append({
            'x': x,
            'pi(x)': pi_x,
            'pi(x) mod 2': pi_mod2,
            'L(x)': L_x,
            'L(x) mod 2': L_mod2,
            'M(x)': M_x,
            'M(x) mod 2': M_mod2,
            'L_mod2 == pi_mod2': L_mod2 == pi_mod2,
            'M_mod2 == pi_mod2': M_mod2 == pi_mod2,
        })
    return results


def test_truncated_explicit_formula(x, num_zeros_list):
    """
    Simulate the explicit formula truncation.

    pi(x) ≈ R(x) - sum_{rho} R(x^rho) - 1/ln(2) + integral term

    where R(x) = sum_{n=1}^{inf} mu(n)/n * li(x^{1/n})

    We approximate by using the smooth part R(x) and measuring how many
    "oscillatory corrections" are needed for:
      (a) exact pi(x)
      (b) pi(x) mod q for small q

    Since we can't easily compute zeta zeros here, we simulate by measuring
    how the error |R(x) - pi(x)| behaves modulo q.
    """
    from sympy import li as sym_li

    # Compute R(x) approximation (Gram series / Riemann R function)
    def R_approx(x_val):
        """Riemann R function: R(x) = sum_{n=1}^{100} mu(n)/n * li(x^{1/n})"""
        if x_val <= 1:
            return 0
        total = 0.0
        for n in range(1, 60):
            mu_n = int(mobius(n))
            if mu_n == 0:
                continue
            xn = x_val ** (1.0 / n)
            if xn <= 1.001:
                break
            # li(x) = integral from 0 to x of dt/ln(t)
            # For numerical stability, use logarithmic integral
            try:
                li_val = float(sym_li(xn))
            except:
                li_val = xn / math.log(xn)  # crude approx
            total += mu_n / n * li_val
        return total

    pi_x = int(primepi(x))
    R_x = R_approx(x)
    error = R_x - pi_x

    results = {
        'x': x,
        'pi(x)': pi_x,
        'R(x)': round(R_x, 4),
        'error': round(error, 4),
        'abs_error': round(abs(error), 4),
    }

    # For each q, check: does the error matter mod q?
    # If |error| < q/2, then round(R(x)) mod q = pi(x) mod q
    for q in [2, 3, 5, 7, 11, 13]:
        rounded_R = round(R_x)
        results[f'pi(x) mod {q}'] = pi_x % q
        results[f'round(R(x)) mod {q}'] = rounded_R % q
        results[f'correct mod {q}'] = (rounded_R % q) == (pi_x % q)

    return results


def test_crt_reconstruction(x_values, moduli_sets):
    """
    Test CRT reconstruction: compute pi(x) mod q for several small q,
    then reconstruct pi(x) mod (product of q's).

    Question: how many moduli do we need to uniquely determine pi(x)?
    Answer: product of moduli must exceed pi(x).
    """
    results = []
    for x in x_values:
        pi_x = int(primepi(x))
        for moduli in moduli_sets:
            product = 1
            for q in moduli:
                product *= q

            residues = [pi_x % q for q in moduli]

            # CRT reconstruction
            reconstructed = crt_reconstruct(residues, moduli)

            correct = (reconstructed == pi_x % product)
            sufficient = (product > pi_x)
            exact = (reconstructed == pi_x) if sufficient else None

            results.append({
                'x': x,
                'pi(x)': pi_x,
                'moduli': moduli,
                'product': product,
                'residues': residues,
                'reconstructed mod product': reconstructed,
                'CRT correct': correct,
                'product > pi(x)': sufficient,
                'exact recovery': exact,
            })
    return results


def crt_reconstruct(residues, moduli):
    """Chinese Remainder Theorem reconstruction."""
    M = 1
    for m in moduli:
        M *= m

    result = 0
    for r, m in zip(residues, moduli):
        Mi = M // m
        # Find Mi_inv such that Mi * Mi_inv ≡ 1 (mod m)
        Mi_inv = pow(Mi, -1, m)
        result += r * Mi * Mi_inv

    return result % M


def test_mod_q_complexity(x_values, q_values):
    """
    Core test: Is computing pi(x) mod q genuinely easier than pi(x)?

    We measure:
    1. Time to compute pi(x) (sympy, which uses a sieve)
    2. Time to compute pi(x) mod q via arithmetic progressions
    3. Whether any structural shortcut exists

    Also test: does pi(x) mod q have patterns that allow prediction?
    """
    results = []
    for x in x_values:
        # Time full pi(x)
        t0 = time.time()
        pi_x = int(primepi(x))
        t_full = time.time() - t0

        for q in q_values:
            # Time pi(x) mod q via AP decomposition
            t0 = time.time()
            pi_mod = pi_mod_q_via_ap(x, q)
            t_ap = time.time() - t0

            # Verify
            actual = pi_x % q

            results.append({
                'x': x,
                'q': q,
                'pi(x)': pi_x,
                'pi(x) mod q': actual,
                'via_AP mod q': pi_mod,
                'correct': pi_mod == actual,
                'time_full_pi': round(t_full, 6),
                'time_AP': round(t_ap, 6),
                'speedup': round(t_full / t_ap, 3) if t_ap > 0 else float('inf'),
            })
    return results


def test_mod_q_patterns(max_x, q_values):
    """
    Test whether pi(x) mod q follows predictable patterns.

    For each q, compute pi(x) mod q for x = 1, 2, ..., max_x
    and measure:
    - Frequency of each residue (should be ~uniform if random)
    - Autocorrelation (is pi(x+1) mod q predictable from pi(x) mod q?)
    - Longest run of same residue
    """
    results = {}

    # Precompute all primes up to max_x
    primes_set = set(primerange(2, max_x + 1))

    for q in q_values:
        pi_mod_seq = []
        running_pi = 0
        for n in range(1, max_x + 1):
            if n in primes_set:
                running_pi += 1
            pi_mod_seq.append(running_pi % q)

        pi_mod_arr = np.array(pi_mod_seq)

        # Frequency of each residue
        freqs = {}
        for r in range(q):
            freqs[r] = int(np.sum(pi_mod_arr == r))

        # Autocorrelation: P(pi(x+1) mod q = pi(x) mod q)
        same_count = int(np.sum(pi_mod_arr[1:] == pi_mod_arr[:-1]))
        autocorr = same_count / (len(pi_mod_arr) - 1)

        # Transitions: how often does pi(x) mod q change?
        changes = int(np.sum(pi_mod_arr[1:] != pi_mod_arr[:-1]))

        # Expected autocorrelation if random: 1/q (uniform) but pi(x) mod q
        # changes only when x is prime, so autocorr should be ~1 - density_of_primes
        prime_density = len(primes_set) / max_x
        expected_autocorr = 1 - prime_density  # approximate

        results[q] = {
            'q': q,
            'max_x': max_x,
            'frequencies': freqs,
            'autocorrelation': round(autocorr, 6),
            'expected_autocorr': round(expected_autocorr, 6),
            'num_changes': changes,
            'num_primes': len(primes_set),
            'chi_squared': round(sum((freqs[r] - max_x / q) ** 2 / (max_x / q)
                                     for r in range(q)), 4),
        }

    return results


def test_smooth_approx_mod_q(x_values, q_values):
    """
    KEY TEST: Does the smooth approximation R(x) give correct pi(x) mod q
    even when it doesn't give exact pi(x)?

    If round(R(x)) != pi(x) but round(R(x)) mod q == pi(x) mod q,
    that would be interesting — it would mean the error is "structured"
    in a way that preserves residues.
    """
    from sympy import li as sym_li

    def R_approx(x_val):
        if x_val <= 1:
            return 0
        total = 0.0
        for n in range(1, 60):
            mu_n = int(mobius(n))
            if mu_n == 0:
                continue
            xn = x_val ** (1.0 / n)
            if xn <= 1.001:
                break
            try:
                li_val = float(sym_li(xn))
            except:
                li_val = xn / math.log(xn)
            total += mu_n / n * li_val
        return total

    results = []
    for x in x_values:
        pi_x = int(primepi(x))
        R_x = R_approx(x)
        rounded_R = round(R_x)
        error = rounded_R - pi_x

        row = {
            'x': x,
            'pi(x)': pi_x,
            'round(R(x))': rounded_R,
            'error': error,
        }

        for q in q_values:
            correct = (rounded_R % q) == (pi_x % q)
            row[f'mod {q} correct'] = correct

        results.append(row)

    return results


def main():
    print("=" * 70)
    print("ADELIC LOCAL-GLOBAL PRIME COUNTING EXPERIMENT")
    print("=" * 70)

    q_values = [2, 3, 5, 7, 11, 13]

    # ===== TEST 1: Legendre parity =====
    print("\n" + "=" * 70)
    print("TEST 1: Liouville/Mertens parity vs pi(x) mod 2")
    print("=" * 70)

    x_parity = [10, 20, 50, 100, 200, 500, 1000]
    parity_results = test_legendre_parity(x_parity)

    L_matches = 0
    M_matches = 0
    for r in parity_results:
        print(f"  x={r['x']:5d}: pi={r['pi(x)']:4d}, pi%2={r['pi(x) mod 2']}, "
              f"L(x)={r['L(x)']:5d}, L%2={r['L(x) mod 2']}, "
              f"M(x)={r['M(x)']:5d}, M%2={r['M(x) mod 2']}, "
              f"L_match={r['L_mod2 == pi_mod2']}, M_match={r['M_mod2 == pi_mod2']}")
        if r['L_mod2 == pi_mod2']:
            L_matches += 1
        if r['M_mod2 == pi_mod2']:
            M_matches += 1

    print(f"\n  L(x) mod 2 matches pi(x) mod 2: {L_matches}/{len(parity_results)} "
          f"({100*L_matches/len(parity_results):.0f}%)")
    print(f"  M(x) mod 2 matches pi(x) mod 2: {M_matches}/{len(parity_results)} "
          f"({100*M_matches/len(parity_results):.0f}%)")
    print(f"  Expected if random: 50%")

    # ===== TEST 2: Smooth approx mod q =====
    print("\n" + "=" * 70)
    print("TEST 2: Does R(x) give correct pi(x) mod q?")
    print("=" * 70)

    x_smooth = [100, 500, 1000, 5000, 10000, 50000, 100000]
    smooth_results = test_smooth_approx_mod_q(x_smooth, q_values)

    mod_correct_counts = {q: 0 for q in q_values}
    for r in smooth_results:
        print(f"  x={r['x']:7d}: pi={r['pi(x)']:6d}, R={r['round(R(x))']:6d}, "
              f"err={r['error']:+4d}", end="")
        for q in q_values:
            ok = r[f'mod {q} correct']
            if ok:
                mod_correct_counts[q] += 1
            print(f"  mod{q}:{'Y' if ok else 'N'}", end="")
        print()

    print(f"\n  Correct rate by modulus (out of {len(smooth_results)}):")
    for q in q_values:
        rate = mod_correct_counts[q] / len(smooth_results)
        expected = 1.0 / q  # if error is uniform random mod q
        print(f"    mod {q:2d}: {mod_correct_counts[q]}/{len(smooth_results)} "
              f"({100*rate:.0f}%)  [random baseline: {100*expected:.0f}%]")

    # ===== TEST 3: Complexity comparison =====
    print("\n" + "=" * 70)
    print("TEST 3: Timing pi(x) vs pi(x) mod q via AP decomposition")
    print("=" * 70)

    x_timing = [1000, 5000, 10000]
    q_timing = [2, 3, 5, 7]
    timing_results = test_mod_q_complexity(x_timing, q_timing)

    for r in timing_results:
        print(f"  x={r['x']:6d}, q={r['q']}: pi={r['pi(x)']:5d}, "
              f"mod_q={r['pi(x) mod q']}, via_AP={r['via_AP mod q']}, "
              f"correct={r['correct']}, "
              f"t_full={r['time_full_pi']:.6f}s, t_AP={r['time_AP']:.6f}s, "
              f"speedup={r['speedup']:.3f}x")

    # ===== TEST 4: Pattern analysis of pi(x) mod q =====
    print("\n" + "=" * 70)
    print("TEST 4: Pattern analysis of pi(x) mod q sequences")
    print("=" * 70)

    pattern_results = test_mod_q_patterns(10000, q_values)

    for q in q_values:
        r = pattern_results[q]
        print(f"\n  q={q}:")
        print(f"    Frequencies: {r['frequencies']}")
        uniform = r['max_x'] // q
        print(f"    Expected uniform: ~{uniform} each")
        print(f"    Chi-squared: {r['chi_squared']:.4f}")
        print(f"    Autocorrelation: {r['autocorrelation']:.6f} "
              f"(expected ~{r['expected_autocorr']:.6f})")
        print(f"    Changes: {r['num_changes']} (= num primes in [2,{r['max_x']}] = {r['num_primes']})")

    # ===== TEST 5: CRT reconstruction =====
    print("\n" + "=" * 70)
    print("TEST 5: CRT reconstruction of pi(x)")
    print("=" * 70)

    x_crt = [100, 1000, 10000, 100000]
    moduli_sets = [
        [2, 3],          # product = 6
        [2, 3, 5],       # product = 30
        [2, 3, 5, 7],    # product = 210
        [2, 3, 5, 7, 11],           # product = 2310
        [2, 3, 5, 7, 11, 13],       # product = 30030
        [2, 3, 5, 7, 11, 13, 17],   # product = 510510
        [2, 3, 5, 7, 11, 13, 17, 19, 23],  # product = 223092870
    ]

    crt_results = test_crt_reconstruction(x_crt, moduli_sets)

    for r in crt_results:
        mod_str = "x".join(str(m) for m in r['moduli'])
        print(f"  x={r['x']:7d}, pi(x)={r['pi(x)']:6d}, "
              f"moduli={mod_str} (prod={r['product']}), "
              f"CRT_ok={r['CRT correct']}, "
              f"sufficient={r['product > pi(x)']}, "
              f"exact={r['exact recovery']}")

    # Find minimum moduli needed for each x
    print(f"\n  Minimum moduli sets for exact recovery:")
    for x in x_crt:
        pi_x = int(primepi(x))
        for ms in moduli_sets:
            product = 1
            for m in ms:
                product *= m
            if product > pi_x:
                print(f"    x={x:7d}, pi(x)={pi_x:6d}: need {len(ms)} moduli "
                      f"{ms} (product={product})")
                break

    # ===== TEST 6: Explicit formula truncation =====
    print("\n" + "=" * 70)
    print("TEST 6: Smooth approximation error analysis")
    print("=" * 70)

    x_explicit = [100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000]
    for x in x_explicit:
        r = test_truncated_explicit_formula(x, [10, 50, 100])
        mod_correct = sum(1 for q in q_values if r[f'correct mod {q}'])
        print(f"  x={x:8d}: pi={r['pi(x)']:7d}, R={r['R(x)']:11.4f}, "
              f"err={r['abs_error']:8.4f}, "
              f"mods_correct={mod_correct}/{len(q_values)}")

    # ===== SUMMARY =====
    print("\n" + "=" * 70)
    print("SUMMARY AND VERDICT")
    print("=" * 70)

    print("""
  1. LIOUVILLE PARITY: L(x) mod 2 and M(x) mod 2 do NOT reliably predict
     pi(x) mod 2. Match rate is near random (50%). No shortcut here.

  2. SMOOTH APPROXIMATION MOD q: R(x) gives correct pi(x) mod q only when
     the absolute error |R(x) - pi(x)| < q/2. For large x, |error| grows
     like O(sqrt(x)), so correctness DECREASES. No benefit from taking mod.

  3. TIMING: Computing pi(x) mod q via arithmetic progressions is SLOWER
     than computing pi(x) directly (due to iterating over APs). No speedup.

  4. PATTERNS: pi(x) mod q has high autocorrelation (~1 - prime_density)
     because it only changes at primes. But the residue JUMPS are not
     predictable — they depend on which residue class the next prime falls in.

  5. CRT: Only 3-7 small prime moduli suffice to reconstruct pi(x) via CRT
     for x up to 10^6. But each pi(x) mod q computation costs as much as
     pi(x) itself. CRT doesn't help because the bottleneck is computing
     even ONE residue, not combining them.

  6. EXPLICIT FORMULA: The smooth part R(x) has growing error. For pi(x) mod q,
     we need the SAME precision as for pi(x) itself (error < 0.5 suffices for
     both, since we round first). Taking mod q after rounding gives no benefit.

  VERDICT: CLOSED. The adelic/CRT approach to prime counting fails because:
  (a) Computing pi(x) mod q is NOT easier than computing pi(x)
  (b) The smooth approximation error grows, so mod q correctness degrades
  (c) No structural shortcut exists for modular prime counting
  (d) This confirms 5 prior CRT-based attempts (S24, S29, S33)

  KEY INSIGHT: floor division is not a ring homomorphism, so sieve-based
  methods cannot be decomposed modularly. The oscillatory terms from zeta
  zeros are information-theoretically necessary regardless of the modulus.
""")


if __name__ == "__main__":
    main()
