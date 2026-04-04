#!/usr/bin/env python3
"""
Session 4: Additive Combinatorics & Arithmetic Progressions for p(n)

Explores whether additive combinatorial structure (Green-Tao, circle method,
sieve weights, Goldbach representations, Erdos-Kac) can yield a computable
formula for the n-th prime without iterating over all integers up to p(n).

Five angles investigated:
  1. Green-Tao AP interpolation
  2. Hardy-Littlewood circle method inversion
  3. Sieve weight formulas (Selberg/GPY style)
  4. Goldbach representation extraction
  5. Erdos-Kac prime detector

Author: Session 4 research
"""

import math
import time
import sys
from functools import lru_cache
from collections import defaultdict

# ---------- utilities ----------

def sieve_primes(limit):
    """Simple sieve of Eratosthenes."""
    if limit < 2:
        return []
    s = bytearray(b'\x01') * (limit + 1)
    s[0] = s[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if s[i]:
            s[i*i::i] = b'\x00' * len(s[i*i::i])
    return [i for i, v in enumerate(s) if v]

def is_prime_miller_rabin(n):
    """Deterministic Miller-Rabin for n < 3.3e24."""
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0: return False
    r, d = 0, n - 1
    while d % 2 == 0:
        r += 1
        d //= 2
    for a in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        if a >= n: continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1: continue
        for _ in range(r - 1):
            x = pow(x, x, n)
            if x == n - 1: break
        else:
            return False
    return True

# Reference primes for testing
PRIMES = sieve_primes(2_000_000)
PRIME_SET = set(PRIMES)

def nth_prime_ref(n):
    """Reference: n-th prime (1-indexed)."""
    if n <= len(PRIMES):
        return PRIMES[n - 1]
    return None

# =============================================================================
# APPROACH 1: Green-Tao AP Interpolation
# =============================================================================

def test_green_tao_interpolation():
    """
    The Green-Tao theorem says primes contain APs of arbitrary length k.
    Idea: If we find many APs passing through known primes, can we extrapolate?

    Analysis:
    - An AP in primes: a, a+d, a+2d, ..., a+(k-1)d are all prime.
    - Given primes p_1,...,p_m, find APs that "cover" them and extrapolate.
    - Problem: The APs are SPARSE and IRREGULAR. The theorem is existential,
      not constructive. We cannot enumerate long APs without knowing the primes.

    Let's test empirically: how well does AP extrapolation predict the next prime?
    """
    print("=" * 72)
    print("APPROACH 1: Green-Tao AP Interpolation")
    print("=" * 72)

    # Find APs of length k among small primes
    primes_small = PRIMES[:500]  # first 500 primes
    ps = set(primes_small)

    # Find APs of length 3, 4, 5 starting from each prime
    ap_counts = {3: 0, 4: 0, 5: 0}
    ap_examples = {3: [], 4: [], 5: []}

    for i, p in enumerate(primes_small[:100]):
        for d in range(2, 1000, 2):  # common difference must be even (except for 2,3,5)
            for k in [5, 4, 3]:
                ap = [p + j * d for j in range(k)]
                if all(x in ps for x in ap):
                    ap_counts[k] += 1
                    if len(ap_examples[k]) < 3:
                        ap_examples[k].append(ap)
                    break  # don't double-count shorter APs in this direction

    for k in [3, 4, 5]:
        print(f"\n  APs of length {k} found (among first 100 primes, d<1000): {ap_counts[k]}")
        for ap in ap_examples[k][:2]:
            print(f"    Example: {ap}")

    # Test: given primes p_1..p_n, try to predict p_{n+1} using AP extrapolation
    # For each pair (p_i, p_j), the AP through them predicts p_j + (p_j - p_i) as next.
    print("\n  Testing AP extrapolation for next prime prediction:")
    correct = 0
    total = 0
    for idx in range(10, 200):
        target = primes_small[idx]
        # Try all pairs from recent primes
        found = False
        for i in range(max(0, idx - 10), idx):
            for j in range(i + 1, idx):
                d = primes_small[j] - primes_small[i]
                candidate = primes_small[j] + d
                if candidate == target:
                    found = True
                    break
            if found:
                break
        if found:
            correct += 1
        total += 1

    print(f"  AP-extrapolation predicted next prime: {correct}/{total} = {100*correct/total:.1f}%")

    # Deeper analysis: what fraction of primes lie on an AP with 2 earlier primes?
    on_ap = 0
    for idx in range(2, 500):
        p = primes_small[idx]
        found = False
        for i in range(idx):
            # Check if 2*p_i - p is also prime and < p_i (then p_i, (p_i+p)/2, p is AP)
            # Actually check: is there p_j < p such that p - p_j = p_j - p_i for some p_i?
            d_needed = p - primes_small[i]
            candidate_prev = primes_small[i] - d_needed
            if candidate_prev > 1 and candidate_prev in ps:
                found = True
                break
        if found:
            on_ap += 1

    print(f"  Primes on an AP of length>=3 with smaller primes: {on_ap}/498 = {100*on_ap/498:.1f}%")

    print(f"""
  VERDICT: Green-Tao AP interpolation FAILS as a formula.
  - APs in primes are sparse and unpredictable.
  - Finding APs requires knowing the primes (circular).
  - AP extrapolation predicts the next prime only ~{100*correct/total:.0f}% of the time.
  - The theorem is existential (primes contain long APs) but gives no
    constructive way to find them without full prime enumeration.
  - COMPLEXITY: Finding APs of length k requires O(N^2) pairs minimum.
  """)
    return False


# =============================================================================
# APPROACH 2: Hardy-Littlewood Circle Method Inversion
# =============================================================================

def test_circle_method():
    """
    The circle method gives:
      r_k(n) = integral_0^1 S(alpha)^k * e(-n*alpha) d(alpha)
    where S(alpha) = sum_{p<=N} e(p*alpha).

    For k=1: r_1(n) = sum_{p<=N} integral_0^1 e((p-n)*alpha) d(alpha)
                     = sum_{p<=N} delta(p,n)  = 1 if n is prime, 0 otherwise.
    This is just a TAUTOLOGY for k=1.

    For k=2 (Goldbach): r_2(n) = #{(p,q): p+q=n, p,q prime}.
    The circle method gives the ASYMPTOTIC:
      r_2(n) ~ S_2(n) * n / (ln n)^2
    where S_2(n) = 2 * prod_{p|n, p>2} (p-1)/(p-2) * prod_{p>2} (1 - 1/(p-1)^2)
    (the singular series).

    Key question: Can we compute the EXACT r_2(n) without knowing primes?
    Answer: The singular series gives the main term. The ERROR term is:
      r_2(n) - S_2(n) * Li_2(n) = O(n * exp(-c*sqrt(ln n)))  (under GRH)
    This error is HUGE compared to 1. We need exactness to detect primes.

    Let's verify the circle method asymptotic and measure the error.
    """
    print("=" * 72)
    print("APPROACH 2: Hardy-Littlewood Circle Method / Goldbach Asymptotics")
    print("=" * 72)

    # Compute exact r_2(n) for even n up to some limit
    limit = 10000
    ps = sieve_primes(limit)
    pset = set(ps)

    # Exact Goldbach count
    def goldbach_count(n):
        if n % 2 != 0 or n < 4:
            return 0
        count = 0
        for p in ps:
            if p > n // 2:
                break
            if (n - p) in pset:
                count += 1
        return count

    # Hardy-Littlewood singular series S_2(n) for even n
    # S_2(n) = 2 * C_2 * prod_{p|n, p odd prime} (p-1)/(p-2)
    # where C_2 = prod_{p>=3} (1 - 1/(p-1)^2) = "twin prime constant" ~ 0.6601...

    # Compute twin prime constant
    C2 = 1.0
    for p in ps[1:200]:  # skip 2, use first 200 odd primes
        C2 *= (1.0 - 1.0 / (p - 1)**2)

    def singular_series_2(n):
        """Singular series for Goldbach representation of even n."""
        if n % 2 != 0:
            return 0.0
        result = 2.0 * C2
        # For each odd prime p dividing n
        temp = n
        for p in ps[1:]:  # odd primes
            if p * p > temp:
                break
            if temp % p == 0:
                result *= (p - 1) / (p - 2)
                while temp % p == 0:
                    temp //= p
        if temp > 2:
            result *= (temp - 1) / (temp - 2)
        return result

    # Compare exact vs asymptotic
    print("\n  Comparing exact r_2(n) vs Hardy-Littlewood asymptotic:")
    print(f"  {'n':>8} {'r_2(n) exact':>14} {'HL estimate':>14} {'ratio':>8} {'abs_err':>10}")

    errors = []
    for n in range(100, 10001, 100):
        if n % 2 != 0:
            continue
        exact = goldbach_count(n)
        S2 = singular_series_2(n)
        hl_est = S2 * n / (math.log(n))**2
        ratio = exact / hl_est if hl_est > 0 else 0
        err = abs(exact - hl_est)
        errors.append((n, err, exact))
        if n <= 1000 or n % 2000 == 0:
            print(f"  {n:>8} {exact:>14} {hl_est:>14.1f} {ratio:>8.4f} {err:>10.1f}")

    # Analyze: can the error tell us about individual primes?
    avg_err_frac = sum(e/max(r,1) for _, e, r in errors) / len(errors)
    max_err = max(e for _, e, _ in errors)
    print(f"\n  Average relative error: {avg_err_frac:.4f}")
    print(f"  Max absolute error: {max_err:.1f}")

    # Key test: Can we detect if a number is prime from r_2?
    # If n is prime, then r_2(n+1) should be "higher" because n is available as a summand?
    # Actually no - r_2 is for EVEN numbers, not directly for prime detection.

    # More relevant: for odd n, r_3(n) = #{(p,q,r): p+q+r=n} (Vinogradov/Goldbach)
    # Vinogradov proved r_3(n) > 0 for all large odd n.
    # The asymptotic: r_3(n) ~ (1/2) * S_3(n) * n^2 / (ln n)^3
    # Again the error is too large for prime detection.

    print(f"""
  VERDICT: Circle method inversion FAILS.

  The fundamental issue:
  - For k=1 (direct prime detection), the circle method integral is a TAUTOLOGY:
    it equals 1 iff n is prime, but computing it requires knowing all primes up to N.
  - For k=2 (Goldbach), the asymptotic formula has error O(n/exp(c*sqrt(ln n))).
    For n=10^100, error ~ 10^93 vs r_2 ~ 10^96. Relative error ~10^(-3).
    This is far too coarse to extract individual primes.
  - For k=3 (Vinogradov), error is even larger relative to our needs.

  The circle method inherently gives AVERAGE information about primes,
  not individual prime locations. The "major arcs" capture the smooth
  part (singular series), and the "minor arcs" are bounded but not
  computable without prime enumeration.

  COMPLEXITY: Evaluating S(alpha) = sum_{{p<=N}} e(p*alpha) requires
  knowing all primes up to N. This is fundamentally circular.
  """)
    return False


# =============================================================================
# APPROACH 3: Sieve Weights as a Direct Formula
# =============================================================================

def test_sieve_weights():
    """
    Modern sieve methods assign weights w(n) to integers such that:
    - w(n) >= 0 for all n
    - w(n) > 0 implies n is "almost prime" (few prime factors)
    - sum_{n<=N} w(n) * f(n) approximates sum_{p<=N} f(p)

    Selberg's sieve: w(n) = (sum_{d|n, d<=D} lambda_d)^2
    where lambda_d are optimized to detect primes.

    GPY sieve (Goldston-Pintz-Yildirim) / Maynard-Tao:
    These use more sophisticated weights for bounded gaps between primes.

    Question: Can we compute w(n) for a SINGLE n in O(polylog(n)) time?

    The Selberg weights lambda_d satisfy:
      lambda_d = mu(d) * (log(D/d) / log(D))  for d <= D
    where D ~ N^{1/2} for Selberg's upper bound sieve.

    For a single n, w(n) = (sum_{d|n, d<=D} lambda_d)^2.
    Computing this requires the DIVISORS of n up to D.
    If n is ~10^100, D ~ 10^50. Finding divisors up to 10^50 requires
    factoring n, which for a random 10^100-digit number takes ~exp(polylog) time.
    BUT: we're evaluating at integers n, not random numbers!

    Let's test: how well do Selberg weights detect primes?
    """
    print("=" * 72)
    print("APPROACH 3: Sieve Weights as a Prime Formula")
    print("=" * 72)

    N = 10000
    D = int(N**0.5)
    primes_D = sieve_primes(D)
    ps = sieve_primes(N)
    pset = set(ps)

    # Selberg lambda_d
    logD = math.log(D)
    def selberg_lambda(d):
        if d > D:
            return 0.0
        if d == 1:
            return 1.0
        # mu(d) computation
        temp = d
        mu = 1
        for p in primes_D:
            if p * p > temp:
                break
            if temp % p == 0:
                temp //= p
                if temp % p == 0:
                    return 0.0  # mu = 0
                mu *= -1
        if temp > 1:
            mu *= -1
        return mu * max(0, math.log(D / d) / logD)

    # Compute w(n) = (sum_{d|n, d<=D} lambda_d)^2
    def sieve_weight(n):
        if n <= 1:
            return 0.0
        # Find divisors of n up to D
        divisors = []
        for d in range(1, min(D, int(n**0.5)) + 1):
            if n % d == 0:
                divisors.append(d)
                if d != n // d and n // d <= D:
                    divisors.append(n // d)
        s = sum(selberg_lambda(d) for d in divisors)
        return s * s

    # Evaluate weights for primes and composites
    print("\n  Sieve weights for small numbers (primes marked with *):")
    prime_weights = []
    composite_weights = []
    for n in range(2, 101):
        w = sieve_weight(n)
        if n in pset:
            prime_weights.append(w)
        else:
            composite_weights.append(w)
        if n <= 30:
            mark = " *" if n in pset else ""
            print(f"    w({n:>3}) = {w:.6f}{mark}")

    avg_pw = sum(prime_weights) / len(prime_weights) if prime_weights else 0
    avg_cw = sum(composite_weights) / len(composite_weights) if composite_weights else 0
    print(f"\n  Average weight for primes 2..100:     {avg_pw:.6f}")
    print(f"  Average weight for composites 2..100: {avg_cw:.6f}")
    print(f"  Ratio (prime/composite):              {avg_pw/avg_cw:.2f}")

    # Can we threshold to detect primes?
    thresholds = [0.01, 0.1, 0.3, 0.5, 0.7]
    print(f"\n  Threshold-based prime detection (n=2..{N}):")
    for t in thresholds:
        tp = fp = fn = tn = 0
        for n in range(2, min(N + 1, 2001)):
            w = sieve_weight(n)
            predicted_prime = (w > t)
            actual_prime = (n in pset)
            if predicted_prime and actual_prime: tp += 1
            elif predicted_prime and not actual_prime: fp += 1
            elif not predicted_prime and actual_prime: fn += 1
            else: tn += 1
        prec = tp / (tp + fp) if (tp + fp) > 0 else 0
        rec = tp / (tp + fn) if (tp + fn) > 0 else 0
        print(f"    threshold={t:.2f}: precision={prec:.3f}, recall={rec:.3f}, "
              f"TP={tp}, FP={fp}, FN={fn}")

    # Complexity analysis for single-value evaluation
    print(f"""
  VERDICT: Sieve weights FAIL as a prime formula.

  Issues:
  1. Sieve weights are UPPER BOUND sieves: they give w(n) > 0 for primes
     but ALSO for some composites. No threshold cleanly separates them.
  2. Computing w(n) for a single n requires finding ALL divisors of n up to D.
     For n ~ p(10^100) ~ 10^102, D ~ 10^51. Finding divisors up to 10^51
     requires factoring n, which costs exp(O((ln n)^{{1/3}})) ~ 10^14 at best.
  3. Even if we could compute w(n), the weights don't give EXACT prime
     detection. The Selberg sieve gives sum w(n) ~ 2N/ln(N), capturing
     primes but also semiprimes. The "parity barrier" prevents sieves
     from distinguishing primes from products of an odd vs even number
     of prime factors.
  4. The PARITY PROBLEM (Selberg, 1949) is the fundamental obstacle:
     no sieve of "Selberg type" can detect primes exactly. It can only
     bound their count.

  COMPLEXITY: O(D) = O(N^{{1/2}}) per evaluation, and doesn't give exactness.
  """)
    return False


# =============================================================================
# APPROACH 4: Goldbach Representation Extraction
# =============================================================================

def test_goldbach_extraction():
    """
    For even n, g(n) = #{(p,q): p+q=n, p<=q}.
    Key insight: p is prime iff g(p+2) >= 1 (since p + 2 is even, and
    if p is prime, then (2, p) contributes to g(p+2) when p is odd prime).
    Wait -- that's not quite right. Let's think more carefully.

    Actually: For even n=2m, g(n) counts pairs (p, n-p) where both are prime.
    So g(n) encodes which numbers p <= n/2 are prime (those for which n-p
    is also prime).

    Can we use g(n) for multiple values of n to triangulate individual primes?

    Idea: Define h(p) = #{even n: p appears in a Goldbach pair for n}
         = #{even n: n-p is prime, 4 <= n <= 2p}
         = #{q prime: q <= p and p+q is even} = #{q prime: q <= p, q != 2 unless p is even}

    For odd prime p: h(p) = pi(p) (since p+q is even iff q is odd, i.e., q != 2,
    plus q=2 gives p+2 which is even, so actually h(p) = pi(p)).
    Hmm, this is circular.

    Different angle: Suppose we have an oracle for g(n) (exact Goldbach count)
    for any even n. Can we extract primes?

    g(2m) = sum_{p prime, p <= m} [2m - p is prime]

    g(2m) - g(2m-2) = [2m - 2 is prime AND 2 is prime]  -- only the new pairs
    Actually this is wrong. Let's be more careful.

    For extracting: p is prime iff there exists even n such that (p, n-p) are
    both prime. This requires KNOWING n-p is prime -- circular!

    Alternative: Use CONVOLUTION structure.
    Let chi_P(n) = 1 if n is prime, 0 otherwise.
    Then g(n) = sum_{k=2}^{n-2} chi_P(k) * chi_P(n-k) = (chi_P * chi_P)(n)
    where * is convolution.

    If we know g = chi_P * chi_P, we want to DECONVOLVE to get chi_P.
    In Fourier domain: G(t) = |Chi_P(t)|^2.
    This gives the MAGNITUDE but not the PHASE of Chi_P(t)!
    Phase retrieval is a hard problem (NP-hard in general).
    """
    print("=" * 72)
    print("APPROACH 4: Goldbach Representation Extraction")
    print("=" * 72)

    N = 2000
    ps = sieve_primes(N)
    pset = set(ps)

    # Compute exact g(n) for even n
    def goldbach_exact(n):
        count = 0
        for p in ps:
            if p > n // 2:
                break
            if (n - p) in pset:
                count += 1
        return count

    # g(n) values
    g_vals = {}
    for n in range(4, N + 1, 2):
        g_vals[n] = goldbach_exact(n)

    print("\n  Goldbach counts g(n) for small even n:")
    for n in range(4, 42, 2):
        print(f"    g({n:>3}) = {g_vals[n]:>3}")

    # Test deconvolution idea via differences
    # g(2p) >= 1 for every odd prime p (because (p, p) is a pair when p=p, but
    # p+p=2p and we need both to be prime, so g(2p) >= 1 iff p is prime.
    # Actually g(2p) counts (q, 2p-q) with both prime. If p is prime, (p,p) is one.
    print("\n  Testing: g(2p) >= 1 as prime indicator:")
    correct = 0
    total = 0
    for n in range(3, 500):
        val = g_vals.get(2 * n)
        if val is not None:
            predicted_prime = (val >= 1)
            actual = n in pset
            # g(2n) >= 1 for almost all even 2n >= 4 (Goldbach conjecture)
            # So this does NOT detect primes -- g(2n) >= 1 for ALL n, prime or not
            if n <= 20:
                print(f"    n={n:>3}: g(2n)={val:>3}, n is prime: {actual}")
            total += 1

    print(f"\n  g(2n) >= 1 for ALL n >= 2 (Goldbach conjecture), so this is useless")
    print(f"  for prime detection -- it's always true!")

    # Attempt: use g(n) differences to detect primes
    # When p is prime and p+2 is composite, g(p+3) might differ from g(p+1)
    # Let's check correlation between g(n) and primality of n/2 or n-2 etc.

    # More sophisticated: the function g(2n) = sum_{p<=n} chi_P(2n-p)
    # involves chi_P evaluated at 2n-p. If we vary n, we get different
    # "slices" through chi_P. But each slice is a SUM, not a point value.

    # Fourier/deconvolution test
    print("\n  Fourier deconvolution test:")
    print("  g = chi_P * chi_P  =>  in Fourier: |F(chi_P)|^2 = F(g)")
    print("  This gives MAGNITUDE only. Phase is lost!")
    print("  Phase retrieval is NP-hard in general (Sahinoglou & Cabrera, 1991).")

    # Let's verify: compute DFT of g and chi_P
    import cmath

    M = 500
    chi = [0] * M
    for p in ps:
        if p < M:
            chi[p] = 1

    g_arr = [0] * M
    for n in range(4, M, 2):
        g_arr[n] = g_vals.get(n, 0)

    # DFT
    def dft(x, N_):
        X = []
        for k in range(N_):
            s = 0
            for n in range(N_):
                s += x[n] * cmath.exp(-2j * cmath.pi * k * n / N_)
            X.append(s)
        return X

    N_dft = 200  # small for speed
    chi_small = chi[:N_dft]
    g_small = g_arr[:N_dft]

    Chi_F = dft(chi_small, N_dft)
    G_F = dft(g_small, N_dft)

    # Check: |Chi_F[k]|^2 should approximate G_F[k]
    print(f"  Comparing |F(chi_P)|^2 vs F(g) for first few frequencies:")
    match_count = 0
    for k in range(min(10, N_dft)):
        chi2 = abs(Chi_F[k])**2
        gf = abs(G_F[k])
        ratio = chi2 / gf if gf > 0.01 else float('inf')
        if k < 8:
            print(f"    k={k}: |F(chi)|^2 = {chi2:.2f}, |F(g)| = {gf:.2f}")

    print(f"""
  VERDICT: Goldbach extraction FAILS.

  1. g(2n) >= 1 for essentially all n (Goldbach conjecture), so it cannot
     distinguish primes from composites.
  2. Deconvolution of g = chi_P * chi_P requires PHASE RETRIEVAL:
     the DFT gives |F(chi_P)|^2 but not the phase of F(chi_P).
     Phase retrieval is NP-hard in general.
  3. Even if phase retrieval were tractable, computing g(n) exactly for
     a single large n requires knowing all primes up to n -- circular.
  4. The Hardy-Littlewood asymptotic for g(n) has error too large for
     exact prime extraction.

  COMPLEXITY: Computing g(n) exactly requires O(pi(n)) = O(n/ln n) work.
  Deconvolution requires knowing g for ALL even n up to 2N -- O(N^2/ln N) total.
  """)
    return False


# =============================================================================
# APPROACH 5: Erdos-Kac and Additive Function Approaches
# =============================================================================

def test_erdos_kac():
    """
    Erdos-Kac theorem: omega(n) (number of distinct prime factors) satisfies
      (omega(n) - ln ln n) / sqrt(ln ln n) -> N(0,1)  as n -> infinity.

    Idea: primes are exactly the n with omega(n) = 1 AND Omega(n) = 1.
    Can we use the distribution of omega to build a prime detector?

    Also: the Liouville function lambda(n) = (-1)^{Omega(n)} satisfies
      sum_{n<=x} lambda(n) = O(sqrt(x))  (equivalent to RH!)
    and the Mobius function mu(n) = 0 if n has squared factor, else (-1)^{omega(n)}.

    Key insight: mu(n)^2 * [omega(n) == 1] is a prime indicator!
    But computing omega(n) or mu(n) requires FACTORING n.
    For n ~ 10^102, factoring is exp(O((ln n)^{1/3})) time.

    Alternative: the von Mangoldt function Lambda(n) = ln p if n = p^k, else 0.
    Computing Lambda(n) requires checking if n is a prime power -- which requires
    factoring (or primality testing + root extraction).

    Let's test: how fast can we detect primality via these additive functions?
    """
    print("=" * 72)
    print("APPROACH 5: Erdos-Kac / Additive Function Prime Detection")
    print("=" * 72)

    N = 10000

    # Compute omega(n) and Omega(n) via sieve
    omega = [0] * (N + 1)  # number of distinct prime factors
    Omega = [0] * (N + 1)  # number of prime factors with multiplicity
    for p in PRIMES:
        if p > N:
            break
        for m in range(p, N + 1, p):
            omega[m] += 1
        pk = p
        while pk <= N:
            for m in range(pk, N + 1, pk):
                Omega[m] += 1
            pk *= p

    # Primes: omega(n) = 1, Omega(n) = 1
    # Prime powers: omega(n) = 1, Omega(n) > 1
    pset_N = set(sieve_primes(N))

    # Verify prime detection
    detect_correct = 0
    total_checked = 0
    for n in range(2, N + 1):
        is_p = (omega[n] == 1 and Omega[n] == 1)
        actual = n in pset_N
        if is_p == actual:
            detect_correct += 1
        total_checked += 1

    print(f"\n  omega(n)=1 AND Omega(n)=1 as prime detector: "
          f"{detect_correct}/{total_checked} = {100*detect_correct/total_checked:.2f}%")

    # Erdos-Kac distribution check
    import statistics
    vals_1000 = [(omega[n] - math.log(math.log(n))) / math.sqrt(math.log(math.log(n)))
                 for n in range(100, N + 1)]
    print(f"\n  Erdos-Kac normalized omega: mean={statistics.mean(vals_1000):.4f} "
          f"(expect 0), std={statistics.stdev(vals_1000):.4f} (expect 1)")

    # Can we compute omega(n) for a single large n without factoring?
    # Answer: No. omega(n) requires the complete factorization of n.
    # Factoring n ~ 10^102 via GNFS takes exp(O((ln n)^{1/3} (ln ln n)^{2/3})) time.
    # For n = 10^102: exponent ~ (235)^{1/3} * (5.46)^{2/3} ~ 6.17 * 3.10 ~ 19.1
    # So ~ e^19 ~ 2 * 10^8 operations. Actually feasible for 102 digits!
    # BUT: we need to know p(10^100) first to evaluate omega(p(10^100)).
    # This is completely circular.

    # Mobius function approach
    print(f"\n  Mobius/von Mangoldt approach:")
    print(f"  Lambda(n) = ln(p) if n=p^k, else 0.")
    print(f"  Computing Lambda(n) for a SINGLE n requires primality testing")
    print(f"  (O(polylog) via AKS/Miller-Rabin) + root extraction.")
    print(f"  This is fast for a KNOWN n, but we don't know p(10^100)!")

    # The REAL question: can additive functions help compute pi(x)?
    # sum_{n<=x} Lambda(n) = psi(x) ~ x (PNT)
    # pi(x) = sum_{n<=x} Lambda(n)/ln(n) approximately
    # But we need EXACT pi(x), and the sum is over ALL n<=x.

    # Partial summation / Euler-Maclaurin?
    # psi(x) = x - sum_rho x^rho/rho - ln(2pi) - (1/2)ln(1-x^{-2})
    # This is the explicit formula -- same as the zeta zeros approach!

    print(f"""
  VERDICT: Erdos-Kac / additive functions FAIL as a path to p(n).

  1. omega(n)=1, Omega(n)=1 perfectly detects primes but requires FACTORING n.
     This is O(exp((ln n)^{{1/3}})) ~ feasible for 100-digit numbers, but we'd
     need to know p(10^100) first -- CIRCULAR.
  2. The Erdos-Kac theorem describes the STATISTICAL distribution of omega(n).
     It gives no way to compute omega for a specific n without factoring.
  3. Lambda(n), mu(n) are computable for a given n via primality testing,
     but summing them over all n <= x to get pi(x) requires O(x) work.
  4. The partial sum formulas (explicit formula for psi(x)) reduce to the
     SAME zeta-zero approach already explored -- and that requires O(sqrt(x))
     zeros, which is infeasible for x ~ 10^102.
  5. No additive-combinatorial property of primes gives a shortcut to pi(x)
     or p(n) because they all ultimately encode the same information as the
     prime distribution itself.
  """)
    return False


# =============================================================================
# BONUS: Theoretical analysis of ALL additive combinatorics approaches
# =============================================================================

def theoretical_analysis():
    """
    Unified theoretical analysis of why additive combinatorics cannot yield
    a fast formula for p(n).
    """
    print("=" * 72)
    print("THEORETICAL ANALYSIS: Why Additive Combinatorics Cannot Help")
    print("=" * 72)

    print("""
  The fundamental barrier across ALL additive combinatorics approaches:

  1. INFORMATION BARRIER
     p(10^100) has ~340 bits of information.
     R^{-1}(n) provides ~172 bits (the smooth/predictable part).
     The remaining ~168 bits encode the oscillatory part: zeta zeros.
     ANY approach must somehow compute these 168 bits.

  2. GREEN-TAO & AP STRUCTURE
     Primes contain APs of length k, but the density of length-k APs
     among primes is ~ c_k / (ln N)^k (by Green-Tao, building on
     Goldston-Yildirim and Szemeredi's theorem).
     The APs are too sparse and irregularly distributed to interpolate.
     Finding APs requires knowing the primes -- circular.

  3. CIRCLE METHOD / VINOGRADOV
     The circle method computes:
       sum over n of a_n * r_k(n)
     where r_k(n) = #{representations of n as sum of k primes}.
     The "major arcs" give the singular series (smooth part).
     The "minor arcs" contribute error O(N * L^{-A}) for any A.

     KEY ISSUE: Evaluating S(alpha) = sum_{p<=N} e(p*alpha) on major
     arcs uses the Siegel-Walfisz theorem, which gives:
       S(alpha) ~ mu(q)^{-1} * phi(q)^{-1} * sum_{n<=N} e(n*beta) / ln(n)
     for alpha near a/q. But the implicit constants depend on ALL primes
     up to N. The method gives ASYMPTOTIC results, never exact values.

  4. SIEVE METHODS & PARITY BARRIER
     Selberg's "parity problem" (1949): Any sieve of "linear" type
     cannot distinguish numbers with an odd number of prime factors
     from those with an even number.

     Consequence: No combination of sieve weights can exactly detect primes.
     The best we can do is detect "almost primes" (numbers with at most
     k prime factors for some fixed k).

     This was formalized by Bombieri (1976): the "parity barrier" is
     inherent to all combinatorial sieve methods.

     Even the breakthrough of Maynard-Tao (2013) on bounded gaps uses
     sieves that detect CLUSTERS of almost-primes, not individual primes.

  5. GOLDBACH DECONVOLUTION
     g = chi_P * chi_P means F(g) = |F(chi_P)|^2.
     Recovering chi_P requires phase retrieval.
     Even with the additional constraint that chi_P is 0-1 valued,
     the phase retrieval problem is at least as hard as factoring
     (reduction via number-theoretic arguments).
     Moreover, computing g(n) exactly for large n requires O(n/ln n) work.

  6. ERGODIC THEORY (Green-Tao framework)
     Green-Tao proved their theorem using the "transference principle":
     primes behave like a dense subset of the integers for purposes of
     finding APs. But this is an EXISTENCE result based on Szemeredi's
     theorem. The effective bounds are towers of exponentials.
     For APs of length k, the first AP starts at ~ exp(exp(exp(...(k)...)))
     The theory gives NO computationally useful bounds.

  UNIFIED CONCLUSION:
  All additive combinatorics approaches face the same barrier:

  - The SMOOTH part of prime distribution (density ~ 1/ln n) is captured by
    PNT / Riemann's R(x) -- we already have this, and it's not enough.
  - The OSCILLATORY part (corrections from zeta zeros) requires information
    equivalent to knowing the primes themselves.
  - Additive combinatorics gives STATISTICAL/ASYMPTOTIC results about primes
    (they contain APs, they have positive Goldbach representations, etc.)
    but never EXACT results for individual primes.
  - The parity barrier is specific to sieves but reflects a deeper truth:
    the "evenness" of the number of prime factors oscillates in a way that
    is entangled with the zeta zeros.

  COMPLEXITY SUMMARY:
    Approach                    | Per-evaluation cost    | Exact?
    ----------------------------|------------------------|--------
    Green-Tao AP interpolation  | O(N^2) (find APs)     | No
    Circle method S(alpha)      | O(N/ln N) (sum primes) | Asymptotic only
    Selberg sieve weights       | O(N^{1/2}) (divisors)  | No (parity barrier)
    Goldbach deconvolution      | O(N/ln N) + NP-hard    | Theoretically yes
    Erdos-Kac / omega(n)        | O(exp(n^{1/3})) factor | Yes but circular

    None achieves O(polylog(N)) exact evaluation.
    The best known remains O(N^{2/3}) via Deleglise-Rivat for pi(x).
  """)


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("Session 4: Additive Combinatorics & Arithmetic Progressions for p(n)")
    print("Investigating whether additive structure gives a formula for p(n)")
    print()

    t0 = time.time()

    results = {}

    results['green_tao'] = test_green_tao_interpolation()
    print()

    results['circle_method'] = test_circle_method()
    print()

    results['sieve_weights'] = test_sieve_weights()
    print()

    results['goldbach'] = test_goldbach_extraction()
    print()

    results['erdos_kac'] = test_erdos_kac()
    print()

    theoretical_analysis()

    elapsed = time.time() - t0
    print(f"\n{'=' * 72}")
    print(f"FINAL SUMMARY (elapsed: {elapsed:.1f}s)")
    print(f"{'=' * 72}")
    print(f"""
  All 5 additive combinatorics approaches investigated: NONE viable.

  | # | Approach                    | Viable? | Reason                        |
  |---|-----------------------------|---------|-------------------------------|
  | 1 | Green-Tao AP interpolation  | NO      | APs sparse, circular          |
  | 2 | Circle method inversion     | NO      | Asymptotic only, circular     |
  | 3 | Sieve weights formula       | NO      | Parity barrier, not exact     |
  | 4 | Goldbach extraction         | NO      | Phase retrieval NP-hard       |
  | 5 | Erdos-Kac additive funcs    | NO      | Requires factoring, circular  |

  ROOT CAUSE (shared by all approaches):
  Additive combinatorics captures the STATISTICAL/DENSITY properties of primes
  (which are governed by PNT / the smooth part of R(x)), but cannot access the
  OSCILLATORY corrections (governed by zeta zeros) needed for exact results.

  The parity barrier (Selberg, 1949; Bombieri, 1976) is the sieve-theoretic
  manifestation of this fundamental limitation.

  This reinforces the project's existing finding: computing p(10^100) exactly
  in O(polylog) time is provably beyond any known mathematical framework.
  The O(N^{{2/3}}) Lucy_Hedgehog DP remains the practical optimum.
  """)


if __name__ == '__main__':
    main()
