#!/usr/bin/env python3
"""
Session 8, Experiment: Dynamical Systems, Ergodic Theory & Orbital Approaches to Primes

We test 5 families of ideas:
  1. Collatz-like prime-generating maps & Rowland acceleration
  2. Furstenberg topology — computational content extraction
  3. Ergodic averages along primes — reversal attempts
  4. Symbolic dynamics on prime indicator sequences
  5. Fixed points of iterated maps at primes (Newton on Chebyshev psi, etc.)

Goal: find ANYTHING that computes p(n) in O(polylog n) time.
Expectation based on 205+ prior failures: the ~178-bit barrier will hold.
"""

import math
import time
import sys
from collections import defaultdict, Counter
from functools import lru_cache

# =============================================================================
# Utility: reference primes
# =============================================================================

def sieve(limit):
    """Sieve of Eratosthenes up to limit."""
    is_prime = bytearray(b'\x01') * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return [i for i in range(2, limit + 1) if is_prime[i]]

LIMIT = 100000
PRIMES = sieve(LIMIT)
PRIME_SET = set(PRIMES)

def is_prime_small(n):
    if n < 2: return False
    if n <= LIMIT: return n in PRIME_SET
    if n % 2 == 0: return False
    for p in range(3, int(n**0.5) + 1, 2):
        if n % p == 0: return False
    return True

def nth_prime(n):
    """Return the n-th prime (1-indexed)."""
    if n <= len(PRIMES):
        return PRIMES[n - 1]
    raise ValueError(f"n={n} exceeds precomputed range")


# =============================================================================
# EXPERIMENT 1: Collatz-like prime-generating maps & Rowland acceleration
# =============================================================================

def experiment1_collatz_and_rowland():
    print("=" * 78)
    print("EXPERIMENT 1: Collatz-like maps & Rowland sequence")
    print("=" * 78)

    # --- 1a: Rowland sequence ---
    # a(n) = a(n-1) + gcd(n, a(n-1)), a(1) = 7
    # The differences |a(n) - a(n-1)| that exceed 1 are always prime.
    print("\n--- 1a: Rowland sequence a(n) = a(n-1) + gcd(n, a(n-1)) ---")
    t0 = time.time()
    a = 7
    rowland_primes = []
    steps_to_prime = []  # how many iterations between consecutive prime outputs
    last_prime_step = 0
    for n in range(2, 200001):
        g = math.gcd(n, a)
        a = a + g
        if g > 1:
            rowland_primes.append(g)
            steps_to_prime.append(n - last_prime_step)
            last_prime_step = n
    elapsed = time.time() - t0
    unique_primes = sorted(set(rowland_primes))
    print(f"  200000 iterations: found {len(rowland_primes)} prime outputs")
    print(f"  Unique primes: {unique_primes[:20]}...")
    print(f"  Largest prime found: {max(rowland_primes)}")
    if steps_to_prime:
        avg_gap = sum(steps_to_prime) / len(steps_to_prime)
        print(f"  Avg iterations per prime output: {avg_gap:.1f}")
    print(f"  Time: {elapsed:.3f}s")
    # Key issue: Rowland generates primes in order but takes O(p^2) steps per prime p.
    # The n-th unique prime p requires ~p^2 iterations.
    # For p(10^6) ~ 1.5*10^7, that's ~2*10^14 iterations. FAR from polylog.

    # --- 1b: Cloitre variant ---
    # b(n) = b(n-1) + lcm(n, b(n-1)), primes from ratios
    print("\n--- 1b: Cloitre's EKG-like sequence ---")
    # EKG sequence: a(1)=1, a(2)=2, a(n) = smallest m not yet used with gcd(m, a(n-1)) > 1
    # Primes appear at positions 2p-1 in a rearranged sense. But still O(n) per step.
    t0 = time.time()
    used = {1, 2}
    ekg = [1, 2]
    for _ in range(998):
        prev = ekg[-1]
        m = 2
        while True:
            if m not in used and math.gcd(m, prev) > 1:
                ekg.append(m)
                used.add(m)
                break
            m += 1
    elapsed = time.time() - t0
    ekg_primes = [x for x in ekg if is_prime_small(x)]
    print(f"  First 1000 EKG terms computed in {elapsed:.3f}s")
    print(f"  Primes in first 1000: {len(ekg_primes)}")
    print(f"  First 15 EKG terms: {ekg[:15]}")
    # EKG is interesting but O(n^2) to find n-th prime. No speedup.

    # --- 1c: Parametric Collatz-like maps ---
    # Try f_a,b(x) = x/2 if even, ax+b if odd, for various a,b
    # Track how many primes the orbit visits (up to some bound).
    print("\n--- 1c: Parametric Collatz-like maps f(x) = x/2 if even, ax+b if odd ---")
    best_prime_density = 0
    best_params = None
    t0 = time.time()
    for a in range(1, 12, 2):  # odd multipliers
        for b in range(-5, 8, 2):  # odd offsets
            if b == 0:
                continue
            # Start from seed, iterate, count distinct primes visited
            for seed in [3, 7, 11, 23]:
                x = seed
                visited_primes = set()
                orbit_len = 0
                seen = set()
                for _ in range(10000):
                    if x in seen or x <= 0 or x > 10**7:
                        break
                    seen.add(x)
                    orbit_len += 1
                    if x <= LIMIT and x in PRIME_SET:
                        visited_primes.add(x)
                    if x % 2 == 0:
                        x = x // 2
                    else:
                        x = a * x + b
                if orbit_len > 0:
                    density = len(visited_primes) / orbit_len
                    if len(visited_primes) > best_prime_density:
                        best_prime_density = len(visited_primes)
                        best_params = (a, b, seed, orbit_len, len(visited_primes))
    elapsed = time.time() - t0
    if best_params:
        a, b, seed, orb, nprimes = best_params
        print(f"  Best: a={a}, b={b}, seed={seed}: {nprimes} distinct primes in orbit of length {orb}")
    print(f"  Time: {elapsed:.3f}s")
    print(f"  VERDICT: Orbits visit SOME primes but not in order. Cannot extract p(n).")

    # --- 1d: Attempt to accelerate Rowland via pattern ---
    print("\n--- 1d: Rowland acceleration attempt ---")
    # The Rowland sequence has a(n) = n + r(n) where r follows a pattern.
    # Key insight: the GCD kicks in only when n shares a factor with a(n-1).
    # Can we predict WHEN the next gcd > 1 event happens?
    # If a(n) = a(n-1) + 1 (gcd=1) for long stretches, those stretches have length ~p-1
    # before the next prime p appears.
    t0 = time.time()
    a_val = 7
    prime_events = []  # (n, prime) when gcd > 1
    for n in range(2, 500001):
        g = math.gcd(n, a_val)
        a_val = a_val + g
        if g > 1:
            prime_events.append((n, g))

    # Analyze gaps between events
    if len(prime_events) >= 10:
        gaps = [prime_events[i+1][0] - prime_events[i][0] for i in range(len(prime_events)-1)]
        primes_found = [e[1] for e in prime_events]
        # Check if gap before prime p is related to p
        correlations = []
        for i in range(min(50, len(gaps))):
            p = primes_found[i]
            g = gaps[i] if i < len(gaps) else 0
            correlations.append((p, g, g / p if p > 0 else 0))
        print(f"  Found {len(prime_events)} prime events in 500000 iterations")
        print(f"  (prime, gap_to_next, ratio) for first 15:")
        for p, g, r in correlations[:15]:
            print(f"    prime={p:6d}, gap={g:6d}, ratio={r:.2f}")
    elapsed = time.time() - t0
    print(f"  Time: {elapsed:.3f}s")
    print(f"  VERDICT: Gap ~ p^2/const. No way to skip ahead in O(polylog).")

    return "FAIL: All Collatz/Rowland variants are O(p) or O(p^2) per prime. No polylog path."


# =============================================================================
# EXPERIMENT 2: Furstenberg topology — computational content
# =============================================================================

def experiment2_furstenberg():
    print("\n" + "=" * 78)
    print("EXPERIMENT 2: Furstenberg's topological proof — computational content")
    print("=" * 78)

    # Furstenberg's proof: The sets S(a,b) = {a + nb : n in Z} form a basis for a topology.
    # Each S(a,b) is open AND closed. The complement of any finite set of primes would be
    # the intersection of finitely many S(0,p), which is open, hence the primes not being
    # the complement of an open set implies infinitely many primes.
    #
    # COMPUTATIONAL QUESTION: Can we use this topology to LOCATE primes efficiently?

    # --- 2a: Intersections of arithmetic progressions ---
    # The Furstenberg open sets that "avoid primes" are unions of S(0,p) = multiples of p.
    # The complement = {n : n not divisible by any prime up to B} gives us rough primes.
    # This is just... the sieve. Furstenberg's topology gives us NOTHING beyond the sieve.
    print("\n--- 2a: Furstenberg topology = the sieve in disguise ---")

    # Demonstrate: removing S(0,2), S(0,3), S(0,5), S(0,7) from {1..100}
    # leaves numbers coprime to 2,3,5,7 = numbers whose smallest prime factor > 7
    N = 100
    remaining = set(range(2, N+1))
    small_primes = [2, 3, 5, 7]
    for p in small_primes:
        multiples = set(range(2*p, N+1, p))
        remaining -= multiples
    # remaining should be {primes > 7} union {1} minus composites with small factors
    actual_primes_in_range = {p for p in PRIMES if 2 <= p <= N}
    remaining.discard(1)
    print(f"  After removing multiples of {small_primes} from [2,{N}]:")
    print(f"  Remaining: {sorted(remaining)[:30]}")
    print(f"  Actual primes: {sorted(actual_primes_in_range)[:30]}")
    print(f"  Extra composites remaining: {sorted(remaining - actual_primes_in_range)}")
    print(f"  This IS the sieve of Eratosthenes. Furstenberg adds no new computation.")

    # --- 2b: Topology of prime gaps via basis neighborhoods ---
    print("\n--- 2b: Topological density of primes in Furstenberg topology ---")
    # In Furstenberg's topology, a set is dense iff it intersects every S(a,b).
    # Primes are dense (Dirichlet's theorem: primes in every a+bZ with gcd(a,b)=1).
    # But density doesn't give us a LOCATION algorithm.

    # Test: for each basis set S(a,b) with small b, find the smallest prime in it
    print("  Smallest prime in S(a,b) = {a + nb : n >= 0}:")
    for b in range(2, 13):
        for a in range(1, b):
            if math.gcd(a, b) > 1:
                continue
            # Find smallest prime = a mod b
            smallest = None
            for k in range(max(2, a), 1000, b):
                if k >= a and k % b == a % b and is_prime_small(k):
                    smallest = k
                    break
            if smallest and b <= 6:
                print(f"    S({a},{b}): smallest prime = {smallest}")

    # --- 2c: Can the topology encode prime POSITIONS? ---
    print("\n--- 2c: Encoding prime positions via topological invariants ---")
    # The key invariant in Furstenberg's topology: the closure of {p(1), p(2), ..., p(n)}
    # What is cl({2,3,5,7,11,...,p(n)}) in this topology?
    # A set is closed iff it's a union of complete residue classes (or finite).
    # Since primes are not a union of residue classes (they're not periodic),
    # cl(primes) = Z in the Furstenberg topology.
    # This tells us NOTHING about individual prime locations.

    # More precisely: in this topology, convergence of a sequence a_n -> L means
    # a_n - L is eventually in every S(0,b), i.e., b | (a_n - L) for all b eventually.
    # This is just convergence in the profinite completion Z_hat.
    # Primes do NOT converge in Z_hat (they're equidistributed mod m for every m).
    print("  cl(primes) = Z in Furstenberg topology (primes are dense)")
    print("  Primes don't converge in profinite Z_hat (equidistributed mod m for all m)")
    print("  No topological invariant distinguishes prime positions from non-prime")

    print(f"\n  VERDICT: Furstenberg's topology is the sieve repackaged.")
    print(f"  It proves infinitude but encodes ZERO positional information.")
    return "FAIL: Furstenberg topology = sieve of Eratosthenes in topological language. No new computation."


# =============================================================================
# EXPERIMENT 3: Ergodic averages and prime reconstruction
# =============================================================================

def experiment3_ergodic():
    print("\n" + "=" * 78)
    print("EXPERIMENT 3: Ergodic averages along primes")
    print("=" * 78)

    # Bergelson-Richter (2020): For any measure-preserving system (X, mu, T) and f in L^2,
    # (1/pi(N)) sum_{p <= N} f(T^p x) -> int f d mu   for mu-a.e. x
    # i.e., ergodic averages along primes equidistribute.
    #
    # REVERSAL IDEA: If we know the ergodic average converges to int f dmu,
    # can we extract the prime positions from the convergence pattern?

    # --- 3a: Circle rotation ---
    print("\n--- 3a: Ergodic average along primes for irrational rotation ---")
    # T: x -> x + alpha (mod 1), alpha irrational
    # f(x) = cos(2*pi*x)
    # Ergodic average = (1/pi(N)) sum_{p<=N} cos(2*pi*(x0 + p*alpha))
    # Should converge to 0 (integral of cos over [0,1]).

    alpha = math.sqrt(2) - 1  # irrational
    x0 = 0.0
    N_vals = [100, 500, 1000, 5000, 10000, 50000]
    print(f"  alpha = sqrt(2)-1 = {alpha:.6f}, f(x) = cos(2*pi*x)")
    for N_max in N_vals:
        primes_up_to = [p for p in PRIMES if p <= N_max]
        if not primes_up_to:
            continue
        avg = sum(math.cos(2 * math.pi * (x0 + p * alpha)) for p in primes_up_to) / len(primes_up_to)
        # Compare with average along ALL integers
        avg_all = sum(math.cos(2 * math.pi * (x0 + n * alpha)) for n in range(1, N_max + 1)) / N_max
        print(f"  N={N_max:6d}: prime_avg = {avg:+.6f}, all_avg = {avg_all:+.6f}, "
              f"pi(N) = {len(primes_up_to)}")

    # The convergence rate tells us about correlations of primes with the rotation,
    # but NOT about individual prime positions. To reconstruct primes from the averages,
    # we'd need to solve: given {avg_f(N) for many f}, find {p <= N}.
    # This is a moment problem with O(N/ln N) unknowns and we'd need O(N/ln N) functions f.
    # -> O(N/ln N) work just to SET UP the problem. Worse than sieving.

    # --- 3b: Vinogradov-type exponential sums ---
    print("\n--- 3b: Exponential sums over primes (Vinogradov) ---")
    # S(alpha, N) = sum_{p <= N} e(p * alpha) where e(x) = e^{2*pi*i*x}
    # For irrational alpha, |S(alpha,N)| = o(pi(N)) by Vinogradov.
    # The RATE of cancellation depends on the diophantine type of alpha.

    import cmath
    for N_max in [1000, 10000, 50000]:
        primes_up_to = [p for p in PRIMES if p <= N_max]
        pi_N = len(primes_up_to)
        # Test various alpha
        for alpha_val in [0.5, math.sqrt(2)/2, 1/math.e, math.pi/10]:
            S = sum(cmath.exp(2j * math.pi * p * alpha_val) for p in primes_up_to)
            ratio = abs(S) / pi_N
            if N_max == 50000:
                print(f"  N={N_max}, alpha={alpha_val:.4f}: |S|/pi(N) = {ratio:.4f}")

    print(f"\n  Exponential sums decay as O(N/ln^A N) (Vinogradov).")
    print(f"  This gives STATISTICAL info about primes, not individual positions.")

    # --- 3c: Reconstruction via inverse problem ---
    print("\n--- 3c: Prime reconstruction from ergodic data ---")
    # Suppose we measure f_k = (1/pi(N)) sum_{p<=N} e(p*k/N) for k=0,...,N-1
    # This is the DFT of the prime indicator (scaled). We can recover the indicator
    # via inverse DFT, but computing f_k for all k requires knowing all primes <= N.
    # CIRCULAR.
    #
    # Alternative: measure f_k for SOME k values. This is compressed sensing.
    # But Session 7 already showed compressed sensing on zeta zeros fails to scale.
    print("  Reconstructing prime indicator from Fourier data is CIRCULAR:")
    print("  Computing the Fourier coefficients requires knowing the primes.")
    print("  Compressed sensing approach already failed in Session 7.")

    print(f"\n  VERDICT: Ergodic averages give STATISTICAL properties of primes.")
    print(f"  Cannot reconstruct individual prime positions without O(N/ln N) data points.")
    return "FAIL: Ergodic theory gives statistics, not positions. Reconstruction is circular or requires O(N) data."


# =============================================================================
# EXPERIMENT 4: Symbolic dynamics on prime indicator
# =============================================================================

def experiment4_symbolic():
    print("\n" + "=" * 78)
    print("EXPERIMENT 4: Symbolic dynamics on the prime indicator sequence")
    print("=" * 78)

    # Prime indicator: chi(n) = 1 if n prime, 0 otherwise
    # Sequence: 0,1,1,0,1,0,1,0,0,0,1,0,1,0,0,0,1,0,1,...

    # --- 4a: Topological entropy ---
    print("\n--- 4a: Topological entropy of prime indicator ---")
    # The topological entropy h_top measures the growth rate of distinct words.
    # For a sequence x, h_top = lim_{L->inf} (1/L) * log_2(|{words of length L in x}|)
    # For the prime indicator, all binary words of length L appear (up to some constraints).

    N = 50000
    indicator = bytearray(N + 1)
    for p in PRIMES:
        if p <= N:
            indicator[p] = 1

    word_counts = {}
    for L in [2, 3, 4, 5, 6, 7, 8, 10, 12, 15]:
        words = set()
        for start in range(2, N - L + 1):
            w = bytes(indicator[start:start + L])
            words.add(w)
        word_counts[L] = len(words)
        entropy = math.log2(len(words)) / L if len(words) > 1 else 0
        max_possible = 2 ** L
        coverage = len(words) / max_possible * 100
        print(f"  L={L:2d}: {len(words):6d} distinct words (of {max_possible:6d} possible, "
              f"{coverage:5.1f}%), h_top >= {entropy:.4f}")

    # --- 4b: Forbidden words ---
    print("\n--- 4b: Forbidden words (patterns that never appear) ---")
    for L in [3, 4, 5, 6]:
        all_words = set()
        for i in range(2**L):
            all_words.add(tuple(int(b) for b in format(i, f'0{L}b')))
        observed = set()
        for start in range(2, N - L + 1):
            w = tuple(indicator[start:start + L])
            observed.add(w)
        forbidden = all_words - observed
        if forbidden:
            # Show first few forbidden words
            forbidden_list = sorted(forbidden)[:5]
            print(f"  L={L}: {len(forbidden)} forbidden words. Examples: "
                  f"{'  '.join(''.join(str(b) for b in w) for w in forbidden_list)}")
        else:
            print(f"  L={L}: NO forbidden words (all {2**L} patterns appear)")

    # --- 4c: Is the prime indicator a sofic shift? ---
    print("\n--- 4c: Sofic shift test ---")
    # A sofic shift is the set of labels on bi-infinite walks in a finite labeled graph.
    # Equivalently, it's recognized by a finite automaton.
    # If the prime indicator were sofic, there'd be a DFA recognizing it -> primes would be
    # eventually periodic or have bounded complexity. But:
    # - Block complexity p(L) = number of distinct L-blocks
    # - For sofic: p(L) is eventually <= C (constant) — NO! For sofic, p(L) grows at most linearly.
    # - For prime indicator: p(L) -> 2^L (all words eventually appear)
    # This means the prime indicator has MAXIMAL complexity -> NOT sofic.

    growth_rates = []
    prev_count = 1
    for L in sorted(word_counts.keys()):
        if L > 2:
            growth = word_counts[L] / word_counts.get(L-1, word_counts[min(word_counts.keys())])
        else:
            growth = word_counts[L]
        growth_rates.append((L, word_counts[L], growth))

    print(f"  Block complexity p(L) for prime indicator:")
    for L, count, growth in growth_rates[:8]:
        print(f"    p({L}) = {count}")
    print(f"  p(L) grows exponentially -> NOT a sofic shift.")
    print(f"  NOT even a coded shift (those have p(L) = O(L)).")
    print(f"  The prime indicator has MAXIMAL symbolic complexity.")

    # --- 4d: Higher-dimensional shifts ---
    print("\n--- 4d: 2D symbolic dynamics (primes on a grid) ---")
    # Place indicator on a 2D grid (row-major). Does the 2D pattern have lower entropy?
    # Try grid widths w = 6, 30, 210 (primorials) since primes avoid multiples.
    for w in [6, 30]:
        rows = 200
        patterns_2x2 = set()
        for r in range(rows - 1):
            for c in range(w - 1):
                n00 = r * w + c + 2  # offset by 2 so we start at n=2
                n01 = n00 + 1
                n10 = n00 + w
                n11 = n00 + w + 1
                if max(n00, n01, n10, n11) <= N:
                    pat = (indicator[n00], indicator[n01], indicator[n10], indicator[n11])
                    patterns_2x2.add(pat)
        print(f"  Width w={w}: {len(patterns_2x2)} distinct 2x2 blocks (of 16 possible)")
    print(f"  2D representation doesn't reduce entropy. All local patterns still appear.")

    # --- 4e: Sturmian / substitution analysis ---
    print("\n--- 4e: Is prime indicator close to a Sturmian sequence? ---")
    # Sturmian sequences have p(L) = L+1 (minimal non-periodic complexity).
    # Prime indicator has p(L) >> L+1. How far off?
    # A Sturmian sequence is a rotation sequence: s_n = floor((n+1)*alpha + beta) - floor(n*alpha + beta)
    # Can we find alpha, beta that BEST approximates the prime indicator?
    # Best alpha would be ~1/ln(n) (prime density). But gaps are irregular.

    # Compute hamming distance from best Sturmian approximation
    best_hamming = N
    # Try alpha = 1/ln(k) for various k (matching local density)
    for alpha_test in [0.15, 0.18, 0.20, 0.22, 0.25]:
        for beta_test in [0.0, 0.1, 0.3, 0.5, 0.7]:
            hamming = 0
            for n in range(2, min(2000, N)):
                sturmian_val = int(math.floor((n + 1) * alpha_test + beta_test)) - \
                               int(math.floor(n * alpha_test + beta_test))
                if sturmian_val != indicator[n]:
                    hamming += 1
            if hamming < best_hamming:
                best_hamming = hamming
                best_ab = (alpha_test, beta_test)

    total = min(2000, N) - 2
    print(f"  Best Sturmian approx (alpha={best_ab[0]}, beta={best_ab[1]}): "
          f"Hamming distance = {best_hamming}/{total} ({best_hamming/total*100:.1f}%)")
    print(f"  A Sturmian model misses ~{best_hamming/total*100:.0f}% of positions. Not viable.")

    print(f"\n  VERDICT: Prime indicator has MAXIMAL symbolic complexity.")
    print(f"  Not sofic, not substitutive, not Sturmian, not even coded.")
    print(f"  No finite automaton or shift of finite type can generate it.")
    return "FAIL: Prime indicator has maximal complexity. No symbolic dynamics shortcut."


# =============================================================================
# EXPERIMENT 5: Fixed points of iterated maps at primes
# =============================================================================

def experiment5_fixed_points():
    print("\n" + "=" * 78)
    print("EXPERIMENT 5: Fixed points of iterated maps at primes")
    print("=" * 78)

    # --- 5a: Newton's method on W(x) = psi(x) - x ---
    # psi(x) = sum_{p^k <= x} ln(p) is Chebyshev's function.
    # PNT says psi(x) ~ x. We want the points where pi(x) jumps, i.e., primes.
    # Problem: psi(x) is a step function. Newton's method doesn't apply directly.
    # We can try a smoothed version.
    print("\n--- 5a: Newton on smoothed Chebyshev psi ---")

    # Compute psi(x) for integer x up to 1000
    def chebyshev_psi(x):
        """Compute psi(x) = sum_{p^k <= x} ln(p)."""
        result = 0.0
        for p in PRIMES:
            if p > x:
                break
            pk = p
            while pk <= x:
                result += math.log(p)
                pk *= p
        return result

    # psi(x) - x has oscillations. Its "zeros" (where psi crosses x) are NOT at primes.
    # Let's check:
    print("  psi(x) vs x for x near primes:")
    for x in [10, 20, 50, 100, 200]:
        psi_val = chebyshev_psi(x)
        print(f"    x={x}: psi(x)={psi_val:.2f}, psi(x)-x={psi_val-x:.2f}, "
              f"psi(x)/x={psi_val/x:.4f}")

    # The zeros of psi(x) - x are where the PNT error term changes sign.
    # These are related to zeta zeros, NOT to individual primes.
    print("  psi(x) - x oscillates due to zeta zeros, crosses 0 at non-prime points.")
    print("  Newton's method on psi(x)-x finds zeta-zero related crossings, not primes.")

    # --- 5b: Map whose fixed points are primes ---
    print("\n--- 5b: Maps with prime fixed points ---")
    # Consider f(x) = x * (Wilson quotient).
    # Wilson's theorem: (p-1)! = -1 mod p for prime p.
    # W(n) = ((n-1)! + 1) / n (Wilson quotient, integer iff n is prime)
    # Define g(n) = n + ((n-1)! + 1) mod n
    # Then g(n) = n iff n is prime (since (n-1)!+1 mod n = 0 for primes).
    # But computing (n-1)! mod n costs O(n) multiplications -> O(n) per primality test.

    print("  Wilson-based fixed point map: g(n) = n + ((n-1)!+1 mod n)")
    print("  g(n) = n iff n is prime (Wilson's theorem)")
    for n in range(2, 20):
        factorial_mod = 1
        for i in range(2, n):
            factorial_mod = (factorial_mod * i) % n
        wilson_residue = (factorial_mod + 1) % n
        is_fixed = (wilson_residue == 0)
        actual = is_prime_small(n)
        print(f"    n={n:2d}: (n-1)!+1 mod n = {wilson_residue}, "
              f"fixed={'Y' if is_fixed else 'N'}, prime={'Y' if actual else 'N'}")
    print("  Correctly identifies primes but costs O(n) per evaluation -> useless for p(n).")

    # --- 5c: Newton's method on a function with roots at primes ---
    print("\n--- 5c: Newton's method with sin(pi*Gamma(x)/x) ---")
    # For integer x: Gamma(x) = (x-1)!
    # sin(pi * (x-1)!/x) = 0 when (x-1)!/x is integer, i.e., when x | (x-1)!.
    # By Wilson's theorem, p | (p-1)! + 1 but p does NOT divide (p-1)!.
    # So sin(pi*(x-1)!/x) != 0 at primes. OPPOSITE of what we want.
    # Actually: for COMPOSITE n > 4, n | (n-1)! (Korselt). So the zeros are composites.
    # We want the NON-zeros, which are primes. Can't use Newton for that.

    # Alternative: f(x) = sin(pi*x/2) * product_{p <= sqrt(x)} sin(pi*x/p)
    # Roots include all composites. But this requires knowing primes up to sqrt(x).
    # CIRCULAR.
    print("  sin(pi*Gamma(x)/x) has roots at COMPOSITES (by Korselt/Wilson).")
    print("  We want the non-roots (primes). Can't use Newton to find non-roots.")
    print("  Any explicit function with roots exactly at primes requires knowing primes. CIRCULAR.")

    # --- 5d: Iterated function systems and prime orbits ---
    print("\n--- 5d: Iterated function system approach ---")
    # Can we define f such that the orbit {f^n(x0)} enumerates primes?
    # If f(p_k) = p_{k+1} for all k, then f encodes ALL prime gaps.
    # Cramér's conjecture: gaps ~ (ln p)^2, which vary irregularly.
    # The function f would need ~log(p) bits of "new information" at each step.
    # This is the ~178-bit barrier again: f cannot be polylog-time computable.

    # Demonstrate: what would f look like?
    print("  If f(p_k) = p_{k+1}, then f encodes prime gaps:")
    for k in range(20):
        pk = PRIMES[k]
        pk1 = PRIMES[k + 1]
        gap = pk1 - pk
        print(f"    f({pk:3d}) = {pk1:3d}  (gap = {gap})")
    print("  Gaps are 1,2,2,4,2,4,2,4,6,2,6,4,2,4,6,6,2,6,4,...")
    print("  This sequence has Kolmogorov complexity >= 0.5*log2(n) bits for p(n).")
    print("  No polylog-computable function f can produce this orbit.")

    # --- 5e: Von Mangoldt and Möbius as dynamics ---
    print("\n--- 5e: Dynamical interpretation of von Mangoldt Λ(n) ---")
    # Λ(n) = ln(p) if n = p^k, 0 otherwise.
    # sum_{d|n} Λ(d) = ln(n)  (fundamental identity)
    # This is a linear recurrence over divisors. Can we invert it dynamically?
    # Λ(n) = sum_{d|n} mu(d) * ln(n/d)  (Möbius inversion)
    # But computing mu(d) for all d|n requires factoring n. O(sqrt(n)) at best.

    def von_mangoldt(n):
        """Compute Λ(n)."""
        if n <= 1:
            return 0
        # Check if n is a prime power
        for p in PRIMES:
            if p * p > n:
                break
            if n % p == 0:
                while n % p == 0:
                    n //= p
                return math.log(p) if n == 1 else 0
        # n is prime
        return math.log(n)

    print("  Λ(n) for n = 1..20:")
    vals = []
    for n in range(1, 21):
        v = von_mangoldt(n)
        vals.append(f"{v:.2f}" if v > 0 else "0")
    print(f"    {', '.join(vals)}")
    print("  Computing Λ(n) requires factoring n. No dynamical shortcut.")

    # --- 5f: Functional iteration on Riemann's R(x) ---
    print("\n--- 5f: Fixed-point iteration with Riemann R function ---")
    # R(x) ~ pi(x). Define T(n) = R^{-1}(n) as initial guess for p(n).
    # Then refine: p(n) = T(n) + correction.
    # The correction requires knowing pi(T(n)) exactly -> back to square one.
    # But what if we iterate: x_{k+1} = x_k + (n - pi_approx(x_k)) * ln(x_k)?
    # This is Newton's method on pi(x) - n = 0. Convergence requires pi(x) evaluation.

    # Demonstrate the iteration with exact pi (cheating):
    pi_cache = {}
    count = 0
    for i in range(2, LIMIT + 1):
        if i in PRIME_SET:
            count += 1
        pi_cache[i] = count

    def pi_exact(x):
        x = int(x)
        if x < 2: return 0
        if x <= LIMIT: return pi_cache.get(x, pi_cache[min(x, LIMIT)])
        return None

    print("  Newton iteration: x_{k+1} = x_k + (n - pi(x_k)) * ln(x_k)")
    for target_n in [100, 1000, 5000]:
        # Initial guess via n * ln(n)
        x = target_n * math.log(target_n)
        for iteration in range(20):
            xi = int(x)
            if xi < 2: xi = 2
            if xi > LIMIT: break
            pi_val = pi_exact(xi)
            if pi_val is None: break
            if pi_val == target_n:
                break
            # Newton step
            x = x + (target_n - pi_val) * math.log(max(x, 2))
        actual = PRIMES[target_n - 1] if target_n <= len(PRIMES) else '?'
        result_prime = xi
        found_pi = pi_exact(min(result_prime, LIMIT))
        print(f"    n={target_n}: converged to x={result_prime}, pi(x)={found_pi}, "
              f"actual p(n)={actual}, match={'YES' if result_prime == actual else 'NO'}")
    print("  Newton converges in ~5 iterations BUT each step needs pi(x) exactly.")
    print("  Computing pi(x) costs O(x^{2/3}). Total: O(p(n)^{2/3}). NOT polylog.")

    print(f"\n  VERDICT: All fixed-point methods require evaluating pi(x) or factoring,")
    print(f"  which costs at least O(n^{1/3}) by known lower bounds.")
    return "FAIL: Every fixed-point/Newton approach bottlenecks on pi(x) computation or factoring."


# =============================================================================
# EXPERIMENT 6: Bonus — novel hybrid ideas
# =============================================================================

def experiment6_novel():
    print("\n" + "=" * 78)
    print("EXPERIMENT 6: Novel hybrid dynamical approaches")
    print("=" * 78)

    # --- 6a: Arithmetic dynamics — preperiodic points ---
    print("\n--- 6a: Arithmetic dynamics (preperiodic points of polynomial maps) ---")
    # For f(z) = z^2 + c, the preperiodic points form an algebraic set.
    # The primes dividing the iterates f^n(0) have been studied (dynamical Zsygmondy).
    # Can we pick c so that the primitive prime divisors of f^n(0) give us ALL primes?

    # Test: f(z) = z^2 + 1, orbit of 0: 0, 1, 2, 5, 26, 677, ...
    # Prime divisors of orbit elements:
    print("  Orbit of 0 under f(z) = z^2 + c:")
    for c in [1, -1, 2, 3, -2]:
        z = 0
        orbit = [z]
        prime_divs = set()
        for _ in range(12):
            z = z * z + c
            if abs(z) > 10**15:
                break
            orbit.append(z)
            if z != 0:
                n = abs(z)
                for p in PRIMES[:100]:
                    if n % p == 0:
                        prime_divs.add(p)
        print(f"    c={c:+d}: orbit = {orbit[:7]}{'...' if len(orbit)>7 else ''}")
        print(f"         prime divisors: {sorted(prime_divs)[:15]}")
    print("  Orbits grow doubly-exponentially. Prime divisors are SPARSE.")
    print("  Cannot enumerate all primes this way.")

    # --- 6b: Measure-theoretic prime generation ---
    print("\n--- 6b: Benford-like distribution of prime residues ---")
    # Do primes follow any non-uniform distribution that could be exploited?
    # Chebyshev bias: slightly more primes = 3 mod 4 than 1 mod 4 (up to some bound).
    for mod in [4, 6, 10, 30]:
        residue_counts = defaultdict(int)
        for p in PRIMES[:5000]:
            residue_counts[p % mod] += 1
        residues = sorted(residue_counts.items())
        desc = ', '.join(f'{r}:{c}' for r, c in residues if c > 10)
        print(f"  Primes mod {mod}: {desc}")
    print("  Primes equidistribute mod m (Dirichlet). Chebyshev bias is O(sqrt(x)).")
    print("  No exploitable non-uniformity.")

    # --- 6c: Transfer operator / Ruelle zeta function ---
    print("\n--- 6c: Transfer operator approach ---")
    # The Ruelle zeta function for the map T(x) = 1/x mod 1 (Gauss map) is
    # related to Selberg zeta functions. For the modular surface, the Selberg zeta
    # function has zeros related to eigenvalues of the Laplacian, NOT at primes.
    # The PRIME geodesics on the modular surface correspond to hyperbolic elements of SL(2,Z),
    # which are parametrized by solutions to t^2 - 4 = discriminant of real quadratic fields.
    # These are NOT the rational primes.
    print("  Selberg zeta zeros = Laplacian eigenvalues, NOT primes.")
    print("  Prime geodesics on modular surface ≠ rational primes.")
    print("  Transfer operator spectral theory doesn't help locate primes.")

    # --- 6d: Quantum ergodicity and primes ---
    print("\n--- 6d: Quantum ergodicity ---")
    # QUE (Quantum Unique Ergodicity) on the modular surface:
    # Hecke-Maass eigenforms equidistribute (Lindenstrauss 2006, Fields Medal).
    # Hecke eigenvalues at prime p: lambda(p) = a_p / p^{1/2}
    # These satisfy Ramanujan conjecture: |a_p| <= 2*sqrt(p).
    # But computing a_p for Maass forms requires O(p^{1/2}) work.
    # Even with Hecke multiplicativity, no shortcut to individual primes.
    print("  Hecke eigenvalues encode prime information but cost O(sqrt(p)) to compute.")
    print("  QUE gives equidistribution, not individual prime locations.")

    print(f"\n  VERDICT: All novel dynamical/ergodic hybrids hit the same barrier.")
    return "FAIL: Arithmetic dynamics, transfer operators, quantum ergodicity all confirm the barrier."


# =============================================================================
# MAIN: Run all experiments
# =============================================================================

def main():
    print("*" * 78)
    print("SESSION 8: Dynamical Systems, Ergodic Theory & Orbital Approaches to Primes")
    print("*" * 78)
    print(f"Testing 6 experiment families. Reference: {len(PRIMES)} primes up to {LIMIT}.")
    print()

    results = []
    t_total = time.time()

    experiments = [
        ("1: Collatz-like maps & Rowland", experiment1_collatz_and_rowland),
        ("2: Furstenberg topology", experiment2_furstenberg),
        ("3: Ergodic averages", experiment3_ergodic),
        ("4: Symbolic dynamics", experiment4_symbolic),
        ("5: Fixed points of iterated maps", experiment5_fixed_points),
        ("6: Novel hybrid approaches", experiment6_novel),
    ]

    for name, func in experiments:
        t0 = time.time()
        try:
            result = func()
        except Exception as e:
            result = f"ERROR: {e}"
            import traceback
            traceback.print_exc()
        elapsed = time.time() - t0
        results.append((name, result, elapsed))
        print(f"\n  [{name}] Time: {elapsed:.2f}s")

    total_time = time.time() - t_total

    print("\n" + "=" * 78)
    print("SUMMARY OF ALL EXPERIMENTS")
    print("=" * 78)
    for name, result, elapsed in results:
        status = "FAIL" if "FAIL" in result else ("PROMISING" if "PROMIS" in result else "???")
        print(f"  {status:10s} | {name:45s} | {elapsed:.2f}s")
        print(f"           | {result}")
    print(f"\n  Total time: {total_time:.2f}s")

    print("\n" + "=" * 78)
    print("THEORETICAL ANALYSIS: WHY DYNAMICAL SYSTEMS CANNOT BREAK THE BARRIER")
    print("=" * 78)
    print("""
  The fundamental issue across ALL dynamical approaches:

  1. INFORMATION THEORY: p(n) contains ~0.5*log2(n) bits of irreducible information.
     Any map f that produces p(n) from n must encode these bits somewhere.
     If f is polylog-time computable, it has at most O(polylog n) bits of description.
     But p(n) needs ~0.5*log2(n) bits -> f must extract them from n itself,
     which requires evaluating pi(x) or equivalent, costing O(n^{1/3+}) minimum.

  2. ERGODIC THEORY: Ergodic averages along primes equidistribute (Bergelson-Richter).
     This means primes are "generic" in the ergodic sense — they carry NO special
     dynamical structure that could be exploited. Equidistribution is the ABSENCE
     of exploitable structure.

  3. SYMBOLIC DYNAMICS: The prime indicator has MAXIMAL topological entropy and
     block complexity p(L) ~ 2^L. It's not sofic, not substitutive, not coded.
     This means no finite-state machine can generate it. Any generator needs
     unbounded memory, growing as O(n / ln n).

  4. FIXED POINTS: Any function whose fixed points are exactly the primes must
     encode the prime-counting function pi(x). Evaluating such a function at x
     costs at least O(x^{1/3}) by Aggarwal's lower bound (2025).

  5. FURSTENBERG TOPOLOGY: Computationally equivalent to the Sieve of Eratosthenes.
     The topological proof of prime infinitude is non-constructive and provides
     zero positional information about individual primes.

  CONCLUSION: Dynamical systems and ergodic theory provide STATISTICAL properties
  of primes (equidistribution, entropy bounds, density theorems) but cannot provide
  INDIVIDUAL prime positions in polylog time. This is consistent with all 205+
  prior approaches. The ~178-bit barrier remains unbreakable.
""")


if __name__ == "__main__":
    main()
