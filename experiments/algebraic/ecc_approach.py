"""
Session 6: Error-Correcting Code Perspective on Primes

RADICAL IDEA: What if we treat the approximate formula (Cipolla/R^{-1})
as a "noisy channel" and the exact prime sequence as the "codeword"?

The "noise" is δ(n) = p(n) - approx(n). If this noise has structure
that allows error correction, we can "decode" the exact prime.

Key question: Is there a SYNDROME that tells us the exact error?

Approach 1: Modular Syndrome Decoding
  For each small prime q, compute approx(n) mod q.
  Compare with the constraint p(n) mod q ≠ 0.
  Can these constraints uniquely determine p(n)?

Approach 2: Nearest Prime Decoding
  Given approx(n), find the k-th prime nearest to it, where k depends
  on systematic bias corrections.

Approach 3: Prime Constellation Matching
  Primes come in patterns (twin primes, triplets, etc.).
  Given the approximate location, which constellation does p(n) belong to?

Approach 4: Wheel Factorization + Position Recovery
  Primes > 5 are all ≡ 1,7,11,13,17,19,23,29 (mod 30).
  This gives us a "wheel" with 8 slots per revolution.
  Can we predict WHICH slot p(n) is in, and which revolution?
"""

import numpy as np
from sympy import prime, nextprime, prevprime, isprime, primepi
import time

def wheel_position(p, wheel_mod=30):
    """Return the position of p in the wheel mod wheel_mod."""
    residues = sorted([r for r in range(1, wheel_mod)
                       if all(r % q != 0 for q in [2, 3, 5])])
    r = p % wheel_mod
    if r in residues:
        return residues.index(r)
    return -1  # Not a valid wheel position (only for p=2,3,5)

def experiment_wheel_prediction():
    """
    Can we predict which wheel position p(n) occupies?

    For wheel mod 30, there are 8 valid residues: 1,7,11,13,17,19,23,29.
    If we could predict the residue, we'd reduce candidates by 8/30 = 73%.
    """
    print("="*70)
    print("EXPERIMENT: WHEEL POSITION PREDICTION")
    print("="*70)

    N = 10000
    wheel_mod = 30
    residues_30 = [1, 7, 11, 13, 17, 19, 23, 29]

    # Compute wheel positions for all primes
    primes_list = [int(prime(n)) for n in range(4, N + 4)]  # Skip 2, 3, 5
    positions = [wheel_position(p) for p in primes_list]

    # Distribution of wheel positions
    from collections import Counter
    pos_counts = Counter(positions)
    print("\nWheel position distribution (mod 30):")
    for i, r in enumerate(residues_30):
        count = pos_counts.get(i, 0)
        print(f"  pos {i} (≡{r:2d} mod 30): {count} ({100*count/len(positions):.1f}%)")

    # Can we predict position from n?
    # Build features
    positions = np.array(positions)
    ns = np.arange(4, N + 4, dtype=np.float64)

    # Simple prediction: use n mod 8 (since there are 8 positions)
    simple_pred = (ns.astype(int) * 3 + 1) % 8
    simple_acc = np.mean(simple_pred == positions)
    print(f"\nRandom baseline accuracy: {100/8:.1f}%")
    print(f"Simple (n*3+1 mod 8) accuracy: {100*simple_acc:.1f}%")

    # Consecutive position transitions
    transitions = np.zeros((8, 8), dtype=int)
    for i in range(len(positions) - 1):
        transitions[positions[i], positions[i+1]] += 1

    # Normalize to get transition probabilities
    trans_probs = transitions / transitions.sum(axis=1, keepdims=True)

    print("\nTransition matrix (rows=from, cols=to):")
    print("     " + "  ".join(f" {r:2d}" for r in residues_30))
    for i, r in enumerate(residues_30):
        row = " ".join(f"{trans_probs[i,j]:.2f}" for j in range(8))
        print(f"  {r:2d}: {row}")

    # Can transition matrix predict next position?
    markov_pred = np.array([np.argmax(trans_probs[positions[i]])
                            for i in range(len(positions)-1)])
    markov_acc = np.mean(markov_pred == positions[1:])
    print(f"\nMarkov chain prediction accuracy: {100*markov_acc:.1f}%")

    # What about using the GAP to predict position?
    gaps = np.diff(primes_list)
    # Gap mod 30 determines position change
    pos_changes = np.diff(positions) % 8
    gap_mod30 = gaps % 30

    # For each gap mod 30, what's the position change?
    print("\nGap mod 30 → position change:")
    for g in sorted(set(gap_mod30[:100])):
        mask = gap_mod30 == g
        if mask.sum() > 10:
            changes = pos_changes[mask[:len(pos_changes)]]
            most_common = Counter(changes).most_common(3)
            print(f"  gap��{g:2d}: {most_common}")

def experiment_nearest_prime_decoding():
    """
    Given R^{-1}(n) as approximate, how often is the nearest prime correct?
    And: can a systematic bias correction improve this?
    """
    print("\n" + "="*70)
    print("EXPERIMENT: NEAREST PRIME DECODING")
    print("="*70)

    N = 5000

    # Cipolla approximation
    correct_nearest = 0
    correct_biased = 0
    correct_k_nearest = {k: 0 for k in range(1, 11)}

    errors = []

    for n in range(2, N + 2):
        p_n = int(prime(n))
        ln_n = np.log(n)
        ln_ln_n = np.log(ln_n) if ln_n > 1 else 0.1

        # Cipolla approximation (3rd order)
        approx = n * (ln_n + ln_ln_n - 1 + (ln_ln_n - 2) / ln_n +
                      (ln_ln_n**2 - 6*ln_ln_n + 11) / (2 * ln_n**2))

        error = p_n - approx
        errors.append(error)

        # Nearest prime to approx
        nearest = int(round(approx))
        if nearest < 2:
            nearest = 2
        np_above = int(nextprime(nearest - 1))
        np_below = int(prevprime(nearest + 1)) if nearest > 2 else 2

        if abs(np_above - approx) < abs(np_below - approx):
            nearest_prime = np_above
        else:
            nearest_prime = np_below

        if nearest_prime == p_n:
            correct_nearest += 1

        # Count how many primes away the true p(n) is from nearest_prime
        # Find primes near approx
        candidates = []
        p = max(2, int(approx) - 100)
        if not isprime(p):
            p = int(nextprime(p))
        while p < approx + 100:
            candidates.append(p)
            p = int(nextprime(p))

        if p_n in candidates:
            idx = candidates.index(p_n)
            # Find which k-nearest it is
            dists = [(abs(c - approx), i) for i, c in enumerate(candidates)]
            dists.sort()
            for k_rank, (d, orig_idx) in enumerate(dists, 1):
                if candidates[orig_idx] == p_n:
                    for k in range(k_rank, 11):
                        correct_k_nearest[k] += 1
                    break

    print(f"\nResults for n=2..{N+1}:")
    print(f"  Nearest prime = p(n): {correct_nearest}/{N} ({100*correct_nearest/N:.1f}%)")
    for k in [1, 2, 3, 5, 10]:
        print(f"  Within {k}-nearest primes: {correct_k_nearest[k]}/{N} ({100*correct_k_nearest[k]/N:.1f}%)")

    # Study the rank (which k-nearest) distribution
    errors = np.array(errors)
    print(f"\nApproximation error statistics:")
    print(f"  Mean: {np.mean(errors):.2f}")
    print(f"  Std: {np.std(errors):.2f}")
    print(f"  Mean/average_gap: {np.mean(errors)/np.mean(np.log(N*np.log(N))):.2f}")

    # Can we correct the bias?
    # The mean error grows with n. Fit it.
    ns = np.arange(2, N + 2, dtype=np.float64)
    # Fit: bias(n) ≈ a * ln(n)^2 / n + b * ln(ln(n))^3 / ln(n) + ...
    ln_ns = np.log(ns)
    ln_ln_ns = np.log(ln_ns)

    X_bias = np.column_stack([
        ln_ln_ns**3 / ln_ns,
        ln_ln_ns**2 / ln_ns,
        ln_ln_ns / ln_ns,
        1 / ln_ns,
        np.ones(N),
    ])

    bias_coeffs, _, _, _ = np.linalg.lstsq(X_bias, errors, rcond=None)
    predicted_bias = X_bias @ bias_coeffs

    corrected_errors = errors - predicted_bias
    print(f"\nAfter bias correction:")
    print(f"  RMS error: {np.sqrt(np.mean(corrected_errors**2)):.2f}")
    correct_after_bias = 0
    for i, n in enumerate(range(2, N + 2)):
        p_n = int(prime(n))
        approx_corrected = round(n * (np.log(n) + np.log(np.log(n)) - 1 +
                                      (np.log(np.log(n)) - 2) / np.log(n) +
                                      (np.log(np.log(n))**2 - 6*np.log(np.log(n)) + 11) / (2 * np.log(n)**2))
                                 + predicted_bias[i])
        nearest = int(nextprime(approx_corrected - 1)) if approx_corrected > 2 else 2
        if nearest == p_n:
            correct_after_bias += 1
    print(f"  Nearest prime after correction: {correct_after_bias}/{N} ({100*correct_after_bias/N:.1f}%)")

def experiment_prime_gap_distribution():
    """
    Study whether the gap distribution at a given scale can help decode.

    Key insight: If we know the DISTRIBUTION of gaps near p(n),
    we can compute the probability of each candidate being p(n).
    The candidate with highest probability is our guess.
    """
    print("\n" + "="*70)
    print("EXPERIMENT: PROBABILISTIC GAP DECODING")
    print("="*70)

    N = 10000
    primes_list = [int(prime(n)) for n in range(1, N + 2)]
    gaps = [primes_list[i+1] - primes_list[i] for i in range(len(primes_list)-1)]

    # Gap distribution
    from collections import Counter
    gap_counts = Counter(gaps)

    print("Most common gaps:")
    for gap, count in gap_counts.most_common(10):
        print(f"  gap={gap}: {count} ({100*count/len(gaps):.1f}%)")

    # Cramér model: gaps ~ Exponential(ln(p))
    # Test goodness-of-fit
    for n_start in [100, 1000, 5000]:
        local_gaps = gaps[n_start:n_start+500]
        local_primes = primes_list[n_start:n_start+500]
        expected_gap = np.mean([np.log(p) for p in local_primes])
        actual_mean = np.mean(local_gaps)
        actual_std = np.std(local_gaps)

        # For exponential, std = mean
        print(f"\n  At n≈{n_start}: E[gap]={expected_gap:.2f}, "
              f"actual mean={actual_mean:.2f}, std={actual_std:.2f}")
        print(f"    std/mean = {actual_std/actual_mean:.3f} "
              f"(exponential would give 1.0)")

    # Key question: given that we know the approximate location,
    # what's the expected number of primes we need to check?
    # If error is ε and average gap is g, we check about ε/g candidates.

    # For Cipolla at n=10000: error std ≈ 50, gap ≈ ln(10^5) ≈ 11.5
    # So ~50/11.5 ≈ 4.3 candidates on average
    print("\n\nExpected candidates to check (error_std / avg_gap):")
    for n in [100, 1000, 10000]:
        p = int(prime(n))
        ln_n = np.log(n)
        ln_ln_n = np.log(ln_n)
        approx = n * (ln_n + ln_ln_n - 1 + (ln_ln_n - 2) / ln_n)

        # Compute error std from nearby points
        local_errors = []
        for k in range(max(2, n-50), n+50):
            pk = int(prime(k))
            lk = np.log(k)
            llk = np.log(lk)
            a = k * (lk + llk - 1 + (llk - 2) / lk)
            local_errors.append(pk - a)

        err_std = np.std(local_errors)
        avg_gap = np.log(p)
        candidates = err_std / avg_gap

        print(f"  n={n}: p(n)≈{p}, error_std≈{err_std:.1f}, "
              f"avg_gap≈{avg_gap:.1f}, candidates≈{candidates:.1f}")

def experiment_ensemble_voting():
    """
    Use MULTIPLE independent approximations and take a "vote"
    on which prime is correct.

    Each approximation gives a different estimate. The true prime
    should be near ALL of them. The intersection of confidence
    intervals might pin down the exact prime.
    """
    print("\n" + "="*70)
    print("EXPERIMENT: ENSEMBLE VOTING")
    print("="*70)

    N = 2000
    correct = 0

    for n in range(10, N + 10):
        p_n = int(prime(n))
        ln_n = np.log(n)
        ln_ln_n = np.log(ln_n)

        # Multiple approximations
        approx1 = n * (ln_n + ln_ln_n - 1)  # Cipolla order 1
        approx2 = n * (ln_n + ln_ln_n - 1 + (ln_ln_n - 2) / ln_n)  # Order 2
        approx3 = n * (ln_n + ln_ln_n - 1 + (ln_ln_n - 2) / ln_n +
                       (ln_ln_n**2 - 6*ln_ln_n + 11) / (2 * ln_n**2))  # Order 3

        # Rosser's bounds: n*(ln(n) + ln(ln(n)) - 1) < p(n) < n*(ln(n) + ln(ln(n)))
        # for n >= 6
        upper = n * (ln_n + ln_ln_n) if n >= 6 else p_n + 10
        lower = n * (ln_n + ln_ln_n - 1) if n >= 6 else max(2, p_n - 10)

        # Average estimate
        avg_approx = (approx1 + approx2 + approx3) / 3

        # Weighted average (weight by expected accuracy)
        # Order 3 is best, give it highest weight
        weighted = 0.1 * approx1 + 0.3 * approx2 + 0.6 * approx3

        # Find primes in [lower, upper]
        candidates = []
        p = max(2, int(lower))
        if not isprime(p):
            p = int(nextprime(p))
        while p <= upper + 10:  # slight margin
            candidates.append(p)
            p = int(nextprime(p))

        # Vote: which candidate is closest to the weighted average?
        if candidates:
            best = min(candidates, key=lambda c: abs(c - weighted))
            if best == p_n:
                correct += 1

    print(f"\nEnsemble voting accuracy (n=10..{N+9}): {correct}/{N} ({100*correct/N:.1f}%)")

def main():
    print("Session 6: Error-Correcting Code Perspective")
    print(f"Started: {time.strftime('%H:%M:%S')}\n")

    t0 = time.time()

    experiment_wheel_prediction()
    experiment_nearest_prime_decoding()
    experiment_prime_gap_distribution()
    experiment_ensemble_voting()

    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"Total time: {elapsed:.1f}s")

if __name__ == '__main__':
    main()
