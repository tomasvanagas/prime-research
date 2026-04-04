"""
Session 4: PRIME GAP ENCODING / PREDICTION

THE LAST CREATIVE IDEA:

p(n) = 2 + Σ_{k=1}^{n-1} g(k)  where g(k) = p(k+1) - p(k) is the k-th gap.

If we could PREDICT g(k) exactly, we'd have p(n) without sieving.
g(k) is always even (for k ≥ 2) and g(k) ~ ln(p(k)) on average.

What if gaps have a COMPRESSIBLE pattern?

Cramér's model: g(k) ≈ Poisson(ln p(k))
But actual gaps have more structure:
  - Gaps are always even (for p > 2)
  - Gaps avoid certain residues mod 6 (gaps ≡ 0 mod 6 are common)
  - The Hardy-Littlewood conjecture predicts exact density of each gap size
  - Consecutive gaps are weakly correlated (Lemke Oliver & Soundararajan, 2016)

IDEA: What if we could express g(n) as:
  g(n) = 2 * round(f(n)) where f(n) is a smooth function?

This would give p(n) = 2 + 2 * Σ round(f(k))
If f is computable in O(polylog), and the rounding is always correct...

But the Prime Number Theorem already tells us f(n) ≈ ln(n * ln(n)) / 2.
The ROUNDING is the hard part. The fluctuations in g(n) around its mean
are the same random walk as δ(n) = p(n) - R^{-1}(n).

Let me test: how predictable are gaps from local information?
"""

import math
import time

def sieve(n):
    s = [True] * (n+1)
    s[0] = s[1] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]:
            for j in range(i*i, n+1, i):
                s[j] = False
    return [i for i, v in enumerate(s) if v]

primes = sieve(500000)
gaps = [primes[i+1] - primes[i] for i in range(len(primes)-1)]

print("=" * 70)
print("PRIME GAP PREDICTION ANALYSIS")
print("=" * 70)

# Test 1: Can consecutive gaps predict the next gap?
print("\n--- Gap autocorrelation ---")
N = len(gaps)
mean_g = sum(gaps) / N

for lag in range(1, 11):
    corr_num = sum((gaps[i] - mean_g) * (gaps[i+lag] - mean_g) for i in range(N-lag))
    corr_den = sum((gaps[i] - mean_g)**2 for i in range(N))
    corr = corr_num / corr_den
    print(f"  lag {lag:2d}: autocorrelation = {corr:.6f}")

# Test 2: Gap prediction from history (linear regression on last K gaps)
print("\n--- Gap prediction accuracy ---")
for K in [1, 2, 3, 5, 10, 20]:
    correct = 0
    total = 0
    for i in range(K, min(N, 10000)):
        # Simple: predict g(i) from average of last K gaps
        predicted = sum(gaps[i-j-1] for j in range(K)) / K
        predicted_gap = 2 * round(predicted / 2)  # even
        if predicted_gap == gaps[i]:
            correct += 1
        total += 1
    print(f"  K={K:2d}: {correct}/{total} ({100*correct/total:.2f}%) exact gap predictions")

# Test 3: Gap prediction using ln(p) model
print("\n--- Gap prediction from ln(p) model ---")
correct = 0
total = 0
for i in range(1, min(N, 10000)):
    p = primes[i]
    expected_gap = math.log(p)
    predicted_gap = 2 * round(expected_gap / 2)
    if predicted_gap == gaps[i]:
        correct += 1
    total += 1
print(f"  ln(p) model: {correct}/{total} ({100*correct/total:.2f}%) exact")

# Test 4: Hardy-Littlewood prediction for specific gap sizes
print("\n--- Hardy-Littlewood singular series ---")
# P(gap = 2k) ≈ 2 * C_2 * Π_{p|k, p>2} (p-1)/(p-2) / ln²(N)
# where C_2 ≈ 0.6601618... is the twin prime constant

# Count gap frequencies
gap_freq = {}
for g in gaps[:50000]:
    gap_freq[g] = gap_freq.get(g, 0) + 1

print(f"  Top 15 most common gaps (first 50000 gaps):")
sorted_gaps = sorted(gap_freq.items(), key=lambda x: -x[1])
for g, count in sorted_gaps[:15]:
    frac = count / 50000
    expected = 2 / (math.log(primes[min(50000, len(primes)-1)]) ** 2) if g == 2 else 0
    print(f"    gap={g:3d}: count={count:5d} ({100*frac:.2f}%)")

# Test 5: CUMULATIVE gap prediction (does the average converge fast?)
print("\n--- Cumulative prediction: p(n) = 2 + Σ predicted_gaps ---")
# Strategy: use g(k) ≈ ln(p̂(k)) where p̂ is the running estimate
p_est = 2.0
max_n = 10000
errors = []
for k in range(1, max_n):
    # Predict gap from current estimate
    g_pred = math.log(max(p_est, 2))
    p_est += g_pred
    actual = primes[k]
    error = abs(p_est - actual)
    if k in [100, 500, 1000, 5000, 10000-1]:
        errors.append((k+1, round(p_est), actual, error))

print(f"  {'n':>6s} {'predicted':>12s} {'actual':>12s} {'error':>10s} {'rel_err':>10s}")
for n, pred, act, err in errors:
    print(f"  {n:6d} {pred:12d} {act:12d} {err:10.1f} {err/act*100:9.4f}%")

# Test 6: Cramér random model - how well does it match?
print("\n--- Cramér random model test ---")
print("If gaps are independent Poisson(ln p), then:")
print("  Var(p(n) - n ln n) should grow as n * (ln n)^2")
import random
random.seed(42)

# Simulate Cramér model
sim_errors = []
for trial in range(100):
    p_sim = 2
    for k in range(1, 1000):
        g = max(1, int(random.expovariate(1.0 / math.log(max(p_sim, 3)))))
        if g % 2 == 1 and p_sim > 2:
            g += 1
        p_sim += g

    n = 1000
    actual_1000 = primes[n-1]
    sim_errors.append(p_sim - actual_1000)

mean_err = sum(sim_errors) / len(sim_errors)
std_err = (sum((e - mean_err)**2 for e in sim_errors) / len(sim_errors))**0.5
print(f"  Cramér model for p(1000): mean_error={mean_err:.1f}, std={std_err:.1f}")
print(f"  Actual p(1000) = {primes[999]}")
print(f"  This std ≈ {std_err:.1f} means ~{std_err:.0f} candidates to check")

# Test 7: WHEEL FACTORIZATION approach to gap prediction
print("\n--- Wheel factorization for gap prediction ---")
# Primes > 5 are ≡ 1,7,11,13,17,19,23,29 mod 30
# So gaps between consecutive primes > 5 are from a restricted set
# Valid gaps mod 30: 2,4,6,8,10,12,14,16,18,20,22,24,26,28 (even only)
# But further constrained by mod 30 positions

# What fraction of gaps can be predicted just from the wheel?
wheel30 = [1,7,11,13,17,19,23,29]
wheel_gaps = []
for i in range(len(wheel30)):
    for j in range(len(wheel30)):
        for skip in range(5):
            g = wheel30[j] - wheel30[i] + 30 * skip
            if g > 0:
                wheel_gaps.append(g)

wheel_gap_set = set(wheel_gaps)
possible_gaps = set(g for g in gaps[3:] if g in wheel_gap_set)  # gaps after p=7
print(f"  Possible gaps from wheel mod 30: {len(wheel_gap_set)} values")
print(f"  Actual distinct gaps seen: {len(set(gaps[3:50000]))}")
print(f"  Gaps consistent with wheel: {sum(1 for g in gaps[3:50000] if g in wheel_gap_set)}/{min(50000-3, len(gaps)-3)}")

# FINAL ANALYSIS
print("\n" + "=" * 70)
print("CONCLUSION: Prime Gap Prediction")
print("=" * 70)
print("""
Gap prediction cannot work for exact p(n) because:

1. AUTOCORRELATION is near-zero (≈ -0.006 to 0.01)
   → Consecutive gaps are essentially INDEPENDENT
   → No sequential prediction possible

2. Best prediction (average of K previous gaps) gives ~20% accuracy
   → Must be EXACTLY right for each of n-1 gaps
   → Probability of getting ALL gaps right: 0.20^(n-1) → 0

3. Cumulative error: even small per-gap errors accumulate
   → After n gaps: error ≈ √n * σ_gap
   → For n = 10^100: error ≈ 10^50 * O(ln n) ≈ 10^52

4. The gap sequence has POSITIVE entropy:
   → Each gap carries ~log₂(ln p) ≈ log₂(ln(n ln n)) bits of information
   → For n = 10^100: each gap carries ~9 bits
   → Total information: n * 9 ≈ 9 × 10^100 bits
   → No polylog-parameter formula can encode this

The prime gap sequence is essentially INCOMPRESSIBLE.
This is another angle on the same fundamental barrier.
""")
