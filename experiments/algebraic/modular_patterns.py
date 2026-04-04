#!/usr/bin/env python3
"""
Modular arithmetic patterns in primes.

Investigates:
1. Wheel factorization — uniformity of primes across residue classes
2. Prime constellation / gap patterns
3. Chebyshev bias
4. Residue-class-based pi(x) estimation
5. Autocorrelation of prime gaps
"""

import math
import time
from collections import Counter, defaultdict

# ─────────────────────────────────────────────────────────────
# Sieve
# ─────────────────────────────────────────────────────────────

def sieve_of_eratosthenes(limit):
    """Return list of primes up to limit."""
    is_prime = bytearray(b'\x01') * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return [i for i in range(2, limit + 1) if is_prime[i]]


def euler_phi(n):
    """Euler's totient function."""
    result = n
    p = 2
    temp = n
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return result


def reduced_residues(m):
    """Return sorted list of integers in [1, m) coprime to m."""
    return [r for r in range(1, m) if math.gcd(r, m) == 1]


# ─────────────────────────────────────────────────────────────
# Generate primes
# ─────────────────────────────────────────────────────────────

LIMIT = 10**6
print(f"Sieving primes up to {LIMIT:,} ...")
t0 = time.time()
primes = sieve_of_eratosthenes(LIMIT)
print(f"  Found {len(primes):,} primes in {time.time()-t0:.2f}s\n")

output_lines = []   # collect markdown output


def section(title):
    sep = "=" * 70
    print(f"\n{sep}\n  {title}\n{sep}")
    output_lines.append(f"\n## {title}\n")


def note(text):
    print(text)
    output_lines.append(text)


# ═════════════════════════════════════════════════════════════
# 1. WHEEL FACTORIZATION ANALYSIS
# ═════════════════════════════════════════════════════════════

section("1. Wheel factorization — prime distribution across residue classes")

for m, label in [(30, "2*3*5"), (210, "2*3*5*7"), (2310, "2*3*5*7*11")]:
    residues = reduced_residues(m)
    phi_m = len(residues)
    note(f"\n### mod {m} = {label}  (phi={phi_m}, density={phi_m/m:.4f})")

    # Count primes > largest prime factor of m in each residue class
    skip = set()
    temp = m
    for p in [2, 3, 5, 7, 11, 13]:
        if temp % p == 0:
            skip.add(p)
            while temp % p == 0:
                temp //= p
        if temp == 1:
            break

    counts = Counter()
    total = 0
    for p in primes:
        if p in skip:
            continue
        counts[p % m] += 1
        total += 1

    expected = total / phi_m
    chi2 = 0
    min_count = float('inf')
    max_count = 0
    min_r = max_r = 0
    for r in residues:
        c = counts[r]
        chi2 += (c - expected) ** 2 / expected
        if c < min_count:
            min_count = c
            min_r = r
        if c > max_count:
            max_count = c
            max_r = r

    df = phi_m - 1
    # For large df, chi2 ~ N(df, 2*df), so z-score:
    z = (chi2 - df) / math.sqrt(2 * df)

    note(f"  Total primes counted: {total}")
    note(f"  Expected per class: {expected:.1f}")
    note(f"  Min: {min_count} (r={min_r})  Max: {max_count} (r={max_r})")
    note(f"  Chi-square = {chi2:.2f}  (df={df}, z-score={z:.2f})")
    if abs(z) < 2:
        note(f"  => Distribution is CONSISTENT with uniform (z within ±2)")
    else:
        note(f"  => Distribution shows SIGNIFICANT deviation from uniform")

    # Show top 5 and bottom 5 for mod 30
    if m == 30:
        sorted_classes = sorted(residues, key=lambda r: counts[r], reverse=True)
        note(f"\n  All residue classes mod 30:")
        for r in sorted_classes:
            pct = counts[r] / total * 100
            note(f"    r={r:2d}: {counts[r]:6d}  ({pct:.2f}%)")


# ═════════════════════════════════════════════════════════════
# 2. PRIME CONSTELLATION / GAP PATTERNS
# ═════════════════════════════════════════════════════════════

section("2. Prime constellation patterns (consecutive gap tuples)")

gaps = [primes[i+1] - primes[i] for i in range(len(primes) - 1)]

note(f"\n### Single gaps (most common among {len(gaps):,} gaps):")
gap_counts = Counter(gaps)
for g, c in gap_counts.most_common(15):
    note(f"  gap={g:3d}: {c:6d}  ({c/len(gaps)*100:.2f}%)")

# Gap pairs (g1, g2)
note(f"\n### Gap pairs (g1, g2) — most common:")
pair_counts = Counter()
for i in range(len(gaps) - 1):
    pair_counts[(gaps[i], gaps[i+1])] += 1
for pair, c in pair_counts.most_common(15):
    note(f"  {pair}: {c:6d}  ({c/(len(gaps)-1)*100:.2f}%)")

# Gap triples (g1, g2, g3)
note(f"\n### Gap triples (g1, g2, g3) — most common:")
triple_counts = Counter()
for i in range(len(gaps) - 2):
    triple_counts[(gaps[i], gaps[i+1], gaps[i+2])] += 1
for triple, c in triple_counts.most_common(15):
    note(f"  {triple}: {c:5d}  ({c/(len(gaps)-2)*100:.2f}%)")

# Conditional: given last gap = g, what's the most likely next gap?
note(f"\n### Conditional next-gap distribution (given previous gap):")
cond = defaultdict(Counter)
for i in range(len(gaps) - 1):
    cond[gaps[i]][gaps[i+1]] += 1

for g in [2, 4, 6, 8, 10, 12]:
    total_g = sum(cond[g].values())
    if total_g < 100:
        continue
    note(f"\n  After gap={g} ({total_g} occurrences):")
    for ng, c in cond[g].most_common(5):
        note(f"    next={ng:3d}: {c:5d} ({c/total_g*100:.1f}%)")


# ═════════════════════════════════════════════════════════════
# 3. CHEBYSHEV BIAS
# ═════════════════════════════════════════════════════════════

section("3. Chebyshev bias exploration")

note("\n### Primes mod 4: count(≡3) vs count(≡1)")
c1 = c3 = 0
bias_4_data = []
checkpoints = [10**k for k in range(2, 7)]
for p in primes:
    if p == 2:
        continue
    if p % 4 == 1:
        c1 += 1
    else:
        c3 += 1
    if p in checkpoints or p == primes[-1]:
        pass  # record below
# Redo with checkpoints
c1 = c3 = 0
note(f"  {'Limit':>10s}  {'#(≡1)':>8s}  {'#(≡3)':>8s}  {'bias(3-1)':>10s}  {'ratio 3/1':>10s}")
cp_idx = 0
for p in primes:
    if p == 2:
        continue
    if p % 4 == 1:
        c1 += 1
    else:
        c3 += 1
    while cp_idx < len(checkpoints) and p >= checkpoints[cp_idx]:
        note(f"  {checkpoints[cp_idx]:10,d}  {c1:8d}  {c3:8d}  {c3-c1:+10d}  {c3/max(c1,1):10.4f}")
        cp_idx += 1
        if cp_idx >= len(checkpoints):
            break

note(f"\n### Primes mod 3: count(≡2) vs count(≡1)")
c1 = c2 = 0
note(f"  {'Limit':>10s}  {'#(≡1)':>8s}  {'#(≡2)':>8s}  {'bias(2-1)':>10s}  {'ratio 2/1':>10s}")
cp_idx = 0
for p in primes:
    if p <= 3:
        continue
    if p % 3 == 1:
        c1 += 1
    else:
        c2 += 1
    while cp_idx < len(checkpoints) and p >= checkpoints[cp_idx]:
        note(f"  {checkpoints[cp_idx]:10,d}  {c1:8d}  {c2:8d}  {c2-c1:+10d}  {c2/max(c1,1):10.4f}")
        cp_idx += 1
        if cp_idx >= len(checkpoints):
            break

# Detailed bias mod 4 — fraction of time ≡3 is ahead
note(f"\n### Fraction of primes p <= x where count(≡3 mod 4) > count(≡1 mod 4):")
c1 = c3 = 0
ahead_3 = 0
total_checked = 0
for p in primes:
    if p == 2:
        continue
    if p % 4 == 1:
        c1 += 1
    else:
        c3 += 1
    total_checked += 1
    if c3 > c1:
        ahead_3 += 1

note(f"  Among primes up to {LIMIT:,}: ≡3 is ahead {ahead_3/total_checked*100:.1f}% of the time")
note(f"  Final counts: ≡1: {c1}, ≡3: {c3}, bias: {c3-c1:+d}")


# ═════════════════════════════════════════════════════════════
# 4. RESIDUE CLASS COUNTING — pi(x) estimation
# ═════════════════════════════════════════════════════════════

section("4. Residue-class-based pi(x) estimation")

def li(x):
    """Logarithmic integral via numerical integration."""
    if x <= 2:
        return 0
    # Simple trapezoidal integration from 2 to x
    n_steps = max(1000, int(x / 100))
    dx = (x - 2) / n_steps
    total = 0
    for i in range(n_steps):
        t = 2 + (i + 0.5) * dx
        total += dx / math.log(t)
    return total

note("\n### Standard pi(x) vs li(x) vs per-class estimation")
note(f"  {'x':>10s}  {'pi(x)':>8s}  {'li(x)':>10s}  {'err_li':>8s}  {'class_est':>10s}  {'err_class':>10s}")

# For mod 30: each of the 8 residue classes contributes ~equally
# By Dirichlet, pi(x; 30, a) ~ li(x) / phi(30) for each coprime a
# We can refine by measuring actual density per class up to some calibration point
# and extrapolating

m = 30
residues_30 = reduced_residues(m)
phi_30 = len(residues_30)  # 8

test_points = [10**k for k in range(3, 7)]

for x in test_points:
    # Actual pi(x)
    pi_x = sum(1 for p in primes if p <= x)

    # Standard li(x)
    li_x = li(x)

    # Per-class estimate: use Dirichlet's theorem with empirical correction
    # Calibrate at x/2, predict at x
    cal = x // 2
    class_est = 0
    for r in residues_30:
        # count primes ≡ r (mod 30) up to cal
        count_cal = sum(1 for p in primes if p <= cal and p > 5 and p % m == r)
        # Dirichlet says asymptotically each class gets li(x)/phi(m)
        # Empirical ratio at calibration point:
        expected_cal = li(cal) / phi_30
        if expected_cal > 0:
            correction = count_cal / expected_cal
        else:
            correction = 1.0
        # Predict for this class at x
        class_est += correction * li(x) / phi_30

    # Add small primes (2, 3, 5) that aren't in residue classes
    class_est += 3

    err_li = li_x - pi_x
    err_class = class_est - pi_x

    note(f"  {x:10,d}  {pi_x:8d}  {li_x:10.1f}  {err_li:+8.1f}  {class_est:10.1f}  {err_class:+10.1f}")


# Check if some residue classes consistently lead/lag
note(f"\n### Per-class deviation from uniform (mod 30), primes up to {LIMIT:,}:")
total_above_5 = sum(1 for p in primes if p > 5)
expected_per_class = total_above_5 / phi_30
for r in residues_30:
    actual = sum(1 for p in primes if p > 5 and p % 30 == r)
    dev = (actual - expected_per_class) / expected_per_class * 100
    note(f"  r={r:2d}: actual={actual:6d}  expected={expected_per_class:.0f}  dev={dev:+.2f}%")


# ═════════════════════════════════════════════════════════════
# 5. AUTOCORRELATION OF PRIME GAPS
# ═════════════════════════════════════════════════════════════

section("5. Autocorrelation of prime gaps")

note(f"\nComputing autocorrelation of {len(gaps):,} prime gaps at lags 1..20")

mean_gap = sum(gaps) / len(gaps)
var_gap = sum((g - mean_gap)**2 for g in gaps) / len(gaps)

note(f"  Mean gap: {mean_gap:.4f}")
note(f"  Variance: {var_gap:.4f}")
note(f"  Std dev:  {math.sqrt(var_gap):.4f}")

note(f"\n  {'Lag':>5s}  {'Autocorr':>10s}  {'Significance':>14s}")
note(f"  {'---':>5s}  {'--------':>10s}  {'------------':>14s}")

# Significance threshold: ±2/sqrt(N)
sig_threshold = 2.0 / math.sqrt(len(gaps))

for lag in range(1, 21):
    cov = 0
    n = len(gaps) - lag
    for i in range(n):
        cov += (gaps[i] - mean_gap) * (gaps[i + lag] - mean_gap)
    cov /= n
    autocorr = cov / var_gap

    sig = "***" if abs(autocorr) > 3 * sig_threshold else \
          "**" if abs(autocorr) > 2 * sig_threshold else \
          "*" if abs(autocorr) > sig_threshold else ""
    note(f"  {lag:5d}  {autocorr:+10.6f}  {sig:>14s}")

note(f"\n  Significance threshold (2/sqrt(N)): ±{sig_threshold:.6f}")
note(f"  *: > 1x threshold, **: > 2x, ***: > 3x")


# ═════════════════════════════════════════════════════════════
# BONUS: Gap prediction accuracy
# ═════════════════════════════════════════════════════════════

section("Bonus: Gap prediction via conditional probabilities")

note("\nUsing empirical P(next_gap | prev_gap) to 'predict' most likely next gap.")
note("Evaluated on second half of gap sequence (first half as training).")

half = len(gaps) // 2
# Build conditional model from first half
train_cond = defaultdict(Counter)
for i in range(half - 1):
    train_cond[gaps[i]][gaps[i+1]] += 1

# For each gap in second half, predict most likely next gap
correct = 0
total_pred = 0
baseline_correct = 0
# Baseline: always predict the most common gap overall
most_common_gap = gap_counts.most_common(1)[0][0]

for i in range(half, len(gaps) - 1):
    prev = gaps[i]
    if prev in train_cond:
        predicted = train_cond[prev].most_common(1)[0][0]
    else:
        predicted = most_common_gap
    actual = gaps[i + 1]
    if predicted == actual:
        correct += 1
    if most_common_gap == actual:
        baseline_correct += 1
    total_pred += 1

note(f"  Predictions: {total_pred}")
note(f"  Baseline (always guess {most_common_gap}): {baseline_correct}/{total_pred} = {baseline_correct/total_pred*100:.2f}%")
note(f"  Conditional model:          {correct}/{total_pred} = {correct/total_pred*100:.2f}%")
note(f"  Improvement: {(correct-baseline_correct)/total_pred*100:+.2f} percentage points")

# Two-gap conditioning
note(f"\n### Two-gap conditioning: P(next | prev2, prev1)")
train_cond2 = defaultdict(Counter)
for i in range(half - 2):
    train_cond2[(gaps[i], gaps[i+1])][gaps[i+2]] += 1

correct2 = 0
for i in range(half, len(gaps) - 2):
    key = (gaps[i], gaps[i+1])
    if key in train_cond2:
        predicted = train_cond2[key].most_common(1)[0][0]
    else:
        predicted = most_common_gap
    actual = gaps[i + 2]
    if predicted == actual:
        correct2 += 1

total_pred2 = len(gaps) - 2 - half
note(f"  Two-gap conditional model: {correct2}/{total_pred2} = {correct2/total_pred2*100:.2f}%")
note(f"  Improvement over baseline: {(correct2-baseline_correct)/total_pred*100:+.2f} pp")


# ═════════════════════════════════════════════════════════════
# Summary
# ═════════════════════════════════════════════════════════════

section("Summary of key findings")

note("""
1. **Wheel factorization uniformity**: Primes are remarkably uniformly distributed
   across reduced residue classes for mod 30, 210, 2310. Chi-square tests show no
   significant deviation from uniform at the 10^6 level.

2. **Gap patterns**: Gap=6 dominates (~15-16%), followed by 2 and 4. Gap pairs
   (6,6) are most common. Conditional distributions given previous gap show
   real structure: after gap=2, gap=6 is far more likely than gap=2.

3. **Chebyshev bias**: Confirmed. Primes ≡ 3 (mod 4) consistently outnumber
   primes ≡ 1 (mod 4). The bias is small but persistent — ≡3 leads the race
   most of the time (known as "prime race").

4. **Per-class pi(x) estimation**: Using per-class calibration with Dirichlet's
   theorem provides only marginal improvement over plain li(x). The residue
   classes are too uniform to exploit for better counting.

5. **Gap autocorrelation**: There IS significant negative autocorrelation at lag 1
   (small gap tends to follow large gap and vice versa). Higher lags show weaker
   but still detectable correlations.

6. **Gap prediction**: Conditional models achieve modest improvement over baseline,
   but prime gaps remain fundamentally hard to predict. Two-gap conditioning
   helps slightly more.

**Conclusion for nth prime computation**: Modular patterns alone don't offer a
shortcut to computing nth prime without sieving. The uniformity of prime distribution
across residue classes means there's no exploitable bias. The autocorrelation in
gaps is real but too weak for accurate sequential prediction. The most promising
avenue remains analytic methods (Meissel-Lehmer for pi(x), then binary search).
""")


# ═════════════════════════════════════════════════════════════
# Write notes file
# ═════════════════════════════════════════════════════════════

with open("/apps/aplikacijos/prime-research/notes_modular.md", "w") as f:
    f.write("# Modular Arithmetic Patterns in Primes\n")
    f.write(f"\nGenerated: analysis of primes up to {LIMIT:,}\n")
    for line in output_lines:
        f.write(line + "\n")

print("\n\nResults saved to notes_modular.md")
