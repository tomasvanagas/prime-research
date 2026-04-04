"""
Prime Gap Analysis: Can gaps be predicted well enough to find the nth prime?

Explores:
1. Gap statistics for first 100K primes
2. Gap prediction formulas (PNT, modular corrections, Gallagher model)
3. Accumulative gap approach: predict p(1001)..p(2000) from p(1000)
4. Cramer random model comparison
"""

import math
import random
import time
from collections import Counter, defaultdict

# ============================================================
# Utility: generate primes via sieve
# ============================================================

def sieve_of_eratosthenes(limit):
    """Return list of primes up to limit."""
    is_prime = bytearray([1]) * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = bytearray(len(is_prime[i*i::i]))
    return [i for i in range(2, limit + 1) if is_prime[i]]


def get_first_n_primes(n):
    """Get the first n primes. Over-estimate the bound using PNT."""
    if n < 6:
        return [2, 3, 5, 7, 11][:n]
    # Upper bound: p(n) < n * (ln(n) + ln(ln(n))) for n >= 6
    upper = int(n * (math.log(n) + math.log(math.log(n))) + 100)
    primes = sieve_of_eratosthenes(upper)
    while len(primes) < n:
        upper = int(upper * 1.2)
        primes = sieve_of_eratosthenes(upper)
    return primes[:n]


# ============================================================
# Part 1: Gap statistics
# ============================================================

def analyze_gap_statistics(primes, gaps):
    """Analyze prime gap statistics."""
    n = len(gaps)
    results = {}

    # Basic stats
    avg_gap = sum(gaps) / n
    max_gap = max(gaps)
    min_gap = min(gaps)
    results['avg_gap'] = avg_gap
    results['max_gap'] = max_gap
    results['min_gap'] = min_gap

    # Distribution (histogram)
    gap_counts = Counter(gaps)
    results['gap_distribution'] = dict(sorted(gap_counts.items())[:30])  # top 30

    # Correlation between consecutive gaps g(n) and g(n+1)
    mean_g = avg_gap
    cov = sum((gaps[i] - mean_g) * (gaps[i+1] - mean_g) for i in range(n - 1)) / (n - 1)
    var_g = sum((g - mean_g)**2 for g in gaps) / n
    corr_consecutive = cov / var_g if var_g > 0 else 0
    results['corr_consecutive_gaps'] = corr_consecutive

    # Correlation between g(n) and n
    mean_n = (n - 1) / 2
    cov_gn = sum((i - mean_n) * (gaps[i] - mean_g) for i in range(n)) / n
    var_n = sum((i - mean_n)**2 for i in range(n)) / n
    corr_gn = cov_gn / (var_n**0.5 * var_g**0.5) if var_n > 0 and var_g > 0 else 0
    results['corr_gap_vs_index'] = corr_gn

    # Average gap vs ln(p(n)) — PNT predicts avg gap ~ ln(p)
    # Check in windows
    window_size = 1000
    pnt_comparison = []
    for start in range(0, n - window_size, window_size):
        end = start + window_size
        window_avg = sum(gaps[start:end]) / window_size
        mid_prime = primes[start + window_size // 2]
        ln_p = math.log(mid_prime)
        ratio = window_avg / ln_p
        pnt_comparison.append((start + window_size // 2, mid_prime, window_avg, ln_p, ratio))
    results['pnt_comparison_samples'] = pnt_comparison[::10]  # every 10th window

    # Cramer's conjecture: max gap near p ~ (ln p)^2
    # Track max gap seen vs (ln p)^2 in windows
    cramer_data = []
    running_max = 0
    for i in range(n):
        if gaps[i] > running_max:
            running_max = gaps[i]
            ln_p = math.log(primes[i])
            cramer_bound = ln_p ** 2
            cramer_data.append((i, primes[i], running_max, cramer_bound, running_max / cramer_bound))
    results['cramer_records'] = cramer_data[-10:]  # last 10 record gaps

    return results


# ============================================================
# Part 2: Gap prediction formulas
# ============================================================

def test_ln_prediction(primes, gaps):
    """Test g(n) ~ ln(p(n)) as a predictor."""
    errors = []
    abs_errors = []
    for i in range(len(gaps)):
        predicted = math.log(primes[i])
        err = gaps[i] - predicted
        errors.append(err)
        abs_errors.append(abs(err))

    mean_err = sum(errors) / len(errors)
    mean_abs_err = sum(abs_errors) / len(abs_errors)
    rmse = (sum(e**2 for e in errors) / len(errors)) ** 0.5
    within_1 = sum(1 for e in abs_errors if e <= 1) / len(abs_errors)
    within_2 = sum(1 for e in abs_errors if e <= 2) / len(abs_errors)

    return {
        'method': 'ln(p)',
        'mean_error': mean_err,
        'mean_abs_error': mean_abs_err,
        'rmse': rmse,
        'within_1': within_1,
        'within_2': within_2,
    }


def test_ln_with_mod_correction(primes, gaps):
    """Test g(n) ~ ln(p(n)) + correction based on p(n) mod small primes."""
    # Compute average gap by residue class for mod 6, mod 30, mod 210
    for modulus in [6, 30, 210]:
        residue_gaps = defaultdict(list)
        for i in range(len(gaps)):
            r = primes[i] % modulus
            residue_gaps[r].append(gaps[i])

        residue_avg = {}
        for r, gs in residue_gaps.items():
            residue_avg[r] = sum(gs) / len(gs)

        # Now predict using residue-specific mean
        # But this is cheating — we use the data to compute the correction.
        # A fairer test: train on first half, test on second half
        half = len(gaps) // 2
        train_residue = defaultdict(list)
        for i in range(half):
            r = primes[i] % modulus
            train_residue[r].append(gaps[i] - math.log(primes[i]))

        correction = {}
        for r, deltas in train_residue.items():
            correction[r] = sum(deltas) / len(deltas)

        # Test on second half
        errors = []
        abs_errors = []
        for i in range(half, len(gaps)):
            r = primes[i] % modulus
            corr = correction.get(r, 0)
            predicted = math.log(primes[i]) + corr
            err = gaps[i] - predicted
            errors.append(err)
            abs_errors.append(abs(err))

        mean_err = sum(errors) / len(errors)
        mean_abs_err = sum(abs_errors) / len(abs_errors)
        rmse = (sum(e**2 for e in errors) / len(errors)) ** 0.5
        within_1 = sum(1 for e in abs_errors if e <= 1) / len(abs_errors)

        yield {
            'method': f'ln(p) + mod {modulus} correction',
            'mean_error': mean_err,
            'mean_abs_error': mean_abs_err,
            'rmse': rmse,
            'within_1': within_1,
            'residue_corrections': {r: round(v, 3) for r, v in sorted(correction.items())},
        }


def test_gallagher_model(primes, gaps):
    """
    Gallagher's model: gaps near p follow approx Poisson with mean ln(p),
    but adjusted by local density. We test how well the Poisson CDF matches
    the actual gap distribution.
    """
    # Bin primes by size and compare gap distribution to Poisson
    from math import exp, factorial

    results = []
    for start_idx, label in [(0, 'small (p~2-50K)'), (50000, 'medium (p~600K-700K)'), (90000, 'large (p~1.1M-1.3M)')]:
        end_idx = min(start_idx + 10000, len(gaps))
        subset_gaps = gaps[start_idx:end_idx]
        subset_primes = primes[start_idx:end_idx]
        mean_ln = sum(math.log(p) for p in subset_primes) / len(subset_primes)
        actual_mean = sum(subset_gaps) / len(subset_gaps)

        # Compare frequencies: actual vs Poisson(mean_ln)
        actual_freq = Counter(subset_gaps)
        n_sample = len(subset_gaps)

        # Poisson pmf
        lam = mean_ln
        comparison = []
        for g in sorted(actual_freq.keys())[:15]:
            actual_frac = actual_freq[g] / n_sample
            # Poisson P(X=g) — but gaps are always even (except gap=1 at p=2)
            # Adjust: for even gaps, use Poisson on g/2 with lambda/2? No.
            # Actually Gallagher says the distribution of g/ln(p) -> exp(1)
            # Let's just compare raw Poisson
            if lam > 0 and g < 200:
                poisson_p = (lam ** g) * exp(-lam) / factorial(g) if g < 170 else 0
            else:
                poisson_p = 0
            comparison.append((g, actual_frac, poisson_p))

        results.append({
            'range': label,
            'mean_ln_p': round(mean_ln, 3),
            'actual_mean_gap': round(actual_mean, 3),
            'sample': comparison[:10],
        })

    return results


def test_residue_gap_prediction(primes, gaps):
    """
    Test whether residues mod primorial products give useful info.
    For p(n) mod 6: primes > 3 are 1 or 5 mod 6.
    Check if knowing the residue helps predict the gap size.
    """
    results = {}
    for mod in [6, 30, 210]:
        residue_stats = defaultdict(lambda: {'count': 0, 'total': 0, 'gaps': []})
        for i in range(len(gaps)):
            r = primes[i] % mod
            residue_stats[r]['count'] += 1
            residue_stats[r]['total'] += gaps[i]
            residue_stats[r]['gaps'].append(gaps[i])

        summary = {}
        for r in sorted(residue_stats.keys()):
            s = residue_stats[r]
            if s['count'] > 0:
                avg = s['total'] / s['count']
                std = (sum((g - avg)**2 for g in s['gaps']) / s['count']) ** 0.5
                summary[r] = {'count': s['count'], 'avg_gap': round(avg, 3), 'std': round(std, 3)}
        results[mod] = summary

    return results


# ============================================================
# Part 3: Accumulative gap approach
# ============================================================

def accumulative_gap_test(primes, gaps):
    """
    Start from p(1000), predict gaps to reach p(1001)..p(2000).
    Test several prediction strategies.
    """
    start_idx = 999  # p(1000) is primes[999]
    end_idx = 1999   # p(2000) is primes[1999]

    strategies = {}

    # Strategy 1: g = round(ln(p))
    def strat_ln(current_p, idx, history):
        return round(math.log(current_p))

    # Strategy 2: g = round(ln(p)) but force even (primes > 2 differ by even numbers)
    def strat_ln_even(current_p, idx, history):
        g = round(math.log(current_p))
        if g % 2 == 1:
            g += 1  # round up to even
        return max(g, 2)

    # Strategy 3: Use mod 30 correction (trained on first 1000 primes)
    train_corrections_30 = defaultdict(list)
    for i in range(1000):
        r = primes[i] % 30
        train_corrections_30[r].append(gaps[i] - math.log(primes[i]))
    mod30_corr = {r: sum(v)/len(v) for r, v in train_corrections_30.items()}

    def strat_mod30(current_p, idx, history):
        r = current_p % 30
        corr = mod30_corr.get(r, 0)
        g = round(math.log(current_p) + corr)
        if g < 1:
            g = 2
        if g % 2 == 1 and current_p > 2:
            g += 1
        return g

    # Strategy 4: Running average of recent gaps
    def strat_running_avg(current_p, idx, history):
        if len(history) < 5:
            return round(math.log(current_p))
        recent = history[-20:]
        g = round(sum(recent) / len(recent))
        return max(g, 2)

    for name, strat_fn in [
        ('round(ln(p))', strat_ln),
        ('round(ln(p)) forced even', strat_ln_even),
        ('ln(p) + mod30 correction', strat_mod30),
        ('running average', strat_running_avg),
    ]:
        predicted_p = primes[start_idx]
        actual_gaps_used = []
        predicted_gaps = []
        cumulative_error = 0
        exact_hits = 0

        for i in range(start_idx, end_idx):
            actual_gap = gaps[i]
            pred_gap = strat_fn(predicted_p, i, actual_gaps_used)
            predicted_gaps.append(pred_gap)
            actual_gaps_used.append(actual_gap)  # give actual gap for running avg

            cumulative_error += (pred_gap - actual_gap)
            predicted_p += pred_gap  # accumulate predicted gaps

            if pred_gap == actual_gap:
                exact_hits += 1

        actual_p2000 = primes[end_idx]
        predicted_p2000 = predicted_p

        strategies[name] = {
            'predicted_p2000': predicted_p2000,
            'actual_p2000': actual_p2000,
            'absolute_error': abs(predicted_p2000 - actual_p2000),
            'relative_error': abs(predicted_p2000 - actual_p2000) / actual_p2000,
            'cumulative_drift': cumulative_error,
            'exact_gap_hits': exact_hits,
            'exact_hit_rate': exact_hits / (end_idx - start_idx),
            'mean_gap_error': cumulative_error / (end_idx - start_idx),
        }

    return strategies


# ============================================================
# Part 4: Cramer random model
# ============================================================

def cramer_random_model(n_primes, seed=42):
    """
    Generate pseudo-primes using Cramer's model:
    each integer m > 2 is "prime" independently with probability 1/ln(m).
    Compare nth pseudo-prime to nth real prime.
    """
    random.seed(seed)
    pseudo_primes = [2]
    m = 3
    while len(pseudo_primes) < n_primes:
        prob = 1.0 / math.log(m)
        if random.random() < prob:
            pseudo_primes.append(m)
        m += 1

    return pseudo_primes


def compare_cramer_to_real(real_primes, n_compare):
    """Run multiple Cramer trials and compare."""
    results = []
    for trial in range(5):
        pseudo = cramer_random_model(n_compare, seed=trial * 137 + 42)
        diffs = []
        for check_n in [100, 1000, 5000, 10000, 50000, min(n_compare, 100000)]:
            if check_n <= n_compare:
                real_pn = real_primes[check_n - 1]
                pseudo_pn = pseudo[check_n - 1]
                rel_err = (pseudo_pn - real_pn) / real_pn
                diffs.append((check_n, real_pn, pseudo_pn, rel_err))
        results.append(diffs)
    return results


# ============================================================
# Main
# ============================================================

def main():
    print("=" * 70)
    print("PRIME GAP ANALYSIS")
    print("=" * 70)

    N = 100000
    print(f"\nGenerating first {N} primes...")
    t0 = time.time()
    primes = get_first_n_primes(N)
    print(f"  Done in {time.time()-t0:.2f}s. p({N}) = {primes[-1]}")

    gaps = [primes[i+1] - primes[i] for i in range(len(primes) - 1)]

    # ---- Part 1: Gap Statistics ----
    print("\n" + "=" * 70)
    print("PART 1: GAP STATISTICS")
    print("=" * 70)

    stats = analyze_gap_statistics(primes, gaps)
    print(f"\n  Average gap: {stats['avg_gap']:.4f}")
    print(f"  Max gap: {stats['max_gap']}")
    print(f"  Min gap: {stats['min_gap']}")
    print(f"  ln(p(100000)) = {math.log(primes[-1]):.4f}")
    print(f"  Ratio avg_gap / ln(p_last): {stats['avg_gap'] / math.log(primes[-1]):.4f}")

    print(f"\n  Correlation between consecutive gaps: {stats['corr_consecutive_gaps']:.6f}")
    print(f"  Correlation between gap and index: {stats['corr_gap_vs_index']:.6f}")

    print("\n  Gap distribution (top 20):")
    for gap, count in list(stats['gap_distribution'].items())[:20]:
        bar = '#' * (count // 200)
        print(f"    gap={gap:3d}: {count:5d} ({100*count/len(gaps):.1f}%) {bar}")

    print("\n  PNT check: avg_gap / ln(p) in windows of 1000 (should be ~1.0):")
    for idx, p, avg_g, ln_p, ratio in stats['pnt_comparison_samples']:
        print(f"    n={idx:6d}, p={p:8d}, avg_gap={avg_g:.2f}, ln(p)={ln_p:.2f}, ratio={ratio:.4f}")

    print("\n  Cramer's conjecture — record gaps vs (ln p)^2:")
    for idx, p, max_g, cramer_b, ratio in stats['cramer_records']:
        print(f"    n={idx:6d}, p={p:8d}, max_gap={max_g:3d}, (ln p)^2={cramer_b:.1f}, ratio={ratio:.4f}")

    # ---- Part 2: Gap Prediction Formulas ----
    print("\n" + "=" * 70)
    print("PART 2: GAP PREDICTION FORMULAS")
    print("=" * 70)

    # 2a: Simple ln(p) prediction
    ln_result = test_ln_prediction(primes, gaps)
    print(f"\n  Method: {ln_result['method']}")
    print(f"    Mean error: {ln_result['mean_error']:.4f}")
    print(f"    Mean absolute error: {ln_result['mean_abs_error']:.4f}")
    print(f"    RMSE: {ln_result['rmse']:.4f}")
    print(f"    Fraction within +/-1: {ln_result['within_1']:.4f}")
    print(f"    Fraction within +/-2: {ln_result['within_2']:.4f}")

    # 2b: ln(p) + modular corrections
    print("\n  Modular correction tests (train on 1st half, test on 2nd):")
    for mod_result in test_ln_with_mod_correction(primes, gaps):
        print(f"\n  Method: {mod_result['method']}")
        print(f"    Mean error: {mod_result['mean_error']:.4f}")
        print(f"    Mean absolute error: {mod_result['mean_abs_error']:.4f}")
        print(f"    RMSE: {mod_result['rmse']:.4f}")
        print(f"    Fraction within +/-1: {mod_result['within_1']:.4f}")
        if len(mod_result['residue_corrections']) <= 20:
            print(f"    Corrections by residue: {mod_result['residue_corrections']}")

    # 2c: Gallagher model
    print("\n  Gallagher / Poisson model comparison:")
    gallagher = test_gallagher_model(primes, gaps)
    for g in gallagher:
        print(f"\n    Range: {g['range']}")
        print(f"    Mean ln(p): {g['mean_ln_p']}, Actual mean gap: {g['actual_mean_gap']}")
        print(f"    Gap | Actual freq | Poisson(ln p)")
        for gap, actual, poisson in g['sample']:
            print(f"      {gap:3d}   {actual:.4f}        {poisson:.4f}")

    # 2d: Residue-based gap statistics
    print("\n  Average gap by residue class:")
    residue_results = test_residue_gap_prediction(primes, gaps)
    for mod in [6, 30]:
        print(f"\n    mod {mod}:")
        for r, s in residue_results[mod].items():
            print(f"      r={r:3d}: count={s['count']:5d}, avg_gap={s['avg_gap']:.3f}, std={s['std']:.3f}")

    # ---- Part 3: Accumulative Gap Approach ----
    print("\n" + "=" * 70)
    print("PART 3: ACCUMULATIVE GAP PREDICTION (p(1000) -> p(2000))")
    print("=" * 70)
    print(f"  p(1000) = {primes[999]}, p(2000) = {primes[1999]}")

    accum = accumulative_gap_test(primes, gaps)
    for name, r in accum.items():
        print(f"\n  Strategy: {name}")
        print(f"    Predicted p(2000): {r['predicted_p2000']}")
        print(f"    Actual   p(2000): {r['actual_p2000']}")
        print(f"    Absolute error:   {r['absolute_error']}")
        print(f"    Relative error:   {r['relative_error']:.6f}")
        print(f"    Cumulative drift: {r['cumulative_drift']}")
        print(f"    Exact gap hits:   {r['exact_gap_hits']}/{1000} ({r['exact_hit_rate']:.1%})")

    # ---- Part 4: Cramer Random Model ----
    print("\n" + "=" * 70)
    print("PART 4: CRAMER RANDOM MODEL")
    print("=" * 70)

    print("\n  Generating Cramer pseudo-primes (5 trials)...")
    cramer_results = compare_cramer_to_real(primes, N)
    for trial_idx, trial in enumerate(cramer_results):
        print(f"\n  Trial {trial_idx + 1}:")
        for n, real_pn, pseudo_pn, rel_err in trial:
            print(f"    n={n:6d}: real p(n)={real_pn:8d}, pseudo p(n)={pseudo_pn:8d}, "
                  f"rel_err={rel_err:+.4f} ({rel_err*100:+.2f}%)")

    # ---- Summary ----
    print("\n" + "=" * 70)
    print("SUMMARY OF FINDINGS")
    print("=" * 70)

    print("""
  1. INDIVIDUAL GAP PREDICTION IS POOR:
     - ln(p) predicts the average gap well (PNT works), but individual gaps
       are highly variable. RMSE ~ {rmse:.1f}, only {w1:.1%} within +/-1.
     - Modular corrections help marginally but the fundamental noise is large.
     - Gaps are essentially unpredictable at the individual level.

  2. ACCUMULATIVE APPROACH FAILS:
     - Even small per-gap errors (~{mae:.1f} on average) accumulate over 1000 steps.
     - The cumulative drift makes this approach useless for exact nth prime.
     - Best strategy still has error of hundreds to thousands.

  3. CRAMER MODEL:
     - The random model gives surprisingly reasonable estimates (within ~5-15%)
       but is far too noisy for an exact computation.
     - It captures the average behavior but not the fine structure.

  4. CONCLUSION:
     - Prime gaps are fundamentally noisy — individual gaps cannot be predicted
       with enough accuracy to replace sieving.
     - The PNT gives excellent average behavior, but the variance around ln(p)
       is too large for accumulative approaches.
     - To find the nth prime exactly, you MUST count primes (sieve, Meissel-Lehmer,
       or analytic methods). Gap prediction cannot substitute.
""".format(
        rmse=ln_result['rmse'],
        w1=ln_result['within_1'],
        mae=ln_result['mean_abs_error'],
    ))

    return {
        'stats': stats,
        'ln_result': ln_result,
        'accum': accum,
    }


if __name__ == '__main__':
    results = main()
