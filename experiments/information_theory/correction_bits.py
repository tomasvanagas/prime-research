#!/usr/bin/env python3
"""
Session 11: Bit-level analysis of the correction term delta(n) = p(n) - round(R^{-1}(n))

Goal: find ANY exploitable pattern in the bits of delta(n) that could allow
      O(polylog) computation instead of O(x^{2/3}) sieving.
"""

import sys, os, time, math, gzip, json
from collections import Counter, defaultdict
import numpy as np
from scipy import stats

import mpmath
from mpmath import mp, mpf, log, li, power, zeta, sqrt, pi, fsum
from sympy import primerange, prime, primepi, factorint, isprime

# ---------------------------------------------------------------------------
# 1. Riemann's R function and its inverse
# ---------------------------------------------------------------------------

def R(x, terms=200):
    """Riemann's R function: R(x) = sum_{k=1}^{terms} mu(k)/k * li(x^{1/k})"""
    # We use the Gram series which converges faster
    mp.dps = 50
    x = mpf(x)
    if x <= 1:
        return mpf(0)
    lnx = log(x)
    s = mpf(1)
    term = mpf(1)
    for k in range(1, terms + 1):
        term *= lnx / k
        s += term / (k * zeta(k + 1))
    return s


def R_inv(n, tol=1e-12, max_iter=100):
    """Inverse of R: given n, find x such that R(x) = n.
    Uses Newton's method with R'(x) = 1/(x*ln(x)) approximately."""
    mp.dps = 50
    n = mpf(n)
    if n <= 0:
        return mpf(2)

    # Initial guess using prime number theorem: p(n) ~ n*ln(n)
    nn = float(n)
    if nn < 2:
        x = mpf(2)
    else:
        guess = nn * math.log(nn) + nn * math.log(math.log(nn)) if nn > 5 else nn * 2.5
        x = mpf(guess)

    for _ in range(max_iter):
        rx = R(x)
        err = rx - n
        if abs(float(err)) < tol:
            break
        # R'(x) ~ 1/(ln(x)) for large x
        deriv = 1 / log(x)
        x = x - err / deriv
        if x < 2:
            x = mpf(2)
    return x


# ---------------------------------------------------------------------------
# 2. Compute delta(n) for a range of n
# ---------------------------------------------------------------------------

def compute_deltas(n_max, sample_points=None):
    """Compute delta(n) = p(n) - round(R^{-1}(n)) for n in range."""
    print(f"Computing deltas for n up to {n_max}...")

    # Get all primes up to the n_max-th prime
    # For n_max = 10^5, p(n_max) ~ 1.3 * 10^6
    upper_bound = int(n_max * (math.log(n_max) + math.log(math.log(n_max)) + 3)) if n_max > 10 else 100
    primes_list = list(primerange(2, upper_bound + 1))

    if len(primes_list) < n_max:
        upper_bound = int(upper_bound * 1.5)
        primes_list = list(primerange(2, upper_bound + 1))

    primes_list = primes_list[:n_max]
    print(f"  Got {len(primes_list)} primes, largest = {primes_list[-1]}")

    if sample_points is None:
        sample_points = list(range(1, min(n_max + 1, len(primes_list) + 1)))

    deltas = {}
    t0 = time.time()
    report_interval = max(1, len(sample_points) // 20)

    for idx, n in enumerate(sample_points):
        if n < 1 or n > len(primes_list):
            continue
        pn = primes_list[n - 1]
        r_inv_n = R_inv(n)
        rounded = int(round(float(r_inv_n)))
        delta = pn - rounded
        deltas[n] = delta

        if (idx + 1) % report_interval == 0:
            elapsed = time.time() - t0
            print(f"  {idx+1}/{len(sample_points)} done ({elapsed:.1f}s) last: n={n}, p(n)={pn}, R^-1={rounded}, delta={delta}")

    elapsed = time.time() - t0
    print(f"  Completed in {elapsed:.1f}s")
    return deltas, primes_list


# ---------------------------------------------------------------------------
# 3. Analysis functions
# ---------------------------------------------------------------------------

def analyze_basic_stats(deltas):
    """Basic statistics of delta values."""
    vals = list(deltas.values())
    ns = sorted(deltas.keys())

    abs_vals = [abs(v) for v in vals]
    signs = [1 if v >= 0 else -1 for v in vals]

    results = {
        'count': len(vals),
        'mean': float(np.mean(vals)),
        'std': float(np.std(vals)),
        'median': float(np.median(vals)),
        'abs_mean': float(np.mean(abs_vals)),
        'abs_max': max(abs_vals),
        'positive_frac': sum(1 for v in vals if v >= 0) / len(vals),
        'zero_count': sum(1 for v in vals if v == 0),
    }

    # Bit-length statistics
    bit_lengths = [v.bit_length() if v != 0 else 0 for v in abs_vals]
    results['mean_bit_length'] = float(np.mean(bit_lengths))
    results['max_bit_length'] = max(bit_lengths)

    # How does |delta| scale with n?
    if len(ns) > 100:
        log_n = [math.log(n) for n in ns if n > 1]
        log_abs_delta = [math.log(abs(deltas[n]) + 1) for n in ns if n > 1]
        if len(log_n) > 2:
            slope, intercept, r, p, se = stats.linregress(log_n, log_abs_delta)
            results['scaling_exponent'] = slope
            results['scaling_r_squared'] = r**2

    return results


def analyze_bit_patterns(deltas):
    """Analyze binary representation of |delta(n)|."""
    ns = sorted(deltas.keys())
    results = {}

    # Extract bit sequences
    # For each delta, get its binary representation (sign + magnitude)
    signs = []
    magnitudes = []
    for n in ns:
        d = deltas[n]
        signs.append(0 if d >= 0 else 1)
        magnitudes.append(abs(d))

    # --- Sign sequence analysis ---
    sign_runs = []
    current_run = 1
    for i in range(1, len(signs)):
        if signs[i] == signs[i-1]:
            current_run += 1
        else:
            sign_runs.append(current_run)
            current_run = 1
    sign_runs.append(current_run)

    results['sign_mean_run_length'] = float(np.mean(sign_runs))
    results['sign_max_run_length'] = max(sign_runs)
    # For random binary: mean run ~ 2

    # Sign autocorrelation
    s_arr = np.array(signs, dtype=float)
    s_arr = s_arr - s_arr.mean()
    if np.std(s_arr) > 0:
        autocorr = np.correlate(s_arr, s_arr, mode='full')
        autocorr = autocorr[len(autocorr)//2:]
        autocorr = autocorr / autocorr[0]
        results['sign_autocorr_lag1'] = float(autocorr[1]) if len(autocorr) > 1 else 0
        results['sign_autocorr_lag2'] = float(autocorr[2]) if len(autocorr) > 2 else 0
        results['sign_autocorr_lag5'] = float(autocorr[5]) if len(autocorr) > 5 else 0
        results['sign_autocorr_lag10'] = float(autocorr[10]) if len(autocorr) > 10 else 0

    # --- Individual bit analysis ---
    # Look at bit k of |delta(n)| across all n
    max_bits = max(m.bit_length() for m in magnitudes if m > 0)
    bit_freqs = {}
    for k in range(max_bits):
        bits_at_k = [(m >> k) & 1 for m in magnitudes]
        ones = sum(bits_at_k)
        total = len(bits_at_k)
        bit_freqs[k] = ones / total
    results['bit_frequencies'] = {k: round(v, 4) for k, v in bit_freqs.items()}

    # --- Consecutive delta correlation ---
    delta_vals = [deltas[n] for n in ns]
    if len(delta_vals) > 2:
        d_arr = np.array(delta_vals, dtype=float)
        # Pearson correlation between consecutive deltas
        r_consec, p_consec = stats.pearsonr(d_arr[:-1], d_arr[1:])
        results['consecutive_pearson_r'] = float(r_consec)
        results['consecutive_pearson_p'] = float(p_consec)

        # Spearman rank correlation
        rho, p_rho = stats.spearmanr(d_arr[:-1], d_arr[1:])
        results['consecutive_spearman_rho'] = float(rho)

    # --- XOR analysis: delta(n) XOR delta(n+1) ---
    xor_bits = []
    for i in range(len(ns) - 1):
        if ns[i+1] == ns[i] + 1:
            x = abs(deltas[ns[i]]) ^ abs(deltas[ns[i+1]])
            xor_bits.append(x.bit_length() if x > 0 else 0)
    if xor_bits:
        results['xor_consecutive_mean_bits'] = float(np.mean(xor_bits))
        results['xor_consecutive_max_bits'] = max(xor_bits)

    return results


def analyze_modular_patterns(deltas):
    """Check if delta(n) mod m has non-uniform distribution for small m."""
    ns = sorted(deltas.keys())
    results = {}

    for m in [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 24, 30]:
        residues = [deltas[n] % m for n in ns]
        counts = Counter(residues)
        # Expected uniform: each residue appears len(ns)/m times
        expected = len(ns) / m
        chi2 = sum((counts.get(r, 0) - expected)**2 / expected for r in range(m))
        # p-value from chi-squared test
        p_value = 1 - stats.chi2.cdf(chi2, m - 1)
        results[f'mod_{m}_chi2'] = round(chi2, 2)
        results[f'mod_{m}_p_value'] = round(p_value, 6)
        results[f'mod_{m}_distribution'] = dict(sorted(counts.items()))

    return results


def analyze_compression(deltas):
    """Estimate Kolmogorov complexity via compression ratio."""
    ns = sorted(deltas.keys())
    results = {}

    # Convert deltas to byte sequences
    delta_bytes = b''
    sign_bytes = b''
    magnitude_bytes = b''

    for n in ns:
        d = deltas[n]
        sign_bytes += (0 if d >= 0 else 1).to_bytes(1, 'big')
        mag = abs(d)
        nbytes = max(1, (mag.bit_length() + 7) // 8)
        magnitude_bytes += mag.to_bytes(nbytes, 'big')
        # Fixed-width encoding for delta
        d_shifted = d + 2**20  # shift to make positive
        delta_bytes += d_shifted.to_bytes(4, 'big')

    raw_size = len(delta_bytes)
    compressed = gzip.compress(delta_bytes, compresslevel=9)
    results['raw_bytes'] = raw_size
    results['compressed_bytes'] = len(compressed)
    results['compression_ratio'] = round(len(compressed) / raw_size, 4)

    # Sign sequence compression
    sign_raw = len(sign_bytes)
    sign_comp = len(gzip.compress(sign_bytes, compresslevel=9))
    results['sign_compression_ratio'] = round(sign_comp / sign_raw, 4)

    # Magnitude compression
    mag_raw = len(magnitude_bytes)
    mag_comp = len(gzip.compress(magnitude_bytes, compresslevel=9))
    results['magnitude_compression_ratio'] = round(mag_comp / mag_raw, 4)

    # Compare with random: generate random deltas with same distribution
    rng = np.random.default_rng(42)
    abs_vals = [abs(deltas[n]) for n in ns]
    random_signs = rng.choice([-1, 1], size=len(ns))
    random_mags = rng.choice(abs_vals, size=len(ns))
    random_deltas = random_signs * random_mags

    rand_bytes = b''
    for d in random_deltas:
        d_shifted = int(d) + 2**20
        rand_bytes += d_shifted.to_bytes(4, 'big')

    rand_comp = len(gzip.compress(rand_bytes, compresslevel=9))
    results['random_compression_ratio'] = round(rand_comp / len(rand_bytes), 4)
    results['compressibility_advantage'] = round(
        (rand_comp / len(rand_bytes)) / (len(compressed) / raw_size), 4
    )

    return results


def analyze_n_features(deltas, primes_list):
    """Try to predict delta(n) from features of n."""
    ns = sorted(deltas.keys())
    results = {}

    # Features of n
    def digit_sum(x):
        return sum(int(d) for d in str(x))

    def num_prime_factors(x):
        if x < 2:
            return 0
        return sum(factorint(x).values())

    def omega(x):
        """Number of distinct prime factors."""
        if x < 2:
            return 0
        return len(factorint(x))

    # Build feature matrix (for a subset to keep it fast)
    subset = [n for n in ns if n <= 5000 and n >= 2]

    features = []
    targets_sign = []
    targets_mag = []
    targets_delta = []

    print("  Building feature matrix...")
    for n in subset:
        d = deltas[n]
        f = [
            n,
            math.log(n),
            n % 2, n % 3, n % 4, n % 5, n % 6,
            digit_sum(n),
            digit_sum(n) % 9,
            num_prime_factors(n),
            omega(n),
            1 if isprime(n) else 0,
            n % 10,
            n % 30,
        ]
        features.append(f)
        targets_sign.append(0 if d >= 0 else 1)
        targets_mag.append(abs(d))
        targets_delta.append(d)

    X = np.array(features, dtype=float)
    y_sign = np.array(targets_sign)
    y_mag = np.array(targets_mag, dtype=float)
    y_delta = np.array(targets_delta, dtype=float)

    # Correlation of each feature with delta
    feature_names = ['n', 'log_n', 'n%2', 'n%3', 'n%4', 'n%5', 'n%6',
                     'digit_sum', 'digit_sum%9', 'num_pf', 'omega', 'is_prime',
                     'n%10', 'n%30']

    correlations = {}
    for i, name in enumerate(feature_names):
        r, p = stats.pearsonr(X[:, i], y_delta)
        correlations[name] = {'r': round(r, 4), 'p': round(p, 6)}
    results['feature_correlations'] = correlations

    # Sign prediction accuracy using majority vote per feature bin
    best_sign_acc = 0.5
    best_sign_feature = None
    for i, name in enumerate(feature_names):
        if name in ['n', 'log_n']:
            continue
        bins = defaultdict(list)
        for j, n in enumerate(subset):
            bins[int(X[j, i])].append(y_sign[j])
        correct = 0
        total = 0
        for b, signs_in_bin in bins.items():
            majority = 1 if sum(signs_in_bin) > len(signs_in_bin)/2 else 0
            correct += sum(1 for s in signs_in_bin if s == majority)
            total += len(signs_in_bin)
        acc = correct / total if total > 0 else 0.5
        if acc > best_sign_acc:
            best_sign_acc = acc
            best_sign_feature = name

    results['best_sign_prediction_accuracy'] = round(best_sign_acc, 4)
    results['best_sign_feature'] = best_sign_feature

    # --- delta(n) vs prime gaps ---
    # Is delta correlated with the local prime gap?
    gaps = []
    delta_at_gap = []
    for n in subset:
        if n >= 2 and n < len(primes_list):
            gap = primes_list[n] - primes_list[n-1]  # gap after p(n)
            gaps.append(gap)
            delta_at_gap.append(deltas[n])

    if len(gaps) > 10:
        r_gap, p_gap = stats.pearsonr(gaps, delta_at_gap)
        results['delta_vs_gap_correlation'] = round(r_gap, 4)
        results['delta_vs_gap_p_value'] = round(p_gap, 6)

    return results


def analyze_difference_sequence(deltas):
    """Analyze delta(n+1) - delta(n) and higher-order differences."""
    ns = sorted(deltas.keys())
    results = {}

    # First differences
    consecutive = [(ns[i], ns[i+1]) for i in range(len(ns)-1) if ns[i+1] == ns[i] + 1]
    first_diff = [deltas[b] - deltas[a] for a, b in consecutive]

    if len(first_diff) > 10:
        results['first_diff_mean'] = round(float(np.mean(first_diff)), 4)
        results['first_diff_std'] = round(float(np.std(first_diff)), 4)
        results['first_diff_median'] = float(np.median(first_diff))

        # Note: first_diff = p(n+1) - p(n) - (round(R^{-1}(n+1)) - round(R^{-1}(n)))
        # = gap(n) - dR
        # So first differences encode the "gap prediction error"

        # Autocorrelation of first differences
        fd = np.array(first_diff, dtype=float)
        fd = fd - fd.mean()
        if np.std(fd) > 0:
            ac = np.correlate(fd, fd, mode='full')
            ac = ac[len(ac)//2:]
            ac = ac / ac[0]
            results['first_diff_autocorr_lag1'] = round(float(ac[1]), 4) if len(ac) > 1 else 0
            results['first_diff_autocorr_lag2'] = round(float(ac[2]), 4) if len(ac) > 2 else 0

    # Second differences
    if len(first_diff) > 10:
        second_diff = [first_diff[i+1] - first_diff[i] for i in range(len(first_diff)-1)]
        results['second_diff_mean'] = round(float(np.mean(second_diff)), 4)
        results['second_diff_std'] = round(float(np.std(second_diff)), 4)

    return results


def analyze_bit_predictability(deltas):
    """For each bit position k, try to predict bit k of |delta(n+1)| from |delta(n)|."""
    ns = sorted(deltas.keys())
    results = {}

    consecutive = [(ns[i], ns[i+1]) for i in range(len(ns)-1) if ns[i+1] == ns[i] + 1]

    max_bits = 15  # analyze up to 15 bits

    for k in range(max_bits):
        bits_curr = [((abs(deltas[a]) >> k) & 1) for a, b in consecutive]
        bits_next = [((abs(deltas[b]) >> k) & 1) for a, b in consecutive]

        # Same-bit prediction: does bit k of delta(n) predict bit k of delta(n+1)?
        agree = sum(1 for c, n in zip(bits_curr, bits_next) if c == n)
        total = len(bits_curr)
        results[f'bit_{k}_same_prediction'] = round(agree / total, 4) if total > 0 else 0.5

        # Cross-bit: does bit 0 of delta(n) predict bit k of delta(n+1)?
        bits0_curr = [((abs(deltas[a]) >> 0) & 1) for a, b in consecutive]
        agree0 = sum(1 for c, n in zip(bits0_curr, bits_next) if c == n)
        results[f'bit_{k}_from_bit0'] = round(agree0 / total, 4) if total > 0 else 0.5

    return results


def analyze_local_structure(deltas, primes_list):
    """Check if delta(n) has local structure — e.g., in windows of 30 (primorial)."""
    ns = sorted(deltas.keys())
    results = {}

    # Delta distribution by n mod 30
    mod30_stats = defaultdict(list)
    for n in ns:
        mod30_stats[n % 30].append(deltas[n])

    mod30_means = {}
    for r in sorted(mod30_stats.keys()):
        vals = mod30_stats[r]
        if len(vals) > 5:
            mod30_means[r] = round(float(np.mean(vals)), 2)
    results['delta_mean_by_n_mod30'] = mod30_means

    # ANOVA: is the mean delta significantly different across n mod 30 groups?
    groups = [mod30_stats[r] for r in sorted(mod30_stats.keys()) if len(mod30_stats[r]) > 5]
    if len(groups) > 2:
        f_stat, p_val = stats.f_oneway(*groups)
        results['anova_n_mod30_F'] = round(f_stat, 4)
        results['anova_n_mod30_p'] = round(p_val, 6)

    # Delta distribution by n mod 6
    mod6_stats = defaultdict(list)
    for n in ns:
        mod6_stats[n % 6].append(deltas[n])
    mod6_means = {}
    for r in sorted(mod6_stats.keys()):
        vals = mod6_stats[r]
        if len(vals) > 5:
            mod6_means[r] = round(float(np.mean(vals)), 2)
    results['delta_mean_by_n_mod6'] = mod6_means

    return results


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("SESSION 11: Bit-level analysis of delta(n) = p(n) - round(R^{-1}(n))")
    print("=" * 70)

    # Phase 1: Compute deltas
    # Start with N=10000 for thorough analysis (scales to more if fast enough)
    N = 10000

    t_start = time.time()
    deltas, primes_list = compute_deltas(N)
    t_compute = time.time() - t_start

    all_results = {'compute_time_s': round(t_compute, 1), 'N': N}

    # Phase 2: Run analyses
    print("\n--- Basic Statistics ---")
    basic = analyze_basic_stats(deltas)
    all_results['basic_stats'] = basic
    for k, v in basic.items():
        print(f"  {k}: {v}")

    print("\n--- Bit Pattern Analysis ---")
    bits = analyze_bit_patterns(deltas)
    all_results['bit_patterns'] = bits
    for k, v in bits.items():
        if k == 'bit_frequencies':
            print(f"  {k}:")
            for bk, bv in sorted(v.items()):
                print(f"    bit {bk}: P(1) = {bv}")
        else:
            print(f"  {k}: {v}")

    print("\n--- Modular Pattern Analysis ---")
    mod = analyze_modular_patterns(deltas)
    all_results['modular_patterns'] = mod
    for k, v in mod.items():
        if 'distribution' not in k:
            print(f"  {k}: {v}")

    print("\n--- Compression Analysis ---")
    comp = analyze_compression(deltas)
    all_results['compression'] = comp
    for k, v in comp.items():
        print(f"  {k}: {v}")

    print("\n--- Feature Correlation Analysis ---")
    feat = analyze_n_features(deltas, primes_list)
    all_results['features'] = feat
    if 'feature_correlations' in feat:
        print("  Feature correlations with delta:")
        for name, info in feat['feature_correlations'].items():
            sig = " ***" if info['p'] < 0.001 else " **" if info['p'] < 0.01 else " *" if info['p'] < 0.05 else ""
            print(f"    {name:15s}: r={info['r']:+.4f}  p={info['p']:.6f}{sig}")
    for k, v in feat.items():
        if k != 'feature_correlations':
            print(f"  {k}: {v}")

    print("\n--- Difference Sequence Analysis ---")
    diff = analyze_difference_sequence(deltas)
    all_results['differences'] = diff
    for k, v in diff.items():
        print(f"  {k}: {v}")

    print("\n--- Bit Predictability Analysis ---")
    bitpred = analyze_bit_predictability(deltas)
    all_results['bit_predictability'] = bitpred
    for k, v in bitpred.items():
        deviation = abs(v - 0.5)
        flag = " <<<" if deviation > 0.02 else ""
        print(f"  {k}: {v}{flag}")

    print("\n--- Local Structure Analysis ---")
    local = analyze_local_structure(deltas, primes_list)
    all_results['local_structure'] = local
    for k, v in local.items():
        print(f"  {k}: {v}")

    # Phase 3: Summary
    print("\n" + "=" * 70)
    print("SUMMARY OF FINDINGS")
    print("=" * 70)

    findings = []

    # Check compression
    if comp.get('compressibility_advantage', 1) > 1.05:
        findings.append(f"Delta sequence is {comp['compressibility_advantage']:.2f}x more compressible than random")
    else:
        findings.append("Delta sequence compresses similarly to random — high Kolmogorov complexity")

    # Check sign autocorrelation
    if 'sign_autocorr_lag1' in bits:
        ac1 = bits['sign_autocorr_lag1']
        if abs(ac1) > 0.02:
            findings.append(f"Sign autocorrelation at lag 1: {ac1:.4f} (non-trivial)")
        else:
            findings.append(f"Sign autocorrelation at lag 1: {ac1:.4f} (essentially random)")

    # Check consecutive correlation
    if 'consecutive_pearson_r' in bits:
        r = bits['consecutive_pearson_r']
        if abs(r) > 0.05:
            findings.append(f"Consecutive delta correlation: r={r:.4f} (meaningful)")
        else:
            findings.append(f"Consecutive delta correlation: r={r:.4f} (weak/none)")

    # Check modular patterns
    for m in [2, 3, 6, 30]:
        p_key = f'mod_{m}_p_value'
        if p_key in mod and mod[p_key] < 0.01:
            findings.append(f"Significant non-uniformity in delta mod {m} (p={mod[p_key]})")

    # Check bit prediction
    for k in range(5):
        key = f'bit_{k}_same_prediction'
        if key in bitpred and abs(bitpred[key] - 0.5) > 0.02:
            findings.append(f"Bit {k} shows predictability: {bitpred[key]:.4f} (vs 0.5 for random)")

    # Check feature correlations
    if 'feature_correlations' in feat:
        for name, info in feat['feature_correlations'].items():
            if info['p'] < 0.001 and abs(info['r']) > 0.05:
                findings.append(f"Feature '{name}' correlates with delta: r={info['r']:.4f}")

    # Check ANOVA
    if 'anova_n_mod30_p' in local and local['anova_n_mod30_p'] < 0.01:
        findings.append(f"ANOVA shows delta depends on n mod 30 (p={local['anova_n_mod30_p']})")

    # Check scaling
    if 'scaling_exponent' in basic:
        findings.append(f"|delta| scales as n^{basic['scaling_exponent']:.3f} (R^2={basic['scaling_r_squared']:.3f})")

    for i, f in enumerate(findings):
        print(f"  {i+1}. {f}")

    all_results['findings'] = findings

    # Save results
    outpath = os.path.join(os.path.dirname(__file__), 'results_session11.json')
    with open(outpath, 'w') as fp:
        # Convert any non-serializable types
        def convert(obj):
            if isinstance(obj, (np.integer,)):
                return int(obj)
            if isinstance(obj, (np.floating,)):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return obj

        json.dump(all_results, fp, indent=2, default=convert)
    print(f"\nResults saved to {outpath}")

    return all_results


if __name__ == '__main__':
    main()
