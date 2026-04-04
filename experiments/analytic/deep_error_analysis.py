"""
DEEP ERROR ANALYSIS: delta(n) = p(n) - R^{-1}(n)
===================================================
Comprehensive structural analysis of the error between
the nth prime and the inverse Riemann R function.

Tests:
  1. Compute delta(n) for n=2..50000
  2. Sign changes and their relation to prime gaps
  3. Periodicity related to primorial numbers
  4. delta(n) mod small numbers structure
  5. Distribution of normalized delta — Gaussian test, Rubinstein-Sarnak
  6. Second difference analysis
  7. Pattern of near-zero crossings
"""

import math
import time
import sys
import os
from collections import defaultdict, Counter

# ============================================================
# SIEVE
# ============================================================

def sieve(n):
    """Sieve of Eratosthenes up to n."""
    s = bytearray(b'\x01') * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n**0.5) + 1):
        if s[i]:
            s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i, v in enumerate(s) if v]

# ============================================================
# HIGH-PRECISION R^{-1} using mpmath
# ============================================================

def compute_R_inv_batch(N_max):
    """Compute R^{-1}(n) for n=1..N_max using fast float R function."""

    # Precompute Mobius values up to 50 (x^(1/k) negligible for k>50 when x<10^7)
    MU = [0] * 51
    MU[1] = 1
    for i in range(1, 51):
        for j in range(2*i, 51, i):
            MU[j] -= MU[i]

    # Active Mobius indices
    mu_indices = [(k, MU[k]) for k in range(1, 51) if MU[k] != 0]

    EULER_GAMMA = 0.5772156649015329

    def _li_fast(x):
        """Logarithmic integral li(x) via Ramanujan series (float)."""
        if x <= 1.0:
            return 0.0
        ln_x = math.log(x)
        r = EULER_GAMMA + math.log(abs(ln_x))
        t = 1.0
        for k in range(1, 100):
            t *= ln_x / k
            contrib = t / k
            r += contrib
            if abs(contrib) < 1e-15 * abs(r):
                break
        return r

    def R_func_fast(x):
        """Riemann R function using float arithmetic."""
        if x <= 1.0:
            return 0.0
        result = 0.0
        for k, mu_k in mu_indices:
            xk = x ** (1.0 / k)
            if xk <= 1.0001:
                break
            result += mu_k / k * _li_fast(xk)
        return result

    def inv_R_fast(n_val):
        """Newton inversion of R(x) = n using float."""
        if n_val <= 1:
            return 2.0
        ln_n = math.log(n_val)
        x = n_val * ln_n
        if n_val > 5:
            x += n_val * math.log(ln_n)

        for _ in range(50):
            rx = R_func_fast(x)
            rpx = 1.0 / math.log(x)
            dx = (n_val - rx) * math.log(x)  # (n - R(x)) / R'(x)
            x += dx
            if abs(dx) < 1e-10:
                break
        return x

    print(f"Computing R^{{-1}}(n) for n=1..{N_max} using fast float Riemann R...")
    t0 = time.time()
    r_inv = [0.0]  # index 0 unused
    for n in range(1, N_max + 1):
        r_inv.append(inv_R_fast(n))
        if n % 10000 == 0:
            elapsed = time.time() - t0
            rate = n / elapsed
            eta = (N_max - n) / rate
            print(f"  n={n:6d} | R^{{-1}}={r_inv[-1]:.4f} | "
                  f"{elapsed:.1f}s elapsed, ~{eta:.0f}s remaining")

    print(f"  Done in {time.time()-t0:.1f}s")

    # Verify accuracy at a few points
    test_primes = {1000: 7919, 5000: 48611, 10000: 104729}
    print("  Accuracy check:")
    for n, pn in test_primes.items():
        if n <= N_max:
            err = pn - r_inv[n]
            print(f"    R^{{-1}}({n}) = {r_inv[n]:.4f}, p({n}) = {pn}, delta = {err:.4f}")

    return r_inv

# ============================================================
# MAIN ANALYSIS
# ============================================================

def main():
    N_MAX = 50000
    SIEVE_LIMIT = 700000  # p(50000) ~ 611953

    print("=" * 70)
    print("DEEP ERROR ANALYSIS: delta(n) = p(n) - R^{-1}(n)")
    print("=" * 70)

    # Step 1: Get primes
    print(f"\n[1] Sieving primes up to {SIEVE_LIMIT}...")
    t0 = time.time()
    primes = sieve(SIEVE_LIMIT)
    print(f"    Found {len(primes)} primes in {time.time()-t0:.2f}s")
    if len(primes) < N_MAX:
        print(f"    ERROR: need {N_MAX} primes, only got {len(primes)}")
        SIEVE_LIMIT = int(SIEVE_LIMIT * 1.5)
        primes = sieve(SIEVE_LIMIT)
        print(f"    Re-sieved to {SIEVE_LIMIT}, got {len(primes)} primes")
    # primes[0]=2, primes[1]=3, ... so primes[n-1] = p(n)

    # Step 2: Compute R^{-1}
    r_inv = compute_R_inv_batch(N_MAX)

    # Step 3: Compute delta(n) = p(n) - R^{-1}(n)
    print(f"\n[2] Computing delta(n) for n=2..{N_MAX}...")
    delta = [0.0] * (N_MAX + 1)  # delta[n] = p(n) - R^{-1}(n)
    p = [0] * (N_MAX + 1)  # p[n] = nth prime
    for n in range(1, N_MAX + 1):
        p[n] = primes[n - 1]
        delta[n] = float(p[n]) - r_inv[n]

    # Prime gaps
    gaps = [0] * (N_MAX + 1)
    for n in range(1, N_MAX):
        gaps[n] = p[n + 1] - p[n]

    # ============================================================
    # TEST 2: SIGN CHANGES OF delta(n)
    # ============================================================
    print("\n" + "=" * 70)
    print("[TEST 2] SIGN CHANGES OF delta(n)")
    print("=" * 70)

    sign_changes = []
    for n in range(3, N_MAX + 1):
        if delta[n] * delta[n - 1] < 0:
            sign_changes.append(n)

    print(f"  Total sign changes in [3, {N_MAX}]: {len(sign_changes)}")
    print(f"  Average interval between sign changes: "
          f"{(N_MAX - 2) / max(len(sign_changes), 1):.2f}")

    # Are sign changes related to prime gaps?
    sc_gaps = [gaps[n - 1] for n in sign_changes if n - 1 >= 1 and n - 1 < N_MAX]
    all_gaps_in_range = [gaps[n] for n in range(2, N_MAX)]
    mean_gap_at_sc = sum(sc_gaps) / max(len(sc_gaps), 1)
    mean_gap_overall = sum(all_gaps_in_range) / max(len(all_gaps_in_range), 1)
    print(f"\n  Mean gap at sign changes: {mean_gap_at_sc:.4f}")
    print(f"  Mean gap overall:         {mean_gap_overall:.4f}")
    print(f"  Ratio: {mean_gap_at_sc / max(mean_gap_overall, 0.001):.4f}")

    # Gap distribution at sign changes vs overall
    sc_gap_counts = Counter(sc_gaps)
    all_gap_counts = Counter(all_gaps_in_range)
    print("\n  Gap distribution comparison (gap: sign_change_frac vs overall_frac):")
    for g in sorted(set(list(sc_gap_counts.keys())[:10] + list(all_gap_counts.keys())[:10])):
        if g == 0:
            continue
        sc_f = sc_gap_counts.get(g, 0) / max(len(sc_gaps), 1)
        all_f = all_gap_counts.get(g, 0) / max(len(all_gaps_in_range), 1)
        if all_f > 0.01 or sc_f > 0.01:
            ratio = sc_f / max(all_f, 1e-10)
            print(f"    gap={g:3d}: at_SC={sc_f:.4f}, overall={all_f:.4f}, ratio={ratio:.3f}")

    # Sign change spacing distribution
    sc_spacings = [sign_changes[i+1] - sign_changes[i] for i in range(len(sign_changes)-1)]
    if sc_spacings:
        print(f"\n  Sign change spacings: mean={sum(sc_spacings)/len(sc_spacings):.2f}, "
              f"min={min(sc_spacings)}, max={max(sc_spacings)}")
        sp_counts = Counter(sc_spacings)
        print("  Most common spacings:")
        for sp, cnt in sp_counts.most_common(10):
            print(f"    spacing={sp}: count={cnt} ({100*cnt/len(sc_spacings):.1f}%)")

    # ============================================================
    # TEST 3: PERIODICITY RELATED TO PRIMORIAL NUMBERS
    # ============================================================
    print("\n" + "=" * 70)
    print("[TEST 3] PERIODICITY vs PRIMORIAL NUMBERS")
    print("=" * 70)

    primorials = [2, 6, 30, 210, 2310, 30030]

    for P in primorials:
        if P > N_MAX // 2:
            break
        # Compute mean delta in each residue class mod P
        residue_sums = defaultdict(lambda: [0.0, 0])
        for n in range(2, N_MAX + 1):
            r = n % P
            residue_sums[r][0] += delta[n]
            residue_sums[r][1] += 1

        residue_means = {}
        for r in range(P):
            if residue_sums[r][1] > 0:
                residue_means[r] = residue_sums[r][0] / residue_sums[r][1]

        if not residue_means:
            continue

        mean_vals = list(residue_means.values())
        overall_mean = sum(mean_vals) / len(mean_vals)
        variance = sum((v - overall_mean)**2 for v in mean_vals) / len(mean_vals)
        std_of_means = math.sqrt(variance) if variance > 0 else 0

        # Compare std of residue means to expected if random
        # If delta were iid, std of means would be ~sigma/sqrt(count)
        delta_vals = [delta[n] for n in range(2, N_MAX + 1)]
        overall_std = math.sqrt(sum(d**2 for d in delta_vals) / len(delta_vals) -
                               (sum(delta_vals) / len(delta_vals))**2)
        avg_count = N_MAX / P
        expected_std = overall_std / math.sqrt(max(avg_count, 1))

        print(f"\n  Primorial P={P}:")
        print(f"    Std of residue-class means: {std_of_means:.6f}")
        print(f"    Expected if random:         {expected_std:.6f}")
        print(f"    Ratio (>1 means structure): {std_of_means / max(expected_std, 1e-15):.3f}")

        # Also check: does delta(n) correlate with n mod P?
        if P <= 30:
            # Show all residue means
            sorted_residues = sorted(residue_means.items())
            top5 = sorted(sorted_residues, key=lambda x: x[1])[:3]
            bot5 = sorted(sorted_residues, key=lambda x: -x[1])[:3]
            print(f"    Lowest  mean residues: {[(r, f'{m:.3f}') for r, m in top5]}")
            print(f"    Highest mean residues: {[(r, f'{m:.3f}') for r, m in bot5]}")

    # ============================================================
    # TEST 4: delta(n) mod SMALL NUMBERS
    # ============================================================
    print("\n" + "=" * 70)
    print("[TEST 4] STRUCTURE IN delta(n) MOD SMALL NUMBERS")
    print("=" * 70)

    # Since delta(n) is real-valued, we look at floor(delta(n)) mod m
    # and also round(delta(n)) mod m
    for m in [2, 3, 4, 5, 6, 7, 10, 12, 30]:
        rounded_mod = Counter()
        for n in range(2, N_MAX + 1):
            rd = int(round(delta[n])) % m
            rounded_mod[rd] += 1

        total = sum(rounded_mod.values())
        expected = total / m
        chi2 = sum((cnt - expected)**2 / expected for cnt in rounded_mod.values())

        # Chi-squared critical value for m-1 degrees of freedom at p=0.01
        # Approximate: for large df, chi2_crit ~ df + 2.33*sqrt(2*df)
        df = m - 1
        chi2_crit = df + 2.33 * math.sqrt(2 * df) if df > 1 else 6.63

        significant = "*** SIGNIFICANT ***" if chi2 > chi2_crit else ""
        print(f"  mod {m:2d}: chi2={chi2:10.2f}, crit={chi2_crit:.2f}, "
              f"df={df} {significant}")
        if chi2 > chi2_crit:
            dist = sorted(rounded_mod.items())
            print(f"         Distribution: {[(r, cnt) for r, cnt in dist]}")

    # Also check: is round(delta(n)) predominantly even or odd?
    even_count = sum(1 for n in range(2, N_MAX + 1) if int(round(delta[n])) % 2 == 0)
    odd_count = N_MAX - 1 - even_count
    print(f"\n  round(delta) parity: even={even_count}, odd={odd_count}, "
          f"ratio={even_count/max(odd_count,1):.4f}")

    # ============================================================
    # TEST 5: DISTRIBUTION ANALYSIS
    # ============================================================
    print("\n" + "=" * 70)
    print("[TEST 5] DISTRIBUTION OF NORMALIZED delta")
    print("=" * 70)

    # Normalization 1: delta(n) / sqrt(p(n))
    norm1 = []
    for n in range(2, N_MAX + 1):
        norm1.append(delta[n] / math.sqrt(p[n]))

    m1 = sum(norm1) / len(norm1)
    v1 = sum((x - m1)**2 for x in norm1) / len(norm1)
    s1 = math.sqrt(v1)
    skew1 = sum((x - m1)**3 for x in norm1) / (len(norm1) * s1**3) if s1 > 0 else 0
    kurt1 = sum((x - m1)**4 for x in norm1) / (len(norm1) * s1**4) - 3 if s1 > 0 else 0

    print(f"  delta/sqrt(p): mean={m1:.6f}, std={s1:.6f}, skew={skew1:.4f}, "
          f"excess_kurt={kurt1:.4f}")
    print(f"    (Gaussian: skew=0, kurt=0)")

    # Normalization 2: delta(n) / sqrt(p(n) * ln(p(n)))
    norm2 = []
    for n in range(2, N_MAX + 1):
        pn = p[n]
        norm2.append(delta[n] / math.sqrt(pn * math.log(pn)))

    m2 = sum(norm2) / len(norm2)
    v2 = sum((x - m2)**2 for x in norm2) / len(norm2)
    s2 = math.sqrt(v2)
    skew2 = sum((x - m2)**3 for x in norm2) / (len(norm2) * s2**3) if s2 > 0 else 0
    kurt2 = sum((x - m2)**4 for x in norm2) / (len(norm2) * s2**4) - 3 if s2 > 0 else 0

    print(f"\n  delta/sqrt(p*ln(p)): mean={m2:.6f}, std={s2:.6f}, skew={skew2:.4f}, "
          f"excess_kurt={kurt2:.4f}")

    # Normalization 3: (delta(n) - mean*sqrt(p)) / (std*sqrt(p))
    # Center and standardize
    std_norm = [(x - m1) / s1 for x in norm1]

    # Kolmogorov-Smirnov test vs Gaussian (manual)
    std_norm_sorted = sorted(std_norm)
    n_pts = len(std_norm_sorted)
    ks_stat = 0
    for i, x in enumerate(std_norm_sorted):
        # CDF of standard normal
        cdf = 0.5 * (1 + math.erf(x / math.sqrt(2)))
        ecdf = (i + 1) / n_pts
        ks_stat = max(ks_stat, abs(ecdf - cdf))

    # KS critical value at alpha=0.05: 1.36/sqrt(n)
    ks_crit = 1.36 / math.sqrt(n_pts)
    print(f"\n  KS test (standardized delta/sqrt(p) vs Gaussian):")
    print(f"    KS stat = {ks_stat:.6f}, critical = {ks_crit:.6f}")
    print(f"    {'REJECT Gaussian' if ks_stat > ks_crit else 'Cannot reject Gaussian'}")

    # Histogram of normalized delta
    print(f"\n  Histogram of delta/sqrt(p) (20 bins):")
    lo, hi = min(norm1), max(norm1)
    nbins = 20
    bin_width = (hi - lo) / nbins
    hist = [0] * nbins
    for x in norm1:
        b = min(int((x - lo) / bin_width), nbins - 1)
        hist[b] += 1
    for i in range(nbins):
        bar = "#" * (hist[i] * 50 // max(hist))
        center = lo + (i + 0.5) * bin_width
        print(f"    [{center:7.3f}] {hist[i]:5d} {bar}")

    # Rubinstein-Sarnak: The distribution should be related to
    # the logarithmic density of {x : pi(x) > li(x)}.
    # The key signature is a BIAS towards positive delta.
    pos_frac = sum(1 for n in range(2, N_MAX + 1) if delta[n] > 0) / (N_MAX - 1)
    print(f"\n  Fraction with delta > 0: {pos_frac:.6f}")
    print(f"  (Rubinstein-Sarnak predicts bias towards positive due to")
    print(f"   Chebyshev's bias in pi(x) - li(x))")

    # ============================================================
    # TEST 6: SECOND DIFFERENCE ANALYSIS
    # ============================================================
    print("\n" + "=" * 70)
    print("[TEST 6] SECOND DIFFERENCE: D2_delta(n) = delta(n+1)-2*delta(n)+delta(n-1)")
    print("=" * 70)

    d2 = [0.0] * (N_MAX + 1)
    for n in range(2, N_MAX):
        d2[n] = delta[n + 1] - 2 * delta[n] + delta[n - 1]

    # Note: D2_delta = (p(n+1)-2p(n)+p(n-1)) - (R^{-1}(n+1)-2R^{-1}(n)+R^{-1}(n-1))
    # The prime part is g(n) - g(n-1) where g(n) = p(n+1)-p(n)
    # The R^{-1} part is smooth, so D2_delta ~ gap_diff + smooth_correction

    d2_vals = [d2[n] for n in range(3, N_MAX)]
    d2_mean = sum(d2_vals) / len(d2_vals)
    d2_var = sum((x - d2_mean)**2 for x in d2_vals) / len(d2_vals)
    d2_std = math.sqrt(d2_var)

    print(f"  D2_delta: mean={d2_mean:.6f}, std={d2_std:.6f}")

    # Compare to gap differences
    gap_diffs = [gaps[n] - gaps[n - 1] for n in range(3, N_MAX)]
    gd_mean = sum(gap_diffs) / len(gap_diffs)
    gd_std = math.sqrt(sum((x - gd_mean)**2 for x in gap_diffs) / len(gap_diffs))
    print(f"  Gap diffs: mean={gd_mean:.6f}, std={gd_std:.6f}")
    print(f"  D2_delta is {'simpler' if d2_std < gd_std else 'NOT simpler'} than gap diffs")

    # Autocorrelation of D2
    d2c = [x - d2_mean for x in d2_vals]
    ac_d2 = [0.0] * 11
    norm_factor = sum(x**2 for x in d2c)
    for lag in range(1, 11):
        ac_d2[lag] = sum(d2c[i] * d2c[i + lag] for i in range(len(d2c) - lag)) / max(norm_factor, 1e-30)
    print(f"  Autocorrelation of D2_delta:")
    for lag in range(1, 11):
        print(f"    lag={lag:2d}: {ac_d2[lag]:.6f}")

    # D2 as integer: since p(n) are integers and R^{-1} is smooth,
    # round(D2_delta) should approximate gap_diff
    d2_as_int = [int(round(d2[n])) for n in range(3, N_MAX)]
    gap_diff_ints = [gaps[n] - gaps[n - 1] for n in range(3, N_MAX)]
    match_count = sum(1 for a, b in zip(d2_as_int, gap_diff_ints) if a == b)
    print(f"\n  round(D2_delta) == gap_diff: {match_count}/{len(d2_as_int)} "
          f"({100*match_count/len(d2_as_int):.2f}%)")

    # ============================================================
    # TEST 7: PATTERN OF NEAR-ZERO CROSSINGS
    # ============================================================
    print("\n" + "=" * 70)
    print("[TEST 7] PATTERN WHERE delta(n) ~ 0 (R^{-1} very accurate)")
    print("=" * 70)

    # Find n where |delta(n)| is small relative to sqrt(p(n))
    rel_delta = [(abs(delta[n]) / math.sqrt(p[n]), n) for n in range(2, N_MAX + 1)]
    rel_delta.sort()

    # Top 100 most accurate
    best = rel_delta[:100]
    best_ns = [n for _, n in best]

    print(f"  Top 20 most accurate R^{{-1}}(n) values:")
    print(f"  {'n':>7s} {'p(n)':>10s} {'R_inv':>14s} {'delta':>12s} "
          f"{'|delta|/sqrt(p)':>16s} {'gap_before':>10s}")
    for i in range(min(20, len(best))):
        rd, n = best[i]
        gb = gaps[n - 1] if n > 1 else 0
        print(f"  {n:7d} {p[n]:10d} {r_inv[n]:14.4f} {delta[n]:12.4f} "
              f"{rd:16.8f} {gb:10d}")

    # Check if these n values have special properties
    best_set = set(best_ns)

    # Are they related to prime gaps?
    best_gaps = [gaps[n - 1] for n in best_ns if n > 1]
    other_gaps = [gaps[n] for n in range(2, N_MAX) if n not in best_set]
    mean_best_gap = sum(best_gaps) / max(len(best_gaps), 1)
    mean_other_gap = sum(other_gaps) / max(len(other_gaps), 1)
    print(f"\n  Mean gap before best-100 n: {mean_best_gap:.4f}")
    print(f"  Mean gap for others:        {mean_other_gap:.4f}")

    # Residues mod small numbers
    print(f"\n  Residues of best-100 n mod small numbers:")
    for m in [2, 3, 6, 10, 30]:
        counts = Counter(n % m for n in best_ns)
        expected = 100 / m
        chi2 = sum((cnt - expected)**2 / expected for cnt in counts.values())
        print(f"    mod {m:2d}: {dict(sorted(counts.items()))} chi2={chi2:.2f}")

    # Spacing of best-100 n
    best_ns_sorted = sorted(best_ns)
    best_spacings = [best_ns_sorted[i+1] - best_ns_sorted[i]
                     for i in range(len(best_ns_sorted) - 1)]
    if best_spacings:
        print(f"\n  Spacings between best n: mean={sum(best_spacings)/len(best_spacings):.1f}, "
              f"min={min(best_spacings)}, max={max(best_spacings)}")

    # Are they near primorial multiples?
    print(f"\n  Distance to nearest primorial multiple:")
    for P in [6, 30, 210, 2310]:
        dists = [min(n % P, P - n % P) for n in best_ns]
        mean_dist = sum(dists) / len(dists)
        # Expected distance if uniform: P/4
        print(f"    P={P:5d}: mean_dist={mean_dist:.2f}, expected_if_random={P/4:.2f}, "
              f"ratio={mean_dist/(P/4):.3f}")

    # ============================================================
    # TEST 8 (BONUS): FOURIER/SPECTRAL ANALYSIS OF delta
    # ============================================================
    print("\n" + "=" * 70)
    print("[TEST 8] SPECTRAL ANALYSIS OF delta(n)")
    print("=" * 70)

    # Compute power spectrum using DFT of centered delta
    delta_vals = [delta[n] for n in range(2, min(N_MAX + 1, 20002))]  # use first 20k
    mean_d = sum(delta_vals) / len(delta_vals)
    centered = [d - mean_d for d in delta_vals]

    # Use manual DFT for a few key frequencies, then FFT-like approach
    N_fft = len(centered)

    # Check specific frequencies related to primorials
    print(f"\n  Power at primorial-related frequencies (N={N_fft}):")
    for P in [2, 6, 30, 210, 2310]:
        if P > N_fft // 2:
            break
        # Frequency = N_fft / P
        freq_idx = N_fft // P
        # Compute DFT at this frequency
        cos_sum = sum(centered[i] * math.cos(2 * math.pi * freq_idx * i / N_fft)
                      for i in range(N_fft))
        sin_sum = sum(centered[i] * math.sin(2 * math.pi * freq_idx * i / N_fft)
                      for i in range(N_fft))
        power = (cos_sum**2 + sin_sum**2) / N_fft**2

        # Compare to average power
        # Quick estimate of average power: var of centered / N
        avg_power = sum(c**2 for c in centered) / N_fft**2

        print(f"    Period {P:5d}: power={power:.6e}, avg_power={avg_power:.6e}, "
              f"ratio={power/max(avg_power,1e-30):.3f}")

    # Also check low frequencies (long-range trends)
    print(f"\n  Power at low frequencies (long periods):")
    for period in [100, 500, 1000, 5000, 10000]:
        freq_idx = N_fft // period
        if freq_idx == 0:
            continue
        cos_sum = sum(centered[i] * math.cos(2 * math.pi * freq_idx * i / N_fft)
                      for i in range(N_fft))
        sin_sum = sum(centered[i] * math.sin(2 * math.pi * freq_idx * i / N_fft)
                      for i in range(N_fft))
        power = (cos_sum**2 + sin_sum**2) / N_fft**2
        avg_power = sum(c**2 for c in centered) / N_fft**2
        print(f"    Period {period:5d}: power={power:.6e}, ratio={power/max(avg_power,1e-30):.3f}")

    # ============================================================
    # TEST 9 (BONUS): CORRELATION WITH PRIME INDEX PROPERTIES
    # ============================================================
    print("\n" + "=" * 70)
    print("[TEST 9] CORRELATION OF delta(n) WITH PROPERTIES OF n")
    print("=" * 70)

    # Is delta correlated with how "smooth" n is?
    def largest_prime_factor(n):
        if n <= 1:
            return 1
        f = 2
        lpf = 1
        while f * f <= n:
            while n % f == 0:
                lpf = f
                n //= f
            f += 1
        if n > 1:
            lpf = n
        return lpf

    def omega(n):
        """Number of distinct prime factors."""
        if n <= 1:
            return 0
        count = 0
        f = 2
        while f * f <= n:
            if n % f == 0:
                count += 1
                while n % f == 0:
                    n //= f
            f += 1
        if n > 1:
            count += 1
        return count

    # Sample to avoid slowness
    sample_size = 5000
    import random
    random.seed(42)
    sample_ns = sorted(random.sample(range(100, N_MAX + 1), sample_size))

    lpf_corr_data = [(largest_prime_factor(n), delta[n] / math.sqrt(p[n]))
                      for n in sample_ns]
    omega_corr_data = [(omega(n), delta[n] / math.sqrt(p[n]))
                        for n in sample_ns]

    # Group by omega
    by_omega = defaultdict(list)
    for om, d in omega_corr_data:
        by_omega[om].append(d)

    print(f"  delta/sqrt(p) grouped by omega(n) [number of prime factors of index]:")
    for om in sorted(by_omega.keys()):
        vals = by_omega[om]
        if len(vals) < 10:
            continue
        m = sum(vals) / len(vals)
        s = math.sqrt(sum((v - m)**2 for v in vals) / len(vals))
        print(f"    omega(n)={om}: count={len(vals):5d}, mean={m:.6f}, std={s:.6f}")

    # Pearson correlation between log(lpf(n)) and normalized delta
    lpf_vals = [math.log(max(d[0], 2)) for d in lpf_corr_data]
    delta_norm_vals = [d[1] for d in lpf_corr_data]
    m_lpf = sum(lpf_vals) / len(lpf_vals)
    m_dn = sum(delta_norm_vals) / len(delta_norm_vals)
    cov = sum((a - m_lpf) * (b - m_dn) for a, b in zip(lpf_vals, delta_norm_vals)) / len(lpf_vals)
    std_lpf = math.sqrt(sum((a - m_lpf)**2 for a in lpf_vals) / len(lpf_vals))
    std_dn = math.sqrt(sum((b - m_dn)**2 for b in delta_norm_vals) / len(delta_norm_vals))
    pearson = cov / (std_lpf * std_dn) if std_lpf > 0 and std_dn > 0 else 0
    print(f"\n  Pearson correlation(log(lpf(n)), delta/sqrt(p)): {pearson:.6f}")

    # ============================================================
    # TEST 10 (BONUS): RUNNING STATISTICS — HOW delta EVOLVES
    # ============================================================
    print("\n" + "=" * 70)
    print("[TEST 10] EVOLUTION OF delta(n) STATISTICS WITH n")
    print("=" * 70)

    windows = [1000, 2000, 5000, 10000, 20000, 30000, 40000, 50000]
    print(f"  {'n_max':>7s} {'mean_d/sqrt(p)':>15s} {'std_d/sqrt(p)':>14s} "
          f"{'skew':>8s} {'kurt':>8s} {'frac_pos':>9s}")
    for w in windows:
        if w > N_MAX:
            break
        vals = [delta[n] / math.sqrt(p[n]) for n in range(2, w + 1)]
        m = sum(vals) / len(vals)
        v = sum((x - m)**2 for x in vals) / len(vals)
        s = math.sqrt(v) if v > 0 else 1e-15
        sk = sum((x - m)**3 for x in vals) / (len(vals) * s**3) if s > 0 else 0
        ku = sum((x - m)**4 for x in vals) / (len(vals) * s**4) - 3 if s > 0 else 0
        fp = sum(1 for x in vals if x > 0) / len(vals)
        print(f"  {w:7d} {m:15.6f} {s:14.6f} {sk:8.4f} {ku:8.4f} {fp:9.6f}")

    # ============================================================
    # TEST 11 (BONUS): FIRST DIFFERENCE — delta(n+1) - delta(n)
    # ============================================================
    print("\n" + "=" * 70)
    print("[TEST 11] FIRST DIFFERENCE: D_delta(n) = delta(n+1) - delta(n)")
    print("=" * 70)

    d1 = [delta[n + 1] - delta[n] for n in range(2, N_MAX)]
    # D_delta = gap(n) - (R^{-1}(n+1) - R^{-1}(n))
    # R^{-1}(n+1) - R^{-1}(n) ~ ln(p(n)) for large n (since R'(x) ~ 1/ln(x))

    d1_mean = sum(d1) / len(d1)
    d1_std = math.sqrt(sum((x - d1_mean)**2 for x in d1) / len(d1))
    print(f"  D_delta: mean={d1_mean:.6f}, std={d1_std:.6f}")

    # What does D_delta look like? It's essentially gap(n) - expected_gap
    # where expected_gap = R^{-1}(n+1) - R^{-1}(n) ~ ln(p(n))
    predicted_gaps = [r_inv[n + 1] - r_inv[n] for n in range(2, N_MAX)]
    actual_gaps = [float(gaps[n]) for n in range(2, N_MAX)]
    gap_errors = [a - p for a, p in zip(actual_gaps, predicted_gaps)]

    ge_mean = sum(gap_errors) / len(gap_errors)
    ge_std = math.sqrt(sum((x - ge_mean)**2 for x in gap_errors) / len(gap_errors))
    print(f"  Gap prediction error (gap - (R_inv(n+1)-R_inv(n))): mean={ge_mean:.6f}, std={ge_std:.6f}")
    print(f"  Predicted gap mean: {sum(predicted_gaps)/len(predicted_gaps):.6f}")
    print(f"  Actual gap mean:    {sum(actual_gaps)/len(actual_gaps):.6f}")

    # ============================================================
    # SUMMARY
    # ============================================================
    print("\n" + "=" * 70)
    print("SUMMARY OF FINDINGS")
    print("=" * 70)
    print(f"""
  1. Sign changes: {len(sign_changes)} in [3,{N_MAX}]
     - Average spacing: {(N_MAX-2)/max(len(sign_changes),1):.1f}
     - Gap distribution at sign changes vs overall: see above

  2. Primorial periodicity: check ratio column above (>2 = significant)

  3. Modular structure: check chi2 values above

  4. Distribution:
     - delta/sqrt(p):        mean={m1:.6f}, std={s1:.6f}, skew={skew1:.4f}, kurt={kurt1:.4f}
     - delta/sqrt(p*ln(p)):  mean={m2:.6f}, std={s2:.6f}, skew={skew2:.4f}, kurt={kurt2:.4f}
     - KS test vs Gaussian:  stat={ks_stat:.6f}, crit={ks_crit:.6f}
     - Fraction positive:    {pos_frac:.6f}

  5. Second difference D2_delta: std={d2_std:.6f} vs gap_diff std={gd_std:.6f}

  6. Near-zero delta: see pattern analysis above

  7. Spectral: check ratio column (>>1 = peak at that frequency)

  KEY QUESTION: Is there exploitable structure for predicting delta(n)?
    """)

if __name__ == "__main__":
    main()
