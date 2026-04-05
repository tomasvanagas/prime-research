"""
Binary string structure and clustering patterns of delta(n) = p(n) - round(R^{-1}(n)).
Analyzes bit-level statistics, clustering, transition matrices, information-theoretic
measures, and run/pattern structure.
"""

import numpy as np
import sys
import os
from collections import Counter, defaultdict
from itertools import product
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from math import log2, log, factorial

OUTPUT_FILE = os.path.join(os.path.dirname(__file__), "binary_clustering_results.txt")
DELTA_FILE = "/apps/aplikacijos/prime-research/experiments/information_theory/kt_delta_values.npy"

out_lines = []

def pr(s=""):
    out_lines.append(s)
    print(s)

def section(title):
    pr()
    pr("=" * 72)
    pr(f"  {title}")
    pr("=" * 72)
    pr()


def to_signed_binary(val, width):
    """Convert signed integer to binary string of given width (two's complement)."""
    if val >= 0:
        return format(val, f'0{width}b')
    else:
        return format((1 << width) + val, f'0{width}b')


def bit_arrays(deltas, width):
    """Convert delta array to bit matrix (N x width), two's complement."""
    N = len(deltas)
    bits = np.zeros((N, width), dtype=np.int8)
    for i, d in enumerate(deltas):
        if d >= 0:
            v = d
        else:
            v = (1 << width) + d
        for b in range(width):
            bits[i, width - 1 - b] = (v >> b) & 1
    return bits


# ============================================================
# EXPERIMENT 1: Binary string analysis
# ============================================================
def experiment_binary(deltas):
    section("EXPERIMENT 1: BINARY STRING ANALYSIS")

    abs_max = max(abs(deltas.min()), abs(deltas.max()))
    width = int(np.ceil(np.log2(abs_max + 1))) + 1  # +1 for sign
    pr(f"Delta range: [{deltas.min()}, {deltas.max()}]")
    pr(f"Binary width (two's complement): {width} bits")
    pr()

    bits = bit_arrays(deltas, width)
    N = len(deltas)

    # 1a. Bit bias: fraction with bit set
    pr("--- Bit bias (fraction with bit=1) ---")
    bit_fracs = bits.mean(axis=0)
    for b in range(width):
        pr(f"  Bit {b:2d} (2^{width-1-b:2d}): {bit_fracs[b]:.6f}  "
           f"(bias from 0.5: {bit_fracs[b]-0.5:+.6f})")

    # Generate random integers with same mean/std for comparison
    rng = np.random.default_rng(42)
    rand_deltas = np.round(rng.normal(deltas.mean(), deltas.std(), N)).astype(int)
    rand_deltas = np.clip(rand_deltas, -(1 << (width - 1)), (1 << (width - 1)) - 1)
    rand_bits = bit_arrays(rand_deltas, width)
    rand_fracs = rand_bits.mean(axis=0)

    pr()
    pr("--- Comparison to Gaussian random with same mean/std ---")
    for b in range(width):
        pr(f"  Bit {b:2d}: delta={bit_fracs[b]:.6f}  random={rand_fracs[b]:.6f}  "
           f"diff={bit_fracs[b]-rand_fracs[b]:+.6f}")

    # 1b. Bit-level entropy per position
    pr()
    pr("--- Bit-level entropy per position (bits, max=1.0) ---")
    for b in range(width):
        p = bit_fracs[b]
        if p == 0 or p == 1:
            H = 0.0
        else:
            H = -p * log2(p) - (1 - p) * log2(1 - p)
        pr(f"  Bit {b:2d}: H = {H:.6f}")

    # 1c. Bit-pair correlations (Pearson)
    pr()
    pr("--- Bit-pair correlations (top 15 by |r|) ---")
    corrs = []
    for i in range(width):
        for j in range(i + 1, width):
            r = np.corrcoef(bits[:, i].astype(float), bits[:, j].astype(float))[0, 1]
            corrs.append((i, j, r))
    corrs.sort(key=lambda x: abs(x[2]), reverse=True)
    for i, j, r in corrs[:15]:
        pr(f"  Bits ({i:2d},{j:2d}): r = {r:+.6f}")

    # 1d. Joint entropy for adjacent bit pairs
    pr()
    pr("--- Joint entropy for adjacent bit pairs ---")
    for b in range(width - 1):
        joint = Counter()
        for row in bits:
            joint[(row[b], row[b + 1])] += 1
        H = 0.0
        for cnt in joint.values():
            p = cnt / N
            if p > 0:
                H -= p * log2(p)
        pr(f"  Bits ({b},{b+1}): H_joint = {H:.6f} bits (indep would be ~2.0)")

    return bits, width


# ============================================================
# EXPERIMENT 2: Delta clustering
# ============================================================
def experiment_clustering(deltas):
    section("EXPERIMENT 2: DELTA CLUSTERING")

    N = len(deltas)
    ns = np.arange(1, N + 1, dtype=float)

    # Normalize for clustering
    X = np.column_stack([
        (ns - ns.mean()) / ns.std(),
        (deltas - deltas.mean()) / deltas.std()
    ])

    for k in [2, 4, 8, 16]:
        km = KMeans(n_clusters=k, n_init=10, random_state=42, max_iter=300)
        labels = km.fit_predict(X)
        sil = silhouette_score(X, labels, sample_size=min(10000, N))
        pr(f"k={k:2d}: silhouette = {sil:.6f}")

        # Cluster sizes
        sizes = np.bincount(labels)
        pr(f"       cluster sizes: {sorted(sizes, reverse=True)}")

        # Check if clusters relate to n mod 6
        if k <= 8:
            pr(f"       n mod 6 distribution per cluster:")
            for c in range(k):
                mask = labels == c
                mods = (np.arange(1, N + 1)[mask]) % 6
                mod_counts = np.bincount(mods, minlength=6)
                mod_frac = mod_counts / mod_counts.sum()
                pr(f"         cluster {c}: " +
                   " ".join(f"{m}:{f:.3f}" for m, f in enumerate(mod_frac)))

        # Check n mod 30
        if k == 4:
            pr(f"       n mod 30 distribution per cluster (top 5 residues):")
            for c in range(k):
                mask = labels == c
                mods = (np.arange(1, N + 1)[mask]) % 30
                mod_counts = np.bincount(mods, minlength=30)
                top5 = np.argsort(mod_counts)[-5:][::-1]
                pr(f"         cluster {c}: " +
                   " ".join(f"{r}:{mod_counts[r]}" for r in top5))
        pr()

    # Chi-squared test: are cluster labels independent of n mod 6?
    pr("--- Chi-squared test: cluster labels (k=4) vs n mod 6 ---")
    km4 = KMeans(n_clusters=4, n_init=10, random_state=42)
    labels4 = km4.fit_predict(X)
    contingency = np.zeros((4, 6), dtype=int)
    for i in range(N):
        contingency[labels4[i], (i + 1) % 6] += 1
    chi2, pval, dof, _ = stats.chi2_contingency(contingency)
    pr(f"  chi2 = {chi2:.2f}, dof = {dof}, p-value = {pval:.6e}")
    pr()


# ============================================================
# EXPERIMENT 3: Transition patterns
# ============================================================
def experiment_transitions(deltas):
    section("EXPERIMENT 3: TRANSITION PATTERNS")

    N = len(deltas)
    n_bins = 30

    # Quantile-based binning
    quantiles = np.percentile(deltas, np.linspace(0, 100, n_bins + 1))
    quantiles[0] -= 1
    bin_indices = np.digitize(deltas, quantiles) - 1
    bin_indices = np.clip(bin_indices, 0, n_bins - 1)

    # Build transition matrix
    T = np.zeros((n_bins, n_bins), dtype=float)
    for i in range(N - 1):
        T[bin_indices[i], bin_indices[i + 1]] += 1

    # Normalize rows
    row_sums = T.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    P = T / row_sums

    # Check doubly stochastic
    col_sums = P.sum(axis=0)
    row_sums_p = P.sum(axis=1)
    pr(f"Transition matrix: {n_bins} x {n_bins} bins (quantile-based)")
    pr(f"Row sums: mean={row_sums_p.mean():.6f}, std={row_sums_p.std():.6f}")
    pr(f"Col sums: mean={col_sums.mean():.6f}, std={col_sums.std():.6f}, "
       f"min={col_sums.min():.6f}, max={col_sums.max():.6f}")
    pr(f"Is approximately doubly stochastic? "
       f"{'YES' if col_sums.std() < 0.05 else 'NO'} (col_sum std={col_sums.std():.6f})")

    # Eigenvalue analysis
    eigenvalues = np.linalg.eigvals(P)
    eigenvalues_sorted = sorted(np.abs(eigenvalues), reverse=True)
    pr()
    pr("--- Eigenvalue spectrum (top 10 by magnitude) ---")
    for i, ev in enumerate(eigenvalues_sorted[:10]):
        pr(f"  lambda_{i}: |lambda| = {ev:.6f}")

    spectral_gap = eigenvalues_sorted[0] - eigenvalues_sorted[1]
    pr(f"Spectral gap: {spectral_gap:.6f}")

    # Compare to random shuffle
    rng = np.random.default_rng(123)
    shuffled = deltas.copy()
    rng.shuffle(shuffled)
    bin_shuf = np.digitize(shuffled, quantiles) - 1
    bin_shuf = np.clip(bin_shuf, 0, n_bins - 1)
    T_shuf = np.zeros((n_bins, n_bins), dtype=float)
    for i in range(N - 1):
        T_shuf[bin_shuf[i], bin_shuf[i + 1]] += 1
    row_sums_shuf = T_shuf.sum(axis=1, keepdims=True)
    row_sums_shuf[row_sums_shuf == 0] = 1
    P_shuf = T_shuf / row_sums_shuf

    ev_shuf = sorted(np.abs(np.linalg.eigvals(P_shuf)), reverse=True)
    pr()
    pr("--- Shuffled (iid baseline) eigenvalues (top 5) ---")
    for i, ev in enumerate(ev_shuf[:5]):
        pr(f"  lambda_{i}: |lambda| = {ev:.6f}")
    pr(f"Shuffled spectral gap: {ev_shuf[0] - ev_shuf[1]:.6f}")

    # Frobenius distance between P and uniform
    uniform = np.ones_like(P) / n_bins
    frob_actual = np.linalg.norm(P - uniform, 'fro')
    frob_shuffled = np.linalg.norm(P_shuf - uniform, 'fro')
    pr()
    pr(f"Frobenius distance from uniform: actual={frob_actual:.6f}, "
       f"shuffled={frob_shuffled:.6f}")

    # Diagonal dominance check
    diag_mean = np.mean(np.diag(P))
    pr(f"Mean diagonal of P: {diag_mean:.6f} (uniform would be {1/n_bins:.6f})")
    pr(f"Diagonal enhancement ratio: {diag_mean * n_bins:.4f}x")


# ============================================================
# EXPERIMENT 4: Information-theoretic measures
# ============================================================
def lempel_ziv_complexity(seq):
    """Lempel-Ziv 76 complexity."""
    s = list(seq)
    n = len(s)
    i, c = 0, 1
    vocab = set()
    l = 1
    while i + l <= n:
        substr = tuple(s[i:i + l])
        if substr in vocab:
            l += 1
        else:
            vocab.add(substr)
            c += 1
            i = i + l
            l = 1
    return c


def approx_entropy(seq, m, r):
    """Approximate entropy (ApEn)."""
    N = len(seq)
    def phi(m_):
        patterns = defaultdict(int)
        for i in range(N - m_ + 1):
            pat = []
            for j in range(m_):
                pat.append(seq[i + j])
            patterns[tuple(pat)] += 1
        total = N - m_ + 1
        return sum((c / total) * log(c / total) for c in patterns.values())
    return phi(m) - phi(m + 1)


def sample_entropy(seq, m, r_frac):
    """Sample entropy using binned sequence."""
    N = len(seq)
    seq_arr = np.array(seq, dtype=float)
    r = r_frac * seq_arr.std()

    def count_matches(template_len):
        count = 0
        templates = []
        for i in range(N - template_len):
            templates.append(seq_arr[i:i + template_len])
        for i in range(len(templates)):
            for j in range(i + 1, len(templates)):
                if np.max(np.abs(templates[i] - templates[j])) <= r:
                    count += 1
        return count

    # For large N, subsample
    if N > 5000:
        idx = np.random.default_rng(42).choice(N - m - 1, size=5000, replace=False)
        sub = seq_arr[idx]
        r = r_frac * sub.std()
        # Use simplified version
        B = 0
        A = 0
        templates_m = [sub[i:i + m] if i + m <= len(sub) else None for i in range(len(sub))]
        # Too slow for full, use approximation via binning
        return _sampen_binned(seq, m, r_frac)

    B = count_matches(m)
    A = count_matches(m + 1)
    if B == 0:
        return float('inf')
    return -log(A / B)


def _sampen_binned(seq, m, r_frac):
    """Approximate SampEn by discretizing."""
    arr = np.array(seq, dtype=float)
    std = arr.std()
    r = r_frac * std
    # Bin the sequence
    binned = np.round(arr / r).astype(int)

    N = len(binned)
    # Count template matches for length m and m+1
    count_m = defaultdict(int)
    count_m1 = defaultdict(int)
    for i in range(N - m):
        count_m[tuple(binned[i:i + m])] += 1
        count_m1[tuple(binned[i:i + m + 1])] += 1

    B = sum(c * (c - 1) // 2 for c in count_m.values())
    A = sum(c * (c - 1) // 2 for c in count_m1.values())

    if B == 0:
        return float('inf')
    return -log(A / B)


def permutation_entropy(seq, order):
    """Permutation entropy."""
    N = len(seq)
    perm_counts = defaultdict(int)
    for i in range(N - order + 1):
        window = seq[i:i + order]
        perm = tuple(np.argsort(window))
        perm_counts[perm] += 1
    total = N - order + 1
    H = 0.0
    for cnt in perm_counts.values():
        p = cnt / total
        if p > 0:
            H -= p * log2(p)
    max_H = log2(factorial(order))
    return H, max_H, len(perm_counts), factorial(order)


def experiment_information(deltas):
    section("EXPERIMENT 4: INFORMATION-THEORETIC MEASURES")

    N = len(deltas)

    # 4a. Lempel-Ziv complexity (on sign sequence for tractability)
    pr("--- Lempel-Ziv complexity ---")
    sign_seq = np.sign(deltas)
    lz_sign = lempel_ziv_complexity(sign_seq[:50000])
    # Random baseline
    rng = np.random.default_rng(42)
    rand_sign = rng.choice([-1, 0, 1], size=50000,
                           p=[np.mean(sign_seq < 0), np.mean(sign_seq == 0), np.mean(sign_seq > 0)])
    lz_rand = lempel_ziv_complexity(rand_sign)
    pr(f"  Sign sequence (N=50000): LZ = {lz_sign}")
    pr(f"  Random baseline:         LZ = {lz_rand}")
    pr(f"  Ratio (delta/random):    {lz_sign/lz_rand:.6f}")

    # LZ on quantized delta (16 bins)
    quantiles16 = np.percentile(deltas, np.linspace(0, 100, 17))
    quantiles16[0] -= 1
    binned16 = np.digitize(deltas, quantiles16) - 1
    binned16 = np.clip(binned16, 0, 15)
    lz_q16 = lempel_ziv_complexity(binned16[:50000])
    rand_q16 = rng.integers(0, 16, size=50000)
    lz_rand16 = lempel_ziv_complexity(rand_q16)
    pr(f"  Quantized (16 bins, N=50000): LZ = {lz_q16}")
    pr(f"  Random 16-symbol baseline:    LZ = {lz_rand16}")
    pr(f"  Ratio: {lz_q16/lz_rand16:.6f}")
    pr()

    # 4b. Approximate entropy
    pr("--- Approximate Entropy (ApEn) ---")
    # Use binned sequence for speed
    for m in [2, 3]:
        apen = approx_entropy(binned16[:20000], m, 0)
        rand_apen = approx_entropy(rng.integers(0, 16, size=20000), m, 0)
        pr(f"  m={m}: ApEn(delta) = {apen:.6f},  ApEn(random) = {rand_apen:.6f}")
    pr()

    # 4c. Sample entropy (binned approximation)
    pr("--- Sample Entropy (SampEn, binned approx) ---")
    for m in [2, 3]:
        for r_frac in [0.2]:
            se = _sampen_binned(deltas[:20000].tolist(), m, r_frac)
            se_rand = _sampen_binned(
                rng.normal(deltas.mean(), deltas.std(), 20000).tolist(), m, r_frac)
            pr(f"  m={m}, r=0.2*std: SampEn(delta) = {se:.6f}, SampEn(random) = {se_rand:.6f}")
    pr()

    # 4d. Permutation entropy
    pr("--- Permutation Entropy ---")
    delta_list = deltas[:50000].tolist()
    for order in [3, 4, 5, 6, 7]:
        H, H_max, n_perms, total_perms = permutation_entropy(delta_list, order)
        pr(f"  Order {order}: H = {H:.6f} / {H_max:.6f} (normalized: {H/H_max:.6f}), "
           f"observed {n_perms}/{total_perms} permutations")


# ============================================================
# EXPERIMENT 5: Runs and patterns
# ============================================================
def experiment_runs(deltas):
    section("EXPERIMENT 5: RUNS AND PATTERNS")

    N = len(deltas)
    signs = np.sign(deltas)

    # 5a. Run lengths (consecutive same-sign)
    pr("--- Run length distribution (same sign) ---")
    runs = []
    current_sign = signs[0]
    current_len = 1
    for i in range(1, N):
        if signs[i] == current_sign:
            current_len += 1
        else:
            if current_sign != 0:  # skip zeros
                runs.append(current_len)
            current_sign = signs[i]
            current_len = 1
    if current_sign != 0:
        runs.append(current_len)

    runs = np.array(runs)
    pr(f"  Total runs: {len(runs)}")
    pr(f"  Mean run length: {runs.mean():.4f}")
    pr(f"  Std run length: {runs.std():.4f}")
    pr(f"  Max run length: {runs.max()}")
    pr(f"  Median run length: {np.median(runs):.1f}")

    # Distribution of run lengths
    run_counts = Counter(runs)
    pr()
    pr("  Run length distribution (top 15):")
    for length in sorted(run_counts.keys())[:15]:
        cnt = run_counts[length]
        pr(f"    length {length:3d}: {cnt:6d} ({cnt/len(runs)*100:.2f}%)")

    # Compare to geometric distribution
    p_switch = 1.0 / runs.mean()
    pr()
    pr(f"  Geometric fit: p = {p_switch:.6f}")
    pr("  Comparison (observed vs geometric):")
    for length in range(1, 11):
        obs = run_counts.get(length, 0) / len(runs)
        geo = p_switch * (1 - p_switch) ** (length - 1)
        pr(f"    length {length}: observed={obs:.6f}, geometric={geo:.6f}, "
           f"ratio={obs/geo:.4f}" if geo > 0 else f"    length {length}: obs={obs:.6f}")

    # KS test against geometric
    from scipy.stats import kstest
    geo_cdf = lambda x: 1 - (1 - p_switch) ** np.floor(x)
    ks_stat, ks_p = kstest(runs, geo_cdf)
    pr(f"  KS test vs geometric: stat={ks_stat:.6f}, p-value={ks_p:.6e}")

    # 5b. Sign trigrams and 4-grams
    pr()
    pr("--- Sign n-gram analysis ---")

    # Map signs to {-, 0, +}
    sign_chars = {-1: '-', 0: '0', 1: '+'}
    sign_str = [sign_chars[s] for s in signs]

    for gram_len in [3, 4]:
        pr(f"\n  {gram_len}-grams:")
        gram_counts = Counter()
        for i in range(N - gram_len + 1):
            gram = tuple(sign_str[i:i + gram_len])
            gram_counts[gram] += 1

        total = N - gram_len + 1
        all_possible = list(product(['-', '0', '+'], repeat=gram_len))
        observed = set(gram_counts.keys())
        missing = [g for g in all_possible if g not in observed]

        pr(f"    Total possible: {len(all_possible)}")
        pr(f"    Observed: {len(observed)}")
        pr(f"    Missing: {len(missing)}")

        if missing and len(missing) <= 20:
            for g in missing:
                pr(f"      Missing: {''.join(g)}")

        # Top and bottom by frequency
        sorted_grams = gram_counts.most_common()
        pr(f"    Most common:")
        for gram, cnt in sorted_grams[:5]:
            pr(f"      {''.join(gram)}: {cnt} ({cnt/total*100:.3f}%)")
        pr(f"    Least common:")
        for gram, cnt in sorted_grams[-5:]:
            pr(f"      {''.join(gram)}: {cnt} ({cnt/total*100:.3f}%)")

    # 5c. Delta-value trigrams (binned)
    pr()
    pr("--- Delta value n-gram forbidden patterns (4 bins: Q1/Q2/Q3/Q4) ---")
    quartiles = np.percentile(deltas, [25, 50, 75])
    qbin = np.digitize(deltas, quartiles)  # 0,1,2,3

    for gram_len in [3, 4]:
        gram_counts = Counter()
        for i in range(N - gram_len + 1):
            gram = tuple(qbin[i:i + gram_len])
            gram_counts[gram] += 1

        total_possible = 4 ** gram_len
        observed = len(gram_counts)
        expected_per = (N - gram_len + 1) / total_possible
        pr(f"  {gram_len}-grams: {observed}/{total_possible} observed "
           f"(expected ~{expected_per:.1f} per pattern)")

        # Chi-squared uniformity test
        counts_arr = np.array([gram_counts.get(g, 0)
                               for g in product(range(4), repeat=gram_len)])
        expected_arr = np.full_like(counts_arr, dtype=float,
                                    fill_value=counts_arr.sum() / len(counts_arr))
        chi2 = np.sum((counts_arr - expected_arr) ** 2 / expected_arr)
        dof = len(counts_arr) - 1
        from scipy.stats import chi2 as chi2_dist
        p_val = 1 - chi2_dist.cdf(chi2, dof)
        pr(f"    Chi-squared uniformity: chi2={chi2:.2f}, dof={dof}, p={p_val:.6e}")


# ============================================================
# MAIN
# ============================================================
def main():
    pr("BINARY STRING STRUCTURE AND CLUSTERING ANALYSIS OF DELTA(n)")
    pr(f"delta(n) = p(n) - round(R^{{-1}}(n))")
    pr(f"Data: {DELTA_FILE}")

    deltas = np.load(DELTA_FILE)
    N = len(deltas)
    pr(f"N = {N}")
    pr(f"Range: [{deltas.min()}, {deltas.max()}]")
    pr(f"Mean: {deltas.mean():.4f}, Std: {deltas.std():.4f}")
    pr(f"Median: {np.median(deltas):.1f}")
    pr(f"Fraction positive: {np.mean(deltas > 0):.6f}")
    pr(f"Fraction zero: {np.mean(deltas == 0):.6f}")
    pr(f"Fraction negative: {np.mean(deltas < 0):.6f}")

    experiment_binary(deltas)
    experiment_clustering(deltas)
    experiment_transitions(deltas)
    experiment_information(deltas)
    experiment_runs(deltas)

    section("SUMMARY")
    pr("See individual sections above for detailed findings.")

    # Write results
    with open(OUTPUT_FILE, 'w') as f:
        f.write('\n'.join(out_lines) + '\n')
    pr()
    pr(f"Results saved to {OUTPUT_FILE}")


if __name__ == '__main__':
    main()
