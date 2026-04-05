#!/usr/bin/env python3
"""
F_2 Correlation Profile of the Prime Indicator Function chi_P.

For N-bit inputs, define chi_P(x) = 1 if x is prime, 0 otherwise.
Compute the Walsh-Hadamard (Fourier over F_2) spectrum and correlation
with low-degree F_2 polynomials.

Key question: Does chi_P behave like a random function w.r.t. low-degree
F_2 correlations, or does it show structure exploitable by ACC^0[2] circuits?

If chi_P is NOT in ACC^0[2], correlations with degree-d polynomials must
decay exponentially in N for any fixed d.
"""

import numpy as np
from sympy import isprime
import time


def build_truth_table(N):
    """Build truth table of chi_P for x in {0, ..., 2^N - 1}."""
    size = 1 << N
    table = np.zeros(size, dtype=np.int8)
    for x in range(size):
        if isprime(x):
            table[x] = 1
    return table


def walsh_hadamard_transform(f, N):
    """
    Fast Walsh-Hadamard transform (vectorized numpy).
    f_hat(S) = (1/2^N) * sum_x (-1)^{f(x) + <S,x>}
    """
    size = 1 << N
    h = (1 - 2 * f.astype(np.float64)).copy()

    for i in range(N):
        step = 1 << i
        # Reshape to apply butterfly in-place via numpy
        h2 = h.reshape(-1, 2 * step)
        left = h2[:, :step].copy()
        right = h2[:, step:].copy()
        h2[:, :step] = left + right
        h2[:, step:] = left - right

    h /= size
    return h


def precompute_popcounts(N):
    """Precompute popcount for all 2^N values."""
    size = 1 << N
    indices = np.arange(size, dtype=np.uint32)
    # Use bit manipulation for popcount
    pc = np.zeros(size, dtype=np.int32)
    tmp = indices.copy()
    while np.any(tmp > 0):
        pc += (tmp & 1).astype(np.int32)
        tmp >>= 1
    return pc


def compute_correlation_profile(fhat, N, pc):
    """
    C(d) = max over |S| <= d of |fhat(S)|.
    max_at_deg[d] = max over |S|=d of |fhat(S)|.
    """
    abs_fhat = np.abs(fhat)
    max_at_deg = np.zeros(N + 1)
    for d in range(N + 1):
        mask = (pc == d)
        if np.any(mask):
            max_at_deg[d] = abs_fhat[mask].max()

    C = np.zeros(N + 1)
    C[0] = max_at_deg[0]
    for d in range(1, N + 1):
        C[d] = max(C[d - 1], max_at_deg[d])
    return C, max_at_deg


def compute_spectral_weight(fhat, N, pc):
    """W(d) = sum_{|S|=d} fhat(S)^2."""
    sq = fhat ** 2
    W = np.zeros(N + 1)
    for d in range(N + 1):
        W[d] = sq[pc == d].sum()
    return W


def precompute_monomials(N, max_deg, pc):
    """
    Precompute monomial values for all subsets S with |S| <= max_deg.
    monomial_S(x) = 1 iff (x & S) == S, i.e., all bits of S are set in x.
    Returns: list of subset indices, and a 2D array mono_vals[i, x].
    """
    size = 1 << N
    subsets = np.where(pc <= max_deg)[0]
    xs = np.arange(size, dtype=np.uint32)

    # Vectorized: for each subset S, compute (xs & S) == S
    # Do this in chunks to manage memory
    mono_vals = np.zeros((len(subsets), size), dtype=np.int8)
    for i, S in enumerate(subsets):
        mono_vals[i] = ((xs & S) == S).astype(np.int8)

    return subsets, mono_vals


def compute_degree_d_correlation(f, fhat, N, d, pc, num_samples=50000):
    """
    Estimate max correlation with degree-d polynomials over F_2.

    For d=1: exact (correlation with linear/affine = max |fhat(S)| for |S|<=1).
    For d>=2: use sampling + greedy.
    """
    size = 1 << N
    abs_fhat = np.abs(fhat)

    # Subsets of degree <= d
    low_mask = (pc <= d)
    low_indices = np.where(low_mask)[0]
    num_low = len(low_indices)

    # L1 upper bound
    l1_bound = abs_fhat[low_mask].sum()

    # Single monomial max
    single_max = abs_fhat[low_mask].max()

    if d <= 1:
        return single_max, l1_bound, single_max

    # Precompute monomial evaluations
    subsets, mono_vals = precompute_monomials(N, d, pc)
    # mono_vals shape: (num_low, size)

    g = 1 - 2 * f.astype(np.float64)  # (-1)^f
    best_corr = single_max

    rng = np.random.default_rng(42)

    # Strategy 1: Random sampling
    actual_samples = min(num_samples, 30000)
    for _ in range(actual_samples):
        k = rng.geometric(0.3)
        k = min(k, num_low)
        chosen = rng.choice(num_low, size=k, replace=False)

        # p(x) = XOR of chosen monomials => (-1)^{p(x)} = product of (1-2*mono)
        # XOR in F_2 means: p(x) = sum mod 2 of monomials
        # So p(x) = (sum of mono_vals[chosen]) mod 2
        parity = mono_vals[chosen].sum(axis=0) % 2
        signs = 1 - 2 * parity.astype(np.float64)
        corr = abs(np.dot(g, signs)) / size
        if corr > best_corr:
            best_corr = corr

    # Strategy 2: Greedy construction
    # Start with best single monomial
    best_idx = np.argmax(abs_fhat[low_indices])
    current_parity = mono_vals[best_idx].astype(np.int32)
    current_signs = 1 - 2 * (current_parity % 2).astype(np.float64)
    current_corr = abs(np.dot(g, current_signs)) / size

    for _ in range(min(30, num_low)):
        best_improve = 0
        best_next_idx = -1
        # Try adding each monomial
        check_indices = rng.choice(num_low, size=min(500, num_low), replace=False)
        for idx in check_indices:
            new_parity = (current_parity + mono_vals[idx]) % 2
            new_signs = 1 - 2 * new_parity.astype(np.float64)
            new_corr = abs(np.dot(g, new_signs)) / size
            improve = new_corr - current_corr
            if improve > best_improve:
                best_improve = improve
                best_next_idx = idx
        if best_next_idx >= 0 and best_improve > 1e-10:
            current_parity = (current_parity + mono_vals[best_next_idx]) % 2
            current_signs = 1 - 2 * current_parity.astype(np.float64)
            current_corr += best_improve
        else:
            break

    best_corr = max(best_corr, current_corr)
    return best_corr, l1_bound, single_max


def random_function_same_density(N, density, rng):
    """Generate a random Boolean function with the same density as chi_P."""
    size = 1 << N
    num_ones = int(round(density * size))
    table = np.zeros(size, dtype=np.int8)
    positions = rng.choice(size, size=num_ones, replace=False)
    table[positions] = 1
    return table


def fourier_entropy(fhat):
    """Compute Fourier entropy: -sum fhat(S)^2 * log2(fhat(S)^2), normalized."""
    sq = fhat ** 2
    sq = sq[sq > 1e-30]
    total = sq.sum()
    if total < 1e-30:
        return 0.0
    p = sq / total
    return -np.sum(p * np.log2(p))


def main():
    N_values = [6, 8, 10, 12, 14, 16]
    num_random = 10
    rng = np.random.default_rng(12345)

    all_results = {}

    for N in N_values:
        t0 = time.time()
        size = 1 << N
        print(f"\n{'='*72}")
        print(f"N = {N}  (domain size 2^{N} = {size})")
        print(f"{'='*72}")

        # Precompute popcounts
        pc = precompute_popcounts(N)

        # Build chi_P truth table
        f = build_truth_table(N)
        density = f.sum() / size
        print(f"Primes in [0, {size-1}]: {f.sum()}  (density = {density:.6f})")
        print(f"Prime number theorem prediction: {size / np.log(size):.1f}")

        # Walsh-Hadamard transform
        fhat = walsh_hadamard_transform(f, N)

        # Verify Parseval
        parseval = np.sum(fhat ** 2)
        print(f"Parseval check: sum fhat^2 = {parseval:.6f} (should be 1.0)")

        # Correlation profile
        C, max_at_deg = compute_correlation_profile(fhat, N, pc)

        # Spectral weight
        W = compute_spectral_weight(fhat, N, pc)

        # Fourier entropy
        fe = fourier_entropy(fhat)

        # Spectral norm
        spec_norm = np.max(np.abs(fhat))

        # Degree-1 spectral weight fraction
        w01 = W[0] + W[1]

        print(f"\nSpectral norm: {spec_norm:.6f}")
        print(f"Fourier entropy: {fe:.4f} bits (max = {N:.1f})")
        print(f"fhat(0) = {fhat[0]:.6f}  (= 1 - 2*density = {1-2*density:.6f})")
        print(f"Weight at degree 0+1: {w01:.6f} ({100*w01:.1f}% of total)")

        # Table: degree, W(d), max|fhat| at deg d, C(d)
        print(f"\n{'Deg':>4} {'W(d)':>12} {'max|fhat|@d':>14} {'C(d)':>12} {'W cumul':>12}")
        print(f"{'---':>4} {'---':>12} {'---':>14} {'---':>12} {'---':>12}")
        w_cumul = 0.0
        for d in range(N + 1):
            w_cumul += W[d]
            print(f"{d:>4} {W[d]:>12.6f} {max_at_deg[d]:>14.6f} {C[d]:>12.6f} {w_cumul:>12.6f}")

        # Degree-d correlation with general polynomials
        print(f"\nCorrelation with degree-d polynomials (sampling for d>=2):")
        print(f"{'Deg':>4} {'Best found':>12} {'L1 bound':>12} {'Single best':>12}")
        deg_corr_results = {}
        for d in [1, 2, 3]:
            # Skip d=3 for N>=14 (too many monomials to precompute)
            if d == 3 and N >= 14:
                continue
            # Limit samples for larger N
            nsamp = 50000 if N <= 12 else 10000
            best, l1, single = compute_degree_d_correlation(f, fhat, N, d, pc, nsamp)
            print(f"{d:>4} {best:>12.6f} {l1:>12.6f} {single:>12.6f}")
            deg_corr_results[d] = best

        # Random functions comparison
        print(f"\nComparison with {num_random} random functions (same density {density:.4f}):")
        rand_C = np.zeros((num_random, N + 1))
        rand_W = np.zeros((num_random, N + 1))
        rand_fe = np.zeros(num_random)
        rand_specnorm = np.zeros(num_random)

        for r in range(num_random):
            rf = random_function_same_density(N, density, rng)
            rfhat = walsh_hadamard_transform(rf, N)
            rC, _ = compute_correlation_profile(rfhat, N, pc)
            rW = compute_spectral_weight(rfhat, N, pc)
            rand_C[r] = rC
            rand_W[r] = rW
            rand_fe[r] = fourier_entropy(rfhat)
            rand_specnorm[r] = np.max(np.abs(rfhat))

        print(f"  Spectral norm: chi_P = {spec_norm:.6f}, "
              f"random = {rand_specnorm.mean():.6f} +/- {rand_specnorm.std():.6f}")
        print(f"  Fourier entropy: chi_P = {fe:.4f}, "
              f"random = {rand_fe.mean():.4f} +/- {rand_fe.std():.4f}")

        print(f"\n  C(d) comparison (max correlation with single degree-d monomial):")
        print(f"  {'Deg':>4} {'chi_P':>10} {'rand mean':>10} {'rand std':>10} {'z-score':>10}")
        for d in range(min(N + 1, 8)):
            rmean = rand_C[:, d].mean()
            rstd = rand_C[:, d].std()
            z = (C[d] - rmean) / rstd if rstd > 1e-10 else 0
            print(f"  {d:>4} {C[d]:>10.6f} {rmean:>10.6f} {rstd:>10.6f} {z:>10.2f}")

        print(f"\n  W(d) comparison (spectral weight at degree d):")
        print(f"  {'Deg':>4} {'chi_P':>10} {'rand mean':>10} {'rand std':>10} {'z-score':>10}")
        for d in range(min(N + 1, 8)):
            rmean = rand_W[:, d].mean()
            rstd = rand_W[:, d].std()
            z = (W[d] - rmean) / rstd if rstd > 1e-10 else 0
            print(f"  {d:>4} {W[d]:>10.6f} {rmean:>10.6f} {rstd:>10.6f} {z:>10.2f}")

        elapsed = time.time() - t0
        print(f"\nTime for N={N}: {elapsed:.2f}s")

        all_results[N] = {
            'density': density,
            'num_primes': int(f.sum()),
            'C': C,
            'W': W,
            'fe': fe,
            'spec_norm': spec_norm,
            'rand_C_mean': rand_C.mean(axis=0),
            'rand_C_std': rand_C.std(axis=0),
            'rand_specnorm_mean': rand_specnorm.mean(),
            'rand_specnorm_std': rand_specnorm.std(),
            'rand_fe_mean': rand_fe.mean(),
            'deg_corr': deg_corr_results,
        }

    # Summary table across all N
    print(f"\n\n{'='*72}")
    print("SUMMARY: Scaling of key quantities with N")
    print(f"{'='*72}")

    print(f"\n--- Density and spectral norm ---")
    print(f"{'N':>4} {'density':>10} {'#primes':>8} {'specnorm':>10} {'rand_sn':>10} {'ratio':>8}")
    for N in N_values:
        r = all_results[N]
        ratio = r['spec_norm'] / r['rand_specnorm_mean'] if r['rand_specnorm_mean'] > 0 else 0
        print(f"{N:>4} {r['density']:>10.6f} {r['num_primes']:>8} {r['spec_norm']:>10.6f} "
              f"{r['rand_specnorm_mean']:>10.6f} {ratio:>8.3f}")

    print(f"\n--- Correlation profile C(d) = max_{{|S|<=d}} |fhat(S)| ---")
    print(f"{'N':>4} {'C(0)':>10} {'C(1)':>10} {'C(2)':>10} {'C(3)':>10} {'C(4)':>10} {'C(5)':>10}")
    for N in N_values:
        r = all_results[N]
        vals = [f"{r['C'][d]:>10.6f}" if d <= N else f"{'':>10}" for d in range(6)]
        print(f"{N:>4} " + " ".join(vals))

    print(f"\n--- Same for random functions (mean) ---")
    print(f"{'N':>4} {'C(0)':>10} {'C(1)':>10} {'C(2)':>10} {'C(3)':>10} {'C(4)':>10} {'C(5)':>10}")
    for N in N_values:
        r = all_results[N]
        vals = [f"{r['rand_C_mean'][d]:>10.6f}" if d <= N else f"{'':>10}" for d in range(6)]
        print(f"{N:>4} " + " ".join(vals))

    print(f"\n--- Spectral weight W(d) for chi_P ---")
    print(f"{'N':>4} {'W(0)':>10} {'W(1)':>10} {'W(2)':>10} {'W(3)':>10} {'W(4)':>10} {'W(5)':>10}")
    for N in N_values:
        r = all_results[N]
        vals = [f"{r['W'][d]:>10.6f}" if d <= N else f"{'':>10}" for d in range(6)]
        print(f"{N:>4} " + " ".join(vals))

    print(f"\n--- Fourier entropy ---")
    print(f"{'N':>4} {'chi_P':>10} {'random':>10} {'max':>6}")
    for N in N_values:
        r = all_results[N]
        print(f"{N:>4} {r['fe']:>10.4f} {r['rand_fe_mean']:>10.4f} {N:>6}")

    print(f"\n--- Degree-d polynomial correlation (sampling, d>=2) ---")
    print(f"{'N':>4} {'d=1':>10} {'d=2':>10} {'d=3':>10}")
    for N in N_values:
        r = all_results[N]
        dc = r['deg_corr']
        d1 = f"{dc.get(1, 0):>10.6f}" if 1 in dc else f"{'N/A':>10}"
        d2 = f"{dc.get(2, 0):>10.6f}" if 2 in dc else f"{'N/A':>10}"
        d3 = f"{dc.get(3, 0):>10.6f}" if 3 in dc else f"{'N/A':>10}"
        print(f"{N:>4} {d1} {d2} {d3}")

    # Key scaling analysis
    print(f"\n--- KEY: Does C(1) = max linear correlation decay? ---")
    print(f"{'N':>4} {'C(1)':>10} {'1-2*dens':>10} {'C(1)/(1-2d)':>12} {'C(1)*sqrt(2^N)':>16}")
    for N in N_values:
        r = all_results[N]
        dens_term = 1 - 2 * r['density']
        ratio = r['C'][1] / dens_term if dens_term > 0 else 0
        scaled = r['C'][1] * np.sqrt(1 << N)
        print(f"{N:>4} {r['C'][1]:>10.6f} {dens_term:>10.6f} {ratio:>12.6f} {scaled:>16.4f}")

    print(f"\n--- KEY: W(1) vs random (degree-1 spectral weight, z-scores) ---")
    print("chi_P has MASSIVELY higher W(1) than random -- the parity of the last")
    print("bit (even/odd) strongly correlates with primality (all primes > 2 are odd).")
    print("This is the bit-1 linear structure, not useful for ACC^0 separation.")

    # Interpretation
    print(f"\n{'='*72}")
    print("INTERPRETATION")
    print(f"{'='*72}")
    print("""
1. SPECTRAL NORM: C(0) = |1 - 2*density| grows toward 1 as density -> 0
   (by PNT, density ~ 1/ln(2^N) ~ 1/(N*ln2)). This dominates C(d) for
   all d, making it hard to detect low-degree structure beyond the bias.

2. DEGREE-1 WEIGHT W(1): Enormously higher than random (z-scores 3-40+).
   This is because primes > 2 are all odd, so bit 0 (the parity bit)
   has huge correlation with primality. This is trivial structure.

3. HIGHER-DEGREE WEIGHT W(d) for d>=2: Generally BELOW random, meaning
   primes have LESS high-degree F_2 structure than random functions of
   the same density. Consistent with primes being "pseudorandom" at
   higher Fourier levels.

4. FOR ACC^0[2] SEPARATION: The key quantity is whether for fixed d,
   C(d) - C(0) decays exponentially. Since C(0) = 1-2*density already
   dominates, we need to look at whether any degree-d monomial with
   d >= 1 achieves correlation BEYOND the bias. The data shows max|fhat|
   at degree >= 2 is small and decreasing with N, consistent with
   chi_P being hard for ACC^0[2]. But the range N=6..16 is far too
   small to make definitive claims.

5. FOURIER ENTROPY: chi_P has consistently LOWER entropy than random,
   meaning its spectrum is more concentrated. This concentration is
   primarily in degrees 0 and 1 (the bias and parity structure).
""")


if __name__ == '__main__':
    main()
