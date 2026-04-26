"""
Random-subset vs first-K zero sampling for psi(x) recovery.

Fresh framing
-------------
Compressed sensing tells us: if a signal is K-sparse in some basis, then
~K log(N) random measurements suffice to recover it (vs N nyquist samples).

The explicit formula is a Fourier-type series in zeros gamma_k. Standard
practice: use the first K zeros (lowest frequencies). This is a *low-pass
filter*. But if the residual has hidden sparsity in the zero basis, then
RANDOM subsets of zeros may give better recovery than the deterministic
"first-K" choice (which suffers from sharp truncation Gibbs phenomena).

Question: Holding K constant, does the average residual differ between
  (a) using zeros[0:K]   (low-pass / contiguous prefix), vs
  (b) using a random size-K subset of zeros[0:M] for large M?
If (b) gives lower error: zeros have correlations exploitable for
  compressed-sensing recovery.
If (a) and (b) give same error: standard low-pass is optimal under the
  assumption of i.i.d. zero contributions; no structural compression.

We test on x in [10^6, 10^6 + 200) at K = 50, holding total pool size
M = 2000.
"""

import math
import numpy as np


def sieve_primes(N):
    sieve = np.ones(N + 1, dtype=bool)
    sieve[0] = sieve[1] = False
    for i in range(2, int(math.isqrt(N)) + 1):
        if sieve[i]:
            sieve[i * i :: i] = False
    return np.nonzero(sieve)[0]


def true_psi_at(N):
    primes = sieve_primes(N)
    psi = 0.0
    for p in primes.tolist():
        lp = math.log(p)
        pk = p
        while pk <= N:
            psi += lp
            pk *= p
    return psi


def load_zeros(path, n):
    zeros = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                zeros.append(float(line))
            except ValueError:
                continue
            if len(zeros) >= n:
                break
    return np.array(zeros, dtype=np.float64)


def psi_partial_subset(x, zeros_subset, scale=1.0):
    """Approximation summing only the zeros in `zeros_subset`, optionally
    rescaling the contribution by `scale` (used to compensate for partial
    sampling: if we sample fraction f of the pool, scale=1/f estimates
    the full sum)."""
    log_x = math.log(x)
    sqrt_x = math.sqrt(x)
    rho_abs = np.sqrt(0.25 + zeros_subset * zeros_subset)
    phi = np.arctan2(zeros_subset, 0.5)
    osc = np.cos(zeros_subset * log_x - phi) / rho_abs
    return x - 2.0 * sqrt_x * osc.sum() * scale - math.log(2.0 * math.pi)


def main():
    M = 2000
    K_first = 50
    K_rand = 50
    N_trials = 200  # number of x values
    N_random_subsets = 50  # random subsets to average over

    zeros = load_zeros(
        "/apps/aplikacijos/prime-research/data/zeta_zeros_2000.txt", M
    )
    print(f"Pool: {len(zeros)} zeros, gamma in [{zeros[0]:.2f}, {zeros[-1]:.2f}]")

    # Build incremental psi_true for x in [base, base + N_trials)
    base = 10**6
    primes_full = sieve_primes(base + N_trials + 10)
    prime_power_log = {}
    for p in primes_full.tolist():
        pk = p
        while pk <= base + N_trials + 10:
            prime_power_log[pk] = math.log(p)
            pk *= p
    psi_true_at_base = true_psi_at(base)

    psi_true = psi_true_at_base
    xs = []
    psi_trues = []
    for delta in range(0, N_trials):
        x = base + delta
        if delta > 0 and x in prime_power_log:
            psi_true += prime_power_log[x]
        xs.append(x)
        psi_trues.append(psi_true)
    psi_trues = np.array(psi_trues)

    # (a) First-K residuals
    zk_first = zeros[:K_first]
    res_first = np.array(
        [psi_partial_subset(x, zk_first) - pt for x, pt in zip(xs, psi_trues)]
    )

    # (b) Random-K residuals (averaged over many random subsets)
    rng = np.random.default_rng(0xC0DE)
    res_rand_collect = []
    for trial in range(N_random_subsets):
        idx = rng.choice(M, size=K_rand, replace=False)
        zk_rand = zeros[idx]
        # Scale: a random size-K subset estimates the full M-sum with factor M/K.
        # But the FULL sum is what we want to compare against truth, not
        # the FIRST-K sum. To compare apples to apples, do NOT rescale —
        # measure the residual of the unscaled sum, then look at the
        # *spread* across random subsets.
        r = np.array(
            [psi_partial_subset(x, zk_rand) - pt for x, pt in zip(xs, psi_trues)]
        )
        res_rand_collect.append(r)
    res_rand = np.array(res_rand_collect)  # shape (N_random_subsets, N_trials)

    # Statistics
    print("\n--- Residual statistics over x in [1e6, 1e6+200) ---")
    print(
        f"First-K=50:    mean={res_first.mean():+8.3f}  "
        f"std={res_first.std():.3f}  "
        f"max|res|={np.max(np.abs(res_first)):.3f}  "
        f"rms={math.sqrt(np.mean(res_first**2)):.3f}"
    )
    rand_means = res_rand.mean(axis=1)
    rand_stds = res_rand.std(axis=1)
    rand_rms = np.sqrt((res_rand**2).mean(axis=1))
    print(
        f"Random-K=50 (averaged over {N_random_subsets} subsets):\n"
        f"  mean of subset-mean:  {rand_means.mean():+8.3f}\n"
        f"  mean of subset-std:   {rand_stds.mean():.3f}\n"
        f"  mean of subset-rms:   {rand_rms.mean():.3f}\n"
        f"  best-subset rms:      {rand_rms.min():.3f}\n"
        f"  worst-subset rms:     {rand_rms.max():.3f}"
    )

    # Compare with FULL pool (M=2000) residual as reference
    res_full = np.array(
        [psi_partial_subset(x, zeros) - pt for x, pt in zip(xs, psi_trues)]
    )
    print(
        f"Full pool M=2000:  rms={math.sqrt(np.mean(res_full**2)):.3f}  "
        f"max|res|={np.max(np.abs(res_full)):.3f}"
    )

    # Test if averaging across random subsets reduces noise (variance reduction)
    avg_rand = res_rand.mean(axis=0)  # average over random subsets
    print(
        f"\nAveraging across {N_random_subsets} random subsets:\n"
        f"  rms of mean-over-subsets: {math.sqrt(np.mean(avg_rand**2)):.3f}"
    )
    print(
        "  (If this matches single-subset rms / sqrt(N_random_subsets), then "
        "subsets are independent estimators — no info gain over first-K).\n"
        "  (If it matches first-K rms, then random gives no advantage.)"
    )

    # Correlation between residual and gamma-spacing structure of pool
    # (test if "lucky" subsets are characterized by their gamma distribution)
    spreads = np.array([np.std(zeros[idx]) for idx in
                        [rng.choice(M, size=K_rand, replace=False)
                         for _ in range(50)]])
    print(f"\nSubset gamma-std varies in [{spreads.min():.1f}, {spreads.max():.1f}]")
    print("  (Consistency check on subset diversity.)\n")


if __name__ == "__main__":
    main()
