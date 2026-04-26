"""F1: Persistent Homology of Detrended Prime Gap Sequence.

Conjecture (TCDP): The number of "long-lived" 0-dim persistence bars in the
Vietoris-Rips filtration of {(n/N, d_n) : n=1..N} is o(N / log N), where
d_n = p_n - n*ln(n) - n*(ln ln n - 1) is the detrended prime sequence.

A "long-lived" bar = persists beyond the noise floor (gap > 95th percentile).

If TCDP holds, primes carry a polylog topological signature; if not, all bars
are noise-scale and TDA gives no compression.

We use 0-dim persistence (== single-linkage dendrogram), which is computable
exactly via scipy without a dedicated TDA library.

We compare against:
- Cramér random primes: each integer k included independently w.p. 1/log(k).
- True primes.
If TCDP-bar-count(true) << TCDP-bar-count(Cramér), TDA is detecting structure.
"""

import sys
import time

import numpy as np
from sympy import primerange
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist


def detrend(primes):
    """Return d_n = p_n - n ln n - n (ln ln n - 1)."""
    n = np.arange(1, len(primes) + 1, dtype=np.float64)
    trend = n * np.log(n + 1) + n * (np.log(np.log(n + 2)) - 1)
    return np.array(primes, dtype=np.float64) - trend


def cramer_pseudoprimes(N_target, seed=0):
    """Generate ~N_target Cramér-model primes."""
    rng = np.random.default_rng(seed)
    # We need integers k such that random < 1/log(k); collect ~N_target
    # Estimate range needed: nth prime ~ n log n, so go up to ~2 * N_target * log(N_target)
    upper = int(2 * N_target * max(np.log(N_target), 2))
    ks = np.arange(2, upper + 1)
    probs = 1.0 / np.log(ks.astype(np.float64))
    keep = rng.random(len(ks)) < probs
    pseudo = ks[keep][:N_target]
    return pseudo.tolist()


def persistence_intervals(cloud):
    """Compute 0-dim persistence intervals via single-linkage hierarchical clustering.

    Returns sorted list of bar lengths (deaths). One bar (the infinite one) is dropped.
    """
    if len(cloud) < 2:
        return np.array([])
    D = pdist(cloud)
    Z = linkage(D, method="single")
    # Z[:,2] = merge distances. Each merge kills one connected component.
    # Bar lengths = merge distances (births are all 0 for 0-dim).
    deaths = Z[:, 2].copy()
    return np.sort(deaths)


def count_long_lived(deaths, noise_threshold_q=0.50):
    """Count bars longer than the q-quantile of all bar lengths.

    A "long-lived" bar is one above the median by default. The number of such
    bars relative to N is the relevant quantity for TCDP.
    """
    if len(deaths) == 0:
        return 0
    threshold = np.quantile(deaths, noise_threshold_q)
    return int(np.sum(deaths > 5 * threshold))  # 5x median = clearly above noise


def main():
    Ns = [200, 500, 1000, 2000]
    print(f"=== F1: Persistent Homology of Detrended Primes ===\n")

    for N in Ns:
        primes = list(primerange(2, 100000))[:N]
        d_true = detrend(primes)

        # Build cloud: normalize n to [0, 1] and d to unit std
        n_norm = np.linspace(0, 1, N)
        d_norm = (d_true - d_true.mean()) / (d_true.std() + 1e-12)
        cloud_true = np.column_stack([n_norm, d_norm])

        deaths_true = persistence_intervals(cloud_true)
        long_true = count_long_lived(deaths_true)
        med_true = np.median(deaths_true) if len(deaths_true) else 0
        max_true = deaths_true[-1] if len(deaths_true) else 0
        log_threshold = N / max(np.log(N), 1.0)

        # Cramér baseline
        pseudo = cramer_pseudoprimes(N, seed=42)
        d_pseudo = detrend(pseudo)
        d_pseudo_norm = (d_pseudo - d_pseudo.mean()) / (d_pseudo.std() + 1e-12)
        cloud_pseudo = np.column_stack([n_norm, d_pseudo_norm])
        deaths_pseudo = persistence_intervals(cloud_pseudo)
        long_pseudo = count_long_lived(deaths_pseudo)

        # IID uniform baseline (pure noise)
        rng = np.random.default_rng(7)
        d_iid = rng.standard_normal(N)
        cloud_iid = np.column_stack([n_norm, d_iid])
        deaths_iid = persistence_intervals(cloud_iid)
        long_iid = count_long_lived(deaths_iid)

        print(f"N = {N}")
        print(f"  True primes:   {long_true:5d} long-lived bars  "
              f"(median={med_true:.4f}, max={max_true:.4f})")
        print(f"  Cramér model:  {long_pseudo:5d} long-lived bars")
        print(f"  IID Gaussian:  {long_iid:5d} long-lived bars")
        print(f"  N/log N target: {log_threshold:.1f}")
        ratio = long_true / log_threshold if log_threshold > 0 else float("inf")
        print(f"  long_true / (N/log N) = {ratio:.3f}")
        if long_true < long_pseudo / 2:
            print(f"  *** True primes have FEWER bars than Cramér => possible structure! ***")
        elif long_true > 2 * long_pseudo:
            print(f"  True primes have MORE bars than Cramér => anti-structure")
        else:
            print(f"  True primes ~ Cramér => no detectable topological compression")
        print()

    # Scaling test
    print(f"\n=== Scaling: long_true(N) growth ===")
    Ns_big = [200, 400, 800, 1600]
    longs = []
    for N in Ns_big:
        primes = list(primerange(2, 100000))[:N]
        d_true = detrend(primes)
        n_norm = np.linspace(0, 1, N)
        d_norm = (d_true - d_true.mean()) / (d_true.std() + 1e-12)
        cloud = np.column_stack([n_norm, d_norm])
        deaths = persistence_intervals(cloud)
        longs.append(count_long_lived(deaths))
        print(f"  N={N:5d}: long_lived={longs[-1]:5d}, ratio to N={longs[-1]/N:.4f}, "
              f"ratio to N/log N={longs[-1]/(N/np.log(N)):.4f}")

    # Fit growth: log(long) = alpha * log(N) + beta
    if len(longs) >= 3 and all(x > 0 for x in longs):
        ln_n = np.log(Ns_big)
        ln_l = np.log(longs)
        alpha, beta = np.polyfit(ln_n, ln_l, 1)
        print(f"\n  Fit: long(N) ~ N^{alpha:.3f} * exp({beta:.3f})")
        if alpha < 0.7:
            print(f"  *** sub-linear scaling => possible compression ***")
        elif alpha > 0.95:
            print(f"  ~linear scaling => no compression (every prime is its own bar)")


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"\nTotal time: {time.time() - t0:.2f} s")
