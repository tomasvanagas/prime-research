"""
Hankel singular-value test for hidden linear-recurrence structure in chi_P.

Idea (fresh perspective, session 58):
   If the prime indicator chi_P(n) satisfies *any* finite-order linear
   recurrence over the reals (or, by extension, any near-low-rank
   approximate recurrence), then the Hankel matrix
        H[i, j]  =  chi_P(i + j)         (0 <= i, j < K)
   has rank r equal to the order of that recurrence (Berlekamp-Massey).

   The prime indicator is *not* expected to obey such a recurrence -- but
   that's also what we expect for a random Bernoulli sequence with the
   same density.  The interesting question is whether real primes are
   either:
     (a) *more* compressible than random (lower effective rank, hidden
         recurrence-like structure), or
     (b) *less* compressible than random (anti-correlated, somehow more
         "spread out" than chance), or
     (c) statistically indistinguishable.

   We compute singular values of H_chi and H_random, look at the decay,
   and report effective ranks at several thresholds.

   This is a different basis than DCT, wavelet, or Walsh -- it directly
   tests for linear-recurrence structure (hidden algebraic generators).

Pass / fail
   FAIL  -> SV spectra of H_chi and H_random are statistically
            indistinguishable (effective ranks within ~1 sigma).
            Adds a 23rd pseudorandomness measure.
   WEAK  -> H_chi has effective rank ~10x lower than H_random, suggesting
            a partial low-rank approximation.
   PASS  -> H_chi has effective rank polylog(K) << K (so primes obey a
            short linear recurrence -- this would be a major result).
"""

import numpy as np


def sieve_to(N: int) -> np.ndarray:
    s = np.ones(N + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(N ** 0.5) + 1):
        if s[p]:
            s[p * p :: p] = False
    return s


def hankel_from_seq(seq: np.ndarray, K: int) -> np.ndarray:
    """H[i, j] = seq[i + j] for 0 <= i, j < K. Need len(seq) >= 2*K - 1."""
    H = np.empty((K, K), dtype=np.float64)
    for i in range(K):
        H[i] = seq[i : i + K]
    return H


def effective_rank(svals: np.ndarray, threshold: float) -> int:
    """Smallest r such that ||sigma_{r:}||^2 / ||sigma||^2 < threshold."""
    s2 = svals ** 2
    total = s2.sum()
    if total == 0:
        return 0
    cum = np.cumsum(s2[::-1])[::-1] / total  # tail energy
    r = int(np.argmax(cum < threshold))
    if cum[-1] >= threshold:
        return len(svals)
    return r


def main():
    K = 600
    N = 2 * K + 50
    print(f"# Hankel SV test for chi_P, K={K}, N={N}")

    chi = sieve_to(N).astype(np.float64)
    print(f"# pi({N}) = {int(chi.sum())}, density = {chi.mean():.4f}")

    H = hankel_from_seq(chi, K)
    print("# computing SVD of H_chi...")
    s_chi = np.linalg.svd(H, compute_uv=False)

    # Naive control: Bernoulli with same density as chi.
    p = float(chi.mean())
    rs_naive = []
    for seed in range(5):
        rng = np.random.default_rng(seed + 100)
        rand_seq = (rng.random(len(chi)) < p).astype(np.float64)
        H_r = hankel_from_seq(rand_seq, K)
        rs_naive.append(np.linalg.svd(H_r, compute_uv=False))
    s_naive_mean = np.mean(rs_naive, axis=0)
    s_naive_std = np.std(rs_naive, axis=0)

    # Sieve-aware control: Bernoulli on integers coprime to 30.
    # Density on the skeleton matches chi's density on that skeleton.
    skeleton = np.array([np.gcd(n, 30) == 1 for n in range(len(chi))], dtype=bool)
    skeleton[2] = skeleton[3] = skeleton[5] = True  # tiny primes always allowed
    chi_on_skel = chi.copy()
    chi_on_skel[~skeleton] = 0.0
    p_skel = float(chi[skeleton].mean()) if skeleton.sum() > 0 else 0.0
    rs_sieve = []
    for seed in range(5):
        rng = np.random.default_rng(seed + 200)
        rand_seq = (rng.random(len(chi)) < p_skel).astype(np.float64)
        rand_seq[~skeleton] = 0.0
        # Force 2,3,5 to mimic chi exactly (trivial small-prime contribution).
        for q in (2, 3, 5):
            if q < len(rand_seq):
                rand_seq[q] = 1.0
        H_r = hankel_from_seq(rand_seq, K)
        rs_sieve.append(np.linalg.svd(H_r, compute_uv=False))
    s_sieve_mean = np.mean(rs_sieve, axis=0)
    s_sieve_std = np.std(rs_sieve, axis=0)
    print(f"# skeleton density = {skeleton.mean():.4f}, "
          f"p on skel = {p_skel:.4f}")

    # Compare effective ranks at thresholds.
    print()
    print("# effective rank at tail-energy threshold")
    print("# threshold   chi_P     naive_mean     sieve_mean")
    for thresh in [1e-1, 1e-2, 1e-3, 1e-4, 1e-6]:
        r_chi = effective_rank(s_chi, thresh)
        rn = float(np.mean([effective_rank(s, thresh) for s in rs_naive]))
        rs_ = float(np.mean([effective_rank(s, thresh) for s in rs_sieve]))
        print(f"#  {thresh:.0e}    {r_chi:5d}     {rn:7.1f}        {rs_:7.1f}")

    # Direct comparison of leading SV ratios vs sieve-aware control.
    print()
    print("# leading singular values: chi_P vs SIEVE-AWARE control")
    print("#  k     s_chi      s_sieve_mean    z_sieve")
    for k in [0, 1, 2, 3, 4, 5, 6, 7, 8, 16, 32, 64, 128, 256]:
        if k < len(s_chi):
            z = (s_chi[k] - s_sieve_mean[k]) / max(s_sieve_std[k], 1e-9)
            print(f"#  {k:4d}  {s_chi[k]:.4e}   {s_sieve_mean[k]:.4e}   {z:+.3f}")

    # Save data.
    np.savez(
        "/apps/aplikacijos/prime-research/experiments/wildcard/hankel_chi_P_session58_data.npz",
        s_chi=s_chi,
        s_naive_mean=s_naive_mean,
        s_naive_std=s_naive_std,
        s_sieve_mean=s_sieve_mean,
        s_sieve_std=s_sieve_std,
    )

    # Verdict relative to sieve control.
    r1_chi = effective_rank(s_chi, 1e-2)
    r1_sieve = float(np.mean([effective_rank(s, 1e-2) for s in rs_sieve]))

    # Check max |z_sieve| in mid-range k = [10, 100] (above leading low-rank
    # structure, below noise floor) -- this is the meaningful regime.
    z_mid = [(s_chi[k] - s_sieve_mean[k]) / max(s_sieve_std[k], 1e-9)
             for k in range(10, 100)]
    z_max = float(np.max(np.abs(z_mid)))
    print()
    print(f"# r_99% chi_P = {r1_chi}, sieve-control = {r1_sieve:.0f}")
    print(f"# max |z| over k in [10, 100] vs sieve control = {z_max:.2f}")
    if z_max < 3.0:
        print("VERDICT: FAIL — after subtracting wheel-sieve structure, chi_P "
              "Hankel is statistically indistinguishable from random.")
    elif z_max < 5.0:
        print(f"VERDICT: WEAK — mild non-random Hankel structure (max z={z_max:.2f}).")
    else:
        print(f"VERDICT: SIGNAL — chi_P Hankel deviates significantly from "
              f"sieve-aware random (max z={z_max:.2f}).")


if __name__ == "__main__":
    main()
