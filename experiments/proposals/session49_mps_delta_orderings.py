#!/usr/bin/env python3
"""
Session 49 / Proposal A — Tensor-train bond dimension of delta(n) under
multiple index orderings.

The hypothesis: although prior tensor-network sieve experiments found
volume-law entanglement under the *natural* ordering, a non-trivial
permutation of the index n may yield low-rank TT structure. If even one
of the orderings tested below produces bond dim D = polylog(N), then
delta(n) is compactly representable and computable in polylog time.

We test these orderings:
  1. identity (1,2,...,N)
  2. bit-reversal of the binary index
  3. Gray-code ordering
  4. Morton (Z-order) of (high half, low half) of n
  5. 2-adic ordering (sort by largest power of 2 dividing n, then by quotient)
  6. sort by R^{-1}(n) (already a smooth coordinate)
  7. random permutation (sanity baseline)
  8. true Gaussian noise (volume-law baseline)

For each ordering, we reshape delta into a 2x2x...x2 tensor (log2(N) legs),
do sequential SVD with epsilon = 1e-3 truncation, and record bond dimensions
at every cut.

Output: experiments/proposals/session49_mps_delta_orderings_results.md
"""

import math
import time
import numpy as np

CACHE = "/apps/aplikacijos/prime-research/experiments/proposals/session49_data.npz"

# -------- Index orderings ----------------------------------------------------

def identity_perm(N):
    return np.arange(N)

def bit_reverse_perm(N):
    L = int(round(math.log2(N)))
    assert 2**L == N
    out = np.zeros(N, dtype=np.int64)
    for i in range(N):
        b = i
        r = 0
        for _ in range(L):
            r = (r << 1) | (b & 1)
            b >>= 1
        out[i] = r
    return out

def gray_perm(N):
    """Map index i to its Gray code reflection: i -> i XOR (i>>1)."""
    return np.array([i ^ (i >> 1) for i in range(N)], dtype=np.int64)

def morton_perm(N):
    """Z-order: interleave bits of (a, b) where i = a + 2^half_lo * b.
    For L = half_lo + half_hi (with half_hi >= half_lo), the result fits in L bits.
    """
    L = int(round(math.log2(N)))
    half_lo = L // 2
    half_hi = L - half_lo
    out = np.zeros(N, dtype=np.int64)
    for i in range(N):
        a = i & ((1 << half_lo) - 1)
        b = i >> half_lo
        m = 0
        # Interleave bits 0..half_lo-1 of a with bits 0..half_lo-1 of b
        for k in range(half_lo):
            m |= ((a >> k) & 1) << (2 * k)
            m |= ((b >> k) & 1) << (2 * k + 1)
        # Append remaining bits of b at top
        if half_hi > half_lo:
            m |= (b >> half_lo) << (2 * half_lo)
        out[i] = m
    # Verify it is a permutation in [0, N)
    assert out.max() < N and out.min() >= 0 and len(set(out.tolist())) == N
    return out

def two_adic_perm(N):
    """Sort by 2-adic valuation, then by quotient. Index i+1 (i=0..N-1)."""
    keys = []
    for i in range(N):
        n = i + 1
        v = 0
        m = n
        while m % 2 == 0:
            v += 1
            m //= 2
        keys.append((v, m, i))
    keys.sort()
    return np.array([k[2] for k in keys], dtype=np.int64)

def sort_by_smooth_perm(N, smooth_values):
    return np.argsort(smooth_values)

def random_perm(N, seed=0):
    rng = np.random.default_rng(seed)
    p = np.arange(N)
    rng.shuffle(p)
    return p

# -------- Tensor train decomposition ----------------------------------------

def tt_bond_dimensions(seq, eps_rel=1e-3):
    """Sequential SVD on a 2-leg-per-site reshape. Returns list of bond dims."""
    N = len(seq)
    L = int(round(math.log2(N)))
    assert 2**L == N
    T = np.array(seq, dtype=np.float64).reshape([2]*L)
    bonds = []
    left_dim = 1
    for site in range(L - 1):
        # Reshape to (left_dim*2, rest)
        mat = T.reshape(left_dim * 2, -1)
        U, S, Vt = np.linalg.svd(mat, full_matrices=False)
        # Truncate at eps_rel * S[0]
        if S[0] < 1e-15:
            r = 1
        else:
            r = int(np.sum(S > eps_rel * S[0]))
            r = max(r, 1)
        bonds.append(r)
        # Continue with the truncated remainder; for analysis we don't actually
        # need to keep U, just shape. Keep S * Vt as the new tensor.
        new_T = np.diag(S[:r]) @ Vt[:r, :]
        T = new_T.reshape([r] + [2]*(L - 1 - site))
        left_dim = r
    return bonds

# -------- Main --------------------------------------------------------------

def main():
    N = 8192  # 2^13
    print(f"Loading cached delta(n) for n=1..{N}", flush=True)
    d = np.load(CACHE)
    primes = d["primes"][:N]
    rinv = d["rinv"][:N]
    delta = d["delta"][:N]
    print(f"delta range: [{delta.min():.2f}, {delta.max():.2f}]  std={delta.std():.3f}", flush=True)

    # Build orderings
    print("Building orderings...", flush=True)
    orderings = {
        "identity": identity_perm(N),
        "bit_reverse": bit_reverse_perm(N),
        "gray": gray_perm(N),
        "morton": morton_perm(N),
        "two_adic": two_adic_perm(N),
        "sort_by_Rinv": sort_by_smooth_perm(N, rinv),
        "random": random_perm(N, seed=42),
    }

    # Baselines
    rng = np.random.default_rng(7)
    gauss_noise = rng.normal(size=N) * float(delta.std())
    smooth_signal = np.cos(2 * np.pi * np.arange(N) / 50)  # rank-1-ish

    # Pre-normalize each sequence (for cleaner SVD comparison)
    def normalize(x):
        x = np.asarray(x, dtype=np.float64)
        n = np.linalg.norm(x)
        return x / n if n > 0 else x

    rows = []
    print("\nTT bond-dim analysis (eps_rel=1e-3):", flush=True)
    print(f"{'ordering':20s}  {'max_bond':>8s}  {'sum_bond':>9s}  {'bonds':<60s}")
    for name, perm in orderings.items():
        seq = normalize(delta[perm])
        bonds = tt_bond_dimensions(seq, eps_rel=1e-3)
        rows.append((name, max(bonds), sum(bonds), bonds))
        print(f"{name:20s}  {max(bonds):8d}  {sum(bonds):9d}  {bonds}")

    seq = normalize(gauss_noise)
    bonds = tt_bond_dimensions(seq, eps_rel=1e-3)
    rows.append(("gaussian_noise", max(bonds), sum(bonds), bonds))
    print(f"{'gaussian_noise':20s}  {max(bonds):8d}  {sum(bonds):9d}  {bonds}")

    seq = normalize(smooth_signal)
    bonds = tt_bond_dimensions(seq, eps_rel=1e-3)
    rows.append(("smooth_cosine", max(bonds), sum(bonds), bonds))
    print(f"{'smooth_cosine':20s}  {max(bonds):8d}  {sum(bonds):9d}  {bonds}")

    # Save results
    out_md = "/apps/aplikacijos/prime-research/experiments/proposals/session49_mps_delta_orderings_results.md"
    L = int(round(math.log2(N)))
    with open(out_md, "w") as f:
        f.write("# session49_mps_delta_orderings — results\n\n")
        f.write(f"N={N} (L={L} legs), eps_rel=1e-3 SVD truncation.\n\n")
        f.write(f"delta(n) = p(n) - R^(-1)(n) for n=1..{N}, std={delta.std():.3f}.\n\n")
        f.write("| ordering | max_bond | sum_bond | bond profile |\n")
        f.write("|---|---|---|---|\n")
        for name, mb, sb, bonds in rows:
            f.write(f"| {name} | {mb} | {sb} | {bonds} |\n")
        f.write("\n## Interpretation\n\n")
        f.write("Volume-law (max_bond ~ N^{1/2} = 90) means TT compression buys\n")
        f.write("nothing for delta. Polylog bond dim under any ordering would be\n")
        f.write("a major positive signal. Smooth cosine baseline shows the rank-2\n")
        f.write("ceiling for a single-frequency signal.\n\n")
        max_natural = max(b for n_, b, _, _ in rows if n_ == "identity")
        max_gauss = max(b for n_, b, _, _ in rows if n_ == "gaussian_noise")
        f.write(f"natural ordering max_bond = {max_natural}, gaussian = {max_gauss}.\n")
        if max_natural < max_gauss * 0.7:
            f.write("VERDICT: delta has notable structure under natural ordering.\n")
        else:
            f.write("VERDICT: natural ordering looks volume-law-like.\n")
        best_ord = min(rows, key=lambda r: r[1])
        if best_ord[0] not in ("smooth_cosine",):
            f.write(f"Best ordering: **{best_ord[0]}** with max_bond={best_ord[1]}.\n")
            if best_ord[1] <= 3 * L:
                f.write("This bond dim is *sub-linear* in N^{1/2} -- worth deeper study.\n")
            else:
                f.write("Bond dim still scales like sqrt(N); proposal A closed.\n")

    print(f"\nWrote {out_md}", flush=True)

if __name__ == "__main__":
    main()
