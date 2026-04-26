"""
Tensor-Train rank of the prime indicator chi_P on [0, 2^L).

Hypothesis: if chi_P has TT-rank polylog L, then pi(2^L) is computable
in polylog time via TT prefix-sum.

Compares against three baselines:
  - random Bernoulli with same density,
  - "n divisible by 7" (low TT-rank, ~7),
  - random binary string (max TT-rank).

Reports per-bipartition rank profile and singular-value tail.
"""
import numpy as np
from sympy import isprime


def tt_decompose(T_flat, dims, eps=1e-10):
    """TT-SVD with relative truncation eps. Returns ranks list."""
    L = len(dims)
    ranks = [1]
    A = T_flat.reshape(dims[0], -1)
    for ell in range(L - 1):
        U, S, Vt = np.linalg.svd(A, full_matrices=False)
        s0 = S[0] if S[0] > 0 else 1.0
        r = int(np.sum(S > eps * s0))
        r = max(r, 1)
        ranks.append(r)
        A = (np.diag(S[:r]) @ Vt[:r])
        # reshape into (r * dims[ell+1], rest)
        if ell + 1 < L - 1:
            A = A.reshape(r * dims[ell + 1], -1)
        else:
            A = A.reshape(r * dims[ell + 1], 1)
    return ranks


def tt_rank_profile(arr, name, L):
    arr = np.asarray(arr, dtype=np.float64)
    assert arr.size == 2 ** L
    ranks = tt_decompose(arr, [2] * L)
    print(f"{name:30s} TT-ranks: {ranks}")
    print(f"{name:30s} max rank = {max(ranks)}, sum = {sum(ranks)}")
    return ranks


def main():
    rng = np.random.default_rng(42)
    for L in [10, 12, 14, 16]:
        N = 2 ** L
        print(f"\n=== L = {L}, N = 2^{L} = {N} ===")

        # Prime indicator
        chi = np.array([1.0 if isprime(i) else 0.0 for i in range(N)])
        density = chi.mean()

        # Baselines
        random_bern = (rng.random(N) < density).astype(np.float64)
        div_by_7 = np.array([1.0 if i % 7 == 0 else 0.0 for i in range(N)])
        random_binary = (rng.random(N) < 0.5).astype(np.float64)

        ranks_prime = tt_rank_profile(chi, "prime indicator chi_P", L)
        ranks_bern = tt_rank_profile(random_bern, f"random Bernoulli p={density:.4f}", L)
        ranks_div7 = tt_rank_profile(div_by_7, "n divisible by 7", L)
        ranks_rand = tt_rank_profile(random_binary, "random binary p=0.5", L)

        # Compression ratio: TT params vs raw
        def tt_params(ranks):
            total = 0
            ranks_full = ranks + [1]
            for ell in range(L):
                total += ranks_full[ell] * 2 * ranks_full[ell + 1]
            return total

        print(f"\nL={L} compression (TT params / 2^L):")
        print(f"  prime:    {tt_params(ranks_prime):>8d} / {N} = {tt_params(ranks_prime)/N:.4f}")
        print(f"  Bernoulli:{tt_params(ranks_bern):>8d} / {N} = {tt_params(ranks_bern)/N:.4f}")
        print(f"  div7:     {tt_params(ranks_div7):>8d} / {N} = {tt_params(ranks_div7)/N:.4f}")
        print(f"  random:   {tt_params(ranks_rand):>8d} / {N} = {tt_params(ranks_rand)/N:.4f}")


if __name__ == "__main__":
    main()
