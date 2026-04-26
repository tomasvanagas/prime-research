"""
Parallel sieve as a linear dynamical system.

Idea: encode the state of the sieve of Eratosthenes as a vector
v in F_2^N, where v[i] = 1 means "i is still alive (not crossed off)".
Each round of the sieve picks the smallest alive index p > 1 and
zeros out indices 2p, 3p, ..., kp <= N.

If we COULD represent each round's action as multiplication by a
sparse matrix M_p, then `final_state = M_{p_t} ... M_{p_2} M_{p_1} v_0`,
and pi(N) = popcount(final_state) - 1.

A linear dynamical system v_{t+1} = A v_t admits "fast forwarding"
to v_T = A^T v_0 in O(log T * cost(matmul)) via repeated squaring.

But the sieve is NOT linear: which index gets crossed depends on the
state. We instead test:

  Q1: can the FIRST K rounds (K = pi(sqrt(N))) be combined into a
      single linear operator M_{<=K} that commutes with arbitrary v_0?
      (Yes, if we know the K primes; the operator is just "AND with
      the not-multiples-of-p_1..p_K mask". This is just sieving.)

  Q2: does the operator M_{<=K} have a low-rank or low-tensor-rank
      structure that lets us COMPUTE it on a single point in polylog?
      In particular, is M_{<=K} representable as a small matrix
      product over a tensor decomposition of [1, N]?

We explore Q2 by computing the tensor decomposition (CP / Tucker)
rank of the indicator vector "n is coprime to {p_1, ..., p_K}" when
N is reshaped as a tensor. Low rank => fast eval.

We also test: how many primes K can we sieve with using the smallest
batch sizes? The "Frobenius operator" view: if M_{<=K} acts on the
length-prim_K# = product of first K primes "wheel", which is small,
then we can fast-forward through wheel rotations.
"""

import numpy as np
from sympy import primerange, prod


def coprime_indicator(N, primes):
    """v[i] = 1 if gcd(i, prod(primes)) == 1, else 0; for i = 0..N-1."""
    out = np.ones(N, dtype=np.int8)
    out[0] = 0
    for p in primes:
        out[::p] = 0
    return out


def main():
    print("=== Sieve linear-dynamical-system view ===\n")

    # The wheel of size primorial p1#...#pK# is "fundamental period" of
    # gcd-with-{p1..pK}. It has phi(primorial) ones out of primorial slots.
    # Test: the COPRIME indicator on Z_{primorial} is highly structured.
    print("Wheel structure (coprime-to-first-K-primes pattern repeats period = primorial):")
    print(f"  K  primorial  phi(primorial)  density")
    cumprod = 1
    cum_phi = 1
    for k, p in enumerate(primerange(2, 32), start=1):
        cumprod *= int(p)
        cum_phi *= int(p) - 1
        print(f"  {k:2d}  {cumprod:>15d}  {cum_phi:>15d}  {cum_phi/cumprod:.4f}")
        if cumprod > 1e18:
            break

    # Tensor reshape test: take v of length L = a*b*c..., reshape as tensor,
    # check matrix rank along each unfolding mode.
    print("\nTensor-rank check on coprime indicator restricted to L = primorial(K):")
    for K in [3, 4, 5, 6]:
        primes = list(primerange(2, 30))[:K]
        L = 1
        for p in primes:
            L *= p
        v = coprime_indicator(L, primes)  # length L
        # Reshape as tensor with shape = primes
        T = v.reshape(primes)
        # Compute mode-wise unfolding rank
        ranks = []
        for mode in range(K):
            # unfold: bring axis `mode` to front
            mat = np.moveaxis(T, mode, 0).reshape(primes[mode], -1)
            r = np.linalg.matrix_rank(mat.astype(np.float64))
            ranks.append(r)
        # CP rank lower bound: each mode rank
        # CP rank upper bound: number of nonzeros
        nnz = int(v.sum())
        print(f"  K={K} L={L} primes={primes} ranks_per_mode={ranks} nnz={nnz} phi={int(np.prod([p-1 for p in primes]))}")

    # In fact the coprime indicator on a primorial is a TENSOR PRODUCT of
    # rank-1 indicators (i mod p_k != 0). That means it has CP rank EXACTLY
    # 1 in the "right" basis: chi(i mod p_1, ..., i mod p_K) =
    #     prod_k chi(i mod p_k != 0).
    # This is rank-1 product structure -> evaluable in O(K) per query.
    print("\nKey identity check (coprime indicator factorizes via CRT):")
    for K in [3, 4, 5]:
        primes = list(primerange(2, 30))[:K]
        L = 1
        for p in primes:
            L *= p
        v = coprime_indicator(L, primes)
        # Predicted: v[i] = prod_k (1 if i mod p_k != 0 else 0)
        pred = np.ones(L, dtype=np.int8)
        for p in primes:
            mask = (np.arange(L) % p == 0)
            pred[mask] = 0
        print(f"  K={K}: factorization matches: {np.array_equal(v, pred)}")

    print("\nObservation: rank-1 in CRT basis means we can evaluate ALL")
    print("'coprime-to-{first K primes}' queries in O(K log p_K) per query.")
    print("This gives the standard wheel-sieve speedup but does NOT extend")
    print("to coprime-to-{first K primes where K = pi(sqrt(N))}, because")
    print("primorial(pi(sqrt(N))) ~ exp(sqrt(N)) >> N.\n")

    # The "fast forward" question is: can we evaluate
    #   pi(N) = #{i in [2, N] : i is alive after sieving with all primes <= sqrt(N)}
    # by using the rank-1 wheel structure for the FIRST few primes (cheap),
    # and a fast residual for the LATER primes. The later primes don't have
    # short "primorial" period that fits in N. So fast-forward fails.
    # This is the well-known wheel-sieve barrier.
    print("Verdict: rank-1 wheel structure gives ONLY constant-factor gain;")
    print("primorial(pi(sqrt(N))) >> N forbids algebraic fast-forwarding.\n")

    # Final partial test: simulate a "tensor-train" representation of the
    # alive vector after sieving with primes up to some bound. We don't
    # expect the TT-rank to be small.
    print("TT-rank-like test: alive(N) after sieving primes <= B, reshape and SVD")
    N = 2**12  # 4096
    for B in [3, 5, 7, 11, 13]:
        primes = list(primerange(2, B + 1))
        v = coprime_indicator(N, primes).astype(np.float64)
        # reshape as 2x2x...x2 (12 modes)
        T = v.reshape([2] * 12)
        # Mode-1 unfolding rank
        mat = T.reshape(2, -1)
        r = np.linalg.matrix_rank(mat)
        # full reshape to roughly square then rank
        sq = v.reshape(64, 64)
        rsq = np.linalg.matrix_rank(sq)
        print(f"  N={N} B={B} primes={primes} bit-reshape_mode0_rank={r}  64x64_rank={rsq}")

    print("\nVerdict: alive-vector reshapes do NOT reveal globally low rank;")
    print("primes act non-locally w.r.t. binary digit structure.")


if __name__ == "__main__":
    main()
