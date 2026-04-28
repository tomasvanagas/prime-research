#!/usr/bin/env python3
"""Direct Fourier coefficient computation at d=20 (N=2^20).

Without doing the full SVD, compute |S_q^k|^2 = |sum_{p' prime <= N}
exp(2 pi i k p' / q)|^2 for k coprime to q. The sum over k gives
N * E_V_q^primitive (the energy of chi_P in V_q^primitive subspace).

Compares with empirical block sums from S82's d=20 JSON.
"""

import json
import math
import numpy as np


def sieve(N):
    is_p = np.ones(N + 1, dtype=bool)
    is_p[:2] = False
    for i in range(2, int(N**0.5) + 1):
        if is_p[i]:
            is_p[i * i :: i] = False
    return is_p


def primes_up_to(N):
    return np.nonzero(sieve(N))[0]


def fourier_coef_at_q(prime_list, q, k, N):
    """Compute S_q^k = sum_{p' prime <= N} exp(2 pi i k * p' / q)."""
    angle = 2 * np.pi * k * prime_list / q
    s = float(np.sum(np.cos(angle))) + 1j * float(np.sum(np.sin(angle)))
    return s


def primitive_v_q_energy_direct(prime_list, N, q):
    """Sum over k coprime to q of |S_q^k|^2 / N. This equals the energy
    of chi_P (as length-N vector with values in {0,1}) in V_q^primitive
    SUBSPACE, NOT yet 2D-reshape-restricted.

    Note: this is the FULL N-vector energy of chi_P in V_q^primitive,
    which is a subspace of R^N. It does NOT depend on the (W^j, W^{d-j})
    reshape of chi_P. The 2D reshape's (W^j, W^{d-j})-rank-2 basis vectors
    span the same V_q^primitive subspace because each character lift is
    a function of (n mod q), which is preserved by reshape.
    """
    energy = 0.0
    by_k = {}
    for k in range(1, q):
        if math.gcd(k, q) != 1:
            continue
        S = fourier_coef_at_q(prime_list, q, k, N)
        contrib = (abs(S) ** 2) / N
        by_k[k] = contrib
        energy += contrib
    return energy, by_k


def main():
    # Compute for d = 14, 16, 18, 20.
    for d in [14, 16, 18, 20]:
        N = 2 ** d
        primes = primes_up_to(N)
        pi_N = len(primes)
        print(f"\n=== d={d} N={N} pi(N)={pi_N} ===")

        for p in [3, 5, 7, 11]:
            e_p, by_k_p = primitive_v_q_energy_direct(primes, N, p)
            e_2p, by_k_2p = primitive_v_q_energy_direct(primes, N, 2 * p)
            e_4p, by_k_4p = primitive_v_q_energy_direct(primes, N, 4 * p) if 4 * p <= 256 else (0.0, {})
            e_8p, _ = primitive_v_q_energy_direct(primes, N, 8 * p) if 8 * p <= 256 else (0.0, {})
            e_16p, _ = primitive_v_q_energy_direct(primes, N, 16 * p) if 16 * p <= 256 else (0.0, {})

            # Predicted main term: pi(N)^2 / ((p-1) N).
            pred_p = pi_N ** 2 / ((p - 1) * N)
            pred_2p = pi_N ** 2 / ((p - 1) * N)
            stack_pred = pred_p + pred_2p

            stack_emp = e_p + e_2p + e_4p + e_8p + e_16p

            # Compute K = stack * p * log^2(N) / N.
            log2N = math.log(N) ** 2
            K = stack_emp * p * log2N / N
            K_pred = stack_pred * p * log2N / N

            print(f"  p={p}: V_p={e_p:.2f}  V_2p={e_2p:.2f}  V_4p={e_4p:.2f}  "
                  f"V_8p={e_8p:.2f}  V_16p={e_16p:.2f}")
            print(f"      pred V_p + V_2p = {stack_pred:.2f}, empirical V_p + V_2p = {e_p + e_2p:.2f}, "
                  f"ratio = {(e_p + e_2p)/stack_pred:.3f}")
            print(f"      stack pred (V_p + V_2p only): {stack_pred:.2f}, "
                  f"K_pred = {K_pred:.3f}; stack actual (V_p..V_16p): {stack_emp:.2f}, K = {K:.3f}")


if __name__ == "__main__":
    main()
