"""
Spike-fraction signature across arithmetic indicators.

Tests whether the S168/S190 21%-spike fraction is specific to chi_P
(via the principal-character / coprime-support mechanism) or generic
to dense arithmetic indicators.

For each indicator f : [1, N] -> {0, 1}, computes:

    Spike(f, Q, N) = sum_{q sqf, 2 <= q <= Q} ||P_{V_q^prim} f||^2 / ||f||^2

where V_q^prim is the primitive-additive-character subspace at conductor q.

Indicators tested:
    f1 = chi_P              (S168/S190 control, expect ~21% at Q=N^{0.185})
    f2 = chi_{Omega odd}    (Liouville indicator (1 - lambda)/2)
    f3 = chi_{Omega = 2}    (semiprime indicator, similar density to chi_P at large N)
    f4 = chi_{Omega = 3}    (3-almost-prime)
    f5 = mu^2               (squarefree indicator)
    f6 = chi_{Omega even}   (complement of f2)
    f7 = chi_{Omega = 1}    (= chi_P, sanity)

Output: per-(f, Q, d) spike fraction + per-q breakdown.
"""

import json
import math
import sys
from pathlib import Path

import numpy as np


def sieve_omega_mu(N):
    """Return (chi_P, Omega, mu) arrays of length N+1."""
    is_p = np.ones(N + 1, dtype=bool)
    is_p[:2] = False
    Omega = np.zeros(N + 1, dtype=np.int32)
    mu = np.ones(N + 1, dtype=np.int8)
    mu[0] = 0
    # Smallest prime factor sieve
    spf = np.zeros(N + 1, dtype=np.int32)
    for i in range(2, N + 1):
        if spf[i] == 0:
            for j in range(i, N + 1, i):
                if spf[j] == 0:
                    spf[j] = i
    for i in range(4, N + 1):
        if spf[i] != i:  # composite
            is_p[i] = False
    # Omega and mu via prime-power factorisation
    for n in range(2, N + 1):
        x = n
        Om = 0
        sf = True
        while x > 1:
            p = spf[x]
            cnt = 0
            while x % p == 0:
                x //= p
                cnt += 1
            Om += cnt
            if cnt >= 2:
                sf = False
        Omega[n] = Om
        if not sf:
            mu[n] = 0
        elif Om % 2 == 0:
            mu[n] = 1
        else:
            mu[n] = -1
    return is_p, Omega, mu


def squarefree_qs(Qmax):
    is_sqf = [True] * (Qmax + 1)
    p = 2
    while p * p <= Qmax:
        for k in range(p * p, Qmax + 1, p * p):
            is_sqf[k] = False
        p += 1
    return [q for q in range(2, Qmax + 1) if is_sqf[q]]


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def per_q_energy(f, N, q):
    """
    Compute || P_{V_q^prim} f ||^2 for the primitive-character subspace.

    V_q^prim = span { exp(2 pi i a n / q) : a in [1, q-1], gcd(a, q) = 1 }.
    Inner product on [1, N] uses ||e_a^q||^2 = N when N is a multiple of q
    (when not, off by O(1) per a; we accept this finite-q error).

    Energy = (1/N) * sum_{a coprime to q, a != 0} | sum_r F_q[r] e^{-2 pi i a r / q} |^2
    where F_q[r] = sum_{n=1..N, n mod q == r} f(n).
    """
    # Build F_q (size q)
    n_idx = np.arange(1, N + 1)
    Fq = np.bincount(n_idx % q, weights=f[1 : N + 1].astype(np.float64), minlength=q)
    # Discrete Fourier transform on Z/qZ
    Fq_hat = np.fft.fft(Fq)
    # Sum |Fq_hat[a]|^2 for a coprime to q (a in [1, q-1])
    e = 0.0
    for a in range(1, q):
        if gcd(a, q) == 1:
            e += abs(Fq_hat[a]) ** 2
    return e / N


def spike_breakdown(f, N, Q):
    """
    Return dict { q: per_q_energy/||f||^2 } for q sqf in [2, Q],
    plus total spike fraction.
    """
    norm2 = float((f[1 : N + 1] ** 2).sum())
    if norm2 == 0.0:
        return {"per_q": {}, "spike_fraction": 0.0, "norm2": 0.0}
    sqfs = squarefree_qs(Q)
    per_q = {}
    total = 0.0
    for q in sqfs:
        e_q = per_q_energy(f, N, q)
        per_q[q] = e_q / norm2
        total += e_q
    return {"per_q": per_q, "spike_fraction": total / norm2, "norm2": norm2}


def run(d_list, Q_factors, out_dir):
    all_results = {}
    for d in d_list:
        N = 1 << d
        print(f"\n=== d={d}, N={N}", file=sys.stderr)
        is_p, Omega, mu = sieve_omega_mu(N)

        chi_P = is_p.astype(np.int8)
        chi_Om_odd = (Omega % 2 == 1).astype(np.int8)
        chi_Om_odd[0] = 0  # Omega(0) undefined; safe pad
        chi_Om2 = (Omega == 2).astype(np.int8)
        chi_Om3 = (Omega == 3).astype(np.int8)
        chi_mu2 = (mu != 0).astype(np.int8)
        chi_mu2[0] = 0  # exclude n=0
        chi_mu2[1] = 1  # mu(1) = 1, squarefree
        chi_Om_even = ((Omega % 2 == 0) & (np.arange(N + 1) >= 2)).astype(np.int8)

        indicators = {
            "chi_P": chi_P,
            "chi_Omega_odd": chi_Om_odd,
            "chi_Omega_2": chi_Om2,
            "chi_Omega_3": chi_Om3,
            "mu2": chi_mu2,
            "chi_Omega_even": chi_Om_even,
        }

        d_results = {"d": d, "N": N, "Q_targets": {}, "indicator_norm2": {}}
        for name, f in indicators.items():
            d_results["indicator_norm2"][name] = float((f[1 : N + 1] ** 2).sum())

        Q_list = sorted(set([2, 6, 8] + [max(2, round(N ** alpha)) for alpha in Q_factors] + [30]))
        d_results["Q_list"] = Q_list

        for Q in Q_list:
            d_results["Q_targets"][Q] = {}
            for name, f in indicators.items():
                br = spike_breakdown(f, N, Q)
                d_results["Q_targets"][Q][name] = {
                    "spike_fraction": br["spike_fraction"],
                    "norm2": br["norm2"],
                    # Keep per-q at lowest Q only to limit JSON size:
                    "per_q": br["per_q"] if Q == Q_list[0] or Q in (8, 13, 30) else None,
                }
            line = f"d={d:2d} Q={Q:3d} | " + " ".join(
                f"{name}={d_results['Q_targets'][Q][name]['spike_fraction']:.4f}"
                for name in indicators
            )
            print(line, file=sys.stderr)

        all_results[d] = d_results

    out_path = out_dir / "spike_indicator_specificity.json"
    with open(out_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nWrote {out_path}", file=sys.stderr)
    return all_results


def main():
    out_dir = Path(__file__).parent
    d_list = [14, 16, 18]
    Q_factors = [0.185, 0.21]
    run(d_list, Q_factors, out_dir)


if __name__ == "__main__":
    main()
