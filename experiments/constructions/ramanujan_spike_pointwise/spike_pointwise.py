"""
Pointwise Ramanujan-spike approximator T_Q(n) for chi_P.

Construction (project-internal derivation from S166 / S168 / S169):

    P_{V_q^prim} chi_P (n)  ~=  (pi(N) - r(q)) * mu(q) c_q(n) / (N * phi(q))   (squarefree q)

Summing over squarefree q in [1, Q] gives the predicted "spike-block"
representation of chi_P pointwise:

    T_Q(n) := (pi(N)/N) * M_Q(n)
    M_Q(n) := sum_{sqf q <= Q} mu(q) c_q(n) / phi(q)

We use the Hoelder identity for Ramanujan sums to obtain a closed-form
simplification (proved purely from definitions of mu, phi, c_q):

    For q squarefree, gcd(q, n) = d:
        c_q(n) = mu(q/d) * phi(d)
        => mu(q) c_q(n) / phi(q) = mu(d) / phi(q/d)
        => mu(q) c_q(n) / phi(q) = mu(gcd(q,n)) / phi(q/gcd(q,n))

So
    M_Q(n) = sum_{sqf q <= Q} mu(gcd(q,n)) / phi(q/gcd(q,n))

This is the new pointwise scalar field.  By S168/S169:

  ||T_Q||^2 / ||chi_P||^2  -->  Frac(Q, N) * (pi(N) log N / N)
  with Frac(Q, N) = log Q / log N for squarefree-Mertens cumulative,
  predicted ~ 0.21 at Q = N^{0.185} (S169 measurement).

This script:
  1. Builds chi_P, mu, phi at N = 2^d for d in {14, 16, 18}.
  2. Computes T_Q(n) for n in [1, N] across multiple Q.
  3. Measures L2 norm ratio, correlation with chi_P, predictive power.
  4. Reports falsification check: does the L2 ratio match S169 at Q = N^{0.185}?
"""

import math
import json
import os
import sys
from pathlib import Path


def sieve(N):
    """Boolean primality vector chi_P[1..N] (chi_P[0] = 0 padding)."""
    is_p = [False, False] + [True] * (N - 1)
    for p in range(2, int(math.isqrt(N)) + 1):
        if is_p[p]:
            for k in range(p * p, N + 1, p):
                is_p[k] = False
    return is_p  # length N+1


def mobius_table(M):
    """mu[1..M]."""
    mu = [0] * (M + 1)
    mu[1] = 1
    primes = []
    is_p = [True] * (M + 1)
    is_p[0] = is_p[1] = False
    smallest = [0] * (M + 1)
    for i in range(2, M + 1):
        if is_p[i]:
            primes.append(i)
            smallest[i] = i
        for p in primes:
            if p * i > M:
                break
            is_p[p * i] = False
            smallest[p * i] = p
            if i % p == 0:
                break
    # Compute mu via square-free check using smallest prime factor
    for n in range(2, M + 1):
        x = n
        cnt = 0
        sf = True
        while x > 1:
            p = smallest[x]
            c = 0
            while x % p == 0:
                x //= p
                c += 1
            if c >= 2:
                sf = False
                break
            cnt += 1
        mu[n] = 0 if not sf else (-1 if cnt % 2 else 1)
    return mu


def totient_table(M):
    """phi[1..M] via sieve."""
    phi = list(range(M + 1))
    for i in range(2, M + 1):
        if phi[i] == i:  # i is prime
            for j in range(i, M + 1, i):
                phi[j] -= phi[j] // i
    return phi


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def squarefree_list(M, mu):
    return [q for q in range(1, M + 1) if mu[q] != 0]


def compute_M_Q(N, Q, mu, phi, sqf_qs_leQ):
    """
    M_Q(n) = sum_{sqf q <= Q} mu(gcd(q,n)) / phi(q/gcd(q,n))   for n in [1, N].

    Vectorised by q (outer loop): for each squarefree q, walk over n in [1, N]
    in steps that group by d = gcd(q, n).  For squarefree q, d ranges over
    divisors of q.  For each divisor d of q, the residues n with gcd(q,n) = d
    are exactly those n with d | n and gcd(n/d, q/d) = 1.

    Implementation: simpler O(N * #sqf_q) loop with explicit gcd, since
    Q is at most ~N^{0.21} ~ 22 at N=2^18.  Total ops ~ 5e7, fast.
    """
    M = [0.0] * (N + 1)
    for q in sqf_qs_leQ:
        muq = mu[q]
        # mu is +/-1 for squarefree q, never 0.
        for n in range(1, N + 1):
            d = gcd(q, n)
            md = mu[d]  # d divides q, q squarefree => d squarefree => |mu[d]|=1
            phi_q_over_d = phi[q // d]
            # mu(q) c_q(n)/phi(q) == mu(d)/phi(q/d)  by Hoelder identity
            M[n] += md / phi_q_over_d
    return M


def measure(N, d_log2, Q_list):
    is_p = sieve(N)
    M_max = max(Q_list) + 1
    mu = mobius_table(M_max)
    phi = totient_table(M_max)

    pi_N = sum(1 for n in range(1, N + 1) if is_p[n])
    chi_norm2 = pi_N  # ||chi_P||^2 = pi(N)
    mean_chi = pi_N / N

    results = {}
    for Q in Q_list:
        sqf_qs_leQ = [q for q in range(1, Q + 1) if mu[q] != 0]

        M_Q = compute_M_Q(N, Q, mu, phi, sqf_qs_leQ)

        # T_Q(n) := (pi(N)/N) * M_Q(n)
        T_norm2 = 0.0
        T_norm2_no_const = 0.0  # excluding q=1 wheel mode
        inner = 0.0
        # Subtract wheel mode T_1(n) = pi(N)/N * 1, then re-norm.
        # M_Q includes q=1 contribution = +1 (constant).
        # So T_Q without q=1 = T_Q - mean_chi.
        sum_T_p = sum_T_c = 0.0
        sumsq_T_p = sumsq_T_c = 0.0
        n_p = pi_N
        n_c = N - pi_N
        T_vals = [0.0] * (N + 1)
        for n in range(1, N + 1):
            t = mean_chi * M_Q[n]
            T_vals[n] = t
            T_norm2 += t * t
            t0 = t - mean_chi  # subtract wheel-mode constant
            T_norm2_no_const += t0 * t0
            if is_p[n]:
                inner += t
                sum_T_p += t
                sumsq_T_p += t * t
            else:
                sum_T_c += t
                sumsq_T_c += t * t

        # ||chi_P - T_Q||^2 = ||chi_P||^2 - 2<chi_P, T_Q> + ||T_Q||^2
        residual2 = chi_norm2 - 2 * inner + T_norm2

        mean_T_p = sum_T_p / n_p if n_p > 0 else 0.0
        var_T_p = sumsq_T_p / n_p - mean_T_p * mean_T_p if n_p > 0 else 0.0
        mean_T_c = sum_T_c / n_c if n_c > 0 else 0.0
        var_T_c = sumsq_T_c / n_c - mean_T_c * mean_T_c if n_c > 0 else 0.0

        # Pearson correlation (mean-subtracted)
        # cov(chi_P, T_Q) = E[chi_P*T_Q] - E[chi_P]E[T_Q]
        E_chi = pi_N / N
        E_T = (sum_T_p + sum_T_c) / N
        E_chiT = inner / N
        cov_chiT = E_chiT - E_chi * E_T
        var_chi = E_chi * (1 - E_chi)
        var_T = (sumsq_T_p + sumsq_T_c) / N - E_T * E_T
        pearson = cov_chiT / math.sqrt(var_chi * var_T) if var_chi * var_T > 0 else 0.0

        # Predictive accuracy: rank n by T_Q(n).  Top pi(N) ranked entries -
        # how many are actually prime?  This is "precision-at-pi(N)".
        # Sort all n in [1, N] by -T_Q(n).
        order = sorted(range(1, N + 1), key=lambda n: -T_vals[n])
        top_pi = order[:pi_N]
        true_positives = sum(1 for n in top_pi if is_p[n])
        precision_at_piN = true_positives / pi_N

        results[Q] = {
            "Q": Q,
            "pi_N": pi_N,
            "chi_norm2": chi_norm2,
            "T_norm2": T_norm2,
            "T_norm2_no_const": T_norm2_no_const,
            "T_over_chi": T_norm2 / chi_norm2,
            "T_over_chi_no_const": T_norm2_no_const / chi_norm2,
            "inner_chi_T": inner,
            "residual_norm2": residual2,
            "residual_over_chi": residual2 / chi_norm2,
            "raw_correlation": inner / math.sqrt(chi_norm2 * T_norm2)
            if T_norm2 > 0
            else 0.0,
            "pearson_correlation": pearson,
            "mean_T_prime": mean_T_p,
            "mean_T_composite": mean_T_c,
            "predictive_gap": mean_T_p - mean_T_c,
            "var_T_prime": var_T_p,
            "var_T_composite": var_T_c,
            "precision_at_piN": precision_at_piN,
            "baseline_precision": pi_N / N,  # random guess
            "lift_over_baseline": precision_at_piN / (pi_N / N)
            if pi_N > 0
            else 0.0,
            "n_squarefree_q_leQ": len(sqf_qs_leQ),
            "log_Q_over_log_N": math.log(Q) / math.log(N) if Q > 1 else 0.0,
        }
    return {"d": d_log2, "N": N, "pi_N": pi_N, "by_Q": results}


def main():
    out_dir = Path(__file__).parent
    all_results = {}

    Ns = [(14, 1 << 14), (16, 1 << 16), (18, 1 << 18), (20, 1 << 20)]
    for d, N in Ns:
        # S169 prediction: spike fraction at Q ~ N^{0.185}.
        Q_target = max(2, round(N ** 0.185))
        Q_list = sorted(
            set(
                [
                    2,
                    6,
                    Q_target,
                    max(2, round(N ** 0.21)),
                    30,
                    max(2, round(N ** 0.5)),
                ]
            )
        )
        print(
            f"\n=== d={d}, N={N}, pi(N)={'?'} -- Q_targets={Q_list}",
            file=sys.stderr,
        )
        res = measure(N, d, Q_list)
        all_results[d] = res
        for Q, r in res["by_Q"].items():
            print(
                f"d={d:2d} N={N:8d} Q={Q:5d}  "
                f"||T||/pi(N)={r['T_over_chi']:.4f} (no_const={r['T_over_chi_no_const']:.4f})  "
                f"raw={r['raw_correlation']:.3f}  pear={r['pearson_correlation']:.3f}  "
                f"prec@pi={r['precision_at_piN']:.3f}  lift={r['lift_over_baseline']:.2f}",
                file=sys.stderr,
            )

    with open(out_dir / "spike_pointwise_results.json", "w") as f:
        json.dump(all_results, f, indent=2)

    print(json.dumps(all_results, indent=2))


if __name__ == "__main__":
    main()
