#!/usr/bin/env python3
"""S166 squarefree extension — direct Fourier verification.

Theorem (S168, extending S166):
    For squarefree q >= 2 and chi_P : {1,..,N} -> {0,1} the prime indicator,
        E(q, N) := (1/N) sum_{k coprime to q} |sum_{p' prime <= N} omega_q^{k p'}|^2
                = (mu(q))^2 * (pi(N) - r(q))^2 / (phi(q) * N) + R(q, N)
    where r(q) = number of distinct primes dividing q, and the remainder
    R(q, N) is bounded by O(q * Var(q, N) / N) where
        Var(q, N) = sum_{a in (Z/qZ)*} (pi(N; q, a) - (pi(N) - r)/phi(q))^2.

For non-squarefree q (mu(q) = 0): main term VANISHES, only the variance
remainder contributes.

This script verifies the theorem at d in {14, 16, 18, 20} for all
squarefree q in [2, 50] (and also tests a handful of non-squarefree q
to confirm the main-term-vanishes prediction).

Aggregate test: predicts that
    sum_{q squarefree, q <= Q*} E(q, N) ~ (pi(N)^2 / N) * sum_{q sqf, q<=Q*} 1/phi(q)
                                       ~ A * (pi(N)^2 log(Q*) / N)
matches the empirical sum across squarefree q. Then connects to S74's
N^{0.42} spike count via Q* ~ N^{0.21}.
"""

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


def factor(q):
    """Return list of (prime, exponent) for q."""
    f = []
    n = q
    p = 2
    while p * p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            f.append((p, e))
        p += 1
    if n > 1:
        f.append((n, 1))
    return f


def is_squarefree(q):
    if q <= 0:
        return False
    if q == 1:
        return True
    return all(e == 1 for _, e in factor(q))


def mobius(q):
    if q <= 0:
        return 0
    if q == 1:
        return 1
    f = factor(q)
    if any(e > 1 for _, e in f):
        return 0
    return (-1) ** len(f)


def euler_phi(q):
    if q == 1:
        return 1
    res = q
    for p, _ in factor(q):
        res = res // p * (p - 1)
    return res


def num_distinct_prime_factors(q):
    return len(factor(q))


def fourier_energy_at_q(prime_arr, q, N):
    """Compute E(q, N) = (1/N) sum_{k coprime to q} |S_q^k|^2.

    S_q^k = sum_{p' prime <= N} exp(2 pi i k p' / q).
    """
    energy = 0.0
    for k in range(1, q):
        if math.gcd(k, q) != 1:
            continue
        angle = 2 * np.pi * k * prime_arr / q
        c = float(np.sum(np.cos(angle)))
        s = float(np.sum(np.sin(angle)))
        energy += (c * c + s * s) / N
    return energy


def predicted_main(q, pi_N):
    """Main-term prediction: mu(q)^2 * (pi(N) - r(q))^2 / (phi(q) * N).

    Returns the numerator ((pi(N)-r)^2 / phi(q)) — divide by N to get energy.
    For non-squarefree q (mu = 0): returns 0.
    """
    mu = mobius(q)
    if mu == 0:
        return 0.0
    r = num_distinct_prime_factors(q)
    phi_q = euler_phi(q)
    return (pi_N - r) ** 2 / phi_q


def compute_variance(prime_arr, q, pi_N):
    """Compute Var(q, N) = sum_{a in (Z/qZ)*} (pi(N; q, a) - mu_avg)^2
    where mu_avg = (pi(N) - r) / phi(q)."""
    r = num_distinct_prime_factors(q)
    phi_q = euler_phi(q)
    mu_avg = (pi_N - r) / phi_q
    counts = np.zeros(q, dtype=np.int64)
    residues = prime_arr % q
    np.add.at(counts, residues, 1)
    var = 0.0
    for a in range(q):
        if math.gcd(a, q) == 1:
            var += (counts[a] - mu_avg) ** 2
    return var


def main():
    print("# S166 squarefree extension — empirical verification\n")

    # Per-q testing
    ds = [14, 16, 18, 20]

    Q_max = 50
    sqf_qs = [q for q in range(2, Q_max + 1) if is_squarefree(q)]
    nonsqf_qs = [q for q in range(2, Q_max + 1) if not is_squarefree(q)]

    for d in ds:
        N = 2 ** d
        prime_arr = primes_up_to(N)
        pi_N = len(prime_arr)
        print(f"\n=== d={d} N={N} pi(N)={pi_N} ===\n")

        # 1) per-q ratios for squarefree q
        print(f"  squarefree q in [2, {Q_max}]:  E_emp / E_pred (main term only)\n")
        print(f"  {'q':>4} {'phi':>4} {'r':>3} {'mu':>3} {'E_emp/N':>12} "
              f"{'E_pred_main':>12} {'ratio':>8} {'q*Var/N':>10}")
        ratios_sqf = []
        for q in sqf_qs:
            E_emp = fourier_energy_at_q(prime_arr, q, N)
            E_pred_num = predicted_main(q, pi_N)
            E_pred = E_pred_num / N
            var = compute_variance(prime_arr, q, pi_N)
            remainder_bound = q * var / N
            ratio = E_emp / E_pred if E_pred > 0 else float('nan')
            ratios_sqf.append((q, ratio))
            mu = mobius(q)
            r = num_distinct_prime_factors(q)
            phi_q = euler_phi(q)
            print(f"  {q:4d} {phi_q:4d} {r:3d} {mu:+3d} {E_emp:12.4f} "
                  f"{E_pred:12.4f} {ratio:8.5f} {remainder_bound:10.4f}")

        ratios_only = [r for _, r in ratios_sqf if not math.isnan(r)]
        print(f"\n  squarefree ratios: min={min(ratios_only):.4f}, "
              f"max={max(ratios_only):.4f}, median={np.median(ratios_only):.4f}, "
              f"mean={np.mean(ratios_only):.4f}")

        # 2) non-squarefree q: main term should be 0; only variance contributes
        if d == 18:
            print(f"\n  non-squarefree q in [2, {Q_max}] (main term vanishes):\n")
            print(f"  {'q':>4} {'phi':>4} {'E_emp':>10} {'q*Var/N':>10} "
                  f"{'ratio':>8}")
            for q in nonsqf_qs[:15]:
                E_emp = fourier_energy_at_q(prime_arr, q, N)
                var = compute_variance(prime_arr, q, pi_N)
                rem = q * var / N
                ratio = E_emp / rem if rem > 0 else float('nan')
                phi_q = euler_phi(q)
                print(f"  {q:4d} {phi_q:4d} {E_emp:10.4f} {rem:10.4f} {ratio:8.4f}")

        # 3) aggregate: sum over squarefree q <= Q* and compare with
        #    A * pi(N)^2 / N * log(Q*)
        print(f"\n  aggregate: sum E(q, N) for squarefree q <= Q*:\n")
        print(f"  {'Q*':>5} {'sum_emp':>14} {'sum_pred':>14} "
              f"{'pi(N)^2/N * log Q*':>18} {'A_emp':>8}")
        for Q_star in [10, 20, 50, 100, 200]:
            sqf_in_range = [q for q in range(2, Q_star + 1) if is_squarefree(q)]
            if Q_star == 200 and d <= 16:
                # Compute the larger sum using already-computed primes
                pass
            sum_emp = 0.0
            sum_pred = 0.0
            for q in sqf_in_range:
                if Q_star <= 100 or q <= 50:
                    if q <= Q_max:
                        # already in our table — but we computed via fresh call,
                        # which is fine
                        pass
                E_emp = fourier_energy_at_q(prime_arr, q, N)
                E_pred = predicted_main(q, pi_N) / N
                sum_emp += E_emp
                sum_pred += E_pred
            asymp = pi_N ** 2 / N * math.log(Q_star)
            A = sum_pred / asymp
            print(f"  {Q_star:5d} {sum_emp:14.4f} {sum_pred:14.4f} "
                  f"{asymp:18.4f} {A:8.4f}")

    # 4) Wirsing-type asymptotic: sum_{q sqf <= Q} 1/phi(q) ~ B * log Q
    print(f"\n=== Wirsing constant: sum_{{q sqf <= Q}} 1/phi(q) / log(Q) ===\n")
    print(f"  {'Q':>6} {'sum 1/phi(q)':>14} {'log Q':>10} {'ratio (-> B)':>14}")
    cum = 0.0
    for q in range(2, 5001):
        if is_squarefree(q):
            cum += 1.0 / euler_phi(q)
        if q in [10, 50, 100, 200, 500, 1000, 2000, 5000]:
            ratio = cum / math.log(q)
            print(f"  {q:6d} {cum:14.4f} {math.log(q):10.4f} {ratio:14.4f}")


if __name__ == "__main__":
    main()
