"""
BBP-style digit-freeze test for psi(x) via the explicit formula.

Fresh framing
-------------
The Bailey-Borwein-Plouffe formula extracts the k-th hex digit of pi
(the constant) in O(k log k) without computing earlier digits. Its trick
is that the contribution of each summand to the k-th digit decays
geometrically.

Question: does psi(x) (or pi(x)) admit anything analogous?
The explicit formula
    psi(x) = x - sum_rho x^rho/rho - log(2*pi) - (1/2)*log(1 - 1/x^2)
is structurally similar to a Fourier series in log(x) with frequencies
gamma_k. Each zero contributes amplitude sqrt(x)/|rho_k|, decaying as 1/k.
If contributions had fast geometric decay (like BBP), polylog zeros would
fix any digit. If they have slow polynomial decay (like sqrt(x)/K), then
each new digit costs *exponentially* more zeros.

We test empirically:
  - For x in {1e4, 1e5, 1e6, 1e7}, with K zeros in {10, 50, 200, 1000, 2000},
    compute the partial-sum approximation psi_K(x) and the residual
    e_K(x) = psi_K(x) - psi_true(x).
  - Track top stable decimal digit of psi as K grows.
  - Look for anomalous x with super-fast convergence (would suggest
    cancellation structure usable for polylog extraction).

Outcome shapes the verdict:
  (a) e_K ~ sqrt(x) * f(log K)   --> standard sqrt(x) information bound,
                                    no BBP-style shortcut available.
  (b) e_K ~ sqrt(x) / K^a, a>1   --> better than naive, but still polynomial
                                    in zero count vs digits.
  (c) anomalous x with e_K << bound for moderate K --> structural shortcut
                                    candidate.
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


def true_pi_psi(N):
    primes = sieve_primes(N)
    pi_N = int(primes.size)
    psi_N = 0.0
    for p in primes.tolist():
        lp = math.log(p)
        pk = p
        while pk <= N:
            psi_N += lp
            pk *= p
    return pi_N, psi_N


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


def psi_partial(x, zeros_K):
    """
    psi(x) ~ x - 2*sqrt(x) * sum_k cos(gamma_k * log x - phi_k) / |rho_k|
            - log(2*pi)
    where rho_k = 1/2 + i*gamma_k,
          |rho_k| = sqrt(1/4 + gamma_k^2),
          phi_k = atan2(gamma_k, 1/2).
    """
    log_x = math.log(x)
    sqrt_x = math.sqrt(x)
    rho_abs = np.sqrt(0.25 + zeros_K * zeros_K)
    phi = np.arctan2(zeros_K, 0.5)
    osc = np.cos(zeros_K * log_x - phi) / rho_abs
    return x - 2.0 * sqrt_x * osc.sum() - math.log(2.0 * math.pi)


def top_digit_correct(true_val, approx):
    """
    Number of leading decimal digits agreeing between true_val and approx
    when both are written in scientific notation. Cap at 15.
    """
    if true_val == 0:
        return 0
    rel = abs(approx - true_val) / abs(true_val)
    if rel <= 0:
        return 15
    return min(15, max(0, -math.floor(math.log10(rel))))


def main():
    zeros = load_zeros(
        "/apps/aplikacijos/prime-research/data/zeta_zeros_2000.txt", 2000
    )
    print(f"Loaded {len(zeros)} zeros (gamma_max = {zeros[-1]:.2f})")
    print(
        f"Theoretical zero-density at T={zeros[-1]:.0f}: "
        f"~T*log(T)/(2*pi) = "
        f"{zeros[-1] * math.log(zeros[-1]) / (2 * math.pi):.0f}"
    )

    test_xs = [10**4, 10**5, 10**6, 10**7]
    K_values = [10, 50, 200, 1000, 2000]

    print("\n--- Main scan: psi(x) convergence vs K ---")
    print(
        f"{'x':>10} {'K':>6} {'psi_approx':>20} {'residual':>14} "
        f"{'|res|/sqrt(x)':>14} {'digits_OK':>10}"
    )
    for x in test_xs:
        _, psi_true = true_pi_psi(x)
        sqrt_x = math.sqrt(x)
        for K in K_values:
            zk = zeros[:K]
            psi_app = psi_partial(x, zk)
            res = psi_app - psi_true
            digits = top_digit_correct(psi_true, psi_app)
            print(
                f"{x:>10d} {K:>6d} {psi_app:>20.4f} {res:>+14.4f} "
                f"{abs(res) / sqrt_x:>14.4f} {digits:>10d}"
            )

    # Anomaly hunt: scan many x values at fixed x ~ 10^6 and a fixed K,
    # looking for x with anomalously small residual.
    print("\n--- Anomaly hunt: residual distribution at x near 10^6, K=2000 ---")
    K = 2000
    zk = zeros[:K]
    candidates = []
    base = 10**6
    delta_range = 4001
    # Compute psi_true once at x = base + delta_range, then sweep deltas.
    # For each delta, true psi(base+delta) requires re-sieving; instead, we use
    # the fact psi(N+1) - psi(N) = log p if N+1 = p^k for prime p, else 0.
    primes_full = sieve_primes(base + delta_range + 10)
    is_prime_set = set(primes_full.tolist())
    # Build prime-power set up to base+delta_range:
    prime_power_log = {}
    for p in primes_full.tolist():
        if p * p > base + delta_range:
            break
        pk = p
        while pk <= base + delta_range:
            prime_power_log[pk] = math.log(p)
            pk *= p
    # Add primes themselves (k=1):
    for p in primes_full.tolist():
        if p <= base + delta_range:
            prime_power_log[p] = math.log(p)
    # Compute psi_true at base
    _, psi_at_base = true_pi_psi(base)
    psi_true = psi_at_base
    n = base
    residuals = []
    for delta in range(0, delta_range):
        x = base + delta
        if delta > 0:
            # Increment psi by log p if x is a prime power
            if x in prime_power_log:
                psi_true += prime_power_log[x]
        psi_app = psi_partial(x, zk)
        res = psi_app - psi_true
        residuals.append((x, res))
    residuals_arr = np.array([r for _, r in residuals])
    sqrt_x = math.sqrt(base)
    print(f"  N samples: {len(residuals)}, x in [{base}, {base + delta_range - 1}]")
    print(f"  mean residual:        {residuals_arr.mean():+.4f}")
    print(f"  std residual:         {residuals_arr.std():.4f}")
    print(f"  min |residual|:       {np.min(np.abs(residuals_arr)):.4f}")
    print(f"  max |residual|:       {np.max(np.abs(residuals_arr)):.4f}")
    print(f"  sqrt(x):              {sqrt_x:.4f}")
    print(f"  best digit-agreement: {-math.log10(np.min(np.abs(residuals_arr))/base):.2f}")
    # Print top 5 anomalies (smallest |res|)
    abs_res = np.abs(residuals_arr)
    idx = np.argsort(abs_res)[:5]
    print("  top-5 smallest |residual| values:")
    for i in idx:
        x, r = residuals[i]
        print(f"    x={x}, residual={r:+.6f}, |res|/sqrt(x)={abs(r)/sqrt_x:.6f}")

    # Theoretical scaling check: how fast does max-|residual| shrink with K?
    print("\n--- Scaling check: max|residual| over x in [1e6, 1e6+200) vs K ---")
    print(f"{'K':>6} {'max|res|':>14} {'sqrt(x)*log(K)^2':>20} {'ratio':>10}")
    for K in [10, 50, 200, 1000, 2000]:
        zk = zeros[:K]
        psi_at_base_local = psi_at_base
        psi_true_local = psi_at_base_local
        max_res = 0.0
        for delta in range(0, 200):
            x = base + delta
            if delta > 0 and x in prime_power_log:
                psi_true_local += prime_power_log[x]
            psi_app = psi_partial(x, zk)
            r = abs(psi_app - psi_true_local)
            if r > max_res:
                max_res = r
        bound = sqrt_x * math.log(K) ** 2
        print(
            f"{K:>6d} {max_res:>14.4f} {bound:>20.4f} {max_res / bound:>10.4f}"
        )


if __name__ == "__main__":
    main()
