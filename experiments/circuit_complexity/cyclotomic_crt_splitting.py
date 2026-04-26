"""
Cyclotomic CRT splitting for AKS-style polynomial powering.

Question: Does the CRT decomposition Z_n[x]/(x^r - 1) ~= prod_{d|r} Z_n[x]/Phi_d(x)
reduce the maximum matrix-powering dimension below r, enough to threaten the
'AKS in TC^0 via matrix powering' OPEN status?

Setup:
  - For each n in a sample, find the AKS-prescribed r:
      smallest r > 1 with gcd(n, r) = 1 and ord_r(n) > floor(log2(n))^2.
  - Compute the CRT factorization x^r - 1 = prod_{d|r} Phi_d(x).
  - The MAX dimension over factors is max_{d|r} deg(Phi_d) = max_{d|r} phi(d).
  - For TC^0-feasibility, this max-dim must be o(polylog n) (constant depth
    only handles fixed dimension; polylog dim is the genuine open frontier).

Output: For n ranging over primes and composites in [10..10^6], record
  r, log2(n)^2, max_dim/r ratio, dim of largest cyclotomic factor.

Expected obstacle: AKS picks the SMALLEST r passing the ord test. Empirically
this r is usually prime, so phi(r) = r-1 -- the cyclotomic CRT decomposition
saves almost nothing. The 'split' factors collectively still require the SAME
max dimension as the original.
"""

from math import gcd, log2, floor
from sympy import totient, divisors, isprime
import bisect


def aks_r(n: int) -> int:
    """Smallest r > 1 coprime to n with ord_r(n) > floor(log2(n))^2."""
    bound = floor(log2(n)) ** 2
    r = 2
    while True:
        if gcd(n, r) == 1:
            # Compute multiplicative order of n mod r
            k = 1
            v = n % r
            while v != 1:
                v = (v * n) % r
                k += 1
                if k > bound:
                    return r
        r += 1


def aks_max_dim(n: int) -> dict:
    r = aks_r(n)
    divs = divisors(r)
    phis = [int(totient(d)) for d in divs]
    max_dim = max(phis)
    return {
        "n": n,
        "r": r,
        "log2_n_sq": floor(log2(n)) ** 2,
        "max_dim": max_dim,
        "ratio_max_dim_over_r": max_dim / r,
        "num_factors": len(divs),
        "r_is_prime": isprime(r),
    }


def main():
    print(f"{'n':>10} | {'r':>6} | {'log2(n)^2':>10} | {'max_dim':>8} | "
          f"{'max/r':>6} | {'#factors':>8} | {'r prime':>8}")
    print("-" * 72)

    # A representative sample
    sample_ns = []
    # primes near powers of 10
    for k in range(2, 7):
        base = 10 ** k
        # nearest prime above
        m = base + 1
        while not isprime(m):
            m += 1
        sample_ns.append(m)
        # nearest composite (semiprime) above
        m = base + 1
        while isprime(m):
            m += 1
        sample_ns.append(m)
    # a few specific Carmichael-like / Fermat pseudoprimes
    sample_ns.extend([561, 1105, 1729, 2465, 2821, 6601, 8911,
                      1024, 2048, 4096, 65537, 524287])
    sample_ns = sorted(set(sample_ns))

    rows = []
    for n in sample_ns:
        if n < 4:
            continue
        d = aks_max_dim(n)
        rows.append(d)
        print(f"{d['n']:>10} | {d['r']:>6} | {d['log2_n_sq']:>10} | "
              f"{d['max_dim']:>8} | {d['ratio_max_dim_over_r']:>6.3f} | "
              f"{d['num_factors']:>8} | {str(d['r_is_prime']):>8}")

    # Summary statistics
    print("\n--- SUMMARY ---")
    r_prime_count = sum(1 for r in rows if r["r_is_prime"])
    print(f"Sample size: {len(rows)}")
    print(f"r is prime in {r_prime_count} of {len(rows)} cases "
          f"({100*r_prime_count/len(rows):.1f}%)")
    avg_ratio = sum(r["ratio_max_dim_over_r"] for r in rows) / len(rows)
    min_ratio = min(r["ratio_max_dim_over_r"] for r in rows)
    print(f"Average max_dim / r = {avg_ratio:.4f}")
    print(f"Minimum max_dim / r = {min_ratio:.4f}")
    print(f"This means: cyclotomic CRT reduces the largest matrix-powering")
    print(f"dimension by AT MOST a factor of {1/min_ratio:.2f} in this sample.")
    print(f"For TC^0-feasibility we would need max_dim = O(1), i.e. ratio -> 0.")


if __name__ == "__main__":
    main()
