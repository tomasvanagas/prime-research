"""
Proposal C: Test whether Hecke eigenvalues / arithmetic-function fingerprints
can discriminate primes from composites in [1, 10000].

For each n in 1..10000 compute a fingerprint:
  f(n) = (tau(n) mod 691,
          tau(n) mod 5,
          sigma_1(n) mod n,             # this is 1 for primes
          sigma_3(n) mod 691,           # related to E_4 Eisenstein series
          jacobi-like a_p of E: y^2 = x^3 - x mod various small primes,
          ...)

For primes, by Ramanujan, tau(p) ≡ p^11 + 1 (mod 691). Composites generally
do not satisfy this congruence pattern.

We measure: how big a fingerprint do we need before all composites in
[1,10000] differ from every prime's fingerprint?

If the answer is "polylog(N)", then we have an O(polylog) primality oracle.
If we can additionally count primes with given fingerprint cheaply, we get
p(n) in polylog.
"""

from sympy import isprime, divisors, primepi, prime as sym_prime
from sympy.ntheory.modular import crt
from pathlib import Path
import math


def tau_via_recursion(n_max):
    """Compute Ramanujan's tau(n) for n=1..n_max using
    multiplicativity + tau(p^{k+1}) = tau(p)*tau(p^k) - p^11 * tau(p^{k-1}).
    Requires tau(p) for primes p <= n_max.

    Use the eta product expansion: Delta(q) = q * prod_{n>=1}(1-q^n)^24.
    Compute tau(n) by power series expansion to n_max.
    """
    # Build (1 - q^k)^24 product for k=1..ceil(sqrt(n_max))
    # Faster: use Euler pentagonal number theorem for prod(1-q^k):
    # prod_{n>=1}(1-q^n) = sum_{k=-inf}^{inf} (-1)^k q^{k(3k-1)/2}
    # Then raise to 24th power.

    # Step 1: compute prod (1 - q^n) coefficients via pentagonal
    eta = [0] * (n_max + 2)
    eta[0] = 1
    k = 1
    while True:
        e1 = k * (3 * k - 1) // 2
        e2 = k * (3 * k + 1) // 2
        if e1 > n_max and e2 > n_max:
            break
        sign = -1 if (k % 2) else 1
        # The series is sum_{k=-inf..inf} (-1)^k q^{k(3k-1)/2}
        # For k=1: pentagonal numbers 1, 2, 5, 7, 12, 15, ...
        # The k=1 term has sign (-1)^1 = -1; e1=1, e2=2
        # The k=2 term has sign (-1)^2 = +1; e1=5, e2=7
        s = -1 if (k % 2 == 1) else 1
        if e1 <= n_max:
            eta[e1] += s
        if e2 <= n_max:
            eta[e2] += s
        k += 1

    # eta[0..n_max] now has coefficients of prod(1-q^n).
    # Need eta^24 = Delta / q. Use repeated squaring.

    def poly_mul(a, b, n):
        out = [0] * (n + 1)
        for i in range(min(len(a), n + 1)):
            if a[i] == 0:
                continue
            ai = a[i]
            for j in range(min(len(b), n + 1 - i)):
                out[i + j] += ai * b[j]
        return out

    # eta^24 via binary exponentiation
    base = eta[: n_max + 1]
    result = [0] * (n_max + 1)
    result[0] = 1
    exp = 24
    while exp > 0:
        if exp & 1:
            result = poly_mul(result, base, n_max)
        base = poly_mul(base, base, n_max)
        exp >>= 1

    # Now Delta(q) = q * eta^24, so tau(n) = result[n-1] for n >= 1.
    tau = [0] * (n_max + 1)
    for n in range(1, n_max + 1):
        tau[n] = result[n - 1]
    return tau


def sigma(n, k):
    """Sum of k-th powers of divisors."""
    return sum(d**k for d in divisors(n))


def main():
    N = 10000
    print("Computing tau(n) for n=1..N...")
    tau = tau_via_recursion(N)

    primes = [n for n in range(2, N + 1) if isprime(n)]
    composites = [n for n in range(2, N + 1) if not isprime(n)]
    pi_N = len(primes)

    # Fingerprint features. Each is a function n -> small integer.
    def feat_tau_691(n):
        # For primes p, tau(p) ≡ p^11 + 1 (mod 691)
        return (tau[n] - pow(n, 11, 691) - 1) % 691

    def feat_tau_5(n):
        # For primes p coprime to 5, tau(p) ≡ 1 + p (mod 5) (Hardy?)
        # Actually: tau(p) ≡ p + p^4 ≡ p + p^4 (mod 5) by Wilton
        # Use direct test
        return tau[n] % 5

    def feat_sigma1_minus_1_minus_n(n):
        # For primes p: sigma_1(p) = 1 + p. So sigma_1(p) - 1 - p = 0.
        return sigma(n, 1) - 1 - n

    def feat_n_mod_30(n):
        # Mod 30 primes are mostly in {1,7,11,13,17,19,23,29}
        return n % 30

    def feat_fermat_mod7(n):
        # By Fermat little, n^6 ≡ 1 (mod 7) for primes != 7. Composites can fail.
        return pow(n, 6, 7) if math.gcd(n, 7) == 1 else 7

    # Combined fingerprint
    def fp(n, depth):
        feats = [
            feat_tau_691(n),
            feat_tau_5(n),
            feat_sigma1_minus_1_minus_n(n),
            feat_n_mod_30(n),
            feat_fermat_mod7(n),
        ]
        return tuple(feats[:depth])

    # For each fingerprint depth k=1..5:
    # - Build set of fingerprints occurring at primes (only)
    # - Count composites whose fingerprint also appears at some prime (false positives)
    # - Count primes whose fingerprint also appears at some composite (collisions)
    rows = []
    for depth in range(1, 6):
        prime_fps = {fp(p, depth) for p in primes}
        composite_fps_in_prime_set = sum(
            1 for c in composites if fp(c, depth) in prime_fps
        )
        # FN: primes that share a fingerprint with a composite
        composite_fps = {fp(c, depth) for c in composites}
        primes_collide = sum(1 for p in primes if fp(p, depth) in composite_fps)
        rows.append(
            {
                "depth": depth,
                "FP": composite_fps_in_prime_set,
                "FP_rate": composite_fps_in_prime_set / len(composites),
                "primes_collide": primes_collide,
                "primes_collide_rate": primes_collide / len(primes),
                "unique_prime_fps": len(prime_fps),
            }
        )

    # Special test: when feat_tau_691(p) == 0 (Ramanujan congruence),
    # how often does this hold for COMPOSITES too?
    primes_passing_ramanujan = sum(1 for p in primes if feat_tau_691(p) == 0)
    composites_passing_ramanujan = sum(1 for c in composites if feat_tau_691(c) == 0)

    # Write results
    out = Path(__file__).with_suffix("").as_posix() + "_results.md"
    with open(out, "w") as f:
        f.write("# Proposal C — Hecke / arithmetic fingerprint as primality oracle\n\n")
        f.write("## Setup\n\n")
        f.write(
            f"For n in [2, {N}] (containing {pi_N} primes, {len(composites)} composites), "
            "compute fingerprints of varying depth combining:\n\n"
            "1. tau(n) - n^11 - 1 (mod 691) -- the Ramanujan congruence (zero on primes)\n"
            "2. tau(n) mod 5\n"
            "3. sigma_1(n) - 1 - n -- zero iff n is prime or 1\n"
            "4. n mod 30\n"
            "5. n^6 mod 7\n\n"
        )

        f.write("## Ramanujan-congruence test (depth 1)\n\n")
        f.write(
            f"- Primes for which tau(p) ≡ p^11 + 1 (mod 691): "
            f"{primes_passing_ramanujan} / {pi_N} = "
            f"{primes_passing_ramanujan/pi_N*100:.2f}%\n"
        )
        f.write(
            f"- Composites for which tau(c) ≡ c^11 + 1 (mod 691): "
            f"{composites_passing_ramanujan} / {len(composites)} = "
            f"{composites_passing_ramanujan/len(composites)*100:.2f}%\n\n"
        )

        f.write("## Fingerprint discrimination as a function of depth\n\n")
        f.write(
            "*FP* = composites whose full fingerprint also occurs at some prime "
            "(false positives, the metric we want to drive to 0).\n"
            "*Collide* = primes whose fingerprint also occurs at some composite.\n\n"
        )
        f.write("| depth | unique prime fps | composite FP | FP rate | primes colliding | collide rate |\n")
        f.write("|---|---|---|---|---|---|\n")
        for r in rows:
            f.write(
                f"| {r['depth']} | {r['unique_prime_fps']} | {r['FP']} | "
                f"{r['FP_rate']*100:.2f}% | {r['primes_collide']} | "
                f"{r['primes_collide_rate']*100:.2f}% |\n"
            )

        f.write("\n## Interpretation\n\n")
        if rows[-1]["FP"] == 0 and rows[-1]["primes_collide"] == 0:
            f.write(
                "With depth 5, the fingerprint **perfectly separates primes "
                "from composites** in [2, 10000]. This is a polylog primality "
                "oracle for n <= 10000. Path remains open; escalate to larger N.\n"
            )
        elif rows[-1]["FP_rate"] < 0.001:
            f.write(
                f"With depth 5, FP rate is {rows[-1]['FP_rate']*100:.3f}% — "
                "very low but not zero. The fingerprint compresses primality "
                "almost-but-not-quite. Open question: does adding more features "
                "linearly close the gap, or do we asymptote?\n\n"
                "**Verdict:** Path partially open; further features needed.\n"
            )
        else:
            f.write(
                f"With depth 5, FP rate is {rows[-1]['FP_rate']*100:.3f}% — too "
                "high for a primality oracle. The fingerprints chosen do not "
                "separate primes from composites.\n\n"
                "Note: feature 3 (sigma_1(n) - 1 - n) IS a perfect discriminator "
                "by definition (zero iff prime), but its computation requires "
                "factoring n, which is *not* polylog. The interesting question is "
                "whether features 1+2 (cheaply computable Hecke data) approach "
                "perfect discrimination.\n\n"
                "**Verdict:** Path closed for the chosen features. The Ramanujan "
                "congruence is necessary for primes but not sufficient — many "
                "composites also satisfy it.\n"
            )

        f.write("\n## Note on cost of features\n\n")
        f.write(
            "- tau(n) per single n: O(polylog n) is **conjectural**. The known "
            "best is Edixhoven et al.'s O(polylog) algorithm for tau(p) at "
            "primes p (assuming GRH for ell-adic Galois reps). For composites, "
            "tau(n) is computed via multiplicativity from tau(p) at the prime "
            "factors of n — but this requires factoring n, breaking polylog.\n"
            "- sigma_k(n): O(sqrt(n)) by enumerating divisors — not polylog.\n"
            "- n mod m: O(polylog).\n\n"
            "So unless tau(n) for *composite* n can be computed without factoring "
            "(open problem, likely false), the entire fingerprint is not polylog "
            "to evaluate. This is a structural barrier even before considering "
            "the discrimination rate.\n"
        )

    print(f"Wrote {out}")
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
