"""Proposal 26 follow-up — do Fejer-recovery failures cluster by arithmetic features of x?

The original proposal26 experiment showed that Fejer damping gives a constant-factor
improvement over sharp truncation at intermediate T (e.g. T=100: 67% recovery vs sharp 53%
on x in [100, 3000]). The natural follow-up: if Fejer FAILURES concentrate on a
recognisable subset of inputs (smooth-x, bad residues mod small primes, near a
prime jump), a hybrid algorithm could route "easy" x through Fejer at small T and
"hard" x through sharp truncation at full T, giving a real speedup.

This script measures, for x in [100, 3000] and T in {30, 100, 300}, whether Fejer
failure indicator F(x) is predictable from cheap arithmetic features of x:
  - smoothness: B(x) = largest prime factor of x
  - x mod m for m in {2, 3, 5, 6, 30, 210}
  - distance to nearest prime
  - parity, sigma_0(x) (number of divisors)

If any feature has predictive power well above chance, a hybrid scheme is possible.

Run: python3 proposal26_fejer_failure_clustering.py
"""

from pathlib import Path
import time

import mpmath
from sympy import primepi, factorint, divisor_count, nextprime, prevprime
from sympy.ntheory import mobius


ZEROS_FILE = Path(__file__).resolve().parents[2] / "data" / "zeta_zeros_2000.txt"


def load_zeros(n_max: int = 2000):
    out = []
    with open(ZEROS_FILE) as fh:
        for line in fh:
            line = line.strip()
            if line:
                out.append(mpmath.mpf(line))
                if len(out) >= n_max:
                    break
    return out


def riemann_R(x):
    x = mpmath.mpf(x)
    s = mpmath.mpf(0)
    k = 1
    while True:
        root = x ** (mpmath.mpf(1) / k)
        if root <= mpmath.mpf("1.05"):
            break
        mu = mobius(k)
        if mu != 0:
            s += mpmath.mpf(mu) / k * mpmath.li(root)
        k += 1
        if k > 200:
            break
    return s


def li_via_Ei(rho, log_x):
    return mpmath.ei(rho * log_x)


def explicit_pi_fejer(x, zeros, T, sharp=False):
    x = mpmath.mpf(x)
    log_x = mpmath.log(x)
    s = riemann_R(x)
    for gamma in zeros:
        if sharp:
            if gamma > T:
                break
            K = mpmath.mpf(1)
        else:
            if gamma > 2 * T:
                break
            u = gamma / (2 * T)
            K = (mpmath.sin(u) / u) ** 2 if u else mpmath.mpf(1)
        rho = mpmath.mpc(mpmath.mpf("0.5"), gamma)
        s -= 2 * K * li_via_Ei(rho, log_x).real
    return s


def largest_prime_factor(n):
    f = factorint(n)
    return max(f.keys()) if f else 1


def dist_to_nearest_prime(x):
    if x <= 2:
        return 0
    np = nextprime(x - 1)  # smallest prime >= x
    pp = prevprime(x + 1) if x >= 3 else 2  # largest prime <= x
    return min(abs(np - x), abs(x - pp))


def main():
    mpmath.mp.dps = 30
    print("Loading zeros...")
    zeros = load_zeros(2000)
    print(f"  loaded {len(zeros)} zeros")

    xs = list(range(100, 3001, 50))  # 59 sample points for speed
    print(f"\nSampling {len(xs)} x values in [100, 3000].")

    pi_xs = [int(primepi(x)) for x in xs]

    # Cheap arithmetic features
    feats = {
        "lpf": [largest_prime_factor(x) for x in xs],
        "smooth_le_7": [int(largest_prime_factor(x) <= 7) for x in xs],
        "mod6": [x % 6 for x in xs],
        "mod30": [x % 30 for x in xs],
        "is_prime": [int(divisor_count(x) == 2) for x in xs],
        "div_count": [divisor_count(x) for x in xs],
        "near_prime_dist": [dist_to_nearest_prime(x) for x in xs],
        "x_mod_2": [x % 2 for x in xs],
        "x_mod_3": [x % 3 for x in xs],
    }

    for T in [30, 100, 300]:
        print(f"\n=== T = {T} ===")
        t0 = time.time()
        s_sharp = [explicit_pi_fejer(x, zeros, T, sharp=True) for x in xs]
        s_fejer = [explicit_pi_fejer(x, zeros, T, sharp=False) for x in xs]
        sharp_ok = [round(float(s)) == pi for s, pi in zip(s_sharp, pi_xs)]
        fejer_ok = [round(float(s)) == pi for s, pi in zip(s_fejer, pi_xs)]
        sharp_rate = sum(sharp_ok) / len(xs)
        fejer_rate = sum(fejer_ok) / len(xs)
        print(f"  sharp recovery: {sharp_rate:.3f}, fejer recovery: {fejer_rate:.3f}, eval {time.time()-t0:.1f}s")

        # For each feature, partition xs and measure recovery rate per partition.
        # Look for partitions where one bucket has >= 90% recovery and another <= 50% recovery
        # (clear separability) using Fejer.
        for fname, vals in feats.items():
            buckets = {}
            for v, ok in zip(vals, fejer_ok):
                # Bucket continuous features (lpf, near_prime_dist, div_count) into 5 quintiles.
                if fname in ("lpf", "near_prime_dist", "div_count"):
                    pass  # skip per-value bucketing, do quartile analysis below
                else:
                    buckets.setdefault(v, []).append(ok)
            if fname in ("lpf", "near_prime_dist", "div_count"):
                # Quartile split
                sorted_vals = sorted(zip(vals, fejer_ok))
                n = len(sorted_vals)
                q = [sorted_vals[i*n//4:(i+1)*n//4] for i in range(4)]
                rates = [sum(o for _, o in qi)/len(qi) if qi else float("nan") for qi in q]
                spread = max(rates) - min(rates)
                print(f"  feature={fname:>15}  quartile rates={[f'{r:.2f}' for r in rates]}  spread={spread:.3f}")
            else:
                rates = {k: (sum(v) / len(v) if v else float('nan'), len(v)) for k, v in buckets.items()}
                spread = max(r for r, _ in rates.values()) - min(r for r, _ in rates.values())
                summary = ", ".join(f"{k}->{r:.2f}({n})" for k, (r, n) in sorted(rates.items()))
                print(f"  feature={fname:>15}  spread={spread:.3f}  buckets: {summary}")


if __name__ == "__main__":
    main()
