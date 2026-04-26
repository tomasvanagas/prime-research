"""
Continued-fraction expansion of π(x) / li(x).

Idea: This ratio approaches 1 as x → ∞. Khinchin's theorem says the geometric
mean of partial quotients of a "random" real number is K_0 ≈ 2.6854...

If π(x)/li(x) has *bounded* or *eventually periodic* CF, its remainder is
quadratic-irrational-like — possibly with hidden algebraic structure exploitable.

If the partial quotients look exponentially distributed with Khinchin geometric
mean ≈ 2.685, the residual is generic and there's no exploitable Diophantine
structure.

Test:
- Compute π(x) for x = 10, 100, ..., 10^6 (sympy)
- Compute li(x) at high precision (mpmath)
- Compute CF expansion of π(x)/li(x) to ~20 partial quotients
- Report:
    a) sequence of partial quotients
    b) geometric mean
    c) max partial quotient (should be unbounded for random reals; bounded for
       CF of algebraic numbers of degree 2)
    d) statistic distribution

Related: π(x)/li(x) − 1 ~ −log(x)/√x · oscillation, so as x grows the CF
should contain longer initial 0s (since ratio → 1 makes a₀ = 1, then ε CF).
"""

from sympy import primepi
from mpmath import mp, mpf, li, log
import math

mp.dps = 80


def continued_fraction(x_mp, max_terms=30):
    """Return list of partial quotients of x_mp."""
    quotients = []
    x = x_mp
    for _ in range(max_terms):
        a = int(mp.floor(x))
        quotients.append(a)
        frac = x - a
        if abs(frac) < mpf("1e-50"):
            break
        x = 1 / frac
    return quotients


def cf_summary(quotients, label):
    """Print summary statistics for CF."""
    qs = quotients
    # Skip trivial leading 1 if present
    print(f"\n{label}")
    print(f"  Full CF (first 25): {qs[:25]}")
    # statistics on a_1, a_2, ... (drop a_0)
    rest = qs[1:]
    if not rest:
        print("  No nontrivial partial quotients.")
        return
    geom_mean = math.exp(sum(math.log(max(1, a)) for a in rest) / len(rest))
    arith_mean = sum(rest) / len(rest)
    print(f"  partial quotients a_1.. count={len(rest)}, geom_mean={geom_mean:.4f}, arith_mean={arith_mean:.4f}")
    print(f"  max a_i = {max(rest)}, min = {min(rest)}")
    print(f"  Khinchin K_0 ≈ 2.6854 (random reals)")
    # Distribution: count a == 1, 2, 3, large
    from collections import Counter
    c = Counter()
    for a in rest:
        if a == 1: c["1"] += 1
        elif a == 2: c["2"] += 1
        elif a == 3: c["3"] += 1
        elif a < 10: c["4-9"] += 1
        elif a < 100: c["10-99"] += 1
        else: c[">=100"] += 1
    print(f"  distribution: {dict(c)}")


def main():
    test_xs = [10, 100, 1000, 10000, 100000, 1000000]
    print("Continued fraction expansion of π(x) / li(x)")
    print("=" * 60)

    all_qs = []
    for x in test_xs:
        pi_x = int(primepi(x))
        li_x = li(mpf(x))
        ratio = mpf(pi_x) / li_x
        diff = ratio - 1
        print(f"\nx = {x}")
        print(f"  π(x) = {pi_x}")
        print(f"  li(x) = {mp.nstr(li_x, 12)}")
        print(f"  π/li = {mp.nstr(ratio, 25)}")
        print(f"  π/li − 1 = {mp.nstr(diff, 12)}")
        qs = continued_fraction(ratio, max_terms=30)
        cf_summary(qs, f"  CF of π({x})/li({x}):")
        all_qs.extend(qs[1:])  # collect nontrivial parts

    # Aggregate
    print("\n" + "=" * 60)
    print("Aggregate CF statistics across all x:")
    print(f"  total partial quotients: {len(all_qs)}")
    if all_qs:
        gm = math.exp(sum(math.log(max(1, a)) for a in all_qs) / len(all_qs))
        print(f"  combined geometric mean: {gm:.4f}  (Khinchin = 2.685)")
        print(f"  largest:  {sorted(all_qs)[-5:]}")

    # Test 2: For fixed x = 100000, compute CF of (π(x) - li(x))/(sqrt(x)/log(x))
    # The Riemann hypothesis predicts this is bounded; CF should reflect.
    print("\n" + "=" * 60)
    print("Bonus: CF of normalized residual r(x) = (π(x) − li(x)) · log(x) / √x")
    print("Under RH, |r(x)| is bounded; CF behavior reflects 'how irrational'.")
    for x in [10000, 100000, 1000000]:
        pi_x = int(primepi(x))
        li_x = li(mpf(x))
        r = (mpf(pi_x) - li_x) * mp.log(x) / mp.sqrt(x)
        print(f"\n  x = {x}: r(x) = {mp.nstr(r, 20)}")
        # CF of |r|
        qs = continued_fraction(abs(r), max_terms=25)
        cf_summary(qs, f"    CF of |r({x})|:")


if __name__ == "__main__":
    main()
