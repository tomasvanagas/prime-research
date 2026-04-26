"""
Robustness check for A3 finding: do alternative generator sets behave the
same way? If yes, the structural conclusion (spectrum probes omega(n), not
primality) is generator-independent.

Test sets:
  S1 = {2, 3, 5} ∪ inverses  (original)
  S2 = {2, 3, 7} ∪ inverses
  S3 = {2, 3, 5, 7} ∪ inverses (4 base generators, larger degree)
  S4 = {-1, 2, 3, 5} ∪ inverses (include -1 = trivial-ish)

For each, run sweep n in [7..500] coprime to the relevant ideal, count
n_int_eigs vs 2^omega(n). Verify lower bound holds and the spectral
feature does NOT separate primes from prime powers.

Also: run KEY DISCRIMINATION TEST. For EACH n in [7..2000] coprime to 30,
run a primality predictor based on threshold(n_int_eigs) and measure
prime/prime-power error rates. Must be HIGH for prime powers (primality
test fails).
"""

import math
import json
import os
import sys

import numpy as np
from sympy import isprime, factorint
from collections import defaultdict


def coprime_units(n):
    return [a for a in range(1, n) if math.gcd(a, n) == 1]


def ext_gcd(a, b):
    if b == 0:
        return (a, 1, 0)
    g, x1, y1 = ext_gcd(b, a % b)
    return (g, y1, x1 - (a // b) * y1)


def modinv(a, n):
    g, x, _ = ext_gcd(a, n)
    if g != 1:
        return None
    return x % n


def cayley_spectrum(n, generators):
    units = coprime_units(n)
    idx = {u: i for i, u in enumerate(units)}
    m = len(units)
    A = np.zeros((m, m), dtype=np.int8)
    S = set()
    for s in generators:
        if math.gcd(s, n) == 1:
            S.add(s % n)
            inv = modinv(s, n)
            if inv is not None:
                S.add(inv)
    if len(S) == 0:
        return None, S
    for i, u in enumerate(units):
        for s in S:
            j = idx[(u * s) % n]
            A[i, j] = 1
    eigs = np.linalg.eigvalsh(A.astype(np.float64))
    return eigs, S


def n_int_eigs(eigs, tol=1e-6):
    return int(np.sum(np.abs(eigs - np.round(eigs)) < tol))


def main(n_max=1500):
    # We focus on n coprime to 30 (so {2,3,5} all units) AND coprime to 7 (so {2,3,7} all units)
    # Just use n coprime to 2*3*5*7 = 210 to allow all generator sets.
    ns = [n for n in range(11, n_max + 1) if math.gcd(n, 210) == 1]
    print(f"Sweeping {len(ns)} values of n coprime to 210 in [11..{n_max}]")

    gen_sets = {
        "S1_235": (2, 3, 5),
        "S2_237": (2, 3, 7),
        "S3_2357": (2, 3, 5, 7),
    }

    results = []
    for n in ns:
        row = {
            "n": n,
            "is_prime": bool(isprime(n)),
            "factorization": str(factorint(n)),
            "omega": len(factorint(n)),
        }
        for name, gen in gen_sets.items():
            eigs, S = cayley_spectrum(n, gen)
            if eigs is None:
                continue
            row[f"{name}_int"] = n_int_eigs(eigs)
            row[f"{name}_size"] = len(S)
        results.append(row)

    out_dir = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(out_dir, "robustness_features.json"), "w") as f:
        json.dump(results, f, indent=1)

    # Verify lower bound 2^omega for each generator set
    print("\nLower bound check (n_int_eigs >= 2^omega):")
    for name in gen_sets:
        viols = 0
        for r in results:
            if r.get(f"{name}_int", 0) < 2 ** r["omega"]:
                viols += 1
        print(f"  {name}: violations = {viols} / {len(results)}")

    # Group analysis
    print("\nDistribution by omega and is_prime:")
    print(f"{'gen':>10} {'group':>20} {'count':>6} {'min':>6} {'p25':>6} {'med':>6} {'p75':>6} {'max':>6}")
    for name in gen_sets:
        groups = defaultdict(list)
        for r in results:
            key = (r["omega"], r["is_prime"])
            groups[key].append(r.get(f"{name}_int", 0))
        for key in sorted(groups.keys()):
            vals = sorted(groups[key])
            label = f"omega={key[0]},prime={key[1]}"
            n = len(vals)
            mn, p25, med, p75, mx = (vals[0], vals[n//4], vals[n//2], vals[3*n//4], vals[-1])
            print(f"{name:>10} {label:>20} {n:>6d} {mn:>6d} {p25:>6d} {med:>6d} {p75:>6d} {mx:>6d}")

    # KEY DISCRIMINATION TEST: can ANY threshold on n_int_eigs / matrix
    # separate primes from prime powers?
    print("\nDiscrimination test: primes vs prime-powers (omega=1, composite):")
    primes_only = [r for r in results if r["is_prime"]]
    pp_only = [r for r in results if not r["is_prime"] and r["omega"] == 1]
    print(f"  primes: {len(primes_only)}, prime powers: {len(pp_only)}")
    for name in gen_sets:
        if not primes_only or not pp_only:
            continue
        p_ints = sorted([r[f"{name}_int"] for r in primes_only])
        pp_ints = sorted([r[f"{name}_int"] for r in pp_only])
        print(f"  {name}: prime int range [{min(p_ints)}, {max(p_ints)}], pp int range [{min(pp_ints)}, {max(pp_ints)}]; ranges OVERLAP: {min(p_ints) <= max(pp_ints) and min(pp_ints) <= max(p_ints)}")
        # Maximum AUC
        from scipy import stats
        u, p = stats.mannwhitneyu(p_ints, pp_ints)
        auc = u / (len(p_ints) * len(pp_ints))
        auc_disc = max(auc, 1 - auc)
        print(f"    Mann-Whitney AUC for primes vs prime-powers: {auc_disc:.3f} (0.5 = no signal)")

    return results


if __name__ == "__main__":
    n_max = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
    main(n_max)
