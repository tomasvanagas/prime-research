"""
D48 — Bost-Connes endomotive Galois orbits on χ_P-projected KMS_∞ ground states.

Critic-recommended pick from S163 critique. Computes the Galois-orbit
length distribution on equivalence classes of (Z/NZ)* characters under
the χ_P-projection, comparing to random matched-density null and to
the Galois-only character-order distribution.

Falsification criterion in `bc_endomotive_galois_chi_p_results.md`.
"""

import argparse
import json
import sys
from collections import Counter
from itertools import product
from math import gcd, lcm
from pathlib import Path

import numpy as np
from sympy import factorint, primerange
from sympy.ntheory.modular import crt
from sympy.ntheory.residue_ntheory import primitive_root


# ----------------------------------------------------------------------
# (Z/NZ)* cyclic decomposition + discrete log table
# ----------------------------------------------------------------------

def cyclic_decomp(N):
    """Return list of (g, d): G = ⟨g_1⟩ × ... × ⟨g_r⟩ with g_i ∈ (Z/NZ)*
    of order d_i (CRT-lifted from per-prime-power factors)."""
    factors = factorint(N)
    pp_local = []  # list of (m_local, g_local, ord_local)
    for p, k in factors.items():
        m = p ** k
        if p == 2:
            if k == 1:
                continue
            elif k == 2:
                pp_local.append((m, 3, 2))  # generator -1 mod 4
            else:
                pp_local.append((m, m - 1, 2))             # -1 generator, ord 2
                pp_local.append((m, 5, 2 ** (k - 2)))      # 5 generator, ord 2^(k-2)
        else:
            phi_local = (p - 1) * p ** (k - 1)
            g_local = int(primitive_root(m))
            pp_local.append((m, g_local, phi_local))

    gens = []
    for (m, g_local, ord_local) in pp_local:
        n_complement = N // m
        if n_complement == 1:
            g_global = g_local
        else:
            g_global, _ = crt([m, n_complement], [g_local, 1])
            g_global = int(g_global)
        gens.append((g_global, ord_local))
    return gens


def all_dlogs(N, gens):
    """Map a ∈ (Z/NZ)* → (e_1, ..., e_r) such that a ≡ Π g_i^{e_i} mod N."""
    d_list = [d for _, d in gens]
    g_list = [g for g, _ in gens]
    dlog = {}
    for e_tuple in product(*[range(d) for d in d_list]):
        a = 1
        for g, e in zip(g_list, e_tuple):
            a = (a * pow(g, e, N)) % N
        dlog[a] = e_tuple
    return dlog


# ----------------------------------------------------------------------
# Character sums in Z[ζ_M] (exact zero-detection)
# ----------------------------------------------------------------------

def char_exp(k_tuple, e_tuple, M_factors, M):
    """ω_k(a) = ζ_M^{exp}, where exp = Σ_i (M/d_i) k_i e_i mod M."""
    return sum(Mf * k * e for Mf, k, e in zip(M_factors, k_tuple, e_tuple)) % M


def sum_exact_is_zero(exponents, M):
    """Check whether Σ_e ζ_M^e is exactly zero in Z[ζ_M]. Uses cyclotomic
    polynomial reduction: Σ a_e * ζ_M^e in Z[x]/Φ_M(x) is zero iff
    polynomial reduces to 0.

    Implementation: build polynomial p(x) = Σ a_e x^e in Z[x] of degree
    < M, then reduce mod Φ_M(x). Zero iff reduction is zero polynomial.
    """
    from sympy.polys.specialpolys import cyclotomic_poly
    from sympy import Poly, symbols
    if not exponents:
        return True
    x = symbols('x')
    counts = [0] * M
    for e in exponents:
        counts[e] += 1
    if all(c == 0 for c in counts):
        return True
    p = Poly(counts, x)  # coefficients of x^0, x^1, ..., x^{M-1}? Actually sympy Poly takes coeffs in descending order
    # Actually, build properly:
    from sympy import Add, Mul, Pow, Integer
    # Build polynomial Σ counts[i] * x^i
    p_expr = sum(c * x**i for i, c in enumerate(counts) if c != 0)
    if p_expr == 0:
        return True
    p_poly = Poly(p_expr, x, domain='ZZ')
    Phi_M = Poly(cyclotomic_poly(M, x), x, domain='ZZ')
    _, r = p_poly.div(Phi_M)
    return r.is_zero


_ZERO_THRESHOLD = 1e-9


def sum_exact_value(exponents, M, exact_check=False):
    """Return sum Σ ζ_M^e as a complex number AND a flag for exact-zero status.

    Uses a numerical threshold (1e-9) for is_zero by default. The minimum
    non-zero magnitude of Σ a_e ζ_M^e with a_e ∈ N small is bounded below
    by a Kronecker-style estimate that, for our parameter range
    (|exponents| ≤ 343, M ≤ 60), exceeds the threshold by many orders.
    Set exact_check=True to also do exact polynomial reduction in Z[ζ_M]
    (slow); when both are used, the script asserts they agree."""
    if not exponents:
        return 0.0 + 0.0j, True
    counts = np.zeros(M, dtype=int)
    for e in exponents:
        counts[e] += 1
    z = np.exp(2j * np.pi * np.arange(M) / M)
    val = complex(np.dot(counts, z))
    is_zero_num = abs(val) < _ZERO_THRESHOLD
    if exact_check:
        is_zero_alg = sum_exact_is_zero(exponents, M)
        assert is_zero_num == is_zero_alg, (
            f"Numerical/exact disagree at exponents={exponents} M={M}: "
            f"|val|={abs(val)} num={is_zero_num} alg={is_zero_alg}")
        return val, is_zero_alg
    return val, is_zero_num


# ----------------------------------------------------------------------
# Galois orbits on Ĝ under (Z/MZ)*-action
# ----------------------------------------------------------------------

def galois_orbits(d_list, M):
    """Compute Galois orbits on Ĝ ≅ Π_i Z/d_i under σ_a · k = (a k_i mod d_i).
    Returns dict mapping rep → orbit (as set of tuples), and char_to_orbit dict."""
    Galois_group = [a for a in range(1, M + 1) if gcd(a, M) == 1]
    char_tuples = list(product(*[range(d) for d in d_list]))
    char_to_orbit_id = {}
    orbits = []  # list of frozensets

    for kt in char_tuples:
        if kt in char_to_orbit_id:
            continue
        orbit = set()
        for a in Galois_group:
            new_kt = tuple((a * k) % d for k, d in zip(kt, d_list))
            orbit.add(new_kt)
        orbit_id = len(orbits)
        orbits.append(frozenset(orbit))
        for o in orbit:
            char_to_orbit_id[o] = orbit_id

    return orbits, char_to_orbit_id


# ----------------------------------------------------------------------
# Main D48 measurement
# ----------------------------------------------------------------------

def measure_N(N, n_trials=50, seed=42, verbose=True):
    rng = np.random.default_rng(seed)

    gens = cyclic_decomp(N)
    d_list = [d for _, d in gens]
    g_list = [g for g, _ in gens]

    # Edge case: G is trivial?
    if not d_list:
        return None
    M = lcm(*d_list) if len(d_list) > 1 else d_list[0]
    M_factors = [M // d for d in d_list]

    dlog = all_dlogs(N, gens)
    G = sorted(dlog.keys())
    phi_N = len(G)

    P_N = [p for p in primerange(2, N + 1) if gcd(p, N) == 1]
    pi_prime_N = len(P_N)

    char_tuples = list(product(*[range(d) for d in d_list]))
    char_idx = {kt: j for j, kt in enumerate(char_tuples)}

    def compute_S_set(set_of_a):
        """For each character ω, compute exponents [e_1, ..., e_k] for
        ω applied to elements of `set_of_a`. Return list-of-lists indexed
        by char_idx, plus complex S values, plus exact-zero flags."""
        S_complex = np.zeros(len(char_tuples), dtype=complex)
        S_zero = np.zeros(len(char_tuples), dtype=bool)
        for j, kt in enumerate(char_tuples):
            exps = [char_exp(kt, dlog[a], M_factors, M) for a in set_of_a]
            val, is_zero = sum_exact_value(exps, M)
            S_complex[j] = val
            S_zero[j] = is_zero
        return S_complex, S_zero

    if verbose:
        print(f"\n{'='*60}\nN = {N}: |G| = phi(N) = {phi_N}, |P_N| = {pi_prime_N}, M = exp(G) = {M}, "
              f"r = #cyclic factors = {len(gens)}")
        print(f"  d_list = {d_list}, gens = {g_list}")

    S_chi_complex, S_chi_zero = compute_S_set(P_N)
    n_zero_chi = int(S_chi_zero.sum())

    # Galois orbits on full Ĝ
    orbits, char_to_orbit = galois_orbits(d_list, M)
    n_orbits = len(orbits)

    # For each Galois orbit, count: (a) how many characters in orbit have S_chi = 0,
    # (b) the orbit size.
    orbit_zero_counts = np.zeros(n_orbits, dtype=int)
    for j, kt in enumerate(char_tuples):
        if S_chi_zero[j]:
            orbit_zero_counts[char_to_orbit[kt]] += 1

    orbit_sizes = np.array([len(orbits[i]) for i in range(n_orbits)])

    # Equivalence classes under χ_P projection:
    # - All characters with S_chi != 0 → singleton classes
    # - All characters with S_chi == 0 → one merged "zero class"
    # Galois orbit on equivalence classes: orbit O ⊂ Ĝ contributes:
    # - number of distinct singleton classes = (size of O) - (# zero chars in O)
    # - if any zero chars in O, plus 1 for the zero-class image
    # But the merged zero-class is shared across ALL Galois orbits!
    # So globally: # equivalence classes = (# non-zero chars) + (1 if any zeros) = (|Ĝ| - n_zero_chi) + min(n_zero_chi, 1)
    # Galois acts on these: orbit of a non-zero singleton {ω} is its singleton-image orbit (length = # distinct non-zero ω^a).
    # The zero-class is Galois-invariant (its own orbit of length 1).
    # Within a Galois-orbit O of Ĝ, the eq classes from O are:
    #   - one merged-zero element (if any zero in O)
    #   - the non-zero ω^a as singletons; their Galois-orbit length is (size of O) - (# zeros in O), but those are all in the same eq-class orbit.
    # Wait: Galois orbit on eq-classes of {ω}: if ω is a non-zero singleton,
    # the orbit is the set of eq-classes hit by σ_a · ω. Each σ_a · ω is either
    # nonzero (singleton class) or zero (the merged class). So orbit of {ω} =
    # { {ω^a} : ω^a non-zero, a ∈ G_M } ∪ ({zero-class} if any ω^a in Z).
    #
    # For computing orbit-length DISTRIBUTION: each eq-class lies in exactly one
    # orbit. The zero-class is a single eq-class lying in its own orbit (length 1).
    # Singleton non-zero classes lie in orbits whose length is # distinct non-zero ω^a in the Galois-orbit of ω.
    # That is: orbit-length of {ω} (non-zero) = (size of Galois orbit O of ω) - (# zeros in O).
    # ALL singletons in a given Ĝ-orbit O lie in the SAME eq-class orbit, of length (|O| - n_zero(O)).
    # So eq-class orbit-length distribution:
    # - (|O| - n_zero(O)) appears exactly (|O| - n_zero(O)) times for each O with n_zero(O) < |O|
    # - Plus one orbit of length 1 (the zero-class), iff any character is zero globally.

    # Equivalent: histogram of {(|O| - n_zero(O)) for orbits O} weighted by (|O| - n_zero(O)),
    # plus one element of length 1 if n_zero_total > 0.

    eq_orbit_length_distribution = Counter()
    # For each Galois orbit on Ĝ that has at least one nonzero char:
    for i in range(n_orbits):
        nonzero_in_orbit = orbit_sizes[i] - orbit_zero_counts[i]
        if nonzero_in_orbit > 0:
            # contributes 'nonzero_in_orbit' eq-classes, each in an orbit of length 'nonzero_in_orbit'
            eq_orbit_length_distribution[int(nonzero_in_orbit)] += int(nonzero_in_orbit)
    # Plus zero-class if any
    if n_zero_chi > 0:
        eq_orbit_length_distribution[1] += 1

    # Galois-only orbit-length distribution (baseline = (I) prediction)
    galois_only_distribution = Counter()
    for i in range(n_orbits):
        # Each Galois orbit contributes |O| eq-classes (all singletons = unprojected),
        # all in the same Galois orbit of length |O|.
        galois_only_distribution[int(orbit_sizes[i])] += int(orbit_sizes[i])

    # Generate matched-density null distributions
    null_results = []
    for trial in range(n_trials):
        R = sorted(rng.choice(G, size=pi_prime_N, replace=False).tolist())
        S_R, S_R_zero = compute_S_set(R)
        n_zero_R = int(S_R_zero.sum())

        orbit_R_zero_counts = np.zeros(n_orbits, dtype=int)
        for j, kt in enumerate(char_tuples):
            if S_R_zero[j]:
                orbit_R_zero_counts[char_to_orbit[kt]] += 1

        eq_orbit_R = Counter()
        for i in range(n_orbits):
            nonzero_in_orbit = orbit_sizes[i] - orbit_R_zero_counts[i]
            if nonzero_in_orbit > 0:
                eq_orbit_R[int(nonzero_in_orbit)] += int(nonzero_in_orbit)
        if n_zero_R > 0:
            eq_orbit_R[1] += 1

        null_results.append({
            'n_zero': n_zero_R,
            'eq_distrib': dict(eq_orbit_R),
            'mean_S_squared': float(np.mean(np.abs(S_R) ** 2)),
            'max_abs_S': float(np.max(np.abs(S_R))),
        })

    null_zero_counts = np.array([r['n_zero'] for r in null_results])
    null_mean_S_sq = np.array([r['mean_S_squared'] for r in null_results])
    null_max_abs_S = np.array([r['max_abs_S'] for r in null_results])

    chi_mean_S_sq = float(np.mean(np.abs(S_chi_complex) ** 2))
    chi_max_abs_S = float(np.max(np.abs(S_chi_complex)))

    # z-scores
    def z(x, vec):
        m, s = float(vec.mean()), float(vec.std() + 1e-12)
        return (x - m) / s

    z_n_zero = z(n_zero_chi, null_zero_counts)
    z_mean_S_sq = z(chi_mean_S_sq, null_mean_S_sq)
    z_max_abs_S = z(chi_max_abs_S, null_max_abs_S)

    if verbose:
        print(f"  # characters of order d in Ĝ: {sorted(galois_only_distribution.items())}")
        print(f"  # zero S_chi(ω): {n_zero_chi} (z = {z_n_zero:+.2f} vs null mean "
              f"{null_zero_counts.mean():.2f} ± {null_zero_counts.std():.2f})")
        print(f"  mean |S_chi|² = {chi_mean_S_sq:.4f} (z = {z_mean_S_sq:+.2f} vs null mean "
              f"{null_mean_S_sq.mean():.4f} ± {null_mean_S_sq.std():.4f})")
        print(f"  max |S_chi|  = {chi_max_abs_S:.4f} (z = {z_max_abs_S:+.2f} vs null mean "
              f"{null_max_abs_S.mean():.4f} ± {null_max_abs_S.std():.4f})")
        print(f"  Galois-only orbit-length distribution (I-prediction): "
              f"{sorted(galois_only_distribution.items())}")
        print(f"  χ_P eq-class orbit-length distribution:               "
              f"{sorted(dict(eq_orbit_length_distribution).items())}")
        # Show null mean
        all_lengths = set(galois_only_distribution) | set(eq_orbit_length_distribution)
        for r in null_results:
            all_lengths |= set(r['eq_distrib'])
        avg_null = {L: float(np.mean([r['eq_distrib'].get(L, 0) for r in null_results]))
                    for L in sorted(all_lengths)}
        std_null = {L: float(np.std([r['eq_distrib'].get(L, 0) for r in null_results]))
                    for L in sorted(all_lengths)}
        print(f"  Null mean   eq-class orbit-length distribution:       "
              f"{[(L, avg_null[L]) for L in sorted(all_lengths)]}")
        print(f"  Null std    eq-class orbit-length distribution:       "
              f"{[(L, std_null[L]) for L in sorted(all_lengths)]}")
        # Per-bin z-score
        print(f"  Per-bin z-score (chi - null) / null_std:")
        for L in sorted(all_lengths):
            obs = eq_orbit_length_distribution.get(L, 0)
            zL = (obs - avg_null[L]) / (std_null[L] + 1e-9)
            print(f"    L = {L}: chi = {obs}, null = {avg_null[L]:.2f} ± {std_null[L]:.2f}, "
                  f"z = {zL:+.2f}")

    return {
        'N': N,
        'phi_N': phi_N,
        'pi_prime_N': pi_prime_N,
        'M': M,
        'd_list': d_list,
        'g_list': g_list,
        'P_N': P_N,
        'n_zero_chi': n_zero_chi,
        'mean_S_chi_squared': chi_mean_S_sq,
        'max_abs_S_chi': chi_max_abs_S,
        'galois_only_distribution': dict(galois_only_distribution),
        'chi_eq_distribution': dict(eq_orbit_length_distribution),
        'z_n_zero': float(z_n_zero),
        'z_mean_S_squared': float(z_mean_S_sq),
        'z_max_abs_S': float(z_max_abs_S),
        'null_n_zero_mean': float(null_zero_counts.mean()),
        'null_n_zero_std': float(null_zero_counts.std()),
        'null_mean_S_sq_mean': float(null_mean_S_sq.mean()),
        'null_max_abs_S_mean': float(null_max_abs_S.mean()),
        'null_eq_distrib_summaries': null_results,
    }


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--N',
        type=int,
        nargs='*',
        default=[30, 60, 105, 210, 330, 420, 510, 630, 840, 1155, 1260, 2310],
        help='Moduli to test',
    )
    parser.add_argument('--trials', type=int, default=50)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--out', type=str,
                        default='bc_endomotive_galois_chi_p_results.json')
    parser.add_argument('--quiet', action='store_true')
    args = parser.parse_args()

    out_path = Path(__file__).parent / args.out
    results = []

    for N in args.N:
        try:
            res = measure_N(N, n_trials=args.trials, seed=args.seed,
                            verbose=not args.quiet)
            if res is not None:
                # keep null summaries but compress (don't store full eq_distrib lists)
                null_summary = {
                    'null_n_zero_list': [r['n_zero'] for r in res['null_eq_distrib_summaries']],
                    'null_max_abs_S_list': [r['max_abs_S'] for r in res['null_eq_distrib_summaries']],
                    'null_eq_distribs': [r['eq_distrib'] for r in res['null_eq_distrib_summaries']],
                }
                res_serializable = {k: v for k, v in res.items()
                                    if k != 'null_eq_distrib_summaries'}
                res_serializable.update(null_summary)
                results.append(res_serializable)
        except Exception as e:
            print(f"FAILED N={N}: {e}", file=sys.stderr)
            raise

    out_path.write_text(json.dumps(results, indent=2, default=str))
    print(f"\nWritten {out_path}")


if __name__ == '__main__':
    main()
