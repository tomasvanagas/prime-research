"""
D44 — Rédei symbol distribution on admissible prime triples.

Cross-domain import: Mazur 1966 primes-as-knots dictionary;
Morishita 2002 Massey products on prime numbers; Rédei 1939
explicit ternary symbol formula.

Setup:
  * Enumerate primes p < q < r <= N with p ≡ q ≡ r ≡ 1 mod 4 and
    all pairwise Legendre symbols = +1 ("Borromean admissibility").
  * Compute the Rédei symbol [p, q, r] ∈ {±1} via the explicit
    formula (Lemmermeyer 2000, Reciprocity Laws, Ch. 9):

      Find a, b ∈ Z with a² - p b² = q (genus-principal class).
      Let σ ∈ F_r with σ² ≡ p (mod r) (exists since (p/r) = +1).
      Then [p, q, r] = ((a + b σ) mod r / r)  -- Legendre symbol.

  * Compare empirical distribution of [p,q,r] to:
      (i) uniform ±1 (null hypothesis);
      (ii) HL ternary singular series weighted prediction.

Falsification statement (PRE-REGISTERED before run):

  F1 (A-grade, target):  fraction with [p,q,r] = +1, denoted f_+, deviates
                          from 0.5 by more than 5σ_K = 5 / sqrt(K) where
                          K is the number of admissible triples, AND
                          the deviation is correlated with HL S_3 weighting
                          at Pearson r > 0.3 (50 random-permutation
                          significance test, Bonferroni-corrected).
  F2 (B-grade case I):   |f_+ - 0.5| < 5/sqrt(K) at all three N levels —
                          Rédei symbol is unbiased on admissible triples
                          (closes "arithmetic Massey is bilinear-equivalent
                          to HL on average").
  F3 (B-grade case II):  |f_+ - 0.5| > 5/sqrt(K) but f_+ ≠ HL S_3
                          prediction (deviation has a non-HL structural
                          source). Inspect mod-8, mod-16 residue
                          subsets for substructure.
  INC:                    Cannot find norm representation a² - p b² = q
                          for a substantial fraction (> 10%) of admissible
                          (p, q) pairs with (p/q) = +1. Note: this means
                          q is not represented by the principal form
                          x² - p y² but lies in a non-principal genus
                          class. We will skip such pairs (filing
                          INCONCLUSIVE for them) and report the fraction.

Validation:
  Single test case from Wikipedia / Morishita 2012:
    [13, 61, 937] = -1 (Borromean prime triple — pairwise Legendre +1,
                         but jointly linked).
"""

from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import dataclass

from sympy import isprime, primerange
from sympy.ntheory.residue_ntheory import sqrt_mod
from sympy.ntheory.residue_ntheory import is_quad_residue
from sympy.solvers.diophantine.diophantine import diop_DN


# -----------------------------------------------------------
# Number-theoretic primitives
# -----------------------------------------------------------

def jacobi(a: int, n: int) -> int:
    """Jacobi symbol (a/n).  n > 0 odd."""
    a %= n
    result = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                result = -result
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            result = -result
        a %= n
    return result if n == 1 else 0


def legendre(a: int, p: int) -> int:
    """Legendre symbol (a/p), p odd prime.  Returns -1, 0, +1."""
    return jacobi(a, p)


def find_norm_rep(p: int, q: int):
    """
    Find primitive (a, b) with a^2 - p b^2 = q (in Z[sqrt(p)]).

    Returns (a, b) with a, b >= 0, gcd(a, b) = 1, or None.
    """
    sols = diop_DN(p, q)
    if not sols:
        return None
    # diop_DN returns fundamental positive solutions
    best = None
    for (a, b) in sols:
        a, b = abs(int(a)), abs(int(b))
        if math.gcd(a, b) == 1:
            if best is None or (a + b) < (best[0] + best[1]):
                best = (a, b)
    return best


def find_norm_rep_4q(p: int, q: int):
    """
    Find (U, V) with U^2 - p V^2 = 4q AND U == V (mod 2)
    (ie, an element of the maximal order O_p = Z[(1 + sqrt(p))/2]
    of norm q).  Returns (U, V) primitive (gcd(U, V) = 1 when V > 0)
    or None.
    """
    sols = diop_DN(p, 4 * q)
    if not sols:
        return None
    best = None
    for (U, V) in sols:
        U, V = abs(int(U)), abs(int(V))
        if (U % 2) != (V % 2):
            continue
        if best is None or (U + V) < (best[0] + best[1]):
            best = (U, V)
    return best


# -----------------------------------------------------------
# Rédei symbol
# -----------------------------------------------------------

def redei_symbol(p: int, q: int, r: int, max_b: int = 0):
    """
    Compute the Rédei symbol [p, q, r] ∈ {+1, -1, 0}.

    Preconditions: p, q, r distinct odd primes ≡ 1 mod 4 with all
    pairwise Legendre symbols = +1.

    Returns:
      +1 or -1 — well-defined symbol;
       0       — INCONCLUSIVE (q lies outside the principal narrow
                  ideal class of Z[√p], symbol requires algebraic-
                  number machinery beyond the simple norm-rep formula).
       None    — preconditions not met.

    Formula (Lemmermeyer 2000, Reciprocity Laws, Ch. 9.1; specialized
    to maximal order Z[(1 + sqrt(p))/2] for p ≡ 1 mod 4):

      Find α ∈ O_p with N(α) = q, written as α = (U + V√p)/2 with
      U ≡ V (mod 2), so U^2 - p V^2 = 4q.

      Let σ ∈ F_r with σ^2 ≡ p (mod r).  Define
        α_r := (U + V·σ) · 2^{-1}  in F_r.
      Then [p, q, r] = (α_r / r)  -- Legendre symbol in F_r.

    The 'q-rep' formulation a² - p b² = q is the special case
    U = 2a, V = 2b (both even).  We unify both cases via the 4q-rep.
    """
    if r == p or r == q or p == q:
        return None
    if legendre(p, q) != 1 or legendre(p, r) != 1 or legendre(q, r) != 1:
        return None

    # Try the maximal-order representation 4q (covers both
    # integer-rep and half-integer-rep cases).
    rep = find_norm_rep_4q(p, q)
    if rep is None:
        # No element of O_p has norm q ⇒ q-ideal is non-principal in
        # the narrow class group.  Mark INCONCLUSIVE.
        return 0
    U, V = rep

    # σ² ≡ p (mod r);  exists because (p/r) = +1.
    sigma = sqrt_mod(p % r, r)
    if sigma is None:
        return None
    inv2 = pow(2, -1, r)
    val = ((U + V * sigma) * inv2) % r
    if val == 0:
        return None
    return legendre(val, r)


# -----------------------------------------------------------
# Admissible triples
# -----------------------------------------------------------

def primes_1_mod_4(N: int) -> list[int]:
    return [p for p in primerange(5, N + 1) if p % 4 == 1]


def enumerate_admissible_triples(N: int) -> list[tuple[int, int, int]]:
    """
    Triples (p, q, r) with p < q < r ≤ N, all p,q,r ≡ 1 mod 4, all
    pairwise Legendre symbols = +1.
    """
    primes = primes_1_mod_4(N)
    # precompute pairwise Legendre symbols
    legendre_cache: dict[tuple[int, int], int] = {}
    for i, p in enumerate(primes):
        for q in primes[i + 1:]:
            legendre_cache[(p, q)] = legendre(p, q)
            legendre_cache[(q, p)] = legendre(q, p)  # symmetric for both ≡ 1 mod 4

    out: list[tuple[int, int, int]] = []
    for i, p in enumerate(primes):
        for j in range(i + 1, len(primes)):
            q = primes[j]
            if legendre_cache[(p, q)] != 1:
                continue
            for k in range(j + 1, len(primes)):
                r = primes[k]
                if legendre_cache[(p, r)] != 1:
                    continue
                if legendre_cache[(q, r)] != 1:
                    continue
                out.append((p, q, r))
    return out


# -----------------------------------------------------------
# Hardy-Littlewood ternary singular series
# -----------------------------------------------------------

def hl_ternary_singular_series(N: int) -> float:
    """
    Hardy-Littlewood prime k=3-tuple constant for *generic* admissible
    triples (used as a normalisation reference, NOT as a per-triple
    weight).  S_3 = ∏_p (1 - 3/p + 3/p² - νₚ/p³)/(1 - 1/p)³  -- but
    here we just compute the "uniform admissible" product
    ∏_{p ≤ small_cutoff} (1 - 3/p)/(1 - 1/p)³ * correction.

    For our purposes we use the simpler 3-prime constant from
    Hardy-Littlewood:  C_3 = 2 ∏_{p≥3} (1 − 3p − 1)/( (p − 1)^3) ≈ 1.32...
    For uniform admissible triples (no constraints on prime gaps), the
    relevant prediction is just the FRACTION expected to have
    [p,q,r]=+1 if the symbol is governed by HL pair-correlation.
    """
    # Use the well-known Hardy-Littlewood twin-prime constant cubed
    # as a coarse reference; not directly used for per-triple weights.
    return 0.660161815846869573927812110014  # twin-prime constant


# -----------------------------------------------------------
# Statistics
# -----------------------------------------------------------

@dataclass
class StatsResult:
    N: int
    n_primes_1mod4: int
    n_admissible: int
    n_well_defined: int       # excluding INCONCLUSIVE = 0
    n_inconclusive: int        # representation not found
    n_plus: int
    n_minus: int
    f_plus: float
    z_uniform: float           # z-score against null f_+ = 1/2
    p_value_uniform: float     # 2-sided p-value
    elapsed_sec: float


def compute_stats(triples: list[tuple[int, int, int]],
                  N: int, max_b: int) -> StatsResult:
    t0 = time.time()
    n_plus = 0
    n_minus = 0
    n_inconclusive = 0
    primes_set = set(primes_1_mod_4(N))
    for (p, q, r) in triples:
        s = redei_symbol(p, q, r)
        if s == 1:
            n_plus += 1
        elif s == -1:
            n_minus += 1
        elif s == 0:
            n_inconclusive += 1
        # None should not happen given enumeration filter
    K = n_plus + n_minus
    f_plus = n_plus / K if K > 0 else 0.5
    sigma_uniform = 0.5 / math.sqrt(K) if K > 0 else 1.0
    z = (f_plus - 0.5) / sigma_uniform if K > 0 else 0.0
    # 2-sided z->p
    from math import erf
    p_value = math.erfc(abs(z) / math.sqrt(2))
    return StatsResult(
        N=N,
        n_primes_1mod4=len(primes_set),
        n_admissible=len(triples),
        n_well_defined=K,
        n_inconclusive=n_inconclusive,
        n_plus=n_plus,
        n_minus=n_minus,
        f_plus=f_plus,
        z_uniform=z,
        p_value_uniform=p_value,
        elapsed_sec=time.time() - t0,
    )


# -----------------------------------------------------------
# Mod-residue subset analysis (looks for non-HL substructure)
# -----------------------------------------------------------

def mod_subset_analysis(triples: list[tuple[int, int, int]],
                        n_plus_minus: dict[tuple[int, int, int], int],
                        m: int) -> dict[tuple, dict]:
    """
    For each residue triple (a_p, a_q, a_r) mod m with a_i ∈ {1, 5, ...}
    compute fraction with [p,q,r]=+1.  Returns dict keyed by residue triple.
    """
    bucket: dict[tuple, list[int]] = {}
    for tri in triples:
        sym = n_plus_minus.get(tri, 0)
        if sym == 0:
            continue
        key = (tri[0] % m, tri[1] % m, tri[2] % m)
        bucket.setdefault(key, []).append(sym)
    result = {}
    for key, syms in bucket.items():
        n = len(syms)
        np_ = sum(1 for s in syms if s == 1)
        result[key] = {
            "n": n,
            "n_plus": np_,
            "f_plus": np_ / n if n > 0 else None,
        }
    return result


# -----------------------------------------------------------
# Validation
# -----------------------------------------------------------

def validate_borromean():
    """The Wikipedia / Morishita 2012 canonical example."""
    p, q, r = 13, 61, 937
    s = redei_symbol(p, q, r)
    print(f"Validation [13, 61, 937] = {s} (expected -1)")
    if s == -1:
        print("  ✓ Validation PASSED")
        return True
    else:
        print("  ✗ Validation FAILED")
        return False


# -----------------------------------------------------------
# Main
# -----------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=1000,
                        help="prime cutoff")
    parser.add_argument("--max-b", type=int, default=200_000,
                        help="max search bound for norm representation")
    parser.add_argument("--validate-only", action="store_true",
                        help="only validate on (13, 61, 937)")
    parser.add_argument("--out-json", type=str, default=None)
    args = parser.parse_args()

    print("=" * 60)
    print("D44 — Rédei symbol distribution on prime triples")
    print("=" * 60)
    print(f"N = {args.N}, max_b = {args.max_b}")
    print()

    print("Step 1: Validation on canonical Borromean triple (13, 61, 937)")
    valid = validate_borromean()
    print()
    if args.validate_only:
        return

    print(f"Step 2: Enumerate admissible triples up to N = {args.N}")
    t0 = time.time()
    triples = enumerate_admissible_triples(args.N)
    print(f"  found {len(triples)} admissible triples"
          f" in {time.time() - t0:.2f}s")
    print()

    print("Step 3: Compute Rédei symbol distribution")
    # Compute and bucket per-triple
    sym_map: dict[tuple[int, int, int], int] = {}
    n_plus = 0
    n_minus = 0
    n_inc = 0
    t0 = time.time()
    for i, (p, q, r) in enumerate(triples):
        s = redei_symbol(p, q, r, max_b=args.max_b)
        if s == 1:
            n_plus += 1
            sym_map[(p, q, r)] = 1
        elif s == -1:
            n_minus += 1
            sym_map[(p, q, r)] = -1
        elif s == 0:
            n_inc += 1
        if (i + 1) % 500 == 0:
            print(f"  processed {i + 1}/{len(triples)}  "
                  f"(+{n_plus}, -{n_minus}, inc={n_inc})  "
                  f"[{time.time() - t0:.1f}s]")
    elapsed = time.time() - t0
    print(f"  done in {elapsed:.1f}s")
    print()

    K = n_plus + n_minus
    f_plus = n_plus / K if K > 0 else 0.5
    sigma_unif = 0.5 / math.sqrt(K) if K > 0 else 1.0
    z = (f_plus - 0.5) / sigma_unif if K > 0 else 0.0
    p_val = math.erfc(abs(z) / math.sqrt(2))

    print("Step 4: Summary statistics")
    print(f"  admissible triples           : {len(triples)}")
    print(f"  well-defined ([p,q,r] ≠ 0)  : {K}")
    print(f"  inconclusive (no norm rep)  : {n_inc}")
    print(f"  n_plus  ([p,q,r] = +1)       : {n_plus}")
    print(f"  n_minus ([p,q,r] = -1)       : {n_minus}")
    print(f"  f_plus = n_plus / K          : {f_plus:.6f}")
    print(f"  σ_uniform = 0.5/sqrt(K)      : {sigma_unif:.6f}")
    print(f"  z-score (vs uniform 1/2)     : {z:.4f}")
    print(f"  p-value (two-sided)          : {p_val:.6e}")
    print()

    print("Step 5: Mod-residue substructure analysis")
    for m in (8, 12, 24):
        buckets = mod_subset_analysis(triples, sym_map, m)
        # print top-3 most-deviated buckets
        items = []
        for key, info in buckets.items():
            n = info["n"]
            if n < 20:
                continue
            f = info["f_plus"]
            sig = (f - 0.5) / (0.5 / math.sqrt(n))
            items.append((abs(sig), key, info, sig))
        items.sort(reverse=True)
        print(f"  mod {m}: {len(items)} buckets with n ≥ 20")
        for abs_sig, key, info, sig in items[:5]:
            print(f"    {key} : n={info['n']}, "
                  f"f_+={info['f_plus']:.3f}, z={sig:+.2f}")
    print()

    # Save raw results
    out = {
        "N": args.N,
        "max_b": args.max_b,
        "validation_borromean": valid,
        "n_admissible": len(triples),
        "n_well_defined": K,
        "n_inconclusive": n_inc,
        "n_plus": n_plus,
        "n_minus": n_minus,
        "f_plus": f_plus,
        "sigma_uniform": sigma_unif,
        "z_uniform": z,
        "p_value_uniform": p_val,
        "elapsed_sec": elapsed,
        "triples_full": [
            {"p": p, "q": q, "r": r, "sym": sym_map[(p, q, r)]}
            for (p, q, r) in triples if (p, q, r) in sym_map
        ],
    }
    out_path = (args.out_json or
                f"/apps/aplikacijos/prime-research/experiments/algebraic/"
                f"redei_symbol_prime_triples/results_N{args.N}.json")
    with open(out_path, "w") as fh:
        json.dump({k: v for k, v in out.items() if k != "triples_full"}, fh,
                  indent=2)
    print(f"Saved summary to {out_path}")


if __name__ == "__main__":
    main()
