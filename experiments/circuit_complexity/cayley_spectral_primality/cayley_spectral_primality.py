"""
A3 (ATTACK_VECTORS.md §A.A3) — Spectral primality test via custom Cayley graph.

Question: For n coprime to 30, build Cayley graph C_n := Cay((Z/nZ)*, {2,3,5,
2^{-1},3^{-1},5^{-1}}). Does any feature of the Laplacian / adjacency spectrum
of C_n decide primality of n?

Falsification criterion (PRE-DECLARED, before running):
  H_A (positive): there exists a polynomial-time-computable scalar feature
    F(spec(C_n)) such that F(spec(C_n)) is in disjoint ranges for primes
    n and composites n on a held-out validation set [n=600..1000, n coprime
    to 30]. Disjoint = no overlap.
  H_0 (null / closure): every feature we compute either (a) has overlapping
    ranges between primes and composites, or (b) the "feature" is a
    rephrasing of a known primality test (Lucas-like, Fermat-like, order
    test) that already lives in NC^1 / TC^0 / its known complexity class.

Cross-domain ingredient: spectral graph theory of Cayley graphs of finite
abelian groups. For abelian G with character group Ĝ, eigenvalues of
Cay(G,S) are λ_χ = Σ_{s∈S} χ(s) for χ ∈ Ĝ (Babai 1979, "Spectra of Cayley
graphs"). The cyclic-vs-product structure of (Z/nZ)* differs sharply for
primes (cyclic of order n-1) vs composites.

Channelling: Bourgain (spectral gap on Cayley graphs of arithmetic groups,
sum-product mixing arguments) and Lubotzky-Phillips-Sarnak (spectral
features of arithmetic Cayley graphs).

Cross-references:
  E5.1 (BPSW), E5.3 (TC^0 primality requires growing-dim MPOW),
  E5.8 (Brandt diagonalisation closure), E7.10 (AKS modulus orthogonality
  to depth). CLOSED_PATHS lines 354, 383, 384, 385 — but those used
  primes-as-generators (circular). A3 uses fixed generators {2,3,5}.

Sub-experiments:
  1. Build C_n for n in [7..1000] coprime to 30.
  2. Compute adjacency spectrum.
  3. Extract features:
     a. Spectral gap (λ_max - λ_2)
     b. Number of distinct eigenvalues
     c. mult(λ_max)  = [G : ⟨S⟩]
     d. mult(λ_min)
     e. Mean / std / skew of eigenvalue distribution
     f. Number of integer eigenvalues
  4. Plot features vs n, color-coded by primality.
  5. For each feature, compute primal/composite separation:
       max overlap = max{P[F prime] within composite range,
                          P[F composite] within prime range}.
"""

import math
import os
import json
import sys
from collections import Counter

import numpy as np
from sympy import isprime, factorint, totient


def coprime_units(n):
    """Return units mod n: {a in [1,n-1] : gcd(a,n) = 1}."""
    return [a for a in range(1, n) if math.gcd(a, n) == 1]


def modinv(a, n):
    """Modular inverse via extended Euclidean."""
    g, x, _ = ext_gcd(a, n)
    if g != 1:
        raise ValueError(f"no inverse for {a} mod {n}")
    return x % n


def ext_gcd(a, b):
    if b == 0:
        return (a, 1, 0)
    g, x1, y1 = ext_gcd(b, a % b)
    return (g, y1, x1 - (a // b) * y1)


def cayley_adjacency(n, generators=(2, 3, 5)):
    """Build adjacency matrix of Cay((Z/nZ)*, S ∪ S^{-1}) for given generators.

    Returns (A, units) where units is the index → group-element list.
    """
    units = coprime_units(n)
    idx = {u: i for i, u in enumerate(units)}
    m = len(units)
    A = np.zeros((m, m), dtype=np.int8)
    S = set()
    for s in generators:
        if math.gcd(s, n) == 1:
            S.add(s % n)
            S.add(modinv(s, n))
    for i, u in enumerate(units):
        for s in S:
            v = (u * s) % n
            j = idx[v]
            A[i, j] = 1
    return A, units, S


def spectrum_features(n, generators=(2, 3, 5)):
    """Compute eigenvalue features for the Cayley graph of (Z/nZ)*."""
    A, units, S = cayley_adjacency(n, generators)
    m = A.shape[0]
    if m == 0 or len(S) == 0:
        return None
    # symmetric (we added inverses), so eigvalsh
    eigs = np.linalg.eigvalsh(A.astype(np.float64))
    eigs_sorted = np.sort(eigs)[::-1]
    deg = len(S)
    # mult(lambda_max = deg)
    tol = 1e-7
    mult_max = int(np.sum(np.abs(eigs - deg) < tol))
    mult_min = int(np.sum(np.abs(eigs - eigs_sorted[-1]) < tol))
    # spectral gap (deg - lambda_2)
    if m >= 2:
        lambda_2 = float(eigs_sorted[1])
    else:
        lambda_2 = float(eigs_sorted[0])
    gap = deg - lambda_2
    # number of distinct eigenvalues
    rounded = np.round(eigs * 1e6) / 1e6  # tolerance
    n_distinct = len(set(rounded.tolist()))
    # how many are very close to integers
    n_int = int(np.sum(np.abs(eigs - np.round(eigs)) < 1e-6))
    # moments
    mean = float(np.mean(eigs))
    std = float(np.std(eigs))
    skew = float(np.mean((eigs - mean) ** 3) / (std ** 3 + 1e-12))
    kurt = float(np.mean((eigs - mean) ** 4) / (std ** 4 + 1e-12))
    # |lambda| > deg-1 count (near-extremal)
    n_near_top = int(np.sum(eigs > deg - 0.5))
    # |lambda_min|
    lambda_min = float(eigs_sorted[-1])
    # subgroup index ⟨S⟩
    subgroup = generate_subgroup(S, n)
    sg_index = totient(n) // len(subgroup) if len(subgroup) > 0 else 0

    # zero eigenvalue multiplicity
    n_zero = int(np.sum(np.abs(eigs) < 1e-7))
    # near-zero eigenvalues (|λ| < 0.5)
    n_small = int(np.sum(np.abs(eigs) < 0.5))
    # ratio of largest to second largest in absolute value (spectral gap normalized)
    abs_sorted = np.sort(np.abs(eigs))[::-1]
    if m >= 2 and abs_sorted[0] > 1e-9:
        abs_gap_ratio = float(abs_sorted[1] / abs_sorted[0])
    else:
        abs_gap_ratio = 0.0
    # trace of A^4, A^6 (closed walks of length 4, 6)
    eigs_pow_4 = float(np.sum(eigs ** 4))
    eigs_pow_6 = float(np.sum(eigs ** 6))
    # number of eigenvalues that are perfect squares of integers — algebraic structure
    n_perfect_int = int(np.sum(np.abs(eigs ** 2 - np.round(eigs ** 2)) < 1e-6))
    # discrete entropy of |λ_i| / Σ|λ_i|
    abs_eigs = np.abs(eigs)
    if abs_eigs.sum() > 1e-9:
        p_dist = abs_eigs / abs_eigs.sum()
        entropy = float(-np.sum(p_dist * np.log(p_dist + 1e-12)))
    else:
        entropy = 0.0
    # number of unique eigenvalues weighted by (m / mult)
    rounded_2 = np.round(eigs * 1e4) / 1e4
    distinct_with_mult = Counter(rounded_2.tolist())
    max_mult = max(distinct_with_mult.values())
    avg_mult = m / len(distinct_with_mult)

    return {
        "n": n,
        "is_prime": bool(isprime(n)),
        "phi_n": int(totient(n)),
        "factorization": str(factorint(n)),
        "deg": deg,
        "lambda_2": lambda_2,
        "lambda_min": lambda_min,
        "spectral_gap": float(gap),
        "mult_max": mult_max,
        "mult_min": mult_min,
        "n_distinct": n_distinct,
        "n_int_eigs": n_int,
        "mean_eig": mean,
        "std_eig": std,
        "skew_eig": skew,
        "kurt_eig": kurt,
        "n_near_top": n_near_top,
        "sg_index": int(sg_index),
        "matrix_size": m,
        "generators_active": sorted(S),
        "n_zero_eigs": n_zero,
        "n_small_eigs": n_small,
        "abs_gap_ratio": abs_gap_ratio,
        "trace_A4": eigs_pow_4,
        "trace_A6": eigs_pow_6,
        "n_perfect_int": n_perfect_int,
        "spectral_entropy": entropy,
        "max_mult": max_mult,
        "avg_mult": float(avg_mult),
    }


def generate_subgroup(S, n):
    """BFS from identity using generators S in (Z/nZ)*."""
    seen = {1}
    frontier = [1]
    while frontier:
        new = []
        for x in frontier:
            for s in S:
                y = (x * s) % n
                if y not in seen:
                    seen.add(y)
                    new.append(y)
        frontier = new
    return seen


def main(n_max=1000, output_dir=None):
    output_dir = output_dir or os.path.dirname(os.path.abspath(__file__))
    os.makedirs(output_dir, exist_ok=True)

    # n coprime to 30, and n >= 7 (so that {2,3,5} all units)
    ns = [n for n in range(7, n_max + 1) if math.gcd(n, 30) == 1]
    print(f"Computing spectra for {len(ns)} values of n in [7..{n_max}]")
    print(f"  primes: {sum(isprime(n) for n in ns)}")
    print(f"  composites: {sum(not isprime(n) for n in ns)}")

    rows = []
    for i, n in enumerate(ns):
        if i % 50 == 0:
            print(f"  [{i:4d}/{len(ns)}]  n={n:5d}, primes_so_far={sum(r['is_prime'] for r in rows)}")
        feat = spectrum_features(n)
        if feat is not None:
            rows.append(feat)

    out_path = os.path.join(output_dir, "spectrum_features.json")
    with open(out_path, "w") as f:
        json.dump(rows, f, indent=1)
    print(f"\nWrote {len(rows)} rows to {out_path}")

    # quick analysis: separation per feature
    feature_keys = [
        "lambda_2", "spectral_gap", "mult_max", "mult_min", "n_distinct",
        "n_int_eigs", "mean_eig", "std_eig", "skew_eig", "kurt_eig",
        "n_near_top", "sg_index", "lambda_min",
    ]
    primes = [r for r in rows if r["is_prime"]]
    comps = [r for r in rows if not r["is_prime"]]
    print(f"\n{'Feature':<16} {'prime min':>12} {'prime max':>12} {'comp min':>12} {'comp max':>12}  separation?")
    print("-" * 90)
    sep_summary = {}
    for k in feature_keys:
        pvs = np.array([r[k] for r in primes], dtype=np.float64)
        cvs = np.array([r[k] for r in comps], dtype=np.float64)
        p_min, p_max = pvs.min(), pvs.max()
        c_min, c_max = cvs.min(), cvs.max()
        # disjoint?
        disjoint_lo = p_max < c_min
        disjoint_hi = c_max < p_min
        disjoint = disjoint_lo or disjoint_hi
        # overlap fraction (composites in prime range)
        overlap_p_in_c = np.mean((pvs >= c_min) & (pvs <= c_max)) if len(pvs) else 0
        overlap_c_in_p = np.mean((cvs >= p_min) & (cvs <= p_max)) if len(cvs) else 0
        # signed difference
        sep_med = (np.median(pvs) - np.median(cvs)) / (np.std(pvs) + np.std(cvs) + 1e-9)
        marker = " *DISJOINT*" if disjoint else ""
        print(f"{k:<16} {p_min:>12.4f} {p_max:>12.4f} {c_min:>12.4f} {c_max:>12.4f}  med-z={sep_med:+.2f}{marker}")
        sep_summary[k] = {
            "prime_min": float(p_min), "prime_max": float(p_max),
            "comp_min": float(c_min), "comp_max": float(c_max),
            "disjoint": bool(disjoint),
            "overlap_p_in_c": float(overlap_p_in_c),
            "overlap_c_in_p": float(overlap_c_in_p),
            "med_separation_z": float(sep_med),
        }

    sep_path = os.path.join(output_dir, "separation_summary.json")
    with open(sep_path, "w") as f:
        json.dump(sep_summary, f, indent=1)
    print(f"\nSeparation summary written to {sep_path}")

    # also dump per-composite features so we can inspect them
    print("\nComposite-by-composite feature detail:")
    print(f"{'n':>6} {'fac':>20} {'sg_idx':>7} {'mult_max':>9} {'gap':>9} {'lambda_2':>9} {'lambda_min':>11}")
    for r in comps:
        print(f"{r['n']:>6} {r['factorization']:>20} {r['sg_index']:>7d} {r['mult_max']:>9d} {r['spectral_gap']:>9.4f} {r['lambda_2']:>9.4f} {r['lambda_min']:>11.4f}")

    return rows, sep_summary


if __name__ == "__main__":
    n_max = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    main(n_max=n_max)
