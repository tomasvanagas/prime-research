"""
A7 commit thread, S211: plethysm sub-frame attack on f_chi_P^(n).

This script implements GL_n irrep decomposition machinery in SymPy and
attacks the open plethysm sub-frame of A7 (Mulmuley-Sohoni GCT).

Strategy
--------
1. Compute the GL_n character of any polynomial-functor representation
   on V = C^n by evaluating the trace at a generic diagonal element
   diag(x_1, ..., x_n) and expanding in the Schur basis {s_lambda} via
   Vandermonde inversion.

2. Three GL_n-modules of interest at this stage:
   - Sym^k(Sym^d V), the ambient ring of degree-(kd) elements in
     polynomials of degree d (the natural Plethysm coefficient table).
   - T_f := span{X . f : X in gl_n}, the GL_n-tangent space of the
     orbit at f.  As a sub-GL_n-module of Sym^d(V*) (taken degree by
     degree).  For a polynomial f in n variables, X . f means the
     Lie-derivative.
   - The deg-1 part of the coordinate ring of orbit closure,
     V* (Sym^d V) modulo the tangent direction; since Sym^d V is
     irreducible this is a quotient of an irreducible module.

3. Compare f_chi_P^(n) against:
   - matched-support random baselines (already done at orbit-dim level
     in S156; we now extend to plethysm level)
   - elementary symmetric e_d (= prod sum_(|S|=d) prod_S x_i)
   - power sum p_d
   - determinantal "comparison" object (we use det_2 for n=4 below)

The session 1 deliverables:
- Plethysm coefficient table for (n=4, d in {1,2,3}, k in {1,2,3,4}).
- Tangent-space irrep decomposition of f_chi_P^(n), e_d, p_d for n=4
  and matched-support random baselines.
- Honest report of whether anything distinguishes f_chi_P from baseline
  at the plethysm sub-frame.

What would falsify the eventual closure:
- Finding a partition lambda such that the multiplicity of S_lambda in
  C[orbit closure of f_chi_P]_d is non-zero but the same multiplicity
  in C[orbit closure of (random-coef matched-support)]_d is zero.
- Or vice-versa.
- Or finding a partition lambda for which the tangent-space
  multiplicity differs between f_chi_P and matched baseline.

Cross-domain refs (cited in session synthesis):
  Mulmuley-Sohoni 2001 SIAM J Comput 31, 496;
  Burgisser-Ikenmeyer-Panova 2017 FOCS arXiv:1604.06431;
  Macdonald 1995 'Symmetric Functions and Hall Polynomials' OUP;
  Stanley 1999 'Enumerative Combinatorics II' ch. 7 plethysm.
"""

import itertools
import json
import os
import random
import sys
from fractions import Fraction
from functools import lru_cache

import sympy as sp


# ============================================================
# Partition utilities
# ============================================================

def partitions_of(N, max_parts=None, max_part=None):
    """Yield partitions of N as decreasing tuples."""
    if N == 0:
        yield ()
        return
    if max_part is None:
        max_part = N
    def gen(remaining, cap, parts_so_far):
        if remaining == 0:
            yield tuple(parts_so_far)
            return
        upper = min(remaining, cap)
        if max_parts is not None and len(parts_so_far) >= max_parts:
            return
        for p in range(upper, 0, -1):
            yield from gen(remaining - p, p, parts_so_far + [p])
    yield from gen(N, max_part, [])


def transpose_partition(lam):
    """Conjugate partition (reflect Young diagram across diagonal)."""
    if not lam:
        return ()
    L = lam[0]
    return tuple(sum(1 for p in lam if p >= i) for i in range(1, L + 1))


# ============================================================
# Schur polynomial via Jacobi-Trudi (in finite n variables)
# ============================================================

@lru_cache(maxsize=None)
def schur_polynomial_in_n_vars(lam, n, vars_tuple):
    """
    Compute s_lambda(x_1, ..., x_n) using bialternant formula
        s_lam = det(x_i^(lam_j + n - j)) / det(x_i^(n - j))
    Returns a sympy expression.
    """
    if len(lam) > n:
        return sp.Integer(0)
    # Pad lam to length n
    L = list(lam) + [0] * (n - len(lam))
    syms = vars_tuple
    # Numerator: det of [x_i ^ (L_j + n - j)] for i,j = 1..n
    # j is column, i is row, with j-th column exponent L_j + n - j
    num_mat = sp.Matrix(n, n,
        lambda i, j: syms[i] ** (L[j] + n - 1 - j))
    # Denominator: Vandermonde det of [x_i ^ (n - j)]
    den_mat = sp.Matrix(n, n,
        lambda i, j: syms[i] ** (n - 1 - j))
    num = num_mat.det()
    den = den_mat.det()
    return sp.simplify(sp.cancel(num / den))


# ============================================================
# GL_n character expansion: write a symmetric polynomial in n vars
# as a linear combination of Schur polynomials s_lambda(x_1,...,x_n)
# for partitions lambda with at most n parts.
# ============================================================

def expand_in_schur_basis(P, syms, total_degree, max_parts):
    """
    Expand a symmetric polynomial P(x_1,...,x_n) (homogeneous of total
    degree D) in the Schur basis.  Returns dict {lambda: integer multiplicity}.

    Uses Vandermonde-quotient: P = sum_lam c_lam s_lam, so
    P * Vandermonde = sum_lam c_lam det(x_i^(lam_j + n - j)).
    Read off c_lam by extracting the coefficient of the monomial
    x_1^(lam_1+n-1) x_2^(lam_2+n-2) ... x_n^(lam_n+0) in
    antisymmetrized P (= P * V).  But we just compute s_lam in the
    n vars and solve the linear system for the small set of partitions.

    Simpler & robust: enumerate partitions of total_degree with
    at most max_parts parts, compute s_lam, and solve the linear
    system over Q.
    """
    n = len(syms)
    parts_list = [lam for lam in partitions_of(total_degree, max_parts=n)
                  if (max_parts is None or len(lam) <= max_parts)]
    schur_expansions = []
    monomial_dict_target = sp.Poly(P, *syms).as_dict()
    # Build a list of monomial keys appearing across all partitions
    all_keys = set(monomial_dict_target.keys())
    for lam in parts_list:
        s_lam = schur_polynomial_in_n_vars(lam, n, tuple(syms))
        sd = sp.Poly(sp.expand(s_lam), *syms).as_dict() if s_lam != 0 else {}
        all_keys.update(sd.keys())
        schur_expansions.append((lam, sd))
    # Solve linear system: target = sum_lam c_lam * s_lam, in monomial basis.
    keys = sorted(all_keys)
    # Build matrix (rows = monomials, cols = partitions)
    A = sp.zeros(len(keys), len(parts_list))
    b = sp.zeros(len(keys), 1)
    for i, k in enumerate(keys):
        b[i, 0] = monomial_dict_target.get(k, 0)
        for j, (lam, sd) in enumerate(schur_expansions):
            A[i, j] = sd.get(k, 0)
    # Solve A x = b for x (the multiplicities).
    sol = A.solve(b)
    out = {}
    for j, (lam, _) in enumerate(schur_expansions):
        v = sp.nsimplify(sol[j, 0])
        if v != 0:
            out[lam] = v
    return out


# ============================================================
# Sym^k(Sym^d V) plethysm: GL_n character
# ============================================================

def sym_d_monomials(syms, d):
    """List of all monomials of degree d in the given variables."""
    n = len(syms)
    out = []
    for exps in _weak_compositions(d, n):
        m = sp.Integer(1)
        for var, e in zip(syms, exps):
            m *= var ** e
        out.append(m)
    return out


def _weak_compositions(d, n):
    """Tuples of length n of non-neg ints summing to d."""
    if n == 1:
        yield (d,)
        return
    for k in range(d + 1):
        for rest in _weak_compositions(d - k, n - 1):
            yield (k,) + rest


def sym_k_of_sym_d_character(n, d, k):
    """
    GL_n character of Sym^k(Sym^d V) where V = C^n.
    This is the elementary symmetric h_k applied to the monomials of
    Sym^d V (as if those monomials were new variables).
    Returns a sympy expression in x_1,...,x_n.
    """
    syms = sp.symbols(f"x1:{n+1}")
    monos = sym_d_monomials(syms, d)
    # h_k(M) = sum over multisets of size k from M of product
    expr = sp.Integer(0)
    for combo in itertools.combinations_with_replacement(range(len(monos)), k):
        term = sp.Integer(1)
        for idx in combo:
            term *= monos[idx]
        expr += term
    return sp.expand(expr), syms


def plethysm_sym_k_sym_d(n, d, k):
    """
    Decompose Sym^k(Sym^d V_n) = sum_lambda b_lambda S_lambda(V_n)
    Returns dict {lambda: multiplicity}.
    Total degree of each Schur is k*d; partitions must have <= n parts.
    """
    char, syms = sym_k_of_sym_d_character(n, d, k)
    return expand_in_schur_basis(char, syms, k * d, max_parts=n)


# ============================================================
# Tangent-space irrep decomposition
# ============================================================

def homogeneous_components(f_expr, syms):
    """Return dict {degree: homogeneous component of f_expr}."""
    poly = sp.Poly(sp.expand(f_expr), *syms)
    out = {}
    for monom, coef in poly.as_dict().items():
        deg = sum(monom)
        if deg in out:
            out[deg] += coef * sp.prod(s ** e for s, e in zip(syms, monom))
        else:
            out[deg] = coef * sp.prod(s ** e for s, e in zip(syms, monom))
    return out


def gl_n_lie_action(f_expr, syms):
    """
    Return list of polynomials [E_{ij} . f]_{i,j=1..n} where
    E_{ij} . f := x_i * d f / d x_j  (action of gl_n by linear
    substitutions on V*).  These span the tangent space T_f of the
    GL_n orbit.
    """
    n = len(syms)
    out = []
    for i in range(n):
        for j in range(n):
            df = sp.diff(f_expr, syms[j])
            out.append(sp.expand(syms[i] * df))
    return out


def linear_span_basis(polys, syms):
    """
    Return a basis (as list of polys) of the linear span of polys
    over Q in the polynomial ring.
    """
    if not polys:
        return []
    # Collect all monomial keys
    poly_dicts = [sp.Poly(sp.expand(p), *syms).as_dict() for p in polys]
    keys = set()
    for d in poly_dicts:
        keys.update(d.keys())
    keys = sorted(keys)
    # Build matrix with each polynomial as a column
    M = sp.zeros(len(keys), len(poly_dicts))
    for col, d in enumerate(poly_dicts):
        for row, k in enumerate(keys):
            M[row, col] = d.get(k, 0)
    # Reduce to find independent columns
    rref, pivots = M.rref()
    basis_polys = []
    for p_col in pivots:
        d = poly_dicts[p_col]
        poly = sp.Integer(0)
        for k, c in d.items():
            mono = sp.prod(s ** e for s, e in zip(syms, k))
            poly += c * mono
        basis_polys.append(sp.expand(poly))
    return basis_polys


def tangent_space_decomposition(f_expr, syms):
    """
    For each homogeneous degree d in f, compute the tangent space
    T^{(d)} = (gl_n . f) intersect Sym^d V*, and report:
      - tangent_dim
      - ambient_dim = dim Sym^d V*
      - is_full
      - weight_set: set of weights (a_1,...,a_n) of monomials reachable
        in the saturation gl_n . f_d (a monomial-level invariant).
      - missing_weights: weights of Sym^d V not reached by tangent.
    """
    out = {}
    comps = homogeneous_components(f_expr, syms)
    n = len(syms)
    for deg, fd in comps.items():
        if fd == 0:
            continue
        actions = gl_n_lie_action(fd, syms)
        basis = linear_span_basis(actions, syms)
        full_dim = sp.binomial(deg + n - 1, n - 1)
        # Compute torus weight occupation
        weights_seen = set()
        for p in basis:
            d_dict = sp.Poly(sp.expand(p), *syms).as_dict()
            for k, c in d_dict.items():
                if c != 0:
                    weights_seen.add(tuple(k))
        # Full set of weights of Sym^deg V
        full_weights = set(_weak_compositions(deg, n))
        missing = full_weights - weights_seen
        out[deg] = {
            "tangent_dim": len(basis),
            "ambient_dim": int(full_dim),
            "is_full": (len(basis) == int(full_dim)),
            "weights_seen": sorted(weights_seen),
            "n_missing_weights": len(missing),
            "missing_weights": sorted(missing),
        }
    return out


# ============================================================
# Polynomial constructors
# ============================================================

def f_chi_P(n):
    """f_chi_P^(n) = sum_{S subset [n], val(S) prime} prod_{i in S} x_i.
       1-indexed: val(S) = sum_{i in S} 2^(i-1)."""
    PRIMES = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,
              79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,
              157,163,167,173,179,181,191,193,197,199}
    syms = sp.symbols(f"x1:{n+1}")
    f = sp.Integer(0)
    for size in range(1, n + 1):
        for S in itertools.combinations(range(1, n + 1), size):
            val = sum(2 ** (i - 1) for i in S)
            if val in PRIMES:
                term = sp.Integer(1)
                for i in S:
                    term *= syms[i - 1]
                f += term
    return f, syms


def matched_random_polynomial(n, support, coef_range=(1, 7), seed=None):
    """Build poly with same support hypergraph, random coefficients."""
    rng = random.Random(seed)
    syms = sp.symbols(f"x1:{n+1}")
    f = sp.Integer(0)
    lo, hi = coef_range
    for S in support:
        c = rng.randint(lo, hi)
        term = sp.Integer(c)
        for i in S:
            term *= syms[i - 1]
        f += term
    return f, syms


def f_chi_P_support(n):
    """Return list of subsets S (as tuples) such that val(S) is prime."""
    PRIMES = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73}
    out = []
    for size in range(1, n + 1):
        for S in itertools.combinations(range(1, n + 1), size):
            val = sum(2 ** (i - 1) for i in S)
            if val in PRIMES:
                out.append(S)
    return out


# ============================================================
# Driver
# ============================================================

def driver(n=4, k_max=4, d_max=3, n_random=10, seed=42):
    print(f"\n=== A7 plethysm sub-frame, n = {n} ===")
    print()

    # ---- Plethysm table ----
    print(f"--- Sym^k(Sym^d C^{n}) decomposition into GL_{n} irreps S_lambda ---")
    plethysm_table = {}
    # Limit total degree to keep tractable; n^k * deg(d)^? complexity
    # Manually filter for tractability
    for d in range(1, d_max + 1):
        for k in range(1, k_max + 1):
            kd = k * d
            # Heuristic skip: skip cases beyond Schur-basis tractability.
            # At n=4 we can do up to total degree ~10 in <60s.
            # At n=5 limit total degree to ~6.
            if n <= 4 and kd > 10:
                continue
            if n == 5 and kd > 6:
                continue
            if n >= 6 and kd > 4:
                continue
            print(f"\n  Sym^{k}(Sym^{d} C^{n}):  total deg = {kd}")
            try:
                decomp = plethysm_sym_k_sym_d(n, d, k)
            except Exception as e:
                print(f"      [skipped: {e}]")
                continue
            for lam, mult in sorted(decomp.items()):
                print(f"      lambda={str(lam):<20s}  mult={mult}")
            plethysm_table[f"Sym^{k}(Sym^{d})"] = {
                "lam_to_mult": {str(lam): str(mult) for lam, mult in decomp.items()},
            }
    print()

    # ---- Tangent space of f_chi_P^(n) ----
    print(f"--- Tangent-space dim of GL_{n}-orbit of f_chi_P^({n}) ---")
    f, syms = f_chi_P(n)
    print(f"  f_chi_P^({n}) = {f}")
    tangent_chi = tangent_space_decomposition(f, syms)
    print(f"  Tangent-dim per homogeneous degree:")
    for deg, info in sorted(tangent_chi.items()):
        print(f"    deg {deg}: tangent_dim = {info['tangent_dim']}, "
              f"ambient_dim = {info['ambient_dim']}, "
              f"is_full = {info['is_full']}, "
              f"n_missing_weights = {info['n_missing_weights']}")
        if info["missing_weights"] and len(info["missing_weights"]) <= 8:
            print(f"      missing weights: {info['missing_weights']}")

    # ---- Tangent space of e_d, p_d ----
    print(f"\n--- Tangent-space dim of comparison polynomials at n={n} ---")
    syms_n = sp.symbols(f"x1:{n+1}")
    comparison = {}
    # e_d: elementary symmetric
    for d in range(1, n + 1):
        e_d = sp.Integer(0)
        for S in itertools.combinations(range(n), d):
            term = sp.Integer(1)
            for i in S:
                term *= syms_n[i]
            e_d += term
        td = tangent_space_decomposition(e_d, syms_n)
        comparison[f"e_{d}"] = td
        print(f"  e_{d}: {td}")
    # p_d: power sum
    for d in range(1, n + 1):
        p_d = sum(s ** d for s in syms_n)
        td = tangent_space_decomposition(p_d, syms_n)
        comparison[f"p_{d}"] = td
        print(f"  p_{d}: {td}")

    # ---- Random matched-support baselines ----
    print(f"\n--- Matched-support random baselines (n_random={n_random}) ---")
    support = f_chi_P_support(n)
    chi_p_dict = {f"deg{d}": info for d, info in tangent_chi.items()}
    baseline_results = []
    chi_p_weight_sets = {f"deg{d}": set(tuple(w) for w in info["weights_seen"])
                         for d, info in tangent_chi.items()}
    for trial in range(n_random):
        f_rand, syms_r = matched_random_polynomial(n, support, seed=seed + trial)
        td = tangent_space_decomposition(f_rand, syms_r)
        baseline_results.append({f"deg{d}": info for d, info in td.items()})
        if trial < 3 or trial == n_random - 1:
            print(f"  trial {trial}:  f = {f_rand}")
            for deg_lbl in sorted(td.keys()):
                info_b = td[deg_lbl]
                ws_chi = chi_p_weight_sets.get(f"deg{deg_lbl}", set())
                ws_b = set(tuple(w) for w in info_b["weights_seen"])
                same = (ws_chi == ws_b)
                print(f"    deg {deg_lbl}: dim={info_b['tangent_dim']}, "
                      f"weight_set==chi_P? {same}")
    # Aggregate per-degree comparison
    print(f"\n--- Aggregate comparison (z-score per degree) ---")
    all_degrees = set()
    for tr in baseline_results:
        for k in tr:
            all_degrees.add(k)
    summary = {}
    for d_label in sorted(all_degrees):
        chi_dim = chi_p_dict.get(d_label, {}).get("tangent_dim", None)
        baseline_dims = [tr[d_label]["tangent_dim"]
                         for tr in baseline_results if d_label in tr]
        if not baseline_dims or chi_dim is None:
            continue
        m = sum(baseline_dims) / len(baseline_dims)
        var = sum((x - m) ** 2 for x in baseline_dims) / len(baseline_dims)
        sd = var ** 0.5
        z = ((chi_dim - m) / sd) if sd > 0 else 0
        ambient = chi_p_dict[d_label]["ambient_dim"]
        is_full = chi_p_dict[d_label]["is_full"]
        summary[d_label] = {
            "chi_p_tangent_dim": chi_dim,
            "ambient_dim": ambient,
            "chi_p_is_full": is_full,
            "baseline_mean": m,
            "baseline_std": sd,
            "z_score": float(z),
            "baseline_dims": baseline_dims,
        }
        print(f"  {d_label}: chi_P_dim={chi_dim}, ambient={ambient}, "
              f"is_full={is_full}, baseline mean={m:.2f}, std={sd:.3f}, "
              f"z={z:.3f}")
    return {
        "n": n,
        "plethysm_table": plethysm_table,
        "f_chi_P_tangent": chi_p_dict,
        "comparison": {k: {f"deg{d}": v for d, v in info.items()}
                       for k, info in comparison.items()},
        "baseline_summary": summary,
    }


if __name__ == "__main__":
    n_arg = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    k_max = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    d_max = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    n_rand = int(sys.argv[4]) if len(sys.argv) > 4 else 10
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 42

    out = driver(n=n_arg, k_max=k_max, d_max=d_max,
                 n_random=n_rand, seed=seed)

    out_path = os.path.join(
        os.path.dirname(__file__) or ".",
        f"plethysm_subframe_n{n_arg}_results.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nWrote {out_path}")
