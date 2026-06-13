"""
A7 commit thread, S212 (slot 2 of 5): second-order plethysm sub-frame.

This script computes I_2(orbit-closure(f)) for several polynomials of
degree d=3 in n=5 variables, decomposes I_2 (and the quotient
Sym^2(Sym^3 V_5) / I_2) into GL_5 irreps S_lambda for partitions
lambda |= 6, and compares results between:

  - tilde f_chi_P^(4)        : degree-3 homogenisation of f_chi_P^(4)
                                (n=4 chi_P padded by y = x_5)
  - det_2 * y                : deg-2 determinant in 4 vars padded by y
  - perm_2 * y               : deg-2 permanent in 4 vars padded by y
  - matched-support random   : K trials with same monomial support as
                                tilde f_chi_P^(4)
  - x_1^3                    : control (cube; orbit closure = veronese)

The "occurrence obstruction" sub-frame asks: which S_lambda occur in
Sym^2(Sym^3 V_5) / I_2(orbit-closure(f)) for one polynomial but NOT for
another (with same support hypergraph)?  An S_lambda found only in
Sym^2/I_2(f_chi_P) but not in Sym^2/I_2(matched-baseline) would be a
chi_P-specific GL_5-invariant -- the first non-support-determined
plethysm-level invariant.

Method
------
1. Compute the GL_n-action matrix R_3(g) on Sym^3 V_n for random g ∈ GL_n.
2. For each polynomial f, generate M random orbit points w_t = R_3(g_t) f
   (each w_t is a length-N coefficient vector in the monomial basis of
   Sym^3 V_n; here N = 35 for n=5, d=3).
3. Build the M x N(N+1)/2 evaluation matrix E whose row t has entries
   (w_t)_α (w_t)_β for α <= β (a basis of Sym^2(Sym^3 V*) under
   the natural identification of Sym^d V with its dual).
4. ker(E) (right null space) = I_2(orbit-closure(f)) ⊆ Sym^2(Sym^3 V).
5. Each basis pair (alpha, beta) carries the diagonal-torus weight
   alpha + beta. So Sym^2(Sym^3 V) splits into weight-graded components,
   and I_2 inherits this grading. Compute dim(I_2)_gamma weight-by-weight.
6. Read off the GL_n-character of I_2 (and of the quotient
   Sym^2(Sym^3 V) / I_2) and decompose into Schur basis via
   Vandermonde inversion.
7. Compare the resulting S_lambda multiplicities across the test
   polynomials. Any difference is an occurrence-obstruction signal.

Cross-domain refs:
  Mulmuley-Sohoni 2001 SIAM J Comput 31, 496;
  Burgisser-Ikenmeyer-Panova 2017 FOCS arXiv:1604.06431;
  Landsberg 2017 *Geometry and Complexity Theory*, ch. 9-10
    (esp. the catalecticant identification of I_2 for cubes).
"""

import itertools
import json
import os
import random
import sys
import time
from fractions import Fraction
from functools import lru_cache
from math import factorial

import numpy as np
import sympy as sp


# ============================================================
# Partition utilities (re-imported from S211)
# ============================================================

def partitions_of(N, max_parts=None, max_part=None):
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


def weak_compositions(d, n):
    if n == 1:
        yield (d,)
        return
    for k in range(d + 1):
        for rest in weak_compositions(d - k, n - 1):
            yield (k,) + rest


def multinomial(parts):
    """Multinomial of |parts|! / prod(parts[i]!)."""
    n_total = sum(parts)
    out = factorial(n_total)
    for p in parts:
        out //= factorial(p)
    return out


# ============================================================
# Schur via bialternant (small support)
# ============================================================

@lru_cache(maxsize=None)
def schur_polynomial_in_n_vars(lam, n, vars_tuple):
    if len(lam) > n:
        return sp.Integer(0)
    L = list(lam) + [0] * (n - len(lam))
    syms = vars_tuple
    num_mat = sp.Matrix(n, n, lambda i, j: syms[i] ** (L[j] + n - 1 - j))
    den_mat = sp.Matrix(n, n, lambda i, j: syms[i] ** (n - 1 - j))
    num = num_mat.det()
    den = den_mat.det()
    return sp.simplify(sp.cancel(num / den))


def expand_in_schur_basis(P, syms, total_degree, max_parts):
    """Decompose symmetric polynomial P (homogeneous of total_degree) into
    sum of Schur polys s_lambda.  Returns dict {lambda: integer multiplicity}."""
    n = len(syms)
    parts_list = [lam for lam in partitions_of(total_degree, max_parts=n)]
    monomial_dict_target = sp.Poly(P, *syms).as_dict()
    all_keys = set(monomial_dict_target.keys())
    schur_expansions = []
    for lam in parts_list:
        s_lam = schur_polynomial_in_n_vars(lam, n, tuple(syms))
        sd = sp.Poly(sp.expand(s_lam), *syms).as_dict() if s_lam != 0 else {}
        all_keys.update(sd.keys())
        schur_expansions.append((lam, sd))
    keys = sorted(all_keys)
    A = sp.zeros(len(keys), len(parts_list))
    b = sp.zeros(len(keys), 1)
    for i, k in enumerate(keys):
        b[i, 0] = monomial_dict_target.get(k, 0)
        for j, (lam, sd) in enumerate(schur_expansions):
            A[i, j] = sd.get(k, 0)
    sol = A.solve(b)
    out = {}
    for j, (lam, _) in enumerate(schur_expansions):
        v = sp.nsimplify(sol[j, 0])
        if v != 0:
            out[lam] = int(v)
    return out


def schur_dim_in_n_vars(lam, n):
    """Dim of irrep S_lambda of GL_n.  Uses Weyl dim formula."""
    if len(lam) > n:
        return 0
    L = list(lam) + [0] * (n - len(lam))
    num = 1
    den = 1
    for i in range(n):
        for j in range(i + 1, n):
            num *= (L[i] - L[j] + j - i)
            den *= (j - i)
    return num // den


# ============================================================
# Sym^d action by GL_n
# ============================================================

def build_action_structure(n, d, basis_d, basis_dict):
    """For each (beta_idx, alpha_idx), enumerate matrices B in N^{n,n}
    with B.sum(axis=1) = alpha and B.sum(axis=0) = beta, and store
    multinomial coefficients.  Each entry (beta_idx, alpha_idx) gets a list
    of (coef, B) pairs.

    Then R_d(g)[beta, alpha] = sum (coef * prod_{i,j} g[i,j]**B[i,j]).
    """
    structure = {}
    for alpha_idx, alpha in enumerate(basis_d):
        # Each row B[i] is a composition of alpha[i] into n parts
        rows_per_i = [list(weak_compositions(int(alpha[i]), n)) for i in range(n)]
        for B_rows in itertools.product(*rows_per_i):
            B = np.array(B_rows, dtype=int)
            beta = tuple(B.sum(axis=0).tolist())
            if beta not in basis_dict:
                continue
            beta_idx = basis_dict[beta]
            # Multinomial coefficient: prod_i (alpha_i! / prod_j B[i,j]!)
            coef = 1
            for i in range(n):
                row_vals = [int(B[i, j]) for j in range(n)]
                coef *= multinomial(row_vals)
            entry = structure.setdefault((beta_idx, alpha_idx), [])
            entry.append((coef, B))
    return structure


def action_matrix(g, structure, N, n):
    """Compute R_d(g) given precomputed structure."""
    R = np.zeros((N, N), dtype=np.float64)
    for (beta_idx, alpha_idx), terms in structure.items():
        s = 0.0
        for coef, B in terms:
            prod = float(coef)
            for i in range(n):
                for j in range(n):
                    e = int(B[i, j])
                    if e > 0:
                        prod *= float(g[i, j]) ** e
            s += prod
        R[beta_idx, alpha_idx] = s
    return R


# ============================================================
# Polynomials in Sym^d V_n (coefficient vectors)
# ============================================================

def poly_to_coef_vec(f_expr, basis_d, syms):
    """Convert sympy poly to coef vector in basis_d."""
    poly = sp.Poly(sp.expand(f_expr), *syms)
    n_terms = len(basis_d)
    vec = np.zeros(n_terms, dtype=np.float64)
    pd = poly.as_dict()
    for k, c in pd.items():
        if k in [tuple(a) for a in basis_d]:
            pass  # convert to numpy
    # Actually do it properly
    bdict = {tuple(a): i for i, a in enumerate(basis_d)}
    for k, c in pd.items():
        kt = tuple(k)
        if kt in bdict:
            vec[bdict[kt]] = float(c)
        else:
            raise ValueError(f"Monomial {k} not in deg-{sum(k)} basis")
    return vec


def f_chi_P_homogenised(n_inner, d_target, syms_full):
    """f_chi_P^(n_inner) homogenised to degree d_target using the n_inner+1-th
    variable as y.  Returns sympy expression in syms_full = (x_1, ..., x_{n_inner}, y).
    Note: requires d_target >= max degree in f_chi_P, which is n_inner."""
    PRIMES = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73}
    n = n_inner
    y = syms_full[n]  # last var is y
    f = sp.Integer(0)
    for size in range(1, n + 1):
        for S in itertools.combinations(range(1, n + 1), size):
            val = sum(2 ** (i - 1) for i in S)
            if val in PRIMES:
                term = sp.Integer(1)
                for i in S:
                    term *= syms_full[i - 1]
                deg = size
                pad = d_target - deg
                if pad < 0:
                    raise ValueError(f"d_target={d_target} too low for f_chi_P^({n})")
                term *= y ** pad
                f += term
    return f


def f_chi_P_support(n_inner):
    """List of subsets S of [n_inner] with val(S) prime."""
    PRIMES = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73}
    out = []
    for size in range(1, n_inner + 1):
        for S in itertools.combinations(range(1, n_inner + 1), size):
            val = sum(2 ** (i - 1) for i in S)
            if val in PRIMES:
                out.append((S, size))
    return out


def matched_random_homogenised(n_inner, d_target, syms_full, coef_range=(1, 7), seed=42):
    rng = random.Random(seed)
    n = n_inner
    y = syms_full[n]
    support = f_chi_P_support(n_inner)
    f = sp.Integer(0)
    lo, hi = coef_range
    for S, deg in support:
        c = rng.randint(lo, hi)
        term = sp.Integer(c)
        for i in S:
            term *= syms_full[i - 1]
        pad = d_target - deg
        term *= y ** pad
        f += term
    return f


def det_2_padded(syms_full, d_target):
    """det_2 = x_1*x_4 - x_2*x_3, then padded by y to degree d_target."""
    if len(syms_full) < 5:
        raise ValueError("Need >= 5 vars for det_2_padded")
    x = syms_full
    return (x[0] * x[3] - x[1] * x[2]) * (x[4] ** (d_target - 2))


def perm_2_padded(syms_full, d_target):
    """perm_2 = x_1*x_4 + x_2*x_3, padded by y."""
    if len(syms_full) < 5:
        raise ValueError("Need >= 5 vars for perm_2_padded")
    x = syms_full
    return (x[0] * x[3] + x[1] * x[2]) * (x[4] ** (d_target - 2))


def x1_cubed(syms_full):
    return syms_full[0] ** 3


# ============================================================
# Main I_2 computation
# ============================================================

def generate_orbit_samples(f_vec, structure, N, n, M, seed=None, coef_max=2):
    """Generate M random orbit points by applying R_d(g) to f_vec for
    random g ∈ GL_n(Z).  Returns (M, N) numpy array of coef vectors."""
    rng = np.random.default_rng(seed)
    orbit_coefs = np.zeros((M, N), dtype=np.float64)
    n_full_rank = 0
    attempts = 0
    while n_full_rank < M and attempts < 10 * M:
        g = rng.integers(-coef_max, coef_max + 1, size=(n, n)).astype(float)
        if abs(np.linalg.det(g)) < 0.5:
            attempts += 1
            continue
        R = action_matrix(g, structure, N, n)
        w = R @ f_vec
        orbit_coefs[n_full_rank] = w
        n_full_rank += 1
        attempts += 1
    if n_full_rank < M:
        raise RuntimeError(f"Could not generate {M} samples; got {n_full_rank}")
    return orbit_coefs


def compute_I_2_weight_dims(f_vec, structure, N, n, M, basis_d, seed=None,
                             coef_max=2):
    """For each weight gamma ∈ N^n with |gamma| = 2d, compute
       dim (I_2)_gamma = (size of weight-gamma block) - rank(E_gamma)
    where E_gamma[t, k] = w_t[alpha_k] * w_t[beta_k] for orbit samples w_t
    and basis pairs (alpha_k, beta_k) of weight gamma.

    The kernel-per-weight computation is exact in the equivariant setting:
    since I_2 ⊆ Sym^2(Sym^d V) is GL_n-stable (hence T-stable), it splits
    weight-by-weight.  The constraint "P(w) = 0 for all w ∈ orbit-closure"
    decomposes into one constraint per weight component of P (by torus
    equivariance applied to the orbit-closure ideal).

    Returns: (weight_dims_I2 dict, weight_dims_Sym2 dict, dim_I2_total,
              orbit_coefs).
    """
    orbit_coefs = generate_orbit_samples(f_vec, structure, N, n, M,
                                          seed=seed, coef_max=coef_max)
    # Group pairs by weight
    pairs = [(i, j) for i in range(N) for j in range(i, N)]
    weight_to_pairs = {}
    for k, (i, j) in enumerate(pairs):
        w = tuple(int(a + b) for a, b in zip(basis_d[i], basis_d[j]))
        weight_to_pairs.setdefault(w, []).append((k, i, j))
    weight_dims_I2 = {}
    weight_dims_Sym2 = {}
    for w, plist in weight_to_pairs.items():
        block_size = len(plist)
        weight_dims_Sym2[w] = block_size
        # Build E_w (M x block_size)
        E_w = np.zeros((M, block_size), dtype=np.float64)
        for col, (_, i, j) in enumerate(plist):
            E_w[:, col] = orbit_coefs[:, i] * orbit_coefs[:, j]
        # Rank via SVD
        if block_size == 1:
            # rank is 0 if all entries are zero, else 1
            rank = 1 if np.max(np.abs(E_w)) > 1e-10 else 0
        else:
            s = np.linalg.svd(E_w, compute_uv=False)
            scale = s[0] if len(s) > 0 else 1.0
            tol = max(M, block_size) * np.finfo(float).eps * max(scale, 1.0)
            rank = int(np.sum(s > tol))
        weight_dims_I2[w] = block_size - rank
    dim_I2_total = sum(weight_dims_I2.values())
    return weight_dims_I2, weight_dims_Sym2, dim_I2_total, orbit_coefs


# Backward-compat: kept name compute_I_2_kernel for callers but reimplemented
def compute_I_2_kernel(f_vec, structure, N, n, M, seed=None,
                       coef_max=2, return_E=False):
    weight_I2, weight_full, dim_I2, orbit_coefs = compute_I_2_weight_dims(
        f_vec, structure, N, n, M, BASIS_D_REF[0], seed=seed, coef_max=coef_max)
    pairs = [(i, j) for i in range(N) for j in range(i, N)]
    return weight_I2, weight_full, dim_I2, pairs


# Reference holder so backward compat helper finds the basis
BASIS_D_REF = [None]


def weight_dims_of_sym_2_sym_d(basis_d, n):
    """For comparison: dim of Sym^2(Sym^d V) graded by weight."""
    out = {}
    for i, alpha in enumerate(basis_d):
        for j, beta in enumerate(basis_d):
            if j < i:
                continue
            w = tuple(int(a + b) for a, b in zip(alpha, beta))
            out[w] = out.get(w, 0) + 1
    return out


def schur_decompose_from_weights(weight_dim_dict, n, total_deg):
    """Given weight dims (dict {tuple -> int}), build the character
    chi(t_1, ..., t_n) = sum_w dim_w * t^w, then expand in Schur basis.
    Returns: dict {lambda: multiplicity}."""
    syms = sp.symbols(f'x1:{n+1}')
    char = sp.Integer(0)
    for w, dim in weight_dim_dict.items():
        if dim == 0:
            continue
        mono = sp.Integer(1)
        for i, e in enumerate(w):
            mono *= syms[i] ** e
        char += dim * mono
    char = sp.expand(char)
    # Schur basis
    return expand_in_schur_basis(char, syms, total_deg, max_parts=n)


# ============================================================
# Driver
# ============================================================

def run_for_polynomial(name, f_expr, syms_full, basis_d, basis_dict,
                       structure, N, n, M, seed):
    """Run the I_2 computation for one polynomial. Returns summary dict."""
    print(f"\n=== {name} ===")
    print(f"  f = {f_expr}")
    d = sum(basis_d[0])  # degree of monomials in basis_d
    f_vec = poly_to_coef_vec(f_expr, basis_d, syms_full)
    weight_dims_I2, weight_dims_full, dim_I2, _ = compute_I_2_weight_dims(
        f_vec, structure, N, n, M, basis_d, seed=seed)
    N2 = N * (N + 1) // 2
    print(f"  dim Sym^2(Sym^{d} V_{n}) = {N2}")
    print(f"  dim I_2(orbit-closure) = {dim_I2}  (per-weight kernel split)")
    print(f"  dim Sym^2 / I_2       = {N2 - dim_I2}")
    weight_dims_quot = {w: weight_dims_full[w] - weight_dims_I2.get(w, 0)
                        for w in weight_dims_full}
    # Schur decomposition (of I_2 and of quotient)
    print("  Computing Schur decomposition of I_2 and of Sym^2 / I_2 ...")
    decomp_I2 = schur_decompose_from_weights(weight_dims_I2, n, total_deg=2 * d)
    decomp_quot = schur_decompose_from_weights(weight_dims_quot, n, total_deg=2 * d)
    decomp_full = schur_decompose_from_weights(weight_dims_full, n, total_deg=2 * d)
    # Sanity: sum should equal Sym^2(Sym^d) decomp
    print("  S_lambda multiplicities (Sym^2(Sym^3) | I_2 | quotient):")
    all_lams = sorted(set(decomp_full.keys()) | set(decomp_I2.keys()) | set(decomp_quot.keys()))
    for lam in all_lams:
        a = decomp_full.get(lam, 0)
        b = decomp_I2.get(lam, 0)
        c = decomp_quot.get(lam, 0)
        marker = " *" if (b == 0 and c > 0) else ("  " if a == b + c else " !!")
        print(f"    {str(lam):<14s}  full={a}  I2={b}  quot={c}{marker}")
    return {
        "name": name,
        "dim_Sym2": N2,
        "dim_I2": dim_I2,
        "dim_quot": N2 - dim_I2,
        "M_samples": M,
        "weight_dims_I2": {str(w): v for w, v in weight_dims_I2.items()},
        "weight_dims_quot": {str(w): v for w, v in weight_dims_quot.items()},
        "schur_full": {str(lam): v for lam, v in decomp_full.items()},
        "schur_I2": {str(lam): v for lam, v in decomp_I2.items()},
        "schur_quot": {str(lam): v for lam, v in decomp_quot.items()},
    }


def driver(n, d, M, n_baseline, seed):
    print(f"=== A7 commit thread, S212: I_2(orbit-closure) sub-frame, n={n}, d={d} ===\n")

    syms_full = sp.symbols(f"x1:{n+1}")
    basis_d = list(weak_compositions(d, n))
    basis_dict = {tuple(a): i for i, a in enumerate(basis_d)}
    N = len(basis_d)
    print(f"dim Sym^{d} V_{n} = {N}")
    N2 = N * (N + 1) // 2
    print(f"dim Sym^2(Sym^{d} V_{n}) = {N2}")
    print(f"M (random orbit samples) = {M}")
    if M < N2 + 50:
        print(f"  WARNING: M = {M} may be too small (recommend M >= {N2 + 50})")

    # Build action structure
    print(f"\nBuilding GL_{n}-action structure on Sym^{d} V_{n} ...")
    t0 = time.time()
    structure = build_action_structure(n, d, basis_d, basis_dict)
    print(f"  structure has {len(structure)} non-zero (beta, alpha) pairs,")
    n_terms_total = sum(len(v) for v in structure.values())
    print(f"  total {n_terms_total} (coef, B) entries; build time {time.time()-t0:.1f}s")

    # Sym^2(Sym^d) full character (for reference)
    weight_dims_full = weight_dims_of_sym_2_sym_d(basis_d, n)
    print(f"\nSym^2(Sym^{d} V_{n}) weight count: {len(weight_dims_full)} weights, "
          f"sum = {sum(weight_dims_full.values())}")

    # Run for each polynomial
    results = {}

    if n >= 5 and d >= 3:
        # 1. tilde f_chi_P^(4) homogenised (in 5+ vars).
        f1 = f_chi_P_homogenised(4, d, syms_full)
        results["chi_P_n4"] = run_for_polynomial(
            f"tilde f_chi_P^(4) (n=4 chi_P padded by y to deg {d})",
            f1, syms_full, basis_d, basis_dict, structure, N, n, M, seed)
        # 2. det_2 * y^{d-2}
        f2 = det_2_padded(syms_full, d)
        results["det_2_y"] = run_for_polynomial(
            f"det_2 * y^{d-2}",
            f2, syms_full, basis_d, basis_dict, structure, N, n, M, seed + 1)
        # 3. perm_2 * y^{d-2}
        f3 = perm_2_padded(syms_full, d)
        results["perm_2_y"] = run_for_polynomial(
            f"perm_2 * y^{d-2}",
            f3, syms_full, basis_d, basis_dict, structure, N, n, M, seed + 2)
        # 5-7. Matched-support baselines for chi_P
        for k in range(n_baseline):
            f_b = matched_random_homogenised(4, d, syms_full, seed=seed + 10 + k)
            results[f"baseline_{k}"] = run_for_polynomial(
                f"matched-support baseline #{k} (seed {seed + 10 + k})",
                f_b, syms_full, basis_d, basis_dict, structure, N, n, M,
                seed + 100 + k)
    elif n == 4 and d == 3:
        # n=4, d=3: use deg-3 component of f_chi_P^(4) directly (no padding)
        chi_P_d3 = (syms_full[0] * syms_full[1] * syms_full[2]
                    + syms_full[0] * syms_full[1] * syms_full[3]
                    + syms_full[0] * syms_full[2] * syms_full[3])
        results["chi_P_n4_d3"] = run_for_polynomial(
            "f_chi_P^(4) deg-3 component (= x1 x2 x3 + x1 x2 x4 + x1 x3 x4)",
            chi_P_d3, syms_full, basis_d, basis_dict, structure, N, n, M, seed)
        # e_3(x_1..x_4) = elementary symmetric poly deg 3 in 4 vars
        e_3 = sum(sp.prod(syms_full[i] for i in S)
                  for S in itertools.combinations(range(4), 3))
        results["e_3"] = run_for_polynomial(
            "e_3 (elementary symmetric deg 3 in 4 vars)",
            e_3, syms_full, basis_d, basis_dict, structure, N, n, M, seed + 1)
        # p_3(x_1..x_4) = power sum deg 3
        p_3 = sum(s ** 3 for s in syms_full[:4])
        results["p_3"] = run_for_polynomial(
            "p_3 (power sum deg 3 in 4 vars)",
            p_3, syms_full, basis_d, basis_dict, structure, N, n, M, seed + 2)
        # x_1 x_2 x_3 (single Plücker monomial)
        plucker = syms_full[0] * syms_full[1] * syms_full[2]
        results["plucker"] = run_for_polynomial(
            "x_1 x_2 x_3 (single triple monomial)",
            plucker, syms_full, basis_d, basis_dict, structure, N, n, M, seed + 3)
        # baselines: random integer coefs on chi_P_d3 support
        chi_P_support_d3 = [(0,1,2), (0,1,3), (0,2,3)]
        rng = random.Random(seed + 5)
        for k in range(n_baseline):
            local_rng = random.Random(seed + 10 + k)
            f_b = sum(local_rng.randint(1, 7) * sp.prod(syms_full[i] for i in S)
                      for S in chi_P_support_d3)
            results[f"baseline_{k}"] = run_for_polynomial(
                f"matched-support baseline #{k} (seed {seed + 10 + k})",
                f_b, syms_full, basis_d, basis_dict, structure, N, n, M,
                seed + 100 + k)
    else:
        raise ValueError(f"Unsupported (n,d) = ({n},{d})")

    # 4. x_1^3 (control: cubical orbit closure = Veronese)
    f4 = x1_cubed(syms_full)
    results["x1_cubed"] = run_for_polynomial(
        "x_1^3 (control: cube)",
        f4, syms_full, basis_d, basis_dict, structure, N, n, M, seed + 3)

    # ---- Cross-comparison summary ----
    print("\n\n=========== CROSS-COMPARISON ===========")
    print("\n--- dim I_2 per polynomial ---")
    for name, r in results.items():
        print(f"  {name:35s}  dim_I2 = {r['dim_I2']:3d}, "
              f"dim_quot = {r['dim_quot']:3d}")

    # Schur multiplicity comparison
    all_keys = set()
    for r in results.values():
        all_keys.update(r['schur_quot'].keys())
    all_keys_sorted = sorted(all_keys)
    # Find chi_P-style key
    chi_P_key = next((k for k in results if k.startswith("chi_P")), None)
    baseline_keys = [k for k in results if k.startswith("baseline_")]
    other_keys = [k for k in results
                  if not (k.startswith("chi_P") or k.startswith("baseline_"))]

    print("\n--- Schur multiplicity in Sym^2 / I_2 (quotient) ---")
    header_cols = []
    if chi_P_key:
        header_cols.append(("chi_P", chi_P_key))
    for k in other_keys:
        header_cols.append((k[:8], k))
    for k in baseline_keys:
        header_cols.append((k.replace("baseline_", "b"), k))
    header = "  " + f"{'lambda':<16s}" + " | ".join(f"{h:>6s}" for h, _ in header_cols)
    print(header)
    print("  " + "-" * len(header))
    occurrence_signals = []
    for lam_str in all_keys_sorted:
        row_vals = [results[k]['schur_quot'].get(lam_str, 0)
                    for _, k in header_cols]
        chi_P_val = results[chi_P_key]['schur_quot'].get(lam_str, 0) if chi_P_key else 0
        baseline_vals = [results[k]['schur_quot'].get(lam_str, 0)
                         for k in baseline_keys]
        flag = ""
        if chi_P_key and baseline_vals:
            if chi_P_val > 0 and all(b == 0 for b in baseline_vals):
                flag += " [chi_P > 0, all_baseline = 0]"
            if chi_P_val == 0 and any(b > 0 for b in baseline_vals):
                flag += " [chi_P = 0, baseline > 0]"
        if flag:
            occurrence_signals.append((lam_str, flag, dict(zip([h for h, _ in header_cols], row_vals))))
        print(f"  {lam_str:<16s}" + " | ".join(f"{v:6d}" for v in row_vals) + flag)

    print(f"\n--- Occurrence-obstruction candidates (S_λ in I_2 differs between chi_P and matched baseline) ---")
    if occurrence_signals:
        print(f"  Found {len(occurrence_signals)} candidates:")
        for lam_str, flag, vals in occurrence_signals:
            print(f"    lambda={lam_str}: {vals}{flag}")
    else:
        print(f"  NONE found across {n_baseline} baselines + reference polynomials")

    return {
        "n": n, "d": d, "M": M, "n_baseline": n_baseline,
        "dim_Sym_d_V": N, "dim_Sym_2_Sym_d_V": N2,
        "results": results,
        "occurrence_signals": [
            {"lambda": s[0], "flag": s[1], "chi_P": s[2], "det_2_y": s[3],
             "perm_2_y": s[4], "x_1_cubed": s[5], "baselines": s[6]}
            for s in occurrence_signals
        ],
    }


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    d = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    M = int(sys.argv[3]) if len(sys.argv) > 3 else 800
    n_baseline = int(sys.argv[4]) if len(sys.argv) > 4 else 3
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 4242

    out = driver(n=n, d=d, M=M, n_baseline=n_baseline, seed=seed)
    out_path = os.path.join(
        os.path.dirname(__file__) or ".",
        f"plethysm_subframe_I2_n{n}_d{d}_results.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\n[done] Wrote {out_path}")
