"""
A7 commit thread, S213 (slot 3 of 5): third-order plethysm sub-frame.

This script computes I_3(orbit-closure(f)) for several polynomials of
degree d=3 in n=4 variables, decomposes I_3 (and the quotient
Sym^3(Sym^3 V_4) / I_3) into GL_4 irreps S_lambda for partitions
lambda |= 9 (at most 4 parts), and compares results between:

  - chi_P_n4_d3 : f_chi_P^(4) deg-3 component = x_1 x_2 x_3 + x_1 x_2 x_4
                  + x_1 x_3 x_4
  - e_3         : elementary symmetric deg 3 in 4 vars
  - p_3         : power sum deg 3 in 4 vars
  - x_1 x_2 x_3 : single triple monomial (Plücker)
  - matched-support random : K trials with same monomial support as
                              chi_P_n4_d3 ((1,2,3),(1,2,4),(1,3,4))
  - x_1^3       : control (cube; orbit closure = cubic Veronese ν_3(P^3))

The "occurrence obstruction" sub-frame at level k=3 asks: which S_lambda
occur in Sym^3(Sym^3 V_4) / I_3(orbit-closure(f)) for one polynomial but
NOT for another (with same support hypergraph)?  An S_lambda found only
in Sym^3/I_3(chi_P) but not in Sym^3/I_3(matched-baseline) would be a
chi_P-specific GL_4-invariant -- the first non-support-determined
plethysm-level invariant.

Per S211 plethysm machinery:
  Sym^3(Sym^3 V_4) = S_(9) + S_(7,2) + S_(6,3) + S_(5,2,2) + S_(4,4,1)
  with Weyl dims                220 + 540 + 480 + 160      + 140 = 1540.

Method (parallel to S212 I_2 computation):
  1. Compute the GL_n-action matrix R_3(g) on Sym^3 V_n for random g ∈ GL_n.
  2. For each f, generate M random orbit points w_t = R_3(g_t) f.
  3. For each torus weight gamma ∈ N^n with |gamma| = 3d = 9, collect
     the basis triples (alpha, beta, delta) of monomials in Sym^3 V_n
     with alpha + beta + delta = gamma (a basis of weight-gamma block
     of Sym^3(Sym^3 V*)).  Build M x |block| evaluation matrix E_gamma
     whose row t records (w_t)_alpha * (w_t)_beta * (w_t)_delta.
  4. ker(E_gamma) = (I_3)_gamma; sum over weights to get full I_3 + char.
  5. Schur-decompose via S211 Vandermonde-inversion machinery.

Cross-domain refs:
  Mulmuley-Sohoni 2001 SIAM J Comput 31, 496;
  Burgisser-Ikenmeyer-Panova 2017 FOCS arXiv:1604.06431;
  Landsberg 2017 *Geometry and Complexity Theory*, ch. 9-10
    (catalecticant identification of I_2 for cubes);
  Iarrobino-Kanev 1999 *Power Sums, Gorenstein Algebras, and Determinantal
    Loci*, LNM 1721 (catalecticant Hilbert function for cubic Veronese
    in higher k).
"""

import itertools
import json
import os
import random
import sys
import time
from functools import lru_cache
from math import factorial

import numpy as np
import sympy as sp


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
    n_total = sum(parts)
    out = factorial(n_total)
    for p in parts:
        out //= factorial(p)
    return out


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
    n = len(syms)
    parts_list = [lam for lam in partitions_of(total_degree, max_parts=n)]
    monomial_dict_target = sp.Poly(P, *syms).as_dict() if P != 0 else {}
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


def build_action_structure(n, d, basis_d, basis_dict):
    structure = {}
    for alpha_idx, alpha in enumerate(basis_d):
        rows_per_i = [list(weak_compositions(int(alpha[i]), n)) for i in range(n)]
        for B_rows in itertools.product(*rows_per_i):
            B = np.array(B_rows, dtype=int)
            beta = tuple(B.sum(axis=0).tolist())
            if beta not in basis_dict:
                continue
            beta_idx = basis_dict[beta]
            coef = 1
            for i in range(n):
                row_vals = [int(B[i, j]) for j in range(n)]
                coef *= multinomial(row_vals)
            entry = structure.setdefault((beta_idx, alpha_idx), [])
            entry.append((coef, B))
    return structure


def action_matrix(g, structure, N, n):
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


def poly_to_coef_vec(f_expr, basis_d, syms):
    poly = sp.Poly(sp.expand(f_expr), *syms)
    n_terms = len(basis_d)
    vec = np.zeros(n_terms, dtype=np.float64)
    pd = poly.as_dict()
    bdict = {tuple(a): i for i, a in enumerate(basis_d)}
    for k, c in pd.items():
        kt = tuple(k)
        if kt in bdict:
            vec[bdict[kt]] = float(c)
        else:
            raise ValueError(f"Monomial {k} not in deg-{sum(k)} basis")
    return vec


def x1_cubed(syms_full):
    return syms_full[0] ** 3


def chi_P_d3_n4(syms_full):
    """f_chi_P^(4) deg-3 component."""
    return (syms_full[0] * syms_full[1] * syms_full[2]
            + syms_full[0] * syms_full[1] * syms_full[3]
            + syms_full[0] * syms_full[2] * syms_full[3])


def e_3_n4(syms_full):
    return sum(sp.prod(syms_full[i] for i in S)
               for S in itertools.combinations(range(4), 3))


def p_3_n4(syms_full):
    return sum(s ** 3 for s in syms_full[:4])


def plucker_n4(syms_full):
    return syms_full[0] * syms_full[1] * syms_full[2]


def matched_baseline(syms_full, support_d3, seed):
    rng = random.Random(seed)
    return sum(rng.randint(1, 7) * sp.prod(syms_full[i] for i in S)
               for S in support_d3)


def generate_orbit_samples(f_vec, structure, N, n, M, seed=None, coef_max=2):
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


def weight_dims_of_sym_3_sym_d(basis_d, n):
    """Dim of Sym^3(Sym^d V) graded by torus weight."""
    out = {}
    for triple in itertools.combinations_with_replacement(range(len(basis_d)), 3):
        i, j, k = triple
        w = tuple(int(basis_d[i][p] + basis_d[j][p] + basis_d[k][p])
                  for p in range(n))
        out[w] = out.get(w, 0) + 1
    return out


def compute_I_3_weight_dims(f_vec, structure, N, n, M, basis_d, seed=None,
                             coef_max=2, verbose=False):
    """For each weight gamma ∈ N^n with |gamma| = 3d, compute
       dim (I_3)_gamma = (size of weight-gamma block) - rank(E_gamma).

    Returns: (weight_dims_I3, weight_dims_full, dim_I3_total, orbit_coefs).
    """
    t_orbit = time.time()
    orbit_coefs = generate_orbit_samples(f_vec, structure, N, n, M,
                                          seed=seed, coef_max=coef_max)
    if verbose:
        print(f"    [orbit gen: {time.time()-t_orbit:.1f}s]")

    # Group triples by weight
    t_group = time.time()
    weight_to_triples = {}
    for i in range(N):
        for j in range(i, N):
            for k in range(j, N):
                w = tuple(int(basis_d[i][p] + basis_d[j][p] + basis_d[k][p])
                          for p in range(n))
                weight_to_triples.setdefault(w, []).append((i, j, k))
    if verbose:
        print(f"    [weight grouping: {time.time()-t_group:.1f}s, "
              f"n_weights = {len(weight_to_triples)}]")

    weight_dims_I3 = {}
    weight_dims_full = {}
    t_svd = time.time()
    max_block = 0
    for w, plist in weight_to_triples.items():
        block_size = len(plist)
        max_block = max(max_block, block_size)
        weight_dims_full[w] = block_size
        # Build E_w (M x block_size): row t = product over orbit_coefs[t, i, j, k]
        E_w = np.empty((M, block_size), dtype=np.float64)
        for col, (i, j, k) in enumerate(plist):
            E_w[:, col] = orbit_coefs[:, i] * orbit_coefs[:, j] * orbit_coefs[:, k]
        if block_size == 1:
            rank = 1 if np.max(np.abs(E_w)) > 1e-10 else 0
        else:
            s = np.linalg.svd(E_w, compute_uv=False)
            scale = s[0] if len(s) > 0 else 1.0
            tol = max(M, block_size) * np.finfo(float).eps * max(scale, 1.0)
            rank = int(np.sum(s > tol))
        weight_dims_I3[w] = block_size - rank
    if verbose:
        print(f"    [SVD per weight: {time.time()-t_svd:.1f}s, "
              f"max block = {max_block}]")
    dim_I3_total = sum(weight_dims_I3.values())
    return weight_dims_I3, weight_dims_full, dim_I3_total, orbit_coefs


def schur_decompose_from_weights(weight_dim_dict, n, total_deg):
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
    return expand_in_schur_basis(char, syms, total_deg, max_parts=n)


def run_for_polynomial(name, f_expr, syms_full, basis_d, basis_dict,
                       structure, N, n, M, seed, do_schur=True):
    print(f"\n=== {name} ===")
    print(f"  f = {f_expr}")
    d = sum(basis_d[0])
    f_vec = poly_to_coef_vec(f_expr, basis_d, syms_full)
    t0 = time.time()
    weight_dims_I3, weight_dims_full, dim_I3, _ = compute_I_3_weight_dims(
        f_vec, structure, N, n, M, basis_d, seed=seed, verbose=True)
    print(f"  [I_3 compute time: {time.time()-t0:.1f}s]")
    N3 = N * (N + 1) * (N + 2) // 6
    print(f"  dim Sym^3(Sym^{d} V_{n}) = {N3}")
    print(f"  dim I_3(orbit-closure) = {dim_I3}  (per-weight kernel split)")
    print(f"  dim Sym^3 / I_3       = {N3 - dim_I3}")
    weight_dims_quot = {w: weight_dims_full[w] - weight_dims_I3.get(w, 0)
                        for w in weight_dims_full}
    decomp_I3 = {}
    decomp_quot = {}
    decomp_full = {}
    if do_schur:
        print("  Computing Schur decomposition of I_3 and Sym^3 / I_3 ...")
        t0 = time.time()
        decomp_I3 = schur_decompose_from_weights(weight_dims_I3, n, total_deg=3 * d)
        decomp_quot = schur_decompose_from_weights(weight_dims_quot, n, total_deg=3 * d)
        decomp_full = schur_decompose_from_weights(weight_dims_full, n, total_deg=3 * d)
        print(f"  [Schur decomp time: {time.time()-t0:.1f}s]")
        print("  S_lambda multiplicities (Sym^3(Sym^d) | I_3 | quotient):")
        all_lams = sorted(set(decomp_full.keys())
                          | set(decomp_I3.keys())
                          | set(decomp_quot.keys()))
        for lam in all_lams:
            a = decomp_full.get(lam, 0)
            b = decomp_I3.get(lam, 0)
            c = decomp_quot.get(lam, 0)
            marker = "  " if a == b + c else " !!"
            print(f"    {str(lam):<14s}  full={a}  I3={b}  quot={c}{marker}")
    return {
        "name": name,
        "dim_Sym3": N3,
        "dim_I3": dim_I3,
        "dim_quot": N3 - dim_I3,
        "M_samples": M,
        "weight_dims_I3": {str(w): v for w, v in weight_dims_I3.items()},
        "weight_dims_quot": {str(w): v for w, v in weight_dims_quot.items()},
        "schur_full": {str(lam): v for lam, v in decomp_full.items()},
        "schur_I3": {str(lam): v for lam, v in decomp_I3.items()},
        "schur_quot": {str(lam): v for lam, v in decomp_quot.items()},
    }


def driver(n, d, M, n_baseline, seed, do_x1_cubed=True, do_schur=True):
    print(f"=== A7 commit thread, S213: I_3(orbit-closure) sub-frame, "
          f"n={n}, d={d}, M={M} ===\n")

    syms_full = sp.symbols(f"x1:{n+1}")
    basis_d = list(weak_compositions(d, n))
    basis_dict = {tuple(a): i for i, a in enumerate(basis_d)}
    N = len(basis_d)
    print(f"dim Sym^{d} V_{n} = {N}")
    N3 = N * (N + 1) * (N + 2) // 6
    print(f"dim Sym^3(Sym^{d} V_{n}) = {N3}")
    if M < N3 + 50:
        print(f"  WARNING: M = {M} may be too small (recommend M >= {N3 + 50})")

    print(f"\nBuilding GL_{n}-action structure on Sym^{d} V_{n} ...")
    t0 = time.time()
    structure = build_action_structure(n, d, basis_d, basis_dict)
    n_pairs = len(structure)
    n_terms_total = sum(len(v) for v in structure.values())
    print(f"  structure has {n_pairs} non-zero (beta, alpha) pairs,")
    print(f"  total {n_terms_total} (coef, B) entries; build time {time.time()-t0:.1f}s")

    weight_dims_full_ref = weight_dims_of_sym_3_sym_d(basis_d, n)
    print(f"\nSym^3(Sym^{d} V_{n}) weight count: {len(weight_dims_full_ref)} weights, "
          f"sum = {sum(weight_dims_full_ref.values())}")
    if do_schur:
        print("\nSchur decomposition of full Sym^3(Sym^3 V_n) [reference]:")
        ref_decomp = schur_decompose_from_weights(weight_dims_full_ref, n, total_deg=3 * d)
        for lam, m in sorted(ref_decomp.items()):
            print(f"  S_{lam} : mult {m}, dim {schur_dim_in_n_vars(lam, n)}")

    if not (n == 4 and d == 3):
        # Sanity-test mode (n != 4 or d != 3): only run x_1^d Veronese control
        # to verify the framework matches catalecticant theory.
        print(f"\n[SANITY MODE: n={n}, d={d} -- only x_1^{d} Veronese control]")
        results_sanity = {}
        results_sanity["x1_pow_d"] = run_for_polynomial(
            f"x_1^{d} (sanity: cubic Veronese in n={n})",
            syms_full[0] ** d, syms_full, basis_d, basis_dict,
            structure, N, n, M, seed, do_schur=do_schur)
        # Expected: dim I_3 = N3 - dim Sym^{3d} V_n
        from math import comb as _comb
        expected_quot = _comb(3 * d + n - 1, n - 1)
        expected_I3 = N3 - expected_quot
        actual_I3 = results_sanity["x1_pow_d"]["dim_I3"]
        actual_quot = results_sanity["x1_pow_d"]["dim_quot"]
        print(f"\n[SANITY CHECK] dim I_3(x_1^{d}) at n={n}, d={d}:")
        print(f"  expected dim quot = dim Sym^{3*d}(V_{n}) = "
              f"C({3*d}+{n}-1, {n}-1) = {expected_quot}")
        print(f"  expected dim I_3                          = {expected_I3}")
        print(f"  actual   dim I_3                          = {actual_I3}")
        print(f"  actual   dim quot                          = {actual_quot}")
        if actual_I3 == expected_I3:
            print(f"  [SANITY PASS]")
        else:
            print(f"  [SANITY FAIL] discrepancy = {actual_I3 - expected_I3}")
        return {"n": n, "d": d, "M": M, "sanity": True,
                "expected_I3": expected_I3, "actual_I3": actual_I3,
                "results": results_sanity}

    chi_P_support_d3 = [(0, 1, 2), (0, 1, 3), (0, 2, 3)]

    results = {}
    results["chi_P_n4_d3"] = run_for_polynomial(
        "f_chi_P^(4) deg-3 component (= x1 x2 x3 + x1 x2 x4 + x1 x3 x4)",
        chi_P_d3_n4(syms_full), syms_full, basis_d, basis_dict,
        structure, N, n, M, seed, do_schur=do_schur)

    results["e_3"] = run_for_polynomial(
        "e_3 (elementary symmetric deg 3 in 4 vars)",
        e_3_n4(syms_full), syms_full, basis_d, basis_dict,
        structure, N, n, M, seed + 1, do_schur=do_schur)

    results["p_3"] = run_for_polynomial(
        "p_3 (power sum deg 3 in 4 vars)",
        p_3_n4(syms_full), syms_full, basis_d, basis_dict,
        structure, N, n, M, seed + 2, do_schur=do_schur)

    results["plucker"] = run_for_polynomial(
        "x_1 x_2 x_3 (single triple monomial)",
        plucker_n4(syms_full), syms_full, basis_d, basis_dict,
        structure, N, n, M, seed + 3, do_schur=do_schur)

    for k in range(n_baseline):
        f_b = matched_baseline(syms_full, chi_P_support_d3, seed + 10 + k)
        results[f"baseline_{k}"] = run_for_polynomial(
            f"matched-support baseline #{k} (seed {seed + 10 + k})",
            f_b, syms_full, basis_d, basis_dict,
            structure, N, n, M, seed + 100 + k, do_schur=do_schur)

    if do_x1_cubed:
        results["x1_cubed"] = run_for_polynomial(
            "x_1^3 (control: cubic Veronese)",
            x1_cubed(syms_full), syms_full, basis_d, basis_dict,
            structure, N, n, M, seed + 200, do_schur=do_schur)

    # ---- Cross-comparison summary ----
    print("\n\n=========== CROSS-COMPARISON ===========")
    print("\n--- dim I_3 per polynomial ---")
    for name, r in results.items():
        print(f"  {name:35s}  dim_I3 = {r['dim_I3']:5d}, "
              f"dim_quot = {r['dim_quot']:5d}")

    if do_schur:
        all_keys = set()
        for r in results.values():
            all_keys.update(r['schur_quot'].keys())
        all_keys_sorted = sorted(all_keys)
        chi_P_key = "chi_P_n4_d3"
        baseline_keys = [k for k in results if k.startswith("baseline_")]
        other_keys = [k for k in results
                      if not (k.startswith("chi_P") or k.startswith("baseline_"))]

        print("\n--- Schur multiplicity in Sym^3 / I_3 (quotient) ---")
        header_cols = [("chi_P", chi_P_key)]
        for k in other_keys:
            header_cols.append((k[:10], k))
        for k in baseline_keys:
            header_cols.append((k.replace("baseline_", "b"), k))
        header = "  " + f"{'lambda':<16s}" + " | ".join(f"{h:>10s}" for h, _ in header_cols)
        print(header)
        print("  " + "-" * len(header))

        occurrence_signals = []
        for lam_str in all_keys_sorted:
            row_vals = [results[k]['schur_quot'].get(lam_str, 0)
                        for _, k in header_cols]
            chi_P_val = results[chi_P_key]['schur_quot'].get(lam_str, 0)
            baseline_vals = [results[k]['schur_quot'].get(lam_str, 0)
                             for k in baseline_keys]
            flag = ""
            if baseline_vals:
                if chi_P_val > 0 and all(b == 0 for b in baseline_vals):
                    flag += " [chi_P > 0, all_baseline = 0]"
                if chi_P_val == 0 and any(b > 0 for b in baseline_vals):
                    flag += " [chi_P = 0, baseline > 0]"
                if (chi_P_val != 0 and any(b != chi_P_val for b in baseline_vals)
                    and not all(b == 0 for b in baseline_vals)):
                    flag += " [chi_P != some baseline]"
            if flag:
                occurrence_signals.append(
                    (lam_str, flag, dict(zip([h for h, _ in header_cols], row_vals))))
            print(f"  {lam_str:<16s}" + " | ".join(f"{v:10d}" for v in row_vals) + flag)

        print("\n--- Occurrence-obstruction candidates ---")
        if occurrence_signals:
            print(f"  Found {len(occurrence_signals)} candidates:")
            for lam_str, flag, vals in occurrence_signals:
                print(f"    lambda={lam_str}: {vals}{flag}")
        else:
            print(f"  NONE found across {n_baseline} baselines + reference polynomials")

    return {
        "n": n, "d": d, "M": M, "n_baseline": n_baseline,
        "dim_Sym_d_V": N, "dim_Sym_3_Sym_d_V": N3,
        "results": results,
    }


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    d = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    M = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
    n_baseline = int(sys.argv[4]) if len(sys.argv) > 4 else 3
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 4242
    do_schur = bool(int(sys.argv[6])) if len(sys.argv) > 6 else True

    out = driver(n=n, d=d, M=M, n_baseline=n_baseline, seed=seed, do_schur=do_schur)
    out_path = os.path.join(
        os.path.dirname(__file__) or ".",
        f"plethysm_subframe_I3_n{n}_d{d}_results.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\n[done] Wrote {out_path}")
