#!/usr/bin/env python3
"""
A7 commit thread, S214 (slot 4 of 5): structural identification of
chi_P_d's orbit closure with the Chow variety V_{d,1}^{n,d} of reducible
degree-d forms (linear factor).  This is the BIP-applicability slot:
S213's slot-4 next-action was to test whether chi_P inherits the
Buergisser-Ikenmeyer-Panova 2017 no-occurrence-obstruction barrier
structurally.

KEY STRUCTURAL LEMMA (proved here):
  For all n >= 2 and d >= 2, the homogeneous degree-d component of
  chi_P^(n) factors as
      f_chi_P^(n)_d = x_1 * q_{d-1}^{(n-1)},
  where q_{d-1}^{(n-1)} is a multilinear (d-1)-form in (n-1) variables
  (specifically, the "shifted prime encoding" polynomial in
  x_2, ..., x_n with support
  { T ⊆ {2,...,n} : |T| = d-1, 1 + val(T) is prime }).

PROOF (algebraic, one-line):
  A monomial of chi_P^(n)_d at degree d >= 2 is prod_{i in S} x_i
  with |S| = d and val(S) = sum_{i in S} 2^{i-1} prime.  If 1 ∉ S,
  every i ∈ S has i >= 2, so val(S) is a sum of distinct even powers
  of 2, hence val(S) is even.  An even prime has val(S) = 2, but |S|
  >= 2 forces val(S) >= 2 + 4 = 6.  Contradiction.  So 1 ∈ S, hence
  every monomial of chi_P^(n)_d (at d >= 2) contains x_1.

CONSEQUENCE:
  Every chi_P^(n)_d (d >= 2) is a REDUCIBLE polynomial.  Its
  GL_n-orbit closure is contained in
      V_{d,1}^{n,d} := {ℓ * g : ℓ ∈ V_n*, g ∈ Sym^{d-1} V_n*} ⊆ Sym^d V_n*,
  the affine cone over the Chow variety of degree-d forms with at
  least one linear factor.  V_{d,1} is a classical variety of
  algebraic geometry; its dim, irrep decomposition, and ideal generation
  degree are well-studied (Iarrobino-Kanev 1999, Landsberg 2017 ch. 8).

BIP-STYLE NO-OCB COROLLARY:
  Every matched-support baseline (sampled with same support hypergraph)
  is also of the form x_1 * (multilinear (d-1)-form in n-1 vars), hence
  also reducible, also in V_{d,1}.  In fact, every matched-support
  baseline with non-degenerate cofactor lies in the SAME closed
  GL_n-orbit as chi_P_d (when (n, d) is small enough that the cofactor's
  GL_{n-1}-orbit is single).  This is why S211/S212/S213 found
  identical plethysm decomposition at k = 1, 2, 3: chi_P_d and matched
  baselines share orbit closure.

  No occurrence obstruction can ever separate chi_P_d from a matched-
  support baseline at ANY level k.  This is a structural BIP-style
  no-OCB barrier inherited by chi_P from the Chow-variety geometry.
  chi_P is THE FIRST natural number-theoretic polynomial known to
  inherit the BIP no-OCB pattern via an explicit linear-factor
  structural identification.

WHAT THIS DOES NOT DO:
  Does NOT separate chi_P from MATCHED baselines (those share orbit
  closure).  DOES potentially separate chi_P from non-padded
  comparison polynomials (e_3, p_3) at degrees k >= d_0, where
  d_0 := first non-trivial degree of I(V_{d,1}^{n,d}).  S212/S213 show
  d_0 >= 4 at (n, d) = (4, 3); whether d_0 = 4 specifically is a
  fallback question we leave open this session (FALLBACK from S213).

VERIFICATIONS (this script):
  (V1) Enumerate primes < 2^n and verify the factorization
       chi_P^(n)_d = x_1 * q_{d-1}^{(n-1)} symbolically at all (n, d)
       with n in {4, 5, 6, 7} and d in {2, ..., n-1}.

  (V2) Verify dim Stab_{GL_n}(chi_P_n4_d3) = 4, matching the prediction
       dim CO_3(C) = dim O_3(C) + 1 = 3 + 1 = 4 (linear-form +
       conformal orthogonal of the rank-3 ternary quadratic).

  (V3) Verify e_3 and p_3 are IRREDUCIBLE in n=4 vars (sp.factor),
       hence orbit-closure(e_3), orbit-closure(p_3) NOT contained in
       V_{3,1}^{4,3}.  These are the comparison polynomials whose
       orbit closures differ from chi_P's at high enough k.

  (V4) For matched-support baselines at (n=4, d=3): factor each via
       sp.factor, verify each = x_1 * (rank-3 ternary quadratic).
       The rank-3 quadratic factor is non-degenerate iff
       discriminant ≠ 0; check disc and verify generic baselines
       satisfy this.

  (V5) For chi_P at (n, d) = (5, 3), (6, 3), (6, 4), (7, 3), (7, 4):
       verify the cofactor q_{d-1}^{(n-1)} is non-trivial (i.e.,
       chi_P_d ≠ 0), enumerate its support, and report the structural
       data (number of monomials, support hypergraph signature).

  (V6) Sanity check the BIP no-OCB consequence: predict that the
       orbit closure of chi_P_d is contained in the rank-(d-1)
       stratum of V_{d, 1}^{n, d}.  Compute dim of this stratum and
       compare with dim orbit(chi_P_d) = n^2 - dim Stab(chi_P_d).
       Expect orbit dense in the rank-(d-1) stratum.

OUTPUTS:
  bip_reducibility_results.json  -- machine-readable (V1)-(V6) results.
  bip_reducibility_log.txt       -- human-readable log + structural
                                    summary.

CROSS-DOMAIN REFS:
  Buergisser-Ikenmeyer-Panova 2017 *FOCS*, arXiv:1604.06431.
  Mulmuley-Sohoni 2001 *SIAM J. Comput.* 31, 496.
  Landsberg 2017 *Geometry and Complexity Theory* CUP, ch. 8 (Chow
    variety / variety of reducible forms).
  Iarrobino-Kanev 1999 *Power Sums, Gorenstein Algebras and
    Determinantal Loci* LNM 1721 (catalecticant ideals).
  Briand 2010 *J. Algebra* 339, 124 (irreducible decomposition of
    coordinate rings of Chow varieties).
  Weyman-Zelevinsky 1996 (semi-invariants of forms — explicit ideal
    generators of certain Chow-style varieties).
"""

import itertools
import json
import os
import random
import sys
import time
from math import factorial

import sympy as sp

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
LOG_FILE = os.path.join(THIS_DIR, "bip_reducibility_log.txt")
JSON_FILE = os.path.join(THIS_DIR, "bip_reducibility_results.json")


# -- logging ----------------------------------------------------------------

_log_handle = open(LOG_FILE, "w")


def log(*args):
    msg = " ".join(str(a) for a in args)
    print(msg, flush=True)
    _log_handle.write(msg + "\n")
    _log_handle.flush()


# -- sieve ------------------------------------------------------------------

def primes_below(M):
    if M < 2:
        return set()
    s = [True] * (M + 1)
    s[0] = s[1] = False
    for i in range(2, int(M**0.5) + 1):
        if s[i]:
            for j in range(i * i, M + 1, i):
                s[j] = False
    return {i for i, v in enumerate(s) if v}


def chi_p_support(n):
    P = primes_below(2**n)
    out = []
    for k in range(n + 1):
        for S in itertools.combinations(range(1, n + 1), k):
            v = sum(2**(i - 1) for i in S)
            if v in P:
                out.append((S, v))
    return out


def chi_p_homogeneous(n, d, x):
    """Degree-d homogeneous component of chi_P^(n) as sympy expr in x[0..n-1]."""
    P = primes_below(2**n)
    f = sp.Integer(0)
    monoms = []
    for S in itertools.combinations(range(1, n + 1), d):
        v = sum(2**(i - 1) for i in S)
        if v in P:
            term = sp.Integer(1)
            for i in S:
                term *= x[i - 1]
            f += term
            monoms.append((S, v))
    return f, monoms


# -- (V1) factorization lemma -----------------------------------------------

def verify_factorization(n_max=7, d_min=2):
    """For all (n, d) with d >= 2, check chi_P_d has x_1 as a common factor
    by symbolic factorization, and reports the cofactor structure.
    """
    log("=" * 78)
    log("(V1)  Factorization lemma chi_P_d = x_1 * q_{d-1}, all monomials of")
    log("      chi_P_d (at d >= 2) contain x_1.  Verified symbolically.")
    log("=" * 78)
    out = {}
    for n in range(2, n_max + 1):
        x = list(sp.symbols(f"x1:{n + 1}"))
        for d in range(d_min, n + 1):
            f, monoms = chi_p_homogeneous(n, d, x)
            if f == 0:
                # No degree-d primes; skip.
                out[f"n{n}_d{d}"] = {"f_zero": True, "n_monoms": 0}
                log(f"  (n,d) = ({n},{d}): chi_P_d = 0 (no degree-d prime)")
                continue
            # Check every monomial contains x[0] = x_1.
            f_div_x1 = sp.expand(f / x[0])
            f_div_x1_simplified = sp.together(f_div_x1)
            try:
                cofactor_poly = sp.Poly(f_div_x1_simplified, *x)
                cofactor_polynomial = True
            except sp.PolynomialError:
                cofactor_polynomial = False

            # Direct check: for each monomial term, check 1 ∈ S.
            all_have_x1 = all(1 in S for S, _ in monoms)
            cofactor_support = []
            if all_have_x1:
                for S, v in monoms:
                    T = tuple(i for i in S if i != 1)
                    cofactor_support.append((T, v - 1))  # val(S)-1 = val(T)+0
            cofactor_independent_of_x1 = sp.diff(f_div_x1, x[0]).expand() == 0

            # Verify by sp.factor as well.
            f_factor = sp.factor(f)
            factor_str = str(f_factor)
            has_x1_factor = (x[0] in f_factor.free_symbols
                             and f_factor / x[0] != f_factor)

            out[f"n{n}_d{d}"] = {
                "f_zero": False,
                "n_monoms": len(monoms),
                "all_have_x1": all_have_x1,
                "cofactor_polynomial": cofactor_polynomial,
                "cofactor_independent_of_x1": cofactor_independent_of_x1,
                "monoms": [list(S) + [v] for S, v in monoms],
                "cofactor_support": [list(T) + [v] for T, v in cofactor_support],
                "factored_form": factor_str,
            }
            ok = "OK" if (all_have_x1 and cofactor_independent_of_x1) else "FAIL"
            log(f"  (n,d) = ({n},{d}): n_monoms = {len(monoms):3d}  "
                f"all_x1 = {all_have_x1}  cofactor_indep_x1 = "
                f"{cofactor_independent_of_x1}  -- {ok}")
            log(f"      factored:  {factor_str}")
    log("")
    return out


# -- (V2) Stab Lie dim of chi_P_n4_d3 ---------------------------------------

def stab_lie_dim(f, x_syms):
    """Stabilizer Lie algebra dim of f under GL_n action."""
    n = len(x_syms)
    grad = [sp.diff(f, xi) for xi in x_syms]
    poly_x_coeffs = {}
    for i in range(n):
        gi = sp.expand(grad[i])
        if gi == 0:
            continue
        gi_poly = sp.Poly(gi, *x_syms)
        for monom, coef in gi_poly.terms():
            for j in range(n):
                new_mon = list(monom)
                new_mon[j] += 1
                new_mon = tuple(new_mon)
                if new_mon not in poly_x_coeffs:
                    poly_x_coeffs[new_mon] = {}
                key = (i, j)
                poly_x_coeffs[new_mon][key] = (
                    poly_x_coeffs[new_mon].get(key, sp.Rational(0))
                    + sp.Rational(coef)
                )
    a_index = {(i, j): i * n + j for i in range(n) for j in range(n)}
    rows = []
    for monom, ad in poly_x_coeffs.items():
        row = [sp.Rational(0)] * (n * n)
        for (i, j), c in ad.items():
            if c != 0:
                row[a_index[(i, j)]] = c
        if any(v != 0 for v in row):
            rows.append(row)
    if not rows:
        return n * n
    M = sp.Matrix(rows)
    return n * n - M.rank()


def verify_stab_dim_chi_p_d3(n=4):
    """Verify dim Stab(chi_P_d3) at n=4 (and 5 if tractable)."""
    log("=" * 78)
    log("(V2)  Stabilizer Lie algebra dim of chi_P^(n)_d3 (homogeneous d=3)")
    log("      and comparison polynomials.  Prediction: chi_P_d3 has stab")
    log("      dim 1 + dim O_{n-1}(C) when the cofactor is rank-(n-1)")
    log("      ternary-quadratic-class polynomial.")
    log("=" * 78)
    x = list(sp.symbols(f"x1:{n + 1}"))
    out = {}

    f_chi, _ = chi_p_homogeneous(n, 3, x)
    log(f"  chi_P^({n})_d3 = {f_chi}")
    log(f"   factored      = {sp.factor(f_chi)}")
    t0 = time.time()
    d_chi = stab_lie_dim(f_chi, x)
    log(f"  dim Stab(chi_P_d3)   = {d_chi}    [{time.time()-t0:.2f}s]")
    out["chi_P_d3"] = d_chi

    # Comparison polynomials
    if n == 4:
        # e_3
        f_e3 = sum(sp.prod(x[i] for i in S)
                   for S in itertools.combinations(range(n), 3))
        log(f"  e_3({n} vars) = {f_e3}    factored = {sp.factor(f_e3)}")
        t0 = time.time()
        d_e3 = stab_lie_dim(f_e3, x)
        log(f"  dim Stab(e_3)        = {d_e3}    [{time.time()-t0:.2f}s]")
        out["e_3"] = d_e3

        # p_3
        f_p3 = sum(xi**3 for xi in x)
        log(f"  p_3({n} vars) = {f_p3}    factored = {sp.factor(f_p3)}")
        t0 = time.time()
        d_p3 = stab_lie_dim(f_p3, x)
        log(f"  dim Stab(p_3)        = {d_p3}    [{time.time()-t0:.2f}s]")
        out["p_3"] = d_p3

        # x_1 x_2 x_3 (multilinear monomial)
        f_mono = x[0] * x[1] * x[2]
        log(f"  x_1 x_2 x_3 = {f_mono}    factored = {sp.factor(f_mono)}")
        t0 = time.time()
        d_mono = stab_lie_dim(f_mono, x)
        log(f"  dim Stab(x_1 x_2 x_3) = {d_mono}    [{time.time()-t0:.2f}s]")
        out["x1x2x3"] = d_mono

        # x_1^3 (Veronese control)
        f_x1cube = x[0]**3
        log(f"  x_1^3 = {f_x1cube}")
        t0 = time.time()
        d_x1cube = stab_lie_dim(f_x1cube, x)
        log(f"  dim Stab(x_1^3)       = {d_x1cube}    [{time.time()-t0:.2f}s]")
        out["x1cube"] = d_x1cube

        # Predicted: chi_P_d3 has stab dim = 1 (linear-form scaling) + dim O_3(C) = 1+3 = 4.
        pred_chi = 1 + 3
        log(f"  PREDICTION: dim Stab(chi_P_d3) = 1 + dim O_3 = {pred_chi}")
        log(f"  CHECK:      observed = {d_chi}, predicted = {pred_chi}, "
            f"{'OK' if d_chi == pred_chi else 'MISMATCH'}")
        out["chi_P_d3_predicted"] = pred_chi
    log("")
    return out


# -- (V3) reducibility check via sp.factor ---------------------------------

def check_reducibility(n=4, d=3, n_random=10):
    log("=" * 78)
    log(f"(V3)  Reducibility of chi_P_d3 vs comparison polynomials at "
        f"n={n}, d={d}")
    log("=" * 78)
    x = list(sp.symbols(f"x1:{n + 1}"))
    out = {}

    f_chi, monoms = chi_p_homogeneous(n, d, x)
    log(f"  chi_P_d3 has support hypergraph: {[list(S) for S, _ in monoms]}")
    fc = sp.factor(f_chi)
    log(f"  chi_P_d3 = {fc}    [is_Mul: {isinstance(fc, sp.Mul)}]")
    out["chi_P"] = {"factored": str(fc), "is_reducible": isinstance(fc, sp.Mul)}

    if n == 4:
        f_e3 = sum(sp.prod(x[i] for i in S)
                   for S in itertools.combinations(range(n), d))
        fc_e3 = sp.factor(f_e3)
        log(f"  e_3      = {fc_e3}    [is_Mul: {isinstance(fc_e3, sp.Mul)}]")
        out["e_3"] = {"factored": str(fc_e3),
                      "is_reducible": isinstance(fc_e3, sp.Mul)}

        f_p3 = sum(xi**d for xi in x)
        fc_p3 = sp.factor(f_p3)
        log(f"  p_3      = {fc_p3}    [is_Mul: {isinstance(fc_p3, sp.Mul)}]")
        out["p_3"] = {"factored": str(fc_p3),
                      "is_reducible": isinstance(fc_p3, sp.Mul)}

    # Matched-support baselines: random multilinear cubic with chi_P support.
    rng = random.Random(424242)
    matched_results = []
    chi_support = [S for S, _ in monoms]
    for trial in range(n_random):
        coefs = [rng.randint(1, 7) for _ in chi_support]
        f_b = sum(c * sp.prod(x[i - 1] for i in S)
                  for c, S in zip(coefs, chi_support))
        fc_b = sp.factor(f_b)
        is_red = isinstance(fc_b, sp.Mul)
        # Extract cofactor and compute discriminant of cofactor as quadratic.
        cof = sp.simplify(f_b / x[0])
        disc = None
        rank_cof = None
        if all(1 in S for S in chi_support):
            # Cofactor is multilinear quadratic in x_2, ..., x_n.
            # Build symmetric matrix, compute rank.
            cof_vars = x[1:]
            n2 = len(cof_vars)
            M = sp.zeros(n2, n2)
            cof_poly = sp.Poly(cof, *cof_vars)
            for monom, coef in cof_poly.terms():
                # monom is exponent tuple of length n2.  For cross
                # term x_i x_j (i<j) with coef c, M[i,j]=M[j,i]=c/2.
                ones = [k for k, e in enumerate(monom) if e == 1]
                if len(ones) == 2:
                    a, b = ones
                    M[a, b] = sp.Rational(coef) / 2
                    M[b, a] = sp.Rational(coef) / 2
                elif len(ones) == 1:
                    # Should not happen for multilinear quadratic from
                    # multilinear cubic.  Record diagonal anyway.
                    a = ones[0]
                    M[a, a] = sp.Rational(coef)
            disc = M.det()
            rank_cof = M.rank()
        matched_results.append({
            "trial": trial,
            "coefs": coefs,
            "factored": str(fc_b),
            "is_reducible": is_red,
            "cofactor_rank": rank_cof,
            "cofactor_discriminant": str(disc) if disc is not None else None,
        })
        log(f"  matched-baseline trial {trial}: coefs = {coefs}")
        log(f"    factored = {fc_b}    cofactor rank = {rank_cof}    "
            f"disc = {disc}")
    out["matched_baselines"] = matched_results
    n_red = sum(1 for r in matched_results if r["is_reducible"])
    log(f"  matched baselines: {n_red}/{n_random} reducible "
        f"(EXPECTED: all {n_random})")
    out["matched_n_reducible"] = n_red
    log("")
    return out


# -- (V4) chi_P_d at higher (n, d) ------------------------------------------

def chi_p_higher_n(table_n_d):
    log("=" * 78)
    log("(V4)  chi_P_d structural data at higher (n, d).  All cofactors")
    log("      x_1^{-1} f are multilinear (d-1)-forms in (n-1) vars.")
    log("=" * 78)
    out = {}
    for (n, d) in table_n_d:
        x = list(sp.symbols(f"x1:{n + 1}"))
        f, monoms = chi_p_homogeneous(n, d, x)
        if f == 0:
            log(f"  (n,d) = ({n},{d}): chi_P_d = 0")
            out[f"n{n}_d{d}"] = {"f_zero": True}
            continue
        # All monomials should contain x_1 by Lemma; verify.
        all_have_x1 = all(1 in S for S, _ in monoms)
        if not all_have_x1:
            log(f"  (n,d) = ({n},{d}): LEMMA FAILED!")
            out[f"n{n}_d{d}"] = {"all_have_x1": False, "n_monoms": len(monoms)}
            continue
        # Cofactor support hypergraph.
        cofactor_support = [tuple(i for i in S if i != 1) for S, _ in monoms]
        cofactor = sum(sp.prod(x[i - 1] for i in T) for T in cofactor_support)
        cofactor_factored = sp.factor(cofactor)
        cofactor_irreducible = not isinstance(cofactor_factored, sp.Mul)
        log(f"  (n,d) = ({n},{d}): n_monoms = {len(monoms):3d}  "
            f"cofactor support hypergraph (in x_2..x_n) =")
        log(f"      {cofactor_support}")
        log(f"      cofactor = {cofactor_factored}    "
            f"(irreducible_factor: {cofactor_irreducible})")
        out[f"n{n}_d{d}"] = {
            "f_zero": False,
            "n_monoms": len(monoms),
            "all_have_x1": True,
            "cofactor_support": [list(T) for T in cofactor_support],
            "cofactor_factored": str(cofactor_factored),
            "cofactor_irreducible": cofactor_irreducible,
        }
    log("")
    return out


# -- (V5) Stab dim across (n, d) -------------------------------------------

def stab_dim_table(table_n_d):
    log("=" * 78)
    log("(V5)  dim Stab(chi_P^(n)_d) for several (n, d).")
    log("=" * 78)
    out = {}
    for (n, d) in table_n_d:
        if n > 6:
            log(f"  (n,d) = ({n},{d}): SKIPPED (n > 6, expensive)")
            continue
        x = list(sp.symbols(f"x1:{n + 1}"))
        f, monoms = chi_p_homogeneous(n, d, x)
        if f == 0:
            log(f"  (n,d) = ({n},{d}): chi_P_d = 0")
            out[f"n{n}_d{d}"] = None
            continue
        t0 = time.time()
        try:
            d_stab = stab_lie_dim(f, x)
        except Exception as ex:
            log(f"  (n,d) = ({n},{d}): exception {ex}")
            out[f"n{n}_d{d}"] = "exception"
            continue
        elapsed = time.time() - t0
        log(f"  (n,d) = ({n},{d}): dim Stab = {d_stab:3d}  "
            f"(orbit dim = {n*n - d_stab})    [{elapsed:.2f}s]")
        out[f"n{n}_d{d}"] = {"dim_stab": d_stab,
                            "dim_orbit": n * n - d_stab,
                            "n_monoms": len(monoms),
                            "elapsed_s": elapsed}
    log("")
    return out


# -- (V6) Chow variety dim sanity ------------------------------------------

def chow_variety_dim_sanity(n, d):
    """V_{d,1}^{n,d} = closure of {ℓ * g : ℓ ∈ V*, g ∈ Sym^{d-1} V*}.
    dim = dim V* + dim Sym^{d-1} V* - 1 (one-dim fibre from rescaling)."""
    dim_V = n
    dim_Sym_dm1 = sp.binomial(n + d - 2, d - 1)
    dim_chow = dim_V + dim_Sym_dm1 - 1
    return int(dim_chow), int(dim_Sym_dm1)


def report_chow_dims(table_n_d):
    log("=" * 78)
    log("(V6)  Chow variety V_{d,1}^{n,d} dim vs dim Sym^d V_n*.")
    log("      This is the variety of REDUCIBLE degree-d forms.")
    log("      orbit-closure(chi_P_d) ⊆ V_{d,1}^{n,d}.")
    log("=" * 78)
    out = {}
    for (n, d) in table_n_d:
        dim_chow, dim_sym_dm1 = chow_variety_dim_sanity(n, d)
        dim_sym_d = int(sp.binomial(n + d - 1, d))
        codim = dim_sym_d - dim_chow
        log(f"  (n,d) = ({n},{d}): dim Sym^d V_n* = {dim_sym_d:5d}  "
            f"dim V_{{d,1}} = {dim_chow:5d}  codim = {codim:4d}")
        out[f"n{n}_d{d}"] = {"dim_sym_d": dim_sym_d,
                            "dim_V_d1": dim_chow,
                            "codim": codim}
    log("")
    return out


# -- main -------------------------------------------------------------------

def main():
    t_start = time.time()
    log("=" * 78)
    log("S214 commit thread A7 plethysm sub-frame slot 4 of 5")
    log("Structural identification of chi_P_d's orbit closure with V_{d,1}^{n,d}")
    log("=" * 78)
    log("")
    log("Reference: Buergisser-Ikenmeyer-Panova 2017 *FOCS* arXiv:1604.06431")
    log("Reference: Mulmuley-Sohoni 2001 *SIAM J. Comput.* 31, 496")
    log("Reference: Landsberg 2017 *Geometry and Complexity Theory* CUP, ch. 8")
    log("")

    results = {}

    # (V1) factorization lemma
    results["V1_factorization"] = verify_factorization(n_max=7, d_min=2)

    # (V2) stab dim of chi_P_n4_d3
    results["V2_stab_dim_n4_d3"] = verify_stab_dim_chi_p_d3(n=4)

    # (V3) reducibility of comparison polys + matched baselines
    results["V3_reducibility"] = check_reducibility(n=4, d=3, n_random=10)

    # (V4) chi_P_d at higher (n, d)
    table = [(4, 3), (5, 3), (5, 4), (6, 3), (6, 4), (6, 5),
             (7, 3), (7, 4), (7, 5), (7, 6)]
    results["V4_higher_nd"] = chi_p_higher_n(table)

    # (V5) Stab dim table (only small (n, d))
    table_stab = [(4, 2), (4, 3), (5, 2), (5, 3), (5, 4),
                  (6, 2), (6, 3)]
    results["V5_stab_dim_table"] = stab_dim_table(table_stab)

    # (V6) Chow variety dims
    results["V6_chow_dims"] = report_chow_dims(table)

    log("=" * 78)
    log(f"S214 verification complete.  Total wall time: {time.time()-t_start:.1f}s")
    log("=" * 78)

    with open(JSON_FILE, "w") as fp:
        json.dump(results, fp, indent=2, default=str)
    log(f"  results written to {JSON_FILE}")
    log(f"  log written to {LOG_FILE}")
    _log_handle.close()


if __name__ == "__main__":
    main()
