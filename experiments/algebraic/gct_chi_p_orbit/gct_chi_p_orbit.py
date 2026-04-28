#!/usr/bin/env python3
"""
A7 — Geometric Complexity Theory invariants of the prime-encoding polynomial f_chi_P^{(n)}

Computes for the multi-affine polynomial

    f_chi_P^{(n)}(x_1, ..., x_n) := sum_{S subseteq [n], val(S) in PRIMES} prod_{i in S} x_i,
    val(S) := sum_{i in S} 2^{i-1},

the following invariants:
  1. dim of the GL_n-stabilizer Lie algebra (linear algebra in n^2 unknowns)
  2. orbit dim = n^2 - dim Stab Lie
  3. discrete S_n-permutation subgroup fixing f (enumerated by brute force)
  4. diagonal-torus stabilizer dimension
  5. partial-derivative-space dim (= dim span(grad f))
  6. Hessian rank (symbolic + at a random rational point)
  7. matched-support random multi-affine baseline (B-grade comparison via z-score)

References:
  * Mulmuley-Sohoni 2001, SIAM J. Comput. 31, 496-526.
  * Buergisser-Ikenmeyer-Panova 2017, FOCS, arXiv:1604.06431.
  * Landsberg 2017, "Geometry and Complexity Theory" CUP, ch. 9-10.

Sanity checks (in same script):
  * det_2 on 4 variables: dim Stab Lie alg = 2*4 - 1 = 7 (Frobenius 1897).
  * e_2 on 4 variables: dim Stab Lie alg = dim o(4) = 6.

Pre-stated falsification criteria (registered before measurement):
  FAL-1: dim Stab >= n^2/2  =>  collapse-symmetric (mode I).
  FAL-2: |perm group| = n!  =>  symmetric-poly collapse (mode E).
  FAL-3: PDS rank = n       =>  sanity (almost-trivial).
  FAL-4 (REAL): |z(f_n vs random)| > 3  =>  arithmetic stab signal (A direction).
                |z| <= 3  =>  matched baseline => B-grade closure.
  FAL-5: rank H(f_n) < n at random pts => Hessian obstruction (A signal).
"""

import itertools
import json
import random
import sys

import sympy as sp


# ---------------------------------------------------------------------------
# Polynomials
# ---------------------------------------------------------------------------

def sieve_primes(M):
    if M < 2:
        return set()
    s = [True] * (M + 1)
    s[0] = s[1] = False
    for i in range(2, int(M**0.5) + 1):
        if s[i]:
            for j in range(i * i, M + 1, i):
                s[j] = False
    return {i for i, x in enumerate(s) if x}


def chi_p_poly(n, x_syms):
    primes = sieve_primes(2**n)
    f = sp.Integer(0)
    monos = []
    for k in range(n + 1):
        for S in itertools.combinations(range(1, n + 1), k):
            v = sum(2**(i - 1) for i in S)
            if v in primes:
                term = sp.Integer(1)
                for i in S:
                    term *= x_syms[i - 1]
                f += term
                monos.append((S, v))
    return f, monos


def chi_p_support(n):
    primes = sieve_primes(2**n)
    return [S for k in range(n + 1)
            for S in itertools.combinations(range(1, n + 1), k)
            if sum(2**(i - 1) for i in S) in primes]


# ---------------------------------------------------------------------------
# Invariants
# ---------------------------------------------------------------------------

def stab_lie_dim(f, x_syms, return_basis=False):
    """
    Stabilizer Lie algebra of f under GL_n action on x. The infinitesimal
    action of A in Mat_n is (A.f)(x) = sum_{i,j} A_{ij} x_j (df/dx_i).
    Stabilizer = {A : A.f == 0 identically}. Linear constraints in A_{ij}.
    """
    n = len(x_syms)
    grad = [sp.diff(f, xi) for xi in x_syms]

    # Build Af as polynomial in x with coefficients linear in A_ij
    # Avoid creating sympy Symbols for A: use direct collection.
    poly_x_coeffs = {}  # exponent_tuple -> dict[(i,j) -> rational coef]
    for i in range(n):
        gi = sp.expand(grad[i])
        if gi == 0:
            continue
        gi_poly = sp.Poly(gi, *x_syms)
        for monom, coef in gi_poly.terms():
            # multiply by x_j for each j: (j-th variable adds 1 to monom[j])
            for j in range(n):
                new_mon = list(monom)
                new_mon[j] += 1
                new_mon = tuple(new_mon)
                if new_mon not in poly_x_coeffs:
                    poly_x_coeffs[new_mon] = {}
                key = (i, j)
                poly_x_coeffs[new_mon][key] = poly_x_coeffs[new_mon].get(key, sp.Rational(0)) + sp.Rational(coef)

    # Each x-monomial yields one constraint row of length n^2.
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
        if return_basis:
            return n * n, []
        return n * n
    M = sp.Matrix(rows)
    rk = M.rank()
    dim_stab = n * n - rk
    if return_basis:
        nullspace = M.nullspace()
        basis = [sp.Matrix(n, n, lambda i, j, v=v: v[i * n + j]) for v in nullspace]
        return dim_stab, basis
    return dim_stab


def discrete_perm_group(f, x_syms):
    n = len(x_syms)
    fixers = []
    for sigma in itertools.permutations(range(n)):
        sub = {x_syms[i]: x_syms[sigma[i]] for i in range(n)}
        if sp.expand(f.subs(sub, simultaneous=True) - f) == 0:
            fixers.append(sigma)
    return fixers


def diag_torus_stab_dim(support, n):
    if not support:
        return n
    rows = []
    for S in support:
        row = [0] * n
        for i in S:
            row[i - 1] = 1
        rows.append(row)
    return n - sp.Matrix(rows).rank()


def pds_dim(f, x_syms):
    derivs = [sp.expand(sp.diff(f, xi)) for xi in x_syms]
    monoms = set()
    for d in derivs:
        if d == 0:
            continue
        for m in sp.Poly(d, *x_syms).monoms():
            monoms.add(m)
    monom_list = sorted(monoms)
    if not monom_list:
        return 0
    M = []
    for d in derivs:
        if d == 0:
            row = [sp.Rational(0)] * len(monom_list)
        else:
            dp = sp.Poly(d, *x_syms)
            row = [sp.Rational(dp.coeff_monomial(m)) for m in monom_list]
        M.append(row)
    return sp.Matrix(M).rank()


def hessian_ranks(f, x_syms, seed=0):
    n = len(x_syms)
    H = sp.Matrix(n, n, lambda i, j: sp.diff(f, x_syms[i], x_syms[j]))
    sym_rank = H.rank()
    rng = random.Random(seed)
    pt = {x: sp.Rational(rng.randint(1, 50), rng.randint(1, 10)) for x in x_syms}
    pt_rank = H.subs(pt).rank()
    return sym_rank, pt_rank


# ---------------------------------------------------------------------------
# Baselines
# ---------------------------------------------------------------------------

def random_polynomial_with_support(n, support, rng):
    x_syms = sp.symbols(f'x_1:{n + 1}')
    g = sp.Integer(0)
    for S in support:
        c = rng.randint(1, 7)
        term = sp.Integer(c)
        for i in S:
            term *= x_syms[i - 1]
        g += term
    return g, x_syms


def baseline_stab_dims(n, support, n_trials, seed=0):
    rng = random.Random(seed)
    dims = []
    for _ in range(n_trials):
        g, xs = random_polynomial_with_support(n, support, rng)
        dims.append(stab_lie_dim(g, xs))
    return dims


# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------

def sanity_det2():
    """det_2 = a*d - b*c on 4 variables; expected dim Stab Lie = 7."""
    a, b, c, d = sp.symbols('a b c d')
    f = a * d - b * c
    return stab_lie_dim(f, [a, b, c, d])


def sanity_e2_n4():
    """e_2(x1..x4) = sum_{i<j} x_i x_j; expected dim Stab Lie = dim o(4) = 6."""
    xs = sp.symbols('x1:5')
    f = sum(xs[i] * xs[j] for i in range(4) for j in range(i + 1, 4))
    return stab_lie_dim(f, xs)


def sanity_perm2():
    """perm_2 = a*d + b*c on 4 vars; expected dim Stab Lie = 5 (n^2 + n - 1 = 5)."""
    a, b, c, d = sp.symbols('a b c d')
    f = a * d + b * c
    return stab_lie_dim(f, [a, b, c, d])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    log_path = "/apps/aplikacijos/prime-research/experiments/algebraic/gct_chi_p_orbit/gct_chi_p_orbit_log.txt"
    out = open(log_path, "w")

    def L(msg):
        print(msg, flush=True)
        out.write(msg + "\n")
        out.flush()

    L("=" * 78)
    L("A7 — GCT invariants of f_chi_P^{(n)}")
    L("=" * 78)

    # --- sanity checks ---
    # In the action f -> f o g of GL_{n^2} on C[Mat_n], the stabilizer of det_n
    # is the image of (GL_n x GL_n) acting via X -> A X B with det(AB) = 1.
    # The map GL_n x GL_n -> GL_{n^2} has 1-dim kernel (lambda I, lambda^{-1} I),
    # so dim image = 2n^2 - 1, and the codim-1 constraint det(AB) = 1 cuts to
    # dim 2n^2 - 2 = 2*4 - 2 = 6 for n=2.  (The textbook value "2n^2 - 1" refers
    # to the abstract group GL_n x GL_n / scalar before identifying its image
    # in GL_{n^2}.)
    L("\n--- Sanity checks (textbook stabilizer dimensions, image in GL_{n^2}) ---")
    d_det2 = sanity_det2()
    L(f"  det_2 (4 vars):  dim Stab Lie = {d_det2}  (expected 6 = 2n^2 - 2 for n=2)")
    d_e2 = sanity_e2_n4()
    L(f"  e_2  (4 vars):   dim Stab Lie = {d_e2}  (expected 6 = dim o(4))")
    d_perm2 = sanity_perm2()
    L(f"  perm_2 (4 vars): dim Stab Lie = {d_perm2}  (expected 6 = same orbit as det_2 for n=2)")
    sanity_ok = (d_det2 == 6 and d_e2 == 6 and d_perm2 == 6)
    L(f"  sanity OK: {sanity_ok}")
    assert sanity_ok, "Sanity check failed; aborting."

    # --- main loop over n ---
    results = {'sanity': {'det_2': d_det2, 'e_2_4': d_e2, 'perm_2': d_perm2}}
    for n in [4, 5, 6]:
        L(f"\n--- n = {n}  (GL_{n}; n^2 = {n*n} params) ---")
        x_syms = sp.symbols(f'x_1:{n + 1}')
        f, monos = chi_p_poly(n, x_syms)
        L(f"  f = {f}")
        L(f"  #monomials = {len(monos)}")

        L("  computing Stab Lie...")
        dim_stab = stab_lie_dim(f, x_syms)
        orbit_dim = n * n - dim_stab
        L(f"  dim Stab Lie  = {dim_stab}")
        L(f"  dim orbit     = {orbit_dim}  (out of n^2 = {n*n})")

        L("  computing perm group...")
        perms = discrete_perm_group(f, x_syms)
        L(f"  |perm group fixing f| = {len(perms)}  (out of n! = {sp.factorial(n)})")
        if len(perms) <= 5:
            L(f"    fixing perms: {perms}")

        sup = chi_p_support(n)
        dt = diag_torus_stab_dim(sup, n)
        L(f"  diag torus stab dim  = {dt}")

        pdsd = pds_dim(f, x_syms)
        L(f"  partial-deriv space dim = {pdsd}  (out of n = {n})")

        sym_rank, pt_rank = hessian_ranks(f, x_syms)
        L(f"  Hessian rank symbolic  = {sym_rank}")
        L(f"  Hessian rank at rand pt = {pt_rank}  (out of n = {n})")

        # Baseline
        n_trials = 30 if n <= 5 else 12
        L(f"  computing baseline (matched random, {n_trials} trials)...")
        baseline = baseline_stab_dims(n, sup, n_trials, seed=42 + n)
        mu = sum(baseline) / n_trials
        sd = (sum((x - mu)**2 for x in baseline) / n_trials)**0.5
        if sd > 0:
            z = (dim_stab - mu) / sd
        else:
            z = float('inf') if dim_stab != mu else 0.0
        mn, mx = min(baseline), max(baseline)
        L(f"  baseline Stab dim: mean={mu:.3f}, std={sd:.3f}, range=[{mn},{mx}]")
        L(f"  z(f_chi_P vs baseline) = {z:+.3f}")

        results[n] = dict(
            f=str(f),
            n_monomials=len(monos),
            dim_stab=dim_stab,
            orbit_dim=orbit_dim,
            perm_group_size=len(perms),
            diag_torus_dim=dt,
            pds_dim=pdsd,
            hessian_sym_rank=sym_rank,
            hessian_pt_rank=pt_rank,
            baseline=baseline,
            baseline_mean=mu,
            baseline_std=sd,
            z=z if z != float('inf') else None,
        )

    # --- Detailed Lie algebra basis at n=4 ---
    L("\n--- Detailed Stab Lie algebra basis at n = 4 ---")
    n = 4
    xs = sp.symbols(f'x_1:{n + 1}')
    f4, _ = chi_p_poly(n, xs)
    dim4, basis4 = stab_lie_dim(f4, xs, return_basis=True)
    L(f"  dim = {dim4}; #basis matrices = {len(basis4)}")
    for k, B in enumerate(basis4):
        L(f"  basis[{k}]:")
        for i in range(n):
            L("    [" + ", ".join(str(B[i, j]) for j in range(n)) + "]")

    # --- Summary ---
    L("\n" + "=" * 78)
    L("SUMMARY (FAL verdicts)")
    L("=" * 78)
    L(f"{'n':>2}|{'#mon':>4}|{'Stab':>5}|{'orbit':>5}|{'#perm':>5}|{'tor':>3}|{'pds':>3}|{'Hsym':>4}|{'Hpt':>3}|{'base μ±σ':>12}|{'z':>6}")
    for n in [4, 5, 6]:
        r = results[n]
        zstr = f"{r['z']:+.2f}" if r['z'] is not None else "  inf"
        L(f"{n:>2}|{r['n_monomials']:>4}|{r['dim_stab']:>5}|{r['orbit_dim']:>5}|"
          f"{r['perm_group_size']:>5}|{r['diag_torus_dim']:>3}|{r['pds_dim']:>3}|"
          f"{r['hessian_sym_rank']:>4}|{r['hessian_pt_rank']:>3}|"
          f"{r['baseline_mean']:.2f}±{r['baseline_std']:.2f}|"
          f"{zstr}")

    L("\n--- FAL verdicts ---")
    for n in [4, 5, 6]:
        r = results[n]
        L(f"n={n}:")
        thr = n * n / 2
        L(f"  FAL-1 (Stab >= n^2/2 = {thr}): Stab={r['dim_stab']}; {'COLLAPSE' if r['dim_stab'] >= thr else 'ok'}")
        nfact = int(sp.factorial(n))
        L(f"  FAL-2 (perm = n! = {nfact}): perm={r['perm_group_size']}; {'symmetric (E)' if r['perm_group_size'] == nfact else 'not symmetric'}")
        L(f"  FAL-3 sanity (PDS=n): PDS={r['pds_dim']}; {'as expected' if r['pds_dim'] == n else 'DEGENERATE'}")
        z = r['z']
        if z is None:
            v4 = "z = inf (zero variance baseline; dev DETECTED)"
        elif abs(z) > 3:
            v4 = f"z = {z:+.2f} : DEVIATES > 3sigma — A-grade direction"
        elif abs(z) > 1:
            v4 = f"z = {z:+.2f} : mild dev"
        else:
            v4 = f"z = {z:+.2f} : matches baseline (B-grade)"
        L(f"  FAL-4 (z vs random): {v4}")
        L(f"  FAL-5 (Hessian rank): pt_rank={r['hessian_pt_rank']}/{n}; "
          f"{'DEFICIT (signal)' if r['hessian_pt_rank'] < n else 'no obstruction'}")

    # --- save JSON ---
    json_path = "/apps/aplikacijos/prime-research/experiments/algebraic/gct_chi_p_orbit/gct_chi_p_orbit_results.json"
    with open(json_path, "w") as fh:
        # serialise
        save_results = {}
        for k, v in results.items():
            if isinstance(v, dict):
                save_results[k] = {kk: (vv if not isinstance(vv, list) else vv) for kk, vv in v.items()}
            else:
                save_results[k] = v
        json.dump(save_results, fh, indent=2, default=str)
    L(f"\nResults JSON saved: {json_path}")
    L(f"Log saved: {log_path}")

    out.close()


if __name__ == "__main__":
    main()
