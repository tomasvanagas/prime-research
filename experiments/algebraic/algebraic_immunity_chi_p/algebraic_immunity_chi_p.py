#!/usr/bin/env python3
"""
Polynomial-method invariants of chi_P (the prime indicator).

ATTACK: §B1 of ATTACK_VECTORS.md (slice rank / Croot-Lev-Pach polynomial method).

Three structurally distinct invariants are computed; each is the
"polynomial method" answer in a different formal sense:

  (1) Algebraic immunity AI(chi_P) over F_2.
      AI(f) = min{ deg(g) : g != 0 and g*f == 0 over F_2^N (or g*(1+f)==0) }.
      The smallest-degree nontrivial annihilator. Standard cryptographic
      invariant (Courtois-Meier 2003). NEVER computed for chi_P in this
      project. If AI grows sub-linearly in N, primes admit a low-degree
      annihilator -- a polynomial-method shortcut.

  (2) F_p multilinear ANF degree deg_p(chi_P).
      View chi_P as a function (F_p)^k -> F_p via base-p expansion of
      n in [0, p^k). Compute its unique multilinear (per-variable
      degree <= p-1) representation in F_p[x_0, ..., x_{k-1}] /
      <x_i^p - x_i>. Degree = max total degree of nonzero coefficient.

  (3) Slice-rank flattening lower bound (Tao 2016).
      For T : X^k -> F, slice_rank(T) >= max_i rank(flatten_i(T)).
      Combined with greedy upper bounds, brackets slice rank.

For each invariant, controls: random Bernoulli matched-density, Liouville
indicator 1[lambda(n)=-1], Mobius squarefree indicator 1[mu(n)!=0],
"all-zero" baseline.

Outcome interpretation:
  AI(chi_P) sub-linear -> A-grade (polynomial-method shortcut).
  AI(chi_P) = Theta(N) matching random -> B-grade negative-shape edge.
  AI(chi_P) below random by margin but still Theta(N) -> B-grade
    structural identification (primes' density-dependent AI dip).

Cross-domain refs:
  - Croot-Lev-Pach 2017 arXiv:1605.01506
  - Ellenberg-Gijswijt 2017 Annals 185, 339
  - Tao 2016 blog "A symmetric formulation of the Croot-Lev-Pach
    /Ellenberg-Gijswijt capset bound"
  - Courtois-Meier 2003 "Algebraic attacks on stream ciphers with
    linear feedback" Eurocrypt LNCS 2656.
"""

import time
import json
from pathlib import Path
import numpy as np
from sympy import isprime, factorint, mobius
from sympy.ntheory.factor_ import core


# ---------------------------------------------------------------
# 1. ARITHMETIC FUNCTIONS
# ---------------------------------------------------------------

def chi_P(n):
    return 1 if (n >= 2 and isprime(n)) else 0


def liouville_pos(n):
    """1 if lambda(n) = +1 (n has even number of prime factors w/ multiplicity)."""
    if n <= 0:
        return 0
    if n == 1:
        return 1
    f = factorint(n)
    omega = sum(f.values())
    return 1 if omega % 2 == 0 else 0


def mobius_nonzero(n):
    """1 iff mu(n) != 0, i.e., n is squarefree."""
    if n <= 0:
        return 0
    return 1 if mobius(n) != 0 else 0


def build_truth_table(func, N):
    """Compute truth table of func: [0, 2^N) -> {0,1}."""
    return np.fromiter((func(n) for n in range(2 ** N)), dtype=np.int8, count=2 ** N)


def build_truth_table_modp(func, p, k):
    """Compute truth table of func: [0, p^k) -> {0,1}."""
    return np.fromiter((func(n) for n in range(p ** k)), dtype=np.int8, count=p ** k)


# ---------------------------------------------------------------
# 2. ALGEBRAIC IMMUNITY OVER F_2
# ---------------------------------------------------------------

def gauss_elim_F2(M):
    """In-place row-reduction over F_2. Returns rank."""
    rows, cols = M.shape
    pivot_row = 0
    for col in range(cols):
        if pivot_row >= rows:
            break
        idx = None
        for r in range(pivot_row, rows):
            if M[r, col] == 1:
                idx = r
                break
        if idx is None:
            continue
        if idx != pivot_row:
            M[[pivot_row, idx]] = M[[idx, pivot_row]]
        for r in range(rows):
            if r != pivot_row and M[r, col] == 1:
                M[r] ^= M[pivot_row]
        pivot_row += 1
    return pivot_row


def monomials_up_to_degree(N, d):
    """Enumerate index sets S subset of {0,...,N-1} with |S| <= d.
    Returns list of bitmasks (each integer's bits indicate the set).
    """
    out = [0]
    if d == 0:
        return out
    # Subsets of size 1..d
    for size in range(1, d + 1):
        # generate all C(N, size) subsets
        idxs = [0] * size
        idxs[0] = 0
        # iterative combination generation
        # simpler: use numpy
        from itertools import combinations
        for combo in combinations(range(N), size):
            mask = 0
            for c in combo:
                mask |= (1 << c)
            out.append(mask)
    return out


def evaluate_monomials(N, masks):
    """Build a matrix V of shape (2^N, len(masks)) where
    V[x, j] = product over bits set in masks[j] of x's bits.
    Computed efficiently in F_2 (just AND match).
    """
    M = len(masks)
    rng = 2 ** N
    V = np.zeros((rng, M), dtype=np.int8)
    # V[x, j] = 1 iff all bits in masks[j] are set in x  (i.e. x & masks[j] == masks[j])
    # This is the multilinear monomial evaluation: prod over i in S of x_i = AND of those bits
    xs = np.arange(rng, dtype=np.int64)
    for j, m in enumerate(masks):
        V[:, j] = ((xs & m) == m).astype(np.int8)
    return V


def algebraic_immunity_F2(f_tt, N, max_d=None, return_certificate=False):
    """
    Algebraic immunity of f : F_2^N -> F_2 with truth table f_tt.

    AI(f) = min(d : exists g != 0, deg(g) <= d, with g*f == 0  OR  g*(1+f) == 0).

    Returns AI; optionally a certificate (annihilator).
    """
    if max_d is None:
        max_d = N
    rng = 2 ** N
    # Indices where f = 1 and f = 0
    one_idx = np.nonzero(f_tt)[0]
    zero_idx = np.nonzero(1 - f_tt)[0]

    # Edge cases
    if len(one_idx) == 0:
        return 0, "f is identically 0; AI = 0 (g = 1 annihilates f)."
    if len(zero_idx) == 0:
        return 0, "f is identically 1; AI = 0 (g = 1 annihilates 1-f)."

    for d in range(0, max_d + 1):
        masks = monomials_up_to_degree(N, d)
        V = evaluate_monomials(N, masks)
        # Want g : F_2^N -> F_2 with deg(g) <= d, expressed as
        # g(x) = sum over masks j of c_j * V[x, j].
        # Annihilator of f: g(x) = 0 for all x with f(x) = 1.
        # That is: V[one_idx, :] @ c == 0 over F_2.
        # Find non-trivial c.
        for sub_idx, label in [(one_idx, "f"), (zero_idx, "1+f")]:
            A = V[sub_idx].copy()
            # Solve A @ c == 0 over F_2 for non-trivial c.
            # Equivalent: rank(A) < |masks|. Get nullspace.
            # Need columns >= rows for non-trivial null possibility for "small" d.
            # But even rows>=cols, can have non-trivial null if A is rank-deficient.
            #
            # Approach: row-reduce A (over F_2) and check nullspace dim.
            n_rows, n_cols = A.shape
            B = A.copy()
            r = gauss_elim_F2(B)
            null_dim = n_cols - r
            if null_dim > 0:
                # Found a nontrivial annihilator at degree d.
                cert = (label, d, n_cols, n_rows, r)
                return d, cert if return_certificate else f"deg-{d} ann of {label}"
    return max_d + 1, "AI > max_d"


# ---------------------------------------------------------------
# 3. F_p MULTILINEAR ANF DEGREE
# ---------------------------------------------------------------

def base_p_digits(n, p, k):
    """Return [d_0, d_1, ..., d_{k-1}] with n = sum d_i * p^i."""
    out = [0] * k
    for i in range(k):
        out[i] = n % p
        n //= p
    return out


def multilinear_anf_Fp(f_tt, p, k):
    """
    Compute the unique representation of f : F_p^k -> F_p as
    f(x_0, ..., x_{k-1}) = sum over (e_0,...,e_{k-1}) in [0,p-1]^k of
                            c_{e_0,...,e_{k-1}} * prod x_i^{e_i}
    in F_p[x_0,...,x_{k-1}] / <x_i^p - x_i>.

    Returns the coefficient tensor c of shape p^k (flattened) over F_p,
    and the maximum total degree of any nonzero coefficient.
    """
    # Build the evaluation matrix M of shape (p^k, p^k):
    #   M[(d_0,...,d_{k-1}), (e_0,...,e_{k-1})] = prod d_i^{e_i} mod p
    # Then f_tt = M @ c. Solve for c by inverting M (or solve linear system).
    #
    # For p prime, the per-variable Vandermonde V[d, e] = d^e mod p has
    # determinant prod_{0 <= a < b < p} (a - b) (Vandermonde) which is
    # nonzero in F_p, so V is invertible. Then M = V tensor V tensor ... tensor V (k factors).
    # Inverse is V^{-1} tensor ... tensor V^{-1}.

    V = np.zeros((p, p), dtype=np.int64)
    for d in range(p):
        for e in range(p):
            if d == 0 and e == 0:
                V[d, e] = 1  # 0^0 = 1 convention
            else:
                V[d, e] = pow(d, e, p)
    # Compute V^{-1} mod p
    V_inv = mod_matrix_inverse_Fp(V, p)
    # Apply V_inv to each axis (tensor product)
    # Reshape f_tt to shape (p,)*k, then apply V_inv along each axis.
    arr = np.array(f_tt, dtype=np.int64).reshape((p,) * k) % p
    for axis in range(k):
        arr = np.moveaxis(arr, axis, 0)
        shape = arr.shape
        flat = arr.reshape(p, -1)
        flat = (V_inv @ flat) % p
        arr = flat.reshape(shape)
        arr = np.moveaxis(arr, 0, axis)
    coefs = arr % p

    # Compute max total degree of nonzero coefficient.
    max_deg = -1
    for idx in np.argwhere(coefs != 0):
        deg = int(idx.sum())
        if deg > max_deg:
            max_deg = deg

    # Also: number of nonzero coefficients (sparsity).
    nnz = int(np.sum(coefs != 0))
    return coefs, max_deg, nnz


def mod_matrix_inverse_Fp(M, p):
    """Inverse of n x n matrix M over F_p (p prime). Gauss-Jordan."""
    n = M.shape[0]
    A = np.hstack([M.copy() % p, np.eye(n, dtype=np.int64)]) % p
    for col in range(n):
        # Find pivot
        pivot = None
        for r in range(col, n):
            if A[r, col] % p != 0:
                pivot = r
                break
        if pivot is None:
            raise ValueError("Singular matrix mod p")
        if pivot != col:
            A[[col, pivot]] = A[[pivot, col]]
        # Scale pivot row
        inv = pow(int(A[col, col]) % p, p - 2, p)  # Fermat's little theorem inverse
        A[col] = (A[col] * inv) % p
        # Eliminate other rows
        for r in range(n):
            if r != col and A[r, col] % p != 0:
                A[r] = (A[r] - A[r, col] * A[col]) % p
    return A[:, n:] % p


# ---------------------------------------------------------------
# 4. SLICE-RANK BOUNDS (greedy upper + flattening lower)
# ---------------------------------------------------------------

def flatten_rank_F2(T_int, axis):
    """Rank over F_2 of T flattened with given axis as 'row' index."""
    M = np.moveaxis(T_int, axis, 0)
    flat = M.reshape(M.shape[0], -1).astype(np.int8)
    return gauss_elim_F2(flat.copy())


def slice_rank_lower_bound_F2(T_int):
    """Tao 2016 inequality: slice_rank(T) >= max_i rank(flatten_i(T))."""
    return max(flatten_rank_F2(T_int, axis) for axis in range(T_int.ndim))


def slice_rank_upper_bound_F2_greedy(T_int, max_iter=10000):
    """Greedy upper bound: peel off the heaviest slice until residual = 0."""
    if T_int.ndim < 2:
        return 1 if np.any(T_int) else 0
    R = T_int.copy().astype(np.int8)
    n = 0
    while np.any(R) and n < max_iter:
        best_w = -1
        best = None
        for axis in range(R.ndim):
            for j in range(R.shape[axis]):
                # Take the slice fixing axis at j.
                idx = [slice(None)] * R.ndim
                idx[axis] = j
                sl = R[tuple(idx)]
                w = int(np.sum(sl != 0))
                if w > best_w:
                    best_w = w
                    best = (axis, j, sl.copy())
        if best is None or best_w == 0:
            break
        axis, j, sl = best
        idx = [slice(None)] * R.ndim
        idx[axis] = j
        # Subtract over F_2 (XOR)
        R[tuple(idx)] ^= sl
        n += 1
    return n


# ---------------------------------------------------------------
# 5. RANDOM CONTROLS
# ---------------------------------------------------------------

def random_truth_table_matched_density(N, density, rng):
    """Random Bernoulli truth table at given density."""
    return (rng.random(2 ** N) < density).astype(np.int8)


# ---------------------------------------------------------------
# 6. EXPERIMENT DRIVERS
# ---------------------------------------------------------------

def run_AI_F2(N_range, n_random_seeds=10, max_d_factor=1.0):
    """
    Run AI(.) over F_2 for chi_P, Liouville_pos, Mobius_nonzero,
    matched-density random Bernoulli, all-zero / all-one baselines.
    """
    print(f"\n{'=' * 72}")
    print("ALGEBRAIC IMMUNITY OVER F_2")
    print(f"{'=' * 72}")
    rows = []
    for N in N_range:
        print(f"\n--- N = {N} (truth tables of length {2**N}) ---")

        # chi_P
        t0 = time.time()
        tt_chi_p = build_truth_table(chi_P, N)
        tt_lam = build_truth_table(liouville_pos, N)
        tt_mu = build_truth_table(mobius_nonzero, N)
        rho_chi = float(np.mean(tt_chi_p))
        rho_lam = float(np.mean(tt_lam))
        rho_mu = float(np.mean(tt_mu))
        print(f"  density chi_P:        {rho_chi:.4f} ({int(np.sum(tt_chi_p))}/{2**N})")
        print(f"  density Liouville+:   {rho_lam:.4f}")
        print(f"  density Mobius!=0:    {rho_mu:.4f}")
        print(f"  built truth tables in {time.time()-t0:.2f}s")

        max_d = max(2, int(round(N * max_d_factor)))

        for label, tt, rho in [
            ("chi_P", tt_chi_p, rho_chi),
            ("Liouville+", tt_lam, rho_lam),
            ("Mobius!=0", tt_mu, rho_mu),
        ]:
            t0 = time.time()
            ai, cert = algebraic_immunity_F2(tt, N, max_d=max_d)
            t = time.time() - t0
            print(f"  AI({label}, N={N}) = {ai}  ({cert})  [{t:.2f}s]")
            rows.append({"N": N, "function": label, "density": rho, "AI": ai, "cert": str(cert), "time_s": t})

        # Random Bernoulli at chi_P density
        ais_random = []
        for seed in range(n_random_seeds):
            rng = np.random.default_rng(1000 + N * 31 + seed)
            tt_rand = random_truth_table_matched_density(N, rho_chi, rng)
            t0 = time.time()
            ai, _ = algebraic_immunity_F2(tt_rand, N, max_d=max_d)
            t = time.time() - t0
            ais_random.append(ai)
            if seed < 3:
                print(f"  AI(random_seed{seed}, N={N}, rho={rho_chi:.4f}) = {ai}  [{t:.2f}s]")
        ais_random = np.array(ais_random)
        print(f"  AI(random, N={N}, rho={rho_chi:.4f}): mean={ais_random.mean():.2f}, "
              f"std={ais_random.std():.2f}, min={ais_random.min()}, max={ais_random.max()}")
        rows.append({"N": N, "function": "random",
                     "density": rho_chi,
                     "AI_mean": float(ais_random.mean()),
                     "AI_std": float(ais_random.std()),
                     "AI_min": int(ais_random.min()),
                     "AI_max": int(ais_random.max()),
                     "n_seeds": n_random_seeds})

    return rows


def run_Fp_anf_degree(p_k_list):
    """For each (p, k), compute multilinear ANF degree of chi_P, controls."""
    print(f"\n{'=' * 72}")
    print("F_p MULTILINEAR ANF DEGREE")
    print(f"{'=' * 72}")
    rows = []
    for (p, k) in p_k_list:
        N = p ** k
        print(f"\n--- p = {p}, k = {k} (truth table size = {N}) ---")
        # Build truth tables; for non-prime p, skip (we want p prime).
        if not isprime(p):
            print(f"  skip: p = {p} is not prime")
            continue

        tt_chi = build_truth_table_modp(chi_P, p, k)
        tt_lam = build_truth_table_modp(liouville_pos, p, k)
        tt_mu  = build_truth_table_modp(mobius_nonzero, p, k)
        rho = float(np.mean(tt_chi))
        print(f"  density chi_P: {rho:.4f}")

        for label, tt in [("chi_P", tt_chi), ("Liouville+", tt_lam), ("Mobius!=0", tt_mu)]:
            t0 = time.time()
            coefs, deg, nnz = multilinear_anf_Fp(tt, p, k)
            t = time.time() - t0
            max_deg = (p - 1) * k
            print(f"  deg_{p}({label}, k={k}) = {deg} (max possible {max_deg}); nnz coeffs = {nnz} / {p**k}  [{t:.2f}s]")
            rows.append({"p": p, "k": k, "function": label,
                         "deg": int(deg), "max_deg": int(max_deg),
                         "nnz_coefs": int(nnz), "total_coefs": int(p ** k),
                         "density": float(np.mean(tt)), "time_s": t})

        # Random Bernoulli at chi_P density
        n_rand = 5
        rng = np.random.default_rng(100 + p * k)
        degs_rand = []
        nnz_rand = []
        for s in range(n_rand):
            tt_rand = (rng.random(N) < rho).astype(np.int64)
            _, deg, nnz = multilinear_anf_Fp(tt_rand, p, k)
            degs_rand.append(deg)
            nnz_rand.append(nnz)
        print(f"  random rho={rho:.4f}: deg mean={np.mean(degs_rand):.2f} "
              f"std={np.std(degs_rand):.2f} min={min(degs_rand)} max={max(degs_rand)}; "
              f"nnz mean={np.mean(nnz_rand):.1f}")
        rows.append({"p": p, "k": k, "function": "random",
                     "deg_mean": float(np.mean(degs_rand)),
                     "deg_std": float(np.std(degs_rand)),
                     "deg_min": int(min(degs_rand)),
                     "deg_max": int(max(degs_rand)),
                     "nnz_mean": float(np.mean(nnz_rand)),
                     "n_seeds": n_rand,
                     "density": rho})
    return rows


def run_slice_rank_brackets(p, k_max=4, mode="chi_P"):
    """Brackets on slice rank for chi_P viewed as a tensor on (Z/p)^k.

    For each k from 2 to k_max:
      - Lower bound: max flattening rank over F_2
      - Upper bound: greedy decomposition (over F_2)

    Note: for p=2 this is same as the bit-partition tensor.
    """
    print(f"\n{'=' * 72}")
    print(f"SLICE RANK BRACKETS over F_2 (function = {mode}, base p = {p})")
    print(f"{'=' * 72}")
    rows = []
    for k in range(2, k_max + 1):
        N = p ** k
        # Build truth table (on integers [0, N))
        if mode == "chi_P":
            tt = build_truth_table_modp(chi_P, p, k)
        elif mode == "Liouville+":
            tt = build_truth_table_modp(liouville_pos, p, k)
        elif mode == "Mobius!=0":
            tt = build_truth_table_modp(mobius_nonzero, p, k)
        else:
            raise ValueError(mode)
        rho = float(np.mean(tt))

        # Reshape into k-tensor of shape (p, p, ..., p)
        T = np.array(tt, dtype=np.int8).reshape((p,) * k)
        print(f"\n--- k = {k}, shape {T.shape}, density {rho:.4f}, |X|^k = {N} ---")

        # Compute rank lower bound (over F_2)
        t0 = time.time()
        lb = slice_rank_lower_bound_F2(T)
        t_lb = time.time() - t0
        print(f"  flatten LB (F_2): {lb}  [{t_lb:.3f}s]")

        # Compute upper bound (greedy)
        t0 = time.time()
        ub = slice_rank_upper_bound_F2_greedy(T)
        t_ub = time.time() - t0
        print(f"  greedy UB (F_2):  {ub}  [{t_ub:.3f}s]")

        # Random control
        rng = np.random.default_rng(2024 + p * 7 + k)
        rand_lbs, rand_ubs = [], []
        for s in range(3):
            T_rand = (rng.random(N) < rho).astype(np.int8).reshape((p,) * k)
            rand_lbs.append(slice_rank_lower_bound_F2(T_rand))
            rand_ubs.append(slice_rank_upper_bound_F2_greedy(T_rand))
        print(f"  random LB: {rand_lbs}, mean={np.mean(rand_lbs):.1f}")
        print(f"  random UB: {rand_ubs}, mean={np.mean(rand_ubs):.1f}")
        rows.append({"p": p, "k": k, "function": mode, "density": rho,
                     "slice_LB": int(lb), "slice_UB": int(ub),
                     "random_LB_mean": float(np.mean(rand_lbs)),
                     "random_UB_mean": float(np.mean(rand_ubs))})
    return rows


def main():
    print("ATTACK: §B1  (slice rank / Croot-Lev-Pach polynomial method on chi_P)")
    print("Goal: distinguish chi_P from a random matched-density function on")
    print("polynomial-method invariants (algebraic immunity, F_p ANF degree,")
    print("slice rank brackets).")

    out_dir = Path(__file__).parent
    results = {}

    # 1. Algebraic immunity F_2  -- range tuned for tractability.
    # AI computation at N=14 needs 2^14 x sum(C(14,i), i<=7) = 16384 x 9908 binary matrix.
    # Tractable; ~10 seconds per AI call with sympy + numpy.
    AI_rows = run_AI_F2(N_range=range(4, 14), n_random_seeds=8)
    results["AI_F2"] = AI_rows

    # 2. F_p multilinear ANF degree
    p_k = [
        (3, 2), (3, 3), (3, 4), (3, 5),
        (5, 2), (5, 3),
        (7, 2), (7, 3),
        (11, 2),
    ]
    Fp_rows = run_Fp_anf_degree(p_k)
    results["Fp_anf_degree"] = Fp_rows

    # 3. Slice rank brackets, p=2 (bit partition), p=3
    SR_rows = []
    SR_rows += run_slice_rank_brackets(p=2, k_max=10, mode="chi_P")
    SR_rows += run_slice_rank_brackets(p=3, k_max=4, mode="chi_P")
    SR_rows += run_slice_rank_brackets(p=2, k_max=10, mode="Liouville+")
    results["slice_rank_brackets"] = SR_rows

    # Save JSON
    with open(out_dir / "algebraic_immunity_chi_p_data.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nSaved data to {out_dir / 'algebraic_immunity_chi_p_data.json'}")
    print("Done.")


if __name__ == "__main__":
    main()
