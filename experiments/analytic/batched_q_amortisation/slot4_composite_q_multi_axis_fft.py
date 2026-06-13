"""
Thread 6, slot 4: composite-q multi-axis FFT for Dirichlet L-eval +
end-to-end zero finder.

Mission (from .commit_state Thread 6 slot 4 PRIORITY (b) recommendation
in S229):

  Extend slot-2's cyclic-group FFT primitive and slot-3's end-to-end
  Dirichlet-L zero finder to *composite* squarefree odd q. Composite q
  is the practical AP regime (q = 30, 60, 105, 1001 are common conductors)
  but slot 3 only handled prime q because it relied on the cyclic
  structure of (Z/qZ)* via primitive_root().

Mathematics:

  For squarefree odd q = p_1 * p_2 * ... * p_k (distinct odd primes),
  the multiplicative group decomposes via CRT:

    (Z/qZ)*  ≅  (Z/p_1 Z)* × (Z/p_2 Z)* × ... × (Z/p_k Z)*
            ≅  Z/(p_1-1) ⊕ Z/(p_2-1) ⊕ ... ⊕ Z/(p_k-1)

  via primitive roots g_i mod p_i. Each n ∈ (Z/qZ)* corresponds to a
  tuple (k_1, ..., k_k) with n ≡ g_i^{k_i} (mod p_i).

  Characters factor: χ_{j_1,...,j_k}(n) = Π_i ω_i^{j_i k_i} where
  ω_i = exp(2πi / (p_i-1)).

  Multi-axis FFT identity: aggregate W indexed by (k_1,...,k_k) ∈
  ⊕_i Z/(p_i-1)Z gives L-values for ALL φ(q) characters via

    L[j_1,...,j_k]
        = Σ_{(k_1,...,k_k)} Π_i ω_i^{j_i k_i} · W[k_1,...,k_k]
        = phi · ifftn(W, axes=(0,...,k-1))[j_1,...,j_k]

  where phi = Π_i (p_i - 1) = φ(q) is the total scale factor.

What this slot does:

  1. CRT-based multi-axis FFT primitive `l_eval_fft_multiaxis` for
     squarefree odd q with k ≥ 2 prime factors.
  2. Cross-check vs `l_eval_direct_composite` (per-character matmul) at
     q ∈ {15, 35, 105, 1001} to verify the multi-axis identity.
  3. End-to-end zero finder for composite q via the same Hardy Z +
     sign-change machinery as slot 3.
  4. End-to-end π(x; q, a) timing at q ∈ {15, 35, 105, 1001} compared
     to slot-3's prime-q baseline.

What would falsify slot 4 progress:

  (i)   Multi-axis FFT does not produce correctly-aggregated L-values
        (cross-method check fails at small q). Closure mode E.
  (ii)  FFT primitive is FASTER on composite q than on adjacent prime q
        (e.g., q=1001 vs q=1009) by an asymptotic factor. Partial-positive.
  (iii) FFT primitive is the SAME or SLOWER on composite q (BLAS already
        saturates the FLOP rate for the small per-axis transforms).
        Closure mode I — extends slot-3's regime but no asymptotic win.

Cross-domain ingredient:

  CRT-based multi-axis FFT over Z/(p_1-1) ⊕ ... ⊕ Z/(p_k-1). This
  generalises the slot-2 cyclic primitive (Bober-Hiary 2017) and is
  equivalent to `numpy.fft.ifftn` on a multi-dimensional aggregate.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
import time
from typing import Dict, List, Tuple

import numpy as np
from scipy.special import loggamma

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from per_q_a_explicit_formula_profile import euler_phi, primitive_root, get_chi0_zeros  # noqa: E402
from slot2_afe_shared_l_eval import compute_complex_pow_all_n  # noqa: E402


# =====================================================================
# CRT decomposition of (Z/qZ)* for squarefree odd q
# =====================================================================


def factor_squarefree(q: int) -> List[int]:
    """Return prime factors of q (only if q is squarefree). Raise otherwise."""
    factors: List[int] = []
    n = q
    p = 2
    while p * p <= n:
        if n % p == 0:
            count = 0
            while n % p == 0:
                count += 1
                n //= p
            if count > 1:
                raise ValueError(f"q={q} is not squarefree (p={p} appears {count} times)")
            factors.append(p)
        p += 1
    if n > 1:
        factors.append(n)
    return factors


def crt_decomp(q: int) -> Tuple[List[int], List[int], List[int]]:
    """For squarefree q, return (primes, group_orders, primitive_roots).

    primes: list of distinct primes p_i with q = Π p_i.
    group_orders: list of (p_i - 1) (the order of (Z/p_iZ)*).
    primitive_roots: list of g_i with g_i a primitive root mod p_i.

    Requires q odd (no p=2). For p=2 the group is trivial and is omitted
    (it contributes a factor of 1 to the product).
    """
    primes = factor_squarefree(q)
    if 2 in primes:
        raise ValueError(f"q={q} contains p=2; multi-axis FFT for q with p=2 not implemented")
    group_orders = [p - 1 for p in primes]
    prim_roots = [primitive_root(p) for p in primes]
    if any(g is None for g in prim_roots):
        bad = [p for p, g in zip(primes, prim_roots) if g is None]
        raise ValueError(f"primitive root not found for primes {bad}")
    return primes, group_orders, prim_roots


def crt_residue_to_tuple(n: int, primes: List[int], prim_roots: List[int]) -> Tuple[int, ...]:
    """Map n ∈ (Z/qZ)* to tuple (k_1, ..., k_K) where n ≡ g_i^{k_i} (mod p_i)."""
    out: List[int] = []
    for p, g in zip(primes, prim_roots):
        r = n % p
        if r == 0:
            raise ValueError(f"n={n} not coprime to p={p}")
        # discrete log of r base g mod p
        cur = 1
        for k in range(p - 1):
            if cur == r:
                out.append(k)
                break
            cur = (cur * g) % p
        else:
            raise ValueError(f"discrete log failed for r={r} mod p={p}")
    return tuple(out)


def build_crt_log_table(q: int, primes: List[int], prim_roots: List[int]
                         ) -> Dict[int, Tuple[int, ...]]:
    """Build a dict mapping each residue r ∈ [1, q-1] coprime to q to its
    CRT log-tuple (k_1, ..., k_K)."""
    table: Dict[int, Tuple[int, ...]] = {}
    for r in range(1, q):
        if math.gcd(r, q) == 1:
            table[r] = crt_residue_to_tuple(r, primes, prim_roots)
    return table


# =====================================================================
# Multi-axis FFT L-eval primitive
# =====================================================================


def l_eval_multiaxis_fft(q: int, t_grid: np.ndarray, N: int) -> Tuple[np.ndarray, List[int]]:
    """Compute L_main(½+it, χ) = Σ_{n ≤ N, gcd(n,q)=1} χ(n)/n^{½+it} for
    ALL χ mod q at all t in t_grid via multi-axis FFT.

    Returns:
      L_flat: (phi, n_t) complex matrix, characters indexed in row-major
              tuple order (k_1, k_2, ..., k_K) → flat = k_1*Π_{i>1} sz_i + ...
      group_orders: [p_1-1, p_2-1, ..., p_K-1] used to decode.
    """
    primes, group_orders, prim_roots = crt_decomp(q)
    log_table = build_crt_log_table(q, primes, prim_roots)
    phi_total = 1
    for p_i in group_orders:
        phi_total *= p_i
    n_t = len(t_grid)

    cp_all = compute_complex_pow_all_n(t_grid, N)  # (N, n_t)

    # Aggregate W with shape (p_1-1, ..., p_K-1, n_t)
    W_shape = tuple(group_orders) + (n_t,)
    W = np.zeros(W_shape, dtype=complex)

    for n_idx in range(N):
        r = (n_idx + 1) % q
        if r in log_table:
            tup = log_table[r]
            # Index with full slice over t-axis
            slc = tup + (slice(None),)
            W[slc] += cp_all[n_idx]

    # Multi-axis ifft on the first len(group_orders) axes
    axes = tuple(range(len(group_orders)))
    L = phi_total * np.fft.ifftn(W, axes=axes)
    L_flat = L.reshape(phi_total, n_t)
    return L_flat, group_orders


def l_eval_direct_composite(q: int, t_grid: np.ndarray, N: int
                              ) -> Tuple[np.ndarray, List[int], np.ndarray]:
    """Direct per-character L-evaluation for composite q.

    Builds the φ(q) × N_co character matrix manually from CRT character
    indexing, then matmuls. Used as ground truth for cross-checking the
    multi-axis FFT primitive.

    Returns (L: shape (phi, n_t), group_orders, char_indices: shape (phi, K)
    where char_indices[flat] = (j_1, ..., j_K)).
    """
    primes, group_orders, prim_roots = crt_decomp(q)
    log_table = build_crt_log_table(q, primes, prim_roots)
    phi_total = 1
    for p_i in group_orders:
        phi_total *= p_i

    # Build coprime list (1..q-1, gcd=1) and residues for n ≤ N
    n_arr = np.arange(1, N + 1)
    n_residues = (n_arr % q).astype(int)
    coprime_mask = np.array([(int(r) != 0 and r in log_table) for r in n_residues])
    n_co = n_arr[coprime_mask]
    n_co_residues = (n_co % q).astype(int).tolist()

    log_n = np.log(n_co.astype(np.float64))
    inv_sqrt_n = 1.0 / np.sqrt(n_co.astype(np.float64))
    phase = -np.outer(log_n, t_grid)
    cp = inv_sqrt_n[:, None] * np.exp(1j * phase)  # (N_co, n_t)

    # Build char indices: char j ↔ (j_1, ..., j_K) by row-major flatten
    char_indices = np.zeros((phi_total, len(group_orders)), dtype=np.int64)
    flat = 0
    multi_idx = [0] * len(group_orders)
    for flat in range(phi_total):
        rem = flat
        for ax in range(len(group_orders) - 1, -1, -1):
            multi_idx[ax] = rem % group_orders[ax]
            rem //= group_orders[ax]
        char_indices[flat] = multi_idx

    # Build (phi, N_co) char_table[j_flat, ci] = Π_i ω_i^{j_i · k_i(coprimes[ci])}
    # via tuple-based vectorised exponentiation
    K = len(group_orders)
    chi_table = np.ones((phi_total, len(n_co)), dtype=complex)
    for i in range(K):
        omega_i = np.exp(2j * np.pi / group_orders[i])
        # k_i for each coprime residue
        k_i_per_res = np.array([log_table[r][i] for r in n_co_residues], dtype=np.int64)
        # j_i for each character
        j_i_per_char = char_indices[:, i]
        # Phase: omega_i ^ (j_i * k_i mod (p_i-1))
        phase_idx = (j_i_per_char[:, None] * k_i_per_res[None, :]) % group_orders[i]
        chi_table *= np.exp(2j * np.pi * phase_idx / group_orders[i])

    L = chi_table @ cp  # (phi, n_t)
    return L, group_orders, char_indices


# =====================================================================
# Cross-method correctness
# =====================================================================


def check_multi_axis_correctness(q: int, t_grid_size: int = 5,
                                   t_max: float = 30.0, N: int = 80,
                                   verbose: bool = True) -> Tuple[float, float]:
    """Compare multi-axis FFT vs direct per-character matmul at composite q.
    Returns (max abs error, max relative error)."""
    t_grid = np.linspace(2.0, t_max, t_grid_size)

    L_dir, _, _ = l_eval_direct_composite(q, t_grid, N)
    L_fft, _ = l_eval_multiaxis_fft(q, t_grid, N)

    err = float(np.max(np.abs(L_dir - L_fft)))
    max_abs = float(np.max(np.abs(L_dir)))
    rel = err / max_abs if max_abs > 0 else 0.0
    if verbose:
        print(f"  q={q} (factors {factor_squarefree(q)}, phi={euler_phi(q)})  "
              f"N={N} n_t={t_grid_size}  "
              f"max|L| = {max_abs:.3e}  err = {err:.3e}  rel = {rel:.3e}")
    return err, rel


# =====================================================================
# Full AFE for composite q (main + reflected term)
# =====================================================================


def gauss_sum_composite(q: int, char_indices: np.ndarray,
                         primes: List[int], prim_roots: List[int],
                         log_table: Dict[int, Tuple[int, ...]]) -> np.ndarray:
    """Compute τ(χ) = Σ_{a ∈ (Z/qZ)*} χ(a) e^{2πi a/q} for all χ mod q.

    For composite q with characters factoring via CRT, the Gauss sum has
    a known product formula when χ is primitive:

      τ(χ) = Π_i χ_inv_i(q/p_i) · τ(χ_i)

    But for the slot 4 baseline we just compute τ directly via the
    multi-axis FFT identity:

      τ(χ_j) = Σ_a χ_j(a) e^{2πi a/q}
             = phi · ifftn(e_a)[j]

    where e_a is indexed by the CRT-tuple of a.

    Returns shape (phi,) complex.
    """
    primes_l, group_orders, prim_roots_l = crt_decomp(q)
    phi_total = 1
    for p_i in group_orders:
        phi_total *= p_i

    e_array = np.zeros(tuple(group_orders), dtype=complex)
    for r, tup in log_table.items():
        e_array[tup] = np.exp(2j * np.pi * r / q)

    axes = tuple(range(len(group_orders)))
    tau_array = phi_total * np.fft.ifftn(e_array, axes=axes)
    return tau_array.reshape(phi_total)


def chi_at_minus_1_composite(q: int, primes: List[int], prim_roots: List[int],
                               char_indices: np.ndarray, group_orders: List[int]
                               ) -> np.ndarray:
    """Compute χ(-1) ∈ {±1} for all χ mod q (composite case).

    For each odd prime p_i, ω_i = exp(2πi/(p_i-1)) and ω_i^{(p_i-1)/2} = -1.
    Since -1 ≡ g_i^{(p_i-1)/2} (mod p_i), we have

        χ_i(-1) = ω_i^{j_i · (p_i-1)/2} = (-1)^{j_i}.

    Therefore χ(-1) = Π_i (-1)^{j_i} = (-1)^{Σ_i j_i}.
    """
    exp_sum = (char_indices.sum(axis=1)) % 2  # (phi,)
    return np.where(exp_sum == 0, 1.0, -1.0)


def is_primitive_character_composite(q: int, char_indices: np.ndarray,
                                       primes: List[int], group_orders: List[int]
                                       ) -> np.ndarray:
    """A character χ = Π_i χ_i is primitive mod q (q squarefree) iff every
    χ_i is non-principal (χ_i ≠ trivial). For squarefree q this means
    j_i ≠ 0 for ALL i."""
    return np.all(char_indices != 0, axis=1)


def compute_full_l_via_afe_composite(
    q: int, t_grid: np.ndarray, N: int,
    method: str = "fft",
) -> Dict[str, np.ndarray]:
    """Full AFE: L(½+it, χ) = M(t,χ) + W_χ · (q/π)^{-it} · ratio_gamma(t, a) · conj(M(t,χ))
    for every primitive χ mod q. Composite q version.

    For non-primitive χ (some j_i = 0), the L-function L(s, χ) factors
    through L(s, χ') where χ' is the primitive character inducing χ.
    For slot 4 we restrict zero-finding to primitive χ; non-primitive χ
    can be handled by the corresponding lower-conductor primitive
    L-function with extra Euler factors at the missing primes.

    method ∈ {"fft", "direct"} selects the main-term primitive.

    Returns dict with:
      L: (phi, n_t)        full AFE values for ALL characters (incl. non-primitive);
                            for non-primitive χ this just adds main+reflected as
                            usual but the formula isn't exact (only used as test).
      W_chi: (phi,)         root numbers
      parity: (phi,)        0 if χ(-1)=1 else 1
      char_indices: (phi, K)
      is_primitive: (phi,)
    """
    primes, group_orders, prim_roots = crt_decomp(q)
    log_table = build_crt_log_table(q, primes, prim_roots)

    # Compute L_main with chosen method
    if method == "fft":
        M, _ = l_eval_multiaxis_fft(q, t_grid, N)
    elif method == "direct":
        M, _, _ = l_eval_direct_composite(q, t_grid, N)
    else:
        raise ValueError(f"unknown method {method!r}")

    # Get char_indices via direct call (cheap)
    _, _, char_indices = l_eval_direct_composite(q, t_grid[:1], 1)
    phi_total = M.shape[0]

    # Gauss sums + root numbers
    tau = gauss_sum_composite(q, char_indices, primes, prim_roots, log_table)
    chi_minus = chi_at_minus_1_composite(q, primes, prim_roots, char_indices, group_orders)
    parity = np.where(chi_minus > 0, 0, 1).astype(np.int64)
    # Root number W_χ = τ(χ) / (i^a · √q)
    i_pow_a = np.where(parity == 0, 1.0 + 0j, 1j)
    W_chi = tau / (i_pow_a * math.sqrt(q))

    is_prim = is_primitive_character_composite(q, char_indices, primes, group_orders)

    log_q_pi = math.log(q / math.pi)
    qp_factor = np.exp(-1j * log_q_pi * t_grid)  # (n_t,)

    z_plus_even = (0.5 + 1j * t_grid) / 2
    z_plus_odd = (1.5 + 1j * t_grid) / 2
    z_minus_even = (0.5 - 1j * t_grid) / 2
    z_minus_odd = (1.5 - 1j * t_grid) / 2
    rg_even = np.exp(loggamma(z_minus_even) - loggamma(z_plus_even))
    rg_odd = np.exp(loggamma(z_minus_odd) - loggamma(z_plus_odd))
    parity_mask = parity[:, None].astype(bool)
    rg = np.where(parity_mask, rg_odd[None, :], rg_even[None, :])  # (phi, n_t)

    refl_factor = W_chi[:, None] * qp_factor[None, :] * rg
    L = M + refl_factor * np.conj(M)

    return dict(L=L, W_chi=W_chi, parity=parity, char_indices=char_indices,
                is_primitive=is_prim, primes=primes, group_orders=group_orders,
                prim_roots=prim_roots, log_table=log_table)


# =====================================================================
# Hardy Z + zero finding (composite q)
# =====================================================================


def compute_hardy_theta_composite(q: int, parity: np.ndarray, W_chi: np.ndarray,
                                    t_grid: np.ndarray) -> np.ndarray:
    """θ_χ(t) = (t/2) log(q/π) + arg Γ((½+a+it)/2) - arg(W_χ)/2.
    Returns (phi, n_t) real."""
    log_q_pi = math.log(q / math.pi)
    arg_W = np.angle(W_chi)

    z_even = (0.5 + 1j * t_grid) / 2
    z_odd = (1.5 + 1j * t_grid) / 2
    arg_g_even = np.imag(loggamma(z_even))
    arg_g_odd = np.imag(loggamma(z_odd))

    parity_mask = parity[:, None].astype(bool)
    arg_gamma = np.where(parity_mask, arg_g_odd[None, :], arg_g_even[None, :])

    return (t_grid[None, :] / 2.0) * log_q_pi + arg_gamma - arg_W[:, None] / 2.0


def zero_brackets_vectorised(Z: np.ndarray, t_grid: np.ndarray,
                              t_min: float = 2.0) -> Dict[int, List[float]]:
    """Same as slot-3: vectorised sign-change bracketing + linear-interp
    refinement. Z shape (phi, n_t)."""
    n_chars, n_t = Z.shape
    sgn = np.sign(Z)
    cross_mask = (sgn[:, :-1] * sgn[:, 1:]) < 0

    valid_t = t_grid[:-1] >= t_min
    cross_mask &= valid_t[None, :]

    dt = t_grid[1:] - t_grid[:-1]
    dZ = Z[:, 1:] - Z[:, :-1]
    safe_dZ = np.where(np.abs(dZ) > 1e-300, dZ, 1.0)
    t_zero = t_grid[None, :-1] - Z[:, :-1] * dt[None, :] / safe_dZ

    j_idx, i_idx = np.where(cross_mask)
    refined = t_zero[j_idx, i_idx]

    db: Dict[int, List[float]] = {j: [] for j in range(n_chars)}
    if len(j_idx) == 0:
        return db
    sort_order = np.argsort(j_idx, kind="stable")
    j_sorted = j_idx[sort_order]
    refined_sorted = refined[sort_order]
    boundaries = np.searchsorted(j_sorted, np.arange(1, n_chars))
    splits = np.split(refined_sorted, boundaries)
    for j, arr in enumerate(splits):
        if arr.size > 0:
            db[j] = sorted(arr.tolist())
    return db


def afe_truncation_oversized(q: int, t_max: float, oversize: float = 12.0) -> int:
    return max(int(math.ceil(oversize * math.sqrt(q * t_max / (2.0 * math.pi)))), 30)


def find_zeros_all_chars_composite(
    q: int, t_max: float, n_t: int,
    n_zeros_target: int | None = None,
    t_min: float = 2.0,
    method: str = "fft",
    oversize: float = 12.0,
) -> Tuple[Dict[int, List[float]], Dict[str, float]]:
    """Find L-zeros for all characters χ mod q (composite squarefree odd q)
    in (t_min, t_max] via Hardy Z + sign-change bracketing on the multi-axis
    FFT primitive."""
    t0_total = time.perf_counter()

    if t_min < 1.0:
        t_min = 1.0
    t_grid = np.linspace(t_min, t_max, n_t)

    N = afe_truncation_oversized(q, t_max, oversize=oversize)

    t_setup_0 = time.perf_counter()
    primes, group_orders, prim_roots = crt_decomp(q)
    log_table = build_crt_log_table(q, primes, prim_roots)
    t_setup = time.perf_counter() - t_setup_0

    t_afe_0 = time.perf_counter()
    afe_data = compute_full_l_via_afe_composite(q, t_grid, N, method=method)
    t_afe = time.perf_counter() - t_afe_0

    L = afe_data["L"]
    W_chi = afe_data["W_chi"]
    parity = afe_data["parity"]
    is_prim = afe_data["is_primitive"]

    t_theta_0 = time.perf_counter()
    theta = compute_hardy_theta_composite(q, parity, W_chi, t_grid)
    Z = np.real(np.exp(1j * theta) * L)
    t_theta = time.perf_counter() - t_theta_0

    t_bracket_0 = time.perf_counter()
    zeros_db_all = zero_brackets_vectorised(Z, t_grid, t_min=t_min)

    # Restrict primitive characters; for non-primitive, fall back to using
    # the principal character's zeta zeros (only meaningful for full
    # principal character at q; non-primitive non-principal would need
    # the inducing primitive character's zeros, which slot 4 doesn't
    # implement — they're flagged as a falsifier for slot 5).
    zeros_db: Dict[int, List[float]] = {}
    phi_total = L.shape[0]
    for j in range(phi_total):
        if is_prim[j]:
            zeros_db[j] = zeros_db_all[j][:n_zeros_target] if n_zeros_target else zeros_db_all[j]
        elif tuple(afe_data["char_indices"][j]) == tuple([0] * len(group_orders)):
            # Principal character
            K_principal = max(n_zeros_target, 100) if n_zeros_target else 1000
            zeta_zeros = get_chi0_zeros(K_principal)
            zeta_zeros_in_range = [g for g in zeta_zeros if t_min <= g <= t_max]
            zeros_db[j] = zeta_zeros_in_range[:n_zeros_target] if n_zeros_target else zeta_zeros_in_range
        else:
            # Non-principal non-primitive: leave empty (slot 4 limitation)
            zeros_db[j] = []
    t_bracket = time.perf_counter() - t_bracket_0

    t_total = time.perf_counter() - t0_total

    timing = {
        "t_total_s": t_total,
        "t_setup_s": t_setup,
        "t_afe_s": t_afe,
        "t_theta_s": t_theta,
        "t_bracket_s": t_bracket,
        "n_t": n_t,
        "N_afe": N,
        "phi": phi_total,
        "n_primitive": int(is_prim.sum()),
    }
    return zeros_db, timing


# =====================================================================
# Sanity check vs mpmath at small composite q
# =====================================================================


def mpmath_first_zeros_composite(q: int, char_index_tuple: Tuple[int, ...],
                                   t_max: float = 30.0,
                                   n_zeros: int = 5,
                                   N_factor: int = 80) -> List[float]:
    """Use mpmath to compute first L-zeros for χ mod q via Hardy Z
    sign-change bracketing on a fine grid with a much larger AFE
    truncation than slot 4. Slow but accurate; only for small q sanity.
    """
    import mpmath as mp
    mp.mp.dps = 30
    primes, group_orders, prim_roots = crt_decomp(q)
    log_table = build_crt_log_table(q, primes, prim_roots)

    omega = [mp.exp(2j * mp.pi / order) for order in group_orders]

    def chi_val(n: int) -> complex:
        if math.gcd(n, q) != 1:
            return mp.mpc(0)
        tup = log_table[n % q]
        val = mp.mpc(1)
        for i, k in enumerate(tup):
            j = char_index_tuple[i]
            val *= omega[i] ** ((j * k) % group_orders[i])
        return val

    def L_func(s):
        N = int(N_factor * (1 + abs(mp.im(s))))
        total = mp.mpc(0)
        for n in range(1, N + 1):
            cv = chi_val(n)
            if cv != 0:
                total += cv / mp.mpc(n) ** s
        return total

    # Compute χ(-1) directly to get parity
    parity = 0 if chi_val(q - 1).real > 0 else 1

    # Compute root number τ(χ) / (i^a √q)
    tau = mp.mpc(0)
    for r in range(1, q):
        cv = chi_val(r)
        if cv != 0:
            tau += cv * mp.exp(2j * mp.pi * r / q)
    i_pow_a = mp.mpc(1) if parity == 0 else mp.mpc(0, 1)
    W_chi = tau / (i_pow_a * mp.sqrt(q))

    log_q_pi = mp.log(mp.mpf(q) / mp.pi)

    def hardy_z(t):
        s = mp.mpc(mp.mpf("0.5"), t)
        L_val = L_func(s)
        # θ(t) = (t/2) log(q/π) + arg Γ((1/2+a+it)/2) - arg(W_χ)/2
        z_arg = (mp.mpf("0.5") + parity + 1j * t) / 2
        theta = (mp.mpf(t) / 2) * log_q_pi + mp.arg(mp.gamma(z_arg)) - mp.arg(W_chi) / 2
        return float((mp.exp(1j * theta) * L_val).real)

    # Fine grid for sign-change bracketing
    n_grid = 300
    dt = (t_max - 1.0) / (n_grid - 1)
    t_grid = [1.0 + i * dt for i in range(n_grid)]
    Z_grid = [hardy_z(t) for t in t_grid]

    zeros: List[float] = []
    for i in range(n_grid - 1):
        if Z_grid[i] * Z_grid[i + 1] < 0 and t_grid[i] >= 2.0:
            # Linear interp refine
            t_zero = t_grid[i] - Z_grid[i] * (t_grid[i + 1] - t_grid[i]) / (Z_grid[i + 1] - Z_grid[i])
            zeros.append(t_zero)
    return sorted(zeros)[:n_zeros]


def sanity_check_composite_q(q: int, char_index_tuple: Tuple[int, ...],
                                t_max: float = 30.0, n_t: int = 4000) -> None:
    """Cross-check slot-4 zero finder vs mpmath ground truth at composite q."""
    print(f"\n=== Sanity check: q={q} (factors {factor_squarefree(q)}), "
          f"χ_index={char_index_tuple}, t_max={t_max} ===")
    zeros_db, timing = find_zeros_all_chars_composite(q, t_max=t_max, n_t=n_t,
                                                        n_zeros_target=10, t_min=2.0)
    # Find flat index for the requested char
    primes, group_orders, _ = crt_decomp(q)
    flat = 0
    for i, j_i in enumerate(char_index_tuple):
        flat = flat * group_orders[i] + j_i
    fft_zeros = zeros_db.get(flat, [])[:5]
    print(f"  slot-4 zero finder (first 5): {[f'{g:.3f}' for g in fft_zeros]}")
    print(f"  timing: {timing}")
    try:
        mpm_zeros = mpmath_first_zeros_composite(q, char_index_tuple, t_max=t_max, n_zeros=5)
        print(f"  mpmath ground truth (first 5): {[f'{g:.3f}' for g in mpm_zeros]}")
        if fft_zeros and mpm_zeros:
            n_compare = min(len(fft_zeros), len(mpm_zeros))
            errs = [abs(fft_zeros[i] - mpm_zeros[i]) for i in range(n_compare)]
            print(f"  abs differences: {[f'{e:.4f}' for e in errs]}")
            print(f"  max diff: {max(errs):.4f}")
    except Exception as e:
        print(f"  mpmath check failed: {e}")


# =====================================================================
# End-to-end π(x; q, a) timing
# =====================================================================


def psi_chi_partial_sum(zeros: List[float], x: float) -> Tuple[float, float]:
    """psi(x, χ) ≈ -Σ_ρ x^ρ/ρ for ρ = ½ + iγ.

    Returns (real, imag) parts.
    """
    g = np.asarray(zeros, dtype=np.float64)
    if g.size == 0:
        return 0.0, 0.0
    log_x = math.log(x)
    sqrt_x = math.sqrt(x)
    cos_t = np.cos(g * log_x)
    sin_t = np.sin(g * log_x)
    denom = 0.25 + g * g
    inv_re = 0.5 / denom
    inv_im = -g / denom
    term_re = (cos_t * inv_re - sin_t * inv_im) * sqrt_x
    term_im = (cos_t * inv_im + sin_t * inv_re) * sqrt_x
    return float(term_re.sum()), float(term_im.sum())


def explicit_formula_pi_q_a_composite(q: int, a: int, x: float,
                                         zeros_db: Dict[int, List[float]],
                                         char_indices: np.ndarray,
                                         group_orders: List[int],
                                         primes: List[int],
                                         prim_roots: List[int],
                                         log_table: Dict[int, Tuple[int, ...]]) -> float:
    """Approximate π(x; q, a) via explicit formula.

    π(x; q, a) ≈ (1/φ(q)) Σ_χ χ̄(a) π(x, χ).
    """
    if a not in log_table:
        raise ValueError(f"a={a} not coprime to q={q}")
    a_tup = log_table[a]
    log_x = math.log(x)
    phi_total = char_indices.shape[0]

    # Compute χ̄(a) for all χ
    K = len(group_orders)
    chi_bar_a = np.ones(phi_total, dtype=complex)
    for i in range(K):
        omega_i = np.exp(-2j * np.pi / group_orders[i])  # ω^{-1} = conj
        j_i_per_char = char_indices[:, i]
        chi_bar_a *= np.exp(-2j * np.pi * (j_i_per_char * a_tup[i]) / group_orders[i])

    pi_q_a = 0.0
    for j in range(phi_total):
        zeros = zeros_db.get(j, [])
        if not zeros:
            continue
        psi_re, psi_im = psi_chi_partial_sum(zeros, x)
        psi_complex = complex(psi_re, psi_im)
        pi_chi_term = -psi_complex / log_x  # leading approximation
        contribution = (chi_bar_a[j] * pi_chi_term).real
        pi_q_a += contribution / phi_total

    # Add principal-character main term π(x)/φ(q) (placeholder x/log x)
    pi_q_a += (x / log_x) / phi_total
    return pi_q_a


def time_end_to_end_pi_q_a_composite(q: int, a: int, x: float,
                                        K_target: int = 200,
                                        t_max: float | None = None,
                                        n_t: int | None = None,
                                        method: str = "fft",
                                        n_repeats: int = 3) -> Dict[str, float]:
    """End-to-end timing for composite q."""
    primes, group_orders, prim_roots = crt_decomp(q)
    log_table = build_crt_log_table(q, primes, prim_roots)
    phi_total = euler_phi(q)

    if t_max is None:
        T = 50.0
        for _ in range(20):
            T = 2 * math.pi * K_target / max(math.log(q * T / (2 * math.pi)) - 1.0, 1.0)
        t_max = T * 1.1

    if n_t is None:
        mean_spacing = 2.0 * math.pi / max(math.log(q * t_max / (2 * math.pi)), 1.0)
        n_t = int(math.ceil(t_max / (0.3 * mean_spacing)))

    times_zerofind = []
    times_full = []
    pi_qa = 0.0
    char_indices_cached = None
    for _ in range(n_repeats):
        t_zf_0 = time.perf_counter()
        zeros_db, _ = find_zeros_all_chars_composite(q, t_max=t_max, n_t=n_t,
                                                        n_zeros_target=K_target,
                                                        method=method)
        t_zf = time.perf_counter() - t_zf_0
        times_zerofind.append(t_zf)

        if char_indices_cached is None:
            _, _, char_indices_cached = l_eval_direct_composite(q, np.array([2.0]), 1)

        t_full_0 = time.perf_counter()
        pi_qa = explicit_formula_pi_q_a_composite(
            q, a, x, zeros_db, char_indices_cached, group_orders,
            primes, prim_roots, log_table)
        t_full = time.perf_counter() - t_full_0
        times_full.append(t_full)

    n_zeros_per_char = [len(v) for v in zeros_db.values()]
    return {
        "q": q,
        "phi": phi_total,
        "x": x,
        "K_target": K_target,
        "t_max": t_max,
        "n_t": n_t,
        "method": method,
        "t_zerofind_s": float(np.median(times_zerofind)),
        "t_partialsum_s": float(np.median(times_full)),
        "t_total_single_query_s": float(np.median(times_zerofind)) + float(np.median(times_full)),
        "min_zeros_per_char": min(n_zeros_per_char) if n_zeros_per_char else 0,
        "median_zeros_per_char": int(np.median(n_zeros_per_char)) if n_zeros_per_char else 0,
        "max_zeros_per_char": max(n_zeros_per_char) if n_zeros_per_char else 0,
        "pi_q_a_estimate": pi_qa,
    }


# =====================================================================
# Main
# =====================================================================


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true",
                        help="small grid for smoke test")
    parser.add_argument("--no-mpmath", action="store_true",
                        help="skip mpmath ground-truth comparison")
    parser.add_argument("--skip-end-to-end", action="store_true",
                        help="skip end-to-end timing (just verify primitive)")
    args = parser.parse_args()

    print("=" * 80)
    print("Slot 4 — composite-q multi-axis FFT primitive correctness")
    print("=" * 80)

    correctness_rows = []
    for q in [15, 35, 105, 1001]:
        try:
            err, rel = check_multi_axis_correctness(q, t_grid_size=5, t_max=30.0, N=80)
            correctness_rows.append(dict(q=q, err=err, rel_err=rel,
                                          phi=euler_phi(q),
                                          factors=str(factor_squarefree(q))))
        except Exception as e:
            print(f"  q={q}: FAILED with {e}")
            correctness_rows.append(dict(q=q, err=float("nan"), rel_err=float("nan"),
                                          phi=euler_phi(q),
                                          factors=str(factor_squarefree(q))))

    print()
    print("=" * 80)
    print("Slot 4 — sanity check vs mpmath at small composite q")
    print("=" * 80)
    if not args.no_mpmath:
        sanity_check_composite_q(q=15, char_index_tuple=(1, 1), t_max=40.0, n_t=4000)
        sanity_check_composite_q(q=35, char_index_tuple=(1, 1), t_max=30.0, n_t=4000)
    else:
        print("  (skipped per --no-mpmath)")

    if args.skip_end_to_end:
        print("(skipped end-to-end per --skip-end-to-end)")
        return

    print()
    print("=" * 80)
    print("Slot 4 — end-to-end π(x; q, a) timing at composite q")
    print("=" * 80)

    if args.quick:
        configs = [(15, 2, 1e6, 50), (35, 3, 1e6, 50)]
    else:
        configs = [
            (15,    2, 1e6,  50),
            (35,    3, 1e6,  50),
            (105,   2, 1e6,  50),
            (105,   2, 1e6, 200),
            (1001,  2, 1e6, 200),
        ]

    e2e_rows = []
    for q, a, x, K in configs:
        t0 = time.perf_counter()
        try:
            result = time_end_to_end_pi_q_a_composite(q, a, x, K_target=K,
                                                        method="fft", n_repeats=3)
            elapsed = time.perf_counter() - t0
            print(f"  q={result['q']:>5} a={a} x={result['x']:.0e} K={result['K_target']:>4}: "
                  f"t_zf={result['t_zerofind_s']*1e3:7.2f}ms  "
                  f"t_psum={result['t_partialsum_s']*1e3:7.2f}ms  "
                  f"t_single={result['t_total_single_query_s']*1e3:7.2f}ms  "
                  f"phi={result['phi']:>4}  "
                  f"zeros[med]={result['median_zeros_per_char']:>4}  "
                  f"({elapsed:.1f}s wall)")
            e2e_rows.append(result)
        except Exception as e:
            print(f"  q={q} FAILED: {e}")

    out_csv = os.path.join(HERE, "slot4_composite_q_end_to_end.csv")
    if e2e_rows:
        fieldnames = list(e2e_rows[0].keys())
        with open(out_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            for r in e2e_rows:
                w.writerow(r)
        print(f"\nwrote {out_csv}")

    # FFT vs DIRECT comparison at composite q
    if not args.quick:
        print()
        print("=" * 80)
        print("Slot 4 — FFT vs DIRECT (no FFT sharing) at composite q")
        print("=" * 80)
        fvd_rows = []
        for q, K in [(105, 50), (1001, 200)]:
            T = 50.0
            for _ in range(20):
                T = 2 * math.pi * K / max(math.log(q * T / (2 * math.pi)) - 1.0, 1.0)
            t_max = T * 1.1
            mean_spacing = 2.0 * math.pi / max(math.log(q * t_max / (2 * math.pi)), 1.0)
            n_t = int(math.ceil(t_max / (0.3 * mean_spacing)))
            print(f"\n--- q={q} K={K} t_max={t_max:.2f} n_t={n_t} ---")
            t_fft_runs, t_dir_runs = [], []
            for _ in range(3):
                t0 = time.perf_counter()
                _, _ = find_zeros_all_chars_composite(q, t_max=t_max, n_t=n_t,
                                                          n_zeros_target=K, method="fft")
                t_fft_runs.append(time.perf_counter() - t0)
                t0 = time.perf_counter()
                _, _ = find_zeros_all_chars_composite(q, t_max=t_max, n_t=n_t,
                                                          n_zeros_target=K, method="direct")
                t_dir_runs.append(time.perf_counter() - t0)
            t_fft = float(np.median(t_fft_runs))
            t_dir = float(np.median(t_dir_runs))
            speedup = t_dir / t_fft if t_fft > 0 else float("nan")
            fvd_rows.append(dict(q=q, K=K, t_max=t_max, n_t=n_t,
                                  t_fft_s=t_fft, t_direct_s=t_dir,
                                  fft_speedup=speedup,
                                  phi=euler_phi(q),
                                  factors=str(factor_squarefree(q))))
            print(f"  FFT method:    {t_fft*1e3:7.2f} ms")
            print(f"  DIRECT method: {t_dir*1e3:7.2f} ms")
            print(f"  FFT speedup:   {speedup:5.2f}x")

        if fvd_rows:
            fvd_csv = os.path.join(HERE, "slot4_composite_fft_vs_direct.csv")
            with open(fvd_csv, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(fvd_rows[0].keys()))
                w.writeheader()
                for r in fvd_rows:
                    w.writerow(r)
            print(f"\nwrote {fvd_csv}")


if __name__ == "__main__":
    main()
