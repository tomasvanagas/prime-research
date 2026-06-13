"""
Thread 6, slot 3: end-to-end Dirichlet-L zero finder using FFT-shared
AFE + Hardy Z + sign-change bracketing.

Mission (from .commit_state Thread 6 slot 3 PRIORITY (a) +
session228_commit_p1_afe_shared_l_eval_results.md "What would qualify
as slot-3 progress"):
  Build a real zero-finder using the FFT-shared L-eval primitive,
  measure end-to-end per-(q, a) pi computation cost vs slot-1 baseline
  (q=1009, K=200: 188ms single-query, 23ms amortised at M=256).

  If the cost drops by >=3x single-query (consistent with slot-2's
  L-eval speedup propagating end-to-end), that confirms slot-2's
  primitive lifts to a real partial-positive on AP-pi.

The pipeline:

  Step 1 (FFT). Compute the AFE main term M(t, chi) for ALL chi mod q
                at all t in t_grid via slot-2's `phi * ifft(W, axis=0)`
                identity. Single length-phi FFT per t.

                M(t, chi) = sum_{n <= N, gcd(n,q)=1} chi(n) / n^{1/2+it}

  Step 2 (REFL). For primitive chi, the truncated AFE includes a
                 reflected piece. Beautifully, for prime q where chi
                 is primitive non-principal:

                   sum_{n <= N} chibar(n) / n^{1/2-it} = conj(M(t, chi))

                 because chibar(n) = conj(chi(n)) and n^{1/2-it} =
                 conj(n^{1/2+it}). So the reflected piece is

                   W_chi * (q/pi)^{-it} * Gamma_ratio(t, a) * conj(M)

                 where W_chi = tau(chi) / (i^a sqrt(q)) is the root
                 number with |W_chi| = 1, a = 0 if chi even, 1 if odd.
                 NO ADDITIONAL FFT. Pure pointwise multiply.

  Step 3 (Z). Hardy Z function:

                Z_chi(t) = exp(i * theta_chi(t)) * L(1/2 + it, chi)

              with

                theta_chi(t) = (t/2) * log(q/pi)
                             + arg Gamma((1/2 + a + it)/2)
                             - arg(W_chi) / 2

              Z_chi(t) is real on the critical line. Sign changes <=>
              zeros.

  Step 4 (BRACKET). For each non-principal chi_j, scan Z_chi_j over
                    t_grid for sign changes. Linear-interpolate
                    refinement: t_zero = t_i - Z_i (t_{i+1}-t_i)/(Z_{i+1}-Z_i).

  Step 5 (PI). Plug refined gammas into slot 1's partial-sum evaluator
               for psi(x, chi). Combine via orthogonality. Return
               pi(x; q, a).

What would falsify the slot-3 prediction (single-query speedup ~3x or
better):

  (i)  Real zero-finding is dominated by post-FFT scanning / refinement
       per zero (need ~K * phi total brackets), and these costs scale
       linearly with K * phi independent of FFT primitive. Then
       speedup is bounded.

  (ii) AFE truncation error is too large at the t-heights needed
       (T ~ 130 for K=200, q=1009 to capture all the way up the
       critical line). Then accuracy fails and we cannot match
       slot-1's K=200 zero count.

  (iii) Reflected-term phase computation (loggamma) is itself slow
        (scipy.special.loggamma is sequential per-element); mitigation:
        vectorise.

What this slot DOES NOT do (yet):

  - Does not produce final pi(x; q, a) = exact integer (only the
    explicit-formula approximation).
  - Does not test composite-q (PRIORITY b - slot 4) or cross-conductor
    batches (PRIORITY c - slot 5).
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

from per_q_a_explicit_formula_profile import (  # noqa: E402
    euler_phi,
    primitive_root,
    get_chi0_zeros,
    time_per_chi_eval,
    time_orthogonality,
)
from slot2_afe_shared_l_eval import (  # noqa: E402
    coprimes_in_log_g_order,
    compute_complex_pow_all_n,
)


# AFE truncation: balanced N = sqrt(q*t/(2pi)) gives error ~1 (too large).
# Oversize N by a factor; remainder ~ N_balanced/N for smooth-windowed
# AFEs, ~1/sqrt(N) for hard-truncated. Use OVERSIZE = 12 for accurate
# zero detection at small t; slot 3 prefers correctness over speed in
# the L-eval primitive (slot 2 confirmed FFT primitive itself works).
def afe_truncation_balanced(q: int, t_max: float) -> int:
    return max(int(math.ceil(math.sqrt(q * t_max / (2.0 * math.pi)))), 10)


def afe_truncation_oversized(q: int, t_max: float, oversize: float = 12.0) -> int:
    return max(int(math.ceil(oversize * math.sqrt(q * t_max / (2.0 * math.pi)))), 30)


# =====================================================================
# Full AFE: main term + reflected term, all chi mod q at once via FFT
# =====================================================================


def compute_main_term_fft(q: int, t_grid: np.ndarray, N: int) -> np.ndarray:
    """Slot-2's L_main = phi * ifft(W, axis=0) where W is in log-g order.
    Returns shape (phi, n_t)."""
    coprimes_logg, log_g = coprimes_in_log_g_order(q)
    phi = len(coprimes_logg)
    coprime_set = set(coprimes_logg)
    n_t = len(t_grid)

    cp_all = compute_complex_pow_all_n(t_grid, N)  # (N, n_t)
    W = np.zeros((phi, n_t), dtype=complex)
    for n_idx in range(N):
        r = (n_idx + 1) % q
        if r in coprime_set:
            W[log_g[r]] += cp_all[n_idx]
    return phi * np.fft.ifft(W, axis=0)


def compute_main_term_direct(q: int, t_grid: np.ndarray, N: int) -> np.ndarray:
    """Per-character DIRECT matmul reference: chi_at_n @ complex_pow.
    Used as baseline to compare against the FFT primitive contribution.
    Returns shape (phi, n_t)."""
    from slot2_afe_shared_l_eval import char_table_at_residues
    n_arr = np.arange(1, N + 1)
    n_residues = (n_arr % q).astype(int)
    coprime_mask = np.array([math.gcd(int(r), q) == 1 and int(r) != 0 for r in n_residues])
    n_co = n_arr[coprime_mask]
    n_co_residues = (n_co % q).astype(int).tolist()

    log_n = np.log(n_co.astype(np.float64))
    inv_sqrt_n = 1.0 / np.sqrt(n_co.astype(np.float64))
    phase = -np.outer(log_n, t_grid)
    cp = inv_sqrt_n[:, None] * np.exp(1j * phase)

    chi_at_n = char_table_at_residues(q, n_co_residues)
    return chi_at_n @ cp


def compute_gauss_sums_fft(q: int) -> np.ndarray:
    """tau(chi_j) = sum_a chi_j(a) e^{2 pi i a / q}.
    For prime q with primitive root g and j-th character chi_j(g^k) = omega^{jk}:
        tau(chi_j) = sum_k omega^{jk} e^{2 pi i g^k / q} = phi * ifft(e_a)[j]
    where e_a[k] = e^{2 pi i g^k / q}.
    Returns shape (phi,) complex."""
    coprimes_logg, _ = coprimes_in_log_g_order(q)
    phi = len(coprimes_logg)
    e_a = np.exp(2j * np.pi * np.array(coprimes_logg, dtype=np.float64) / q)
    return phi * np.fft.ifft(e_a)


def compute_full_l_via_afe(
    q: int, t_grid: np.ndarray, N: int,
    parities: np.ndarray, W_chi: np.ndarray,
    method: str = "fft",
) -> np.ndarray:
    """Full AFE: L(1/2+it, chi) = M + W_chi * (q/pi)^{-it} * Gamma_ratio(t, a) * conj(M).
    Returns shape (phi, n_t) complex.

    Special note: this formula is for PRIMITIVE chi only. For prime q,
    every non-principal character is primitive. The principal character
    L(s, chi_0) = zeta(s) * (1 - q^{-s}) is handled separately.

    method ∈ {"fft", "direct"} selects the main-term primitive.
    """
    if method == "fft":
        M = compute_main_term_fft(q, t_grid, N)  # (phi, n_t)
    elif method == "direct":
        M = compute_main_term_direct(q, t_grid, N)
    else:
        raise ValueError(f"unknown method {method!r}")

    log_q_pi = math.log(q / math.pi)
    qp_factor = np.exp(-1j * log_q_pi * t_grid)  # (n_t,)

    # Gamma((1/2 - it + a)/2) / Gamma((1/2 + it + a)/2)
    # using exp(loggamma(z_-) - loggamma(z_+))
    z_plus_even = (0.5 + 1j * t_grid) / 2
    z_plus_odd = (1.5 + 1j * t_grid) / 2
    z_minus_even = (0.5 - 1j * t_grid) / 2
    z_minus_odd = (1.5 - 1j * t_grid) / 2
    rg_even = np.exp(loggamma(z_minus_even) - loggamma(z_plus_even))
    rg_odd = np.exp(loggamma(z_minus_odd) - loggamma(z_plus_odd))

    # rg has |rg| = 1 since loggamma(z*) = conj(loggamma(z)) for z not on negative real axis;
    # the imag part of loggamma(z_-) - loggamma(z_+) is twice the negative imag of loggamma(z_+),
    # and the real parts cancel.

    # Per-character ratio_gamma based on parity
    parity_mask = parities[:, None].astype(bool)  # (phi, 1) where True iff odd
    rg = np.where(parity_mask, rg_odd[None, :], rg_even[None, :])  # (phi, n_t)

    refl_factor = W_chi[:, None] * qp_factor[None, :] * rg  # (phi, n_t)
    return M + refl_factor * np.conj(M)


# =====================================================================
# Hardy Z + sign-change zero finder (per character)
# =====================================================================


def compute_hardy_theta(
    q: int, parities: np.ndarray, W_chi: np.ndarray, t_grid: np.ndarray
) -> np.ndarray:
    """theta_chi_j(t) = (t/2)*log(q/pi) + arg Gamma((1/2 + a_j + it)/2) - arg(W_chi_j)/2.
    Returns shape (phi, n_t) real."""
    log_q_pi = math.log(q / math.pi)
    arg_W = np.angle(W_chi)  # (phi,)

    # arg Gamma((1/2 + a + it)/2): vectorise per parity
    z_even = (0.5 + 1j * t_grid) / 2
    z_odd = (1.5 + 1j * t_grid) / 2
    arg_g_even = np.imag(loggamma(z_even))  # (n_t,)
    arg_g_odd = np.imag(loggamma(z_odd))  # (n_t,)

    parity_mask = parities[:, None].astype(bool)  # (phi, 1)
    arg_gamma = np.where(parity_mask, arg_g_odd[None, :], arg_g_even[None, :])

    return (t_grid[None, :] / 2.0) * log_q_pi + arg_gamma - arg_W[:, None] / 2.0


def zero_brackets_linear_refine(
    Z_row: np.ndarray, t_grid: np.ndarray, t_min: float = 0.5
) -> List[float]:
    """For one character's Hardy Z values on a regular t_grid, return the
    list of refined zeros via linear interpolation. Skips the small-t
    region (below t_min) where Hardy Z prefactor diverges and AFE is
    inaccurate."""
    sgn = np.sign(Z_row)
    sign_changes = np.where(sgn[:-1] * sgn[1:] < 0)[0]
    gammas = []
    for i in sign_changes:
        if t_grid[i] < t_min:
            continue
        denom = Z_row[i + 1] - Z_row[i]
        if abs(denom) < 1e-300:
            continue
        t_zero = t_grid[i] - Z_row[i] * (t_grid[i + 1] - t_grid[i]) / denom
        gammas.append(float(t_zero))
    return gammas


def zero_brackets_vectorised(
    Z: np.ndarray, t_grid: np.ndarray, t_min: float = 2.0
) -> Dict[int, List[float]]:
    """Vectorised across all characters. Z shape (n_chars, n_t).
    Returns dict[j -> sorted list of refined gammas].

    Uses numpy bulk operations for sign-change detection + linear
    interpolation; per-character split is the only Python work.
    """
    n_chars, n_t = Z.shape
    sgn = np.sign(Z)
    cross_mask = (sgn[:, :-1] * sgn[:, 1:]) < 0  # (n_chars, n_t-1) bool

    valid_t = t_grid[:-1] >= t_min  # (n_t-1,)
    cross_mask &= valid_t[None, :]

    dt = t_grid[1:] - t_grid[:-1]  # (n_t-1,)
    dZ = Z[:, 1:] - Z[:, :-1]
    safe_dZ = np.where(np.abs(dZ) > 1e-300, dZ, 1.0)
    t_zero = t_grid[None, :-1] - Z[:, :-1] * dt[None, :] / safe_dZ  # (n_chars, n_t-1)

    j_idx, i_idx = np.where(cross_mask)
    refined = t_zero[j_idx, i_idx]  # (n_total_zeros,)

    db: Dict[int, List[float]] = {j: [] for j in range(n_chars)}
    if len(j_idx) == 0:
        return db
    # Group by j -- np.split on sorted indices.
    sort_order = np.argsort(j_idx, kind="stable")
    j_sorted = j_idx[sort_order]
    refined_sorted = refined[sort_order]
    boundaries = np.searchsorted(j_sorted, np.arange(1, n_chars))
    splits = np.split(refined_sorted, boundaries)
    for j, arr in enumerate(splits):
        if arr.size > 0:
            db[j] = sorted(arr.tolist())
    return db


def find_zeros_all_chars(
    q: int, t_max: float, n_t: int,
    n_zeros_target: int | None = None,
    t_min: float = 2.0,
    method: str = "fft",
) -> Tuple[Dict[int, List[float]], Dict[str, float]]:
    """Find zeros of L(1/2 + it, chi_j) for j = 0..phi-1 in (t_min, t_max].

    j = 0 (principal) is handled by reading cached zeta-zeros and
    multiplying by the (1 - q^{-s})-correction-factor (which has its own
    zeros on Re(s) = 0 only -- no impact on critical-line zeros).

    Returns (zeros_db, timing_dict). zeros_db[j] is the sorted list of
    refined gammas for chi_j.
    """
    t0_total = time.perf_counter()

    # Setup
    coprimes_logg, log_g = coprimes_in_log_g_order(q)
    phi = len(coprimes_logg)

    # Parity: chi_j(-1) = (-1)^j for prime q (since log_g(q-1) = phi/2)
    parities = (np.arange(phi) % 2).astype(np.int64)

    t_setup_0 = time.perf_counter()

    # Gauss sums via FFT
    gauss_sums = compute_gauss_sums_fft(q)  # (phi,)
    # Root numbers W_chi = tau(chi) / (i^a * sqrt(q))
    i_pow_a = np.where(parities == 0, 1.0 + 0j, 1j)
    W_chi = gauss_sums / (i_pow_a * math.sqrt(q))  # (phi,) with |W_chi| = 1 for primitive

    t_setup = time.perf_counter() - t_setup_0

    # Build t-grid
    if t_min < 1.0:
        t_min = 1.0
    t_grid = np.linspace(t_min, t_max, n_t)

    # AFE truncation: oversize for accurate Hardy Z values without
    # Riemann-Siegel correction terms. Balanced N would give large error.
    N = afe_truncation_oversized(q, t_max, oversize=12.0)

    t_afe_0 = time.perf_counter()
    L = compute_full_l_via_afe(q, t_grid, N, parities, W_chi, method=method)  # (phi, n_t)
    t_afe = time.perf_counter() - t_afe_0

    t_theta_0 = time.perf_counter()
    theta = compute_hardy_theta(q, parities, W_chi, t_grid)  # (phi, n_t)
    Z = np.real(np.exp(1j * theta) * L)  # (phi, n_t) real Hardy Z
    t_theta = time.perf_counter() - t_theta_0

    # Sign-change bracketing per character (vectorised)
    t_bracket_0 = time.perf_counter()
    zeros_db = zero_brackets_vectorised(Z, t_grid, t_min=t_min)
    # Principal character: load cached zeta zeros (override the AFE-based
    # finder for j=0 since L(s, chi_0) = zeta(s) * (1 - q^{-s})).
    if n_zeros_target is None:
        K_principal = 1000
    else:
        K_principal = max(n_zeros_target, 100)
    zeta_zeros = get_chi0_zeros(K_principal)
    zeta_zeros_in_range = [g for g in zeta_zeros if t_min <= g <= t_max]
    zeros_db[0] = zeta_zeros_in_range[:n_zeros_target] if n_zeros_target else zeta_zeros_in_range
    if n_zeros_target is not None:
        for j in range(1, phi):
            zeros_db[j] = zeros_db[j][:n_zeros_target]
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
        "phi": phi,
    }
    return zeros_db, timing


# =====================================================================
# Sanity: verify zero finder against mpmath ground truth (small q only)
# =====================================================================


def _build_chi_at_n_table(q: int) -> np.ndarray:
    """Returns table[j, n] = chi_j(n) for n in [1, q-1] (with n coprime: phi values),
    indexed in residue order n=1..q-1 (zeros at non-coprime residues)."""
    g = primitive_root(q)
    phi = euler_phi(q)
    log_g_dict: Dict[int, int] = {}
    cur = 1
    for k in range(phi):
        log_g_dict[cur] = k
        cur = (cur * g) % q
    table = np.zeros((phi, q), dtype=complex)
    omega = np.exp(2j * np.pi / phi)
    for n in range(1, q):
        if math.gcd(n, q) == 1:
            k = log_g_dict[n]
            for j in range(phi):
                table[j, n] = omega ** ((j * k) % phi)
    return table


def mpmath_first_zeros(q: int, j: int, n_zeros: int = 3) -> List[float]:
    """Use mpmath to compute first n_zeros L-zeros for chi_j mod q.
    Linear search over a fine t-grid + secant refinement on Hardy Z.
    Slow; only for small q sanity."""
    import mpmath as mp
    mp.mp.dps = 25
    chi_table = _build_chi_at_n_table(q)
    chi_vals = chi_table[j]  # (q,)

    def L_func(s):
        # Approximate via Dirichlet L-series (slow but accurate for small q).
        # Use the AFE-style truncated sum with reasonable N.
        N = int(50 * (1 + abs(mp.im(s))))
        total = mp.mpc(0)
        for n in range(1, N + 1):
            r = n % q
            if r == 0:
                continue
            cv = chi_vals[r]
            if abs(cv) < 1e-15:
                continue
            cv_mp = mp.mpc(cv.real, cv.imag)
            total += cv_mp / mp.mpc(n) ** s
        return total

    # mpmath has L-series functions but constructing a custom char is awkward;
    # do a manual scan + secant.
    t_lo = 1.0
    t_hi = 50.0
    n_grid = 300
    t_grid = [t_lo + (t_hi - t_lo) * i / (n_grid - 1) for i in range(n_grid)]

    # Compute |L|^2 on grid
    abs_L_sq = []
    for t in t_grid:
        s = mp.mpc(0.5, t)
        L_val = L_func(s)
        abs_L_sq.append(float(abs(L_val)) ** 2)

    # Find local minima
    minima = []
    for i in range(1, n_grid - 1):
        if abs_L_sq[i] < abs_L_sq[i - 1] and abs_L_sq[i] < abs_L_sq[i + 1]:
            if abs_L_sq[i] < 0.05 * float(np.median(abs_L_sq)):
                # Refine with secant on Re(L) AND Im(L)... we want |L| = 0
                # Use parabolic min on |L|^2 instead
                a = abs_L_sq[i - 1]
                b = abs_L_sq[i]
                c = abs_L_sq[i + 1]
                denom = a - 2 * b + c
                if abs(denom) < 1e-300:
                    minima.append(t_grid[i])
                else:
                    delta = 0.5 * (a - c) / denom
                    minima.append(t_grid[i] + delta * (t_grid[1] - t_grid[0]))
    return sorted(minima)[:n_zeros]


def sanity_check_small_q(q: int = 11, j: int = 1, t_max: float = 50.0,
                          n_t: int = 4000) -> None:
    """Compare slot-3 zero finder vs mpmath ground truth at small q."""
    print(f"\n=== Sanity check: q={q}, chi_{j}, t_max={t_max} ===")
    zeros_db, timing = find_zeros_all_chars(q, t_max=t_max, n_t=n_t,
                                             n_zeros_target=10, t_min=2.0)
    fft_zeros = zeros_db[j][:5]
    print(f"  FFT-based zero finder (chi_{j}, first 5): {[f'{g:.3f}' for g in fft_zeros]}")
    print(f"  timing: {timing}")
    try:
        mpm_zeros = mpmath_first_zeros(q, j, n_zeros=5)
        print(f"  mpmath ground truth (chi_{j}, first 5): {[f'{g:.3f}' for g in mpm_zeros]}")
        if fft_zeros and mpm_zeros:
            n_compare = min(len(fft_zeros), len(mpm_zeros))
            errs = [abs(fft_zeros[i] - mpm_zeros[i]) for i in range(n_compare)]
            print(f"  abs differences (first {n_compare}): "
                  f"{[f'{e:.4f}' for e in errs]}")
            print(f"  max diff: {max(errs):.4f}")
    except Exception as e:
        print(f"  mpmath check failed: {e}")


# =====================================================================
# End-to-end: pi(x; q, a) via explicit formula given zeros + partial sums
# =====================================================================


def psi_chi_partial_sum(zeros: List[float], x: float) -> Tuple[float, float]:
    """Compute psi(x, chi) ~ -sum_rho x^rho / rho where rho = 1/2 + i*gamma.

    Returns (real_part, imag_part). For chi non-principal, psi(x, chi)
    is generally complex; for the principal character it is real.

    x^rho / rho = x^{1/2} * (cos(gamma * log x) + i sin(...)) / (1/2 + i gamma)
                = x^{1/2} * (cos+i sin) * (1/2 - i gamma) / (1/4 + gamma^2)
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
    # NOTE: real Dirichlet-L zeros come in pairs (gamma, -gamma) for chi_0;
    # but for non-principal chi the zeros are NOT symmetric around 0
    # unless chi is real (Legendre, etc.). For complex chi, we should
    # ALSO sum over zeros at conjugate gamma for chi-bar.
    # In the slot-3 simplification we only count gammas > 0; the
    # chi/chi-bar symmetry pair is handled at the orthogonality stage.
    return float(term_re.sum()), float(term_im.sum())


def explicit_formula_pi_q_a(q: int, a: int, x: float, zeros_db: Dict[int, List[float]],
                              chi_table: np.ndarray, residues_coprime: List[int]) -> float:
    """Approximate pi(x; q, a) via the explicit formula.

    pi(x; q, a) ~ (1/phi(q)) * sum_chi chibar(a) * pi(x, chi)

    where pi(x, chi) is derived from psi(x, chi) by Mobius unraveling
    (we use the leading-order approximation pi(x, chi) ~ psi(x, chi) / log x
    + Riemann-prime-counting corrections). For slot-3 we use the
    leading approximation as in slot 1.
    """
    phi = len(residues_coprime)
    a_idx = residues_coprime.index(a)
    log_x = math.log(x)
    sqrt_x = math.sqrt(x)

    pi_q_a = 0.0
    # Leading term from chi_0: pi(x; q, a) ~ pi(x) / phi(q)
    # Approximate pi(x) ~ Li(x) ~ x/log(x); slot-1 used the explicit-formula
    # approximation. Use the same leading R(x) approximation.

    # Sum over characters of chibar(a) * (-sum_rho x^rho/rho) / log x
    for j in range(phi):
        zeros = zeros_db.get(j, [])
        if not zeros:
            continue
        psi_re, psi_im = psi_chi_partial_sum(zeros, x)
        # chibar(a) = conj(chi_table[j, a_idx])
        cba = np.conj(chi_table[j, a_idx])
        # pi(x, chi) ~ -psi(x, chi) / log x for non-principal
        psi_complex = complex(psi_re, psi_im)
        pi_chi_term = -psi_complex / log_x
        contribution = (cba * pi_chi_term).real
        pi_q_a += contribution / phi

    # Add the principal-character "main term" pi(x)/phi(q): use Riemann R(x)
    # approximation. For slot-3 cost measurement, just include a placeholder
    # using x/log(x) -- the cost is in the per-chi loop above, which is
    # what we are profiling.
    pi_q_a += (x / log_x) / phi
    return pi_q_a


def time_end_to_end_pi_q_a(q: int, a: int, x: float,
                             K_target: int = 200,
                             t_max: float | None = None,
                             n_t: int | None = None,
                             n_repeats: int = 3) -> Dict[str, float]:
    """Measure end-to-end pi(x; q, a) cost using slot-3's FFT-Hardy-Z
    zero finder followed by slot-1's partial-sum evaluator.

    Default t_max chosen to capture ~K_target zeros per character via
    Riemann-von Mangoldt density inversion.
    """
    phi = euler_phi(q)
    if t_max is None:
        # N_chi(T) ~ T/(2pi) * (log(qT/(2pi)) - 1)
        # Solve T/(2pi) * (log(qT/(2pi)) - 1) = K via simple iteration
        T = 50.0
        for _ in range(20):
            T = 2 * math.pi * K_target / max(math.log(q * T / (2 * math.pi)) - 1.0, 1.0)
        t_max = T * 1.1  # safety margin

    if n_t is None:
        # Mean spacing at height T_max ~ 2pi / log(qT_max/(2pi)); we want grid
        # spacing <= 0.3 * mean spacing to avoid missing brackets.
        mean_spacing = 2.0 * math.pi / max(math.log(q * t_max / (2 * math.pi)), 1.0)
        n_t = int(math.ceil(t_max / (0.3 * mean_spacing)))

    # Build chi table once
    chi_table = np.zeros((phi, phi), dtype=complex)
    coprimes_sorted = [n for n in range(1, q) if math.gcd(n, q) == 1]
    g = primitive_root(q)
    log_g_d: Dict[int, int] = {}
    cur = 1
    for k in range(phi):
        log_g_d[cur] = k
        cur = (cur * g) % q
    omega = np.exp(2j * np.pi / phi)
    for ci, n in enumerate(coprimes_sorted):
        k = log_g_d[n]
        for j in range(phi):
            chi_table[j, ci] = omega ** ((j * k) % phi)

    # Time zero finder
    times_zerofind = []
    times_full = []
    for _ in range(n_repeats):
        t_zf_0 = time.perf_counter()
        zeros_db, _ = find_zeros_all_chars(q, t_max=t_max, n_t=n_t,
                                             n_zeros_target=K_target)
        t_zf = time.perf_counter() - t_zf_0
        times_zerofind.append(t_zf)

        # Partial-sum eval at this x (per-char, all chars)
        t_full_0 = time.perf_counter()
        pi_qa = explicit_formula_pi_q_a(q, a, x, zeros_db, chi_table, coprimes_sorted)
        t_full = time.perf_counter() - t_full_0
        times_full.append(t_full)

    # Count actual zeros found
    n_zeros_per_char = [len(zeros_db[j]) for j in zeros_db]

    return {
        "q": q,
        "phi": phi,
        "x": x,
        "K_target": K_target,
        "t_max": t_max,
        "n_t": n_t,
        "t_zerofind_s": float(np.median(times_zerofind)),
        "t_partialsum_s": float(np.median(times_full)),
        "t_total_single_query_s": float(np.median(times_zerofind)) + float(np.median(times_full)),
        "min_zeros_per_char": min(n_zeros_per_char) if n_zeros_per_char else 0,
        "median_zeros_per_char": int(np.median(n_zeros_per_char)) if n_zeros_per_char else 0,
        "max_zeros_per_char": max(n_zeros_per_char) if n_zeros_per_char else 0,
        "pi_q_a_estimate": pi_qa,
    }


# =====================================================================
# Main driver
# =====================================================================


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true",
                        help="small grid for smoke test")
    parser.add_argument("--no-mpmath", action="store_true",
                        help="skip mpmath ground-truth comparison")
    args = parser.parse_args()

    print("=" * 80)
    print("Slot 3 - Sanity check: Hardy Z zero finder vs mpmath ground truth")
    print("=" * 80)
    if not args.no_mpmath:
        sanity_check_small_q(q=11, j=1, t_max=50.0, n_t=2000)
        sanity_check_small_q(q=11, j=5, t_max=50.0, n_t=2000)
        sanity_check_small_q(q=31, j=1, t_max=40.0, n_t=2000)
    else:
        print("  (skipped per --no-mpmath)")

    print()
    print("=" * 80)
    print("Slot 3 - Internal cost decomposition (timing dict)")
    print("=" * 80)
    for q in ([101, 1009] if not args.quick else [101]):
        for K_target in [50, 200] if not args.quick else [50]:
            # Choose t_max heuristically
            T = 50.0
            for _ in range(20):
                T = 2 * math.pi * K_target / max(math.log(q * T / (2 * math.pi)) - 1.0, 1.0)
            t_max = T * 1.1
            mean_spacing = 2.0 * math.pi / max(math.log(q * t_max / (2 * math.pi)), 1.0)
            n_t = int(math.ceil(t_max / (0.3 * mean_spacing)))
            print(f"\n--- q={q} K_target={K_target} t_max={t_max:.2f} n_t={n_t} ---")
            zeros_db, timing = find_zeros_all_chars(q, t_max=t_max, n_t=n_t,
                                                      n_zeros_target=K_target)
            n_zeros_per_char = [len(zeros_db[j]) for j in zeros_db]
            print(f"  timing: {timing}")
            print(f"  zeros: min={min(n_zeros_per_char)} med={int(np.median(n_zeros_per_char))} max={max(n_zeros_per_char)}")

    print()
    print("=" * 80)
    print("Slot 3 - End-to-end pi(x; q, a) cost: slot-3 vs slot-1 baseline")
    print("=" * 80)

    rows = []
    if args.quick:
        configs = [(101, 2, 1e6, 50)]
    else:
        configs = [
            (101, 2, 1e6, 50),
            (101, 2, 1e6, 200),
            (1009, 2, 1e6, 50),
            (1009, 2, 1e6, 200),
        ]

    for q, a, x, K in configs:
        t0_setup = time.perf_counter()
        result = time_end_to_end_pi_q_a(q, a, x, K_target=K, n_repeats=3)
        elapsed = time.perf_counter() - t0_setup
        print(f"  q={result['q']:>5} a={a} x={result['x']:.0e} K={result['K_target']:>4}: "
              f"t_zf={result['t_zerofind_s']*1e3:7.2f}ms  "
              f"t_psum={result['t_partialsum_s']*1e3:7.2f}ms  "
              f"t_single={result['t_total_single_query_s']*1e3:7.2f}ms  "
              f"zeros[med]={result['median_zeros_per_char']:>4}  "
              f"pi(x;q,a)~{result['pi_q_a_estimate']:.1f}  "
              f"({elapsed:.1f}s wall)")
        rows.append(result)

    out_csv = os.path.join(HERE, "slot3_end_to_end.csv")
    if rows:
        fieldnames = list(rows[0].keys())
        with open(out_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            for r in rows:
                w.writerow(r)
        print(f"\nwrote {out_csv}")

    # FFT vs DIRECT (no FFT sharing) at q=1009, K=200 -- the within-slot-3
    # comparison that isolates the FFT primitive's contribution to the
    # end-to-end pipeline. The "DIRECT" method uses per-character matmul
    # for L-eval (slot 2's optimised baseline) instead of cyclic FFT.
    if not args.quick:
        print()
        print("=" * 80)
        print("Slot 3 - FFT vs DIRECT (no FFT sharing) -- contribution of primitive")
        print("=" * 80)
        comparisons = []
        for q, K in [(101, 200), (1009, 200)]:
            T = 50.0
            for _ in range(20):
                T = 2 * math.pi * K / max(math.log(q * T / (2 * math.pi)) - 1.0, 1.0)
            t_max = T * 1.1
            mean_spacing = 2.0 * math.pi / max(math.log(q * t_max / (2 * math.pi)), 1.0)
            n_t = int(math.ceil(t_max / (0.3 * mean_spacing)))
            print(f"\n--- q={q} K={K} t_max={t_max:.2f} n_t={n_t} ---")
            t_fft_runs = []
            t_dir_runs = []
            for _ in range(3):
                t0 = time.perf_counter()
                _, _ = find_zeros_all_chars(q, t_max=t_max, n_t=n_t,
                                              n_zeros_target=K, method="fft")
                t_fft_runs.append(time.perf_counter() - t0)
                t0 = time.perf_counter()
                _, _ = find_zeros_all_chars(q, t_max=t_max, n_t=n_t,
                                              n_zeros_target=K, method="direct")
                t_dir_runs.append(time.perf_counter() - t0)
            t_fft = float(np.median(t_fft_runs))
            t_dir = float(np.median(t_dir_runs))
            speedup = t_dir / t_fft if t_fft > 0 else float("nan")
            comparisons.append(dict(q=q, K=K, t_max=t_max, n_t=n_t,
                                     t_fft_s=t_fft, t_direct_s=t_dir,
                                     fft_speedup=speedup))
            print(f"  FFT method:    {t_fft*1e3:7.2f} ms")
            print(f"  DIRECT method: {t_dir*1e3:7.2f} ms")
            print(f"  FFT speedup:   {speedup:5.2f}x  (slot 2 measured 1.13-1.79x on L-eval primitive alone)")

        cmp_csv = os.path.join(HERE, "slot3_fft_vs_direct.csv")
        if comparisons:
            with open(cmp_csv, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(comparisons[0].keys()))
                w.writeheader()
                for r in comparisons:
                    w.writerow(r)
            print(f"\nwrote {cmp_csv}")

    # Headline comparison vs slot 1 baseline (q=1009, K=200)
    print()
    print("=" * 80)
    print("Slot 3 - Headline: comparison vs slot-1 baseline at q=1009, K=200, x=1e6")
    print("=" * 80)
    print(f"  slot-1 baseline (synthetic zeros):")
    print(f"    T_setup_q          = 165 ms")
    print(f"    T_full_per_x       =  22.7 ms")
    print(f"    single-query total = 187.7 ms")
    print(f"    amortised at M=256 =  23.4 ms")
    slot3_q1009_K200 = next((r for r in rows if r["q"] == 1009 and r["K_target"] == 200), None)
    if slot3_q1009_K200:
        zf = slot3_q1009_K200["t_zerofind_s"] * 1e3
        ps = slot3_q1009_K200["t_partialsum_s"] * 1e3
        single = slot3_q1009_K200["t_total_single_query_s"] * 1e3
        amort = zf / 256 + ps
        print(f"  slot-3 measured (real zeros via FFT-Hardy-Z):")
        print(f"    T_setup_q          = {zf:.1f} ms")
        print(f"    T_full_per_x       = {ps:.1f} ms")
        print(f"    single-query total = {single:.1f} ms ({187.7/single:.2f}x vs slot-1)")
        print(f"    amortised at M=256 = {amort:.1f} ms ({23.4/amort:.2f}x vs slot-1)")
        print()
        if single < 187.7 / 3:
            print("  >>> Slot-3 single-query speedup >= 3x: PARTIAL POSITIVE LIFTS <<<")
        elif single < 187.7:
            print(f"  >>> Slot-3 single-query speedup ~{187.7/single:.1f}x: bounded benefit. <<<")
        else:
            print("  >>> Slot-3 single-query speedup < 1x: FFT primitive does NOT lift end-to-end. <<<")


if __name__ == "__main__":
    main()
