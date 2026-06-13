"""
Thread 6, slot 2: AFE-shared L-value evaluation across characters mod q.

Mission (from .commit_state recommended_next_action, slot-1 falsifier α):
  Implement Approximate Functional Equation with shared 1/n^{1/2+it}
  table at fixed t. Measure per-character cost reduction vs slot-1
  baseline (22.7ms per query asymptote at q=1009, K=200, x=1e6).

The truncated AFE main term:

    L(½+it, χ) ≈ Σ_{n ≤ N} χ(n) / n^{½+it}

(plus a reflected sum with χ̄ and a gamma factor, which is structurally
identical for sharing purposes — adds a 2× constant).

Three methods compared (same output, different cost):

  M_DIRECT:   chi_at_n @ complex_pow          shape (n_chars, N_co) @ (N_co, n_t)
              cost ≈ phi(q) · N_co · n_t complex MACs (BLAS)

  M_AGGREGATE: aggregate W[r, j] = Σ_{n ≡ r (q)} 1/n^{½+it_j} (cost: N · n_t)
              then char_table @ W              shape (n_chars, phi) @ (phi, n_t)
              cost ≈ phi(q)² · n_t MACs

  M_FFT:       same aggregate W in log-g order (cost: N · n_t)
              then DFT over (Z/qZ)*: L = phi · ifft(W, axis=0)
              cost ≈ phi(q) · log(phi(q)) · n_t complex ops

When N_co > phi: M_AGGREGATE strictly better than M_DIRECT (smaller matmul).
When N_co < phi: M_AGGREGATE strictly worse than M_DIRECT (larger matmul).
M_FFT: theoretically beats both at scale, by factor (phi · N_co)/(N + phi·log phi).

For our regime (q in {1e2, 1e3, 1e4}, t_max in {50, 200}):
  q=101,   N=phi=100  → AGGREGATE = DIRECT, FFT might or might not win
  q=1009,  phi=1008, N≈130 → DIRECT < AGGREGATE; FFT vs DIRECT uncertain
  q=10007, phi=10006, N≈400 → DIRECT < AGGREGATE; FFT theory wins ~30×

Cross-domain ingredient (CROSS_DOMAIN_TECHNIQUES.md §10):
  FFT-based Dirichlet character transform (Bober-Hiary 2017 multi-character
  L-zero algorithm idea, generalised). Uses cyclic structure of (Z/qZ)*
  for prime q.

What would falsify the slot 2 prediction:
  (i) If FFT-based aggregation is strictly slower than direct BLAS
      matmul for q ≤ 1e4, then the AFE shared-table approach has no
      operational benefit and the slot-1 linear baseline holds. Closure
      mode E for option (α).
  (ii) If FFT is faster but the speedup is constant in q (e.g., always
       ~2x), then the benefit is a constant factor and not a sub-linear
       amortisation. Still a partial-positive but bounded.
  (iii) If FFT speedup grows with q (e.g., follows phi/log(phi) ratio),
        then real partial-positive: per-character L-eval cost amortises
        sub-linearly in phi(q). A-grade-shaped result for slot 2.

What this slot does NOT yet do:
  - Does not reduce the per-query partial-sum eval cost (slot 1's
    `T_per_chi_eval`) which remains O(K) per character.
  - Does not address whether the speedup propagates to FULL π(x;q,a)
    computation (would need the partial-sum step too, or zero-finding).
  - Slot 2's deliverable is the cost reduction on the L-eval primitive
    plus the architectural pattern. Slot 3 extends to zero-finding.
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

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from per_q_a_explicit_formula_profile import (  # noqa: E402
    euler_phi,
    primitive_root,
    all_characters_table,
)


# =====================================================================
# Setup utilities
# =====================================================================


def coprimes_in_log_g_order(q: int) -> Tuple[List[int], Dict[int, int]]:
    """Return (coprimes ordered as g^0, g^1, ..., g^{phi-1}, dict r → log_g(r))."""
    g = primitive_root(q)
    if g is None:
        raise ValueError(f"q={q} has non-cyclic (Z/qZ)*; not handled in slot 2")
    phi = euler_phi(q)
    seq: List[int] = []
    log_g: Dict[int, int] = {}
    cur = 1
    for k in range(phi):
        seq.append(cur)
        log_g[cur] = k
        cur = (cur * g) % q
    return seq, log_g


def char_table_at_residues(q: int, residue_list: List[int]) -> np.ndarray:
    """Build char_table[j, c] = χ_j(residue_list[c]) for j ∈ [0, φ).
    Vectorised: avoids the φ² Python loop in `all_characters_table`.
    For prime q only (cyclic (Z/qZ)*).
    """
    g = primitive_root(q)
    if g is None:
        raise ValueError(f"q={q} has non-cyclic (Z/qZ)*")
    phi = euler_phi(q)
    log_g: Dict[int, int] = {}
    cur = 1
    for k in range(phi):
        log_g[cur] = k
        cur = (cur * g) % q
    indices = np.array([log_g[r] for r in residue_list], dtype=np.int64)
    j_grid = np.arange(phi, dtype=np.int64)
    # phase_idx[j, c] = (j * indices[c]) mod phi
    phase_idx = (j_grid[:, None] * indices[None, :]) % phi
    return np.exp(2j * np.pi * phase_idx / phi)


def afe_truncation(q: int, t_max: float) -> int:
    """AFE truncation: N ≈ √(q · t_max / (2π))."""
    return max(int(math.ceil(math.sqrt(q * t_max / (2.0 * math.pi)))), 10)


def compute_complex_pow_all_n(t_grid: np.ndarray, N: int) -> np.ndarray:
    """Returns shape (N, n_t) complex matrix: cp[n-1, j] = 1/n^{1/2+it_j}.
    Single shared cos/sin computation."""
    n_arr = np.arange(1, N + 1, dtype=np.float64)
    log_n = np.log(n_arr)
    inv_sqrt_n = 1.0 / np.sqrt(n_arr)
    phase = -np.outer(log_n, t_grid)  # (N, n_t)
    return inv_sqrt_n[:, None] * np.exp(1j * phase)


# =====================================================================
# Three L-evaluation methods
# =====================================================================


def l_eval_direct(q: int, t_grid: np.ndarray, N: int) -> np.ndarray:
    """Method DIRECT: char_table @ complex_pow_at_coprimes.
    Output: L[k, j] = Σ_{n ≤ N, gcd(n,q)=1} χ_k(n) / n^{1/2+it_j}.

    Uses `char_table_at_residues` to build only the φ × N_co columns we
    need, avoiding the φ² baseline that all_characters_table builds."""
    n_arr = np.arange(1, N + 1)
    n_residues = (n_arr % q).astype(int)
    coprime_mask = np.array([math.gcd(int(r), q) == 1 and int(r) != 0 for r in n_residues])
    n_co = n_arr[coprime_mask]
    n_co_residues = (n_co % q).astype(int).tolist()

    log_n = np.log(n_co.astype(np.float64))
    inv_sqrt_n = 1.0 / np.sqrt(n_co.astype(np.float64))
    phase = -np.outer(log_n, t_grid)
    cp = inv_sqrt_n[:, None] * np.exp(1j * phase)  # (N_co, n_t)

    chi_at_n = char_table_at_residues(q, n_co_residues)  # (n_chars, N_co)
    return chi_at_n @ cp  # (n_chars, n_t)


def l_eval_aggregate_matmul(q: int, t_grid: np.ndarray, N: int) -> np.ndarray:
    """Method AGGREGATE: aggregate W[r, j] over n with n ≡ r (mod q),
    then char_table @ W. Uses full coprime list (φ × φ table — only run
    for small q where memory permits)."""
    coprimes = [n for n in range(1, q) if math.gcd(n, q) == 1]
    phi = len(coprimes)
    coprime_set = set(coprimes)
    coprime_idx_map = {c: i for i, c in enumerate(coprimes)}
    n_t = len(t_grid)

    cp_all = compute_complex_pow_all_n(t_grid, N)  # (N, n_t)

    W = np.zeros((phi, n_t), dtype=complex)
    for n_idx in range(N):
        r = (n_idx + 1) % q
        if r in coprime_set:
            W[coprime_idx_map[r]] += cp_all[n_idx]

    char_table = char_table_at_residues(q, coprimes)  # (n_chars, phi)
    return char_table @ W  # (n_chars, n_t)


def l_eval_fft(q: int, t_grid: np.ndarray, N: int) -> np.ndarray:
    """Method FFT: aggregate W in log-g order, then DFT over (Z/qZ)*.

    Indexing identity: char_table from all_characters_table is
        table[j, ci] = ω^(j · log_g(coprimes[ci]))   ω = exp(2πi/phi)
    so for any vector A indexed by coprimes:
        Σ_ci χ_j(coprimes[ci]) · A[ci] = Σ_k ω^(j · k) · A_logg[k]
    where A_logg is A reindexed to log-g order. Multiplying by phi · ifft
    along the residue axis gives this for all j at once.
    """
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

    # L[j, t] = Σ_k ω^(j·k) W[k, t] = phi · ifft(W along axis 0)[j, t]
    return phi * np.fft.ifft(W, axis=0)  # (phi, n_t)


# =====================================================================
# Cross-method correctness check
# =====================================================================


def check_methods_agree(q: int, N: int = 100, t_grid_size: int = 5,
                         t_max: float = 50.0, verbose: bool = True) -> Tuple[float, float]:
    """Run all three methods on the same input and report max pairwise error."""
    t_grid = np.linspace(2.0, t_max, t_grid_size)
    L_dir = l_eval_direct(q, t_grid, N)
    L_agg = l_eval_aggregate_matmul(q, t_grid, N)
    L_fft = l_eval_fft(q, t_grid, N)

    err_agg = float(np.max(np.abs(L_dir - L_agg)))
    err_fft = float(np.max(np.abs(L_dir - L_fft)))
    if verbose:
        max_abs = float(np.max(np.abs(L_dir)))
        print(f"  q={q} N={N} n_t={t_grid_size} t_max={t_max}")
        print(f"    max|L| ~ {max_abs:.3e}")
        print(f"    max|L_dir - L_agg| = {err_agg:.3e} (rel {err_agg/max_abs:.3e})")
        print(f"    max|L_dir - L_fft| = {err_fft:.3e} (rel {err_fft/max_abs:.3e})")
    return err_agg, err_fft


# =====================================================================
# Profile
# =====================================================================


def time_call(fn, *args, n_warmup: int = 1, n_repeats: int = 3) -> Tuple[float, np.ndarray]:
    """Time fn(*args). Returns (median wall-clock, last result)."""
    for _ in range(n_warmup):
        result = fn(*args)
    times: List[float] = []
    for _ in range(n_repeats):
        t0 = time.perf_counter()
        result = fn(*args)
        elapsed = time.perf_counter() - t0
        times.append(elapsed)
    return float(np.median(times)), result


def profile(qs: List[int], t_max_grid: List[float],
            n_t_grid: List[int], out_csv: str,
            n_repeats: int = 3,
            max_q_for_aggregate: int = 5000) -> None:
    """For q > max_q_for_aggregate, skip AGGREGATE_MATMUL (φ² memory)."""
    print(f"profile: qs={qs}, t_max_grid={t_max_grid}, n_t_grid={n_t_grid}")
    rows: List[Dict[str, float]] = []
    for q in qs:
        phi = euler_phi(q)
        run_agg = q <= max_q_for_aggregate
        for t_max in t_max_grid:
            N = afe_truncation(q, t_max)
            for n_t in n_t_grid:
                t_grid = np.linspace(2.0, t_max, n_t)

                t_dir, L_dir = time_call(l_eval_direct, q, t_grid, N, n_repeats=n_repeats)
                t_fft, L_fft = time_call(l_eval_fft, q, t_grid, N, n_repeats=n_repeats)
                if run_agg:
                    t_agg, L_agg = time_call(l_eval_aggregate_matmul, q, t_grid, N,
                                              n_repeats=n_repeats)
                    err_agg = float(np.max(np.abs(L_dir - L_agg)))
                    speedup_agg = t_dir / t_agg
                else:
                    t_agg, err_agg, speedup_agg = float("nan"), float("nan"), float("nan")

                err_fft = float(np.max(np.abs(L_dir - L_fft)))
                max_abs = float(np.max(np.abs(L_dir)))
                rel_err_fft = err_fft / max_abs if max_abs > 0 else 0.0

                speedup_fft = t_dir / t_fft

                n_l_evals = phi * n_t
                row = dict(
                    q=q, phi=phi, t_max=t_max, n_t=n_t, N_afe=N,
                    n_l_evals=n_l_evals,
                    t_direct_s=t_dir,
                    t_aggregate_s=t_agg,
                    t_fft_s=t_fft,
                    speedup_agg_vs_direct=speedup_agg,
                    speedup_fft_vs_direct=speedup_fft,
                    per_eval_direct_us=t_dir / n_l_evals * 1e6,
                    per_eval_fft_us=t_fft / n_l_evals * 1e6,
                    err_agg=err_agg,
                    err_fft=err_fft,
                    rel_err_fft=rel_err_fft,
                )
                rows.append(row)
                agg_str = (f"AGG={t_agg*1e3:7.2f}ms ({speedup_agg:5.2f}x)"
                           if run_agg else "AGG=    skip")
                print(
                    f"  q={q:>5} φ={phi:>5} t_max={t_max:>5.0f} n_t={n_t:>4} N={N:>4}: "
                    f"DIR={t_dir*1e3:7.2f}ms  {agg_str}  "
                    f"FFT={t_fft*1e3:7.2f}ms ({speedup_fft:5.2f}x)  "
                    f"err_fft={err_fft:.1e} (rel {rel_err_fft:.1e})"
                )

    if rows:
        fieldnames = list(rows[0].keys())
        with open(out_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            for r in rows:
                w.writerow(r)
        print(f"\nwrote {out_csv} ({len(rows)} rows)")


# =====================================================================
# Speedup-vs-q scaling profile (the slot-2 headline measurement)
# =====================================================================


def profile_scaling(qs: List[int], t_max: float, n_t: int,
                    out_csv: str, n_repeats: int = 3,
                    max_q_for_aggregate: int = 5000) -> None:
    """Hold (t_max, n_t) fixed, scan q. Measures whether speedup_fft scales
    with phi(q) (a-grade-shaped) or stays constant (bounded partial-positive)."""
    print(f"\n=== Scaling profile: t_max={t_max}, n_t={n_t} ===")
    rows: List[Dict[str, float]] = []
    for q in qs:
        phi = euler_phi(q)
        N = afe_truncation(q, t_max)
        t_grid = np.linspace(2.0, t_max, n_t)

        t_dir, L_dir = time_call(l_eval_direct, q, t_grid, N, n_repeats=n_repeats)
        t_fft, L_fft = time_call(l_eval_fft, q, t_grid, N, n_repeats=n_repeats)
        if q <= max_q_for_aggregate:
            t_agg, _ = time_call(l_eval_aggregate_matmul, q, t_grid, N, n_repeats=n_repeats)
            speedup_agg = t_dir / t_agg
        else:
            t_agg, speedup_agg = float("nan"), float("nan")

        err_fft = float(np.max(np.abs(L_dir - L_fft)))
        max_abs = float(np.max(np.abs(L_dir)))

        speedup_fft = t_dir / t_fft
        # Theoretical FFT speedup (asymptotic): (phi · N_co) / (N + phi · log phi)
        N_co = sum(1 for n in range(1, N + 1) if math.gcd(n, q) == 1)
        theoretical = (phi * N_co) / (N + phi * math.log2(phi))

        row = dict(
            q=q, phi=phi, N_afe=N, N_co=N_co, t_max=t_max, n_t=n_t,
            t_direct_ms=t_dir * 1e3,
            t_aggregate_ms=t_agg * 1e3,
            t_fft_ms=t_fft * 1e3,
            speedup_agg=speedup_agg,
            speedup_fft=speedup_fft,
            theoretical_speedup_fft=theoretical,
            err_fft=err_fft,
            rel_err_fft=err_fft / max_abs if max_abs > 0 else 0.0,
        )
        rows.append(row)
        print(
            f"  q={q:>5} φ={phi:>5} N={N:>4} N_co={N_co:>4}: "
            f"DIR={t_dir*1e3:7.2f}ms  AGG={t_agg*1e3:7.2f}ms ({speedup_agg:5.2f}x)  "
            f"FFT={t_fft*1e3:7.2f}ms ({speedup_fft:5.2f}x)  theory={theoretical:5.2f}x  "
            f"err={err_fft:.2e}"
        )

    if rows:
        fieldnames = list(rows[0].keys())
        with open(out_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            for r in rows:
                w.writerow(r)
        print(f"  wrote {out_csv}")


# =====================================================================
# End-to-end: zero-counting via |L|² local minima for one character
# (sanity check that the L-evaluator is producing meaningful values)
# =====================================================================


def count_zeros_via_min(L_values_one_char: np.ndarray, t_grid: np.ndarray,
                        threshold_factor: float = 0.1) -> int:
    """Count local minima of |L|² below threshold_factor × median(|L|²).
    This is a coarse zero-count proxy on the critical line."""
    abs2 = np.abs(L_values_one_char) ** 2
    median_val = float(np.median(abs2))
    threshold = threshold_factor * median_val
    n_zeros = 0
    for i in range(1, len(t_grid) - 1):
        if abs2[i] < threshold and abs2[i] <= abs2[i - 1] and abs2[i] <= abs2[i + 1]:
            n_zeros += 1
    return n_zeros


def riemann_von_mangoldt_count(q: int, T: float) -> float:
    """N_χ(T) ≈ (T / 2π) · (log(qT/(2π)) − 1)."""
    return T / (2 * math.pi) * (math.log(q * T / (2 * math.pi)) - 1.0)


def zero_count_sanity(q: int, t_max: float = 100.0, n_t: int = 8000) -> None:
    """For one non-principal character mod q, count approximate zeros via
    local minima of |L|²; compare to Riemann–von Mangoldt prediction."""
    print(f"\n=== Zero-count sanity at q={q}, t_max={t_max}, n_t={n_t} ===")
    N = afe_truncation(q, t_max)
    # Use a finer truncation for accuracy
    N = max(N, 5 * N)
    t_grid = np.linspace(2.0, t_max, n_t)
    L = l_eval_direct(q, t_grid, N)
    expected = riemann_von_mangoldt_count(q, t_max)

    # Per-character zero count from local minima
    counts = []
    for k in range(1, min(L.shape[0], 6)):  # skip principal (j=0)
        nz = count_zeros_via_min(L[k], t_grid)
        counts.append(nz)
    print(f"  N_afe={N}  expected per non-principal χ ≈ {expected:.1f} zeros in [0, {t_max}]")
    print(f"  measured (j=1..5): {counts}")
    print(f"  median: {sorted(counts)[len(counts)//2]}")


# =====================================================================
# Main
# =====================================================================


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true",
                        help="small grid for smoke test")
    args = parser.parse_args()

    print("=" * 80)
    print("Slot 2 — Cross-method correctness check (q=11, q=31, q=101)")
    print("=" * 80)
    for q in [11, 31, 101]:
        check_methods_agree(q, N=100, t_grid_size=20, t_max=50.0)
    print()

    print("=" * 80)
    print("Slot 2 — Cost profile across (q, t_max, n_t)")
    print("=" * 80)
    if args.quick:
        qs = [11, 101, 1009]
        t_max_grid = [50.0]
        n_t_grid = [200, 800]
    else:
        # Trimmed: q=10007 with large n_t blows out runtime (BLAS matmul
        # dominated). Two representative n_t values, two t_max values.
        qs = [101, 1009, 10007]
        t_max_grid = [50.0, 200.0]
        n_t_grid = [400, 800]

    out_csv = os.path.join(HERE, "slot2_l_eval_profile.csv")
    profile(qs, t_max_grid, n_t_grid, out_csv)

    # Scaling profile: vary q at fixed (t_max, n_t)
    print()
    print("=" * 80)
    print("Slot 2 — Scaling profile: speedup_fft vs q (headline)")
    print("=" * 80)
    if args.quick:
        scaling_qs = [11, 31, 101, 251, 503, 1009]
        scaling_csv = os.path.join(HERE, "slot2_scaling_quick.csv")
    else:
        scaling_qs = [101, 251, 503, 1009, 2003, 5003]
        scaling_csv = os.path.join(HERE, "slot2_scaling.csv")
    profile_scaling(scaling_qs, t_max=100.0, n_t=400, out_csv=scaling_csv)

    # Sanity check: zero count for small q
    print()
    print("=" * 80)
    print("Slot 2 — Zero-count sanity (verifying L-eval is meaningful)")
    print("=" * 80)
    zero_count_sanity(q=11, t_max=80.0, n_t=4000)
    zero_count_sanity(q=31, t_max=80.0, n_t=4000)


if __name__ == "__main__":
    main()
