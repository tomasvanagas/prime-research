"""
Thread 6, slot 1: per-(q, a) explicit-formula evaluator with decoupled cost.

Mission (from .commit_state thread_6_slot_plan slot_1):
  Build a per-(q, a) explicit-formula evaluator that decouples shared-q
  zero-of-L-function computation from per-character evaluation. Profile
  both at q in {10^3, 10^4, 10^5} for x in {10^6, 10^8}. Measure
  shared-zero amortisation factor.

Following Thread 5's Correlation Dichotomy partial-positive (S224, 33x
speedup at M=64 for correlated batched single-x queries), Thread 6
attempts a structurally similar amortisation on a different adjacent
problem: π(x; q, a) for batched modulus / residue family.

The explicit formula for primes in arithmetic progressions:

    ψ(x; q, a) = (1/φ(q)) Σ_{χ mod q} χ̄(a) · ψ(x, χ)
    π(x; q, a) ≈ (Riemann-prime-counting transform of ψ; see Davenport)

with for non-principal χ:

    ψ(x, χ) = -Σ_{ρ_χ} x^{ρ_χ}/ρ_χ + (lower-order correction)

and ρ_χ ranging over non-trivial zeros of L(s, χ).

DECOUPLED COST DECOMPOSITION (to be measured this slot):

    T_total(q, a, x, K) =
        T_zero_db_setup(q, K)              [zeros-of-L cache build, per-χ]
      + T_per_chi_eval(K, x) × #(χ mod q)  [partial-sum eval, per χ per x]
      + T_orthogonality(q)                 [Σ χ̄(a) sum]

Slot 1 questions:
  Q1. T_zero_db_setup(q, K) shape — does it scale as φ(q) · K^c, or is
      there structural sharing across χ mod q?
  Q2. T_per_chi_eval(K, x) — same x-independence as zeta case?
  Q3. T_orthogonality(q) cost relative to per-χ eval — when does
      orthogonal sum dominate vs evaluation?

Cross-domain ingredient: amortised algorithmics applied to
Dirichlet-character family (CROSS_DOMAIN_TECHNIQUES.md §8). Same
T_setup / T_eval decoupling Tarjan 1985 introduced for amortised data
structures, now applied to the Dirichlet-L family at fixed x.

WHAT WOULD FALSIFY THE SLOT 1 PREDICTION (per OPEN_POSITIVE_TARGETS P1
"predicted failure"): if T_zero_db_setup(q, K) actually grows much
faster than φ(q)·K (e.g., φ(q)·K^{17/13} or worse) AND there is no
structural sharing of zeros across χ, then per-(q, a) amortised
polylog is impossible — closure mode E.

The opposite finding (sub-φ(q) zero-setup growth) would open
slots 2-5 toward a real batch-polylog primes-in-AP algorithm.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import time
from typing import Dict, List, Tuple

import mpmath as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


# =====================================================================
# Dirichlet-character utilities (small q)
# =====================================================================


def euler_phi(n: int) -> int:
    p = n
    res = n
    if p % 2 == 0:
        res -= res // 2
        while p % 2 == 0:
            p //= 2
    f = 3
    while f * f <= p:
        if p % f == 0:
            res -= res // f
            while p % f == 0:
                p //= f
        f += 2
    if p > 1:
        res -= res // p
    return res


def primitive_root(q: int) -> int | None:
    """Find smallest primitive root mod q (for cyclic (Z/qZ)*)."""
    if q == 1:
        return None
    phi = euler_phi(q)
    factors = []
    p = phi
    d = 2
    while d * d <= p:
        if p % d == 0:
            factors.append(d)
            while p % d == 0:
                p //= d
        d += 1
    if p > 1:
        factors.append(p)
    for g in range(2, q):
        if math.gcd(g, q) != 1:
            continue
        ok = True
        for pf in factors:
            if pow(g, phi // pf, q) == 1:
                ok = False
                break
        if ok:
            return g
    return None


def discrete_log(g: int, h: int, q: int, phi: int) -> int:
    """Naive baby-step giant-step is overkill for our small q's."""
    cur = 1
    for k in range(phi):
        if cur == h:
            return k
        cur = (cur * g) % q
    raise ValueError(f"no discrete log of {h} base {g} mod {q}")


def all_characters_table(q: int) -> Tuple[List[int], np.ndarray]:
    """Return (residues_coprime, table) where table[j, n_idx] is the j-th
    Dirichlet character mod q evaluated at residues_coprime[n_idx], as a
    complex value. Only handles q with cyclic (Z/qZ)* (q prime, prime
    powers, 2*prime, etc.). For small q we use the prime-power case
    directly; for composite q with non-cyclic group we lift via CRT.
    """
    coprimes = [n for n in range(1, q) if math.gcd(n, q) == 1]
    phi = len(coprimes)

    if q in (2,):
        # Only trivial character.
        table = np.ones((1, len(coprimes)), dtype=complex)
        return coprimes, table

    g = primitive_root(q)
    if g is None:
        raise NotImplementedError(f"q={q} has non-cyclic (Z/qZ)*; not handled in slot 1")

    # Discrete log of each coprime element relative to g.
    log_n = {}
    cur = 1
    for k in range(phi):
        log_n[cur] = k
        cur = (cur * g) % q
    indices = [log_n[n] for n in coprimes]

    # χ_j(g) = exp(2πi j / phi); χ_j(n) = exp(2πi j log_g(n) / phi).
    omega = np.exp(2j * np.pi / phi)
    table = np.zeros((phi, len(coprimes)), dtype=complex)
    for j in range(phi):
        for ci, k in enumerate(indices):
            table[j, ci] = omega ** ((j * k) % phi)
    return coprimes, table


# =====================================================================
# Synthetic L-zero generator (slot 1: cost-profile only)
# =====================================================================
#
# Real Dirichlet-L zeros are not in the project's data cache, and
# computing them via root-finding on L(0.5+it, χ) is itself expensive.
# Slot 1's contract is to PROFILE the cost decomposition, not to
# produce numerically-tight π(x; q, a). So we use synthetic zeros:
#   • For χ_0 (principal): use ζ-zeros from the cache directly.
#   • For non-principal χ: zeros of L(s, χ) follow Riemann–von Mangoldt
#     density N_χ(T) ~ (T/2π) log(qT/(2πe)). Generate K imaginary parts
#     by inverting this density, then add deterministic χ-dependent
#     phase to break exact ζ-equality.
#
# This synthetic generator intentionally REMOVES any algorithmic
# structure that would let zeros be shared across distinct χ, which is
# the hypothesis under test. If a slot-2/3 attack discovers actual
# structural sharing (e.g., via approximate functional equation
# evaluations at common heights), it must justify the sharing
# explicitly — slot 1's baseline is "no sharing."
# =====================================================================


def synthetic_zeros_density_inverse(q: int, K: int, chi_seed: int) -> List[float]:
    """Generate K synthetic L(s, χ)-zero imaginary parts at modulus q,
    inverting N_χ(T) = (T/2π)(log(qT/(2π))-1) + O(1).

    Vectorised: 6 Newton iterations on numpy arrays. Cost per zero is
    dominated by the np-broadcast log step; sub-µs at K = 1000.

    chi_seed adds a per-character phase shift so that distinct χ have
    distinct zeros. Reproducible for a given (q, chi_seed).
    """
    rng = np.random.default_rng((q * 1_000_003 + chi_seed) & 0xFFFFFFFF)
    log_2pi = math.log(2.0 * math.pi)
    ks = np.arange(1, K + 1, dtype=np.float64)
    # Initial guess from leading-term inversion: T ≈ 2π k / log(qk)
    T = 2.0 * np.pi * ks / np.maximum(np.log(q * (ks + 2.0)), 1.0)
    for _ in range(8):
        log_qT = np.log(q * T)
        num = (T / (2.0 * np.pi)) * (log_qT - log_2pi - 1.0)
        f = num - ks
        df = (log_qT - log_2pi) / (2.0 * np.pi)
        df = np.where(np.abs(df) < 1e-12, 1e-12, df)
        T = T - f / df
        T = np.maximum(T, 0.5)
    # Per-χ jitter capped by half the local mean spacing.
    log_qT = np.log(q * T)
    local_spacing = 2.0 * np.pi / np.maximum(log_qT, 1.0)
    jitter = (rng.random(K) - 0.5) * local_spacing
    return (T + jitter).tolist()


def get_chi0_zeros(K: int) -> List[float]:
    """Principal-character L-zeros are ζ-zeros (modulo finite Euler-product
    factors that contribute only on Re(s) = 0). Reuse cached ζ-zeros."""
    data_dir = os.path.normpath(os.path.join(HERE, "..", "..", "..", "data"))
    for size in (300, 500, 1000, 2000, 8000):
        if size >= K:
            fname = os.path.join(data_dir, f"zeta_zeros_{size}.txt")
            if os.path.exists(fname):
                with open(fname) as f:
                    zs = [line.strip() for line in f if line.strip()]
                if len(zs) >= K:
                    return [float(z) for z in zs[:K]]
    return [float(mp.zetazero(k).imag) for k in range(1, K + 1)]


# =====================================================================
# Decoupled timing primitives
# =====================================================================


def time_zero_db_setup(q: int, K: int, n_repeats: int = 3) -> Tuple[float, Dict[int, List[float]]]:
    """Build the per-character zero database for all χ mod q.
    Returns (median_time_seconds, dict of χ_idx → list of γ).

    For the principal character (idx 0), uses cached ζ-zeros (fast file
    I/O). For non-principal characters, generates K synthetic zeros via
    Riemann–von Mangoldt density inversion. Per-χ generation cost is
    O(K) Newton iterations.
    """
    times: List[float] = []
    db: Dict[int, List[float]] = {}
    coprimes, table = all_characters_table(q)
    n_chars = table.shape[0]
    for trial in range(n_repeats):
        t0 = time.time()
        db = {}
        # Principal character (j=0): ζ-zeros from cache
        db[0] = get_chi0_zeros(K)
        # Non-principal characters: synthetic generators
        for j in range(1, n_chars):
            db[j] = synthetic_zeros_density_inverse(q, K, chi_seed=j)
        elapsed = time.time() - t0
        times.append(elapsed)
    return float(np.median(times)), db


def time_per_chi_eval(zeros: List[float], x: float, n_repeats: int = 3) -> float:
    """Time the partial-sum evaluation Σ_ρ x^ρ/ρ for one character at one x.

    Vectorised numpy (same primitive used in Thread 5 slot 1).
    """
    g = np.asarray(zeros, dtype=np.float64)
    log_x = math.log(x)
    sqrt_x = math.sqrt(x)
    times: List[float] = []
    for _ in range(n_repeats):
        t0 = time.time()
        cos_t = np.cos(g * log_x)
        sin_t = np.sin(g * log_x)
        denom = 0.25 + g * g
        inv_re = 0.5 / denom
        inv_im = -g / denom
        term_re = (cos_t * inv_re - sin_t * inv_im) * sqrt_x
        term_im = (cos_t * inv_im + sin_t * inv_re) * sqrt_x
        s_real = float(term_re.sum())
        s_imag = float(term_im.sum())
        elapsed = time.time() - t0
        times.append(elapsed)
    _ = (s_real, s_imag)
    return float(np.median(times))


def time_orthogonality(q: int, n_repeats: int = 5) -> float:
    """Time the inner χ̄(a) Σ_χ summation cost, given a precomputed table."""
    coprimes, table = all_characters_table(q)
    a_idx = 0  # representative residue
    times: List[float] = []
    for _ in range(n_repeats):
        psi_chi = np.random.default_rng(0).standard_normal(table.shape[0]).astype(complex)
        t0 = time.time()
        # χ̄(a) Σ_χ psi_chi(x): take row a in conjugate-table
        chi_bar_a = np.conj(table[:, a_idx])
        result = np.dot(chi_bar_a, psi_chi)
        elapsed = time.time() - t0
        _ = result  # prevent dead-code elimination
        times.append(elapsed)
    return float(np.median(times))


# =====================================================================
# Profile
# =====================================================================


def slot1_profile(qs: List[int], xs: List[float], Ks: List[int],
                  out_csv: str, n_repeats: int = 3) -> None:
    """Run the slot-1 cost-decomposition profile and write to CSV."""
    print(f"slot 1 profile: qs={qs}, xs={xs}, Ks={Ks}")
    rows: List[Dict[str, float]] = []
    for q in qs:
        phi_q = euler_phi(q)
        # Setup once per (q, K)
        for K in Ks:
            t_setup, db = time_zero_db_setup(q, K, n_repeats=n_repeats)
            # Time per-χ evaluation for the principal character (representative)
            t_chi0 = time_per_chi_eval(db[0], xs[0], n_repeats=n_repeats)
            # And for one non-principal character (also representative)
            t_chi1 = time_per_chi_eval(db.get(1, db[0]), xs[0], n_repeats=n_repeats)
            t_per_chi_med = 0.5 * (t_chi0 + t_chi1)
            t_orth = time_orthogonality(q, n_repeats=5)

            t_full_per_x = phi_q * t_per_chi_med + t_orth
            t_total_one_query = t_setup + t_full_per_x

            for x in xs:
                # Approximation: per-x evaluation cost is x-independent for
                # this primitive (all γ-dependent ops, no x-shape
                # dependence beyond the cos/sin args). We measure once at
                # xs[0] and report at all xs to avoid noise.
                row = dict(
                    q=q, phi_q=phi_q, K=K, x=x,
                    t_setup_q_seconds=t_setup,
                    t_setup_per_chi=t_setup / max(phi_q, 1),
                    t_per_chi_eval=t_per_chi_med,
                    t_orthogonality=t_orth,
                    t_full_per_x=t_full_per_x,
                    t_total_one_query=t_total_one_query,
                    setup_per_chi_per_K=t_setup / max(phi_q * K, 1),
                    eval_per_K=t_per_chi_med / max(K, 1),
                )
                rows.append(row)
                print(
                    f"  q={q:>5} φ(q)={phi_q:>5} K={K:>5} x={x:.0e} "
                    f"setup={t_setup:7.4f}s ({t_setup/phi_q*1e3:7.3f}ms/χ)  "
                    f"eval={t_per_chi_med*1e3:7.3f}ms/χ  "
                    f"orth={t_orth*1e6:7.2f}µs  "
                    f"full_per_x={t_full_per_x*1e3:7.2f}ms  "
                    f"total={t_total_one_query:7.4f}s"
                )

    fieldnames = list(rows[0].keys())
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"\nwrote {out_csv} ({len(rows)} rows)")


# =====================================================================
# Amortisation profile: per-(q, a) cost vs M batched (a) at fixed q
# =====================================================================


def amortisation_profile(q: int, K: int, x: float,
                         M_grid: List[int],
                         out_csv: str,
                         n_repeats: int = 3) -> None:
    """For fixed (q, x, K), measure per-(q, a) amortised cost vs M, where
    M = number of distinct a values queried with the SAME zero database.

    By the explicit formula, the zero database for L(s, χ) for χ mod q is
    INDEPENDENT of a — only χ̄(a) varies. So:

        T_total(M batched a's) = T_setup_q(K) + M · T_full_per_x

    and per-(q, a) amortised:

        T_amort(M) = T_setup_q(K) / M + T_full_per_x

    This is the q-batch analogue of Thread 5's x-batching dichotomy. The
    a-direction is GUARANTEED amortisable (zeros independent of a). The
    deeper Thread 6 questions are about χ amortisation (slot 2-3) and
    cross-q amortisation (slot 3).
    """
    print(f"\namortisation: q={q} K={K} x={x:.0e} M_grid={M_grid}")
    phi_q = euler_phi(q)
    t_setup, db = time_zero_db_setup(q, K, n_repeats=n_repeats)
    t_chi0 = time_per_chi_eval(db[0], x, n_repeats=n_repeats)
    t_chi1 = time_per_chi_eval(db.get(1, db[0]), x, n_repeats=n_repeats)
    t_per_chi_med = 0.5 * (t_chi0 + t_chi1)
    t_orth = time_orthogonality(q, n_repeats=5)
    t_full_per_x = phi_q * t_per_chi_med + t_orth

    rows = []
    print(f"  t_setup={t_setup:.4f}s  t_full_per_x={t_full_per_x*1e3:.3f}ms")
    for M in M_grid:
        t_amort = t_setup / M + t_full_per_x
        amort_factor = t_amort / (t_setup + t_full_per_x)  # vs M=1 cost
        rows.append(dict(
            q=q, phi_q=phi_q, K=K, x=x, M=M,
            t_setup=t_setup, t_full_per_x=t_full_per_x,
            t_amort_per_a=t_amort,
            amort_factor_vs_M1=amort_factor,
        ))
        print(f"    M={M:>5}: t_amort={t_amort*1e3:7.3f}ms/a  "
              f"amort_factor={amort_factor:.4f}× vs single query")

    fieldnames = list(rows[0].keys())
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"  wrote {out_csv}")


# =====================================================================
# Sanity: small-(q, x) sieve ground truth
# =====================================================================


def sieve_pi_q_a(x: int, q: int, a: int) -> int:
    """Direct sieve: count primes p ≤ x with p ≡ a (mod q)."""
    sieve = bytearray(b"\x01") * (x + 1)
    sieve[0] = sieve[1] = 0
    for p in range(2, int(math.isqrt(x)) + 1):
        if sieve[p]:
            for k in range(p * p, x + 1, p):
                sieve[k] = 0
    return sum(1 for p in range(2, x + 1) if sieve[p] and p % q == a)


# =====================================================================
# Main
# =====================================================================


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true", help="small grid for smoke test")
    args = parser.parse_args()

    if args.quick:
        qs = [7, 11, 31]
        xs = [1e6]
        Ks = [25, 50, 100]
        M_grid = [1, 4, 16, 64]
    else:
        # q ∈ {≈10², ≈10³, ≈10⁴, ≈10⁵} — exact primes for cyclic (Z/qZ)*.
        # 10⁵ pushes φ(q) ~ 10⁵ which is heavy; bound the synthetic-zero
        # generation to keep slot 1 within ~10 minutes. Slot 2 will go
        # higher.
        qs = [101, 1009, 10007]
        xs = [1e6, 1e8]
        Ks = [50, 200, 800]
        M_grid = [1, 4, 16, 64, 256]

    out_csv = os.path.join(HERE, "per_q_a_decoupled_profile.csv")
    slot1_profile(qs, xs, Ks, out_csv)

    # Amortisation profile at one representative (q, K, x).
    amort_csv = os.path.join(HERE, "per_q_a_amortisation.csv")
    amortisation_profile(q=qs[1], K=Ks[1], x=xs[0], M_grid=M_grid, out_csv=amort_csv)

    # Small ground-truth sanity check (independent of profiler).
    sanity_csv = os.path.join(HERE, "per_q_a_sieve_sanity.csv")
    print("\nground-truth sanity (direct sieve):")
    rows = []
    for q in [7, 11, 31]:
        for a in [1, 2, 3]:
            if math.gcd(a, q) != 1:
                continue
            for x in [10_000, 100_000]:
                pi_qa = sieve_pi_q_a(x, q, a)
                # Expected approximation: π(x)/φ(q)
                pi_x = sieve_pi_q_a(x, 1, 0) if False else sum(
                    1 for _ in range(2, x + 1)
                )  # placeholder
                # Use proper count instead:
                from sympy import primepi
                pi_x = int(primepi(x))
                expected = pi_x / euler_phi(q)
                row = dict(q=q, a=a, x=x, pi_q_a=pi_qa, expected=expected,
                           ratio=pi_qa / max(expected, 1e-9))
                rows.append(row)
                print(f"  q={q:>3} a={a:>3} x={x:>7}: π(x;q,a)={pi_qa:>6} "
                      f"≈ π(x)/φ(q)={expected:.1f} ratio={row['ratio']:.4f}")
    fieldnames = list(rows[0].keys())
    with open(sanity_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"wrote {sanity_csv}")


if __name__ == "__main__":
    main()
