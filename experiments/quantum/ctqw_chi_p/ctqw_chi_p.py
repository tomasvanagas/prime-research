"""ATTACK_VECTORS §D.D5 — Continuous-time quantum walk (CTQW) on the
divisor graph, with prime-indicator target state.

Cross-domain import: continuous-time quantum walks (Childs 2009 PRL 102,
180501 = arXiv:0806.1972; Childs-Cleve-Deotto-Farhi-Gutmann-Spielman 2003
STOC = arXiv:quant-ph/0209131). H = A_G (adjacency-matrix Hamiltonian),
unitary evolution U(t) = exp(-i H t), amplitude
P(t) = |<v_target | U(t) | v_seed>|^2.

D4 (S80) closed Szegedy walks on these graphs: spectral gap of the
discriminant matrix is 1/poly(x), so Szegedy mixing is poly(x) — no
polylog π(x) extraction. CTQW is known to admit *exponential* speedup
over Szegedy on glued binary trees (Childs et al. 2003), driven by
spectral DENSITY and eigenvector overlap rather than gap. So the D4
closure does not transfer directly: we measure the CTQW amplitude on
the prime subset of D_x, and compare to (a) the equilibrium baseline
π(x)/x (a uniform target gives this trivially) and (b) controls that
randomise either the divisor structure or the primality indicator.

Specifically:
  • H = A_{D_x} where D_x is the divisor graph on {1, ..., x}: edges
    (m, n) iff m|n and m != n.
  • v_seed = |1>.
  • v_prime = (1 / sqrt(π(x))) sum_{p ≤ x prime} |p>.
  • P(t) = |<v_prime | exp(-iHt) | 1>|^2.
  • Sweep t ∈ [0, 100], find max P(t) and its peak location.

Controls:
  C1 (target randomisation): replace v_prime with a random Bernoulli
      indicator of size π(x). 200 seeds, report z-score.
  C2 (target randomisation, more stringent): match Cramér + odd-parity
      density (Cramér 1/log n + force odd, plus single even '2').
  C3 (graph randomisation): keep v_prime; replace the divisor graph with
      a random regular graph with matching mean degree. 200 seeds.

A-grade success: identify a (graph, seed) where P(t*) is polylog-
detectable, AND the cost of evaluating |1>(t*) in the prime basis is
polylog (Hamiltonian simulation overhead).
B-grade success: structurally extend D4 (Szegedy closure) to CTQW via a
clean lemma "Szegedy gap << 1/polylog ⇒ CTQW mixing >> polylog" on
the divisor graph family — promotes E7.13 to cover both walk types.
B-grade negative shape (the most likely outcome): max-amplitude
matches the trivial π(x)/x equilibrium within ±2σ across all
controls, closing CTQW the same way Szegedy was closed but via a
spectral-density argument rather than a gap one.

Falsification statement: if max_t P(t) > c * π(x)/x for some non-trivial
constant c > 1.5 with z-score > 3σ stable across x ∈ {64, ..., 512}
AND eigendecomposition reveals an isolated band-edge cluster, the D4
closure does NOT transfer to CTQW. Otherwise, the closure transfers.

Edges: cites E7.13 (D4 Szegedy closure), composes with E7.16 (Friedman
density-and-parity-matched controls).
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from sympy import isprime, primerange


# ---------------------------------------------------------------------------
# Graph constructors
# ---------------------------------------------------------------------------

def divisor_graph(x: int) -> np.ndarray:
    """A_{D_x}: vertices [1..x], edges (m, n) iff m|n and m != n."""
    A = np.zeros((x, x), dtype=np.float64)
    for m in range(1, x + 1):
        for k in range(2, x // m + 1):
            j = m * k
            A[m - 1, j - 1] = 1.0
            A[j - 1, m - 1] = 1.0
    return A


def coprime_graph(x: int) -> np.ndarray:
    """A_{C_x}: vertices [1..x], edges (m, n) iff gcd(m, n) = 1 and m != n."""
    A = np.zeros((x, x), dtype=np.float64)
    for m in range(1, x + 1):
        for n in range(m + 1, x + 1):
            if math.gcd(m, n) == 1:
                A[m - 1, n - 1] = 1.0
                A[n - 1, m - 1] = 1.0
    return A


def random_regular_graph_adj(N: int, mean_degree: float, rng: np.random.Generator) -> np.ndarray:
    """Erdős-Rényi G(N, p) with p = mean_degree/(N-1). Returns symmetric float adjacency."""
    p = mean_degree / max(1, N - 1)
    p = min(max(p, 0.0), 1.0)
    M = rng.random((N, N)) < p
    M = np.triu(M, 1)
    A = (M | M.T).astype(np.float64)
    np.fill_diagonal(A, 0.0)
    return A


# ---------------------------------------------------------------------------
# CTQW core
# ---------------------------------------------------------------------------

def ctqw_amplitude_curve(H: np.ndarray, v_seed: np.ndarray, v_target: np.ndarray,
                         t_grid: np.ndarray) -> np.ndarray:
    """Compute P(t) = |<v_target | exp(-iHt) | v_seed>|^2 across t_grid.

    Eigendecomposition: H = U diag(lam) U^T; <v_t | exp(-iHt) | v_s> =
    sum_k <v_t|u_k><u_k|v_s> exp(-i lam_k t).
    """
    lam, U = eigh(H)
    return ctqw_amplitude_curve_eig(lam, U, v_seed, v_target, t_grid)


def ctqw_amplitude_curve_eig(lam: np.ndarray, U: np.ndarray,
                             v_seed: np.ndarray, v_target: np.ndarray,
                             t_grid: np.ndarray,
                             phase_cache: np.ndarray | None = None) -> np.ndarray:
    """Same as ctqw_amplitude_curve but with precomputed eigh and optional phase cache."""
    a = U.T @ v_seed
    b = U.T @ v_target
    coef = np.conj(b) * a
    if phase_cache is None:
        phase_cache = np.exp(-1j * np.outer(t_grid, lam))
    amp = phase_cache @ coef
    return np.abs(amp) ** 2


def prime_indicator(x: int) -> np.ndarray:
    v = np.zeros(x, dtype=np.float64)
    cnt = 0
    for p in primerange(2, x + 1):
        v[p - 1] = 1.0
        cnt += 1
    if cnt > 0:
        v /= math.sqrt(cnt)
    return v


def random_subset_indicator(x: int, k: int, rng: np.random.Generator) -> np.ndarray:
    """Uniform random subset of size k from [2..x]; vertex 1 always excluded
    (matches the prime support [2..x])."""
    pool = np.arange(2, x + 1)  # 1-indexed
    sel = rng.choice(pool, size=k, replace=False)
    v = np.zeros(x, dtype=np.float64)
    for s in sel:
        v[s - 1] = 1.0
    v /= math.sqrt(k)
    return v


def cramer_odd_indicator(x: int, target_count: int, rng: np.random.Generator) -> np.ndarray:
    """Cramér 1/log n + odd-parity matched indicator on [2..x].

    Bernoulli with probability 1/log n at each n ≥ 2 (forced odd, except n=2
    is included with prob 1 to match parity to {primes}: exactly one even
    member, the rest odd). Then *force* the count to target_count (resample
    if mismatch up to 50 retries; else trim/extend uniformly within odd pool).
    """
    for _ in range(50):
        chosen = []
        # Match parity: include 2 (the single even prime) deterministically.
        chosen.append(2)
        for n in range(3, x + 1):
            if n % 2 == 0:
                continue
            p = 1.0 / math.log(n) if n >= 3 else 0.5
            if rng.random() < p:
                chosen.append(n)
        if len(chosen) == target_count:
            break
    if len(chosen) != target_count:
        odds = [n for n in range(3, x + 1, 2)]
        rng.shuffle(odds)
        chosen = [2] + odds[: target_count - 1]
    v = np.zeros(x, dtype=np.float64)
    for s in chosen:
        v[s - 1] = 1.0
    v /= math.sqrt(len(chosen))
    return v


# ---------------------------------------------------------------------------
# Experiment driver
# ---------------------------------------------------------------------------

def run_one(x: int, t_grid: np.ndarray, n_seeds: int, rng: np.random.Generator,
            graph_kind: str = "divisor"):
    if graph_kind == "divisor":
        A = divisor_graph(x)
    elif graph_kind == "coprime":
        A = coprime_graph(x)
    else:
        raise ValueError(graph_kind)

    H = A.copy()
    e1 = np.zeros(x, dtype=np.float64)
    e1[0] = 1.0

    v_p = prime_indicator(x)
    pi_x = sum(1 for _ in primerange(2, x + 1))

    # One eigendecomposition for the main graph; reuse for all controls C1/C2
    lam, U = eigh(H)
    phase_cache = np.exp(-1j * np.outer(t_grid, lam))

    # Primary measurement
    P_prime = ctqw_amplitude_curve_eig(lam, U, e1, v_p, t_grid, phase_cache)
    max_prime = float(np.max(P_prime))
    tstar_prime = float(t_grid[int(np.argmax(P_prime))])
    mean_prime = float(np.mean(P_prime))
    equilib = pi_x / x  # classical equilibrium

    # Spectral analysis
    overlap_v_p_eig = U.T @ v_p   # |<u_k | v_p>|
    # top-5 eigenvectors by overlap with v_p
    abs_overlap = np.abs(overlap_v_p_eig)
    top_idx = np.argsort(abs_overlap)[::-1][:10]
    top_eig_data = [
        {"k": int(i), "lambda": float(lam[i]), "abs_overlap_v_p": float(abs_overlap[i])}
        for i in top_idx
    ]
    band_edge_lam = [float(lam[0]), float(lam[1]), float(lam[-2]), float(lam[-1])]

    # Controls (reuse cached eigh of H)
    def baseline_run(indicator_fn, n_seeds_):
        peaks, means = [], []
        for s in range(n_seeds_):
            v = indicator_fn(s)
            P = ctqw_amplitude_curve_eig(lam, U, e1, v, t_grid, phase_cache)
            peaks.append(float(np.max(P)))
            means.append(float(np.mean(P)))
        peaks = np.array(peaks)
        means = np.array(means)
        return {
            "peak_mean": float(peaks.mean()),
            "peak_std": float(peaks.std(ddof=1) if len(peaks) > 1 else 0.0),
            "mean_mean": float(means.mean()),
            "mean_std": float(means.std(ddof=1) if len(means) > 1 else 0.0),
        }

    rng_b = np.random.default_rng(int(rng.integers(1 << 30)))
    C1 = baseline_run(lambda s, rng_b=rng_b: random_subset_indicator(x, pi_x, rng_b), n_seeds)
    rng_c = np.random.default_rng(int(rng.integers(1 << 30)))
    C2 = baseline_run(lambda s, rng_c=rng_c: cramer_odd_indicator(x, pi_x, rng_c), n_seeds)

    # C3: random graph control (different graph, primes target).
    mean_deg_div = float(A.sum() / x)
    rng_g = np.random.default_rng(int(rng.integers(1 << 30)))
    g_peaks = []
    g_means = []
    for _ in range(n_seeds):
        Ag = random_regular_graph_adj(x, mean_deg_div, rng_g)
        Pg = ctqw_amplitude_curve(Ag, e1, v_p, t_grid)
        g_peaks.append(float(np.max(Pg)))
        g_means.append(float(np.mean(Pg)))
    g_peaks = np.array(g_peaks)
    g_means = np.array(g_means)
    C3 = {
        "peak_mean": float(g_peaks.mean()),
        "peak_std": float(g_peaks.std(ddof=1) if len(g_peaks) > 1 else 0.0),
        "mean_mean": float(g_means.mean()),
        "mean_std": float(g_means.std(ddof=1) if len(g_means) > 1 else 0.0),
    }

    def zscore(measured, ctrl):
        mu = ctrl["peak_mean"]
        sd = ctrl["peak_std"] if ctrl["peak_std"] > 0 else 1e-300
        return float((measured - mu) / sd)

    def zscore_mean(measured, ctrl):
        mu = ctrl["mean_mean"]
        sd = ctrl["mean_std"] if ctrl["mean_std"] > 0 else 1e-300
        return float((measured - mu) / sd)

    return {
        "x": x,
        "graph_kind": graph_kind,
        "pi_x": pi_x,
        "equilib_pi_x_over_x": float(equilib),
        "max_amp_primes": max_prime,
        "t_star": tstar_prime,
        "mean_amp_primes": mean_prime,
        "ratio_max_over_equilib": float(max_prime / equilib),
        "ratio_mean_over_equilib": float(mean_prime / equilib),
        "C1_random_subset": C1,
        "C2_cramer_odd": C2,
        "C3_random_graph": C3,
        "z_peak_C1": zscore(max_prime, C1),
        "z_peak_C2": zscore(max_prime, C2),
        "z_peak_C3": zscore(max_prime, C3),
        "z_mean_C1": zscore_mean(mean_prime, C1),
        "z_mean_C2": zscore_mean(mean_prime, C2),
        "z_mean_C3": zscore_mean(mean_prime, C3),
        "band_edge_lam": band_edge_lam,
        "top10_eigvec_overlap_v_p": top_eig_data,
        "mean_degree": mean_deg_div,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xs", type=int, nargs="+", default=[64, 128, 256, 512])
    ap.add_argument("--T-max", type=float, default=100.0)
    ap.add_argument("--T-pts", type=int, default=4001)
    ap.add_argument("--seeds", type=int, default=100)
    ap.add_argument("--graph", type=str, default="divisor",
                    choices=["divisor", "coprime"])
    ap.add_argument("--out", type=str, default="ctqw_chi_p_results.json")
    ap.add_argument("--seed", type=int, default=20260427)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    t_grid = np.linspace(0.0, args.T_max, args.T_pts)

    out = {
        "args": vars(args),
        "T_max": args.T_max,
        "T_pts": args.T_pts,
        "seeds": args.seeds,
        "graph": args.graph,
        "results": [],
    }
    for x in args.xs:
        t0 = time.time()
        r = run_one(x, t_grid, args.seeds, rng, graph_kind=args.graph)
        r["wall_seconds"] = time.time() - t0
        print(f"x = {x:>4d}  pi_x = {r['pi_x']:>4d}  "
              f"max_amp = {r['max_amp_primes']:.5f}  "
              f"equilib = {r['equilib_pi_x_over_x']:.5f}  "
              f"ratio = {r['ratio_max_over_equilib']:.3f}  "
              f"z(C1) = {r['z_peak_C1']:+.2f}  "
              f"z(C2) = {r['z_peak_C2']:+.2f}  "
              f"z(C3) = {r['z_peak_C3']:+.2f}  "
              f"({r['wall_seconds']:.1f}s)")
        out["results"].append(r)

    Path(args.out).write_text(json.dumps(out, indent=2))
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
