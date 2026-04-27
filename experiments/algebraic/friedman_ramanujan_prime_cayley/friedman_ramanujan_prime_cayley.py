"""
D20 (ATTACK_VECTORS.md §D.D20) — Friedman / Ramanujan spectral gap of the
prime-generated abelian Cayley graph Cay(Z/NZ, S_N) where
S_N = {±p mod N : p prime, p < N^c}.

QUESTION (verbatim from D20):
For N prime, the abelian Cayley graph G_N = Cay(Z/NZ, S_N) of degree
d = |S_N| ≈ 2 N^c / (c log N) is diagonalised by characters: its
non-trivial eigenvalues are λ_k = Σ_{p < N^c} 2 cos(2π pk/N) for
k = 1, ..., N-1. Let λ_2(G_N) := max_{k ≠ 0} |λ_k| and define the
**Ramanujan ratio** r_N := λ_2 / (2 √(d-1)).

  - r_N ≈ 1   ⇔ Cayley graph is Ramanujan-typical (matches the
                Friedman 2008 random regular graph spectral gap bound).
  - r_N < 1  ⇔ super-Ramanujan: primes are a *better* expander than
                random regular graphs.
  - r_N > 1  ⇔ sub-Ramanujan: primes generate a *worse* expander than
                random regular graphs.

Note: Alon-Boppana gives λ_2 ≥ 2√(d-1) − o(1) as a deterministic lower
bound. So in the limit r_N ≥ 1 − o(1) for all sequences. The interesting
quantitative test is whether r_N(prime) sits at the Friedman-typical
value or above it (sub-Ramanujan).

Companion control: for each (N, c), draw L = 100 RANDOM SUBSETS of Z/NZ
of size matched to |{p : p < N^c}|, symmetrised, and compute their
λ_2 / (2 √(d-1)). Compare prime to the empirical random-subset
distribution via z-score Z_N = (r_N(prime) − mean(r_N(rand))) / σ(r_N(rand)).

PRE-DECLARED FALSIFICATION CRITERION (before running):

  A-grade: |Z_N| > 5 across all 5 N values for at least one c, with
           consistent sign (pure super- or sub-Ramanujan), AND the
           gap |r_N(prime) - r_N(random)| GROWS or stays nonzero as
           N → ∞ (slope on log-log plot has correct sign).

  B-grade (i): r_N(prime) ≈ r_N(random) within ±2σ across all (N, c)
               cells — primes are a Ramanujan-typical generating set,
               adds a 39th pseudorandomness measure (spectral gap of
               primes-as-generators).
  B-grade (ii): r_N(prime) deviates from r_N(random) at single (N, c)
               but with c-or-N-dependence inconsistent with a clean
               structural law.

  F-grade: experiment runs but produces only confirmation that
           λ_0 = |S| (trivial, equals π(N^c)), with no separation
           statistic on λ_2.

Cross-domain ingredient (PROMOTES "Random regular graph spectral gap
(Friedman)" §1 from PROPOSED to USED): Friedman 2008 Memoirs AMS 195
arXiv:cs/0405020. The non-trivial λ_k are exactly DFT magnitudes of the
indicator vector 1_{S_N}; the abelian Cayley spectrum is the discrete
Fourier transform of the symmetric generating set.

Channelling: BOURGAIN. The supremum is exactly the prime exponential
sum |Σ_{p ≤ M} e^{2π i p k / N}|; classical Vinogradov bounds give a
log-loss version of Ramanujan, but the empirical CONSTANT (Ramanujan
ratio) is open in the published literature.

Distinction from prior closures (verified before running):
  - CLOSED_PATHS line 356 (Cay(Z/xZ, primes), circular λ_0 = π(x)):
    closes the TRIVIAL eigenvalue. D20 is about the SECOND eigenvalue.
  - CLOSED_PATHS line 387 (GCD/coprimality graph spectrum = Ramanujan
    sums = Meissel-Lehmer ops): different graph, different question.
  - CLOSED_PATHS line 752 / E7.12 (S79, A3 fixed-generator
    Cay((Z/nZ)*, {2,3,5})): different group, different question
    (primality decision, not spectral gap saturation).
  - CLOSED_PATHS line 754 / E7.13 candidate (S80, D4 Szegedy walks
    on similar graphs): about quantum mixing time, not about classical
    second eigenvalue saturation of Friedman's random regular bound.

Files produced: results JSON + accompanying *_results.md.
"""

from __future__ import annotations

import json
import math
import os
import time
from dataclasses import dataclass, asdict
from typing import Sequence

import numpy as np
from sympy import isprime, primerange, sieve

OUT_DIR = os.path.dirname(os.path.abspath(__file__))


def primes_below(M: int) -> np.ndarray:
    """Return numpy array of all primes p with p < M."""
    return np.array(list(primerange(2, M)), dtype=np.int64)


def cayley_lambda2_via_fft(N: int, S_pos: np.ndarray) -> tuple[float, float, float, np.ndarray]:
    """
    Compute λ_k spectrum of the SYMMETRIC generating set
    S = S_pos ∪ -S_pos (mod N) via FFT on the {0,1}-indicator of S.

    Since S is symmetric, the eigenvalues are real. We return:
      λ_max     = max_{all k} λ_k    (= |S| at k = 0)
      λ_2       = max_{k != 0} |λ_k|
      λ_2_minor = max_{k : N/4 ≤ k ≤ 3N/4} |λ_k|
                  ("minor-arc" band — frequencies bounded away from 0,
                   the regime where Vinogradov's prime-exp-sum bound
                   applies; Friedman's Ramanujan reference is correct
                   here, not at k near 0).
      |λ| spectrum
    """
    indicator = np.zeros(N, dtype=np.float64)
    s_neg = (-S_pos) % N
    indicator[S_pos % N] = 1.0
    indicator[s_neg] = 1.0

    fft_vals = np.fft.fft(indicator)
    re = fft_vals.real
    abs_lambda = np.abs(re)

    lambda_max = abs_lambda[0]
    lambda_2 = float(abs_lambda[1:].max())

    # minor-arc band: N/4 ≤ k ≤ 3N/4
    klo, khi = N // 4, (3 * N) // 4
    lambda_2_minor = float(abs_lambda[klo:khi + 1].max())
    return float(lambda_max), float(lambda_2), float(lambda_2_minor), abs_lambda


def random_subset_lambda2(N: int, half_size: int, n_seeds: int, rng: np.random.Generator,
                          support_low: int = 1, support_high: int | None = None,
                          odd_only: bool = False, coprime_to: int | None = None) -> np.ndarray:
    """
    Random symmetric Cayley control: pick half_size elements uniformly
    at random WITHOUT replacement from [support_low, support_high), build
    S = S_pos ∪ -S_pos. Return array of λ_2 values across n_seeds.

    Two control modes:
      - support_high = None  (default = N): UNIFORM random subset of Z/NZ.
                              This is the "Friedman" reference.
      - support_high = M:    SUPPORT-MATCHED random subset of [2, M).
                              Isolates primes-vs-random within the same
                              concentrated support window — controls
                              away the trivial near-zero-frequency
                              alignment artefact (all generators small
                              ⇒ FFT spike near k = 0).
    """
    if support_high is None:
        support_high = N
    pool = np.arange(support_low, support_high, dtype=np.int64)
    if odd_only:
        pool = pool[pool % 2 == 1]
    if coprime_to is not None and coprime_to > 1:
        # gcd-mask via vectorised modular-arithmetic small-prime filter
        mask = np.ones(pool.shape, dtype=bool)
        # filter by each prime factor of coprime_to
        m = coprime_to
        for p in (2, 3, 5, 7, 11, 13):
            if m % p == 0:
                mask &= (pool % p != 0)
                while m % p == 0:
                    m //= p
        if m > 1:  # un-filtered remaining factor
            mask &= (np.gcd(pool, m) == 1)
        pool = pool[mask]
    if pool.size < half_size:
        # not enough elements to draw without replacement
        out = np.full(n_seeds, np.nan, dtype=np.float64)
        out_minor = np.full(n_seeds, np.nan, dtype=np.float64)
        return out, out_minor
    out = np.empty(n_seeds, dtype=np.float64)
    out_minor = np.empty(n_seeds, dtype=np.float64)
    for i in range(n_seeds):
        chosen = rng.choice(pool, size=half_size, replace=False)
        _, l2, l2_minor, _ = cayley_lambda2_via_fft(N, chosen)
        out[i] = l2
        out_minor[i] = l2_minor
    # stash minor-arc λ_2 in attribute on the array (numpy array can hold an attribute via .out)
    return out, out_minor


@dataclass
class Cell:
    N: int
    c: float
    M: int                      # = N^c (cutoff for primes)
    n_pos: int                  # |{p prime : p < M}|
    d: int                      # graph degree = 2 n_pos (assuming p != -p mod N)
    has_self_inverse: bool      # True if some p satisfies 2p ≡ 0 (mod N)
    lambda_max_prime: float     # should equal d
    lambda_2_prime: float
    lambda_2_minor_prime: float # restricted to k in [N/4, 3N/4]
    arg_k2_prime: int           # frequency k achieving λ_2 (full)
    ramanujan_bound: float      # 2 √(d-1)
    r_prime: float              # λ_2 / (2√(d-1))
    r_prime_minor: float        # λ_2_minor / (2√(d-1))
    # B1 = uniform random subset of Z/NZ \ {0}: full λ_2
    rand_uniform_lambda_2_mean: float
    rand_uniform_lambda_2_std: float
    rand_uniform_r_mean: float
    rand_uniform_r_std: float
    z_score_uniform_r: float
    # B2 = SUPPORT-MATCHED random subset of [2, M)
    rand_supp_lambda_2_mean: float
    rand_supp_lambda_2_std: float
    rand_supp_r_mean: float
    rand_supp_r_std: float
    z_score_supp_r: float
    # MINOR-ARC band r-statistics
    rand_uniform_r_minor_mean: float
    rand_uniform_r_minor_std: float
    z_score_uniform_r_minor: float
    rand_supp_r_minor_mean: float
    rand_supp_r_minor_std: float
    z_score_supp_r_minor: float
    # B3 = parity-matched: random subset of ODD ints in [3, M)
    rand_odd_lambda_2_mean: float
    rand_odd_lambda_2_std: float
    rand_odd_r_mean: float
    rand_odd_r_std: float
    z_score_odd_r: float
    rand_odd_r_minor_mean: float
    rand_odd_r_minor_std: float
    z_score_odd_r_minor: float
    # B4 = HL-W-trick-matched: random subset of [3, M), coprime to W=6
    rand_w6_r_mean: float
    rand_w6_r_std: float
    z_score_w6_r: float
    rand_w6_r_minor_mean: float
    rand_w6_r_minor_std: float
    z_score_w6_r_minor: float
    n_seeds: int
    runtime_seconds: float


def run_cell(N: int, c: float, n_seeds: int, rng: np.random.Generator, verbose: bool = True) -> Cell:
    t0 = time.time()
    M = max(2, int(round(N ** c)))
    S_pos = primes_below(M)
    # require S_pos to be coprime with N — for N prime this is automatic
    # (the only prime that could collide with N is N itself which is < N^c only if c >= 1).
    if not isprime(N):
        # fall back: drop any p that equals N (shouldn't happen since p < M < N for c < 1).
        S_pos = S_pos[S_pos < N]
    n_pos = int(S_pos.size)
    # check whether any 2p ≡ 0 (mod N): possible only if N is even, but N is prime and odd > 2.
    has_self_inverse = bool(((2 * S_pos) % N == 0).any())
    if has_self_inverse:
        # then |S| < 2 n_pos; rare edge case for our N values, but track it.
        s_set = set()
        for p in S_pos.tolist():
            s_set.add(p % N)
            s_set.add((-p) % N)
        d = len(s_set)
    else:
        d = 2 * n_pos
    lambda_max, lambda_2, lambda_2_minor, abs_spec = cayley_lambda2_via_fft(N, S_pos)
    arg_k2 = int(1 + abs_spec[1:].argmax())
    ram_bound = 2.0 * math.sqrt(max(d - 1, 1))
    r_prime = lambda_2 / ram_bound
    r_prime_minor = lambda_2_minor / ram_bound

    # B1: uniform random subset of Z/NZ \ {0}
    rand_uniform_l2, rand_uniform_l2_minor = random_subset_lambda2(
        N, n_pos, n_seeds, rng, support_low=1, support_high=N
    )
    rand_uniform_r = rand_uniform_l2 / ram_bound
    rand_uniform_r_minor = rand_uniform_l2_minor / ram_bound

    # B2: support-matched random subset of [2, M)
    rand_supp_l2, rand_supp_l2_minor = random_subset_lambda2(
        N, n_pos, n_seeds, rng, support_low=2, support_high=M
    )
    rand_supp_r = rand_supp_l2 / ram_bound
    rand_supp_r_minor = rand_supp_l2_minor / ram_bound

    # B3: parity-matched random subset of ODD ints in [3, M).
    # Primes (except 2) are all odd; this control kills the parity spike.
    rand_odd_l2, rand_odd_l2_minor = random_subset_lambda2(
        N, n_pos, n_seeds, rng, support_low=3, support_high=M, odd_only=True
    )
    rand_odd_r = rand_odd_l2 / ram_bound
    rand_odd_r_minor = rand_odd_l2_minor / ram_bound

    # B4: HL-W-trick-matched random subset of [3, M) coprime to W=6.
    # Tests whether the residual prime-vs-B3 deviation IS the next layer
    # of the Hardy-Littlewood sieve (mod-3 component).
    rand_w6_l2, rand_w6_l2_minor = random_subset_lambda2(
        N, n_pos, n_seeds, rng, support_low=3, support_high=M, coprime_to=6
    )
    rand_w6_r = rand_w6_l2 / ram_bound
    rand_w6_r_minor = rand_w6_l2_minor / ram_bound

    def _safe_z(x: float, mean: float, std: float) -> float:
        if std <= 0 or math.isnan(std):
            return float("nan")
        return float((x - mean) / std)

    rand_supp_valid = not np.isnan(rand_supp_l2).all()
    rand_odd_valid = not np.isnan(rand_odd_l2).all()
    rand_w6_valid = not np.isnan(rand_w6_l2).all()

    cell = Cell(
        N=N,
        c=c,
        M=M,
        n_pos=n_pos,
        d=d,
        has_self_inverse=has_self_inverse,
        lambda_max_prime=lambda_max,
        lambda_2_prime=lambda_2,
        lambda_2_minor_prime=lambda_2_minor,
        arg_k2_prime=arg_k2,
        ramanujan_bound=ram_bound,
        r_prime=r_prime,
        r_prime_minor=r_prime_minor,
        rand_uniform_lambda_2_mean=float(rand_uniform_l2.mean()),
        rand_uniform_lambda_2_std=float(rand_uniform_l2.std(ddof=1)),
        rand_uniform_r_mean=float(rand_uniform_r.mean()),
        rand_uniform_r_std=float(rand_uniform_r.std(ddof=1)),
        z_score_uniform_r=_safe_z(r_prime, rand_uniform_r.mean(), rand_uniform_r.std(ddof=1)),
        rand_supp_lambda_2_mean=float(rand_supp_l2.mean()) if rand_supp_valid else float("nan"),
        rand_supp_lambda_2_std=float(rand_supp_l2.std(ddof=1)) if rand_supp_valid else float("nan"),
        rand_supp_r_mean=float(rand_supp_r.mean()) if rand_supp_valid else float("nan"),
        rand_supp_r_std=float(rand_supp_r.std(ddof=1)) if rand_supp_valid else float("nan"),
        z_score_supp_r=_safe_z(r_prime, rand_supp_r.mean(), rand_supp_r.std(ddof=1))
        if rand_supp_valid else float("nan"),
        rand_uniform_r_minor_mean=float(rand_uniform_r_minor.mean()),
        rand_uniform_r_minor_std=float(rand_uniform_r_minor.std(ddof=1)),
        z_score_uniform_r_minor=_safe_z(r_prime_minor,
                                        rand_uniform_r_minor.mean(),
                                        rand_uniform_r_minor.std(ddof=1)),
        rand_supp_r_minor_mean=float(rand_supp_r_minor.mean()) if rand_supp_valid else float("nan"),
        rand_supp_r_minor_std=float(rand_supp_r_minor.std(ddof=1)) if rand_supp_valid else float("nan"),
        z_score_supp_r_minor=_safe_z(r_prime_minor,
                                     rand_supp_r_minor.mean(),
                                     rand_supp_r_minor.std(ddof=1)) if rand_supp_valid else float("nan"),
        rand_odd_lambda_2_mean=float(rand_odd_l2.mean()) if rand_odd_valid else float("nan"),
        rand_odd_lambda_2_std=float(rand_odd_l2.std(ddof=1)) if rand_odd_valid else float("nan"),
        rand_odd_r_mean=float(rand_odd_r.mean()) if rand_odd_valid else float("nan"),
        rand_odd_r_std=float(rand_odd_r.std(ddof=1)) if rand_odd_valid else float("nan"),
        z_score_odd_r=_safe_z(r_prime, rand_odd_r.mean(), rand_odd_r.std(ddof=1))
        if rand_odd_valid else float("nan"),
        rand_odd_r_minor_mean=float(rand_odd_r_minor.mean()) if rand_odd_valid else float("nan"),
        rand_odd_r_minor_std=float(rand_odd_r_minor.std(ddof=1)) if rand_odd_valid else float("nan"),
        z_score_odd_r_minor=_safe_z(r_prime_minor,
                                    rand_odd_r_minor.mean(),
                                    rand_odd_r_minor.std(ddof=1)) if rand_odd_valid else float("nan"),
        rand_w6_r_mean=float(rand_w6_r.mean()) if rand_w6_valid else float("nan"),
        rand_w6_r_std=float(rand_w6_r.std(ddof=1)) if rand_w6_valid else float("nan"),
        z_score_w6_r=_safe_z(r_prime, rand_w6_r.mean(), rand_w6_r.std(ddof=1))
        if rand_w6_valid else float("nan"),
        rand_w6_r_minor_mean=float(rand_w6_r_minor.mean()) if rand_w6_valid else float("nan"),
        rand_w6_r_minor_std=float(rand_w6_r_minor.std(ddof=1)) if rand_w6_valid else float("nan"),
        z_score_w6_r_minor=_safe_z(r_prime_minor,
                                   rand_w6_r_minor.mean(),
                                   rand_w6_r_minor.std(ddof=1)) if rand_w6_valid else float("nan"),
        n_seeds=n_seeds,
        runtime_seconds=float(time.time() - t0),
    )
    if verbose:
        print(
            f"[N={N:>6d}, c={c:.3f}] d={d:>4d}  "
            f"r_p={r_prime:.3f}  r_p^min={r_prime_minor:.3f}  "
            f"Z_odd(full)={cell.z_score_odd_r:+.2f}  Z_odd(min)={cell.z_score_odd_r_minor:+.2f}  "
            f"Z_w6(full)={cell.z_score_w6_r:+.2f}  Z_w6(min)={cell.z_score_w6_r_minor:+.2f}  "
            f"({cell.runtime_seconds:.1f}s)"
        )
    return cell


def diagnose_p2_artefact(N: int, c: float, n_seeds: int = 100,
                         rng: np.random.Generator | None = None) -> dict:
    """
    Diagnostic: check whether the prime-vs-B3 (odd-only) deviation is
    entirely explained by the single even prime p=2 contributing
    (-1)^2 = +1 at the parity frequency, breaking the all-odd alignment.

    Compares:
      - prime set (with p=2)
      - prime set MINUS p=2 (odd primes only)
    against B3 = random subset of ODD ints in [3, M).
    """
    if rng is None:
        rng = np.random.default_rng(20260427 + N)
    M = max(2, int(round(N ** c)))
    primes_full = primes_below(M)
    primes_no_2 = primes_full[primes_full > 2]
    n_pos_full = int(primes_full.size)
    n_pos_no2 = int(primes_no_2.size)

    _, l2_full, l2_min_full, _ = cayley_lambda2_via_fft(N, primes_full)
    _, l2_no2, l2_min_no2, _ = cayley_lambda2_via_fft(N, primes_no_2)

    d_full = 2 * n_pos_full
    d_no2 = 2 * n_pos_no2
    r_full = l2_full / (2.0 * math.sqrt(max(d_full - 1, 1)))
    r_full_min = l2_min_full / (2.0 * math.sqrt(max(d_full - 1, 1)))
    r_no2 = l2_no2 / (2.0 * math.sqrt(max(d_no2 - 1, 1)))
    r_no2_min = l2_min_no2 / (2.0 * math.sqrt(max(d_no2 - 1, 1)))

    # B3 controls matched to each cardinality
    odd_l2_full, odd_l2_min_full = random_subset_lambda2(
        N, n_pos_full, n_seeds, rng, support_low=3, support_high=M, odd_only=True
    )
    odd_l2_no2, odd_l2_min_no2 = random_subset_lambda2(
        N, n_pos_no2, n_seeds, rng, support_low=3, support_high=M, odd_only=True
    )
    odd_r_full = odd_l2_full / (2.0 * math.sqrt(max(d_full - 1, 1)))
    odd_r_full_min = odd_l2_min_full / (2.0 * math.sqrt(max(d_full - 1, 1)))
    odd_r_no2 = odd_l2_no2 / (2.0 * math.sqrt(max(d_no2 - 1, 1)))
    odd_r_no2_min = odd_l2_min_no2 / (2.0 * math.sqrt(max(d_no2 - 1, 1)))

    def safe_z(x, m, s):
        if s <= 0 or math.isnan(s):
            return float("nan")
        return float((x - m) / s)

    return {
        "N": N,
        "c": c,
        "M": M,
        "n_pos_full": n_pos_full,
        "n_pos_no2": n_pos_no2,
        "r_full": r_full,
        "r_full_min": r_full_min,
        "r_no2": r_no2,
        "r_no2_min": r_no2_min,
        "z_full_vs_odd": safe_z(r_full, odd_r_full.mean(), odd_r_full.std(ddof=1)),
        "z_full_min_vs_odd": safe_z(r_full_min, odd_r_full_min.mean(),
                                    odd_r_full_min.std(ddof=1)),
        "z_no2_vs_odd": safe_z(r_no2, odd_r_no2.mean(), odd_r_no2.std(ddof=1)),
        "z_no2_min_vs_odd": safe_z(r_no2_min, odd_r_no2_min.mean(),
                                   odd_r_no2_min.std(ddof=1)),
    }


def main(out_path: str | None = None):
    rng = np.random.default_rng(20260427)
    # Five primes N. Picked from D20's concrete first step.
    N_values = [509, 1009, 4001, 16001, 65537]
    # Two density exponents.
    c_values = [0.5, 2.0 / 3.0]
    n_seeds = 100

    cells: list[dict] = []
    for c in c_values:
        for N in N_values:
            cell = run_cell(N, c, n_seeds=n_seeds, rng=rng, verbose=True)
            cells.append(asdict(cell))

    summary = {
        "experiment": "friedman_ramanujan_prime_cayley",
        "attack_vector": "D20",
        "cross_domain": "Friedman 2008 Memoirs AMS 195 (random regular graph spectral gap) "
        "+ Bourgain Vinogradov prime exponential sums",
        "channelled_mathematician": "Bourgain",
        "N_values": N_values,
        "c_values": c_values,
        "n_seeds": n_seeds,
        "rng_seed": 20260427,
        "cells": cells,
    }

    if out_path is None:
        out_path = os.path.join(OUT_DIR, "friedman_ramanujan_prime_cayley.json")
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nwrote {out_path}")

    # full-band summary
    print("\n=== Summary FULL band (max over all k != 0) ===")
    hdr = (f"{'c':>5}  {'N':>6}  {'d':>5}  {'k_2*':>6}  {'r_prime':>8}"
           f"  {'r_unif':>14}  {'Z_unif':>7}  {'r_supp':>14}  {'Z_supp':>7}")
    print(hdr)
    for c in cells:
        print(
            f"{c['c']:>5.3f}  {c['N']:>6d}  {c['d']:>5d}  {c['arg_k2_prime']:>6d}  "
            f"{c['r_prime']:>8.4f}  "
            f"{c['rand_uniform_r_mean']:>7.4f}±{c['rand_uniform_r_std']:.4f}  "
            f"{c['z_score_uniform_r']:>+7.2f}  "
            f"{c['rand_supp_r_mean']:>7.4f}±{c['rand_supp_r_std']:.4f}  "
            f"{c['z_score_supp_r']:>+7.2f}"
        )

    print("\n=== Summary MINOR-ARC band (max over k in [N/4, 3N/4]) ===")
    print(f"{'c':>5}  {'N':>6}  {'d':>5}  {'r_p^min':>8}"
          f"  {'r_supp^min':>14}  {'Z_supp':>7}  {'r_odd^min':>14}  {'Z_odd':>7}")
    for c in cells:
        print(
            f"{c['c']:>5.3f}  {c['N']:>6d}  {c['d']:>5d}  "
            f"{c['r_prime_minor']:>8.4f}  "
            f"{c['rand_supp_r_minor_mean']:>7.4f}±{c['rand_supp_r_minor_std']:.4f}  "
            f"{c['z_score_supp_r_minor']:>+7.2f}  "
            f"{c['rand_odd_r_minor_mean']:>7.4f}±{c['rand_odd_r_minor_std']:.4f}  "
            f"{c['z_score_odd_r_minor']:>+7.2f}"
        )

    print("\n=== Sieve-cascade Z-scores FULL band: B3 = odd-only, B4 = coprime to 6 ===")
    print(f"{'c':>5}  {'N':>6}  {'d':>5}  {'r_p':>8}"
          f"  {'r_odd':>14}  {'Z_odd':>7}  {'r_w6':>14}  {'Z_w6':>7}")
    for c in cells:
        print(
            f"{c['c']:>5.3f}  {c['N']:>6d}  {c['d']:>5d}  "
            f"{c['r_prime']:>8.4f}  "
            f"{c['rand_odd_r_mean']:>7.4f}±{c['rand_odd_r_std']:.4f}  "
            f"{c['z_score_odd_r']:>+7.2f}  "
            f"{c['rand_w6_r_mean']:>7.4f}±{c['rand_w6_r_std']:.4f}  "
            f"{c['z_score_w6_r']:>+7.2f}"
        )

    print("\n=== Sieve-cascade Z-scores MINOR-ARC band ===")
    print(f"{'c':>5}  {'N':>6}  {'d':>5}  {'r_p^min':>8}"
          f"  {'r_odd^min':>14}  {'Z_odd':>7}  {'r_w6^min':>14}  {'Z_w6':>7}")
    for c in cells:
        print(
            f"{c['c']:>5.3f}  {c['N']:>6d}  {c['d']:>5d}  "
            f"{c['r_prime_minor']:>8.4f}  "
            f"{c['rand_odd_r_minor_mean']:>7.4f}±{c['rand_odd_r_minor_std']:.4f}  "
            f"{c['z_score_odd_r_minor']:>+7.2f}  "
            f"{c['rand_w6_r_minor_mean']:>7.4f}±{c['rand_w6_r_minor_std']:.4f}  "
            f"{c['z_score_w6_r_minor']:>+7.2f}"
        )

    # === Diagnostic: is the prime-vs-B3 deviation entirely the p=2 artefact?
    print("\n=== Diagnostic: prime vs prime-without-2, B3 odd-only control ===")
    print(f"{'c':>5}  {'N':>6}  {'r_full^min':>10}  {'Z_full^min':>10}  "
          f"{'r_no2^min':>10}  {'Z_no2^min':>10}")
    diag = []
    for c in c_values:
        for N in N_values:
            d = diagnose_p2_artefact(N, c, n_seeds=n_seeds, rng=rng)
            diag.append(d)
            print(
                f"{c:>5.3f}  {N:>6d}  {d['r_full_min']:>10.4f}  "
                f"{d['z_full_min_vs_odd']:>+10.2f}  "
                f"{d['r_no2_min']:>10.4f}  {d['z_no2_min_vs_odd']:>+10.2f}"
            )
    summary["diagnostic_p2_artefact"] = diag
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2)

    return summary


if __name__ == "__main__":
    main()
