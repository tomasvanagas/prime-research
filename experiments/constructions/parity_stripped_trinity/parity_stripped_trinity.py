#!/usr/bin/env python3
"""
Parity-stripped trinity (S239, paradigm-shift mode).

Builds chi_tilde_P(n) := chi_P(n) - alpha_2 * (-1)^n where
alpha_2 = mean_n chi_P(n) (-1)^n; this is the orthogonal projection of
chi_P onto {x : <x, (-1)^.> = 0} along the parity vector.

Measures three quantities for both chi_P and chi_tilde_P:
  - Mahler measure log m(f_N)            (E2.20 territory)
  - L^infty norm and arg-max location    (E2.21 territory)
  - Toeplitz 4th spectral moment m_4     (E2.31 territory)

Compares to matched-density Bernoulli baselines.

Run:
    python parity_stripped_trinity.py --N_log2 14,16,18 --bdj_N 500,1000,1500
"""
from __future__ import annotations

import argparse
import json
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np


# -- Sieve ----------------------------------------------------------------

def chi_p_indicator(N: int) -> np.ndarray:
    """Return chi_P(1..N) as a length-N float64 array (1.0 at primes, 0.0 else)."""
    sieve = np.ones(N + 1, dtype=np.bool_)
    sieve[0] = False
    if N >= 1:
        sieve[1] = False
    for p in range(2, int(N ** 0.5) + 1):
        if sieve[p]:
            sieve[p * p :: p] = False
    chi = np.zeros(N, dtype=np.float64)
    chi[:] = sieve[1 : N + 1].astype(np.float64)
    return chi


def parity_vec(N: int) -> np.ndarray:
    """(-1)^n for n=1..N as float64."""
    v = np.empty(N, dtype=np.float64)
    v[0::2] = -1.0  # n=1, 3, 5, ...
    v[1::2] = +1.0  # n=2, 4, 6, ...
    return v


def parity_strip(x: np.ndarray) -> tuple[np.ndarray, float]:
    """Subtract the (-1)^n component of x. Returns (x_tilde, alpha)."""
    N = x.size
    pv = parity_vec(N)
    alpha = float(np.dot(x, pv) / N)
    return x - alpha * pv, alpha


# -- Mahler via Jensen-FFT ------------------------------------------------

def log_mahler(coeffs: np.ndarray, M_factor: int = 4, M_min: int = 4096) -> float:
    """log m(f) where f(z) = sum_{n=0..N-1} coeffs[n] z^n."""
    N = coeffs.size
    M = max(N * M_factor, M_min)
    pad = np.zeros(M, dtype=coeffs.dtype)
    pad[:N] = coeffs
    F = np.abs(np.fft.fft(pad))
    F = np.maximum(F, 1e-300)
    return float(np.mean(np.log(F)))


# -- L^infty + arg-max ----------------------------------------------------

def linfty_with_argmax(coeffs: np.ndarray, M_factor: int = 4, M_min: int = 4096) -> tuple[float, int, int]:
    """|f(e^{2pi i k / M})| at M = max(N * M_factor, M_min). Returns (Linfty, k_argmax, M)."""
    N = coeffs.size
    M = max(N * M_factor, M_min)
    pad = np.zeros(M, dtype=coeffs.dtype)
    pad[:N] = coeffs
    F = np.abs(np.fft.fft(pad))
    k = int(np.argmax(F))
    return float(F[k]), k, M


def evaluate_at_root_of_unity(coeffs: np.ndarray, a: int, q: int) -> complex:
    """f(e^{2 pi i a / q}) by direct sum (small q only, exact phase)."""
    N = coeffs.size
    n = np.arange(1, N + 1)
    return complex(np.sum(coeffs * np.exp(2j * np.pi * a * n / q)))


# -- Toeplitz m_4 ---------------------------------------------------------

def toeplitz_m4(eps: np.ndarray) -> tuple[float, float]:
    """4th spectral moment of T_{ij} = eps[|i-j|+1].

    Returns (m_4, lambda_max_scaled) using BDJ scaling lambda / sqrt(N)."""
    N = eps.size
    # Build symmetric Toeplitz matrix with first row eps[0..N-1]
    # T[i,j] = eps[|i-j|]   (note: abs index is 0..N-1)
    # We use eps directly without offset, matching S204's convention.
    idx = np.abs(np.arange(N)[:, None] - np.arange(N)[None, :])
    T = eps[idx]
    eigs = np.linalg.eigvalsh(T)
    eigs_scaled = eigs / np.sqrt(N)
    m4 = float(np.mean(eigs_scaled ** 4))
    lambda_max = float(np.max(np.abs(eigs_scaled)))
    return m4, lambda_max


def standardise(x: np.ndarray) -> np.ndarray:
    p = float(np.mean(x))
    var = float(np.mean((x - p) ** 2))
    if var <= 0:
        raise ValueError("Zero variance — cannot standardise.")
    return (x - p) / np.sqrt(var)


# -- Bernoulli baseline ----------------------------------------------------

def bernoulli_match_density(N: int, target_count: int, rng: np.random.Generator) -> np.ndarray:
    """Random 0/1 sequence of length N with EXACTLY target_count ones."""
    arr = np.zeros(N, dtype=np.float64)
    idx = rng.choice(N, size=target_count, replace=False)
    arr[idx] = 1.0
    return arr


# -- Main experiment ------------------------------------------------------

@dataclass
class FFTResult:
    N: int
    log_m_chi: float
    log_m_chi_tilde: float
    linfty_chi: float
    linfty_chi_tilde: float
    argmax_q_chi: int
    argmax_a_chi: int
    argmax_q_tilde: int
    argmax_a_tilde: int
    M_grid: int
    sqrt_pi_N: float
    pi_N: int
    alpha_2: float
    f_at_minus_one_chi: complex
    f_at_minus_one_tilde: complex
    log_m_bern_mean: float
    log_m_bern_std: float
    log_m_bern_strip_mean: float
    log_m_bern_strip_std: float
    bern_replicates: int


@dataclass
class BDJResult:
    N: int
    m4_chi: float
    m4_chi_tilde: float
    lambda_max_chi: float
    lambda_max_chi_tilde: float
    m4_bern_mean: float
    m4_bern_std: float
    m4_bern_strip_mean: float
    m4_bern_strip_std: float
    bern_replicates: int


def closest_rational(k: int, M: int, q_max: int = 64) -> tuple[int, int, float]:
    """Best (a, q) with q <= q_max approximating k/M. Returns (a, q, |k/M - a/q|)."""
    target = k / M
    best = (0, 1, abs(target))
    for q in range(1, q_max + 1):
        a = int(round(target * q))
        a %= q
        diff = abs(target - a / q)
        if diff < best[2]:
            best = (a, q, diff)
    return best


def run_fft_block(N: int, n_bern: int, M_factor: int, rng: np.random.Generator) -> FFTResult:
    chi = chi_p_indicator(N)
    pi_N = int(chi.sum())
    chi_tilde, alpha_2 = parity_strip(chi)

    log_m_chi = log_mahler(chi, M_factor=M_factor)
    log_m_tilde = log_mahler(chi_tilde, M_factor=M_factor)
    linfty_chi, k_chi, M = linfty_with_argmax(chi, M_factor=M_factor)
    linfty_tilde, k_tilde, _ = linfty_with_argmax(chi_tilde, M_factor=M_factor)

    a_chi, q_chi, _ = closest_rational(k_chi, M)
    a_tilde, q_tilde, _ = closest_rational(k_tilde, M)

    f_minus_one_chi = evaluate_at_root_of_unity(chi, 1, 2)
    f_minus_one_tilde = evaluate_at_root_of_unity(chi_tilde, 1, 2)

    bern_log_m = []
    bern_log_m_strip = []
    for _ in range(n_bern):
        b = bernoulli_match_density(N, pi_N, rng)
        bern_log_m.append(log_mahler(b, M_factor=M_factor))
        b_tilde, _ = parity_strip(b)
        bern_log_m_strip.append(log_mahler(b_tilde, M_factor=M_factor))

    return FFTResult(
        N=N,
        log_m_chi=log_m_chi,
        log_m_chi_tilde=log_m_tilde,
        linfty_chi=linfty_chi,
        linfty_chi_tilde=linfty_tilde,
        argmax_q_chi=q_chi,
        argmax_a_chi=a_chi,
        argmax_q_tilde=q_tilde,
        argmax_a_tilde=a_tilde,
        M_grid=M,
        sqrt_pi_N=float(np.sqrt(pi_N)),
        pi_N=pi_N,
        alpha_2=alpha_2,
        f_at_minus_one_chi=complex(f_minus_one_chi),
        f_at_minus_one_tilde=complex(f_minus_one_tilde),
        log_m_bern_mean=float(np.mean(bern_log_m)),
        log_m_bern_std=float(np.std(bern_log_m, ddof=1) if len(bern_log_m) > 1 else 0.0),
        log_m_bern_strip_mean=float(np.mean(bern_log_m_strip)),
        log_m_bern_strip_std=float(np.std(bern_log_m_strip, ddof=1) if len(bern_log_m_strip) > 1 else 0.0),
        bern_replicates=n_bern,
    )


def run_bdj_block(N: int, n_bern: int, rng: np.random.Generator) -> BDJResult:
    chi = chi_p_indicator(N)
    pi_N = int(chi.sum())
    eps_chi = standardise(chi)
    chi_tilde, _ = parity_strip(chi)
    eps_tilde = standardise(chi_tilde)

    m4_chi, lmax_chi = toeplitz_m4(eps_chi)
    m4_tilde, lmax_tilde = toeplitz_m4(eps_tilde)

    bern_m4 = []
    bern_m4_strip = []
    for _ in range(n_bern):
        b = bernoulli_match_density(N, pi_N, rng)
        bern_m4.append(toeplitz_m4(standardise(b))[0])
        b_tilde, _ = parity_strip(b)
        bern_m4_strip.append(toeplitz_m4(standardise(b_tilde))[0])

    return BDJResult(
        N=N,
        m4_chi=m4_chi,
        m4_chi_tilde=m4_tilde,
        lambda_max_chi=lmax_chi,
        lambda_max_chi_tilde=lmax_tilde,
        m4_bern_mean=float(np.mean(bern_m4)),
        m4_bern_std=float(np.std(bern_m4, ddof=1) if len(bern_m4) > 1 else 0.0),
        m4_bern_strip_mean=float(np.mean(bern_m4_strip)),
        m4_bern_strip_std=float(np.std(bern_m4_strip, ddof=1) if len(bern_m4_strip) > 1 else 0.0),
        bern_replicates=n_bern,
    )


# -- Driver ---------------------------------------------------------------

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--N_log2", default="14,16,18", help="comma-list of log2 N for FFT block")
    p.add_argument("--bdj_N", default="500,1000,1500", help="comma-list of N for Toeplitz block")
    p.add_argument("--n_bern_fft", type=int, default=50)
    p.add_argument("--n_bern_bdj", type=int, default=20)
    p.add_argument("--M_factor", type=int, default=4)
    p.add_argument("--seed", type=int, default=20260430)
    p.add_argument("--out_json", default="parity_stripped_trinity_results.json")
    args = p.parse_args()

    rng = np.random.default_rng(args.seed)

    log2_list = [int(s.strip()) for s in args.N_log2.split(",") if s.strip()]
    bdj_list = [int(s.strip()) for s in args.bdj_N.split(",") if s.strip()]

    fft_results: list[FFTResult] = []
    print(f"\n=== FFT block (Mahler + L^infty) ===")
    for d in log2_list:
        N = 2 ** d
        t0 = time.time()
        res = run_fft_block(N, args.n_bern_fft, args.M_factor, rng)
        dt = time.time() - t0
        fft_results.append(res)
        print(
            f"  N=2^{d}={N}: pi(N)={res.pi_N}, alpha_2={res.alpha_2:.6f}, "
            f"|f(-1)|={abs(res.f_at_minus_one_chi):.2f} -> |f̃(-1)|={abs(res.f_at_minus_one_tilde):.2e}, "
            f"log m(f)={res.log_m_chi:.4f}, log m(f̃)={res.log_m_chi_tilde:.4f}, "
            f"BERN log m={res.log_m_bern_mean:.4f}±{res.log_m_bern_std:.4f}, "
            f"BERN-strip log m={res.log_m_bern_strip_mean:.4f}±{res.log_m_bern_strip_std:.4f} "
            f"[{dt:.1f}s]"
        )
        print(
            f"           L_inf(f)={res.linfty_chi:.2f} at q={res.argmax_q_chi} (sqrt pi(N)={res.sqrt_pi_N:.2f}); "
            f"L_inf(f̃)={res.linfty_chi_tilde:.2f} at q={res.argmax_q_tilde}"
        )

    bdj_results: list[BDJResult] = []
    print(f"\n=== BDJ block (Toeplitz m_4) ===")
    for N in bdj_list:
        t0 = time.time()
        res = run_bdj_block(N, args.n_bern_bdj, rng)
        dt = time.time() - t0
        bdj_results.append(res)
        print(
            f"  N={N}: m_4(chi_P)={res.m4_chi:.3f}, m_4(chĩ_P)={res.m4_chi_tilde:.3f}, "
            f"lambda_max scaled chi_P={res.lambda_max_chi:.3f} -> chĩ_P={res.lambda_max_chi_tilde:.3f}, "
            f"BERN m_4={res.m4_bern_mean:.3f}±{res.m4_bern_std:.3f}, "
            f"BERN-strip m_4={res.m4_bern_strip_mean:.3f}±{res.m4_bern_strip_std:.3f} [{dt:.1f}s]"
        )

    out = {
        "params": {
            "N_log2": log2_list,
            "bdj_N": bdj_list,
            "n_bern_fft": args.n_bern_fft,
            "n_bern_bdj": args.n_bern_bdj,
            "M_factor": args.M_factor,
            "seed": args.seed,
        },
        "fft_block": [
            {
                **{k: v for k, v in r.__dict__.items() if k not in ("f_at_minus_one_chi", "f_at_minus_one_tilde")},
                "f_at_minus_one_chi_re": r.f_at_minus_one_chi.real,
                "f_at_minus_one_chi_im": r.f_at_minus_one_chi.imag,
                "f_at_minus_one_tilde_re": r.f_at_minus_one_tilde.real,
                "f_at_minus_one_tilde_im": r.f_at_minus_one_tilde.imag,
            }
            for r in fft_results
        ],
        "bdj_block": [r.__dict__ for r in bdj_results],
    }
    out_path = Path(__file__).parent / args.out_json
    out_path.write_text(json.dumps(out, indent=2, default=str))
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
