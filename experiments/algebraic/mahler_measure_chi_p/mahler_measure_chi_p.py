"""
D10 — Mahler measure of the prime indicator polynomial f_N(z) = Σ_{n=1}^N χ_P(n) z^n.

Cross-domain technique: Mahler measure / Lehmer-Boyd / Weil height.
Refs: Smyth 2008 survey (CUP); Lehmer 1933 Ann. Math. 34;
Wikipedia: Mahler measure; Wikipedia: Lehmer's conjecture.

Goal: estimate log m(f_N) = ∫₀¹ log|f_N(e^{2π i θ})| dθ via Jensen-FFT,
compare to matched-density Bernoulli, Liouville, square-free baselines.

Falsifier (PRE-REGISTERED — F3 form):

  F1 (Lehmer-typical):  PRIMES log m / log N → 1/2 within sample noise of
                        density-matched Bernoulli baseline at all N.
                        => 38th pseudorandomness measure (B-grade).
  F2 (cyclotomic):      log m(f_N) = O((log N)^c) at N ≥ 2^14
                        => A-grade — first polylog representation of χ_P.
  F3 (intermediate):    log m(f_N) / log N converges to α ≠ 1/2 with
                        |z| > 3σ vs Bernoulli baseline AND vs deterministic
                        Liouville/squarefree controls.
                        => B-grade negative-shape / new structural fact.

Edge candidates: composes E2.13 (Gowers/HL), E2.14 (Anderson),
E2.15 (algebraic immunity), E2.16 (DPP failure), E2.17 (PH).

Usage:
  python mahler_measure_chi_p.py --quick           # tiny smoke test
  python mahler_measure_chi_p.py --main            # the headline run
  python mahler_measure_chi_p.py --cyclotomic 256  # factor f_N over Q[z]
"""
from __future__ import annotations

import argparse
import json
import os
import time
from dataclasses import dataclass, asdict
from pathlib import Path

import math
import numpy as np


# -----------------------------------------------------------------------------
# Indicator vectors
# -----------------------------------------------------------------------------

def sieve_primes(n: int) -> np.ndarray:
    """Boolean array of length n+1 with sieve[i] = (i is prime), i ≤ n."""
    if n < 2:
        return np.zeros(n + 1, dtype=bool)
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(math.isqrt(n)) + 1):
        if s[p]:
            s[p * p :: p] = False
    return s


def chi_prime_indicator(N: int) -> np.ndarray:
    """χ_P(n) for n ∈ [1, N]. Returns float64 array of length N."""
    s = sieve_primes(N)
    return s[1 : N + 1].astype(np.float64)


def liouville_minus(N: int) -> np.ndarray:
    """Indicator of {n : λ(n) = -1}, n ∈ [1, N]. λ(n) = (-1)^Ω(n)."""
    s = sieve_primes(N)
    omega = np.zeros(N + 1, dtype=np.int32)
    primes = np.where(s)[0]
    for p in primes:
        k = p
        while k <= N:
            cnt = 0
            kk = k
            while kk % p == 0:
                cnt += 1
                kk //= p
            omega[k] += cnt
            k += p
    lam = (-1) ** omega
    return (lam[1 : N + 1] == -1).astype(np.float64)


def squarefree_indicator(N: int) -> np.ndarray:
    """μ²(n) for n ∈ [1, N]: 1 if n is square-free, 0 otherwise."""
    sq = np.ones(N + 1, dtype=bool)
    sq[0] = False
    for p in range(2, int(math.isqrt(N)) + 1):
        sq[p * p :: p * p] = False
    return sq[1 : N + 1].astype(np.float64)


def random_bernoulli(N: int, density: float, rng: np.random.Generator) -> np.ndarray:
    """iid Bernoulli(density) of length N."""
    return (rng.random(N) < density).astype(np.float64)


def random_matched_count(N: int, k: int, rng: np.random.Generator) -> np.ndarray:
    """Random 0/1 vector of length N with EXACTLY k ones (sampled without
    replacement). Closer to a "matched-cardinality" baseline than iid."""
    v = np.zeros(N, dtype=np.float64)
    idx = rng.choice(N, size=k, replace=False)
    v[idx] = 1.0
    return v


# -----------------------------------------------------------------------------
# Mahler measure via Jensen-FFT
# -----------------------------------------------------------------------------

def log_mahler_fft(coeffs: np.ndarray, M: int = 1 << 18, eps_floor: float = 1e-30) -> dict:
    """Estimate log m(f) = ∫₀¹ log |f(e^{2π i θ})| dθ via M-point FFT.

    Pads `coeffs` (which represent f(z) = Σ c_n z^n with c_n the i-th entry,
    z = e^{2π i k / M}) to length M, runs FFT, then averages log|f|.

    Returns dict with log_m, max_log, min_log, deg, M, density, n_pad_zeros,
    zero_amplitude_count.
    """
    coeffs = np.asarray(coeffs, dtype=np.float64)
    deg = len(coeffs)
    if deg > M:
        raise ValueError("M must be >= len(coeffs); pick larger M")
    padded = np.zeros(M, dtype=np.float64)
    padded[:deg] = coeffs
    fz = np.fft.fft(padded)
    abs_fz = np.abs(fz)
    floored = np.maximum(abs_fz, eps_floor)
    log_abs = np.log(floored)
    log_m = float(np.mean(log_abs))
    return dict(
        log_m=log_m,
        max_log_abs=float(np.max(log_abs)),
        min_log_abs=float(np.min(log_abs)),
        max_abs=float(np.max(abs_fz)),
        min_abs=float(np.min(abs_fz)),
        deg=deg,
        M=M,
        density=float(np.mean(coeffs)),
        zero_amplitude_count=int(np.sum(abs_fz < eps_floor)),
    )


# -----------------------------------------------------------------------------
# Cyclotomic factorisation (small N) via sympy
# -----------------------------------------------------------------------------

def cyclotomic_factorisation(N: int, max_print: int = 20) -> dict:
    """Factor f_N(z) over Q[z] with sympy and report cyclotomic share."""
    import sympy
    from sympy import Poly, Symbol, ZZ
    from sympy.polys.specialpolys import cyclotomic_poly

    z = Symbol("z")
    chi = chi_prime_indicator(N)
    expr = sum(int(chi[n - 1]) * z ** n for n in range(1, N + 1))
    poly = Poly(expr, z, domain=ZZ)
    factors = poly.factor_list()
    leading = factors[0]
    fact_list = factors[1]

    # Collect cyclotomic indices
    cyclo_set: list[int] = []
    other_factors: list[tuple] = []
    for fac, mult in fact_list:
        d = fac.degree()
        if d == 0:
            continue
        # Test if fac is cyclotomic Φ_n for some n
        # Φ_n has degree φ(n); only need to check n ≤ small bound.
        matched = False
        for n in range(1, max(2 * d + 4, 200)):
            if sympy.totient(n) != d:
                continue
            if Poly(cyclotomic_poly(n, z), z, domain=ZZ) == fac:
                cyclo_set.append(int(n))
                matched = True
                break
        if not matched:
            # store first 10 coeffs of fac for inspection
            coeff_list = [int(c) for c in fac.all_coeffs()]
            other_factors.append((d, int(mult), coeff_list[: min(10, len(coeff_list))]))

    return dict(
        N=N,
        leading_coefficient=int(leading),
        cyclotomic_factor_indices=sorted(cyclo_set),
        n_cyclotomic_factors=len(cyclo_set),
        n_noncyclo_factors=len(other_factors),
        noncyclo_summary=other_factors[:max_print],
        total_degree=int(poly.degree()),
        total_terms=int(np.sum(chi != 0)),
    )


# -----------------------------------------------------------------------------
# Main experiments
# -----------------------------------------------------------------------------

@dataclass
class Run:
    label: str
    N: int
    M: int
    log_m: float
    log_m_per_logN: float
    density: float
    notes: str = ""


def run_main(out_dir: Path, N_list: list[int], M_fft: int, n_controls: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    rows: list[Run] = []
    raw: dict = {"N_list": N_list, "M_fft": M_fft, "n_controls": n_controls, "seed": seed}

    for N in N_list:
        logN = float(np.log(N))
        # primes
        chi = chi_prime_indicator(N)
        density = float(np.mean(chi))
        k_primes = int(np.sum(chi))
        t0 = time.time()
        meas = log_mahler_fft(chi, M=M_fft)
        meas["t_sec"] = time.time() - t0
        rows.append(Run("PRIMES", N, M_fft, meas["log_m"], meas["log_m"] / logN, density))
        raw[f"primes_N{N}"] = meas

        # Liouville-minus
        lam = liouville_minus(N)
        meas_lam = log_mahler_fft(lam, M=M_fft)
        rows.append(Run("LIOUVILLE", N, M_fft, meas_lam["log_m"], meas_lam["log_m"] / logN, float(np.mean(lam))))
        raw[f"liouville_N{N}"] = meas_lam

        # Square-free
        mu2 = squarefree_indicator(N)
        meas_mu = log_mahler_fft(mu2, M=M_fft)
        rows.append(Run("SQFREE", N, M_fft, meas_mu["log_m"], meas_mu["log_m"] / logN, float(np.mean(mu2))))
        raw[f"sqfree_N{N}"] = meas_mu

        # Bernoulli ensemble matched on density
        bern = []
        for _ in range(n_controls):
            v = random_bernoulli(N, density, rng)
            mb = log_mahler_fft(v, M=M_fft)
            bern.append(mb["log_m"])
        bern = np.array(bern)
        raw[f"bernoulli_N{N}"] = dict(mean=float(bern.mean()), std=float(bern.std(ddof=1)),
                                       median=float(np.median(bern)),
                                       q05=float(np.quantile(bern, 0.05)),
                                       q95=float(np.quantile(bern, 0.95)),
                                       min=float(bern.min()), max=float(bern.max()),
                                       n=n_controls)

        # Matched-cardinality random ensemble (exactly k_primes ones)
        match = []
        for _ in range(n_controls):
            v = random_matched_count(N, k_primes, rng)
            mb = log_mahler_fft(v, M=M_fft)
            match.append(mb["log_m"])
        match = np.array(match)
        raw[f"matched_card_N{N}"] = dict(mean=float(match.mean()), std=float(match.std(ddof=1)),
                                          median=float(np.median(match)),
                                          q05=float(np.quantile(match, 0.05)),
                                          q95=float(np.quantile(match, 0.95)),
                                          min=float(match.min()), max=float(match.max()),
                                          n=n_controls)

        # Z-scores: (PRIMES - mean(baseline)) / std(baseline)
        z_b = (meas["log_m"] - bern.mean()) / bern.std(ddof=1) if bern.std(ddof=1) > 0 else 0.0
        z_m = (meas["log_m"] - match.mean()) / match.std(ddof=1) if match.std(ddof=1) > 0 else 0.0
        z_l = (meas["log_m"] - meas_lam["log_m"])  # Liouville is deterministic → diff
        raw[f"z_scores_N{N}"] = dict(
            primes_log_m=meas["log_m"],
            bernoulli_mean=float(bern.mean()),
            bernoulli_std=float(bern.std(ddof=1)),
            matched_mean=float(match.mean()),
            matched_std=float(match.std(ddof=1)),
            z_primes_vs_bernoulli=z_b,
            z_primes_vs_matched=z_m,
            primes_minus_liouville_log_m=z_l,
            primes_minus_sqfree_log_m=meas["log_m"] - meas_mu["log_m"],
        )

        print(f"  N={N:7d}: PRIMES log m = {meas['log_m']:+.5f}, "
              f"BERN = {bern.mean():+.5f}±{bern.std(ddof=1):.4f}, "
              f"MATCH = {match.mean():+.5f}±{match.std(ddof=1):.4f}, "
              f"LIOUV = {meas_lam['log_m']:+.5f}, "
              f"SQFREE = {meas_mu['log_m']:+.5f}, "
              f"z(B)={z_b:+.2f}, z(M)={z_m:+.2f}")

    raw["rows"] = [asdict(r) for r in rows]
    return raw


def fit_alpha(N_list: list[int], log_m_values: list[float]) -> dict:
    """Linear regression of log(m) vs log(log(N)) and log(N) — extract α
    such that log m ≈ α · log N + β."""
    x = np.log(N_list)
    y = np.array(log_m_values)
    # OLS slope on log N
    n = len(x)
    sx = x.sum()
    sy = y.sum()
    sxx = (x * x).sum()
    sxy = (x * y).sum()
    slope = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    intercept = (sy - slope * sx) / n
    resid = y - (slope * x + intercept)
    rss = float(np.sum(resid * resid))
    tss = float(np.sum((y - y.mean()) ** 2))
    return dict(slope=float(slope), intercept=float(intercept),
                rss=rss, tss=tss, r2=1 - rss / tss if tss > 0 else float("nan"),
                fit="log m ≈ α · log N + β with α = slope")


# -----------------------------------------------------------------------------
# Cyclotomic / verification utility
# -----------------------------------------------------------------------------

def jensen_via_roots_smallN(N: int) -> dict:
    """For small N, also compute log m exactly from roots (mpmath polyroots)."""
    import mpmath
    mpmath.mp.dps = 40
    chi = chi_prime_indicator(N)
    coeffs = chi.tolist()  # constant ... z^N (coeffs[0] = c_0, coeffs[N] = c_N)
    # Build polynomial f(z) = Σ chi[n-1] z^n for n=1..N → leading = χ_P(N)
    poly_high_to_low = [chi[N - 1 - i] if (N - i) >= 1 and (N - i) <= N else 0.0
                        for i in range(N + 1)]
    # Above expression: index i in poly_high_to_low corresponds to z^{N-i}
    # but f has z^n for n=1..N, so coefficient on z^{N-i} = chi[N-i-1] when i<N.
    # Correctly:
    poly_high_to_low = [0.0] * (N + 1)
    for n in range(1, N + 1):
        poly_high_to_low[N - n] = float(chi[n - 1])
    # mpmath wants leading coefficient first
    # Drop leading zeros (deg < N when χ_P(N) = 0)
    while len(poly_high_to_low) > 1 and poly_high_to_low[0] == 0.0:
        poly_high_to_low.pop(0)
    if len(poly_high_to_low) <= 1:
        return dict(N=N, log_m_roots=float("nan"), n_roots=0,
                    roots_outside_unit_circle=0)
    leading = poly_high_to_low[0]
    roots = mpmath.polyroots(poly_high_to_low, maxsteps=200, extraprec=200)
    abs_roots = [float(abs(r)) for r in roots]
    out = sum(np.log(max(1.0, ar)) for ar in abs_roots)
    log_m = float(np.log(abs(leading)) + out)
    return dict(N=N, leading=float(leading), log_m_roots=log_m,
                n_roots=len(roots),
                roots_outside_unit_circle=int(sum(1 for ar in abs_roots if ar > 1.0)),
                max_root_modulus=float(max(abs_roots)),
                min_root_modulus=float(min(abs_roots)))


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(Path(__file__).with_name("results")), help="output JSON dir")
    ap.add_argument("--quick", action="store_true", help="tiny smoke run")
    ap.add_argument("--main", action="store_true", help="headline run")
    ap.add_argument("--cyclotomic", type=int, default=0, help="factor f_N over Q[z] for given N")
    ap.add_argument("--roots", type=int, default=0, help="cross-check log m via roots for small N")
    ap.add_argument("--scale", action="store_true", help="extra-large N scaling run")
    ap.add_argument("--N", type=int, nargs="*", default=None)
    ap.add_argument("--M", type=int, default=1 << 18)
    ap.add_argument("--controls", type=int, default=100)
    ap.add_argument("--seed", type=int, default=20260427)
    args = ap.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.quick:
        raw = run_main(out_dir, [256, 1024, 4096], 1 << 14, 30, args.seed)
        fit = fit_alpha([256, 1024, 4096], [raw[f"primes_N{N}"]["log_m"] for N in [256, 1024, 4096]])
        raw["fit_primes"] = fit
        raw["fit_bernoulli_mean"] = fit_alpha([256, 1024, 4096], [raw[f"bernoulli_N{N}"]["mean"] for N in [256, 1024, 4096]])
        with open(out_dir / "quick.json", "w") as fh:
            json.dump(raw, fh, indent=2)
        print("\nFit α (PRIMES):", raw["fit_primes"])
        print("Fit α (BERN mean):", raw["fit_bernoulli_mean"])

    if args.main:
        N_list = args.N or [1024, 4096, 16384, 65536]
        raw = run_main(out_dir, N_list, args.M, args.controls, args.seed)
        for tag in ["primes", "liouville", "sqfree"]:
            raw[f"fit_{tag}"] = fit_alpha(
                N_list, [raw[f"{tag}_N{N}"]["log_m"] for N in N_list]
            )
        raw["fit_bernoulli_mean"] = fit_alpha(
            N_list, [raw[f"bernoulli_N{N}"]["mean"] for N in N_list]
        )
        raw["fit_matched_mean"] = fit_alpha(
            N_list, [raw[f"matched_card_N{N}"]["mean"] for N in N_list]
        )
        with open(out_dir / "main.json", "w") as fh:
            json.dump(raw, fh, indent=2)
        print("\nFit α — PRIMES :", raw["fit_primes"])
        print("Fit α — LIOUV  :", raw["fit_liouville"])
        print("Fit α — SQFREE :", raw["fit_sqfree"])
        print("Fit α — BERN   :", raw["fit_bernoulli_mean"])
        print("Fit α — MATCH  :", raw["fit_matched_mean"])

    if args.scale:
        # Extra-large scaling without controls (for trend confirmation)
        N_list = args.N or [2 ** k for k in range(10, 19)]
        rows = []
        for N in N_list:
            chi = chi_prime_indicator(N)
            M = max(args.M, 4 * N)
            meas = log_mahler_fft(chi, M=M)
            rows.append(dict(N=N, M=M, log_m=meas["log_m"], density=meas["density"]))
            print(f"  N={N:8d}  log m = {meas['log_m']:+.6f}, density = {meas['density']:.5f}")
        fit = fit_alpha([r["N"] for r in rows], [r["log_m"] for r in rows])
        with open(out_dir / "scale.json", "w") as fh:
            json.dump(dict(rows=rows, fit=fit), fh, indent=2)
        print("Fit α (PRIMES extended):", fit)

    if args.cyclotomic > 0:
        fac = cyclotomic_factorisation(args.cyclotomic)
        with open(out_dir / f"cyclo_N{args.cyclotomic}.json", "w") as fh:
            json.dump(fac, fh, indent=2)
        print(f"\nFactorisation of f_{args.cyclotomic}(z):")
        print(f"  cyclotomic factors at indices: {fac['cyclotomic_factor_indices']}")
        print(f"  # non-cyclotomic factors: {fac['n_noncyclo_factors']}")
        for d, mult, head in fac["noncyclo_summary"]:
            print(f"    deg={d:4d} (mult={mult}) — leading coeffs: {head[:6]}...")

    if args.roots > 0:
        ro = jensen_via_roots_smallN(args.roots)
        with open(out_dir / f"roots_N{args.roots}.json", "w") as fh:
            json.dump(ro, fh, indent=2)
        print(f"\nlog m(f_{args.roots}) via roots: {ro['log_m_roots']:.6f}")
        print(f"  outside unit circle: {ro['roots_outside_unit_circle']} of {ro['n_roots']}")
        print(f"  max |root|: {ro['max_root_modulus']:.4f}")


if __name__ == "__main__":
    main()
