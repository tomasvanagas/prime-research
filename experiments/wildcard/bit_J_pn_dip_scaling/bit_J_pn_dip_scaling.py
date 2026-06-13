#!/usr/bin/env python3
"""
F1.a.i — Dip-depth scaling law of the bit-J RH-shadow valley.

Refines E1.3 (per-bit difficulty of p(n)) via the cross-modulus
universality discovered in S199.  The S199 closure flagged a closed-form
prediction for the relative dip depth `rel(m) := ag_Li(J*(m)) · m`,
where `J*(m) := ⌊log_m(p_N) / 2⌋` is the m-adic half-conductor scale.
This script tests that closed form across the dense modulus grid
m ∈ {2, 3, ..., 30}.

Pre-stated falsifiers in bit_J_pn_dip_scaling_results.md.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
SQRT2 = math.sqrt(2.0)
SQRT_TWO_PI = math.sqrt(2.0 * math.pi)


def norm_cdf(x: float, mu: float, sigma: float) -> float:
    return 0.5 * (1.0 + math.erf((x - mu) / (sigma * SQRT2)))


def euler_phi(m: int) -> int:
    n = m
    phi = m
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            phi -= phi // p
        p += 1
    if n > 1:
        phi -= phi // n
    return phi


def radical(m: int) -> int:
    """Squarefree radical of m."""
    n = m
    r = 1
    p = 2
    while p * p <= n:
        if n % p == 0:
            r *= p
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        r *= n
    return r


def sieve_primes(limit: int) -> np.ndarray:
    if limit < 2:
        return np.array([], dtype=np.uint64)
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = False
    return np.nonzero(sieve)[0].astype(np.uint64)


def li_inverse_vectorised(n_arr: np.ndarray) -> np.ndarray:
    """Newton iteration on the asymptotic Li series. Reused from S146/S199."""
    n = n_arr.astype(np.float64)
    safe = np.maximum(n, 6.0)
    L0 = np.log(safe)
    LL = np.log(np.maximum(L0, 1.5))
    x = safe * (L0 + LL - 1.0)
    x = np.maximum(x, 2.0)
    for _ in range(40):
        x = np.maximum(x, 2.0)
        Lx = np.log(x)
        s = np.ones_like(x)
        term = np.ones_like(x)
        for k in range(1, 10):
            term *= k / Lx
            s += term
        Li = x / Lx * s
        Li_prime = 1.0 / Lx
        delta = (n - Li) / Li_prime
        x = x + delta
        x = np.maximum(x, 2.0)
    small = np.where(n_arr < 6)[0]
    base = np.array([2.0, 3.0, 5.0, 7.0, 11.0])
    for i in small:
        x[i] = base[int(n_arr[i]) - 1]
    return x


def empirical_agreement(primes: np.ndarray, li_round: np.ndarray, m: int, J: int) -> float:
    if J == 0:
        d_p = primes % m
        d_l = li_round % m
    else:
        m_J = m ** J
        d_p = (primes // m_J) % m
        d_l = (li_round // m_J) % m
    return float((d_p == d_l).mean())


def gaussian_y_prediction(mu_e: float, sigma_e: float, m: int, J: int) -> float:
    """Closed form: Y := (r + e)/m^J approximated as Gaussian with
    mu_Y = mu_e/m^J + 1/2, sigma_Y² = sigma_e²/m^(2J) + 1/12.

    P[floor(Y) ≡ 0 mod m] = sum_k [Φ(km+1; mu_Y, sigma_Y) - Φ(km; mu_Y, sigma_Y)].
    """
    m_J = m ** J
    mu_Y = mu_e / m_J + 0.5
    var_Y = (sigma_e / m_J) ** 2 + 1.0 / 12.0
    sigma_Y = math.sqrt(var_Y)
    span = 8.0 * sigma_Y + 1.0
    k_lo = int(math.floor((mu_Y - span) / m)) - 1
    k_hi = int(math.ceil((mu_Y + span) / m)) + 1
    p = 0.0
    for k in range(k_lo, k_hi + 1):
        p += norm_cdf(k * m + 1.0, mu_Y, sigma_Y) - norm_cdf(k * m, mu_Y, sigma_Y)
    return p


def empirical_uniform_r_prediction(e_arr: np.ndarray, m: int, J: int) -> float:
    """Empirical-e + uniform-r: per-n analytic average over r ~ Unif(0, m^J).

    For each n, conditional Pr[floor((r + e_n)/m^J) ≡ 0 mod m] is
    [q_n ≡ 0 mod m] · (1 - frac_n) + [q_n ≡ m-1 mod m] · frac_n,
    where q_n = e_n // m^J, frac_n = (e_n mod m^J) / m^J.
    """
    m_J = m ** J
    q = (e_arr // m_J).astype(np.int64)
    r_e = (e_arr % m_J).astype(np.float64)
    frac = r_e / m_J
    q_mod = (q % m).astype(np.int64)
    contrib = np.zeros_like(frac)
    mask0 = q_mod == 0
    mask_m1 = q_mod == (m - 1)
    contrib[mask0] += 1.0 - frac[mask0]
    contrib[mask_m1] += frac[mask_m1]
    return float(contrib.mean())


def gaussian_uniform_r_prediction(mu_e: float, sigma_e: float, m: int, J: int) -> float:
    """Gaussian-e + uniform-r: integrate per-e formula against Gaussian.

    P[delta mod m == 0]
      = Σ_{k ≡ 0   mod m} ∫_0^1 (1-u) φ(km^J + u·m^J; mu_e, sigma_e) m^J du
      + Σ_{k ≡ m-1 mod m} ∫_0^1   u  φ(km^J + u·m^J; mu_e, sigma_e) m^J du.

    Uses Simpson's rule with 21 nodes per unit interval.
    """
    m_J = m ** J
    span = 7.0 * sigma_e + m_J
    k_lo = int(math.floor((mu_e - span) / m_J)) - 1
    k_hi = int(math.ceil((mu_e + span) / m_J)) + 1
    n_quad = 21
    us = np.linspace(0.0, 1.0, n_quad)
    simpson_w = np.ones(n_quad)
    simpson_w[1:-1:2] = 4.0
    simpson_w[2:-1:2] = 2.0
    simpson_factor = (1.0 / (n_quad - 1)) / 3.0
    one_minus_u = 1.0 - us
    p = 0.0
    inv_norm = 1.0 / (sigma_e * SQRT_TWO_PI)
    for k in range(k_lo, k_hi + 1):
        kmod = k % m
        if kmod == 0 or kmod == (m - 1):
            xs = (k * m_J) + us * m_J
            phi_xs = inv_norm * np.exp(-((xs - mu_e) ** 2) / (2.0 * sigma_e ** 2))
            if kmod == 0:
                integrand = one_minus_u * phi_xs * m_J
            else:
                integrand = us * phi_xs * m_J
            integral = simpson_factor * float(np.sum(simpson_w * integrand))
            p += integral
    return p


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=200_000_000)
    parser.add_argument(
        "--moduli",
        type=str,
        default=",".join(str(m) for m in range(2, 31)),
    )
    parser.add_argument(
        "--out",
        type=str,
        default=str(HERE / "bit_J_pn_dip_scaling_results.json"),
    )
    args = parser.parse_args()

    moduli = [int(x) for x in args.moduli.split(",")]
    print(f"Sieving primes up to {args.limit:,}...", flush=True)
    t0 = time.time()
    primes = sieve_primes(args.limit)
    N = len(primes)
    print(
        f"  Sieved {N:,} primes in {time.time() - t0:.1f}s. "
        f"Last prime = {int(primes[-1]):,}",
        flush=True,
    )

    n_arr = np.arange(1, N + 1, dtype=np.uint64)
    print("  Computing Li^{-1}(n)...", flush=True)
    t0 = time.time()
    li_pred = np.round(li_inverse_vectorised(n_arr)).astype(np.uint64)
    print(f"  Li^{{-1}} done in {time.time() - t0:.1f}s", flush=True)

    e_arr = primes.astype(np.int64) - li_pred.astype(np.int64)
    mu_e = float(e_arr.mean())
    sigma_e = float(e_arr.std())
    median_e = float(np.median(e_arr))
    skew_e = (
        float(np.mean((e_arr - mu_e) ** 3) / sigma_e ** 3) if sigma_e > 0 else 0.0
    )
    pos_frac = float((e_arr > 0).mean())
    e_min = int(e_arr.min())
    e_max = int(e_arr.max())
    print("\nEmpirical e_n stats:", flush=True)
    print(
        f"  N = {N:,}, mean = {mu_e:.2f}, std = {sigma_e:.2f}, "
        f"median = {median_e:.2f}, min = {e_min}, max = {e_max}",
        flush=True,
    )
    print(f"  skew = {skew_e:.3f}, P(e > 0) = {pos_frac:.4f}", flush=True)

    last_prime = int(primes[-1])
    results = {
        "limit": int(args.limit),
        "last_prime": last_prime,
        "N": N,
        "moduli": moduli,
        "e_stats": {
            "mean": mu_e,
            "std": sigma_e,
            "median": median_e,
            "skew": skew_e,
            "pos_frac": pos_frac,
            "min": e_min,
            "max": e_max,
        },
        "rows": [],
    }

    print(
        "\n=== Per-modulus dip-scaling table (J* = floor(log_m(p_N)/2)) ===",
        flush=True,
    )
    header = (
        f"{'m':>3} {'φ(m)':>4} {'rad':>3} {'J*':>3} "
        f"{'mu_Y':>7} {'sig_Y':>7} "
        f"{'ag_emp':>10} {'rel_emp':>7} "
        f"{'rel_pY':>7} {'rel_pR':>7} {'rel_pGR':>7} "
        f"{'errY%':>6} {'errR%':>6} {'errGR%':>6}"
    )
    print(header, flush=True)

    for m in moduli:
        log_m_p = math.log(last_prime) / math.log(m)
        J_star = int(log_m_p / 2)
        m_J = m ** J_star
        mu_Y = mu_e / m_J + 0.5
        sig_Y = math.sqrt((sigma_e / m_J) ** 2 + 1.0 / 12.0)
        phi_m = euler_phi(m)
        rad_m = radical(m)
        ag_emp = empirical_agreement(primes, li_pred, m, J_star)
        baseline = 1.0 / m
        rel_emp = ag_emp / baseline
        ag_pY = gaussian_y_prediction(mu_e, sigma_e, m, J_star)
        rel_pY = ag_pY / baseline
        ag_pR = empirical_uniform_r_prediction(e_arr, m, J_star)
        rel_pR = ag_pR / baseline
        ag_pGR = gaussian_uniform_r_prediction(mu_e, sigma_e, m, J_star)
        rel_pGR = ag_pGR / baseline
        # Relative-error percentages (signed); reference = empirical.
        eps = 1e-12
        err_pY = (rel_pY - rel_emp) / max(abs(rel_emp), eps) * 100.0
        err_pR = (rel_pR - rel_emp) / max(abs(rel_emp), eps) * 100.0
        err_pGR = (rel_pGR - rel_emp) / max(abs(rel_emp), eps) * 100.0
        row = {
            "m": m,
            "phi_m": phi_m,
            "radical_m": rad_m,
            "J_star": J_star,
            "log_m_p": log_m_p,
            "m_J_star": m_J,
            "mu_Y": mu_Y,
            "sigma_Y": sig_Y,
            "ag_emp": ag_emp,
            "rel_emp": rel_emp,
            "ag_pred_gaussian_Y": ag_pY,
            "rel_pred_gaussian_Y": rel_pY,
            "ag_pred_empirical_R": ag_pR,
            "rel_pred_empirical_R": rel_pR,
            "ag_pred_gaussian_R": ag_pGR,
            "rel_pred_gaussian_R": rel_pGR,
            "rel_err_pct_Y": err_pY,
            "rel_err_pct_R": err_pR,
            "rel_err_pct_GR": err_pGR,
        }
        results["rows"].append(row)
        print(
            f"{m:>3} {phi_m:>4} {rad_m:>3} {J_star:>3} "
            f"{mu_Y:>7.3f} {sig_Y:>7.3f} "
            f"{ag_emp:>10.6f} {rel_emp:>7.4f} "
            f"{rel_pY:>7.4f} {rel_pR:>7.4f} {rel_pGR:>7.4f} "
            f"{err_pY:>+6.1f} {err_pR:>+6.1f} {err_pGR:>+6.1f}",
            flush=True,
        )

    # Aggregate F4 statistic — mean abs err of rel_pred_R vs rel_pred_Y.
    err_Y = np.mean([abs(r["rel_pred_gaussian_Y"] - r["rel_emp"]) for r in results["rows"]])
    err_R = np.mean([abs(r["rel_pred_empirical_R"] - r["rel_emp"]) for r in results["rows"]])
    err_GR = np.mean([abs(r["rel_pred_gaussian_R"] - r["rel_emp"]) for r in results["rows"]])
    results["mean_abs_err"] = {
        "rel_pred_gaussian_Y": float(err_Y),
        "rel_pred_empirical_R": float(err_R),
        "rel_pred_gaussian_R": float(err_GR),
    }
    print(
        f"\n=== Aggregate mean |Δrel| ===\n"
        f"  Gaussian-Y    : {err_Y:.4f}\n"
        f"  Empirical-R   : {err_R:.4f}\n"
        f"  Gaussian-R    : {err_GR:.4f}",
        flush=True,
    )

    with open(args.out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
