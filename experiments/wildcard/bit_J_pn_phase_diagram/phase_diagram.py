#!/usr/bin/env python3
"""
F1.a.i.γ — Peak-vs-dip phase diagram of the bit-J RH-shadow agreement.

Refines E1.3 by re-parameterising rel_emp(m) = ag_Li(m, J*) · m in the
two-parameter (α, σ_norm) phase plane:
  α(m)        = (μ_Y(m) mod m) / m  ∈ [0, 1)
  σ_norm(m)   = σ_Y(m) / m         ∈ (0, ∞)

with μ_Y, σ_Y from the Gaussian-Y formula (S218).

Three closed-form candidates compared against rel_emp(m):
  (Y)   Gaussian-Y in original parameterisation (S218 baseline).
  (PE)  Phase-Exact: same Gaussian-Y formula in (α, σ_norm, 1/m) form.
  (PL)  Phase-Lim:   wrapped Gaussian density (1/m → 0 limit).
  (R)   Empirical-r: empirical e + uniform r (S218 essentially-exact baseline).

Pre-stated falsifiers in phase_diagram_results.md.
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


def norm_cdf(x: float, mu: float = 0.0, sigma: float = 1.0) -> float:
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
    for i in range(2, int(limit ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = False
    return np.nonzero(sieve)[0].astype(np.uint64)


def li_inverse_vectorised(n_arr: np.ndarray) -> np.ndarray:
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
    """rel = m · Σ_k [Φ(km+1; mu_Y, sigma_Y) − Φ(km; mu_Y, sigma_Y)]
    with Y := (r + e)/m^J approximated as Gaussian with
    mu_Y = mu_e/m^J + 1/2, sigma_Y² = sigma_e²/m^(2J) + 1/12.
    Matches S218's gaussian_y_prediction.
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
    """Empirical-e + uniform-r — essentially exact in S218."""
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


def phase_exact_prediction(alpha: float, sigma_norm: float, m: int) -> float:
    """P[⌊Y⌋ ≡ 0 mod m] in (α, σ_norm, 1/m) form — exact Gaussian-Y re-param.

    P = Σ_j [Φ_std((j − α)/σ_n + 1/(m σ_n)) − Φ_std((j − α)/σ_n)].

    Returns the probability (multiply by m to get rel). Equivalent to
    gaussian_y_prediction at the same (mu_Y, sigma_Y).
    """
    inv_msig = 1.0 / (m * sigma_norm)
    span = 8.0 * sigma_norm + 1.0
    j_lo = int(math.floor(-span)) - 1
    j_hi = int(math.ceil(span)) + 1
    p = 0.0
    for j in range(j_lo, j_hi + 1):
        z = (j - alpha) / sigma_norm
        p += norm_cdf(z + inv_msig) - norm_cdf(z)
    return p


def phase_lim_prediction(alpha: float, sigma_norm: float) -> float:
    """rel_phase_lim(α, σ_n) = Σ_j (1/(σ_n √2π)) exp(−(j − α)²/(2σ_n²))

    This is the wrapped-Gaussian density at offset α with σ_n.
    """
    span = 8.0 * sigma_norm + 1.0
    j_lo = int(math.floor(-span)) - 1
    j_hi = int(math.ceil(span)) + 1
    inv_norm = 1.0 / (sigma_norm * SQRT_TWO_PI)
    p = 0.0
    for j in range(j_lo, j_hi + 1):
        z = (j - alpha) / sigma_norm
        p += math.exp(-0.5 * z * z)
    return inv_norm * p


def phase_lim_poisson(alpha: float, sigma_norm: float, n_terms: int = 50) -> float:
    """Poisson-summation form: rel_phase_lim(α, σ_n) =
       1 + 2 Σ_{n=1}^∞ exp(−2π² n² σ_n²) cos(2π n α).
    Used as a sanity check on the wrapped-Gaussian sum above.
    """
    s = 1.0
    for n in range(1, n_terms + 1):
        s += 2.0 * math.exp(-2.0 * math.pi ** 2 * n * n * sigma_norm * sigma_norm) * math.cos(2.0 * math.pi * n * alpha)
    return s


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=200_000_000)
    parser.add_argument(
        "--moduli",
        type=str,
        default=",".join(
            str(m) for m in (
                list(range(2, 101))
                + [110, 120, 140, 170, 200, 250, 300, 400, 500, 700, 1000, 1500]
            )
        ),
    )
    parser.add_argument(
        "--out",
        type=str,
        default=str(HERE / "phase_diagram_results.json"),
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
    print(
        f"\nEmpirical e_n: mean = {mu_e:.2f}, std = {sigma_e:.2f}, "
        f"min = {int(e_arr.min())}, max = {int(e_arr.max())}",
        flush=True,
    )

    last_prime = int(primes[-1])
    results = {
        "limit": int(args.limit),
        "last_prime": last_prime,
        "N": N,
        "moduli": moduli,
        "e_stats": {"mean": mu_e, "std": sigma_e, "min": int(e_arr.min()), "max": int(e_arr.max())},
        "rows": [],
    }

    print("\n=== Per-modulus phase-diagram table ===", flush=True)
    print(
        f"{'m':>3} {'J*':>3} {'mu_Y':>8} {'sig_Y':>7} "
        f"{'alpha':>6} {'sigN':>7} "
        f"{'rel_emp':>8} {'rel_pY':>8} {'rel_pPE':>8} {'rel_pPL':>8} {'rel_pR':>8} "
        f"{'errPL%':>7}",
        flush=True,
    )

    for m in moduli:
        log_m_p = math.log(last_prime) / math.log(m)
        J_star = int(log_m_p / 2)
        m_J = m ** J_star
        mu_Y = mu_e / m_J + 0.5
        var_Y = (sigma_e / m_J) ** 2 + 1.0 / 12.0
        sig_Y = math.sqrt(var_Y)
        # Phase coordinates.
        alpha = (mu_Y % m) / m  # in [0, 1)
        sigma_norm = sig_Y / m
        ag_emp = empirical_agreement(primes, li_pred, m, J_star)
        baseline = 1.0 / m
        rel_emp = ag_emp / baseline
        rel_pY = gaussian_y_prediction(mu_e, sigma_e, m, J_star) / baseline
        rel_pPE = phase_exact_prediction(alpha, sigma_norm, m) / baseline  # P / (1/m) = m·P
        rel_pPL = phase_lim_prediction(alpha, sigma_norm)  # already in rel units
        rel_pR = empirical_uniform_r_prediction(e_arr, m, J_star) / baseline
        eps = 1e-12
        err_pPL = (rel_pPL - rel_emp) / max(abs(rel_emp), eps) * 100.0
        row = {
            "m": m,
            "phi_m": euler_phi(m),
            "radical_m": radical(m),
            "J_star": J_star,
            "mu_Y": mu_Y,
            "sigma_Y": sig_Y,
            "alpha": alpha,
            "sigma_norm": sigma_norm,
            "ag_emp": ag_emp,
            "rel_emp": rel_emp,
            "rel_pred_gaussian_Y": rel_pY,
            "rel_pred_phase_exact": rel_pPE,
            "rel_pred_phase_lim": rel_pPL,
            "rel_pred_empirical_R": rel_pR,
            "rel_err_pct_phase_lim": err_pPL,
        }
        results["rows"].append(row)
        print(
            f"{m:>3} {J_star:>3} {mu_Y:>8.3f} {sig_Y:>7.3f} "
            f"{alpha:>6.3f} {sigma_norm:>7.3f} "
            f"{rel_emp:>8.4f} {rel_pY:>8.4f} {rel_pPE:>8.4f} {rel_pPL:>8.4f} {rel_pR:>8.4f} "
            f"{err_pPL:>+7.1f}",
            flush=True,
        )

    # --- F-verdict aggregation ---
    rows = results["rows"]

    err_Y = np.mean([abs(r["rel_pred_gaussian_Y"] - r["rel_emp"]) for r in rows])
    err_PE = np.mean([abs(r["rel_pred_phase_exact"] - r["rel_emp"]) for r in rows])
    err_PL = np.mean([abs(r["rel_pred_phase_lim"] - r["rel_emp"]) for r in rows])
    err_R = np.mean([abs(r["rel_pred_empirical_R"] - r["rel_emp"]) for r in rows])
    results["mean_abs_err"] = {
        "rel_pred_gaussian_Y": float(err_Y),
        "rel_pred_phase_exact": float(err_PE),
        "rel_pred_phase_lim": float(err_PL),
        "rel_pred_empirical_R": float(err_R),
    }
    print(
        f"\n=== Aggregate mean |Δrel| (vs rel_emp) over m ∈ {{2..100}} ===\n"
        f"  Gaussian-Y      : {err_Y:.4f}\n"
        f"  Phase-Exact     : {err_PE:.4f}\n"
        f"  Phase-Lim (1/m→0): {err_PL:.4f}\n"
        f"  Empirical-r     : {err_R:.4f}",
        flush=True,
    )

    # F_α — wrapped-Gaussian-lim accuracy at small σ_norm.
    small_sig = [r for r in rows if r["sigma_norm"] < 0.5]
    if small_sig:
        f_alpha_passes = sum(
            1 for r in small_sig
            if abs(r["rel_pred_phase_lim"] - r["rel_emp"]) / max(r["rel_emp"], 1.0) <= 0.20
        )
        f_alpha_pass_rate = f_alpha_passes / len(small_sig)
    else:
        f_alpha_pass_rate = float("nan")
    results["F_alpha_pass_rate"] = f_alpha_pass_rate
    print(f"\nF_α (rel_phase_lim ≤20% err for σ_norm < 0.5): {f_alpha_pass_rate*100:.1f}% of {len(small_sig)} cells")

    # F_β — U-shape on σ_norm < 1.0 sub-population.
    sig_lt1 = [r for r in rows if r["sigma_norm"] < 1.0]
    decile_means = {}
    decile_counts = {}
    for r in sig_lt1:
        d = min(int(r["alpha"] * 10), 9)
        decile_means.setdefault(d, []).append(r["rel_emp"])
        decile_counts.setdefault(d, 0)
        decile_counts[d] += 1
    decile_avg = {d: float(np.mean(v)) for d, v in decile_means.items()}
    if decile_avg:
        max_d = max(decile_avg, key=decile_avg.get)
        min_d = min(decile_avg, key=decile_avg.get)
        f_beta_holds = (max_d in (0, 9)) and (min_d in (4, 5))
    else:
        max_d = min_d = -1
        f_beta_holds = False
    results["F_beta_decile_means"] = decile_avg
    results["F_beta_decile_counts"] = decile_counts
    results["F_beta_max_decile"] = max_d
    results["F_beta_min_decile"] = min_d
    results["F_beta_holds"] = bool(f_beta_holds)
    print(f"F_β (decile max ∈ {{0,9}} AND min ∈ {{4,5}}): max-decile={max_d}, min-decile={min_d} → {'HOLDS' if f_beta_holds else 'FAILS'}")
    print(f"     α-decile means: {decile_avg}")

    # F_γ — flat for σ_norm ≥ 2.0.
    flat_set = [r for r in rows if r["sigma_norm"] >= 2.0]
    if flat_set:
        flat_passes = sum(1 for r in flat_set if abs(r["rel_emp"] - 1.0) <= 0.5)
        f_gamma_pass_rate = flat_passes / len(flat_set)
    else:
        f_gamma_pass_rate = float("nan")
    results["F_gamma_pass_rate"] = f_gamma_pass_rate
    print(f"F_γ (|rel_emp − 1.0| ≤ 0.5 for σ_norm ≥ 2.0): {f_gamma_pass_rate*100:.1f}% of {len(flat_set)} cells")

    # F_δ — symmetry α ↔ 1−α.
    sig_lt1_sorted = sorted(sig_lt1, key=lambda r: (r["alpha"], r["sigma_norm"]))
    paired_passes = 0
    paired_total = 0
    for r in sig_lt1:
        # Find best mirror partner.
        target_alpha = 1.0 - r["alpha"]
        target_sig = r["sigma_norm"]
        best = None
        best_d = float("inf")
        for r2 in sig_lt1:
            if r2["m"] == r["m"]:
                continue
            d = abs(r2["alpha"] - target_alpha) + abs(r2["sigma_norm"] - target_sig)
            if d < best_d:
                best_d = d
                best = r2
        if best is None or best_d > 0.15:
            continue
        paired_total += 1
        denom = max(r["rel_emp"], best["rel_emp"], eps)
        if abs(r["rel_emp"] - best["rel_emp"]) / denom < 0.30:
            paired_passes += 1
    f_delta_rate = (paired_passes / paired_total) if paired_total > 0 else float("nan")
    results["F_delta_paired_passes"] = paired_passes
    results["F_delta_paired_total"] = paired_total
    results["F_delta_pass_rate"] = f_delta_rate
    print(f"F_δ (mirror-pair |Δrel|/max < 0.30 with d ≤ 0.15): "
          f"{paired_passes}/{paired_total} = {f_delta_rate*100 if paired_total else float('nan'):.1f}% — "
          f"{'HOLDS' if (paired_total and f_delta_rate >= 0.60) else 'FAILS / N/A'}")

    # F_ε — peak localisation.
    peak_top3 = sorted(rows, key=lambda r: r["rel_emp"], reverse=True)[:3]
    f_eps_holds = all(
        ((r["alpha"] < 0.10 or r["alpha"] > 0.90) and r["sigma_norm"] < 0.5)
        for r in peak_top3
    )
    results["F_eps_top3"] = [
        {"m": r["m"], "rel_emp": r["rel_emp"], "alpha": r["alpha"], "sigma_norm": r["sigma_norm"]}
        for r in peak_top3
    ]
    results["F_eps_holds"] = bool(f_eps_holds)
    print(f"F_ε (top-3 peaks have α < 0.10 or > 0.90 AND σ_norm < 0.5):")
    for r in peak_top3:
        print(f"     m={r['m']:>3}  rel={r['rel_emp']:>7.3f}  α={r['alpha']:.3f}  σn={r['sigma_norm']:.3f}")
    print(f"     → {'HOLDS' if f_eps_holds else 'FAILS'}")

    # F_ζ — trough localisation (excluding flat cells).
    non_flat = [r for r in rows if r["sigma_norm"] < 2.0]
    trough_bot3 = sorted(non_flat, key=lambda r: r["rel_emp"])[:3]
    f_zeta_holds = all(
        ((0.30 <= r["alpha"] <= 0.70) or r["sigma_norm"] > 1.0)
        for r in trough_bot3
    )
    results["F_zeta_bot3"] = [
        {"m": r["m"], "rel_emp": r["rel_emp"], "alpha": r["alpha"], "sigma_norm": r["sigma_norm"]}
        for r in trough_bot3
    ]
    results["F_zeta_holds"] = bool(f_zeta_holds)
    print(f"F_ζ (bot-3 troughs (σn<2) have α ∈ [0.30, 0.70] OR σ_norm > 1.0):")
    for r in trough_bot3:
        print(f"     m={r['m']:>3}  rel={r['rel_emp']:>7.4f}  α={r['alpha']:.3f}  σn={r['sigma_norm']:.3f}")
    print(f"     → {'HOLDS' if f_zeta_holds else 'FAILS'}")

    # F_η — wrapped-Gaussian-lim vs Gaussian-Y mean abs error within 30%.
    f_eta_holds = abs(err_Y - err_PL) / max(err_Y, eps) <= 0.30
    results["F_eta_holds"] = bool(f_eta_holds)
    print(f"F_η (Phase-Lim vs Gaussian-Y err agreement ≤ 30%): "
          f"err_Y={err_Y:.4f}, err_PL={err_PL:.4f}, "
          f"|Δ|/err_Y = {abs(err_Y - err_PL)/err_Y*100:.1f}% → "
          f"{'HOLDS' if f_eta_holds else 'FAILS'}")

    # Sanity: Poisson-form vs direct on a couple of cells.
    sanity = []
    for r in rows[:5]:
        ps = phase_lim_poisson(r["alpha"], r["sigma_norm"])
        sanity.append({"m": r["m"], "phase_lim": r["rel_pred_phase_lim"], "poisson": ps,
                       "diff": ps - r["rel_pred_phase_lim"]})
    results["poisson_sanity"] = sanity

    with open(args.out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
