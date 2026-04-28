"""D43.c — K-truncated explicit-formula residual roughness.

Hypothesis: as K → ∞, the residual
   R_K(x) := π(x) − Li(x) + Σ_{k≤K} 2 Re Li(x^{ρ_k})
becomes Hölder-α with α* < 1/2 (KPZ-class roughness localised to high zeros).
Sign convention: the explicit formula is π(x) ≈ Li(x) − Σ 2 Re Li(x^{ρ_k}) − log 2,
so π(x) − Li(x) ≈ −Σ 2 Re Li(x^{ρ_k}) − log 2. To *remove* the zero
contribution we ADD the sum: R_K = (π − Li) + Σ. The challenge spec
(NOVELTY_CHALLENGES.md §D43.c) had the sign reversed; corrected here.
S157 (E2.27) measured α(D_0) ≈ 0.85 — far above KPZ. We test whether the
smoothness is uniformly distributed in γ or concentrated in low frequencies.

Methodology: see `d43c_k_truncated_residual_results.md` for pre-stated
falsifiers F1-F6.

Pre-states:
  F1 — α(D_K) ≤ 0.5 at K=1000 (A-grade KPZ residual signal).
  F2 — α(D_K) ≥ 0.75 for all K (B-grade smoothness-preserved refinement).
  F3 — α drops below 0.5 at some K* finite (A-grade crossover).
  F4 — α(D_0) ∈ [0.75, 0.95] (sanity vs S157).
  F5 — Var(R_K) monotone decreasing.
  F6 — |Δ_asymp − Δ_mpmath| / |Δ_mpmath| ≤ 1e-6 at K=1, x_mid.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np
import pywt

HERE = Path(__file__).parent
PROJECT_ROOT = HERE.parent.parent.parent
ZEROS_FILE = PROJECT_ROOT / "data" / "zeta_zeros_8000.txt"

# Reuse S157 utilities.
sys.path.insert(0, str(PROJECT_ROOT / "experiments" / "analytic" / "kpz_pi_li_d43"))
from d43_kpz_pi_li import sieve_pi_table, kpz_grid, wavelet_holder, cramer_pi  # noqa: E402


def wavelet_holder_band(D, levels_keep, wavelet="db4", max_level=None):
    """Hölder fit restricted to the SUBSET of pywt levels in `levels_keep`.

    Useful for excluding scales that have been zeroed-out by explicit-formula
    truncation (band-pass filtering). Returns dict with alpha, r², and the
    levels actually fit.
    """
    n = len(D)
    if max_level is None:
        max_level = pywt.dwt_max_level(n, pywt.Wavelet(wavelet).dec_len)
        max_level = min(max_level, 12)
    coeffs = pywt.wavedec(D, wavelet, level=max_level)
    raw = []
    for level_idx in range(1, len(coeffs)):
        if level_idx not in levels_keep:
            continue
        c = coeffs[level_idx]
        E = float(np.mean(c ** 2))
        raw.append((level_idx, E, len(c)))
    if len(raw) < 3:
        return None
    levels = np.array([li for li, _, _ in raw], dtype=np.float64)
    log_E = np.array([math.log2(max(E, 1e-300)) for _, E, _ in raw])
    A = np.vstack([levels, np.ones_like(levels)]).T
    sol, _, _, _ = np.linalg.lstsq(A, log_E, rcond=None)
    slope, intercept = sol
    pred = A @ sol
    ss_res = ((log_E - pred) ** 2).sum()
    ss_tot = ((log_E - log_E.mean()) ** 2).sum()
    r2 = 1.0 - ss_res / max(ss_tot, 1e-20)
    alpha = -float(slope) / 2.0
    return {
        "slope": float(slope), "intercept": float(intercept),
        "alpha_holder": alpha, "r2": float(r2),
        "level_energies": [(int(li), e, int(c)) for li, e, c in raw],
        "levels_used": [int(l) for l in levels],
    }


def load_zeros(path: Path, K_max: int) -> np.ndarray:
    """Load first K_max imaginary parts of nontrivial Riemann zeros."""
    out = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            try:
                out.append(float(s))
            except ValueError:
                continue
            if len(out) >= K_max:
                break
    if len(out) < K_max:
        raise RuntimeError(f"Only {len(out)} zeros available, need {K_max}")
    return np.asarray(out, dtype=np.float64)


def li_x_rho_asymp(log_x: np.ndarray, rho: complex,
                   n_terms: int = 6) -> np.ndarray:
    """Asymptotic series for Li(x^ρ) = Ei(ρ log x).

    Ei(z) ~ e^z / z · Σ_{n=0..N} n!/z^n  for |z| → ∞ (Olver §2.2 / A&S 5.1.51).

    Returns a complex numpy array.
    """
    z = rho * log_x          # complex array
    inv_z = 1.0 / z
    series = np.zeros_like(z)
    z_pow = np.ones_like(z)
    fact = 1.0
    for n in range(n_terms):
        if n > 0:
            fact *= n
        series += fact * z_pow
        z_pow = z_pow * inv_z
    return np.exp(z) * inv_z * series


def li_x_rho_mpmath(x: float, rho: complex) -> complex:
    """Reference: Li(x^rho) via mpmath.ei."""
    z = mp.mpc(rho.real, rho.imag) * mp.log(x)
    return complex(mp.ei(z))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--logX", type=int, default=22)
    parser.add_argument("--K_max", type=int, default=2000)
    parser.add_argument("--K_grid", type=str,
                        default="0,1,5,10,50,200,1000,2000")
    parser.add_argument("--n_terms", type=int, default=6,
                        help="Asymptotic series order for Li(x^rho).")
    parser.add_argument("--outdir", type=str, default=str(HERE))
    args = parser.parse_args()

    K_list = [int(s) for s in args.K_grid.split(",")]
    K_max = max(K_list + [args.K_max])

    X = 1 << args.logX
    print(f"[D43.c] X = 2^{args.logX} = {X}, K_max = {K_max}")

    t0 = time.time()
    pi_tab = sieve_pi_table(X)
    print(f"[sieve] pi({X}) = {pi_tab[X]}, {time.time()-t0:.1f}s")

    xs, step = kpz_grid(X)
    print(f"[grid] step={step}, n={len(xs)}")

    log_xs = np.log(xs.astype(np.float64))
    sqrt_xs = np.sqrt(xs.astype(np.float64))

    t0 = time.time()
    li_vals = np.array([float(mp.li(int(xv))) for xv in xs])
    print(f"[Li ref] {time.time()-t0:.1f}s")

    pi_vals = pi_tab[xs].astype(np.float64)
    pre = pi_vals - li_vals  # the "residual before zero subtraction"

    # Cramér control prime model.
    print("[cramer] sampling Bernoulli(1/log n) prime model")
    pi_C_tab = cramer_pi(X, seed=314159)
    pre_C = pi_C_tab[xs].astype(np.float64) - li_vals

    gammas = load_zeros(ZEROS_FILE, K_max)
    print(f"[zeros] loaded {len(gammas)}; γ_1={gammas[0]:.4f}, "
          f"γ_{len(gammas)}={gammas[-1]:.4f}")

    # F6: validate asymptotic vs mpmath at multiple (K, x) pairs.
    # We only need 2*Re(Ei) so we test the real part. The iπ Stokes-term
    # is purely imaginary and contributes nothing to 2 Re Ei (it cancels
    # between the conjugate pair).
    F6_samples = []
    for K_idx in [0, 9, 99, 999, min(len(gammas) - 1, 7999)]:
        for x_pos in [len(xs) // 4, len(xs) // 2, 3 * len(xs) // 4]:
            rho_test = complex(0.5, gammas[K_idx])
            x_test = float(xs[x_pos])
            log_x_test = math.log(x_test)
            asymp = float(np.real(li_x_rho_asymp(
                np.array([log_x_test]), rho_test, n_terms=args.n_terms)[0]))
            ref = float(mp.re(li_x_rho_mpmath(x_test, rho_test)))
            err = abs(asymp - ref)
            rel_err = err / max(abs(ref), 1e-300)
            F6_samples.append({
                "K_idx": K_idx, "gamma": gammas[K_idx], "x": x_test,
                "asymp_re": asymp, "ref_re": ref,
                "abs_err": err, "rel_err": rel_err,
            })
    rel_err_max = max(s["rel_err"] for s in F6_samples)
    abs_err_max = max(s["abs_err"] for s in F6_samples)
    print(f"[F6] over {len(F6_samples)} (K, x) samples:")
    for s in F6_samples:
        print(f"     γ={s['gamma']:9.3f}, x={s['x']:.0f}: "
              f"asymp={s['asymp_re']:+.4e}, ref={s['ref_re']:+.4e}, "
              f"rel_err={s['rel_err']:.2e}")
    print(f"[F6] worst rel_err on Re = {rel_err_max:.3e}, "
          f"worst abs_err = {abs_err_max:.3e}")
    rel_err = rel_err_max
    x_mid_idx = len(xs) // 2

    # Cumulatively ADD zero contributions to `pre = π − Li`.
    # Sign: π − Li ≈ −Σ 2 Re Li(x^ρ_k) − log 2, so adding cancels.
    cumulative = pre.copy()
    cumulative_C = pre_C.copy()
    R_dict: dict[int, np.ndarray] = {}
    R_C_dict: dict[int, np.ndarray] = {}
    var_dict: dict[int, float] = {}
    var_C_dict: dict[int, float] = {}
    if 0 in K_list:
        R_dict[0] = pre.copy()
        R_C_dict[0] = pre_C.copy()
        var_dict[0] = float(np.var(pre, ddof=1))
        var_C_dict[0] = float(np.var(pre_C, ddof=1))

    t0 = time.time()
    K_targets = sorted(set(K_list))
    if 0 in K_targets:
        K_targets.remove(0)
    K_done = 0
    for K_target in K_targets:
        for k in range(K_done, K_target):
            rho = complex(0.5, gammas[k])
            li_rho = li_x_rho_asymp(log_xs, rho, n_terms=args.n_terms)
            two_re = 2.0 * np.real(li_rho)
            # ADD to π's residual (explicit formula sign: see header).
            cumulative = cumulative + two_re
            # SAME procedure on Cramér: should NOT cancel since Cramér
            # has no explicit-formula structure.
            cumulative_C = cumulative_C + two_re
        K_done = K_target
        R_dict[K_target] = cumulative.copy()
        R_C_dict[K_target] = cumulative_C.copy()
        var_dict[K_target] = float(np.var(cumulative, ddof=1))
        var_C_dict[K_target] = float(np.var(cumulative_C, ddof=1))
        print(f"[K={K_target:5d}] var(R_K)={var_dict[K_target]:.4e}, "
              f"var(C_K)={var_C_dict[K_target]:.4e}, "
              f"elapsed {time.time()-t0:.1f}s")

    # Compute D_K and Hölder for each K (full fit + fine-scale-only fit).
    # Fine-scale levels: those NOT subtracted by the K-truncation.
    # We map K to a "subtraction cutoff level" via γ_K and the wavelet
    # frequency band per level.
    # In u-space, du ≈ X^{-2/3}; level_idx j_pywt has scale 2^{j_pywt}·du,
    # corresponding to γ_band ≈ 2π / (2^{j_pywt} · du).
    # γ < γ_K is "below subtraction frequency"; we want levels with
    # γ_band > γ_K (i.e., FINER scales than what we just subtracted).
    n_grid = len(xs)
    du = float(np.log(xs[-1] / xs[0]) / (n_grid - 1))  # mean du
    max_level_pywt = pywt.dwt_max_level(
        n_grid, pywt.Wavelet("db4").dec_len)
    max_level_pywt = min(max_level_pywt, 12)

    def level_to_gamma(j_pywt: int) -> float:
        """Approx γ-band centre for pywt level j_pywt (1=coarsest)."""
        scale = 2 ** j_pywt
        return 2 * math.pi / (scale * du)

    def fine_levels_for_K(K: int) -> list[int]:
        """Levels whose γ-band is ABOVE γ_K (not yet subtracted)."""
        if K == 0:
            return list(range(1, max_level_pywt + 1))
        gamma_K = float(gammas[K - 1])
        keep = []
        for j in range(1, max_level_pywt + 1):
            if level_to_gamma(j) > gamma_K:
                keep.append(j)
        return keep

    holder_dict: dict[int, dict] = {}
    holder_C_dict: dict[int, dict] = {}
    holder_fine_dict: dict[int, dict] = {}
    holder_fine_C_dict: dict[int, dict] = {}

    for K in sorted(R_dict):
        R = R_dict[K]
        R_C = R_C_dict[K]
        D_K = R * log_xs / sqrt_xs
        D_C_K = R_C * log_xs / sqrt_xs
        h = wavelet_holder(D_K, wavelet="db4")
        h_C = wavelet_holder(D_C_K, wavelet="db4")

        levels_keep = fine_levels_for_K(K)
        h_fine = wavelet_holder_band(D_K, set(levels_keep), wavelet="db4")
        h_fine_C = wavelet_holder_band(D_C_K, set(levels_keep), wavelet="db4")

        holder_dict[K] = {
            "alpha": float(h["alpha_holder"]),
            "r2": float(h["r2"]),
            "slope": float(h["slope"]),
            "D_mean": float(D_K.mean()),
            "D_std": float(D_K.std(ddof=1)),
            "Var(R_K)": var_dict[K],
        }
        holder_C_dict[K] = {
            "alpha": float(h_C["alpha_holder"]),
            "r2": float(h_C["r2"]),
            "D_std": float(D_C_K.std(ddof=1)),
            "Var(R_K)": var_C_dict[K],
        }
        holder_fine_dict[K] = (
            None if h_fine is None
            else {"alpha": float(h_fine["alpha_holder"]),
                  "r2": float(h_fine["r2"]),
                  "levels": h_fine["levels_used"],
                  "n_levels": len(h_fine["levels_used"])}
        )
        holder_fine_C_dict[K] = (
            None if h_fine_C is None
            else {"alpha": float(h_fine_C["alpha_holder"]),
                  "r2": float(h_fine_C["r2"]),
                  "levels": h_fine_C["levels_used"],
                  "n_levels": len(h_fine_C["levels_used"])}
        )

        fine_str = (f"FINE α={h_fine['alpha_holder']:+.4f} "
                    f"r²={h_fine['r2']:.4f} ({len(h_fine['levels_used'])} lvls)"
                    if h_fine is not None else "FINE n/a")
        print(f"[holder] K={K:5d} π:  FULL α={h['alpha_holder']:+.4f} r²={h['r2']:.4f}  "
              f"{fine_str}")
        print(f"          K={K:5d} C:  FULL α={h_C['alpha_holder']:+.4f} r²={h_C['r2']:.4f}  "
              f"FINE α={(h_fine_C['alpha_holder'] if h_fine_C else float('nan')):+.4f} "
              f"r²={(h_fine_C['r2'] if h_fine_C else float('nan')):.4f}")

    # Falsifier verdict
    alphas = {K: holder_dict[K]["alpha"] for K in sorted(holder_dict)}
    alphas_C = {K: holder_C_dict[K]["alpha"] for K in sorted(holder_C_dict)}
    alphas_fine = {K: (holder_fine_dict[K]["alpha"]
                       if holder_fine_dict[K] is not None else float("nan"))
                   for K in sorted(holder_fine_dict)}
    alphas_fine_C = {K: (holder_fine_C_dict[K]["alpha"]
                         if holder_fine_C_dict[K] is not None else float("nan"))
                     for K in sorted(holder_fine_C_dict)}

    K_max_done = max(alphas)
    F1_full = alphas.get(1000, alphas[K_max_done]) <= 0.50
    F2 = all(a >= 0.75 for a in alphas.values())
    K_below = [K for K, a in sorted(alphas.items()) if a <= 0.50]
    F3_full = bool(K_below) and all(alphas[K] <= 0.55 for K in alphas
                                     if K >= K_below[0])

    # Genuine-roughness criterion: fine-scale α with valid fit.
    # Use only K values where fine-scale fit is well-conditioned (r² > 0.5)
    # and Cramér control fine-scale α stays high.
    fine_valid = {K: holder_fine_dict[K] for K in holder_fine_dict
                  if holder_fine_dict[K] is not None
                  and holder_fine_dict[K]["r2"] > 0.5}
    F1_fine = any(holder_fine_dict[K] is not None
                   and holder_fine_dict[K]["alpha"] <= 0.5
                   and holder_fine_dict[K]["r2"] > 0.5
                  for K in [1000, 2000, 4000] if K in holder_fine_dict)

    # Cramér control: if subtraction is just noise, fine-scale α_C should
    # be similar to baseline (or even more positive). If π's residual
    # SPECIFICALLY drops below Cramér's, that's structural.
    pi_below_cramer_fine = any(
        holder_fine_dict[K] is not None and holder_fine_C_dict[K] is not None
        and holder_fine_dict[K]["alpha"] < holder_fine_C_dict[K]["alpha"] - 0.2
        and holder_fine_dict[K]["r2"] > 0.5
        for K in alphas_fine if K > 0
    )

    F4 = 0.75 <= alphas[0] <= 0.95
    var_seq = [var_dict[K] for K in sorted(var_dict)]
    F5 = all(var_seq[i] >= var_seq[i + 1] for i in range(len(var_seq) - 1))
    F6 = rel_err <= 1e-6

    out = {
        "params": {
            "logX": args.logX, "X": X, "step": step, "n_grid": len(xs),
            "K_grid": K_list, "n_terms": args.n_terms,
            "n_zeros_loaded": len(gammas),
            "du_mean": du,
            "max_level_pywt": int(max_level_pywt),
        },
        "F6_asymptotic_rel_err_max": float(rel_err),
        "F6_asymptotic_abs_err_max": float(abs_err_max),
        "F6_samples": F6_samples,
        "by_K_pi": {str(K): holder_dict[K] for K in sorted(holder_dict)},
        "by_K_cramer": {str(K): holder_C_dict[K] for K in sorted(holder_C_dict)},
        "by_K_pi_fine": {str(K): holder_fine_dict[K] for K in sorted(holder_fine_dict)},
        "by_K_cramer_fine": {str(K): holder_fine_C_dict[K] for K in sorted(holder_fine_C_dict)},
        "alpha_pi_by_K": {str(K): float(a) for K, a in alphas.items()},
        "alpha_cramer_by_K": {str(K): float(a) for K, a in alphas_C.items()},
        "alpha_fine_pi_by_K": {str(K): float(a) for K, a in alphas_fine.items()},
        "alpha_fine_cramer_by_K": {str(K): float(a) for K, a in alphas_fine_C.items()},
        "var_pi_by_K": {str(K): float(var_dict[K]) for K in sorted(var_dict)},
        "var_cramer_by_K": {str(K): float(var_C_dict[K]) for K in sorted(var_C_dict)},
        "F1_full_alpha_drops_to_kpz": bool(F1_full),
        "F1_fine_genuine_roughness": bool(F1_fine),
        "pi_below_cramer_in_fine": bool(pi_below_cramer_fine),
        "F2_smoothness_preserved_full": bool(F2),
        "F3_full_crossover": bool(F3_full),
        "F4_calibration": bool(F4),
        "F5_var_monotone": bool(F5),
        "F6_asymptotic_accurate": bool(F6),
    }
    out_path = Path(args.outdir) / "d43c_k_truncated_residual_results.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"[done] wrote {out_path}")

    print()
    print("=== FALSIFIER VERDICT ===")
    print(f"  F1 FULL (α≤0.50 at K=1000):    {F1_full}")
    print(f"  F1 FINE (genuine α<0.5 r²>0.5): {F1_fine}")
    print(f"  pi α_fine < cramer α_fine -0.2: {pi_below_cramer_fine}")
    print(f"  F2 (α_full≥0.75 for all K):    {F2}")
    print(f"  F3 FULL (crossover):           {F3_full}")
    print(f"  F4 (α_0 calibrated):           {F4}")
    print(f"  F5 (var(R_K) monotone):        {F5}")
    print(f"  F6 (asymptotic OK):            {F6}")
    print()
    print("  α-sequence — π   (full / fine) vs Cramér (full / fine):")
    for K in sorted(alphas):
        a = alphas[K]
        a_C = alphas_C[K]
        a_f = alphas_fine.get(K, float("nan"))
        a_f_C = alphas_fine_C.get(K, float("nan"))
        var_ratio = var_dict[K] / var_dict[0] if 0 in var_dict else 1.0
        var_ratio_C = var_C_dict[K] / var_C_dict[0] if 0 in var_C_dict else 1.0
        print(f"    K={K:5d}: π α={a:+.3f} | π α_fine={a_f:+.3f}    "
              f"C α={a_C:+.3f} | C α_fine={a_f_C:+.3f}    "
              f"var(R)/var(R_0)={var_ratio:.3f}, "
              f"var(C)/var(C_0)={var_ratio_C:.3f}")


if __name__ == "__main__":
    main()
