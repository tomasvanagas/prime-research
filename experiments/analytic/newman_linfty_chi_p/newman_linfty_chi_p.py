"""
ATTACK_VECTORS §D27 — Newman / Erdős L^infty-flatness of the chi_P-coefficient
polynomial on |z|=1.

Polynomial: f_N(z) = Σ_{n=2}^{N} chi_P(n) z^n
Quantity:   R_N := ||f_N||_∞ / √π(N), where ||f||_∞ = max_{|z|=1} |f(z)|.

Three falsifiers (per ATTACK_VECTORS):
  (a)  R_N(prime) - R_B = O(1) → B-grade closure (7th orthogonal Fourier-side measure)
  (b)  R_N(prime) → c < R_RS = √2 → A-grade (super-Rudin-Shapiro flat)
  (c)  R_N(prime) >> R_B with growth → I-mode partial closure (HL singular series imprint)

Computation: FFT-oversample at M = 16N points on |z|=1; verify at smallest N
with M = 32N. Compare to ensembles: Bernoulli matched-density, signed primes
(Littlewood), Liouville, Rudin-Shapiro at matched length π(N).

Usage:
  python newman_linfty_chi_p.py --Nlist 14,16,18 --n_seeds 100
"""

import argparse
import json
import time
from math import gcd as _gcd
from pathlib import Path

import numpy as np
from sympy import sieve


def primes_up_to(N: int) -> np.ndarray:
    """Return numpy array of primes p in [2, N]."""
    return np.array(list(sieve.primerange(2, N + 1)), dtype=np.int64)


def chi_p_indicator(N: int) -> np.ndarray:
    """Return binary array v[0..N], v[n] = 1 iff n is prime."""
    v = np.zeros(N + 1, dtype=np.int8)
    for p in sieve.primerange(2, N + 1):
        v[p] = 1
    return v


def liouville_array(N: int) -> np.ndarray:
    """Return Liouville lambda(n) ∈ {-1, +1} for n in [0, N], with lambda(0)=0."""
    # lambda(n) = (-1)^Omega(n), Omega = number of prime factors with multiplicity.
    Omega = np.zeros(N + 1, dtype=np.int32)
    for p in sieve.primerange(2, N + 1):
        pk = p
        while pk <= N:
            Omega[pk::pk] += 1
            pk *= p
    lam = np.where(np.arange(N + 1) >= 1, (-1) ** Omega, 0)
    return lam.astype(np.float64)


def rudin_shapiro(M: int) -> np.ndarray:
    """Rudin-Shapiro sequence of length M. Returns ±1 array, P_k coefficients
    where the recursive RS pair P_{k+1}=P_k+x^{2^k}Q_k, Q_{k+1}=P_k-x^{2^k}Q_k.
    Truncate / extend by 0 to length M.
    Key property: ||P_k||_inf <= sqrt(2 * 2^k).
    """
    # Build P, Q until length >= M, then truncate.
    P = np.array([1], dtype=np.int8)
    Q = np.array([1], dtype=np.int8)
    while len(P) < M:
        P_new = np.concatenate([P, Q])
        Q_new = np.concatenate([P, -Q])
        P, Q = P_new, Q_new
    return P[:M].astype(np.float64)


def linf_via_fft(coef: np.ndarray, oversample: int = 16) -> tuple[float, int]:
    """Compute max_{|z|=1} |Σ coef[n] z^n| via FFT on (oversample * len(coef)) pts.
    Returns (max_abs, argmax_index).
    coef: real-valued coefficient array (degree = len-1).
    """
    n = len(coef)
    M = oversample * n
    F = np.fft.fft(coef, n=M)
    abs_F = np.abs(F)
    j = int(np.argmax(abs_F))
    return float(abs_F[j]), j


def linf_off_dc(coef: np.ndarray, oversample: int = 16,
                exclude_radius: int = 1) -> tuple[float, int]:
    """Off-DC L^infty: max_{|z|=1, z != 1} |sum coef[n] z^n|, computed via FFT.
    Excludes a radius of `exclude_radius` indices around k=0 (and the symmetric k=M-1
    for real coef the conjugate is at the same magnitude, no exclusion needed there).
    For 0/1 coefficient polynomials, the DC = sum(coef); excluding it isolates the
    Newman-flatness-style L^infty deviation away from constructive-interference at z=1.
    """
    n = len(coef)
    M = oversample * n
    F = np.fft.fft(coef, n=M)
    abs_F = np.abs(F)
    # Mask k in [0, exclude_radius] and [M-exclude_radius, M-1].
    mask = np.ones(M, dtype=bool)
    mask[:exclude_radius + 1] = False
    mask[M - exclude_radius:] = False
    j_local = int(np.argmax(abs_F[mask]))
    # Map back to original index.
    valid_idx = np.where(mask)[0]
    j = int(valid_idx[j_local])
    return float(abs_F[j]), j


def centered_linf(coef: np.ndarray, oversample: int = 16) -> tuple[float, int]:
    """Compute ||g||_inf where g(z) = sum_{n}(coef[n] - mean) z^n; n>=2 only on coef.
    Equivalent to subtracting mean*indicator(n>=2). Returns (max_abs, argmax)."""
    if len(coef) <= 2:
        return 0.0, 0
    nonzero_support = np.arange(len(coef)) >= 2
    p_hat = float(coef[nonzero_support].mean())
    centered = coef.copy()
    centered[nonzero_support] = centered[nonzero_support] - p_hat
    return linf_via_fft(centered, oversample=oversample)


def major_arc_indices(M: int, Q_max: int, half_arc_radius: int = 0) -> np.ndarray:
    """Return sorted unique FFT indices k corresponding to major arcs:
    k ≈ M * a/q for 1 <= q <= Q_max, 0 <= a < q with gcd(a,q)=1, plus q=1 (k=0).
    Each major-arc point is widened by ±half_arc_radius indices.
    """
    from math import gcd
    s = set()
    for q in range(1, Q_max + 1):
        for a in range(0, q):
            if q == 1 or gcd(a, q) == 1:
                k0 = (a * M) // q  # nearest index
                # Widen
                for r in range(-half_arc_radius, half_arc_radius + 1):
                    s.add((k0 + r) % M)
    return np.array(sorted(s), dtype=np.int64)


def closest_rational(freq_over_M: float, Q_max: int) -> tuple[int, int, float]:
    """Find the closest rational a/q with q <= Q_max to a given frequency in [0,1).
    Returns (a, q, |freq - a/q|)."""
    from math import gcd
    best = (0, 1, 1.0)
    for q in range(1, Q_max + 1):
        for a in range(0, q + 1):
            if q == 1 or gcd(a, q) == 1:
                d = min(abs(freq_over_M - a / q), abs(freq_over_M - 1 + a / q))
                if d < best[2]:
                    best = (a, q, d)
    return best


def minor_arc_linf(coef: np.ndarray, oversample: int = 16, Q_max: int = 16,
                   half_arc_radius: int = 1) -> tuple[float, int, dict]:
    """L^infty over minor arcs: exclude FFT indices in a neighborhood of every
    rational point a/q with q <= Q_max. Also report the absolute values
    AT each major arc point, broken down by q.
    Returns (max_abs_minor, argmax_index_minor, major_arc_breakdown)
    where major_arc_breakdown[q] = list of (a, k_index, |F[k]|).
    """
    from math import gcd
    n = len(coef)
    M = oversample * n
    F = np.fft.fft(coef, n=M)
    abs_F = np.abs(F)

    # Build major-arc index set with widening
    excluded = major_arc_indices(M, Q_max, half_arc_radius=half_arc_radius)
    mask = np.ones(M, dtype=bool)
    mask[excluded] = False

    # Argmax over minor arcs
    valid_idx = np.where(mask)[0]
    j_local = int(np.argmax(abs_F[valid_idx]))
    j = int(valid_idx[j_local])
    max_minor = float(abs_F[j])

    # Major-arc breakdown (without widening — just the central point)
    breakdown = {}
    for q in range(1, Q_max + 1):
        rows = []
        for a in range(0, q):
            if q == 1 or gcd(a, q) == 1:
                k0 = (a * M) // q
                rows.append({"a": a, "k": int(k0), "abs": float(abs_F[k0])})
        breakdown[q] = rows
    return max_minor, j, breakdown


def run_oversample_check(coef: np.ndarray, oversamples=(8, 16, 32, 64)) -> dict:
    """Check ||f||_inf at multiple oversample factors for one polynomial.
    Returns dict mapping oversample factor → max_abs.
    """
    results = {}
    for o in oversamples:
        m, _ = linf_via_fft(coef, oversample=o)
        results[o] = m
    return results


def primes_polynomial_coef(N: int) -> tuple[np.ndarray, int]:
    """Coefficients of f_N(z) = sum_{n=2..N} chi_P(n) z^n, length N+1, plus weight π(N)."""
    v = chi_p_indicator(N)
    weight = int(v.sum())
    return v.astype(np.float64), weight


def bernoulli_matched(N: int, weight: int, rng: np.random.Generator) -> np.ndarray:
    """Random binary array of length N+1; b[n] = 1 with probability weight/(N-1) for n>=2."""
    b = np.zeros(N + 1, dtype=np.float64)
    p = weight / (N - 1)
    b[2:] = (rng.random(N - 1) < p).astype(np.float64)
    return b


def signed_primes(prime_indicator: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Random ε_n ∈ {±1} on the prime support; 0 elsewhere."""
    out = prime_indicator.copy()
    mask = out > 0
    n_primes = int(mask.sum())
    out[mask] = rng.choice([-1.0, 1.0], size=n_primes)
    return out


def rudin_shapiro_at_weight(weight: int) -> np.ndarray:
    """Rudin-Shapiro sequence of length weight (so ||P||_inf / sqrt(weight) ≈ sqrt(2))."""
    return rudin_shapiro(weight)


def liouville_indicator_window(N: int) -> tuple[np.ndarray, int]:
    """f_λ_N(z) = sum_{n=2..N} λ(n) z^n. Coefficient array length N+1.
    Weight (variance proxy) is N-1 since λ(n) ∈ ±1, all n>=2.
    """
    lam = liouville_array(N)
    coef = np.zeros(N + 1, dtype=np.float64)
    coef[2:] = lam[2:]
    return coef, N - 1


def _ensemble_stats(values: list[float]) -> dict:
    arr = np.array(values, dtype=np.float64)
    return {
        "mean": float(arr.mean()),
        "std": float(arr.std(ddof=1)) if len(arr) > 1 else 0.0,
        "min": float(arr.min()),
        "max": float(arr.max()),
        "median": float(np.median(arr)),
        "n": int(len(arr)),
    }


def screen(N_log2: int, n_seeds: int, oversample: int = 16, rng_seed: int = 0) -> dict:
    """Run the full screening at N = 2^N_log2.
    Three notions:
      bare:   R = ||f||_inf / sqrt(weight)               (DC-dominated for 0/1)
      offdc:  R_off = max_{k!=0} |F[k]| / sqrt(weight)  (Newman flatness off DC)
      cent:   R_cent = ||g||_inf / sqrt(weight) where g = f − p_hat (centered)
    """
    N = 2**N_log2
    out = {"N": N, "log2_N": N_log2, "oversample": oversample}

    t0 = time.time()
    coef_p, weight_p = primes_polynomial_coef(N)
    out["weight_primes"] = weight_p

    # half_arc: ~oversample is the central Dirichlet peak; use 4×oversample to
    # suppress nearby side lobes. The total exclusion is ≲ Σφ(q)·(2·half_arc+1).
    half_arc = 4 * oversample
    # Total FFT length.
    M_fft = oversample * (N + 1)
    # Q_max chosen so total exclusion ≤ M_fft/4.
    # Σ_{q=1..Q} φ(q) ≈ 3Q²/π² ≈ 0.3 Q²
    safe_cap = max(2, int(np.sqrt(M_fft / (8 * (2 * half_arc + 1)) / 0.3)))
    Q_max = min(64, safe_cap)
    out["Q_max_used"] = Q_max
    out["safe_cap"] = safe_cap

    # Primes: bare, off-DC, centered, minor-arc
    bare_p, j_bare_p = linf_via_fft(coef_p, oversample=oversample)
    off_p, j_off_p = linf_off_dc(coef_p, oversample=oversample, exclude_radius=1)
    cent_p, j_cent_p = centered_linf(coef_p, oversample=oversample)
    minor_p, j_minor_p, breakdown_p = minor_arc_linf(coef_p, oversample=oversample,
                                                     Q_max=Q_max, half_arc_radius=half_arc)
    out["primes"] = {
        "linf_bare": bare_p,
        "linf_offdc": off_p,
        "linf_centered": cent_p,
        "linf_minor_arc": minor_p,
        "argmax_offdc_index": j_off_p,
        "argmax_offdc_freq_over_M": j_off_p / (oversample * (N + 1)),
        "argmax_minor_index": j_minor_p,
        "argmax_minor_freq_over_M": j_minor_p / (oversample * (N + 1)),
        "R_bare": bare_p / np.sqrt(weight_p),
        "R_offdc": off_p / np.sqrt(weight_p),
        "R_centered": cent_p / np.sqrt(weight_p),
        "R_minor_arc": minor_p / np.sqrt(weight_p),
        "major_arc_breakdown": {str(q): rows for q, rows in breakdown_p.items()},
    }
    out["t_primes_sec"] = time.time() - t0

    # Bernoulli matched-density ensemble (uniform support)
    rng = np.random.default_rng(rng_seed)
    bare_B_list, off_B_list, cent_B_list, minor_B_list = [], [], [], []
    t0 = time.time()
    for s in range(n_seeds):
        coef_b = bernoulli_matched(N, weight_p, rng)
        wb = int(coef_b.sum()) or 1
        bare_b, _ = linf_via_fft(coef_b, oversample=oversample)
        off_b, _ = linf_off_dc(coef_b, oversample=oversample, exclude_radius=1)
        cent_b, _ = centered_linf(coef_b, oversample=oversample)
        minor_b, _, _ = minor_arc_linf(coef_b, oversample=oversample,
                                       Q_max=Q_max, half_arc_radius=half_arc)
        bare_B_list.append(bare_b / np.sqrt(wb))
        off_B_list.append(off_b / np.sqrt(wb))
        cent_B_list.append(cent_b / np.sqrt(wb))
        minor_B_list.append(minor_b / np.sqrt(wb))
    out["bernoulli"] = {
        "R_bare": _ensemble_stats(bare_B_list),
        "R_offdc": _ensemble_stats(off_B_list),
        "R_centered": _ensemble_stats(cent_B_list),
        "R_minor_arc": _ensemble_stats(minor_B_list),
    }
    out["t_bernoulli_sec"] = time.time() - t0

    # Bernoulli matched-density on ODD support only (parity-matched control).
    rng = np.random.default_rng(rng_seed + 7)
    odd_count = (N - 1) // 2  # n in {3,5,...,N or N-1}
    p_on_odd = (weight_p - 1) / odd_count  # density of primes on odd support, excluding p=2
    minor_BO_list, off_BO_list, cent_BO_list = [], [], []
    t0 = time.time()
    for s in range(n_seeds):
        coef_bo = np.zeros(N + 1, dtype=np.float64)
        odd_mask = np.zeros(N + 1, dtype=bool)
        odd_mask[3::2] = True  # n = 3, 5, ..., up to <= N
        coef_bo[odd_mask] = (rng.random(odd_mask.sum()) < p_on_odd).astype(np.float64)
        coef_bo[2] = 1.0  # always include n=2 like primes
        wbo = int(coef_bo.sum()) or 1
        off_bo, _ = linf_off_dc(coef_bo, oversample=oversample, exclude_radius=1)
        cent_bo, _ = centered_linf(coef_bo, oversample=oversample)
        minor_bo, _, _ = minor_arc_linf(coef_bo, oversample=oversample,
                                        Q_max=Q_max, half_arc_radius=half_arc)
        off_BO_list.append(off_bo / np.sqrt(wbo))
        cent_BO_list.append(cent_bo / np.sqrt(wbo))
        minor_BO_list.append(minor_bo / np.sqrt(wbo))
    out["bernoulli_odd_support"] = {
        "p_on_odd": p_on_odd,
        "R_offdc": _ensemble_stats(off_BO_list),
        "R_centered": _ensemble_stats(cent_BO_list),
        "R_minor_arc": _ensemble_stats(minor_BO_list),
    }
    out["t_bernoulli_odd_sec"] = time.time() - t0

    # Signed primes (Littlewood) ensemble
    rng = np.random.default_rng(rng_seed + 1)
    R_L_list = []
    t0 = time.time()
    for s in range(n_seeds):
        coef_l = signed_primes(coef_p, rng)
        norm_l, _ = linf_via_fft(coef_l, oversample=oversample)
        R_L_list.append(norm_l / np.sqrt(weight_p))
    out["signed_primes"] = _ensemble_stats(R_L_list)
    out["t_signed_primes_sec"] = time.time() - t0

    # Liouville
    t0 = time.time()
    coef_lam, weight_lam = liouville_indicator_window(N)
    norm_lam, j_lam = linf_via_fft(coef_lam, oversample=oversample)
    out["liouville"] = {
        "linf": norm_lam,
        "argmax_index": j_lam,
        "argmax_freq_over_M": j_lam / (oversample * (N + 1)),
        "R_bare": norm_lam / np.sqrt(weight_lam),  # all signs ±1, no DC
        "weight": weight_lam,
    }
    out["t_liouville_sec"] = time.time() - t0

    # Rudin-Shapiro at matched weight π(N)
    t0 = time.time()
    coef_rs = rudin_shapiro_at_weight(weight_p)
    norm_rs, _ = linf_via_fft(coef_rs, oversample=oversample)
    out["rudin_shapiro"] = {
        "linf": norm_rs,
        "weight": weight_p,
        "R": norm_rs / np.sqrt(weight_p),
        "theoretical_bound_sqrt2": float(np.sqrt(2)),
    }
    out["t_rs_sec"] = time.time() - t0

    # Z-scores
    def _z(emp, ens):
        if ens["std"] > 0:
            return float((emp - ens["mean"]) / ens["std"])
        return float("nan")

    out["z_primes_vs_bernoulli_offdc"] = _z(out["primes"]["R_offdc"], out["bernoulli"]["R_offdc"])
    out["z_primes_vs_bernoulli_centered"] = _z(out["primes"]["R_centered"], out["bernoulli"]["R_centered"])
    out["z_primes_vs_bernoulli_bare"] = _z(out["primes"]["R_bare"], out["bernoulli"]["R_bare"])
    out["z_primes_vs_bernoulli_minor"] = _z(out["primes"]["R_minor_arc"], out["bernoulli"]["R_minor_arc"])
    out["z_primes_vs_signed_primes"] = _z(out["primes"]["R_offdc"], out["signed_primes"])
    out["z_primes_vs_bernoulli_odd_minor"] = _z(out["primes"]["R_minor_arc"], out["bernoulli_odd_support"]["R_minor_arc"])
    out["z_primes_vs_bernoulli_odd_centered"] = _z(out["primes"]["R_centered"], out["bernoulli_odd_support"]["R_centered"])
    out["z_primes_vs_bernoulli_odd_offdc"] = _z(out["primes"]["R_offdc"], out["bernoulli_odd_support"]["R_offdc"])

    # Salem-Zygmund analytic baseline: ±1 polys give ||f||_inf ≈ sqrt(2N log N), so
    # R = ||f||_inf / sqrt(N) ≈ sqrt(2 log N). For 0/1 binary off-DC, similar (with
    # variance of each coef = p(1-p) instead of 1).
    out["salem_zygmund_R_pm1"] = float(np.sqrt(2 * np.log(N + 1)))

    # Q_max scan: confirm minor-arc R is monotonically decreasing as Q_max grows
    # (more major arcs excluded → less arithmetic structure remains).
    # We compare primes vs uniform Bernoulli at each Q_max.
    out["q_max_scan"] = {}
    n_scan_seeds = max(10, n_seeds // 3)
    Q_list_full = [1, 2, 4, 8, 16, 32, 64, 128]
    Q_list = [q for q in Q_list_full if q <= safe_cap]
    for Q_scan in Q_list:
        prime_minor_R, _, _ = minor_arc_linf(coef_p, oversample=oversample,
                                             Q_max=Q_scan, half_arc_radius=half_arc)
        prime_R = prime_minor_R / np.sqrt(weight_p)
        bern_R_list = []
        rng_inner = np.random.default_rng(rng_seed + 200 + Q_scan)
        for s in range(n_scan_seeds):
            coef_b = bernoulli_matched(N, weight_p, rng_inner)
            wb = int(coef_b.sum()) or 1
            bm, _, _ = minor_arc_linf(coef_b, oversample=oversample,
                                      Q_max=Q_scan, half_arc_radius=half_arc)
            bern_R_list.append(bm / np.sqrt(wb))
        bern_stats = _ensemble_stats(bern_R_list)
        out["q_max_scan"][Q_scan] = {
            "prime_R_minor": prime_R,
            "bern_R_minor": bern_stats,
            "z": _z(prime_R, bern_stats),
        }

    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Nlist", type=str, default="14,16,18",
                    help="Comma-separated log2(N) values")
    ap.add_argument("--n_seeds", type=int, default=50)
    ap.add_argument("--oversample", type=int, default=16)
    ap.add_argument("--rng_seed", type=int, default=42)
    ap.add_argument("--outdir", type=str,
                    default=str(Path(__file__).parent))
    ap.add_argument("--oversample_check_at", type=int, default=14,
                    help="At log2(N) = this value, run multi-oversample validation.")
    args = ap.parse_args()

    Ns = [int(x) for x in args.Nlist.split(",")]
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    all_results = []
    for log2N in Ns:
        print(f"=== N = 2^{log2N} = {2**log2N} ===", flush=True)
        res = screen(log2N, n_seeds=args.n_seeds,
                     oversample=args.oversample, rng_seed=args.rng_seed)
        all_results.append(res)
        p, b, bo = res["primes"], res["bernoulli"], res["bernoulli_odd_support"]
        print(f"  weight (=π(N)) = {res['weight_primes']}", flush=True)
        print(f"  primes:        R_bare={p['R_bare']:.4f}  R_offdc={p['R_offdc']:.4f}  R_cent={p['R_centered']:.4f}  R_minor={p['R_minor_arc']:.4f}", flush=True)
        print(f"  bern_uniform:  R_bare={b['R_bare']['mean']:.4f}±{b['R_bare']['std']:.4f}  R_offdc={b['R_offdc']['mean']:.4f}±{b['R_offdc']['std']:.4f}  R_cent={b['R_centered']['mean']:.4f}±{b['R_centered']['std']:.4f}  R_minor={b['R_minor_arc']['mean']:.4f}±{b['R_minor_arc']['std']:.4f}", flush=True)
        print(f"  bern_odd_supp: R_offdc={bo['R_offdc']['mean']:.4f}±{bo['R_offdc']['std']:.4f}  R_cent={bo['R_centered']['mean']:.4f}±{bo['R_centered']['std']:.4f}  R_minor={bo['R_minor_arc']['mean']:.4f}±{bo['R_minor_arc']['std']:.4f}", flush=True)
        print(f"  signed_primes (Littlewood): R = {res['signed_primes']['mean']:.4f} ± {res['signed_primes']['std']:.4f}", flush=True)
        print(f"  liouville:    R = {res['liouville']['R_bare']:.4f}", flush=True)
        print(f"  rudin_shapiro:R = {res['rudin_shapiro']['R']:.4f}  (sqrt(2)={float(np.sqrt(2)):.4f})", flush=True)
        print(f"  Salem-Zyg(±1):R ≈ {res['salem_zygmund_R_pm1']:.4f}", flush=True)
        print(f"  Z[primes vs uniform-Bern]:   bare={res['z_primes_vs_bernoulli_bare']:+.2f}  offdc={res['z_primes_vs_bernoulli_offdc']:+.2f}  cent={res['z_primes_vs_bernoulli_centered']:+.2f}  minor={res['z_primes_vs_bernoulli_minor']:+.2f}", flush=True)
        print(f"  Z[primes vs odd-Bern]:       offdc={res['z_primes_vs_bernoulli_odd_offdc']:+.2f}  cent={res['z_primes_vs_bernoulli_odd_centered']:+.2f}  minor={res['z_primes_vs_bernoulli_odd_minor']:+.2f}", flush=True)
        print(f"  primes argmax (off-DC) freq/M = {p['argmax_offdc_freq_over_M']:.6f}  primes argmax (minor) freq/M = {p['argmax_minor_freq_over_M']:.6f}", flush=True)
        a, q, d = closest_rational(p['argmax_minor_freq_over_M'], Q_max=200)
        print(f"  primes argmax (minor) closest rational a/q = {a}/{q} (distance {d:.6f})", flush=True)
        # Major-arc breakdown for primes
        Q_used = res["Q_max_used"]
        print(f"  Q_max used (main minor-arc): {Q_used}; safe_cap = {res['safe_cap']}", flush=True)
        print(f"  Major arc R values for primes (|F[a/q]| / sqrt(weight)):", flush=True)
        for q in range(1, Q_used + 1):
            rows = p["major_arc_breakdown"].get(str(q), [])
            R_vals = [r["abs"] / np.sqrt(res["weight_primes"]) for r in rows]
            R_str = " ".join(f"{r:.3f}" for r in R_vals)
            phi_q = sum(1 for a in range(q) if (q == 1 or _gcd(a, q) == 1))
            theory = np.sqrt(res["weight_primes"]) / max(phi_q, 1)
            print(f"     q={q:2d}  φ(q)={phi_q:2d}  HL_pred={theory:.3f}: {R_str}", flush=True)
        print(f"  liouville argmax freq/M = {res['liouville']['argmax_freq_over_M']:.6f}", flush=True)
        print(f"  Q_max scan (R_minor for primes vs uniform-Bern):", flush=True)
        for Q_scan, qres in res["q_max_scan"].items():
            print(f"     Q_max={Q_scan:4d}: prime={qres['prime_R_minor']:.4f}  bern={qres['bern_R_minor']['mean']:.4f}±{qres['bern_R_minor']['std']:.4f}  z={qres['z']:+.2f}", flush=True)
        print(f"  times: primes {res['t_primes_sec']:.1f}s; bernoulli {res['t_bernoulli_sec']:.1f}s; signed {res['t_signed_primes_sec']:.1f}s; lambda {res['t_liouville_sec']:.1f}s", flush=True)

    # Oversample validation at smallest tested N
    if args.oversample_check_at in Ns:
        print(f"--- Oversample validation at N = 2^{args.oversample_check_at} ---", flush=True)
        N = 2**args.oversample_check_at
        coef_p, _ = primes_polynomial_coef(N)
        ovs = [8, 16, 32, 64, 128]
        oversamp_results = {}
        for o in ovs:
            m, _ = linf_via_fft(coef_p, oversample=o)
            oversamp_results[o] = m
            print(f"  oversample={o}: ||f||_inf = {m:.6f}", flush=True)
        all_results.append({"oversample_validation": oversamp_results,
                            "N_check": N})

    out_path = outdir / "results.json"
    with open(out_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
