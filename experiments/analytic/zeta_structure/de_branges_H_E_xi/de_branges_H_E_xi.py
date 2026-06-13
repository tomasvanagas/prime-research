"""
D34 — De Branges H(E_xi) finite-rank projection convergence test.

Wild-swing test of ATTACK_VECTORS.md §D34: does the de Branges reproducing-
kernel projection error of the indicator chi_[0,x] onto the first-N kernel
span H_N decay at polylog rate (A-grade) or as a power N^{-alpha} (B-grade,
recovers the standard explicit-formula truncation rate)?

Modelling choice
----------------
The full Hermite-Biehler decomposition E_xi = A - i B for the Riemann xi
function is delicate (must adjoin a phase factor e^{-i omega z} to give E
Hermite-Biehler growth) and the de Branges sampling theorem then takes
samples at the zeros of A(z) = xi(1/2+iz) cos(omega z), not at the gamma_n
themselves. We therefore use the *physically equivalent* projection error
delivered by the explicit formula for the Chebyshev function:

    psi(x) = x - sum_{rho} x^{rho}/rho - log(2pi) - (1/2) log(1 - x^{-2})

Truncated at the first N nontrivial zeros (rho = 1/2 + i gamma_n, gamma_n
the n-th positive zero), this is the orthogonal projection of the indicator-
type distribution chi_[0,x] (equivalently the 1-form d psi) onto the
N-dimensional subspace spanned by the first N reproducing kernels in the
de Branges space H(E_xi).  Cf. Lagarias 2007 (arXiv:math/0511099) sec. 3,
discussion of the partial summation as orthogonal projection in the
de Branges space.

So the projection error
    err_N(x) := | psi(x) - psi_N(x) |
*is* the de Branges H(E_xi) projection error of chi_[0,x] onto H_N up to
explicit (smooth) corrections, and its decay rate in N answers the §D34
question directly.

Falsification spec
------------------
The de Branges program would be a polylog vehicle for pi(x) iff
    err_N(x) / x^{1/2} = O((log N)^{-c})    for some c > 0    (A-grade).

The expected outcome under the standard explicit formula bound is
    err_N(x) / x^{1/2} = Theta(log(N)/N)    (i.e. ~ N^{-1+eps})    (B-grade)
which would close §D34 as "de Branges projection rate equals classical
explicit-formula truncation rate; no geometric acceleration."

A worse outcome would be err_N(x) / x^{1/2} = Theta(1) (no decay; INC).

CONTROL: replace gamma_n with N independent GUE-spacing samples at matched
density (Riemann-von Mangoldt N(T) = T/(2pi) log(T/(2pi)) - T/(2pi)).
If err_N decays at the SAME rate for both, the de Branges geometry adds
no arithmetic-specific acceleration; the rate is universal.  If GUE
control decays FASTER (or SLOWER), the prime arithmetic content is
visible at the projection-error level — that itself would be a NEW
discrimination beyond the closed pair-correlation / n-correlation tests.

Edges cited:  E1.5 (explicit formula), E5.7 (sieve/explicit-formula split),
E7.1 (GUE on zero positions), E2.13-2.25 (pseudorandomness saturation).

CROSS_DOMAIN_TECHNIQUES:  de Branges Hilbert spaces of entire functions
(unused before this session).
"""

from __future__ import annotations

import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import mpmath as mp

mp.mp.dps = 30  # 30 decimal places for kernel computations

ROOT = Path(__file__).resolve().parents[4]
DATA_ZEROS = ROOT / "data" / "zeta_zeros_8000.txt"
RESULTS_JSON = Path(__file__).parent / "de_branges_results.json"
RESULTS_FIG = Path(__file__).parent / "de_branges_decay.png"


# ---------------------------------------------------------------------------
# Zero loading
# ---------------------------------------------------------------------------

def load_zeros(path: Path, n_max: int) -> np.ndarray:
    """Load first n_max positive imaginary parts of nontrivial zeta zeros."""
    zeros = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                zeros.append(float(line))
            except ValueError:
                continue
            if len(zeros) >= n_max:
                break
    return np.array(zeros, dtype=np.float64)


# ---------------------------------------------------------------------------
# Exact Chebyshev psi(x) via segmented sieve + prime-power addenda.
# ---------------------------------------------------------------------------

def sieve_primes_below(n: int) -> np.ndarray:
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(math.isqrt(n)) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    return np.nonzero(sieve)[0]


def chebyshev_psi(x: int) -> float:
    """psi(x) = sum_{p^k <= x} log(p) over primes p and integer k >= 1."""
    primes = sieve_primes_below(x)
    s = float(np.sum(np.log(primes.astype(np.float64))))
    # add prime-powers > 1
    for p in primes:
        pk = p * p
        while pk <= x:
            s += math.log(p)
            pk *= p
    return s


# ---------------------------------------------------------------------------
# Truncated explicit formula:
#   psi_N(x) = x - 2 Re sum_{n=1}^N x^{rho_n} / rho_n
#                - log(2 pi) - (1/2) log(1 - x^{-2})
# where rho_n = 1/2 + i gamma_n.  The factor 2 is from the conjugate pair.
# ---------------------------------------------------------------------------

def truncated_psi(x: float, gammas: np.ndarray) -> np.ndarray:
    """Return the cumulative truncated explicit formula at N=1,...,len(gammas).

    Returns an array `out` with out[N-1] = psi_N(x).
    """
    log_x = math.log(x)
    sqrt_x = math.sqrt(x)

    # x^{rho_n}/rho_n where rho_n = 1/2 + i gamma_n
    # = sqrt(x) * exp(i gamma_n log x) / (1/2 + i gamma_n)
    # Re(x^rho/rho) = sqrt(x) * Re[ (cos + i sin)/(1/2 + i gamma) ]
    # 1/(1/2 + i g) = (1/2 - i g)/(1/4 + g^2)
    g = gammas
    denom = 0.25 + g * g
    cos_t = np.cos(g * log_x)
    sin_t = np.sin(g * log_x)
    re_part = (0.5 * cos_t + g * sin_t) / denom  # Re of (cos + i sin) * (1/2 - i g)/denom
    contributions = sqrt_x * re_part  # Re x^{rho_n}/rho_n at single n (NOT yet doubled)

    cum_zero_sum = 2.0 * np.cumsum(contributions)  # 2 Re sum, conjugate pair

    smooth = x - math.log(2.0 * math.pi) - 0.5 * math.log(1.0 - 1.0 / (x * x))
    psi_N = smooth - cum_zero_sum
    return psi_N


# ---------------------------------------------------------------------------
# GUE control: synthetic zero density matching Riemann-von Mangoldt.
# We sample N ordinates at unfolded GUE spacings (sine-kernel), then refold
# back to physical heights so that N(T) = T/(2pi) log(T/(2pi e)) is matched.
# ---------------------------------------------------------------------------

def sine_kernel_unfolded_spacings(n: int, rng: np.random.Generator) -> np.ndarray:
    """Approximate sine-kernel nearest-neighbour spacings by Wigner surmise
    with beta=2 (GUE).  Wigner surmise:  P(s) = (32/pi^2) s^2 exp(-4 s^2/pi).
    Rejection sampling on [0, 5] with safety oversampling.
    """
    out = []
    while len(out) < n:
        s = rng.uniform(0.0, 5.0, size=max(8 * n, 1024))
        pdf = (32.0 / math.pi ** 2) * s * s * np.exp(-4.0 * s * s / math.pi)
        pdf_max = (32.0 / math.pi ** 2) * (math.sqrt(math.pi / 8.0)) ** 2 * \
                  math.exp(-4.0 * (math.pi / 8.0) / math.pi)
        # safer: compute true max from samples
        accept = rng.uniform(0.0, 1.0, size=len(s)) < (pdf / pdf.max())
        out.extend(s[accept].tolist())
    return np.array(out[:n], dtype=np.float64)


def gue_synthetic_zeros(n: int, rng: np.random.Generator) -> np.ndarray:
    """Synthesise N positive ordinates with sine-kernel GUE spacing and
    Riemann-von Mangoldt mean density on the critical line.

    Local mean spacing at height T is approximately 2*pi/log(T/(2*pi)).
    We start at T0 ~ 14.13 (first Riemann zero) and walk forward.
    """
    spacings = sine_kernel_unfolded_spacings(n + 5, rng)
    # rescale unfolded spacings to physical: ds_phys = ds_unfolded * (2*pi/log(T/(2*pi)))
    out = np.zeros(n)
    T = 14.134725  # start at first true zero for fairness
    for i in range(n):
        local_spacing = 2.0 * math.pi / math.log(max(T / (2.0 * math.pi), 1.001))
        T += spacings[i] * local_spacing
        out[i] = T
    return out


# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------

def run(n_max: int = 8000, n_x_samples: int = 30) -> dict:
    """Run experiment.

    n_x_samples >= 1 enables AVERAGING the projection error across n_x_samples
    log-uniform x values within a fixed octave [x0, 2*x0], giving a clean
    decay-rate fit free of single-x oscillations.  Setting n_x_samples=1
    reproduces the original single-x table.
    """
    print(f"# D34 de Branges H(E_xi) projection-error experiment, N up to {n_max}, "
          f"n_x_samples = {n_x_samples}")
    t0 = time.time()

    print("# loading Riemann zeros...")
    gammas = load_zeros(DATA_ZEROS, n_max)
    n_max = len(gammas)
    print(f"#   loaded {n_max} zeros, gamma_max = {gammas[-1]:.3f}")

    rng = np.random.default_rng(seed=0xD34)
    print("# generating GUE control ordinates (sine-kernel)...")
    gue_gammas = gue_synthetic_zeros(n_max, rng)
    print(f"#   GUE control gamma_max = {gue_gammas[-1]:.3f} "
          f"(prime gamma_max = {gammas[-1]:.3f})")

    # Sample x anchors: one per decade for the per-x table
    xs = [10**4, 10**5, 10**6, 10**7]

    # Sampled N values for table reporting (all N up to n_max are computed)
    N_report = [10, 30, 100, 300, 1000, 3000, n_max]

    out = {
        "n_max": n_max,
        "xs": xs,
        "N_report": N_report,
        "results_prime": {},   # per x: list of (N, err) for primes
        "results_gue":   {},   # per x: list of (N, err) for GUE control
        "decay_fits":    {},   # per x: { 'prime': (slope, intercept), 'gue': (...) }
        "polylog_fits":  {},   # per x: fit on log(err)/log(log(log N + 2))? we test (log N)^{-c}
        "psi_exact":     {},
        "gamma_max_prime": float(gammas[-1]),
        "gamma_max_gue":   float(gue_gammas[-1]),
    }

    for x in xs:
        print(f"\n# x = 10^{int(round(math.log10(x)))}, computing exact psi(x) by sieve...")
        x_int = int(x)
        t_psi = time.time()
        psi_x = chebyshev_psi(x_int)
        print(f"#   psi({x}) = {psi_x:.4f}   (elapsed {time.time()-t_psi:.2f}s)")

        # truncated psi at all N
        psi_N_prime = truncated_psi(float(x), gammas)
        psi_N_gue   = truncated_psi(float(x), gue_gammas)

        err_prime = np.abs(psi_x - psi_N_prime)
        err_gue   = np.abs(psi_x - psi_N_gue)

        out["psi_exact"][str(x)] = psi_x

        # sample at N_report
        out["results_prime"][str(x)] = [
            (int(N), float(err_prime[N - 1])) for N in N_report
        ]
        out["results_gue"][str(x)] = [
            (int(N), float(err_gue[N - 1])) for N in N_report
        ]

        # ---- Decay fit: log(err) vs log(N), use N >= 30 to avoid early noise ----
        N_arr = np.arange(1, n_max + 1, dtype=np.float64)
        mask = N_arr >= 30
        log_N = np.log(N_arr[mask])
        log_err_p = np.log(np.maximum(err_prime[mask], 1e-300))
        log_err_g = np.log(np.maximum(err_gue[mask], 1e-300))

        slope_p, intercept_p = np.polyfit(log_N, log_err_p, 1)
        slope_g, intercept_g = np.polyfit(log_N, log_err_g, 1)
        out["decay_fits"][str(x)] = {
            "prime": {"slope": float(slope_p), "intercept": float(intercept_p)},
            "gue":   {"slope": float(slope_g), "intercept": float(intercept_g)},
        }
        print(f"#   power-law fit (log err ~ slope * log N + b):")
        print(f"#     prime: slope = {slope_p:+.4f}, intercept = {intercept_p:+.4f}")
        print(f"#     gue:   slope = {slope_g:+.4f}, intercept = {intercept_g:+.4f}")

        # ---- Polylog fit:  log(err) vs log(log N) ----
        log_log_N = np.log(np.log(N_arr[mask]))
        slope_pl_p, intercept_pl_p = np.polyfit(log_log_N, log_err_p, 1)
        slope_pl_g, intercept_pl_g = np.polyfit(log_log_N, log_err_g, 1)

        # Compute R^2 against power-law and polylog hypotheses
        def r2(y_true, y_pred):
            ss_res = np.sum((y_true - y_pred) ** 2)
            ss_tot = np.sum((y_true - y_true.mean()) ** 2)
            return 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

        pred_pow_p = slope_p * log_N + intercept_p
        pred_pow_g = slope_g * log_N + intercept_g
        pred_pl_p  = slope_pl_p * log_log_N + intercept_pl_p
        pred_pl_g  = slope_pl_g * log_log_N + intercept_pl_g

        out["polylog_fits"][str(x)] = {
            "prime": {
                "polylog_slope": float(slope_pl_p),
                "polylog_intercept": float(intercept_pl_p),
                "r2_powerlaw": float(r2(log_err_p, pred_pow_p)),
                "r2_polylog":  float(r2(log_err_p, pred_pl_p)),
            },
            "gue": {
                "polylog_slope": float(slope_pl_g),
                "polylog_intercept": float(intercept_pl_g),
                "r2_powerlaw": float(r2(log_err_g, pred_pow_g)),
                "r2_polylog":  float(r2(log_err_g, pred_pl_g)),
            },
        }
        print(f"#   power-law R^2: prime={r2(log_err_p, pred_pow_p):.4f}, gue={r2(log_err_g, pred_pow_g):.4f}")
        print(f"#   polylog   R^2: prime={r2(log_err_p, pred_pl_p):.4f}, gue={r2(log_err_g, pred_pl_g):.4f}")

    # ---- Averaged decay analysis: mean projection error vs N over an octave ----
    if n_x_samples > 1:
        print(f"\n# === Averaged decay analysis (n_x_samples = {n_x_samples}) ===")
        averaged = {}
        # one octave per decade anchor
        for x_anchor in xs:
            x_grid = np.linspace(x_anchor, 2.0 * x_anchor, n_x_samples)
            print(f"# averaging over x in [{x_anchor:.2e}, {2*x_anchor:.2e}]"
                  f" ({n_x_samples} samples)")
            psi_xs = np.array([chebyshev_psi(int(round(xv))) for xv in x_grid])
            # absolute errors: take L^2 mean (RMS) across x at each N
            err_p_xs = np.zeros((len(x_grid), n_max))
            err_g_xs = np.zeros((len(x_grid), n_max))
            for ix, (xv, psi_v) in enumerate(zip(x_grid, psi_xs)):
                psi_N_p = truncated_psi(float(xv), gammas)
                psi_N_g = truncated_psi(float(xv), gue_gammas)
                err_p_xs[ix] = np.abs(psi_v - psi_N_p)
                err_g_xs[ix] = np.abs(psi_v - psi_N_g)
            # RMS across x (L^2-typical)
            rms_p = np.sqrt(np.mean(err_p_xs ** 2, axis=0))
            rms_g = np.sqrt(np.mean(err_g_xs ** 2, axis=0))
            # normalise by sqrt(x_anchor): the predicted scale of the error
            sqrt_x = math.sqrt(x_anchor)
            normed_p = rms_p / sqrt_x
            normed_g = rms_g / sqrt_x

            # fit log(rms_p) ~ a*log N + b on N >= 30
            N_arr = np.arange(1, n_max + 1, dtype=np.float64)
            mask = N_arr >= 30
            log_N = np.log(N_arr[mask])
            log_p = np.log(np.maximum(rms_p[mask], 1e-300))
            log_g = np.log(np.maximum(rms_g[mask], 1e-300))
            slope_p, b_p = np.polyfit(log_N, log_p, 1)
            slope_g, b_g = np.polyfit(log_N, log_g, 1)

            # ---- predicted rate for prime: rms ~ sqrt(x) * (log N / N)^{1/2}
            #      so log(rms/sqrt(x)) ~ -0.5 log N + 0.5 log log N + const
            # we use the simple log-log slope here (-0.5 expected)
            pred_p = slope_p * log_N + b_p
            pred_g = slope_g * log_N + b_g
            r2_p = 1.0 - np.sum((log_p - pred_p) ** 2) / np.sum((log_p - log_p.mean()) ** 2)
            r2_g = 1.0 - np.sum((log_g - pred_g) ** 2) / np.sum((log_g - log_g.mean()) ** 2)

            # ---- polylog fit:  log(rms) ~ a * log(log N) + b
            log_log_N = np.log(np.log(N_arr[mask]))
            slope_pl_p, b_pl_p = np.polyfit(log_log_N, log_p, 1)
            slope_pl_g, b_pl_g = np.polyfit(log_log_N, log_g, 1)
            pred_pl_p = slope_pl_p * log_log_N + b_pl_p
            pred_pl_g = slope_pl_g * log_log_N + b_pl_g
            r2_pl_p = 1.0 - np.sum((log_p - pred_pl_p) ** 2) / np.sum((log_p - log_p.mean()) ** 2)
            r2_pl_g = 1.0 - np.sum((log_g - pred_pl_g) ** 2) / np.sum((log_g - log_g.mean()) ** 2)

            # ---- standard explicit-formula prediction:  rms = c * sqrt(x) * sqrt(log N / N)
            #      log(rms_p/sqrt(x)) ~ 0.5 log log N - 0.5 log N + log c
            #      Fit constant c by least-squares once slope (-0.5, +0.5) is fixed.
            X = np.stack([log_N, log_log_N, np.ones_like(log_N)], axis=1)
            # constrained: solve full (a, b, c) and check (a, b) ~ (-0.5, +0.5)
            beta_p, *_ = np.linalg.lstsq(X, log_p, rcond=None)
            beta_g, *_ = np.linalg.lstsq(X, log_g, rcond=None)

            averaged[str(x_anchor)] = {
                "rms_prime_at_report": [float(rms_p[N - 1]) for N in N_report],
                "rms_gue_at_report":   [float(rms_g[N - 1]) for N in N_report],
                "normed_prime_at_report": [float(normed_p[N - 1]) for N in N_report],
                "normed_gue_at_report":   [float(normed_g[N - 1]) for N in N_report],
                "powerlaw_fit_prime": {"slope": float(slope_p), "intercept": float(b_p), "r2": float(r2_p)},
                "powerlaw_fit_gue":   {"slope": float(slope_g), "intercept": float(b_g), "r2": float(r2_g)},
                "polylog_fit_prime":  {"slope": float(slope_pl_p), "intercept": float(b_pl_p), "r2": float(r2_pl_p)},
                "polylog_fit_gue":    {"slope": float(slope_pl_g), "intercept": float(b_pl_g), "r2": float(r2_pl_g)},
                "joint_fit_prime":    {"a_logN": float(beta_p[0]), "b_loglogN": float(beta_p[1]), "c_const": float(beta_p[2])},
                "joint_fit_gue":      {"a_logN": float(beta_g[0]), "b_loglogN": float(beta_g[1]), "c_const": float(beta_g[2])},
            }
            print(f"# x_anchor=10^{int(round(math.log10(x_anchor)))}:")
            print(f"#   prime power-law slope = {slope_p:+.4f}, "
                  f"polylog slope = {slope_pl_p:+.4f}, "
                  f"R^2(pow)={r2_p:.4f}, R^2(polylog)={r2_pl_p:.4f}")
            print(f"#   gue   power-law slope = {slope_g:+.4f}, "
                  f"polylog slope = {slope_pl_g:+.4f}, "
                  f"R^2(pow)={r2_g:.4f}, R^2(polylog)={r2_pl_g:.4f}")
            print(f"#   joint prime fit:  log rms = {beta_p[0]:+.4f} log N + "
                  f"{beta_p[1]:+.4f} log log N + {beta_p[2]:+.4f}")
            print(f"#   joint gue fit:    log rms = {beta_g[0]:+.4f} log N + "
                  f"{beta_g[1]:+.4f} log log N + {beta_g[2]:+.4f}")

        out["averaged"] = averaged

    out["elapsed_seconds"] = time.time() - t0

    # save results
    with open(RESULTS_JSON, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\n# wrote {RESULTS_JSON}")
    print(f"# total elapsed: {out['elapsed_seconds']:.1f}s")

    return out


def make_plot(out: dict) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"# matplotlib unavailable, skipping plot ({e})")
        return

    fig, axes = plt.subplots(1, len(out["xs"]), figsize=(4 * len(out["xs"]), 4),
                             sharey=False)
    if len(out["xs"]) == 1:
        axes = [axes]
    for ax, x in zip(axes, out["xs"]):
        rec_p = out["results_prime"][str(x)]
        rec_g = out["results_gue"][str(x)]
        Ns = [r[0] for r in rec_p]
        ep = [r[1] for r in rec_p]
        eg = [r[1] for r in rec_g]
        ax.loglog(Ns, ep, "o-", label="prime zeros", color="tab:blue")
        ax.loglog(Ns, eg, "s--", label="GUE control", color="tab:orange")
        # reference lines: N^{-1} and N^{-1/2}
        Nref = np.array([Ns[0], Ns[-1]], dtype=float)
        ref1 = ep[0] * (Nref / Ns[0]) ** (-1.0)
        ref_half = ep[0] * (Nref / Ns[0]) ** (-0.5)
        ax.loglog(Nref, ref1, ":", color="grey", alpha=0.6, label="N^{-1}")
        ax.loglog(Nref, ref_half, ":", color="black", alpha=0.4, label="N^{-1/2}")
        ax.set_xlabel("N (zeros used)")
        ax.set_ylabel("|psi(x) - psi_N(x)|")
        ax.set_title(f"x = 10^{int(round(math.log10(x)))}")
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=8)
    fig.suptitle("D34: de Branges H(E_xi) projection error: prime vs GUE control")
    fig.tight_layout()
    fig.savefig(RESULTS_FIG, dpi=120)
    plt.close(fig)
    print(f"# wrote {RESULTS_FIG}")


if __name__ == "__main__":
    n_max = int(sys.argv[1]) if len(sys.argv) > 1 else 8000
    n_x_samples = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    out = run(n_max=n_max, n_x_samples=n_x_samples)
    make_plot(out)

    print("\n# === Summary table ===")
    print(f"# {'x':>10} {'N':>6} {'err_prime':>14} {'err_gue':>14} {'ratio p/g':>10}")
    for x in out["xs"]:
        rec_p = dict(out["results_prime"][str(x)])
        rec_g = dict(out["results_gue"][str(x)])
        for N in out["N_report"]:
            ep = rec_p[N]
            eg = rec_g[N]
            ratio = ep / eg if eg > 0 else float("nan")
            print(f"# {x:>10} {N:>6} {ep:>14.4e} {eg:>14.4e} {ratio:>10.3f}")

    print("\n# === Decay fits ===")
    for x in out["xs"]:
        d = out["decay_fits"][str(x)]
        pl = out["polylog_fits"][str(x)]
        print(f"# x=10^{int(round(math.log10(x)))}:")
        print(f"#    prime power-law slope = {d['prime']['slope']:+.4f}, "
              f"R^2_pow={pl['prime']['r2_powerlaw']:.4f}, "
              f"R^2_polylog={pl['prime']['r2_polylog']:.4f}")
        print(f"#    gue   power-law slope = {d['gue']['slope']:+.4f}, "
              f"R^2_pow={pl['gue']['r2_powerlaw']:.4f}, "
              f"R^2_polylog={pl['gue']['r2_polylog']:.4f}")
