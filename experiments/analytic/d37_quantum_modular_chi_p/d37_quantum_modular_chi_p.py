#!/usr/bin/env python3
"""
D37 — Zagier (2010) quantum-modular cocycle defect of f_N(z) at rationals.

The prime generating polynomial is

    f_N(z) := sum_{n <= N} chi_P(n) z^n

We evaluate phi_N(x) := f_N(e^{2 pi i x}) / pi(N) on Q, then test whether
the cocycle defect under the SL_2(Z) generator S = ((0, -1), (1, 0))

    h_S^k(x) := phi_N(x) - x^{-k} * phi_N(-1/x)

extrapolates to a smooth function on R \\ {0} as N -> infinity, for some
weight k in {0, 1/2, 1, 3/2, 2}.

We also test the lower-triangular generator gamma = ((1, 0), (1, 1)) so
that we exercise a non-trivial denominator change beyond x -> x+1
(the T generator's cocycle is trivially zero by 1-periodicity).

Hardy-Littlewood / Vinogradov prediction at major arc a/q with
gcd(a, q) = 1:

    f_N(e^{2 pi i a/q}) = (mu(q)/phi(q)) pi(N) + O(N exp(-c sqrt log N))

(unconditional) or O(N^{1/2 + epsilon}) on GRH. The leading factor
mu(q)/phi(q) depends only on the denominator q and is multiplicative-
arithmetic. Quantum-modular smoothness requires the cocycle to be
samples of a C^infty function; we test against this.

Falsification criterion (pre-stated):
- (E) Maximum-residual polynomial fit error of h_S^k(a/q)/pi(N) on the
  rational grid is at the HL-imprint scale |mu(q)/phi(q)| (non-trivial
  jumps), not at the noise floor 1/sqrt(pi(N)). Expected outcome:
  cocycle is dominated by HL imprint, NOT C^infty.
- (I) Some weight k* yields polynomial residual at the noise-floor
  level, signalling latent quantum-modular structure.
- (A-grade) The empirical cocycle h_S^k* converges to an explicit
  smooth function H: Q -> C as N grows (relative residual O(N^{-1/2})).

We compute on q in {2, 3, 4, 5, 6, 7, 8, 9, 10, 12} and
N in {2^14, 2^16, 2^18, 2^20}.

Cross-domain reference:
- Zagier 2010 "Quantum modular forms" Clay Math. Proc. 11, 659.
- Bringmann-Folsom-Ono-Rhoades 2017 AMS Coll. 64 ch. 21.
- Hardy-Littlewood 1923 Acta Math. 44 (singular series).

EDGES this composes against / cites:
- E2.21 (S138): |f_N(e^{2pi i a/q})| matches mu(q)^2/phi(q) HL imprint.
- E2.20 (S134): Mahler measure deficit of f_N is a constant -0.307 nat.
- E2.13: Gowers U^k matches HL singular series.

Cites Forman/Joswig/etc. NOT applicable; Zagier framework only.
"""
from __future__ import annotations

import json
import math
import sys
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np


# ---------- prime sieve --------------------------------------------------

def sieve_chi_P(N: int) -> np.ndarray:
    """Return chi_P[0..N], a 0/1 array with chi_P[n] = 1 iff n is prime."""
    chi = np.ones(N + 1, dtype=np.int8)
    chi[0] = 0
    chi[1] = 0
    for p in range(2, int(N ** 0.5) + 1):
        if chi[p]:
            chi[p * p :: p] = 0
    return chi


def mu(q: int) -> int:
    """Möbius mu."""
    if q == 1:
        return 1
    n = q
    result = 1
    p = 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            if n % p == 0:
                return 0
            result = -result
        p += 1
    if n > 1:
        result = -result
    return result


def phi(q: int) -> int:
    """Euler totient."""
    if q == 1:
        return 1
    n = q
    result = q
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return result


def hl_imprint(q: int) -> float:
    """Hardy-Littlewood leading factor mu(q)/phi(q)."""
    return mu(q) / phi(q)


# ---------- evaluation of f_N at e^{2 pi i x} ---------------------------

def eval_fN(chi: np.ndarray, x: float) -> complex:
    """
    Compute f_N(e^{2 pi i x}) = sum_{n <= N} chi_P(n) e^{2 pi i n x}.

    For x rational with denominator <= ~1000, the floating-point
    error is controlled (we accumulate primes only, so the sum has
    pi(N) terms with magnitudes 1).
    """
    n_arr = np.arange(len(chi))
    prime_n = n_arr[chi == 1]  # primes up to N
    return complex(np.sum(np.exp(2j * np.pi * x * prime_n)))


# ---------- SL_2(Z) action on Q -----------------------------------------

def S_action(x: float) -> float:
    """S = ((0,-1),(1,0)) acts as x -> -1/x."""
    return -1.0 / x


def gamma11_action(x: float) -> float:
    """gamma = ((1,0),(1,1)) acts as x -> x/(x+1)."""
    return x / (x + 1.0)


def S_factor(x: float, k: float) -> complex:
    """For S = ((0,-1),(1,0)), cz+d = x; cocycle factor (cz+d)^{-k}."""
    if x == 0:
        return float("inf")
    if x > 0:
        return x ** (-k)
    # negative x: take principal branch
    return (abs(x) ** (-k)) * np.exp(-1j * np.pi * k)


def gamma11_factor(x: float, k: float) -> complex:
    """For gamma = ((1,0),(1,1)), cz+d = x+1; cocycle factor (x+1)^{-k}."""
    if x + 1 == 0:
        return float("inf")
    return (x + 1) ** (-k) if (x + 1) > 0 else (abs(x + 1) ** (-k)) * np.exp(
        -1j * np.pi * k
    )


# ---------- test rationals ----------------------------------------------

def test_rationals(qs: list[int]) -> list[Fraction]:
    """Return all reduced a/q with q in qs, 1 <= a < q, gcd(a,q)=1."""
    out: list[Fraction] = []
    for q in qs:
        for a in range(1, q):
            if math.gcd(a, q) == 1:
                out.append(Fraction(a, q))
    return out


# ---------- cocycle computation -----------------------------------------

@dataclass
class CocycleRow:
    a: int
    q: int
    x: float
    Sx: float            # S(x) = -1/x
    fN_x: complex
    fN_Sx: complex
    pi_N: int
    weight_k: float
    factor_S: complex    # x^{-k} or principal branch
    cocycle_h: complex   # f_N(x) - x^{-k} f_N(Sx), normalised by pi_N
    hl_x: float          # mu(q)/phi(q)
    hl_Sx: float         # mu(q')/phi(q') where q' = denom(S(x) reduced)


def compute_cocycle_gamma11(
    chi: np.ndarray, ratls: list[Fraction], k: float
) -> list[CocycleRow]:
    """
    Compute cocycle defect under gamma = ((1, 0), (1, 1)),
    which acts on x by gamma(x) = x / (x + 1) and has cocycle factor
    (cz + d)^{-k} = (x + 1)^{-k}.
    """
    pi_N = int(chi.sum())
    rows: list[CocycleRow] = []
    cache: dict[tuple[int, int], complex] = {}

    def f_at(frac: Fraction) -> complex:
        a, q = frac.numerator, frac.denominator
        a_mod = a % q
        key = (a_mod, q)
        if key in cache:
            return cache[key]
        v = eval_fN(chi, a_mod / q)
        cache[key] = v
        return v

    for x_frac in ratls:
        a, q = x_frac.numerator, x_frac.denominator
        x = float(x_frac)
        # gamma(x) = x/(x+1) = a/(a+q) -- still in lowest terms since gcd(a, a+q) = gcd(a, q) = 1
        gamma_x_frac = Fraction(a, a + q)
        gamma_x = float(gamma_x_frac)
        q_prime = gamma_x_frac.denominator
        fN_x_val = f_at(x_frac)
        fN_gx_val = f_at(gamma_x_frac)
        factor = gamma11_factor(x, k)
        h = (fN_x_val - factor * fN_gx_val) / pi_N
        rows.append(
            CocycleRow(
                a=a,
                q=q,
                x=x,
                Sx=gamma_x,
                fN_x=fN_x_val,
                fN_Sx=fN_gx_val,
                pi_N=pi_N,
                weight_k=k,
                factor_S=factor,
                cocycle_h=h,
                hl_x=hl_imprint(q),
                hl_Sx=hl_imprint(q_prime),
            )
        )
    return rows


def hl_predicted_cocycle_general(rows: list[CocycleRow], factor_fn) -> dict:
    """HL match for any generator with cocycle factor factor_fn(x, k)."""
    h_emp_re = np.array([r.cocycle_h.real for r in rows])
    h_emp_im = np.array([r.cocycle_h.imag for r in rows])
    k = rows[0].weight_k
    h_pred_re, h_pred_im = [], []
    for r in rows:
        f = factor_fn(r.x, k)
        pred = r.hl_x - f * r.hl_Sx
        h_pred_re.append(pred.real if isinstance(pred, complex) else pred)
        h_pred_im.append(pred.imag if isinstance(pred, complex) else 0.0)
    h_pred_re = np.array(h_pred_re)
    h_pred_im = np.array(h_pred_im)
    rss_re = float(np.sum((h_emp_re - h_pred_re) ** 2))
    rss_im = float(np.sum((h_emp_im - h_pred_im) ** 2))
    tss_re = float(np.sum((h_emp_re - h_emp_re.mean()) ** 2))
    tss_im = float(np.sum((h_emp_im - h_emp_im.mean()) ** 2))
    return {
        "rss_re_vs_HL": rss_re,
        "rss_im_vs_HL": rss_im,
        "rel_residual_re_vs_HL": math.sqrt(rss_re / tss_re) if tss_re > 0 else 0.0,
        "rel_residual_im_vs_HL": math.sqrt(rss_im / tss_im) if tss_im > 0 else 0.0,
        "max_abs_emp_re": float(np.max(np.abs(h_emp_re))),
        "max_abs_emp_im": float(np.max(np.abs(h_emp_im))),
    }


def compute_cocycle(
    chi: np.ndarray, ratls: list[Fraction], k: float
) -> list[CocycleRow]:
    pi_N = int(chi.sum())
    rows: list[CocycleRow] = []
    # cache evaluations of f_N at every rational point we need
    cache: dict[tuple[int, int], complex] = {}

    def f_at(frac: Fraction) -> complex:
        # Reduce to [0, 1) for cache key
        a, q = frac.numerator, frac.denominator
        a_mod = a % q
        key = (a_mod, q)
        if key in cache:
            return cache[key]
        v = eval_fN(chi, a_mod / q)
        cache[key] = v
        return v

    for x_frac in ratls:
        a, q = x_frac.numerator, x_frac.denominator
        x = float(x_frac)
        # S action: -1/x as a Fraction
        Sx_frac = Fraction(-q, a)  # = -q/a
        Sx_reduced = Sx_frac - Sx_frac.__floor__()  # in [0, 1)
        Sx = float(Sx_reduced)
        q_prime = Sx_reduced.denominator
        fN_x_val = f_at(x_frac)
        fN_Sx_val = f_at(Sx_reduced)
        factor = S_factor(x, k)
        h = (fN_x_val - factor * fN_Sx_val) / pi_N
        rows.append(
            CocycleRow(
                a=a,
                q=q,
                x=x,
                Sx=Sx,
                fN_x=fN_x_val,
                fN_Sx=fN_Sx_val,
                pi_N=pi_N,
                weight_k=k,
                factor_S=factor,
                cocycle_h=h,
                hl_x=hl_imprint(q),
                hl_Sx=hl_imprint(q_prime),
            )
        )
    return rows


# ---------- smoothness test --------------------------------------------

def smoothness_residual(rows: list[CocycleRow], deg: int = 4) -> dict:
    """
    Test whether real / imag part of h_S(a/q) is well-fit by a degree-d
    polynomial in x = a/q.  The residual is the L^2 norm divided by L^2
    norm of the data; small residual = candidate smooth function.
    """
    xs = np.array([r.x for r in rows], dtype=float)
    h_re = np.array([r.cocycle_h.real for r in rows])
    h_im = np.array([r.cocycle_h.imag for r in rows])
    # use 'unique' x ordering for fit
    order = np.argsort(xs)
    xs_s = xs[order]
    h_re_s = h_re[order]
    h_im_s = h_im[order]
    if len(xs_s) <= deg + 1:
        return {"deg": deg, "residual_re": None, "residual_im": None}
    p_re = np.polyfit(xs_s, h_re_s, deg)
    p_im = np.polyfit(xs_s, h_im_s, deg)
    pred_re = np.polyval(p_re, xs_s)
    pred_im = np.polyval(p_im, xs_s)
    rss_re = float(np.sum((h_re_s - pred_re) ** 2))
    rss_im = float(np.sum((h_im_s - pred_im) ** 2))
    tss_re = float(np.sum((h_re_s - h_re_s.mean()) ** 2))
    tss_im = float(np.sum((h_im_s - h_im_s.mean()) ** 2))
    rel_re = math.sqrt(rss_re / tss_re) if tss_re > 0 else 0.0
    rel_im = math.sqrt(rss_im / tss_im) if tss_im > 0 else 0.0
    return {
        "deg": deg,
        "rss_re": rss_re,
        "rss_im": rss_im,
        "tss_re": tss_re,
        "tss_im": tss_im,
        "rel_residual_re": rel_re,
        "rel_residual_im": rel_im,
        "max_abs_re": float(np.max(np.abs(h_re_s - pred_re))),
        "max_abs_im": float(np.max(np.abs(h_im_s - pred_im))),
    }


def hl_predicted_cocycle(rows: list[CocycleRow]) -> dict:
    """
    Compare empirical h_S^k(a/q) to the HL leading prediction:
      h_HL(a/q) := mu(q)/phi(q) - (a/q)^{-k} * mu(q')/phi(q')
    Reports L^2 closeness.
    """
    h_emp_re = np.array([r.cocycle_h.real for r in rows])
    h_emp_im = np.array([r.cocycle_h.imag for r in rows])
    k = rows[0].weight_k
    h_pred_re = []
    h_pred_im = []
    for r in rows:
        hl_x = r.hl_x
        hl_Sx = r.hl_Sx
        f = S_factor(r.x, k)
        pred = hl_x - f * hl_Sx
        h_pred_re.append(pred.real if isinstance(pred, complex) else pred)
        h_pred_im.append(pred.imag if isinstance(pred, complex) else 0.0)
    h_pred_re = np.array(h_pred_re)
    h_pred_im = np.array(h_pred_im)
    rss_re = float(np.sum((h_emp_re - h_pred_re) ** 2))
    rss_im = float(np.sum((h_emp_im - h_pred_im) ** 2))
    tss_re = float(np.sum((h_emp_re - h_emp_re.mean()) ** 2))
    tss_im = float(np.sum((h_emp_im - h_emp_im.mean()) ** 2))
    return {
        "rss_re_vs_HL": rss_re,
        "rss_im_vs_HL": rss_im,
        "tss_re": tss_re,
        "tss_im": tss_im,
        "rel_residual_re_vs_HL": math.sqrt(rss_re / tss_re) if tss_re > 0 else 0.0,
        "rel_residual_im_vs_HL": math.sqrt(rss_im / tss_im) if tss_im > 0 else 0.0,
        "max_abs_emp_re": float(np.max(np.abs(h_emp_re))),
        "max_abs_emp_im": float(np.max(np.abs(h_emp_im))),
        "max_abs_pred_re": float(np.max(np.abs(h_pred_re))),
        "max_abs_pred_im": float(np.max(np.abs(h_pred_im))),
    }


# ---------- bernoulli null ----------------------------------------------

def bernoulli_null(N: int, density: float, seed: int = 42) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return (rng.random(N + 1) < density).astype(np.int8)


# ---------- main --------------------------------------------------------

def main(N: int, qmax: int, weights: list[float], out_dir: Path) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"# D37 — quantum modular cocycle of chi_P, N = {N:,}, qmax = {qmax}")
    t0 = time.time()
    chi = sieve_chi_P(N)
    pi_N = int(chi.sum())
    print(f"  pi(N) = {pi_N}, sieve time {time.time() - t0:.1f}s")

    qs = list(range(2, qmax + 1))
    ratls = test_rationals(qs)
    print(f"  testing {len(ratls)} reduced rationals a/q with q in [2, {qmax}]")

    summary: dict = {
        "N": N,
        "pi_N": pi_N,
        "qmax": qmax,
        "n_rationals": len(ratls),
        "by_weight": {},
    }

    for k in weights:
        t0 = time.time()
        rows = compute_cocycle(chi, ratls, k)
        smooth = smoothness_residual(rows, deg=4)
        smooth6 = smoothness_residual(rows, deg=6)
        hl_match = hl_predicted_cocycle(rows)
        summary["by_weight"][f"k={k}"] = {
            "smooth_deg4": smooth,
            "smooth_deg6": smooth6,
            "hl_match": hl_match,
            "n_rows": len(rows),
            "compute_time_s": time.time() - t0,
        }
        print(
            f"  k={k}: deg-4 rel-resid (re,im) = ("
            f"{smooth['rel_residual_re']:.3f},"
            f"{smooth['rel_residual_im']:.3f}), "
            f"vs HL rel-resid (re,im) = ("
            f"{hl_match['rel_residual_re_vs_HL']:.3f},"
            f"{hl_match['rel_residual_im_vs_HL']:.3f}), "
            f"|h|_max (re,im) = ("
            f"{hl_match['max_abs_emp_re']:.3f},"
            f"{hl_match['max_abs_emp_im']:.3f})"
        )

        # raw rows dump
        rows_dump = [
            {
                "a": r.a,
                "q": r.q,
                "x": r.x,
                "Sx": r.Sx,
                "h_re": r.cocycle_h.real,
                "h_im": r.cocycle_h.imag,
                "hl_x": r.hl_x,
                "hl_Sx": r.hl_Sx,
            }
            for r in rows
        ]
        with (out_dir / f"rows_k{k}_N{N}.json").open("w") as fh:
            json.dump(rows_dump, fh, indent=2)

    # Test second generator gamma11 = ((1,0),(1,1)) at k=1 to verify
    # HL imprint locks the cocycle for non-S generators too.
    print(f"  Second generator gamma = ((1,0),(1,1)) at k=1.0:")
    rows_g = compute_cocycle_gamma11(chi, ratls, k=1.0)
    smooth_g = smoothness_residual(rows_g, deg=4)
    hl_g = hl_predicted_cocycle_general(rows_g, factor_fn=gamma11_factor)
    summary["gamma11_k1"] = {"smooth_deg4": smooth_g, "hl_match": hl_g}
    print(
        f"    deg-4 rel-resid (re,im) = ("
        f"{smooth_g['rel_residual_re']:.3f},"
        f"{smooth_g['rel_residual_im']:.3f}), "
        f"vs HL rel-resid (re,im) = ("
        f"{hl_g['rel_residual_re_vs_HL']:.3f},"
        f"{hl_g['rel_residual_im_vs_HL']:.3f})"
    )

    # Also run Bernoulli null at all weights for parity:
    print(f"  Bernoulli null:")
    chi_null = bernoulli_null(N, density=pi_N / N, seed=42)
    summary["bernoulli_null"] = {}
    for k in weights:
        rows_null = compute_cocycle(chi_null, ratls, k=k)
        smooth_null = smoothness_residual(rows_null, deg=4)
        hl_null = hl_predicted_cocycle(rows_null)
        summary["bernoulli_null"][f"k={k}"] = {
            "smooth_deg4": smooth_null,
            "hl_match": hl_null,
        }
        print(
            f"    k={k}: deg-4 rel-resid (re,im) = ("
            f"{smooth_null['rel_residual_re']:.3f},"
            f"{smooth_null['rel_residual_im']:.3f}), "
            f"vs HL rel-resid (re,im) = ("
            f"{hl_null['rel_residual_re_vs_HL']:.3f},"
            f"{hl_null['rel_residual_im_vs_HL']:.3f}), "
            f"|h|_max (re,im) = ("
            f"{hl_null['max_abs_emp_re']:.3f},"
            f"{hl_null['max_abs_emp_im']:.3f})"
        )

    with (out_dir / f"summary_N{N}.json").open("w") as fh:
        json.dump(summary, fh, indent=2)
    return summary


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, default=2 ** 16)
    p.add_argument("--qmax", type=int, default=12)
    p.add_argument("--out", type=str, default="results")
    p.add_argument(
        "--weights",
        type=str,
        default="0,0.5,1,1.5,2",
        help="comma-separated weight list",
    )
    args = p.parse_args()
    weights = [float(w) for w in args.weights.split(",")]
    out_dir = Path(__file__).parent / args.out
    main(N=args.N, qmax=args.qmax, weights=weights, out_dir=out_dir)
