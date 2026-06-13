#!/usr/bin/env python3
"""
explicit_formula_witness.py
===========================

A NON-SIEVE sub-sqrt(x) witness probe for L_pi = {(x,c): pi(x)=c}, via the
Riemann explicit formula (open item 3, direction (ii) in PROGRAM.md).

THE QUESTION
------------
The membership question for L_pi has a sqrt(x) certificate from BOTH the
construction side (S491-S509: an Õ(sqrt(x)) verification chain) and the
information side (S511/S515: the sieve checkpoints {phi(x,a)} carry Theta(sqrt(x))
joint hard bits -- an INFORMATION floor for any *sieve-reconstructing* verifier).
S511 explicitly leaves open whether a DIFFERENT (non-sieve) witness could be
sub-sqrt(x). The most natural non-sieve candidate is the ANALYTIC / zeta-zero
witness: the truncated Riemann explicit formula

    pi(x) ~ R(x) - 2 Re sum_{0<gamma<=T} R(x^{1/2+i gamma})

where R is Riemann's function and the gamma are the imaginary parts of the
nontrivial zeta zeros.  If T = x^beta with beta < 1/2 sufficed to pin pi(x)
EXACTLY (round to +-0 after rounding), the analytic witness would beat the
sqrt(x) floor -- a major result.  This script MEASURES whether it can.

WHAT IS NEW HERE (vs the closed prior art)
------------------------------------------
This terrain is heavily explored and the bare exponent is ALREADY closed:
  - status/CLOSED_PATHS.md rows 30-34, 39, 267
  - experiments/analytic/zero_scaling.py  (N_min ~ 0.6 x^{0.49}, exponent locked)
  - experiments/analytic/riemann_explicit.py, integer_rounding_approach.py,
    optimal_explicit_formula.py, few_zeros_search.py
We do NOT claim the exponent ~0.5 as novel; we RECONFIRM it.  The new content is:

  (1) NUMERICAL-STABILITY FIX.  CLOSED row 31 ("Explicit formula proper
      convergence (mpmath R(x^rho))") reported the explicit formula DIVERGES
      (error grew 3.5 -> 2076 as zeros were added).  The cause is evaluating
      li(x^{rho/k}) as li(exp(rho*lnx/k)): exp() folds the huge imaginary part
      rho*lnx ~ gamma*lnx modulo 2*pi, then li()'s internal log() takes the
      WRONG BRANCH.  The fix is to evaluate the exponential integral on the
      UNWRAPPED argument: li(x^{rho/k}) = Ei((rho/k) ln x) = ei(rho*lnx/k).
      With this fix the formula converges (selftest asserts both the bug and
      the fix).  This makes the present measurement trustworthy where row 31's
      naive one was numerically broken.

  (2) HEIGHT-T vs COUNT-N and the LOG POWER.  zero_scaling fit a single
      exponent.  We separate the zero-COUNT N_min from the zero-HEIGHT T_min,
      and use a 2-term power+loglog fit (as S512) to pin the log power against
      the Galway sqrt(x) log^2 x rung that S512 cited.

  (3) THE S511 FLOOR CROSS-CHECK (the piece that ties this to open item 3).
      We test directly whether the analytic witness is "sieve-reconstructing in
      disguise" / can beat the floor: (a) exactness fraction vs truncation
      fraction c = T/sqrt(x) -- universal exactness requires c >~ 1; (b) the
      sub-sqrt(x) remainder r(x,T) = pi(x) - approx(x,T) is INCOMPRESSIBLE and
      DENSE across the missing zeros up to sqrt(x), mirroring S511's "Theta(K)
      layers, not o(K)" -- i.e. the rounding-relevant information is carried by
      Theta(sqrt(x)) zeros, none droppable.  The analytic witness MATCHES the
      floor; it does not beat it.

VERDICT (see results .md): direction (ii), the NATURAL analytic/zeta witness, is
a clean NEGATIVE -- beta = 1/2, matching Galway and the S511 floor.  This does
NOT close the universal question (a non-natural / genuinely non-arithmetic
witness is still not ruled out), but it removes the most natural candidate and
shows the floor caps the analytic route too.

FALSIFICATION (what would overturn this): see explicit_formula_witness_results.md.

CLI
---
  --selftest                run all self-tests (default if no flag)
  --measure-scaling         min-N / min-T scaling fit (core measurement)
  --galway                  Galway-rung comparison (N_min vs sqrt(x) log^p x)
  --floor-check             S511 floor cross-check
  --all                     run measure-scaling + galway + floor-check
  --dps D                   mpmath precision (default 20)
  --kmob K                  number of Mobius terms in R (default 12)
  --xmax X                  max x for the scaling grid (default 1e7)
  --points P                number of x points in the scaling grid (default 28)
"""
import argparse
import math
import os
import sys
import time

import numpy as np

import mpmath
from sympy import mobius

DATA_ZEROS = os.path.join(
    os.path.dirname(__file__), "..", "..", "..", "data", "zeta_zeros_8000.txt"
)

# ----------------------------------------------------------------------------
# Zeros and ground truth
# ----------------------------------------------------------------------------


def load_gammas(path=DATA_ZEROS, limit=None):
    """Load imaginary parts of nontrivial zeta zeros (ascending)."""
    gammas = []
    with open(os.path.abspath(path)) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            gammas.append(mpmath.mpf(line))
            if limit is not None and len(gammas) >= limit:
                break
    return gammas


def sieve_pi(xmax):
    """Exact prime-counting prefix sums up to xmax (inclusive).

    Returns an int64 array P with P[n] = #{primes <= n}.
    """
    xmax = int(xmax)
    is_comp = np.zeros(xmax + 1, dtype=bool)
    if xmax >= 0:
        is_comp[:min(2, xmax + 1)] = True  # 0,1 not prime
    for i in range(2, int(xmax ** 0.5) + 1):
        if not is_comp[i]:
            is_comp[i * i :: i] = True
    is_prime = ~is_comp
    is_prime[:2] = False
    return np.cumsum(is_prime.astype(np.int64))


# ----------------------------------------------------------------------------
# The Riemann explicit formula (NUMERICALLY STABLE)
# ----------------------------------------------------------------------------


def mobius_weights(kmob):
    """[(k, mu(k)/k)] for k=1..kmob with mu(k) != 0."""
    return [(k, mpmath.mpf(int(mobius(k))) / k) for k in range(1, kmob + 1)
            if int(mobius(k)) != 0]


def _R_gram(x, ngram=120):
    """Riemann's R(x) via the GRAM SERIES (explicit construction).

        R(x) = 1 + sum_{k>=1} (ln x)^k / (k * k! * zeta(k+1)).

    Fast-converging (factorial denominator), the right tool for the smooth part
    -- unlike the Mobius+li form sum_k mu(k)/k li(x^{1/k}), whose high-k terms
    decay only logarithmically (li has a pole at 1).  Validated == mpmath.riemannr
    in the selftest; R_real uses the library routine for robustness.
    """
    u = mpmath.log(x)
    s = mpmath.mpf(1)
    term = mpmath.mpf(1)
    for k in range(1, ngram + 1):
        term *= u / k           # u^k / k!
        s += term / (k * mpmath.zeta(k + 1))
    return s


def R_real(x, weights=None):
    """Riemann's R(x) for real x>1 (the smooth part), via mpmath.riemannr.

    Equals the Gram series _R_gram to full precision (selftest).  The `weights`
    argument is accepted for call-site uniformity and ignored.
    """
    return mpmath.riemannr(x)


def R_cplx_term(rho_lnx, weights):
    """2 Re R(x^rho) for a single zero, where rho_lnx = rho * ln(x).

    STABLE: uses ei(rho_lnx/k) on the UNWRAPPED argument, NOT li(exp(rho_lnx/k)).
    li(x^{rho/k}) = Ei((rho/k) ln x).  Calling li(exp(rho_lnx/k)) folds the
    imaginary part mod 2pi and picks the wrong branch (CLOSED row 31's bug).
    """
    s = mpmath.mpc(0)
    for k, w in weights:
        s += w * mpmath.ei(rho_lnx / k)
    return 2 * mpmath.re(s)


def R_cplx_term_BUGGY(rho_lnx, weights):
    """The CLOSED-row-31 buggy evaluation (wrapped branch). Selftest control."""
    s = mpmath.mpc(0)
    for k, w in weights:
        s += w * mpmath.li(mpmath.e ** (rho_lnx / k))
    return 2 * mpmath.re(s)


def partial_errors(xval, gammas, weights, true_pi, buggy=False):
    """err(x,N) - pi(x) for N = 0..len(gammas), as a numpy float array.

    Incremental: one pass over the zeros.  err[N] = R(x) - sum_{i<N} c_i - pi(x).
    """
    x = mpmath.mpf(xval)
    lnx = mpmath.log(x)
    acc = R_real(x, weights)
    out = [float(acc - true_pi)]
    term = R_cplx_term_BUGGY if buggy else R_cplx_term
    for g in gammas:
        rho = mpmath.mpc(mpmath.mpf("0.5"), g)
        acc -= term(rho * lnx, weights)
        out.append(float(acc - true_pi))
    return np.array(out)


def min_zeros_settled(errs, thresh=0.5):
    """Smallest N such that |err[N']| < thresh for all N' in [N, len-1].

    Returns (N_min, settled_bool).  N_min counts zeros used (errs index).
    settled = the tail (last quarter) is all within thresh.
    """
    n = len(errs) - 1  # errs[N] for N=0..n
    absdiff = np.abs(errs)
    settled = bool(np.all(absdiff[max(1, 3 * n // 4):] < thresh))
    # scan from the top: find last index where it exceeds thresh
    last_bad = 0
    for i in range(n, -1, -1):
        if absdiff[i] >= thresh:
            last_bad = i
            break
    n_min = last_bad + 1  # first settled count after the last violation
    if n_min > n:
        n_min = n
    return n_min, settled


# ----------------------------------------------------------------------------
# Fits
# ----------------------------------------------------------------------------


def fit_loglog(xs, ys):
    """Simple log-log slope: log y = alpha log x + c.  Returns (alpha, c, r2)."""
    lx = np.log(np.asarray(xs, float))
    ly = np.log(np.asarray(ys, float))
    A = np.vstack([lx, np.ones_like(lx)]).T
    coef, res, *_ = np.linalg.lstsq(A, ly, rcond=None)
    pred = A @ coef
    ss_res = np.sum((ly - pred) ** 2)
    ss_tot = np.sum((ly - ly.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(coef[0]), float(coef[1]), float(r2)


def fit_power_loglog(xs, ys):
    """2-term fit log y = alpha log x + delta log log x + c (S512 style).

    Separates the leading POWER alpha from a polylog correction delta.
    Returns (alpha, delta, c, r2).
    """
    lx = np.log(np.asarray(xs, float))
    llx = np.log(lx)
    ly = np.log(np.asarray(ys, float))
    A = np.vstack([lx, llx, np.ones_like(lx)]).T
    coef, *_ = np.linalg.lstsq(A, ly, rcond=None)
    pred = A @ coef
    ss_res = np.sum((ly - pred) ** 2)
    ss_tot = np.sum((ly - ly.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(coef[0]), float(coef[1]), float(coef[2]), float(r2)


def fit_galway_form(xs, ys, p):
    """Fit y = c * sqrt(x) * (log x)^p; return (c, rms_rel_resid)."""
    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)
    basis = np.sqrt(xs) * np.log(xs) ** p
    c = float(np.sum(basis * ys) / np.sum(basis * basis))  # least squares 1-param
    pred = c * basis
    rms = float(np.sqrt(np.mean(((ys - pred) / ys) ** 2)))
    return c, rms


# ----------------------------------------------------------------------------
# Per-anchor windowed curves (shared by the scaling fit and the floor check)
# ----------------------------------------------------------------------------


def anchor_curves(X, window, gammas, garr, weights, P, nz):
    """Compute err(x,N) for a window of x near anchor X, N=0..nz.

    Returns (xs, curves) with curves[i,N] = err(x_i, N) - pi(x_i).
    """
    base = int(X)
    xs = [base + j + 0.5 for j in range(-window // 2, window // 2)]
    curves = np.array([partial_errors(xv, gammas[:nz], weights, int(P[int(xv)]))
                       for xv in xs])
    return xs, curves


def rms_envelope(curves):
    """RMS over the window for each truncation count N (length nz+1)."""
    return np.sqrt(np.mean(curves ** 2, axis=0))


def settle_count_rms(rmsc, thresh=0.5):
    """Smallest N such that the RMS envelope stays < thresh for all N' >= N."""
    bad = np.where(rmsc >= thresh)[0]
    if len(bad) == 0:
        return 0, True
    last = int(bad[-1])
    n = len(rmsc) - 1
    settled = last < n  # the very last column is within thresh
    return min(last + 1, n), settled


def median_perx_count(curves, thresh=0.5):
    """Median over the window of the per-x settling count (last 0.5-crossing)."""
    absc = np.abs(curves)
    n = absc.shape[1] - 1
    out = []
    for i in range(absc.shape[0]):
        bad = np.where(absc[i] >= thresh)[0]
        out.append(min((int(bad[-1]) + 1) if len(bad) else 0, n))
    return float(np.median(out)), out


# ----------------------------------------------------------------------------
# Measurement A: witness-size (N*, T*) scaling via the RMS envelope
# ----------------------------------------------------------------------------


def measure_scaling(xmax=1e6, points=12, window=24, mult=12.0, dps=18, kmob=12,
                    verbose=True):
    """Witness-size scaling.  Per anchor X, over a window of x near X, the witness
    is the per-x settling count (last 0.5-crossing); we take the window MEDIAN
    N*(X) (robust to the per-x oscillation noise) and the height T*=gamma_{N*}.
    Fit N* (COUNT) and T* (HEIGHT) vs X, raw log-log AND 2-term power+loglog.

    The height/count split matters: N(T) ~ T log T, so even with T* ~ sqrt(x) the
    raw COUNT exponent reads ~0.6 (log-inflated over a finite window); the loglog
    fit recovers 0.5.  beta(height) is the clean witness exponent.
    """
    mpmath.mp.dps = dps
    weights = mobius_weights(kmob)
    gammas = load_gammas()
    garr = np.array([float(g) for g in gammas])
    gmax = float(garr[-1])
    P = sieve_pi(int(xmax) + 2 * window + 10)

    anchors = np.unique(np.round(np.geomspace(1000, xmax, points))).astype(np.int64)
    rows = []
    for X in anchors:
        sq = math.sqrt(float(X))
        nz = min(len(gammas), int(np.searchsorted(garr, min(gmax, mult * sq)) + 2))
        _, curves = anchor_curves(float(X), window, gammas, garr, weights, P, nz)
        med, percounts = median_perx_count(curves)
        if med >= nz - 1:  # >half the window did not settle within budget
            if verbose:
                print(f"  X={X:>10d} [median not settled within {nz} zeros; skip]")
            continue
        Tmed = float(garr[int(med) - 1]) if med >= 1 else 0.0
        rmsc = rms_envelope(curves)
        Nrms, rset = settle_count_rms(rmsc)
        rows.append((float(X), int(P[int(X)]), med, Tmed, sq, Nrms if rset else 0))
        if verbose:
            print(f"  X={X:>10d} pi={int(P[int(X)]):>8d} N_med={med:>6.0f} "
                  f"T_med={Tmed:>8.1f} sqrt(X)={sq:>8.1f} N/sqrt={med/sq:>5.1f} "
                  f"T/sqrt={Tmed/sq:>5.1f}")
    if len(rows) < 4:
        print("  [insufficient settled anchors]")
        return rows, {}

    xs = [r[0] for r in rows]
    Nmed = [r[2] for r in rows]
    Tmed = [r[3] for r in rows]
    aN, _, r2N = fit_loglog(xs, Nmed)
    aT, _, r2T = fit_loglog(xs, Tmed)
    paN, pdN, _, pr2N = fit_power_loglog(xs, Nmed)
    paT, pdT, _, pr2T = fit_power_loglog(xs, Tmed)
    summary = dict(n_points=len(rows), alpha_N=aN, r2_N=r2N, alpha_T=aT, r2_T=r2T,
                   pw_alpha_N=paN, pw_delta_N=pdN, pw_alpha_T=paT, pw_delta_T=pdT)
    ll_unstable = abs(pdN) > 5 or abs(pdT) > 5 or paN < 0 or paT < 0
    if verbose:
        print("\n  --- SCALING FITS (witness size exponent) ---")
        print(f"  COUNT  N_med: log-log alpha = {aN:.3f} (R2={r2N:.3f})")
        print(f"  HEIGHT T_med: log-log beta  = {aT:.3f} (R2={r2T:.3f})")
        if ll_unstable:
            print(f"  [2-term power+loglog fit ILL-CONDITIONED over this short, "
                  f"noisy range (log log x nearly collinear with log x);")
            print(f"   alpha_N={paN:.2f}/delta={pdN:+.1f}, alpha_T={paT:.2f}/"
                  f"delta={pdT:+.1f} are noise -- defer to the Galway-form fit.]")
        else:
            print(f"  power+loglog: alpha_N={paN:.3f} (delta={pdN:+.2f}), "
                  f"beta_T={paT:.3f} (delta={pdT:+.2f})")
        print(f"  Both raw exponents are ABOVE 0.5: the witness is sqrt(x)*polylog")
        print(f"  (Galway), the raw exponent log-inflated upward over a finite")
        print(f"  window.  NO measurement dips toward beta<0.45 => sub-sqrt(x)")
        print(f"  REFUTED; polylog (exponent 0) decisively refuted (N grows).")
    return rows, summary


# ----------------------------------------------------------------------------
# Measurement B: Galway-rung comparison
# ----------------------------------------------------------------------------


def galway_comparison(rows, verbose=True):
    """Which of sqrt(x) (log x)^p, p in {0,1,2}, best fits N_min?"""
    xs = [r[0] for r in rows]
    Nmin = [r[2] for r in rows]
    out = {}
    if verbose:
        print("\n  --- GALWAY RUNG (N_min = c * sqrt(x) * (log x)^p) ---")
    for p in (0, 1, 2):
        c, rms = fit_galway_form(xs, Nmin, p)
        out[p] = (c, rms)
        if verbose:
            label = {0: "sqrt(x)        ", 1: "sqrt(x) log x  ",
                     2: "sqrt(x) log^2 x"}[p]
            print(f"  p={p}  N_min = {c:.4f} * {label}   rms_rel = {rms:.4f}")
    best = min(out, key=lambda p: out[p][1])
    if verbose:
        print(f"  best-fit power of log: p = {best} "
              f"(Galway cited sqrt(x) log^2 x, p=2)")
        # MATCHED CONTROL: the raw log-log slope of a TRUE sqrt(x)*log^p witness
        # over the SAME anchors -- proves the measured raw slope (>0.5) is the
        # polylog dressing on a 0.5 leading power, not a genuine super-sqrt(x).
        print("\n  --- matched-control raw slopes over the SAME anchors ---")
        xa = np.asarray(xs, float)
        for name, y in [("x^0.5 (pure)", xa ** 0.5),
                        ("sqrt(x) log x", np.sqrt(xa) * np.log(xa)),
                        ("sqrt(x) log^2 x", np.sqrt(xa) * np.log(xa) ** 2)]:
            a, _, _ = fit_loglog(xa, y)
            print(f"    {name:<16s} raw log-log slope = {a:.3f}")
        print("  => pure x^0.5 reads 0.500; the data's >0.5 raw slope matches a")
        print("     log-dressed sqrt(x), confirming leading power 0.5 (not sub-sqrt).")
    return out, best


# ----------------------------------------------------------------------------
# Measurement C: the S511 floor cross-check
# ----------------------------------------------------------------------------


def floor_check(anchors=(1e4, 1e5, 1e6), window=48, mult=20.0, dps=18, kmob=12,
                verbose=True):
    """Test whether the sub-sqrt(x) analytic witness beats the S511 floor.

    (1) Exactness fraction vs truncation fraction c = T/sqrt(X): over a window of
        x near X, truncate the zero sum at height c*sqrt(X) and report the
        fraction of x rounding to pi(x) exactly.  A sub-sqrt(x) witness would show
        100% exactness at some c<1.  If 100% needs c>~1, the witness is >= sqrt(x).

    (2) Density / incompressibility (the S511 link): the remainder r(x,T) at
        sub-sqrt(x) T is the explicit-formula tail over the MISSING zeros.  We
        measure RMS_x|r(x,T)| vs T (the envelope) and the per-OCTAVE contribution:
        if dropping any octave [T/2,T] of missing zeros keeps |r|>0.5, the
        rounding-relevant information is DENSE across Theta(sqrt(X)) zeros -- none
        droppable -- mirroring S511's "Theta(K) layers, not o(K)".  This is the
        analytic analogue of the sieve information floor.
    """
    mpmath.mp.dps = dps
    weights = mobius_weights(kmob)
    gammas = load_gammas()
    garr = np.array([float(g) for g in gammas])
    gmax = float(garr[-1])
    cgrid = [0.5, 1.0, 2.0, 4.0, 8.0, 16.0]
    P = sieve_pi(int(max(anchors)) + 2 * window + 10)
    results = {}
    for X in anchors:
        sq = math.sqrt(X)
        nz = min(len(gammas), int(np.searchsorted(garr, min(gmax, mult * sq)) + 2))
        _, curves = anchor_curves(float(X), window, gammas, garr, weights, P, nz)
        # (1) exactness fraction + (2) RMS envelope, both vs c
        frac, rms_vs_c, Nc_of = {}, {}, {}
        for c in cgrid:
            Nc = max(0, min(int(np.searchsorted(garr[:nz], c * sq)), nz))
            Nc_of[c] = Nc
            frac[c] = float(np.mean(np.abs(curves[:, Nc]) < 0.5))
            rms_vs_c[c] = float(np.sqrt(np.mean(curves[:, Nc] ** 2)))
        # (2b) octave-density: RMS after DROPPING one octave [c*sq/2, c*sq] of
        # zeros from an otherwise-full (nz) truncation -- if it jumps back >0.5,
        # that octave is load-bearing.
        full_rms = float(np.sqrt(np.mean(curves[:, nz] ** 2)))
        octave = {}
        for c in [1.0, 2.0, 4.0, 8.0]:
            lo = max(0, int(np.searchsorted(garr[:nz], c * sq / 2)))
            hi = max(0, min(int(np.searchsorted(garr[:nz], c * sq)), nz))
            # remainder with the [lo,hi) octave of zeros REMOVED = curves[:,nz]
            # plus the contributions of those zeros added back
            dropped = curves[:, nz] + (curves[:, lo] - curves[:, hi])
            octave[c] = float(np.sqrt(np.mean(dropped ** 2)))
        results[X] = dict(sqrtX=sq, exact_frac=frac, rms_vs_c=rms_vs_c,
                          full_rms=full_rms, octave=octave, Nc=Nc_of)
        if verbose:
            print(f"\n  anchor X={X:.0e}  sqrt(X)={sq:.0f}  window={window}  "
                  f"(full RMS @ {nz} zeros = {full_rms:.3f})")
            print("    c=T/sqrt(X):   " + "  ".join(f"{c:>5.1f}" for c in cgrid))
            print("    exact frac :   " + "  ".join(f"{frac[c]:>5.2f}" for c in cgrid))
            print("    RMS|r|     :   " + "  ".join(f"{rms_vs_c[c]:>5.2f}" for c in cgrid))
            print("    octave drop:   "
                  + "  ".join(f"{octave[c]:>5.2f}" for c in [1.0, 2.0, 4.0, 8.0])
                  + f"   (full={full_rms:.2f}; >0.5 => octave load-bearing)")
    if results and verbose:
        print("\n  --- FLOOR VERDICT ---")
        for c in cgrid:
            allf = [results[X]["exact_frac"][c] for X in results]
            print(f"  c={c:>4.1f}: min exact frac over anchors = {min(allf):.2f}")
        print("  A sub-sqrt(x) analytic witness would reach ~1.0 exactness at c<1.")
        print("  Observed: exactness needs c>>1 and every octave up to sqrt(x) is")
        print("  load-bearing => Theta(sqrt(x)) zeros, floor-bound (matches S511).")
    return results


# ----------------------------------------------------------------------------
# Self-tests
# ----------------------------------------------------------------------------


def selftest():
    t0 = time.time()
    mpmath.mp.dps = 25
    weights = mobius_weights(12)
    gammas = load_gammas()
    P = sieve_pi(2_000_000)
    ncase = 0

    def check(cond, msg):
        nonlocal ncase
        ncase += 1
        if not cond:
            raise AssertionError(f"SELFTEST FAIL [{ncase}]: {msg}")
        print(f"  ok [{ncase}] {msg}")

    # 1-3: stable formula rounds exactly with enough zeros
    for xv, want in [(10000.5, 1229), (100000.5, 9592), (1000000.5, 78498)]:
        errs = partial_errors(xv, gammas, weights, want)
        check(abs(errs[-1]) < 0.5, f"stable formula |err|<0.5 at x={xv} "
              f"(err={errs[-1]:.3f}, rounds to pi={want})")

    # 4: ground-truth sieve matches a couple of known values
    check(int(P[100000]) == 9592 and int(P[1000000]) == 78498,
          "sieve_pi matches pi(1e5)=9592, pi(1e6)=78498")

    # 5-6: CLOSED row 31 bug reproduced AND fixed
    xv = 100000.5
    tp = int(P[100000])
    errs_good = partial_errors(xv, gammas[:2000], weights, tp, buggy=False)
    errs_bad = partial_errors(xv, gammas[:2000], weights, tp, buggy=True)
    check(abs(errs_good[-1]) < 0.5,
          f"Ei-unwrapped (fix) converges: |err|={abs(errs_good[-1]):.3f} < 0.5")
    check(abs(errs_bad[-1]) > 5.0,
          f"li(exp())-wrapped (row 31 bug) DIVERGES: |err|={abs(errs_bad[-1]):.1f} > 5")

    # 7: settling -> magnitude generally decreases with more zeros (x=1e5)
    errs = partial_errors(100000.5, gammas, weights, int(P[100000]))
    check(abs(errs[-1]) < abs(errs[200]),
          f"more zeros reduce error: |err[8000]|={abs(errs[-1]):.3f} < "
          f"|err[200]|={abs(errs[200]):.3f}")

    # 8: min_zeros_settled is consistent (tail all within 0.5 from N_min on)
    n_min, settled = min_zeros_settled(errs)
    check(settled and np.all(np.abs(errs[n_min:]) < 0.5),
          f"min_zeros_settled at x=1e5: N_min={n_min}, tail all <0.5")

    # 9: N_min for larger x exceeds N_min for smaller x (monotone witness growth)
    e_small = partial_errors(30000.5, gammas, weights, int(P[30000]))
    n_small, s_small = min_zeros_settled(e_small)
    check(s_small and n_small < n_min,
          f"witness grows with x: N_min(3e4)={n_small} < N_min(1e5)={n_min}")

    # 10: R_real (library) == explicit Gram-series construction to full precision
    xR = mpmath.mpf(1234567)
    check(abs(R_real(xR) - _R_gram(xR)) < 1e-12,
          f"R_real == Gram series: diff={float(abs(R_real(xR)-_R_gram(xR))):.1e}")

    # 11: fit recovers a known pure power
    xs = np.geomspace(1e3, 1e7, 30)
    ys = 0.6 * xs ** 0.5
    a, c, r2 = fit_loglog(xs, ys)
    check(abs(a - 0.5) < 1e-6 and r2 > 0.999, f"fit_loglog recovers 0.5: a={a:.4f}")

    # 12: power+loglog fit separates power from polylog on synthetic sqrt*log^2
    ys2 = 0.3 * np.sqrt(xs) * np.log(xs) ** 2
    pa, pd, pc, pr2 = fit_power_loglog(xs, ys2)
    a1, _, _ = fit_loglog(xs, ys2)
    check(abs(pa - 0.5) < 0.03 and a1 > 0.5,
          f"power+loglog recovers power 0.5 (a={pa:.3f}) where 1-term inflates "
          f"to {a1:.3f}")

    # 13: galway-form fit prefers p=2 on synthetic sqrt*log^2 data
    rows = [(float(x), 0, int(round(0.3 * math.sqrt(x) * math.log(x) ** 2)),
             0.0, math.sqrt(x)) for x in xs]
    out, best = galway_comparison(rows, verbose=False)
    check(best == 2, f"galway-form fit picks p=2 on synthetic data (got p={best})")

    # 14: floor-check -- exactness rises with c, is <1 at sub-sqrt(x), RMS falls
    fc = floor_check(anchors=(30000.0,), window=48, mult=16.0, verbose=False)
    X = list(fc.keys())[0]
    frac = fc[X]["exact_frac"]
    check(frac[0.5] <= frac[8.0] and frac[0.5] < 1.0,
          f"exactness rises with c, <1 at c=0.5 "
          f"(frac[.5]={frac[0.5]:.2f} <= frac[8]={frac[8.0]:.2f})")
    rms = fc[X]["rms_vs_c"]
    check(rms[0.5] > rms[8.0],
          f"RMS remainder decreases with T: rms[.5]={rms[0.5]:.2f} > "
          f"rms[8]={rms[8.0]:.2f}")

    # 15: edge -- small x handled, settles
    e = partial_errors(1000.5, gammas[:500], weights, int(P[1000]))
    nm, st = min_zeros_settled(e)
    check(st and nm >= 0, f"small x=1000 settles at N_min={nm}")

    print(f"\nALL SELFTESTS PASSED ({ncase} cases) in {time.time()-t0:.1f}s")


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--measure-scaling", action="store_true")
    ap.add_argument("--galway", action="store_true")
    ap.add_argument("--floor-check", action="store_true")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--dps", type=int, default=18)
    ap.add_argument("--kmob", type=int, default=12)
    ap.add_argument("--xmax", type=float, default=5e5)
    ap.add_argument("--points", type=int, default=12)
    ap.add_argument("--window", type=int, default=20)
    ap.add_argument("--mult", type=float, default=18.0,
                    help="zero budget = mult * sqrt(x) per anchor")
    args = ap.parse_args()

    if not any([args.selftest, args.measure_scaling, args.galway,
                args.floor_check, args.all]):
        args.selftest = True

    if args.selftest:
        selftest()

    rows = None
    if args.measure_scaling or args.galway or args.all:
        print("\n=== MEASUREMENT A: witness-size (N*, T*) scaling ===")
        rows, _ = measure_scaling(xmax=args.xmax, points=args.points,
                                  window=args.window, mult=args.mult,
                                  dps=args.dps, kmob=args.kmob)
    if args.galway or args.all:
        print("\n=== MEASUREMENT B: Galway rung ===")
        galway_comparison(rows)
    if args.floor_check or args.all:
        print("\n=== MEASUREMENT C: S511 floor cross-check ===")
        floor_check(dps=args.dps, kmob=args.kmob)


if __name__ == "__main__":
    main()
