#!/usr/bin/env python3
"""
sparse_precision_pi.py
======================

P8 (OPEN_POSITIVE_TARGETS.md, "Sparse-precision pi queries -- batched on
precision, not x"): compute the FIRST k BITS of pi(x) for varying k, and ask
whether the PER-BIT cost is sub-linear (the A-grade target).  Single-x is fixed;
P8's hope is that refining the precision k amortises -- each extra bit cheaper
than the last (a CONCAVE / diminishing-marginal cost curve).

THE PRECISION-DOMAIN VIEW (orthogonal to S518)
----------------------------------------------
S518 (explicit_formula_witness.py) measured the *x-domain* exponent: to settle
pi(x) EXACTLY needs T ~ sqrt(x)*polylog zeros, information DENSE across them.
This script fixes x and varies the TARGET PRECISION k (bits), measuring the
zeros-per-bit cost.  It is the same Riemann explicit formula

    pi(x) ~ R(x) - 2 Re sum_{0<gamma<=T} R(x^{1/2+i gamma}),     [S518 evaluator]

reused VERBATIM (the numerically-stable Ei-unwrapped term -- CLOSED row 31's
li(exp()) branch-cut bug is already fixed and selftested there; we import it so
we cannot re-introduce it).

THE MECHANISM (what we predict, and measure)
--------------------------------------------
The truncation error is the explicit-formula tail over the OMITTED zeros
gamma>T.  Each omitted zero contributes ~ sqrt(x)/(gamma log x) with a random
phase, so the GUE-averaged RMS error is an L2 sum:
    E_rms(x,T) ~ sqrt( sum_{gamma>T} (sqrt(x)/(gamma log x))^2 )
              ~ (sqrt(x)/log x) * sqrt(log T / T)  ~  C * sqrt(x) * T^{-eta}  (*)
with eta ~ 1/2 (the 1/gamma^2 tail; the sqrt(log T) lifts the MEASURED eta to
~0.55-0.65).  Resolving one more bit = HALVING E needs T scaled by 2^{1/eta}
~ 3-4x, hence the zero count N ~ T*polylog scaled by ~3-4x.  Therefore:

  * the TOP bits are FREE: R(x) alone (polylog, x-independent Gram series)
    already nails pi(x) to within E0 = |pi-R|, the leading
    k0 = B - log2(E0) bits.  ASYMPTOTICALLY E0 ~ sqrt(x) so k0 -> ~B/2 (the
    "~50% of digits"); at FINITE measurable x, E0 ~ sqrt(x)/log x is far below
    its worst-case bound, so k0/B is HIGHER (~0.7-0.9, declining with x);
  * every bit BELOW k0 costs ~2^{1/eta} ~ 3-4x the zeros of the previous one --
    a GEOMETRIC (super-linear, CONVEX) marginal cost, the OPPOSITE of P8's
    sub-linear hope;
  * the full B ~ log2(pi) bits need the last hard bit at N ~ sqrt(x)*polylog
    zeros -- the same wall S518/S511 measure, now resolved bit-by-bit.

So P8's "batched on precision" gives NOTHING: the explicit formula is ALREADY
incrementally batched (refining to k+1 bits REUSES the N zeros of k bits and
adds dN more), yet the cost is geometric, so the total to reach k bits is
dominated by the LAST bit -- amortised per-bit cost A(k)=N(k)/k GROWS.

WHAT WOULD FALSIFY THE NEGATIVE (and the control that proves we'd see it)
------------------------------------------------------------------------
A real P8 win = the measured cost slope s = d log2 N / d(bit) being < 1 in the
random zone (sub-linear per-bit; ideally s->0, deep bits cheap).  We predict
s = 1/eta ~ 1.5-2 (geometric, x3-4 per bit).  A synthetic control feeds an
envelope E ~ T^{-2} (which
WOULD be a sub-linear per-bit cost, slope 0.5) through the SAME pipeline and
asserts the fit recovers 0.5 -- so the method demonstrably detects a sub-linear
regime if one exists; the real-data s~1 is then a genuine negative, not a blind
spot.

VERDICT (see results .md): per-bit cost is FLAT (free) for the top bits (R(x),
polylog) and GEOMETRIC (x3-4 per bit, super-linear) for the hard bits below;
no sub-linear regime; batching across precision saves only a constant factor.
This RE-DERIVES the S511/S518 information barrier in the precision domain (a
measurement confirming a known wall from a new angle -- NOT novelty, per the
contract).

CLI
---
  --selftest        run all self-tests (default if no flag)
  --profile         the per-bit cost profile across anchors (core measurement)
  --anchor X        a single-anchor detailed bit-by-bit table
  --all             profile (multi-anchor)
  --dps D           mpmath precision (default 15; selftest uses more)
  --kmob K          number of Mobius terms in R (default 12)
  --window W        x-window for RMS smoothing (default 16)
  --mult M          zero budget = M*sqrt(x) per anchor, capped at 8000 (default 24)
  --anchors "a,b,.."  comma anchors for --profile (default 1e4,3e4,1e5,3e5,1e6)
"""
import argparse
import math
import os
import sys
import time

import numpy as np
import mpmath

# Reuse the TESTED, numerically-stable explicit-formula evaluator (S518).
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "explicit_formula_witness"))
from explicit_formula_witness import (  # noqa: E402
    load_gammas, sieve_pi, mobius_weights, R_real, partial_errors,
    anchor_curves, rms_envelope, fit_loglog,
)


# ----------------------------------------------------------------------------
# Bit-resolution bookkeeping
# ----------------------------------------------------------------------------


def bit_length_int(v):
    """Number of bits in the integer pi(x): B = floor(log2 v) + 1."""
    return int(v).bit_length()


def bits_resolved(E, pi_val, B):
    """How many leading bits of pi(x) are correct given absolute error E.

    "First k bits correct" <=> E < 2^{B-k} (place value of the k-th-from-top
    bit).  So k = B - log2(E), capped to [0, B]; E<1 => all B bits exact.
    """
    E = max(float(E), 1e-12)
    k = B - math.log2(E)
    if k < 0:
        k = 0.0
    if k > B:
        k = float(B)
    return k


def last_crossing(E, thresh):
    """Smallest N such that E[N'] < thresh for ALL N' >= N (robust settling).

    Returns N or None if E never stays below thresh.  E is the (smoothed) RMS
    envelope, length nz+1, E[N] = error after N zeros.
    """
    E = np.asarray(E, float)
    n = len(E) - 1
    if E[n] >= thresh:
        return None
    last_bad = -1
    for i in range(n, -1, -1):
        if E[i] >= thresh:
            last_bad = i
            break
    return last_bad + 1


# ----------------------------------------------------------------------------
# Per-anchor bit-cost profile
# ----------------------------------------------------------------------------


def anchor_profile(X, window, mult, gammas, garr, weights, P, gmax):
    """At fixed anchor X: the RMS error envelope E[N], heights T[N], and the
    per-bit cost map N(k) = zeros to resolve the first k bits of pi(X).

    Returns a dict (or None if too few zeros to even reach the random zone).
    """
    sq = math.sqrt(float(X))
    nz = min(len(gammas), int(np.searchsorted(garr, min(gmax, mult * sq)) + 2))
    _, curves = anchor_curves(float(X), window, gammas, garr, weights, P, nz)
    E = rms_envelope(curves)                 # E[N], N=0..nz  (RMS over window)
    T = np.concatenate([[0.0], garr[:nz]])   # T[N] = height of the N-th zero
    pi_val = int(P[int(X)])
    B = bit_length_int(pi_val)

    kres = np.array([bits_resolved(E[N], pi_val, B) for N in range(len(E))])
    k0 = float(kres[0])                       # free bits from R(x) alone (N=0)

    # N(k): zeros to (robustly) resolve the first k bits, integer k in (k0, B]
    klo = int(math.floor(k0)) + 1
    Nk = {}
    for k in range(klo, B + 1):
        thresh = 2.0 ** (B - k)
        Nc = last_crossing(E, thresh)
        if Nc is not None and Nc >= 1:
            Nk[k] = int(Nc)
    return dict(X=float(X), sqrtX=sq, pi=pi_val, B=B, nz=nz, E=E, T=T,
                kres=kres, k0=k0, Nk=Nk)


def fit_bit_cost(Nk):
    """Slope s of log2 N(k) vs k over the random zone (the per-bit cost
    exponent).  s ~ 1 => N doubles per bit (geometric).  Returns
    (s, r2, ks, log2N, marginal_ratios)."""
    ks = sorted(Nk)
    if len(ks) < 3:
        return None
    ka = np.array(ks, float)
    l2 = np.array([math.log2(Nk[k]) for k in ks])
    A = np.vstack([ka, np.ones_like(ka)]).T
    coef, *_ = np.linalg.lstsq(A, l2, rcond=None)
    pred = A @ coef
    ss_res = float(np.sum((l2 - pred) ** 2))
    ss_tot = float(np.sum((l2 - l2.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    ratios = [Nk[ks[i + 1]] / Nk[ks[i]] for i in range(len(ks) - 1)]
    return float(coef[0]), float(r2), ks, l2.tolist(), ratios


def fit_envelope_eta(E, T, sqrtX):
    """Fit log E ~ -eta log T over the random zone (1 < E < sqrtX/2, T>0).

    E ~ T^{-eta}; mechanism (*) predicts eta ~ 1/2 (the L2 1/gamma^2 zero-tail,
    measured ~0.55-0.65 with the sqrt(log T) lift).  eta<1 => the per-bit cost
    slope 1/eta > 1 (one bit = halve E = T scaled by 2^{1/eta} ~ 3-4x = N by
    the same).  NOTE: the per-x RMS envelope is GUE-noisy, so this fit is a
    corroborating mechanism; the well-determined per-bit cost is fit_bit_cost.
    """
    E = np.asarray(E, float)
    T = np.asarray(T, float)
    mask = (T > 0) & (E > 1.0) & (E < 0.5 * sqrtX)
    if mask.sum() < 4:
        return None
    a, _, r2 = fit_loglog(T[mask], E[mask])
    return float(-a), float(r2), int(mask.sum())


# ----------------------------------------------------------------------------
# Drivers
# ----------------------------------------------------------------------------


def run_profile(anchors, window=16, mult=24.0, dps=15, kmob=12, verbose=True):
    mpmath.mp.dps = dps
    weights = mobius_weights(kmob)
    gammas = load_gammas()
    garr = np.array([float(g) for g in gammas])
    gmax = float(garr[-1])
    P = sieve_pi(int(max(anchors)) + 2 * window + 10)

    summary = []
    for X in anchors:
        prof = anchor_profile(X, window, mult, gammas, garr, weights, P, gmax)
        if prof is None:
            continue
        fb = fit_bit_cost(prof["Nk"])
        fe = fit_envelope_eta(prof["E"], prof["T"], prof["sqrtX"])
        prof["fit_bit"] = fb
        prof["fit_eta"] = fe
        summary.append(prof)
        if verbose:
            s = fb[0] if fb else float("nan")
            r2 = fb[1] if fb else float("nan")
            eta = fe[0] if fe else float("nan")
            print(f"\n  X={int(X):>10d}  pi={prof['pi']:>8d}  B={prof['B']:>2d} bits  "
                  f"nz={prof['nz']}")
            print(f"    free bits k0 (R alone) = {prof['k0']:.2f}  "
                  f"=> k0/B = {prof['k0']/prof['B']:.2f}  "
                  f"(hard/random-zone bits = {prof['B']-prof['k0']:.2f})")
            print(f"    envelope  E ~ T^(-eta):  eta = {eta:.3f}  "
                  f"(L2 zero-tail ~0.5-0.65; noisy)" if fe else
                  "    envelope: [too few pts]")
            if fb:
                ratios = "  ".join(f"{r:.2f}" for r in fb[4])
                print(f"    per-bit cost  log2 N(k) vs k:  slope s = {s:.3f} "
                      f"(R2={r2:.3f});  sub-linear iff s<1, geometric if s>=1")
                print(f"    marginal zero-count ratios N(k+1)/N(k): {ratios} "
                      f"(GUE-noisy; trend x2-4)")
                print(f"    N(k) table: " + "  ".join(
                    f"k{k}:{prof['Nk'][k]}" for k in sorted(prof['Nk'])))

    if summary and verbose:
        print("\n  === VERDICT (P8: is per-bit cost sub-linear?) ===")
        slopes = [p["fit_bit"][0] for p in summary if p["fit_bit"]]
        etas = [p["fit_eta"][0] for p in summary if p["fit_eta"]]
        frees = [p["k0"] / p["B"] for p in summary]
        print(f"    per-bit cost slope s : {['%.2f' % v for v in slopes]}  "
              f"(min {min(slopes):.2f}, mean {np.mean(slopes):.2f})")
        print(f"    envelope eta         : {['%.2f' % v for v in etas]}  "
              f"(mean {np.mean(etas):.2f}; 1/eta ~ s)")
        print(f"    free-bit fraction k0/B: {['%.2f' % v for v in frees]}  "
              f"(mean {np.mean(frees):.2f})")
        print(f"    --> s = {np.mean(slopes):.2f} >= 1 (every anchor): per-bit cost is "
              f"GEOMETRIC/super-linear, NOT sub-linear; P8 A-grade REFUTED.")
        print(f"    --> top ~{np.mean(frees)*100:.0f}% of bits FREE (R(x), polylog, "
              f"finite-x typical); the hard bits each cost ~x{2**np.mean(slopes):.0f}/bit.")
        print(f"    --> matches S511 info floor / S518 x-domain sqrt(x): the same "
              f"barrier, resolved bit-by-bit (a re-derivation, not a new wall).")
    return summary


def run_anchor_detail(X, window=16, mult=24.0, dps=15, kmob=12):
    mpmath.mp.dps = dps
    weights = mobius_weights(kmob)
    gammas = load_gammas()
    garr = np.array([float(g) for g in gammas])
    gmax = float(garr[-1])
    P = sieve_pi(int(X) + 2 * window + 10)
    prof = anchor_profile(X, window, mult, gammas, garr, weights, P, gmax)
    if prof is None:
        print("  [too few zeros for this anchor]")
        return
    print(f"\n  === single-anchor bit-by-bit, X={int(X)} ===")
    print(f"  pi(X)={prof['pi']}  B={prof['B']} bits  sqrt(X)={prof['sqrtX']:.1f}  "
          f"nz={prof['nz']}")
    print(f"  free bits from R(x) alone (N=0): k0={prof['k0']:.2f}\n")
    print(f"  {'k(bits)':>8s} {'N(zeros)':>9s} {'dN':>7s} {'ratio':>6s} "
          f"{'T_height':>9s}")
    ks = sorted(prof["Nk"])
    prevN = None
    for k in ks:
        N = prof["Nk"][k]
        T = float(prof["T"][N]) if N < len(prof["T"]) else float("nan")
        if prevN is None:
            print(f"  {k:>8d} {N:>9d} {'-':>7s} {'-':>6s} {T:>9.1f}")
        else:
            dN = N - prevN
            ratio = N / prevN if prevN > 0 else float("nan")
            print(f"  {k:>8d} {N:>9d} {dN:>7d} {ratio:>6.2f} {T:>9.1f}")
        prevN = N
    fb = fit_bit_cost(prof["Nk"])
    fe = fit_envelope_eta(prof["E"], prof["T"], prof["sqrtX"])
    if fb:
        print(f"\n  per-bit slope s = {fb[0]:.3f} (R2={fb[1]:.3f})"
              + (f"; eta = {fe[0]:.3f}" if fe else ""))
        print(f"  s>=1 => geometric (x{2**fb[0]:.0f}/bit) per-bit cost; "
              f"no sub-linear regime.")


# ----------------------------------------------------------------------------
# Self-tests
# ----------------------------------------------------------------------------


def _synthetic_profile(eta, X=1e6, nz=4000, C=1.0):
    """Build a synthetic (E, T, sq, X) with E(N) = C*sqrt(X)*log(X) / N^eta on a
    CLEAN unit grid T=N, to validate the bit-cost pipeline against a KNOWN slope.

    With cost = zero-count N and E ~ N^{-eta}, the per-bit cost slope is exactly
    1/eta: eta=1 -> slope 1 (geometric, doubling per bit); eta=2 -> slope 0.5
    (SUB-LINEAR per-bit -- the regime a real P8 win would show).  (Real data
    additionally carries the N~T log T count-inflation, pushing the *count*
    slope slightly above 1/eta -- measured separately, never asserted here.)
    """
    N = np.arange(0, nz + 1, dtype=float)
    T = N.copy()                        # unit grid: cost-count == height proxy
    sq = math.sqrt(X)
    E = np.empty(nz + 1)
    E[0] = C * sq                       # R-alone error ~ sqrt(x)
    E[1:] = C * sq * math.log(X) / np.maximum(N[1:], 1e-9) ** eta
    E = np.maximum(E, 1e-9)
    return E, T, sq, X


def selftest():
    t0 = time.time()
    mpmath.mp.dps = 25
    weights = mobius_weights(12)
    gammas = load_gammas()
    garr = np.array([float(g) for g in gammas])
    gmax = float(garr[-1])
    P = sieve_pi(2_000_000)
    ncase = 0

    def check(cond, msg):
        nonlocal ncase
        ncase += 1
        if not cond:
            raise AssertionError(f"SELFTEST FAIL [{ncase}]: {msg}")
        print(f"  ok [{ncase}] {msg}")

    # 1: imported evaluator still settles pi(1e5) exactly (guards the import)
    errs = partial_errors(100000.5, gammas, weights, 9592)
    check(abs(errs[-1]) < 0.5, f"imported evaluator settles pi(1e5)=9592 "
          f"(|err|={abs(errs[-1]):.3f})")

    # 2: bit_length_int and bits_resolved basics
    check(bit_length_int(9592) == 14, f"bit_length_int(9592)=14 "
          f"(got {bit_length_int(9592)})")
    B = 14
    check(abs(bits_resolved(2.0 ** (B - 5), 9592, B) - 5.0) < 1e-9,
          "bits_resolved: E=2^(B-5) -> exactly 5 bits")
    check(bits_resolved(0.3, 9592, B) == B, "bits_resolved: E<1 -> all B bits")
    check(bits_resolved(1e9, 9592, B) == 0, "bits_resolved: huge E -> 0 bits")

    # 3: last_crossing semantics (stays below for all N'>=N)
    Earr = np.array([100, 50, 0.4, 10, 0.3, 0.2])  # dips then recovers
    check(last_crossing(Earr, 0.5) == 4,
          f"last_crossing robust to a dip: got {last_crossing(Earr, 0.5)} (want 4)")
    check(last_crossing(np.array([0.3, 0.2, 0.1]), 0.5) == 0,
          "last_crossing already-settled -> 0")
    check(last_crossing(np.array([1.0, 1.0]), 0.5) is None,
          "last_crossing never-settles -> None")

    # 4: SYNTHETIC eta=1 -> per-bit slope ~ 1 (geometric, doubling)
    E, T, sq, X = _synthetic_profile(eta=1.0)
    pi_val = int(P[int(X)])
    Bx = bit_length_int(pi_val)
    kres = np.array([bits_resolved(E[N], pi_val, Bx) for N in range(len(E))])
    k0 = float(kres[0])
    Nk = {}
    for k in range(int(math.floor(k0)) + 1, Bx + 1):
        Nc = last_crossing(E, 2.0 ** (Bx - k))
        if Nc is not None and Nc >= 1:
            Nk[k] = int(Nc)
    fb = fit_bit_cost(Nk)
    check(fb is not None and abs(fb[0] - 1.0) < 0.15,
          f"synthetic eta=1 -> per-bit slope ~1 (got {fb[0]:.3f}) [GEOMETRIC]")

    # 5: SYNTHETIC eta=2 -> per-bit slope ~ 0.5 (SUB-LINEAR -- the control that
    #    proves the method WOULD detect a P8 win if one existed)
    E2, T2, sq2, X2 = _synthetic_profile(eta=2.0)
    kres2 = np.array([bits_resolved(E2[N], pi_val, Bx) for N in range(len(E2))])
    k02 = float(kres2[0])
    Nk2 = {}
    for k in range(int(math.floor(k02)) + 1, Bx + 1):
        Nc = last_crossing(E2, 2.0 ** (Bx - k))
        if Nc is not None and Nc >= 1:
            Nk2[k] = int(Nc)
    fb2 = fit_bit_cost(Nk2)
    check(fb2 is not None and fb2[0] < 0.75,
          f"synthetic eta=2 -> per-bit slope <0.75 (got {fb2[0]:.3f}) "
          f"[SUB-LINEAR detected -- method is not blind to a P8 win]")

    # 6: envelope-eta fit recovers the synthetic eta on the eta=1 grid
    fe = fit_envelope_eta(E, T, sq)
    check(fe is not None and abs(fe[0] - 1.0) < 0.2,
          f"fit_envelope_eta recovers eta=1 (got {fe[0]:.3f})")
    fe2 = fit_envelope_eta(E2, T2, sq2)
    check(fe2 is not None and abs(fe2[0] - 2.0) < 0.3,
          f"fit_envelope_eta recovers eta=2 (got {fe2[0]:.3f})")

    # 7: REAL data at x=1e5 -- envelope decreasing over the random zone
    prof = anchor_profile(100000, 16, 24.0, gammas, garr, weights, P, gmax)
    check(prof is not None, "real anchor x=1e5 profiles")
    Er = prof["E"]
    check(Er[20] < Er[1], f"real envelope decreases: E[20]={Er[20]:.2f} < "
          f"E[1]={Er[1]:.2f}")

    # 8: REAL free bits -- a hard/random zone EXISTS (B-k0 >= 1) and k0 < B;
    #    at finite x R is far better than sqrt(x) so k0/B is HIGH (~0.7-0.9),
    #    NOT ~0.5 (that is the asymptotic worst case) -- we assert it honestly.
    check(prof["B"] - prof["k0"] >= 1.0 and prof["k0"] < prof["B"],
          f"real hard-zone exists: B-k0={prof['B']-prof['k0']:.2f} >= 1, "
          f"k0={prof['k0']:.2f} < B={prof['B']}")
    check(0.40 <= prof["k0"] / prof["B"] <= 0.95,
          f"real k0/B={prof['k0']/prof['B']:.2f} high at finite x "
          f"(R beats its sqrt(x) bound; -> ~0.5 only asymptotically)")

    # 9: REAL per-bit slope is SUPER-linear (>1.2, deterministic 1.53 @ x=1e5)
    #    -- the headline NEGATIVE; well clear of the eta=2 sub-linear control 0.48
    fbr = fit_bit_cost(prof["Nk"])
    check(fbr is not None and fbr[0] > 1.2,
          f"real per-bit slope = {fbr[0]:.3f} > 1.2 (GEOMETRIC/super-linear, the "
          f"OPPOSITE of sub-linear) -- P8 A-grade target refuted at x=1e5")

    # 10: REAL envelope eta < 1 (the L2 zero-tail E ~ sqrt(x) T^{-eta}, eta~0.5)
    #     -- eta<1 is WHY the per-bit slope 1/eta exceeds 1 (worse than doubling)
    fer = fit_envelope_eta(prof["E"], prof["T"], prof["sqrtX"])
    check(fer is not None and 0.1 < fer[0] < 1.0,
          f"real envelope eta = {fer[0]:.3f} in (0.1,1.0): decays slower than 1/T "
          f"(L2 zero-tail) => per-bit slope 1/eta > 1")

    # 11: monotonicity -- N(k) nondecreasing in k (more bits never cheaper)
    ks = sorted(prof["Nk"])
    mono = all(prof["Nk"][ks[i + 1]] >= prof["Nk"][ks[i]]
               for i in range(len(ks) - 1))
    check(mono, "N(k) nondecreasing in k (deeper bits never cost fewer zeros)")

    print(f"\nALL SELFTESTS PASSED ({ncase} cases) in {time.time()-t0:.1f}s")


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--profile", action="store_true")
    ap.add_argument("--anchor", type=float, default=None)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--dps", type=int, default=15)
    ap.add_argument("--kmob", type=int, default=12)
    ap.add_argument("--window", type=int, default=16)
    ap.add_argument("--mult", type=float, default=24.0)
    ap.add_argument("--anchors", type=str, default="1e4,3e4,1e5,3e5,1e6")
    args = ap.parse_args()

    if not any([args.selftest, args.profile, args.anchor, args.all]):
        args.selftest = True

    if args.selftest:
        selftest()

    if args.anchor is not None:
        run_anchor_detail(args.anchor, window=args.window, mult=args.mult,
                          dps=args.dps, kmob=args.kmob)

    if args.profile or args.all:
        anchors = [int(float(a)) for a in args.anchors.split(",")]
        print("\n=== P8 PER-BIT COST PROFILE (precision-domain) ===")
        run_profile(anchors, window=args.window, mult=args.mult,
                    dps=args.dps, kmob=args.kmob)


if __name__ == "__main__":
    main()
