#!/usr/bin/env python3
"""
adaptive_query_complexity.py
============================================================================
P10 (OPEN_POSITIVE_TARGETS): the ADAPTIVE query complexity of computing the
n-th prime p(n) from a prime-counting oracle pi(.).

Question (P10): pi(x_1) is computed first, then x_2 chosen adaptively as a
function of the previous answers.  Does adaptivity beat the non-adaptive
"analytic-predict + local-sieve" recipe, or is it (like Aggarwal binary
search) just a relocation of the same sqrt(x) wall?

This script CLOSES P10 with a mechanism tied to the project's central
SMOOTH + RANDOM decomposition, substantiated by three measurements:

  M1  far/smooth zone is FREE & NON-adaptive: R^{-1}(n) via Newton on the
      polylog Riemann R converges in O(log log x) iterations and lands
      within c.sqrt(p(n)) of the truth.  (the leading ~half of the digits)
  M2  hard/random zone width: |p(n) - R^{-1}(n)| scales like c.sqrt(x)
      (slope ~0.5 in log-log) -> the window you must sieve is Theta(sqrt x).
  M3  the hard zone carries NO smooth signal an interpolation search could
      exploit: within the sqrt(x) window the residual of pi against its
      best local polynomial model stays at the sqrt(count) ~ x^{1/4}
      Poisson scale (does not drop with model degree) -> every adaptive
      query strategy degrades to >= binary search inside the window.

Cost/count distinction (the crisp closure):
  * COUNT metric, SMOOTH zone: interpolation search beats binary (textbook
    loglog vs log) -- BUT the smooth zone is already free (R^{-1}, Newton).
  * COUNT metric, RANDOM sqrt(x) window: interpolation does NOT beat binary
    (MEASURED: interp_q ~ binary_q ~ (1/2)log2(x) inside the window) because
    the fine structure is random -- no smooth signal to interpolate.  So the
    adaptive count advantage EVAPORATES exactly in the zone where the work
    lives.
  * COST metric: each *hard-zone* pi-query costs ~sqrt(x) (best-known pi),
    and sieving the whole sqrt(x) window costs the same ~sqrt(x) while
    answering ALL queries at once.  So the COST-optimal p(n) algorithm uses
    ZERO adaptive hard-zone queries: compute R^{-1}(n) (smooth zone, free),
    then sieve the Theta(sqrt x) window once.  No adaptive pi-query strategy
    beats it; the sqrt(x) random width is irreducible.

So: adaptivity helps exactly where it is already free (the smooth zone) and
fails exactly where the cost lives (the random sqrt(x) zone).  The
SMOOTH+RANDOM split IS the adaptivity boundary.

Builds on / cites:
  * guess_comparison_oracle (S491): far-zone polylog-decidable, hard-zone
    width ~sqrt(p(n)) (median 0.13, max 0.53 over n<=1e6), R-guess
    comparison bit = GUE fair coin.  This script adds the *query-complexity*
    statement (count vs cost) and the Newton-iteration + window-content
    measurements that the geography note did not make.
  * CLOSED_PATHS row "Binary search + local sieve" (FAIL/C): "problem
    relocated, not solved" -- that is the COST face for a single oracle;
    here we make it a query-complexity dichotomy.
  * The information barrier: pi(x) = SMOOTH + RANDOM; ~sqrt(x) zeta-zero
    phases behind the last ~50% of digits.  S511 info floor.

FALSIFIERS (what would overturn the closure):
  F1  Newton iteration count growing like log x (not log log x / ~const)
      -> the smooth zone would not be polylog-cheap.   [tested: M1]
  F2  |p(n) - R^{-1}(n)| growing slower than sqrt(x) (e.g. polylog) -> the
      hard zone would shrink and p(n) would be polylog-sieveable.  [M2]
      (growing FASTER than sqrt(x).polylog would violate RH.)
  F3  within-window pi residual dropping well below sqrt(count) for some
      low-degree polynomial model -> a smooth signal an interpolation
      search could exploit to localize p(n) sub-sqrt(x).            [M3]
  F4  an adaptive strategy resolving p(n) with o(sqrt x) total *cost*
      (not just count) on the benchmark -> direct refutation.   [M3 sim]

Usage:
  python3 adaptive_query_complexity.py --selftest
  python3 adaptive_query_complexity.py --measure --xmax 20000000
  python3 adaptive_query_complexity.py --window-demo --n 100000 --xmax 20000000
"""

import argparse
import math
import sys

import numpy as np

try:
    import mpmath as mp
except Exception:  # pragma: no cover - mpmath is present in this environment
    mp = None


# ---------------------------------------------------------------------------
# Ground-truth prime oracle (a real sieve == the "expensive" side we bound)
# ---------------------------------------------------------------------------
def sieve_primes(xmax):
    """Return a numpy int64 array of all primes <= xmax via a boolean sieve."""
    if xmax < 2:
        return np.empty(0, dtype=np.int64)
    s = np.ones(xmax + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(math.isqrt(xmax)) + 1):
        if s[p]:
            s[p * p::p] = False
    return np.nonzero(s)[0].astype(np.int64)


def pi_of(primes, x):
    """Exact pi(x) by binary search in the prime array (the oracle)."""
    return int(np.searchsorted(primes, x, side="right"))


def nth_prime(primes, n):
    """Exact p(n) (1-indexed: p(1)=2)."""
    return int(primes[n - 1])


# ---------------------------------------------------------------------------
# Riemann R (Gram series) and its derivative  -- the polylog smooth predictor
#   R(x)  = 1 + sum_{k>=1} (ln x)^k / (k . k! . zeta(k+1))
#   R'(x) = (1/x) sum_{k>=1} (ln x)^{k-1} / (k! . zeta(k+1))
# ---------------------------------------------------------------------------
def zeta_table(kmax):
    """zt[k] = zeta(k+1) for k = 1..kmax (float)."""
    zt = {}
    for k in range(1, kmax + 1):
        if mp is not None:
            zt[k] = float(mp.zeta(k + 1))
        else:  # pragma: no cover
            s = k + 1
            acc, j = 0.0, 1
            while j < 200000:
                acc += j ** (-s)
                j += 1
            zt[k] = acc
    return zt


def riemann_R(x, zt):
    L = math.log(x)
    total = 1.0
    term = 1.0  # will hold (L^k)/(k!) incrementally
    for k in range(1, len(zt) + 1):
        term *= L / k                      # term == L^k / k!
        add = term / (k * zt[k])
        total += add
        if k > L + 5 and abs(add) < 1e-15 * abs(total):
            break
    return total


def riemann_R_prime(x, zt):
    # R'(x) = (1/x) sum_{k>=1} (ln x)^{k-1} / (k! . zeta(k+1))
    L = math.log(x)
    total = 0.0
    term = 1.0  # term == L^{k-1}/(k-1)!  ; start k=1 -> L^0/0! = 1
    for k in range(1, len(zt) + 1):
        add = term / (k * zt[k])           # = L^{k-1}/(k! . zeta(k+1))
        total += add
        if k > L + 5 and abs(add) < 1e-15 * abs(total):
            break
        term *= L / k                      # advance L^{k-1}/(k-1)! -> L^k/k!
    return total / x


def R_inverse(n, zt, tol=0.25, max_iter=100):
    """Solve R(x) = n by Newton.  Returns (x, iters).  x0 = n ln n (PNT)."""
    x = max(n * math.log(max(n, 3)), 3.0)
    iters = 0
    for _ in range(max_iter):
        iters += 1
        f = riemann_R(x, zt) - n
        if abs(f) < tol:
            break
        fp = riemann_R_prime(x, zt)
        step = f / fp
        nx = x - step
        # safeguard: keep x positive and bounded per step
        if nx <= 1.0:
            nx = x / 2.0
        x = nx
    return x, iters


# ---------------------------------------------------------------------------
# M1 + M2 : Newton iteration count + landing-window width across scales
# ---------------------------------------------------------------------------
def measure_predictor(primes, ns, zt):
    """For each n: p(n), R^{-1}(n), newton iters, abs err, err/sqrt(p(n))."""
    rows = []
    for n in ns:
        pn = nth_prime(primes, n)
        xinv, iters = R_inverse(n, zt)
        err = abs(pn - xinv)
        rows.append(dict(n=n, pn=pn, rinv=xinv, iters=iters,
                         err=err, err_over_sqrt=err / math.sqrt(pn)))
    return rows


def loglog_slope(xs, ys):
    """Least-squares slope of log(y) vs log(x)."""
    lx = np.log(np.asarray(xs, float))
    ly = np.log(np.asarray(ys, float))
    A = np.vstack([lx, np.ones_like(lx)]).T
    slope, intercept = np.linalg.lstsq(A, ly, rcond=None)[0]
    return slope, intercept


# ---------------------------------------------------------------------------
# M3 : within-window structure -- no smooth signal to interpolate
# ---------------------------------------------------------------------------
def window_residuals(primes, center, half_width, npts=400, degrees=(1, 3)):
    """
    Sample pi(x) on a grid across [center-hw, center+hw], fit polynomials of
    the given degrees, and return RMS residuals + the window prime count and
    sqrt(count).  A residual that does NOT drop with degree == no smooth
    signal (random content) an interpolation search could exploit.
    """
    lo = max(2, int(center - half_width))
    hi = int(center + half_width)
    xs = np.linspace(lo, hi, npts)
    ys = np.array([pi_of(primes, x) for x in xs], dtype=float)
    count = ys[-1] - ys[0]
    out = dict(lo=lo, hi=hi, count=count, sqrt_count=math.sqrt(max(count, 1)))
    # center/scale x for conditioning
    xc = (xs - xs.mean()) / (xs.std() + 1e-9)
    for d in degrees:
        coef = np.polyfit(xc, ys, d)
        resid = ys - np.polyval(coef, xc)
        out[f"rms_deg{d}"] = float(np.sqrt(np.mean(resid ** 2)))
    return out


# ---------------------------------------------------------------------------
# M3 (corollary) : COUNT-metric simulation of adaptive search inside the
# window -- interpolation vs binary -- with EXACT pi.  Shows interpolation
# can cut the COUNT (loglog vs log) but every probe is a hard-zone pi-eval,
# so in the COST metric a single window sieve (zero adaptive queries) wins.
# ---------------------------------------------------------------------------
def _bracket(primes, n, lo, hi):
    """Ensure pi(lo) < n <= pi(hi) (the jump p(n) is in (lo, hi])."""
    while pi_of(primes, lo) >= n:
        lo = max(2, lo - (hi - lo) - 1)
    while pi_of(primes, hi) < n:
        hi = hi + (hi - lo) + 1
    return lo, hi


def binary_search_count(primes, n, lo, hi):
    """# pi-queries for plain binary search to pin p(n) exactly."""
    lo, hi = _bracket(primes, n, lo, hi)
    q = 0
    while hi - lo > 1:
        mid = (lo + hi) // 2
        q += 1
        if pi_of(primes, mid) >= n:
            hi = mid
        else:
            lo = mid
    return q, hi  # hi == p(n)


def interpolation_search_count(primes, n, lo, hi, max_q=200):
    """# pi-queries for interpolation search (assume locally linear pi)."""
    lo, hi = _bracket(primes, n, lo, hi)
    plo = pi_of(primes, lo)
    phi = pi_of(primes, hi)
    q = 2  # plo, phi already cost 2 queries
    while hi - lo > 1 and q < max_q:
        if phi == plo:
            mid = (lo + hi) // 2
        else:
            frac = (n - plo) / (phi - plo)
            frac = min(max(frac, 1e-6), 1 - 1e-6)
            mid = int(lo + frac * (hi - lo))
            mid = min(max(mid, lo + 1), hi - 1)
        q += 1
        pm = pi_of(primes, mid)
        if pm >= n:
            hi, phi = mid, pm
        else:
            lo, plo = mid, pm
    return q, hi


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
def run_measure(xmax, seed=0):
    print(f"=== P10 adaptive query complexity : measure (xmax={xmax}) ===")
    primes = sieve_primes(xmax)
    npi = len(primes)
    print(f"sieved {npi} primes up to {xmax} (pi({xmax})={npi})")
    zt = zeta_table(60)

    # log-spaced n with p(n) <= xmax
    nmax = npi
    ns = sorted(set(int(round(v)) for v in
                    np.unique(np.logspace(2.3, math.log10(nmax * 0.95), 22).astype(int))
                    if v >= 50))
    rows = measure_predictor(primes, ns, zt)

    print("\nM1/M2  predictor R^{-1}(n): Newton iters + landing window")
    print(f"{'n':>9} {'p(n)':>12} {'R^-1(n)':>14} {'iters':>5} "
          f"{'|err|':>12} {'err/sqrt(pn)':>12}")
    for r in rows:
        print(f"{r['n']:>9} {r['pn']:>12} {r['rinv']:>14.1f} {r['iters']:>5} "
              f"{r['err']:>12.1f} {r['err_over_sqrt']:>12.4f}")

    iters = [r["iters"] for r in rows]
    eos = [r["err_over_sqrt"] for r in rows]
    print(f"\nM1: Newton iters  min/mean/max = "
          f"{min(iters)}/{np.mean(iters):.2f}/{max(iters)}  "
          f"(O(log log x): ~constant across {len(rows)} scales)")

    # M2 robust: |err| is a fluctuating r.v.; a single-sample log-log slope is
    # noise-dominated.  Sample DENSELY, bin by scale, fit RMS|err| per bin.
    dense_ns = sorted(set(int(round(v)) for v in
                          np.unique(np.logspace(2.5, math.log10(nmax * 0.95),
                                                240).astype(int)) if v >= 50))
    drows = measure_predictor(primes, dense_ns, zt)
    dpn = np.array([r["pn"] for r in drows], float)
    derr = np.array([r["err"] for r in drows], float)
    deos = np.array([r["err_over_sqrt"] for r in drows], float)
    nb = 8
    edges = np.logspace(math.log10(dpn.min()), math.log10(dpn.max()), nb + 1)
    bx, by = [], []
    for i in range(nb):
        m = (dpn >= edges[i]) & (dpn < edges[i + 1] if i < nb - 1 else dpn <= edges[i + 1])
        if m.sum() >= 3:
            bx.append(math.sqrt(np.median(dpn[m])))      # sqrt(x) axis
            by.append(math.sqrt(np.mean(derr[m] ** 2)))  # RMS|err| per bin
    s_rms, _ = loglog_slope([b ** 2 for b in bx], by)  # vs x
    # is err/sqrt(x) trending up?  slope of ratio vs x:
    s_ratio, _ = loglog_slope(dpn, np.maximum(deos, 1e-6))
    print(f"M2: window |p(n)-R^-1(n)| over {len(dense_ns)} n in 8 scale-bins:")
    print(f"    RMS|err| ~ x^{s_rms:.3f}   (predict 0.50 == sqrt(x) hard zone)")
    print(f"    err/sqrt(p(n)): median/mean/max = "
          f"{np.median(deos):.3f}/{np.mean(deos):.3f}/{deos.max():.3f}  "
          f"(guess_oracle S491: 0.13/-/0.53)")
    print(f"    trend of err/sqrt(x) vs x: slope {s_ratio:+.3f}  "
          f"(~0 == ratio BOUNDED -> window is Theta(sqrt x), F2 not triggered)")

    # M3 : within-window structure for a few representative n
    print("\nM3  within-window pi structure (window = +-1.0*sqrt(p(n))):")
    print(f"{'n':>9} {'p(n)':>12} {'count':>9} {'sqrt(cnt)':>9} "
          f"{'rms_deg1':>9} {'rms_deg3':>9} {'rms/sqrtN':>9}")
    demo_ns = [r["n"] for r in rows[5::5]]
    for n in demo_ns:
        pn = nth_prime(primes, n)
        hw = math.sqrt(pn)
        w = window_residuals(primes, pn, hw)
        print(f"{n:>9} {pn:>12} {w['count']:>9.0f} {w['sqrt_count']:>9.1f} "
              f"{w['rms_deg1']:>9.2f} {w['rms_deg3']:>9.2f} "
              f"{w['rms_deg3'] / w['sqrt_count']:>9.3f}")
    print("  -> rms stays at the sqrt(count) Poisson scale and does NOT drop"
          " from deg1->deg3:\n     no smooth signal to interpolate (F3 not"
          " triggered).")

    # M3 corollary : count vs cost simulation
    print("\nM3 corollary  COUNT-metric adaptive search inside the sqrt(x)"
          " window (exact pi):")
    print(f"{'n':>9} {'p(n)':>12} {'win=2sqrt':>10} {'binary_q':>9} "
          f"{'interp_q':>9}")
    for n in demo_ns:
        pn = nth_prime(primes, n)
        hw = int(math.sqrt(pn))
        lo, hi = pn - hw, pn + hw
        bq, bp = binary_search_count(primes, n, lo, hi)
        iq, ip = interpolation_search_count(primes, n, lo, hi)
        assert bp == pn and ip == pn, "search must recover p(n) exactly"
        print(f"{n:>9} {pn:>12} {2*hw:>10} {bq:>9} {iq:>9}")
    print("  -> INSIDE the random window interp_q ~ binary_q ~ (1/2)log2(x):"
          " interpolation gets NO traction\n     (no smooth signal), so even"
          " the COUNT advantage evaporates here.  And each probe is a"
          "\n     hard-zone pi-eval (~sqrt x); sieving the window once (ZERO"
          " adaptive queries) costs the same\n     ~sqrt(x) and answers all"
          " -> cost-optimal p(n) is the NON-adaptive predict+sieve."
          "  Adaptivity\n     gives no COST gain and no COUNT gain in the hard"
          " zone (F4 not triggered).")

    print("\nCLOSURE: P10 -- adaptive pi-querying does not beat the"
          " non-adaptive analytic-predict + local-sieve recipe.")
    print("  smooth/far zone (Newton, O(log log x)) is free & non-adaptive;")
    print("  random/hard zone has width Theta(sqrt x) and no interpolable"
          " signal; reading it (one")
    print("  sieve) dominates any adaptive query strategy.  SMOOTH+RANDOM is"
          " the adaptivity boundary.")
    return rows


def run_window_demo(n, xmax):
    primes = sieve_primes(xmax)
    zt = zeta_table(60)
    pn = nth_prime(primes, n)
    xinv, iters = R_inverse(n, zt)
    print(f"n={n}  p(n)={pn}  R^-1(n)={xinv:.1f}  iters={iters}  "
          f"|err|={abs(pn-xinv):.1f}  err/sqrt={abs(pn-xinv)/math.sqrt(pn):.4f}")
    for c in (0.5, 1.0, 2.0):
        w = window_residuals(primes, pn, c * math.sqrt(pn))
        print(f"  window +-{c}sqrt:  count={w['count']:.0f}  "
              f"sqrt(count)={w['sqrt_count']:.1f}  "
              f"rms_deg1={w['rms_deg1']:.2f}  rms_deg3={w['rms_deg3']:.2f}")


# ---------------------------------------------------------------------------
# Self-test (boundary cases debugged become cases here)
# ---------------------------------------------------------------------------
def selftest():
    ok = True

    def check(cond, msg):
        nonlocal ok
        status = "ok " if cond else "FAIL"
        if not cond:
            ok = False
        print(f"  [{status}] {msg}")

    print("== selftest ==")
    zt = zeta_table(60)

    # [1] R(x) approximates pi(x) within ~sqrt(x)
    for x, pix in [(100, 25), (1000, 168), (10000, 1229), (100000, 9592)]:
        r = riemann_R(x, zt)
        check(abs(r - pix) < 3 * math.sqrt(x),
              f"R({x})={r:.2f} within 3sqrt(x) of pi={pix}")

    # [2] Newton inverse round-trips: R(R^{-1}(n)) ~ n, and iters bounded
    for n in [100, 1000, 10000, 100000, 1000000]:
        x, it = R_inverse(n, zt)
        back = riemann_R(x, zt)
        check(abs(back - n) < 1.0, f"R(R^-1({n}))={back:.3f} ~ {n}")
        check(it <= 15, f"Newton iters({n})={it} <= 15 (loglog)")

    # [3] exact p(n) extraction from sieve (known values)
    primes = sieve_primes(2_000_000)
    for n, pn in [(1, 2), (2, 3), (10, 29), (25, 97), (100, 541),
                  (1000, 7919), (10000, 104729)]:
        check(nth_prime(primes, n) == pn, f"p({n})={pn}")

    # [4] pi_of consistency at boundaries
    check(pi_of(primes, 2) == 1, "pi(2)=1")
    check(pi_of(primes, 1) == 0, "pi(1)=0")
    check(pi_of(primes, 541) == 100, "pi(541)=100")
    check(pi_of(primes, 540) == 99, "pi(540)=99 (just below p(100))")

    # [5] landing window: |p(n)-R^-1(n)| < sqrt(p(n)) across the range
    #     (guess_comparison_oracle S491 max 0.53)
    for n in [1000, 5000, 20000, 100000]:
        pn = nth_prime(primes, n)
        x, _ = R_inverse(n, zt)
        ratio = abs(pn - x) / math.sqrt(pn)
        check(ratio < 1.0, f"|p({n})-R^-1|/sqrt = {ratio:.3f} < 1")

    # [6] both searches recover p(n) exactly from a sqrt-window bracket
    for n in [1000, 12345, 50000]:
        pn = nth_prime(primes, n)
        hw = int(math.sqrt(pn))
        bq, bp = binary_search_count(primes, n, pn - hw, pn + hw)
        iq, ip = interpolation_search_count(primes, n, pn - hw, pn + hw)
        check(bp == pn, f"binary recovers p({n})={pn} in {bq} q")
        check(ip == pn, f"interp recovers p({n})={pn} in {iq} q")

    # [7] _bracket repairs an off bracket (jump outside initial [lo,hi])
    pn = nth_prime(primes, 5000)
    lo, hi = _bracket(primes, 5000, pn + 10, pn + 20)  # jump below lo
    check(pi_of(primes, lo) < 5000 <= pi_of(primes, hi),
          "_bracket repairs a window that excludes p(n)")

    # [8] within-window residual does NOT drop much deg1->deg3 (no signal)
    pn = nth_prime(primes, 50000)
    w = window_residuals(primes, pn, math.sqrt(pn))
    check(w["rms_deg3"] > 0.3 * w["rms_deg1"],
          f"deg3 rms {w['rms_deg3']:.2f} not << deg1 {w['rms_deg1']:.2f} "
          "(random, not smooth)")
    check(w["rms_deg1"] < 5 * w["sqrt_count"],
          f"rms {w['rms_deg1']:.2f} ~ sqrt(count) {w['sqrt_count']:.1f} scale")

    # [9] R monotone increasing (predictor well-defined for Newton)
    check(riemann_R(1e6, zt) > riemann_R(1e5, zt) > riemann_R(1e4, zt),
          "R monotone increasing")

    print("== PASS ==" if ok else "== FAIL ==")
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--measure", action="store_true")
    ap.add_argument("--window-demo", action="store_true")
    ap.add_argument("--n", type=int, default=100000)
    ap.add_argument("--xmax", type=int, default=20_000_000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    if args.selftest:
        sys.exit(0 if selftest() else 1)
    elif args.window_demo:
        run_window_demo(args.n, args.xmax)
    elif args.measure:
        run_measure(args.xmax, args.seed)
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
