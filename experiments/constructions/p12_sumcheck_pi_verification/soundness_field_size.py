#!/usr/bin/env python3
"""
soundness_field_size.py -- the FIELD-SIZE requirement of the pi(x) verification
chain, and which constraint (correctness vs soundness) binds at large x.

The compressed pi(x) verification chain (compressed_layer.run_chain) is a
polynomial IOP over a prime field F_q, made non-interactive by Fiat-Shamir.
Two DISTINCT lower bounds on q govern whether it works at a given x:

  (1) CORRECTNESS / integer-range.  The trace identities (U*a + R - X = 0, the
      bit recompositions) carry integers up to ~X; if q <= X a wrong quotient
      can ALIAS the truth mod q (a different Lv-bit u'' with u''*a+r'' = X + k*q).
      Empirically q > X suffices (n=8: q=257 > X=255 accepts honestly).
          ==> log2 q  >  log2 X = n = log2 x   (+ a tiny measured slack).
      [This face is already CLOSED: status/CLOSED_PATHS.md row 498, the F_{q^2}
       refutation.  We only re-measure the slack here, to place the crossover.]

  (2) SOUNDNESS / Schwartz-Zippel.  Every random-challenge round binds a
      polynomial identity of degree <= d; a cheating prover false-accepts that
      round with prob <= d/q (DeMillo-Lipton-Schwartz-Zippel).  Over the whole
      chain there are R(x) such field challenges, carrying a degree budget
      B(x) = sum of per-round degrees <= R*d_max.  Union bound:
          total field-soundness error  <=  B(x) / q.
      For lambda-bit security (error < 2^-lambda):
          log2 q  >  lambda + log2 B(x).
      [NOT previously measured.  Row 498 notes the element count "only shrinks
       Schwartz-Zippel error, orthogonal to the integer-range hole" -- this
       script quantifies that orthogonal axis.]

What this script MEASURES (all RAM-light, cache-resident -- safe to run beside
the detached large-x reach):
  --budget      R(x), B(x), d_max via deterministic instrumentation of the REAL
                run_chain (a counting rng + a sum-check degree tally), fit the
                exponent in x.  Expect R,B = Theta(sqrt x * polylog).
  --sz-atomic   the per-round 1/q law on the deg-3 trace zero-test: build the
                witness ONCE, corrupt it, re-verify with many fresh challenge
                sets at several primes q>X, confirm accept-rate ~ c/q.
  --sz-chain    end-to-end composition: run the REAL chain with a cheat over
                many seeds at small q>X, confirm the false-accept rate is
                <= B/q and tracks 1/q (no soundness AMPLIFICATION past the
                union bound).  This is the first empirical chain-soundness
                audit at non-negligible error (prior cheat tests are single-shot
                at BIG_Q where the error is ~0 and a 10x inflation is invisible).
  --extrapolate the two log2 q curves, their crossover, and the p(10^100) need.

HEADLINE (see _results.md): the soundness budget is B = K * c_BK with K =
pi(sqrt x) layers (Theta(sqrt x / log x)) and c_BK = B/K ~ linear in log x, so
B = Theta(sqrt x) (the 1/log x in the prime count CANCELS the log x in the
per-layer round budget; measured single-power alpha ~ 0.55).  Hence soundness
needs log2 q ~ lambda + n/2 - log2 n + O(log n) while correctness needs
log2 q ~ n.  These CROSS at n ~ 2*lambda (x ~ 10^63 for lambda=100): below the
crossover the lambda-bit security floor binds; ABOVE it (incl. the p(10^100)
target, n~332) the integer-range/correctness requirement binds, q ~ 333 bits,
and that field is automatically MORE than large enough for 100-bit soundness
(q_soundness ~ 272 bits < 333).

FALSIFIERS (each would overturn a stated claim):
  F1  B(x) exponent in x reads clearly > 0.5 (e.g. > 0.6 with the 2-term fit)
      -> soundness budget is NOT sqrt(x)*polylog; the crossover analysis is wrong.
  F2  the atomic accept-rate does NOT halve when q doubles (slope of log(rate)
      vs log(q) far from -1) -> the per-round Schwartz-Zippel model is wrong in
      this implementation (a soundness bug or correlated challenges).
  F3  a chain cheat's measured false-accept rate EXCEEDS B/q at some q>X
      -> the union bound is violated (soundness amplification / a skipped check).
  F4  honest run_chain REJECTS at some prime q in (X, big] -> correctness needs
      q strictly larger than X (the q>X slack is more than the measured ~1 bit),
      moving the crossover.
"""

import argparse
import math
import sys
import time

import numpy as np

import compressed_prover_mult_trace as cpmt
import compressed_layer as CL

Q = cpmt.DEFAULT_Q
BIG_Q = cpmt.BIG_Q


# --------------------------------------------------------------------------
# tiny prime helper (avoid a hard sympy dep in the hot path; sympy only for
# generating the small probe primes once)
# --------------------------------------------------------------------------

def _is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    i = 3
    while i * i <= n:
        if n % i == 0:
            return False
        i += 2
    return n != 1


def next_prime(n):
    """Smallest prime > n."""
    m = n + 1
    if m <= 2:
        return 2
    if m % 2 == 0:
        m += 1
    while not _is_prime(m):
        m += 2
    return m


# --------------------------------------------------------------------------
# (1) deterministic instrumentation: count field challenges R and degree
#     budget B of the REAL run_chain.
# --------------------------------------------------------------------------

class CountRng:
    """Wraps np.random.default_rng(seed); counts every .integers() call (each
    is exactly one scalar field challenge in this protocol) and forwards
    everything else verbatim, so the wrapped chain verdict is bit-identical."""

    def __init__(self, seed):
        self._rng = np.random.default_rng(seed)
        self.n_chal = 0

    def integers(self, *a, **k):
        self.n_chal += 1
        return self._rng.integers(*a, **k)

    def __getattr__(self, name):
        return getattr(self._rng, name)


def patch_sumcheck(acc):
    """Replace `sumcheck` in every module that bound it with a tally wrapper
    accumulating (calls, rounds, degree-budget, max degree).  Returns a
    restore() callable.  Mirrors prover_opcount_scaling.py's per-module patch
    (a `from X import sumcheck` makes a private binding in each module)."""
    orig = cpmt.sumcheck

    def wrapped(claim, tables, terms, degree, rng, q=cpmt.DEFAULT_Q):
        nvars = int(round(math.log2(len(next(iter(tables.values()))))))
        acc["calls"] += 1
        acc["rounds"] += nvars
        acc["budget"] += nvars * degree           # per-round Schwartz-Zippel deg
        acc["dmax"] = max(acc["dmax"], degree)
        return orig(claim, tables, terms, degree, rng, q)

    patched = []
    for mod in list(sys.modules.values()):
        if mod is None:
            continue
        try:
            if getattr(mod, "sumcheck", None) is orig:
                setattr(mod, "sumcheck", wrapped)
                patched.append(mod)
        except Exception:
            pass

    def restore():
        for mod in patched:
            setattr(mod, "sumcheck", orig)

    return restore


def measure_budget(n, q=BIG_Q, seed=1, **cfg):
    """Run honest run_chain at x=2^n-1 over field q and return the field-
    soundness accounting: R (challenges), rounds, B (degree budget), d_max,
    plus the chain's own structural numbers.  Deterministic in (n,q,cfg)."""
    x = (1 << n) - 1
    acc = {"calls": 0, "rounds": 0, "budget": 0, "dmax": 0}
    restore = patch_sumcheck(acc)
    try:
        crng = CountRng(seed)
        res = CL.run_chain(x, crng, q=q, **cfg)
    finally:
        restore()
    R = crng.n_chal
    # non-sum-check challenges (evaluation points: tstar, z, tau, rho, ...) each
    # bind a multilinear identity at degree 1 -- add them to the budget at deg 1.
    nonsc = R - acc["rounds"]
    B = acc["budget"] + max(0, nonsc) * 1
    return dict(n=n, x=x, q=q, accepted=bool(res["accepted"]),
                claimed=res["claimed"], layers=res["layers"],
                R=R, rounds=acc["rounds"], sc_calls=acc["calls"],
                B=B, dmax=acc["dmax"], comm=res["comm"])


# --------------------------------------------------------------------------
# (2a) atomic Schwartz-Zippel: the deg-3 trace zero-test, witness built once.
# --------------------------------------------------------------------------

def _corrupt_constraint(W, q, kind):
    """Apply a trace corruption to a freshly built witness (mirrors the cheat
    branches of run_primitive) so verify_constraints must reject it."""
    Lv, D = W["Lv"], W["D"]
    d = 1 % D
    if kind == "u_consistent":          # wrong quotient, bits kept consistent
        nu = (int(W["u"][d]) + 1)
        W["tabs"]["U"] = W["tabs"]["U"].copy(); W["tabs"]["U"][d] = nu % q
        for j in range(Lv):
            W["tabs"][f"Ub{j}"] = W["tabs"][f"Ub{j}"].copy()
            W["tabs"][f"Ub{j}"][d] = (nu >> (Lv - 1 - j)) & 1
    elif kind == "u_value":             # quotient word off, bits untouched
        W["tabs"]["U"] = W["tabs"]["U"].copy()
        W["tabs"]["U"][d] = (int(W["tabs"]["U"][d]) + 1) % q
    elif kind == "r_value":             # remainder corrupted
        W["tabs"]["R"] = W["tabs"]["R"].copy()
        W["tabs"]["R"][d] = (int(W["tabs"]["R"][d]) + 1) % q
    elif kind == "nonbit":              # a 'bit' set to 2 (range violation)
        W["tabs"]["Ub0"] = W["tabs"]["Ub0"].copy()
        W["tabs"]["Ub0"][0] = 2
    else:
        raise ValueError(kind)


def sz_atomic(n=10, p=3, m=None, dstart=1, kind="u_value", qs=None,
              n_seeds=4000, base_seed=0):
    """Build+corrupt the trace witness ONCE per q, re-run verify_constraints
    with n_seeds fresh challenge sets, return per-q accept rate.  Confirms the
    deg-3 zero-test false-accepts at rate ~ c/q (Schwartz-Zippel)."""
    X = (1 << n) - 1
    if m is None:
        m = max(2, (n + 1) // 2)            # ~ log2 D, D ~ sqrt(X)
    if qs is None:
        qs = [next_prime(X), next_prime(2 * X), next_prime(4 * X),
              next_prime(8 * X)]
    out = []
    for q in qs:
        W = cpmt.build_witness(X, p, m, dstart, q)
        _corrupt_constraint(W, q, kind)
        acc = 0
        stats = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
        for s in range(n_seeds):
            ok = cpmt.verify_constraints(W, np.random.default_rng(base_seed + s),
                                         stats, q)
            acc += int(ok)
        out.append(dict(q=q, rate=acc / n_seeds, accepts=acc, n=n_seeds))
    return dict(n=n, X=X, m=m, kind=kind, rows=out)


# --------------------------------------------------------------------------
# (2b) chain-level composition: false-accept of a real cheat over many seeds.
# --------------------------------------------------------------------------

def sz_chain(n=8, label="corrupt_layer", chain_kw=None, qs=None, n_seeds=1500,
             base_seed=10000, **cfg):
    """Run the REAL run_chain with the lie in `chain_kw` (e.g.
    {'corrupt_layer':1} or {'cheat':'delta_pi'}) over n_seeds fresh challenge
    seeds at each prime q>X; return per-q false-accept rate + the union-bound
    B/q.  Two regimes show up: an output-value lie (delta_pi) is caught
    DETERMINISTICALLY at the final base equality (false-accept ~ 0, no q
    dependence), while a self-consistent INTERNAL lie (corrupt_layer) is caught
    only by Schwartz-Zippel (false-accept ~ c/q) -- the latter sets the field."""
    if chain_kw is None:
        chain_kw = {"corrupt_layer": 1}
    x = (1 << n) - 1
    if qs is None:
        qs = [next_prime(x), next_prime(4 * x)]
    # the degree budget B is field-size-independent (structural); take it once.
    bud = measure_budget(n, q=qs[0], **cfg)
    out = []
    for q in qs:
        acc = 0
        for s in range(n_seeds):
            r = CL.run_chain(x, np.random.default_rng(base_seed + s),
                             q=q, **chain_kw, **cfg)
            acc += int(r["accepted"])
        out.append(dict(q=q, false_accept=acc / n_seeds, accepts=acc,
                        n=n_seeds, union_bound=bud["B"] / q))
    return dict(n=n, x=x, label=label, B=bud["B"], rows=out)


# --------------------------------------------------------------------------
# fitting + extrapolation
# --------------------------------------------------------------------------

def fit_power_log(ns, ys):
    """Fit log2(y) = alpha*n + delta*log2(n) + c  (separates the power-of-x
    exponent alpha from a polylog dressing -- the S531/S532b lesson).  Returns
    (alpha, delta, c, rms_resid)."""
    n = np.asarray(ns, float)
    ly = np.log2(np.asarray(ys, float))
    A = np.column_stack([n, np.log2(n), np.ones_like(n)])
    coef, *_ = np.linalg.lstsq(A, ly, rcond=None)
    resid = ly - A @ coef
    rms = float(np.sqrt(np.mean(resid ** 2)))
    return float(coef[0]), float(coef[1]), float(coef[2]), rms


def fit_simple_exponent(ns, ys):
    """alpha in y ~ (2^n)^alpha = single-power slope of log2(y) vs n."""
    n = np.asarray(ns, float)
    ly = np.log2(np.asarray(ys, float))
    A = np.column_stack([n, np.ones_like(n)])
    coef, *_ = np.linalg.lstsq(A, ly, rcond=None)
    return float(coef[0]), float(coef[1])


def _log2_pi_sqrt_x(n):
    """log2 of K = pi(sqrt(x)) at x=2^n, PNT leading term sqrt(x)/ln(sqrt(x)).
    sqrt(x)=2^(n/2); ln(2^(n/2)) = (n/2)*ln2.  This is the layer count K that
    dominates the soundness budget B -- the structural sqrt(x) backbone."""
    half = n / 2.0
    return half - math.log2(half * math.log(2.0))


def extrapolate(budget_rows, lam=100, target_log10=100, corr_slack_bits=1.0):
    """Derive log2 q for soundness vs correctness, the crossover, and the
    p(10^target_log10) requirement.

    KEY (S511/S515 lesson): a small-n single-power fit of B(x) reads alpha < 0.5
    because B = Theta(K * polylog) with K = pi(sqrt x) = Theta(sqrt x / log x) --
    the 1/log x in the PRIME COUNT itself deflates the finite-n exponent.  So we
    DO NOT extrapolate B by its power fit (which would underestimate B at large
    x); we anchor on the EXACT pi(sqrt x) backbone (measured `layers` == K) and
    a measured polylog factor c_BK = B / K, which grows ~linearly in n (= log x).
    """
    ns = [r["n"] for r in budget_rows]
    Bs = [r["B"] for r in budget_rows]
    Ks = [r["layers"] for r in budget_rows]          # exact pi(sqrt x)
    a, d, c, rms = fit_power_log(ns, Bs)              # reported diagnostic only
    a1, _ = fit_simple_exponent(ns, Bs)
    # c_BK = B/K, fit linear in n (polylog):  c_BK ~ p1*n + p0
    cbk = [B / K for B, K in zip(Bs, Ks)]
    p1, p0 = np.polyfit(np.asarray(ns, float), np.asarray(cbk, float), 1)
    p1, p0 = float(p1), float(p0)

    def log2B(n):                       # log2 K(n) + log2 c_BK(n)
        return _log2_pi_sqrt_x(n) + math.log2(max(1.0, p1 * n + p0))

    def log2q_sound(n):
        return lam + log2B(n)

    def log2q_corr(n):
        return n + corr_slack_bits

    lo, hi = 8.0, 100000.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if log2q_sound(mid) - log2q_corr(mid) > 0:
            lo = mid          # soundness still binds (needs bigger q)
        else:
            hi = mid
    n_cross = 0.5 * (lo + hi)
    x_cross_log10 = n_cross * math.log10(2)

    n_target = target_log10 / math.log10(2)
    return dict(alpha=a, delta=d, c=c, rms=rms, alpha_simple=a1,
                cbk_p1=p1, cbk_p0=p0, cbk_at_target=p1 * n_target + p0,
                lam=lam, corr_slack_bits=corr_slack_bits,
                n_cross=n_cross, x_cross_log10=x_cross_log10,
                n_target=n_target,
                q_sound_target_bits=log2q_sound(n_target),
                q_corr_target_bits=log2q_corr(n_target),
                binds_target=("correctness" if log2q_corr(n_target) >
                              log2q_sound(n_target) else "soundness"),
                log2B_target=log2B(n_target),
                log2K_target=_log2_pi_sqrt_x(n_target))


# --------------------------------------------------------------------------
# correctness slack: smallest prime q for which honest run_chain accepts.
# --------------------------------------------------------------------------

def min_honest_q(n, seed=5, **cfg):
    """Smallest prime q (scanning up from just below X) for which the honest
    chain accepts -- measures the correctness slack over X."""
    x = (1 << n) - 1
    q = next_prime(x - 1)            # first prime >= x (or just below)
    tries = 0
    while tries < 40:
        ok = CL.run_chain(x, np.random.default_rng(seed), q=q, **cfg)["accepted"]
        if ok:
            return q
        q = next_prime(q)
        tries += 1
    return None


# --------------------------------------------------------------------------
# drivers
# --------------------------------------------------------------------------

def run_budget(ns, q=BIG_Q, cfg=None, label="default"):
    cfg = cfg or {}
    print(f"== soundness budget of run_chain ({label} config, field "
          f"{'BIG_Q' if q == BIG_Q else q}) ==")
    print(f"{'n':>3} {'x':>12} {'K=layers':>9} {'R(chal)':>9} {'rounds':>7} "
          f"{'B(budget)':>10} {'dmax':>5} {'B/K':>8} {'accept':>7}")
    rows = []
    for n in ns:
        r = measure_budget(n, q=q, **cfg)
        rows.append(r)
        print(f"{r['n']:>3} {r['x']:>12} {r['layers']:>9} {r['R']:>9} "
              f"{r['rounds']:>7} {r['B']:>10} {r['dmax']:>5} "
              f"{r['B'] / r['layers']:>8.1f} {str(r['accepted']):>7}")
    if len(rows) >= 3:
        a, d, c, rms = fit_power_log([r["n"] for r in rows],
                                     [r["B"] for r in rows])
        a1, _ = fit_simple_exponent([r["n"] for r in rows],
                                    [r["B"] for r in rows])
        aR, dR, cR, _ = fit_power_log([r["n"] for r in rows],
                                      [r["R"] for r in rows])
        print(f"  B(x) 2-term fit: alpha={a:.3f} (power of x), delta={d:+.2f} "
              f"(polylog), rms={rms:.4f};  single-power alpha={a1:.3f}")
        print(f"  R(x) 2-term fit: alpha={aR:.3f}  "
              f"=> R,B ~ x^{a:.2f} * polylog  (sqrt(x) <=> alpha=0.5)")
    return rows


def run_sz_atomic(n=10, kinds=("u_value", "u_consistent", "r_value"),
                  n_seeds=4000):
    print(f"== atomic Schwartz-Zippel: deg-3 trace zero-test, n={n}, "
          f"{n_seeds} seeds/q ==")
    for kind in kinds:
        res = sz_atomic(n=n, kind=kind, n_seeds=n_seeds)
        qs = [r["q"] for r in res["rows"]]
        rates = [r["rate"] for r in res["rows"]]
        print(f"  cheat={kind:12} X={res['X']} m={res['m']}")
        for r in res["rows"]:
            print(f"     q={r['q']:>8}  accept_rate={r['rate']:.5f} "
                  f"({r['accepts']}/{r['n']})  rate*q={r['rate'] * r['q']:.2f}")
        # slope of log(rate) vs log(q): Schwartz-Zippel => -1
        nz = [(q, rt) for q, rt in zip(qs, rates) if rt > 0]
        if len(nz) >= 2:
            lq = np.log([q for q, _ in nz])
            lr = np.log([rt for _, rt in nz])
            slope = float(np.polyfit(lq, lr, 1)[0])
            print(f"     => log(rate) vs log(q) slope = {slope:.2f}  "
                  f"(Schwartz-Zippel predicts -1)")


def run_sz_chain(n=8, n_seeds=1500):
    print(f"== chain composition: real run_chain false-accept, n={n}, "
          f"{n_seeds} seeds/q ==")
    lies = [("delta_pi (output-value lie)", {"cheat": "delta_pi"}),
            ("corrupt_layer=1 (self-consistent internal lie)",
             {"corrupt_layer": 1})]
    qs = [next_prime((1 << n) - 1), next_prime(4 * ((1 << n) - 1)),
          next_prime(16 * ((1 << n) - 1))]
    for label, kw in lies:
        res = sz_chain(n=n, label=label, chain_kw=kw, qs=qs, n_seeds=n_seeds)
        print(f"  {label}:   B={res['B']}")
        for r in res["rows"]:
            ok = r["false_accept"] <= r["union_bound"] + 1e-9
            print(f"     q={r['q']:>8}  false_accept={r['false_accept']:.5f} "
                  f"({r['accepts']}/{r['n']})  rate*q={r['false_accept'] * r['q']:.2f}"
                  f"  union B/q={r['union_bound']:.4f}  within_bound={ok}")


def run_extrapolate(ns=(8, 10, 12, 14, 16), q=BIG_Q, lam=100):
    rows = [measure_budget(n, q=q) for n in ns]
    # measure the correctness slack at a couple of small n
    slacks = []
    for n in (8, 10, 12):
        mq = min_honest_q(n)
        if mq is not None:
            slacks.append(math.log2(mq) - n)
    slack = max(1.0, max(slacks) if slacks else 1.0)
    ex = extrapolate(rows, lam=lam, corr_slack_bits=slack)
    print(f"== field-size extrapolation (lambda={lam} bits) ==")
    print(f"  small-n single-power fit of B: alpha={ex['alpha_simple']:.3f} "
          f"(DEFLATED < 0.5 by the pi(sqrt x)=Theta(sqrt x/log x) 1/log x; "
          f"S511/S515 PNT discount) -- NOT used for extrapolation")
    print(f"  structural: B = c_BK * K, K=pi(sqrt x) exact; "
          f"c_BK = {ex['cbk_p1']:.2f}*n {ex['cbk_p0']:+.1f} (polylog, ~linear in "
          f"log x); c_BK(target)={ex['cbk_at_target']:.0f}")
    print(f"  measured correctness slack over X: ~{slack:.2f} bits "
          f"(min honest log2 q - n)")
    print(f"  soundness curve:   log2 q > {lam} + log2 K(x) + log2 c_BK "
          f"~ {lam} + n/2 - log2 n + polylog")
    print(f"  correctness curve: log2 q > n + {slack:.1f}")
    print(f"  CROSSOVER at n ~ {ex['n_cross']:.0f}  "
          f"(x ~ 10^{ex['x_cross_log10']:.0f})")
    print(f"  at p(10^100): n={ex['n_target']:.0f}, "
          f"log2 K(sqrt x)={ex['log2K_target']:.0f}, log2 B={ex['log2B_target']:.0f};  "
          f"q_soundness ~ {ex['q_sound_target_bits']:.0f} bits, "
          f"q_correctness ~ {ex['q_corr_target_bits']:.0f} bits "
          f"=> {ex['binds_target'].upper()} BINDS")
    return ex


# --------------------------------------------------------------------------
# selftest
# --------------------------------------------------------------------------

def selftest():
    print("== selftest ==")
    # [1] next_prime
    assert next_prime(255) == 257 and next_prime(256) == 257
    assert next_prime(1) == 2 and next_prime(2) == 3
    assert _is_prime(2 ** 31 - 1) and not _is_prime(2 ** 31 - 2)
    print("  [1] prime helper OK")

    # [2] CountRng is a verbatim drop-in: same verdict as a plain rng, and the
    #     counter equals the number of integers() calls.
    x = (1 << 8) - 1
    base = CL.run_chain(x, np.random.default_rng(3), q=257)
    crng = CountRng(3)
    patched_acc = {"calls": 0, "rounds": 0, "budget": 0, "dmax": 0}
    restore = patch_sumcheck(patched_acc)
    try:
        counted = CL.run_chain(x, crng, q=257)
    finally:
        restore()
    assert base["accepted"] == counted["accepted"]
    assert base["claimed"] == counted["claimed"]
    assert base["comm"] == counted["comm"], (base["comm"], counted["comm"])
    assert crng.n_chal > 0 and patched_acc["calls"] > 0
    assert patched_acc["rounds"] <= crng.n_chal       # sumcheck rounds subset
    assert patched_acc["dmax"] >= 2
    print(f"  [2] CountRng + sumcheck patch bit-identical verdict; "
          f"R={crng.n_chal} rounds={patched_acc['rounds']} "
          f"dmax={patched_acc['dmax']} OK")

    # [3] sumcheck restore: after restore, the symbol is the original again.
    assert cpmt.sumcheck.__name__ == "sumcheck"
    r2 = measure_budget(8, q=257)
    assert r2["accepted"] and r2["R"] > 0 and r2["B"] >= r2["rounds"]
    print(f"  [3] restore OK; measure_budget(n=8) R={r2['R']} B={r2['B']}")

    # [4] budget is field-independent (structural): same R,B,dmax over q & BIG_Q.
    a = measure_budget(8, q=257)
    b = measure_budget(8, q=BIG_Q)
    assert (a["R"], a["B"], a["dmax"], a["rounds"]) == \
           (b["R"], b["B"], b["dmax"], b["rounds"]), (a, b)
    print(f"  [4] R,B,dmax field-independent (q=257 == BIG_Q): "
          f"R={a['R']} B={a['B']}")

    # [5] budget grows with n and B/sqrt(x) is bounded (sqrt-x scaling).
    rs = [measure_budget(n, q=BIG_Q) for n in (8, 10, 12)]
    assert rs[0]["B"] < rs[1]["B"] < rs[2]["B"]
    ratios = [r["B"] / math.sqrt(r["x"]) for r in rs]
    # B/sqrt(x) should not blow up (it carries only polylog); allow a generous
    # band for the small-n discrete K=pi(sqrt x) lumps.
    assert max(ratios) / min(ratios) < 6.0, ratios
    print(f"  [5] B grows; B/sqrt(x) bounded ratios={[f'{v:.2f}' for v in ratios]}")

    # [6] atomic Schwartz-Zippel: honest accepts always; corrupted accept-rate
    #     is in (0,1) at q just over X and ROUGHLY halves when q doubles.
    X = (1 << 8) - 1
    q1, q2 = next_prime(X), next_prime(2 * X)
    W = cpmt.build_witness(X, 3, 4, 1, q1)
    st = {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}
    assert all(cpmt.verify_constraints(W, np.random.default_rng(s), st, q1)
               for s in range(50)), "honest constraint must always accept"
    r1 = sz_atomic(n=8, kind="u_value", qs=[q1, q2], n_seeds=3000)
    a1, a2 = r1["rows"][0]["rate"], r1["rows"][1]["rate"]
    assert 0.0 < a1 < 0.5, a1
    assert a2 < a1, (a1, a2)                  # bigger field -> fewer false accepts
    # rate*q roughly constant (Schwartz-Zippel): within 2x across the doubling
    rq1, rq2 = a1 * q1, a2 * q2
    assert 0.4 < rq1 / rq2 < 2.5, (rq1, rq2)
    print(f"  [6] atomic SZ: rate(q={q1})={a1:.4f} rate(q={q2})={a2:.4f}; "
          f"rate*q {rq1:.1f} vs {rq2:.1f} (~const) OK")

    # [7] a nonbit corruption (range violation) is caught with rate ~0 even at
    #     small q (it fails the bit-range part with overwhelming prob).
    Wn = cpmt.build_witness(X, 3, 4, 1, q1)
    _corrupt_constraint(Wn, q1, "nonbit")
    acc = sum(cpmt.verify_constraints(Wn, np.random.default_rng(s), st, q1)
              for s in range(300))
    assert acc <= 300 * 0.2, acc
    print(f"  [7] nonbit corruption rejected (accepts {acc}/300) OK")

    # [8] chain composition, two regimes:
    #   (a) self-consistent INTERNAL lie (corrupt_layer): caught only by
    #       Schwartz-Zippel -> false-accept in (0,1), <= B/q, drops with q.
    #   (b) output-value lie (delta_pi): caught DETERMINISTICALLY at the final
    #       base equality -> false-accept ~ 0 at every q (no field dependence).
    qsc = [next_prime(X), next_prime(8 * X)]
    rc = sz_chain(n=8, label="corrupt_layer", chain_kw={"corrupt_layer": 1},
                  qs=qsc, n_seeds=600)
    for row in rc["rows"]:
        assert row["false_accept"] <= row["union_bound"] + 1e-9, row
    assert 0.0 < rc["rows"][0]["false_accept"] < 0.5, rc["rows"][0]
    assert rc["rows"][1]["false_accept"] <= rc["rows"][0]["false_accept"] + 0.01
    rd = sz_chain(n=8, label="delta_pi", chain_kw={"cheat": "delta_pi"},
                  qs=qsc, n_seeds=300)
    assert all(r["false_accept"] == 0.0 for r in rd["rows"]), rd["rows"]
    print(f"  [8] corrupt_layer within B/q & drops with q "
          f"{[f'{r['false_accept']:.4f}' for r in rc['rows']]}; "
          f"delta_pi caught deterministically (0 false-accepts) OK")

    # [9] extrapolation anchored on pi(sqrt x): B = K*c_BK with K=pi(sqrt x)
    #     exact and c_BK ~ linear in n (= log x), so B = Theta(sqrt x) (the
    #     1/log x in K cancels the log x in c_BK; single-power alpha ~ 0.5).
    #     Crossover near 2*lambda; the p(10^100) target binds CORRECTNESS.
    rows = [measure_budget(n, q=BIG_Q) for n in (8, 10, 12, 14)]
    cbk = [r["B"] / r["layers"] for r in rows]
    assert cbk[0] < cbk[-1], cbk                     # c_BK grows (polylog in x)
    assert cbk[-1] / cbk[0] < (rows[-1]["n"] / rows[0]["n"]) * 1.5, cbk  # ~linear
    ex = extrapolate(rows, lam=100, corr_slack_bits=1.0)
    assert ex["n_cross"] > 100, ex["n_cross"]        # crossover near 2*lambda
    assert ex["binds_target"] == "correctness", ex["binds_target"]
    assert 0.45 < ex["alpha_simple"] < 0.70, ex["alpha_simple"]   # ~ sqrt x
    assert ex["q_corr_target_bits"] > ex["q_sound_target_bits"], ex
    print(f"  [9] c_BK polylog {[f'{v:.0f}' for v in cbk]}; B~sqrt(x) alpha="
          f"{ex['alpha_simple']:.3f}; crossover n~"
          f"{ex['n_cross']:.0f}; p(10^100) q_corr={ex['q_corr_target_bits']:.0f}"
          f">q_sound={ex['q_sound_target_bits']:.0f} -> binds "
          f"{ex['binds_target']} OK")

    print("ALL SELFTESTS PASSED")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--budget", action="store_true",
                    help="measure R(x), B(x), d_max over --ns")
    ap.add_argument("--sz-atomic", action="store_true",
                    help="per-round Schwartz-Zippel law on the trace zero-test")
    ap.add_argument("--sz-chain", action="store_true",
                    help="end-to-end chain composition false-accept audit")
    ap.add_argument("--extrapolate", action="store_true",
                    help="log2 q curves + crossover + p(10^100) requirement")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--ns", type=str, default="8,10,12,14,16")
    ap.add_argument("--seeds", type=int, default=4000)
    ap.add_argument("--chain-seeds", type=int, default=1500)
    ap.add_argument("--lam", type=int, default=100)
    args = ap.parse_args()

    ns = [int(v) for v in args.ns.split(",")]

    if args.selftest:
        selftest()
        return
    did = False
    if args.budget or args.all:
        run_budget(ns, q=BIG_Q, label="default"); print(); did = True
    if args.sz_atomic or args.all:
        run_sz_atomic(n=ns[1] if len(ns) > 1 else ns[0],
                      n_seeds=args.seeds); print(); did = True
    if args.sz_chain or args.all:
        run_sz_chain(n=ns[0], n_seeds=args.chain_seeds); print(); did = True
    if args.extrapolate or args.all:
        run_extrapolate(tuple(ns), lam=args.lam); print(); did = True
    if not did:
        selftest()


if __name__ == "__main__":
    main()
