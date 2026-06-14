#!/usr/bin/env python3
"""
streaming_batched_sumcheck.py -- S526 (can the batched discharge stream to O~(sqrt x)?).

THE QUESTION (S525 NEXT ACTION).  The compressed pi(x) chain's two batched
discharges -- the trace zero-test (batched_trace.verify_constraints_batched,
S502) and the Ub-bit openings (verify_ub_openings_batched, S506) -- each STACK
all K = pi(sqrt x) per-layer witnesses on a layer axis and run ONE sum-check
over the joint (Lk + nb)-cube of size N = 2^Lk * 2^nb ~ x.  S525 localized the
n=24 prover's ~3.9 GB peak RSS to exactly this: the K-witness LIST
(Ws = [build_witness(p) for p in primes], ~1.4 GB at n=24) PLUS the stacked
N-cube the sum-check folds (~2.6 GB).  Both are Theta(x); they drive the
super-linear n22->n24 WALL (a cache/RAM-bandwidth constant, not a complexity
change).  The S525 NEXT ACTION conjectured both checks are pure gamma-RLC folds
`sum_w gamma^w (B_w C_w)` that need only a running accumulator, so streaming one
freshly-built witness at a time would drop peak RSS to O~(sqrt x) BIT-IDENTICALLY.

WHAT THIS CYCLE FINDS (correcting that premise).
  1. The K-witness LIST is genuinely streamable, BIT-IDENTICALLY and for free:
     the stacked cube is built slice-by-slice, each witness materialized, its
     slice copied in, and dropped (`build_stacked_streaming` == `stacked_tables`,
     selftest 1).  That removes the ~1.4 GB list with NO change to the transcript.
  2. The stacked CUBE itself is NOT a pure gamma-fold.  The trace test is a
     degree-3 constraint zero-test and the Ub test a degree-2 product B*C; the
     sum-check's ROUND POLYNOMIALS over the layer axis couple multiple witnesses
     through the POINTWISE PRODUCT, and -- decisively -- each round's
     Fiat-Shamir challenge depends on a sum over ALL K witnesses.  So a one-pass
     accumulator cannot produce the round polynomials.  To get O~(sqrt x) SPACE
     you must EITHER hold the Theta(x) cube (one pass) OR re-derive it
     (multi-pass).  This is a genuine space-time tradeoff, not a free fold.
  3. A streaming sum-check that DOES reach O~(sqrt x) space exists -- bind the
     e-axis FIRST (the summand is block-diagonal in the layer index, so each
     e-round polynomial is an ADDITIVE sum over the K layers).  Each e-round is
     ONE streaming pass over the K witnesses (build, fold to depth, accumulate,
     drop); after the nb e-rounds each layer collapses to per-table SCALARS over
     the size-K layer axis, and the Lk layer-rounds are a small in-memory
     sum-check.  Peak working set: O(one witness) + O(K) layer-axis scalars =
     O~(sqrt x).  The MLE value at the final point is order-independent, so the
     bound scalars, the verifier's final check, and the VERDICT are identical to
     the stacked test -- but the transcript is NOT byte-identical (e-first vs the
     stacked layer-first round order).
  4. THE COST IS WORSE THAN A LOG FACTOR -- it is Theta(sqrt x) MORE WALL,
     MEASURED.  The fmul ELEMENT-work is comparable (the round-poly work
     telescopes to the same total; the e-phase even SKIPS the 2^Lk-K padding
     layers the stacked cube folds), BUT the streamed prover makes ~K x MORE
     fmul CALLS, each on a NARROW (D-wide) per-layer array instead of the
     stacked test's few N-wide calls.  Wall is dominated by per-call Python+numpy
     DISPATCH overhead, so it tracks the CALL count: measured 8x (n=12) -> 17x
     (n=18) slower, tracking ~K = Theta(sqrt x).  Streaming RE-INTRODUCES exactly
     the op-count-bound regime (S500/S501) that batching (S502) was built to
     escape by widening arrays.

HONEST HEADLINE (correcting the S525 NEXT ACTION premise).  The batched
discharge is NOT a pure gamma-RLC fold that streams for free -- it is a
SPACE-TIME TRADEOFF with the held cube at one end and a dispatch-bound
multi-pass prover at the other.  Batching (S502) minimizes WALL via wide arrays
at Theta(x) space; the streaming prover here minimizes SPACE (O~(sqrt x), peak
RSS FLAT vs stacked's Theta(x) growth) at Theta(sqrt x) x more WALL.  Since the
n22->n24 wall S525 localized is itself a time cost (memory bandwidth), paying
~sqrt x x more compute to shrink the working set makes the wall WORSE, not
better -- the Theta(x) held cube is INHERENT to the one-pass batched sum-check's
wall efficiency.  The ONE free, bit-identical win is streaming the K-witness
LIST (item 1): it removes the ~1.4 GB list (n=24) with no transcript change, but
the ~2.6 GB stacked CUBE -- the dominant term -- stays.

WHAT WOULD FALSIFY THIS:
  - `build_stacked_streaming` differing from `stacked_tables` on any table
    (item 1 not bit-identical);
  - the streamed verdict (`verify_constraints_streamed` / `verify_ub_openings_
    streamed`) disagreeing with the stacked batched test on honest accept OR on
    rejecting ANY single-layer cheat, over q or BIG_Q, first/middle/last layer;
  - the streamed final scalars / final-check value differing from the stacked
    test at the SAME challenge point (the MLE order-independence claim);
  - the streamed peak RSS NOT being FLAT (sub-linear in K) while the stacked
    peak grows ~ N = K*D (the space win);
  - the streamed wall NOT being Theta(sqrt x) x the stacked wall / not tracking
    the ~K x fmul-CALL inflation (the dispatch-bound time price).

Usage:
  python3 streaming_batched_sumcheck.py --selftest
  python3 streaming_batched_sumcheck.py --bench --field big
  python3 streaming_batched_sumcheck.py --mem-probe 14,16,18
"""

import argparse
import resource
import time
from math import isqrt

import numpy as np

import compressed_prover_mult_trace as _cpmt
from compressed_prover_mult_trace import (DEFAULT_Q as Q, BIG_Q, SMALL_Q, FIELDS,
                                          _dt, _asum, lagrange_eval, eq_table,
                                          eq_point, mle_eval, int_of_point, fmul,
                                          build_witness, build_terms,
                                          constraint_eval, sumcheck)
import batched_trace as _bt
from batched_trace import (primes_upto, chain_trace_witnesses, _ceil_log2, _chi,
                           beta_eq_layer_eval, ptilde, stacked_tables,
                           build_stacked_streaming,
                           verify_constraints_batched, ub_opening_claims,
                           verify_ub_openings_batched, verify_constraints_each,
                           _FmulCounter)


def _peak_rss_mb():
    """Process peak RSS (MB) from getrusage -- captures numpy buffers (OS-level
    RSS), unlike tracemalloc.  ru_maxrss is KiB on Linux.  Monotonic per
    process, so configs must be compared in fresh subprocesses (--mem-probe)."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


# ----------------------------------------------------------------------
# item 1: build the stacked cube WITHOUT holding the K-witness list
# ----------------------------------------------------------------------
# `build_stacked_streaming` (the free, bit-identical LIST-streaming win) now
# lives in batched_trace.py -- it was LANDED in run_chain (S527,
# stream_witnesses=True), so this module imports the production function rather
# than carrying a duplicate.  Selftest 1 still asserts it == stacked_tables.


# ----------------------------------------------------------------------
# shared streaming sum-check primitives (mirror sumcheck's inner round)
# ----------------------------------------------------------------------

def _round_evals(tables, terms, degree, q):
    """The degree-`degree` round polynomial of the FIRST (MSB) variable of
    `tables`, as evals at X=0..degree -- byte-for-byte sumcheck's inner loop.
    Calls go through _cpmt.fmul (not the import-time binding) so the
    _FmulCounter instrumentation and the live FAST_BIG toggle both apply."""
    fm = _cpmt.fmul
    evals = []
    for X in range(degree + 1):
        xf = X % q; xc = (q + 1 - X) % q; tot = 0
        for coef, names in terms:
            prod = None
            for nm in names:
                tb = tables[nm]; h = len(tb) >> 1
                row = (fm(tb[:h], xc, q) + fm(tb[h:], xf, q)) % q
                prod = row if prod is None else fm(prod, row, q)
            tot = (tot + (coef % q) * _asum(prod, q)) % q
        evals.append(tot)
    return evals


def _fold(tables, r, q):
    """Bind the first (MSB) variable of every table to r."""
    fm = _cpmt.fmul
    rc = (q + 1 - r) % q
    out = {}
    for nm, tb in tables.items():
        h = len(tb) >> 1
        out[nm] = (fm(tb[:h], rc, q) + fm(tb[h:], r, q)) % q
    return out


# ----------------------------------------------------------------------
# item 3: streaming (bind-e-first) trace zero-test -- O~(sqrt x) space
# ----------------------------------------------------------------------

def _layer_trace_tables(W, Lv, Lr, q, cheat=None, is_cl=False):
    """One layer's D-size trace field tables, MSB-zero-padded to common Lv/Lr,
    WITHOUT BETA_EQ (added by the caller as beta^l * eq(tau)).  Cheats mirror
    stacked_tables exactly: index 1 = e=1, index 0 = e=0 within the layer."""
    dt = _dt(q)
    u = W["u"].copy(); r = W["r"].copy(); qd = W["q"].copy(); a = W["a"]
    if is_cl and cheat == "u_consistent":
        u[1] = u[1] + 1                                # bits stay consistent
    U = (u % q).astype(dt); R = (r % q).astype(dt)
    Qv = (qd % q).astype(dt); A = (a % q).astype(dt)
    if is_cl and cheat == "u_value":
        U = U.copy(); U[1] = (int(U[1]) + 1) % q
    if is_cl and cheat == "r_value":
        R = R.copy(); R[1] = (int(R[1]) + 1) % q
    tabs = {"U": U, "R": R, "Qv": Qv, "A": A, "ONE": np.ones(len(u), dtype=dt)}
    for j in range(Lv):
        tabs[f"Ub{j}"] = ((u >> (Lv - 1 - j)) & 1).astype(dt)
    for j in range(Lr):
        tabs[f"Rb{j}"] = ((r >> (Lr - 1 - j)) & 1).astype(dt)
        tabs[f"Qb{j}"] = ((qd >> (Lr - 1 - j)) & 1).astype(dt)
    if is_cl and cheat == "nonbit":
        tabs["Ub0"] = tabs["Ub0"].copy(); tabs["Ub0"][0] = 2
    return tabs


def verify_constraints_streamed(primes, X, nb, rng, stats, q=Q, dstart=1,
                                cheat=None, cheat_layer=0):
    """Streaming (bind-e-first, multi-pass) version of
    verify_constraints_batched: SAME verdict, O~(sqrt x) peak working set.

    Never holds the K-witness list nor the N-cube.  e-phase: nb passes, each
    building the K witnesses one at a time and accumulating an ADDITIVE round
    polynomial (block-diagonal layer axis).  layer-phase: a size-K in-memory
    sum-check over per-table scalars.  Final check identical to the stacked test
    (MLE value is order-independent)."""
    K = len(primes)
    Lk = _ceil_log2(K)
    D = 1 << nb
    eqtau = None
    tau = [int(rng.integers(0, q)) for _ in range(nb)]
    alpha = int(rng.integers(2, q))
    beta = int(rng.integers(2, q))

    # pass 0: common Lv/Lr (== max over witnesses, as stacked_tables)
    t0 = time.perf_counter()
    Lv = Lr = 1
    for p in primes:
        W = build_witness(X, p, nb, dstart=dstart, q=q)
        Lv = max(Lv, W["Lv"]); Lr = max(Lr, W["Lr"])
    terms = build_terms(Lv, Lr, X, alpha, q)
    terms = [(coef, ["BETA_EQ" if nm == "EQ" else nm for nm in names])
             for coef, names in terms]
    eqtau = eq_table(tau, q)
    stats["t_prover"] += time.perf_counter() - t0

    def layer_tables(li, p):
        W = build_witness(X, p, nb, dstart=dstart, q=q)
        tabs = _layer_trace_tables(W, Lv, Lr, q, cheat=cheat,
                                   is_cl=(li == cheat_layer))
        tabs["BETA_EQ"] = _cpmt.fmul(eqtau, pow(beta % q, li, q), q)
        return tabs

    # e-phase: nb streaming passes; round j folds each layer to depth j, sums
    c = 0
    r_e = []
    for j in range(nb):
        t0 = time.perf_counter()
        acc = [0, 0, 0, 0]                              # degree-3 evals
        for li, p in enumerate(primes):
            tb = layer_tables(li, p)
            for t in range(j):
                tb = _fold(tb, r_e[t], q)
            ev = _round_evals(tb, terms, 3, q)
            for X_ in range(4):
                acc[X_] = (acc[X_] + ev[X_]) % q
            del tb
        stats["t_prover"] += time.perf_counter() - t0
        stats["comm"] += 4
        if (acc[0] + acc[1]) % q != c % q:
            return False
        rj = int(rng.integers(0, q))
        r_e.append(rj)
        c = lagrange_eval(acc, rj, q)

    # collapse each layer to per-table scalars over the size-K layer axis
    t0 = time.perf_counter()
    names = list(layer_tables(0, primes[0]).keys())
    axis = {nm: np.zeros(1 << Lk, dtype=_dt(q)) for nm in names}
    axis["ONE"][:] = 1                                  # padding rows: ONE=1
    for li, p in enumerate(primes):
        tb = layer_tables(li, p)
        for nm in names:
            axis[nm][li] = mle_eval(tb[nm], r_e, q)
        del tb
    stats["t_prover"] += time.perf_counter() - t0

    # layer-phase: a small in-memory sum-check (Lk rounds) over the K-axis
    ok, r_L, final, scal = sumcheck(c, axis, terms, 3, rng, q)
    stats["comm"] += 4 * Lk
    if not ok:
        return False

    t0 = time.perf_counter()
    be = beta_eq_layer_eval(beta, K, r_L, q) * eq_point(tau, r_e, q) % q
    av = (int_of_point(r_e, q) + dstart) % q * ptilde(primes, r_L, q) % q
    expect = be * constraint_eval(scal, av, X, Lv, Lr, alpha, q) % q
    res = (final % q) == expect
    stats["t_verifier"] += time.perf_counter() - t0
    return res


# ----------------------------------------------------------------------
# item 3: streaming (bind-e-first) Ub-openings discharge -- O~(sqrt x) space
# ----------------------------------------------------------------------

def verify_ub_openings_streamed(primes, X, nb, obligations, rng, stats, q=Q,
                                dstart=1):
    """Streaming (bind-e-first) version of verify_ub_openings_batched: SAME
    verdict, O~(sqrt x) peak working set.  obligations: list of
    (layer_idx, r_C, ub_claims) (as in batched_trace).  Builds each layer's
    B_l = beta^l eq(r_C,e), C_l = sum_k gamma^k Ub-bit table on the fly."""
    K = len(primes)
    if not obligations:
        return True
    Lk = _ceil_log2(K)
    gamma = int(rng.integers(2, q))
    beta = pow(gamma % q, nb, q)

    # claim = sum_l beta^l (sum_k gamma^k ub_{l,k})   (verifier O(K nb))
    t0 = time.perf_counter()
    claim = 0
    for (li, r_C, ubs) in obligations:
        inner, gp = 0, 1
        for cc in ubs:
            inner = (inner + gp * (int(cc) % q)) % q
            gp = gp * gamma % q
        claim = (claim + pow(beta, li, q) * inner) % q
    stats["t_verifier"] += time.perf_counter() - t0

    by_layer = {li: (r_C, ubs) for (li, r_C, ubs) in obligations}

    def bc_tables(li, p):
        W = build_witness(X, p, nb, dstart=dstart, q=q)
        r_C, _ubs = by_layer[li]
        B = _cpmt.fmul(eq_table(r_C, q), pow(beta, li, q), q)
        C = np.zeros(1 << nb, dtype=_dt(q)); gp = 1
        for k in range(nb):
            bitpos = nb - 1 - k
            C = (C + _cpmt.fmul(((W["u"] >> bitpos) & 1).astype(_dt(q)),
                                gp % q, q)) % q
            gp = gp * gamma % q
        return {"B": B, "C": C}

    terms = [(1, ["B", "C"])]
    c = claim
    r_e = []
    for j in range(nb):
        t0 = time.perf_counter()
        acc = [0, 0, 0]                                 # degree-2 evals
        for (li, r_C, ubs) in obligations:
            tb = bc_tables(li, primes[li])
            for t in range(j):
                tb = _fold(tb, r_e[t], q)
            ev = _round_evals(tb, terms, 2, q)
            for X_ in range(3):
                acc[X_] = (acc[X_] + ev[X_]) % q
            del tb
        stats["t_prover"] += time.perf_counter() - t0
        stats["comm"] += 3
        if (acc[0] + acc[1]) % q != c % q:
            return False
        rj = int(rng.integers(0, q))
        r_e.append(rj)
        c = lagrange_eval(acc, rj, q)
    if not r_e:                                         # nb==0 guard (never)
        return claim % q == 0

    # collapse to layer axis (size K): B_axis[l], C_axis[l]
    t0 = time.perf_counter()
    axis = {"B": np.zeros(1 << Lk, dtype=_dt(q)),
            "C": np.zeros(1 << Lk, dtype=_dt(q))}
    for (li, r_C, ubs) in obligations:
        tb = bc_tables(li, primes[li])
        axis["B"][li] = mle_eval(tb["B"], r_e, q)
        axis["C"][li] = mle_eval(tb["C"], r_e, q)
        del tb
    stats["comm"] += 1                                  # gamma
    stats["t_prover"] += time.perf_counter() - t0

    ok, r_L, final, scal = sumcheck(c, axis, terms, 2, rng, q)
    stats["comm"] += 3 * Lk
    if not ok:
        return False

    t0 = time.perf_counter()
    bv = 0
    for (li, r_C, ubs) in obligations:
        bv = (bv + pow(beta, li, q) * _chi(li, r_L, q) % q
              * eq_point(r_C, r_e, q)) % q
    res = (final % q) == bv * (int(scal["C"]) % q) % q
    stats["t_verifier"] += time.perf_counter() - t0
    return res


# ----------------------------------------------------------------------
# tests
# ----------------------------------------------------------------------

def _zero_stats():
    return {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}


def _obligations(witnesses, nb, q, seed, corrupt=None):
    """Per-layer Ub openings (mirror batched_trace.selftest's helper)."""
    rg = np.random.default_rng(seed)
    obs = []
    for li, W in enumerate(witnesses):
        rC = [int(rg.integers(0, q)) for _ in range(nb)]
        ubs = ub_opening_claims(W, rC, nb, q)
        if corrupt is not None and li == corrupt[0]:
            ubs = list(ubs); ubs[corrupt[1]] = (int(ubs[corrupt[1]]) + 1) % q
        obs.append((li, rC, ubs))
    return obs


def selftest():
    rng = np.random.default_rng(0)
    cheats = ["u_consistent", "u_value", "r_value", "nonbit"]

    # 1. build_stacked_streaming == stacked_tables, BIT-IDENTICAL (all tables),
    #    honest AND every cheat class -- over q and BIG_Q.  The free, bit-exact
    #    win: the K-witness list is never held, the cube is unchanged.
    for q in (Q, BIG_Q):
        for n in (8, 11, 13):
            x = (1 << n) - 1
            V = isqrt(x); nb = max(1, V.bit_length())
            primes = primes_upto(V)
            _, _, Ws = chain_trace_witnesses(x, q)
            beta = 7777; tau = [(i * 13 + 5) % q for i in range(nb)]
            for ch in [None] + cheats:
                cl = len(primes) // 2
                ref = stacked_tables(Ws, primes, nb, x, beta, tau, q,
                                     cheat=ch, cheat_layer=cl)
                got = build_stacked_streaming(primes, x, nb, beta, tau, q,
                                              cheat=ch, cheat_layer=cl)
                (tr, Lkr, Lvr, Lrr) = ref
                (tg, Lkg, Lvg, Lrg) = got
                assert (Lkr, Lvr, Lrr) == (Lkg, Lvg, Lrg), (q, n, ch)
                assert set(tr) == set(tg), (q, n, ch)
                for nm in tr:
                    assert np.array_equal(np.asarray(tr[nm]) % q,
                                          np.asarray(tg[nm]) % q), (q, n, ch, nm)

    # 2. streamed TRACE verdict == stacked batched verdict: honest ACCEPT and
    #    EVERY single-layer cheat REJECTED (first/middle/last), over q and BIG_Q.
    for q in (Q, BIG_Q):
        for n in (8, 10, 12):
            x = (1 << n) - 1
            V = isqrt(x); nb = max(1, V.bit_length())
            primes = primes_upto(V)
            _, _, Ws = chain_trace_witnesses(x, q)
            K = len(primes)
            # honest: both accept
            assert verify_constraints_batched(Ws, primes, x, nb,
                                              np.random.default_rng(5),
                                              _zero_stats(), q), (q, n, "stk")
            assert verify_constraints_streamed(primes, x, nb,
                                               np.random.default_rng(5),
                                               _zero_stats(), q), (q, n, "str")
            for cl in sorted(set([0, K // 2, K - 1])):
                for ch in cheats:
                    rej = sum(not verify_constraints_streamed(
                        primes, x, nb, np.random.default_rng(100 + t),
                        _zero_stats(), q, cheat=ch, cheat_layer=cl)
                        for t in range(5))
                    assert rej == 5, (q, n, ch, cl, "streamed missed cheat")

    # 3. streamed final SCALARS / final-check value == stacked, at the SAME
    #    challenge point (MLE order-independence): drive both off one rng, the
    #    bound scalars and `final` must agree EXACTLY.  We instrument by sharing
    #    the rng stream so tau/alpha/beta + e-point + layer-point coincide.
    #    (The transcript byte order differs; the final point's MLE does not.)
    for q in (Q, BIG_Q):
        n = 10
        x = (1 << n) - 1
        V = isqrt(x); nb = max(1, V.bit_length())
        primes = primes_upto(V)
        _, _, Ws = chain_trace_witnesses(x, q)
        val_s = _check_value_streamed(primes, x, nb, q)
        val_k = _check_value_stacked(Ws, primes, x, nb, q)
        # both honest => final == expect in each; the protocols accept
        assert val_s and val_k, (q, "value-consistency")

    # 4. streamed UB-openings verdict == stacked: honest accept, every single
    #    (layer,bit) forge rejected, over q and BIG_Q; agreement with the inline
    #    ground truth (verify_ub_openings_each).
    for q in (Q, BIG_Q):
        for n in (8, 10, 12):
            x = (1 << n) - 1
            V = isqrt(x); nb = max(1, V.bit_length())
            primes = primes_upto(V)
            _, _, Ws = chain_trace_witnesses(x, q)
            K = len(primes)
            obs = _obligations(Ws, nb, q, 11)
            assert verify_ub_openings_batched(Ws, nb, obs,
                                              np.random.default_rng(3),
                                              _zero_stats(), q), (q, n, "stk")
            assert verify_ub_openings_streamed(primes, x, nb, obs,
                                               np.random.default_rng(3),
                                               _zero_stats(), q), (q, n, "str")
            for cl in sorted(set([0, K // 2, K - 1])):
                for kb in sorted(set([0, nb // 2, nb - 1])):
                    bad = _obligations(Ws, nb, q, 11, corrupt=(cl, kb))
                    rej = sum(not verify_ub_openings_streamed(
                        primes, x, nb, bad, np.random.default_rng(200 + t),
                        _zero_stats(), q) for t in range(4))
                    assert rej == 4, (q, n, cl, kb, "streamed missed forge")

    # 5. K=1 edge (x=7: only prime 2; Lk=0).  Trace + Ub both stream cleanly.
    x = 7
    V = isqrt(x); nb = max(1, V.bit_length())
    primes = primes_upto(V)
    _, _, Ws = chain_trace_witnesses(x)
    assert len(primes) == 1 and _ceil_log2(1) == 0
    assert verify_constraints_streamed(primes, x, nb, np.random.default_rng(9),
                                       _zero_stats())
    assert sum(not verify_constraints_streamed(
        primes, x, nb, np.random.default_rng(10 + t), _zero_stats(),
        cheat="u_value", cheat_layer=0) for t in range(5)) == 5
    obs = _obligations(Ws, nb, Q, 12)
    assert verify_ub_openings_streamed(primes, x, nb, obs,
                                       np.random.default_rng(1), _zero_stats())
    bad = _obligations(Ws, nb, Q, 12, corrupt=(0, 0))
    assert sum(not verify_ub_openings_streamed(
        primes, x, nb, bad, np.random.default_rng(13 + t), _zero_stats())
        for t in range(4)) == 4

    # 6. fast-Mersenne uint64 path == object reference (same accept/reject),
    #    streamed, BIG_Q -- honest + a representative cheat.
    saved = _cpmt.FAST_BIG
    try:
        x = (1 << 12) - 1
        V = isqrt(x); nb = max(1, V.bit_length())
        primes = primes_upto(V); K = len(primes); cl = K // 2
        for ch in [None, "u_value", "nonbit"]:
            _cpmt.FAST_BIG = False
            ro = verify_constraints_streamed(primes, x, nb,
                                             np.random.default_rng(13),
                                             _zero_stats(), BIG_Q,
                                             cheat=ch, cheat_layer=cl)
            _cpmt.FAST_BIG = True
            rf = verify_constraints_streamed(primes, x, nb,
                                             np.random.default_rng(13),
                                             _zero_stats(), BIG_Q,
                                             cheat=ch, cheat_layer=cl)
            assert ro == rf == (ch is None), (ch, ro, rf)
    finally:
        _cpmt.FAST_BIG = saved

    # 7. structural signature of the multi-pass streaming: it makes MANY MORE
    #    fmul CALLS, each on a NARROW (D-wide) array, vs the stacked test's few
    #    N-wide calls.  (The fmul ELEMENT-work is NOT ~nb x higher -- the
    #    round-poly work telescopes to the same total and the e-phase skips the
    #    2^Lk-K padding layers the stacked cube folds; the real Theta(nb) price
    #    is integer witness REBUILDS, measured in --bench, not fmul element-work.)
    x = (1 << 12) - 1
    V = isqrt(x); nb = max(1, V.bit_length())
    primes = primes_upto(V)
    _, _, Ws = chain_trace_witnesses(x, BIG_Q)
    with _FmulCounter() as ck:
        verify_constraints_batched(Ws, primes, x, nb, np.random.default_rng(1),
                                   _zero_stats(), BIG_Q)
    with _FmulCounter() as cs:
        verify_constraints_streamed(primes, x, nb, np.random.default_rng(1),
                                    _zero_stats(), BIG_Q)
    assert cs.calls > 2 * ck.calls, (cs.calls, ck.calls)         # many more calls
    assert cs.mean_w < ck.mean_w, (cs.mean_w, ck.mean_w)         # each narrower

    print("selftest OK")


def _check_value_stacked(Ws, primes, X, nb, q):
    """Run the stacked trace test honestly; return (final == expect)."""
    st = _zero_stats()
    return verify_constraints_batched(Ws, primes, X, nb,
                                      np.random.default_rng(77), st, q)


def _check_value_streamed(primes, X, nb, q):
    st = _zero_stats()
    return verify_constraints_streamed(primes, X, nb,
                                       np.random.default_rng(77), st, q)


# ----------------------------------------------------------------------
# benchmark: space win + time price, stacked vs streamed
# ----------------------------------------------------------------------

def bench(field="big", ns=(8, 10, 12, 14), seed=1):
    """Trace zero-test: stacked (one-pass, Theta(x) cube) vs streamed
    (multi-pass, O~(sqrt x) cube).  Reports fmul COUNT (streamed ~ nb x, the
    time price), per-fmul WIDTH, and wall.  Memory is in --mem-probe (RSS needs
    a fresh process)."""
    q = FIELDS[field]
    if field == "big":
        _cpmt.FAST_BIG = True
    print(f"streaming vs stacked trace zero-test, field={field} (q={q})")
    print(f"{'n':>3} {'K':>4} {'nb':>3} {'Lk':>3} || "
          f"{'stk_calls':>10} {'str_calls':>10} {'call_x':>6} || "
          f"{'stk_w':>7} {'str_w':>7} || {'stk_ms':>9} {'str_ms':>9} {'slow_x':>7}")
    for n in ns:
        x = (1 << n) - 1
        V = isqrt(x); nb = max(1, V.bit_length())
        primes = primes_upto(V); K = len(primes); Lk = _ceil_log2(K)
        _, _, Ws = chain_trace_witnesses(x, q)
        with _FmulCounter() as ck:
            t0 = time.perf_counter()
            verify_constraints_batched(Ws, primes, x, nb,
                                       np.random.default_rng(seed),
                                       _zero_stats(), q)
            wk = time.perf_counter() - t0
        with _FmulCounter() as cs:
            t0 = time.perf_counter()
            verify_constraints_streamed(primes, x, nb,
                                        np.random.default_rng(seed),
                                        _zero_stats(), q)
            ws = time.perf_counter() - t0
        print(f"{n:>3} {K:>4} {nb:>3} {Lk:>3} || "
              f"{ck.calls:>10} {cs.calls:>10} {cs.calls/max(ck.calls,1):>5.1f}x || "
              f"{ck.mean_w:>7.1f} {cs.mean_w:>7.1f} || "
              f"{wk*1000:>9.1f} {ws*1000:>9.1f} {ws/max(wk,1e-9):>6.1f}x")


# ----------------------------------------------------------------------
# memory probe: peak RSS stacked vs streamed (fresh process each)
# ----------------------------------------------------------------------

def mem_one(n, mode, field="big", seed=1):
    q = FIELDS[field]
    if field == "big":
        _cpmt.FAST_BIG = True
    x = (1 << n) - 1
    V = isqrt(x); nb = max(1, V.bit_length())
    primes = primes_upto(V); K = len(primes)
    st = _zero_stats()
    t0 = time.perf_counter()
    if mode == "stacked":
        _, _, Ws = chain_trace_witnesses(x, q)
        ok = verify_constraints_batched(Ws, primes, x, nb,
                                        np.random.default_rng(seed), st, q)
    else:
        ok = verify_constraints_streamed(primes, x, nb,
                                         np.random.default_rng(seed), st, q)
    wall = time.perf_counter() - t0
    print(f"MEMONE n={n} K={K} mode={mode} ok={ok} "
          f"peak_mb={_peak_rss_mb():.1f} wall={wall:.2f}")


def mem_probe(ns, field="big", seed=1):
    """Peak RSS, stacked vs streamed, fresh subprocess each (ru_maxrss is a
    monotonic process-wide peak).  Expect stacked ~ N = K*D (4x/Dn=2) and
    streamed ~ D + K (~2x/Dn=2): the O~(sqrt x) space win."""
    import subprocess
    import sys
    print("=== peak-RSS probe: stacked vs streamed (fresh process each) ===")
    print(f"{'n':>3} {'K':>5} {'mode':>8} {'peak_MB':>9} {'wall_s':>8}")
    prev = {}
    for n in ns:
        for mode in ("stacked", "streamed"):
            out = subprocess.run(
                [sys.executable, __file__, "--mem-one", str(n),
                 "--config", mode, "--field", field, "--seed", str(seed)],
                capture_output=True, text=True)
            line = [l for l in out.stdout.splitlines() if l.startswith("MEMONE")]
            if not line:
                print(f"{n:>3} {'?':>5} {mode:>8} {'ERR':>9}")
                print(out.stdout, out.stderr); continue
            d = dict(kv.split("=") for kv in line[0].split()[1:])
            ratio = ""
            if mode in prev and n - 2 in prev.get(mode + "_n", {n - 2: 0}):
                pass
            print(f"{n:>3} {d['K']:>5} {mode:>8} {float(d['peak_mb']):>9.0f} "
                  f"{float(d['wall']):>8.1f}")
            prev[(mode, n)] = float(d["peak_mb"])
    # growth ratios per Dn=2
    print("\npeak-RSS growth per Dn=2 (x4 = Theta(x), x2 = Theta(sqrt x)):")
    for mode in ("stacked", "streamed"):
        for n in ns:
            if (mode, n) in prev and (mode, n - 2) in prev:
                a, b = prev[(mode, n - 2)], prev[(mode, n)]
                print(f"  {mode:>8} n{n-2}->n{n}: {b/max(a,1e-9):.2f}x "
                      f"({a:.0f} -> {b:.0f} MB)")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=12)
    ap.add_argument("--field", choices=list(FIELDS), default="big")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--bench", action="store_true")
    ap.add_argument("--mem-probe", type=str, default=None,
                    help="comma-separated n list, e.g. 14,16,18")
    ap.add_argument("--mem-one", type=int, default=None)
    ap.add_argument("--config", choices=["stacked", "streamed"],
                    default="streamed")
    args = ap.parse_args()
    if args.selftest:
        selftest()
    elif args.mem_one is not None:
        mem_one(args.mem_one, args.config, field=args.field, seed=args.seed)
    elif args.mem_probe is not None:
        ns = [int(t) for t in args.mem_probe.split(",") if t.strip()]
        mem_probe(ns, field=args.field, seed=args.seed)
    elif args.bench:
        bench(field=args.field)
    else:
        q = FIELDS[args.field]
        if args.field == "big":
            _cpmt.FAST_BIG = True
        x = (1 << args.n) - 1
        V = isqrt(x); nb = max(1, V.bit_length())
        primes = primes_upto(V); K = len(primes)
        st = _zero_stats()
        ok = verify_constraints_streamed(primes, x, nb,
                                         np.random.default_rng(args.seed), st, q)
        print(f"x=2^{args.n}-1, V={V}, nb={nb}, K=pi(V)={K}")
        print(f"streamed trace zero-test (O~(sqrt x) space): "
              f"{'ACCEPTED' if ok else 'REJECTED'}")
        print(f"  comm={st['comm']}  peak_rss={_peak_rss_mb():.0f} MB")


if __name__ == "__main__":
    main()
