#!/usr/bin/env python3
"""
leaf_open.py -- S505, the multilinear-MLE polynomial-commitment OPENING primitive.

THE PROBLEM IT SOLVES.  Every protocol in this directory ultimately bottoms out
in a "leaf opening": the verifier holds a claim "S~(pt) == claimed" about a
committed O(sqrt x)-size table S, and discharges it by a DIRECT MLE evaluation
mle_eval(S, pt) -- an O(2^nb) = O(sqrt x) loop.  In the K = pi(sqrt x)-layer
compressed chain (compressed_layer.run_chain) those direct evaluations are the
LAST sqrt x-/p-linear term left in the verifier: with K layers each doing O(1)
leaf opens, the leaf-opening cost alone is O(K sqrt x) ~ Õ(x) -- the only thing
keeping the per-layer verifier from being honestly polylog.

WHAT THIS MODULE BUILDS.  The standard "evaluate a committed multilinear at a
point via ONE sum-check against the eq(point, .) table" reduction:

    S~(pt) = sum_w S[w] * eq~(pt, w)          (the MLE evaluation identity)

is itself a degree-2 sum-check claim.  open_eval runs that sum-check.  Its effect
is a REDUCTION, not a free lunch: it converts the direct claim "S~(pt)=claimed"
into (a) a polylog sum-check transcript [verifier O(nb), comm O(nb)] PLUS (b) a
single RESIDUAL claim "S~(r)=residual" at a FRESH random point r.  The O(2^nb)
work moves entirely to the prover's round polynomials and to that one residual.

The residual is the whole point.  In ISOLATION the caller must still discharge it
(mle_eval -- which is why the standalone is, by itself, NO asymptotic win: you've
moved a sqrt x eval from pt to r).  But the residual is exactly a claim about the
SAME committed table at a NEW point, i.e. precisely the shape of a chain claim --
so in compressed_layer.run_chain(pcs=True) it is THREADED as the carried claim and
discharged transitively at the chain's base (the S_0 closes), never per-layer.
That is the "sum-check-delegated eq-fold whose own leaf is the next layer's claim"
of the program's NEXT ACTION, and it is UNCONDITIONAL (no commitment, no crypto):
soundness rides the existing GKR layer sum-checks down to the closed base.

open_batch folds k claims about ONE table into a SINGLE residual via a random-
linear-combination of eq tables (powers of one gamma) + one degree-2 sum-check.
It is the drop-in replacement for line_batch_pair / batch_on_table.  Costs of the
folding it replaces (chaining k claims by (k-1) line restrictions):
  - VERIFIER: (k-1) closing mle_evals = O(k * 2^nb)   -> open_batch O(k nb) (NO 2^nb)
  - PROVER:   (k-1)(nb+1) line-point evals = O(k nb 2^nb) -> open_batch O(k 2^nb)
so open_batch sheds a 2^nb factor on the verifier and an nb factor on the prover,
emitting ONE residual instead of closing.

SOUNDNESS (Schwartz-Zippel / LFKN sum-check).  If claimed != S~(pt) then the
sum-check's round-1 identity g(0)+g(1)=claimed fails (the honest summand sums to
the TRUE value), so a wrong claim is rejected outright; a prover who instead lies
in the round polynomials is caught with probability >= 1 - 2*nb/q (degree 2, nb
rounds).  For open_batch a wrong c_i survives the gamma-combine with probability
<= (k-1)/q.  The residual scal["S"] is the prover's IMPLIED value S~(r); it is
NOT trusted here -- the caller discharges it (a true close, or the chain base).
The eq scalar IS grounded: the verifier recomputes eq~(pt, r) in O(nb) and checks
the prover's folded eq-table value against it, so the prover cannot fake the eq
side of the product.

WHAT WOULD FALSIFY THIS.  open_eval/open_batch disagreeing with mle_eval on an
honest prover; accepting a wrong claimed value or a wrong batched claim above the
field bound; the residual not equaling the true S~(r) for an honest prover; the
verifier work (excluding the prover's table fold) scaling with 2^nb instead of nb;
or run_chain(pcs=True) changing the chain verdict (honest accept / cheat reject)
vs pcs=False.

Usage:
  python3 leaf_open.py --selftest
  python3 leaf_open.py --bench
"""

import argparse
import time

import numpy as np

from compressed_prover_mult_trace import (DEFAULT_Q as Q, BIG_Q, SMALL_Q,
                                           FIELDS, _dt, _asum, lagrange_eval,
                                           eq_table, eq_point, mle_eval, sumcheck,
                                           fmul)


def _stats():
    return {"t_prover": 0.0, "t_verifier": 0.0, "comm": 0}


def open_eval(S, pt, claimed, rng, stats=None, q=Q):
    """Reduce the opening claim  S~(pt) == claimed  to a single RESIDUAL opening
    S~(r) == residual  at a fresh random point r, via ONE degree-2 sum-check of
    claimed = sum_w S[w] eq~(pt, w).

    Verifier work O(n) (n = len(pt)) + the sum-check round checks O(n); NO O(2^n)
    table evaluation on the verifier side.  The O(2^n) cost lives in the prover's
    eq table + round polynomials and in the returned residual (to be discharged by
    the caller -- a downstream chain claim, the chain base, or a real commitment).

    Returns dict(ok, r, residual)."""
    if stats is None:
        stats = _stats()
    n = len(pt)
    S = np.asarray(S)
    t0 = time.perf_counter()
    EQ = eq_table(pt, q)                                   # prover O(2^n)
    stats["t_prover"] += time.perf_counter() - t0
    ok, r, final, scal = sumcheck(claimed, {"S": S.copy(), "EQ": EQ},
                                  [(1, ["S", "EQ"])], 2, rng, q)
    stats["comm"] += 3 * n                                 # n rounds x (deg+1)=3
    if not ok:
        return dict(ok=False, r=r, residual=None)
    t0 = time.perf_counter()
    eqr = eq_point(pt, r, q)                                # verifier O(n)
    ok = (int(scal["EQ"]) % q == eqr                       # eq grounded to public pt
          and (final % q) == int(scal["S"]) * eqr % q)     # sum-check final binding
    stats["t_verifier"] += time.perf_counter() - t0
    return dict(ok=ok, r=r, residual=int(scal["S"]) % q)


def open_batch(S, claims, rng, stats=None, q=Q):
    """Fold k opening claims {(pt_i, c_i): S~(pt_i)==c_i} on the SAME table S into
    ONE residual  S~(r) == residual  via a random-linear-combination eq-fold +
    one degree-2 sum-check.  gamma <- F_q public coin; combined claim sum_i
    gamma^i c_i over the combined eq table EQ_comb[w] = sum_i gamma^i eq~(pt_i, w).

    Verifier O(k n) (no 2^n term); the line-restriction folding it replaces pays
    (k-1) closing mle_evals = O(k 2^n) on the verifier and O(k n 2^n) on the
    prover.  Returns dict(ok, r, residual)."""
    if stats is None:
        stats = _stats()
    if len(claims) == 1:
        return open_eval(S, claims[0][0], claims[0][1], rng, stats, q)
    n = len(claims[0][0])
    S = np.asarray(S)
    gamma = int(rng.integers(0, q))
    # combined claim (verifier O(k))
    t0 = time.perf_counter()
    comb, pw = 0, 1
    for (_, c) in claims:
        comb = (comb + pw * (int(c) % q)) % q
        pw = pw * gamma % q
    stats["t_verifier"] += time.perf_counter() - t0
    # combined eq table (prover O(k 2^n))
    t0 = time.perf_counter()
    EQ = np.zeros(1 << n, dtype=_dt(q))
    pw = 1
    for (pt, _) in claims:
        EQ = (EQ + fmul(eq_table(pt, q), pw % q, q)) % q
        pw = pw * gamma % q
    stats["t_prover"] += time.perf_counter() - t0
    ok, r, final, scal = sumcheck(comb, {"S": S.copy(), "EQ": EQ},
                                  [(1, ["S", "EQ"])], 2, rng, q)
    stats["comm"] += 1 + 3 * n                             # gamma + n rounds x 3
    if not ok:
        return dict(ok=False, r=r, residual=None)
    # verifier recomputes EQ_comb~(r) = sum_i gamma^i eq~(pt_i, r)  (O(k n))
    t0 = time.perf_counter()
    eqr, pw = 0, 1
    for (pt, _) in claims:
        eqr = (eqr + pw * eq_point(pt, r, q)) % q
        pw = pw * gamma % q
    ok = (int(scal["EQ"]) % q == eqr
          and (final % q) == int(scal["S"]) * eqr % q)
    stats["t_verifier"] += time.perf_counter() - t0
    return dict(ok=ok, r=r, residual=int(scal["S"]) % q)


def close_eval(S, pt, claimed, q=Q):
    """The trivial O(2^n) ground-truth close: mle_eval(S, pt) == claimed.  Used at
    the chain BASE (S_0, a one-time O(sqrt x) within budget) and in tests to
    confirm a residual / opening is the true MLE value.  This is the operation a
    real commitment opening would replace; open_eval REDUCES a non-base claim to
    a residual instead of paying this per layer."""
    return int(mle_eval(np.asarray(S), pt, q)) % q == int(claimed) % q


# ----------------------------------------------------------------------
# tests
# ----------------------------------------------------------------------

def selftest():
    rng = np.random.default_rng(0)

    for q in (Q, BIG_Q, SMALL_Q):
        # 1. open_eval AGREES with mle_eval on an honest prover, and its residual
        #    is the TRUE S~(r).  Across n and random tables.
        for n in (1, 2, 4, 6, 8):
            for _ in range(6):
                S = (rng.integers(0, q, size=1 << n)).astype(_dt(q))
                pt = [int(rng.integers(0, q)) for _ in range(n)]
                true = mle_eval(S, pt, q)
                st = _stats()
                res = open_eval(S, pt, true, rng, st, q)
                assert res["ok"], (q, n, "honest open rejected")
                assert close_eval(S, res["r"], res["residual"], q), \
                    (q, n, "residual != true S~(r)")
                # the residual at r equals the direct MLE there
                assert res["residual"] == mle_eval(S, res["r"], q), (q, n)
                # verifier comm is O(n), NOT O(2^n)
                assert st["comm"] == 3 * n, (q, n, st["comm"])

            # 2. a WRONG claimed value is rejected (round-1 identity fails)
            for _ in range(6):
                S = (rng.integers(0, q, size=1 << n)).astype(_dt(q))
                pt = [int(rng.integers(0, q)) for _ in range(n)]
                bad = (int(mle_eval(S, pt, q)) + 1) % q
                assert not open_eval(S, pt, bad, rng, _stats(), q)["ok"], \
                    (q, n, "wrong claim accepted")

        # 3. open_batch folds k claims into ONE residual, agrees with mle_eval,
        #    residual is the true S~(r); verifier comm O(n) (+1 for gamma),
        #    independent of k.
        for n in (2, 4, 6):
            S = (rng.integers(0, q, size=1 << n)).astype(_dt(q))
            for k in (1, 2, 3, 5):
                pts = [[int(rng.integers(0, q)) for _ in range(n)]
                       for _ in range(k)]
                claims = [(pt, mle_eval(S, pt, q)) for pt in pts]
                st = _stats()
                res = open_batch(S, claims, rng, st, q)
                assert res["ok"], (q, n, k, "honest batch rejected")
                assert close_eval(S, res["r"], res["residual"], q), (q, n, k)
                if k > 1:
                    assert st["comm"] == 1 + 3 * n, (q, n, k, st["comm"])

                # 4. corrupting ANY single claim value is rejected
                for j in range(k):
                    bad = list(claims)
                    bad[j] = (bad[j][0], (int(bad[j][1]) + 1) % q)
                    assert not open_batch(S, bad, rng, _stats(), q)["ok"], \
                        (q, n, k, j, "corrupt batch claim accepted")

        # 5. all-zero table: S~(pt)=0 everywhere; honest accepts, residual 0
        for n in (2, 5):
            S = np.zeros(1 << n, dtype=_dt(q))
            pt = [int(rng.integers(0, q)) for _ in range(n)]
            res = open_eval(S, pt, 0, rng, _stats(), q)
            assert res["ok"] and res["residual"] == 0, (q, n, "zero table")
            assert not open_eval(S, pt, 1, rng, _stats(), q)["ok"], (q, n)

    # 6. agreement with the chain's mle_eval semantics on a REAL committed table
    #    (a compressed Lucy small-side layer), the actual use site.
    import compressed_layer as _cl
    for n in (8, 10, 12):
        x = (1 << n) - 1
        V = _cl.isqrt(x); nb = max(1, V.bit_length())
        _, sm, lg = _cl.compressed_lucy(x)
        for tab in (sm[len(sm) // 2], lg[len(lg) // 2]):
            S = _cl.pad(tab, nb)
            pts = [[int(rng.integers(0, Q)) for _ in range(nb)] for _ in range(3)]
            claims = [(pt, mle_eval(S, pt)) for pt in pts]
            res = open_batch(S, claims, rng, _stats(), Q)
            assert res["ok"] and close_eval(S, res["r"], res["residual"]), (n,)

    print("selftest OK")


def bench(q=Q):
    """The headline of the primitive: the open_eval VERIFIER is O(nb), flat in the
    table size 2^nb, while the direct mle_eval close it replaces is O(2^nb).  We
    report both, plus comm.  (The prover's eq-table + round-poly work IS O(2^nb) --
    that cost is the prover's, Õ(x) over the chain; the win is the VERIFIER term,
    which is what scales the whole-chain verifier from Õ(x) to Õ(sqrt x) once the
    residuals thread to the base.)  Falsifier: open_eval t_verifier growing with
    2^nb rather than ~nb."""
    rng = np.random.default_rng(7)
    print(f"q = {q};  open_eval VERIFIER (sum-check, O(nb)) vs direct mle_eval "
          f"close (O(2^nb))")
    print(f"{'nb':>3} {'2^nb':>9} {'open_tV_us':>11} {'mle_close_us':>13} "
          f"{'speedup':>8} {'open_comm':>10}")
    for nb in (8, 10, 12, 14, 16, 18):
        trials = 20 if nb <= 14 else 6
        S = (rng.integers(0, q, size=1 << nb)).astype(_dt(q))
        t_open = t_mle = 0.0
        comm = 0
        for t in range(trials):
            pt = [int(rng.integers(0, q)) for _ in range(nb)]
            true = mle_eval(S, pt, q)
            st = _stats()
            res = open_eval(S, pt, true, np.random.default_rng(100 + t), st, q)
            assert res["ok"]
            t_open += st["t_verifier"]
            comm = st["comm"]
            t0 = time.perf_counter()
            _ = mle_eval(S, pt, q)                          # the close it replaces
            t_mle += time.perf_counter() - t0
        to = t_open / trials * 1e6
        tm = t_mle / trials * 1e6
        print(f"{nb:>3} {1 << nb:>9} {to:>11.2f} {tm:>13.2f} "
              f"{tm / max(to, 1e-9):>7.1f}x {comm:>10}")


def bench_batch(q=Q, nb=14):
    """open_batch (one sum-check, residual threaded) vs the line-restriction fold
    it replaces: folding k same-table claims by sequential line_batch_pair costs
    (k-1) direct mle_eval closes = O(k 2^nb) on the VERIFIER; open_batch is O(k nb)
    + ONE residual.  Reports the verifier cost of each vs k at fixed nb."""
    import compressed_layer as _cl
    rng = np.random.default_rng(7)
    S = (rng.integers(0, q, size=1 << nb)).astype(_dt(q))
    print(f"q = {q}, nb = {nb} (2^nb = {1 << nb});  open_batch vs line-restriction "
          f"folding, VERIFIER cost vs #claims k")
    print(f"{'k':>3} {'open_tV_us':>11} {'open_comm':>10} || "
          f"{'lines_tV_us':>12} {'lines_comm':>11} {'speedup':>8}")
    for k in (2, 3, 5, 8):
        pts = [[int(rng.integers(0, q)) for _ in range(nb)] for _ in range(k)]
        claims = [(pt, mle_eval(S, pt, q)) for pt in pts]
        st = _stats()
        res = open_batch(S, claims, np.random.default_rng(1), st, q)
        assert res["ok"]
        # baseline: sequential line_batch_pair (the current batch_on_table path)
        sb = _stats()
        pt0, v0 = claims[0]
        for (pt1, v1) in claims[1:]:
            pt0, v0, okp = _cl.line_batch_pair(S, pt0, v0, pt1, v1, nb,
                                               np.random.default_rng(2), sb, q)
            assert okp
        print(f"{k:>3} {st['t_verifier']*1e6:>11.2f} {st['comm']:>10} || "
              f"{sb['t_verifier']*1e6:>12.2f} {sb['comm']:>11} "
              f"{sb['t_verifier']/max(st['t_verifier'],1e-9):>7.1f}x")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--bench", action="store_true")
    ap.add_argument("--bench-batch", action="store_true")
    ap.add_argument("--field", choices=list(FIELDS), default="q")
    args = ap.parse_args()
    if args.selftest:
        selftest()
    elif args.bench:
        bench(FIELDS[args.field])
    elif args.bench_batch:
        bench_batch(FIELDS[args.field])
    else:
        selftest()
