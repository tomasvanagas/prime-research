#!/usr/bin/env python3
"""
chain_profile.py — S501, attribute the compressed-chain field wall-clock.

The S500 correction left one concrete question open: the full compressed
delegated+structured pi(x) chain (compressed_layer.run_chain) is "op-count-
bound on many small per-layer cubes" -- so the lever for a fast large-x run is
reducing the COUNT of small numpy field-ops, not the per-multiply cost.  To
know WHICH small-cube kernel to batch we need the wall-clock attributed across
kernels.  This is that measurement: a READ-ONLY profiler over run_chain.

Two passes, because the two facts we need come from different instruments:

  (1) cProfile pass -> per-function (ncalls, tottime, cumtime).  cProfile does
      exclusive (self) timing at the C level with low overhead, so tottime is a
      clean non-double-counted partition of the wall.  We also pull the
      caller-breakdown of the two hot functions (fmul, sumcheck) so the
      function-level cost maps back onto the logical kernels (large_reduce /
      small_reduce / verify_trace_region / verify_affine_region /
      inner_verify_div / chain_layer_reduce_structured).

  (2) size pass -> per-primitive mean / max ARRAY SIZE (the piece cProfile
      cannot give).  A light monkeypatch that only accumulates a size counter
      (no perf_counter) wraps fmul / sumcheck / eq_table / mle_eval / _asum,
      so it barely perturbs the run and we DON'T use it for timing.  sumcheck
      calls are bucketed by (nvars, degree, n_terms) -- a fingerprint that
      separates the layer cubes (nvars=nb) from the inner wiring cubes
      (nvars=l~log p) from the trace zero-test (many terms).

CORRECTION of the NEXT-ACTION premise.  PROGRAM.md asked this profiler to
assert "the profile sums match run_chain's t_prover+t_verifier".  That is
WRONG and this file documents why: run_chain's sumcheck() calls are NOT inside
any t_prover/t_verifier accounting region (sumcheck never touches stats), and
the sum-check loop is ~92% of the wall.  So t_prover+t_verifier accounts for
only ~7% of run_chain's wall -- the reported per-layer stats massively
undercount the true cost.  The honest check this file makes instead:
faithfulness (instrumented run == reference: claimed pi, comm, accept), the two
instruments agreeing on call counts, and the kernel table covering ~all the
wall.  The 7% gap is reported as a headline finding.

What would falsify the results here: the instrumented run disagreeing with the
reference run_chain (claimed pi / comm / accept); the cProfile ncalls and the
size-pass call counts disagreeing for the same primitive; the ranked kernel
table failing to cover the bulk of the wall (a hidden cost); or fmul NOT being
the dominant tottime with a small mean array size (which would refute the
"op-count-bound on small arrays" reading).

Usage:
  python3 chain_profile.py --selftest
  python3 chain_profile.py --n 16 --field big          # the headline profile
  python3 chain_profile.py --n 14 --field q            # uint64 demo field
  python3 chain_profile.py --n 16 --field big --rows 18
"""

import argparse
import cProfile
import io
import pstats
import time
from collections import defaultdict
from contextlib import contextmanager
from math import isqrt

import numpy as np

import compressed_prover_mult_trace as cpmt
import compressed_layer as cl
import lucy_dp_delegated_wiring as ldw
from compressed_prover_mult_trace import DEFAULT_Q as Q, BIG_Q, FIELDS


# ----------------------------------------------------------------------
# the run under test: the production-config compressed chain
# ----------------------------------------------------------------------

def run(n, q, delegate=True, structured=True, seed=5):
    """One full compressed pi(x) chain, x = 2^n - 1.  Returns run_chain's dict
    plus measured wall time."""
    x = (1 << n) - 1
    t0 = time.perf_counter()
    res = cl.run_chain(x, np.random.default_rng(seed), delegate=delegate,
                       structured=structured, q=q)
    res = dict(res)
    res["wall"] = time.perf_counter() - t0
    res["x"] = x
    res["V"] = isqrt(x)
    res["nb"] = max(1, isqrt(x).bit_length())
    return res


# ----------------------------------------------------------------------
# pass 1 — cProfile: per-function tottime / cumtime / ncalls
# ----------------------------------------------------------------------

# functions we name in the ranked table (basename-agnostic: matched by funcname,
# summed over lines if a name appears more than once).
_NAMED = [
    "fmul", "sumcheck", "_sumcheck_rec", "lagrange_eval", "_asum", "eq_table",
    "mle_eval", "eq_point", "aff_rel_eval", "ge_const_eval", "chain_layer_reduce_structured",
    "inner_chain_vectors", "inner_r_table", "_modmul_arr", "verify_constraints",
    "verify_trace_region", "verify_affine_region", "verify_waff_value",
    "large_reduce", "small_reduce", "inner_verify_div", "line_batch_pair",
    "batch_on_table", "build_witness", "build_terms", "constraint_eval",
    "run_chain", "compressed_lucy", "pad",
]


def cprofile_run(n, q, **kw):
    pr = cProfile.Profile()
    pr.enable()
    res = cl.run_chain((1 << n) - 1, np.random.default_rng(kw.get("seed", 5)),
                       delegate=kw.get("delegate", True),
                       structured=kw.get("structured", True), q=q)
    pr.disable()
    ps = pstats.Stats(pr)
    return ps, res


def func_rows(ps):
    """Collapse pstats by funcname (over its filename:lineno keys).  Returns
    name -> dict(ncalls, tottime, cumtime)."""
    agg = defaultdict(lambda: dict(ncalls=0, tottime=0.0, cumtime=0.0))
    for (fn, lineno, name), (cc, nc, tt, ct, callers) in ps.stats.items():
        a = agg[name]
        a["ncalls"] += nc
        a["tottime"] += tt
        # cumtime double counts on recursion; take the max line's cumtime as a
        # representative (recursive funcs flagged separately if needed).
        a["cumtime"] = max(a["cumtime"], ct)
    return agg


def callers_of(ps, target):
    """Return {caller_funcname: ncalls} for direct callers of `target`."""
    out = defaultdict(int)
    for (fn, lineno, name), (cc, nc, tt, ct, callers) in ps.stats.items():
        if name != target:
            continue
        for (cfn, cl_, cname), callinfo in callers.items():
            ncc = callinfo[0] if isinstance(callinfo, tuple) else callinfo
            out[cname] += ncc
    return out


# ----------------------------------------------------------------------
# pass 2 — size instrumentation (mean / max array size per primitive)
# ----------------------------------------------------------------------

def _arr_size(a):
    return int(getattr(a, "size", 1))


class SizeRec:
    """Accumulates per-primitive call count + array-size stats with negligible
    overhead (no timing).  sumcheck is also bucketed by (nvars,degree,nterms)."""

    def __init__(self):
        self.n = defaultdict(int)
        self.size_sum = defaultdict(int)
        self.size_max = defaultdict(int)
        self.sc_buckets = defaultdict(lambda: [0, 0])   # key -> [count, size_sum]

    def add(self, name, size):
        self.n[name] += 1
        if size is not None:
            self.size_sum[name] += size
            if size > self.size_max[name]:
                self.size_max[name] = size

    def add_sumcheck(self, tables, terms, degree):
        first = next(iter(tables.values()))
        size = _arr_size(first)
        nvars = int(round(np.log2(len(first))))
        self.add("sumcheck", size)
        key = (nvars, degree, len(terms))
        b = self.sc_buckets[key]
        b[0] += 1
        b[1] += size


@contextmanager
def patch_sizes(rec):
    """Monkeypatch the size-recording wrappers across every module that holds a
    binding to the primitive (compressed_layer / lucy_dp_delegated_wiring import
    the names directly, so each binding must be patched)."""
    modules = [cpmt, cl, ldw]
    o_fmul = cpmt.fmul
    o_asum = cpmt._asum
    o_eqt = cpmt.eq_table
    o_mle = cpmt.mle_eval
    o_sc = cpmt.sumcheck

    def w_fmul(a, b, q):
        rec.add("fmul", max(_arr_size(a), _arr_size(b)))
        return o_fmul(a, b, q)

    def w_asum(arr, q):
        rec.add("_asum", _arr_size(arr))
        return o_asum(arr, q)

    def w_eqt(z, q=Q):
        rec.add("eq_table", 1 << len(z))
        return o_eqt(z, q)

    def w_mle(table, pt, q=Q):
        rec.add("mle_eval", _arr_size(table))
        return o_mle(table, pt, q)

    def w_sc(claim, tables, terms, degree, rng, q=Q):
        rec.add_sumcheck(tables, terms, degree)
        return o_sc(claim, tables, terms, degree, rng, q)

    saved = []
    for m in modules:
        for nm, w in [("fmul", w_fmul), ("_asum", w_asum), ("eq_table", w_eqt),
                      ("mle_eval", w_mle), ("sumcheck", w_sc)]:
            if hasattr(m, nm):
                saved.append((m, nm, getattr(m, nm)))
                setattr(m, nm, w)
    try:
        yield rec
    finally:
        for m, nm, orig in saved:
            setattr(m, nm, orig)


def size_run(n, q, **kw):
    rec = SizeRec()
    with patch_sizes(rec):
        res = cl.run_chain((1 << n) - 1,
                           np.random.default_rng(kw.get("seed", 5)),
                           delegate=kw.get("delegate", True),
                           structured=kw.get("structured", True), q=q)
    return rec, res


# ----------------------------------------------------------------------
# report
# ----------------------------------------------------------------------

def report(n, q, rows=16, delegate=True, structured=True, seed=5):
    fld = ("big=2^61-1 [object]" if q == BIG_Q and not cpmt.fast_big(q)
           else "big=2^61-1 [uint64 Mersenne]" if q == BIG_Q
           else f"q={q}")
    x = (1 << n) - 1
    V = isqrt(x); nb = max(1, V.bit_length())
    K = len(cl.compressed_lucy(x)[0])
    print(f"=== compressed chain profile: x = 2^{n}-1, V = {V}, nb = {nb}, "
          f"K = {K} layers ===")
    print(f"field {fld};  mode = "
          f"{'delegate' if delegate else 'automaton'}"
          f"{'+structured' if (delegate and structured) else ''}")

    ref = run(n, q, delegate, structured, seed)
    print(f"reference run_chain: claimed pi = {ref['claimed']} "
          f"(sieve {cl.sieve_pi(x)}), accepted = {ref['accepted']}, "
          f"wall = {ref['wall']*1000:.0f} ms")
    acct = (ref["t_prover"] + ref["t_verifier"])
    print(f"  reported t_prover = {ref['t_prover']*1000:.1f} ms, "
          f"t_verifier = {ref['t_verifier']*1000:.1f} ms  "
          f"-> SUM = {acct*1000:.1f} ms = {100*acct/ref['wall']:.1f}% of wall")
    print(f"  *** the other {100*(1-acct/ref['wall']):.0f}% is the sum-check "
          f"loop, which never touches stats (NEXT-ACTION premise corrected) ***")

    ps, _ = cprofile_run(n, q, delegate=delegate, structured=structured,
                         seed=seed)
    total = ps.total_tt
    agg = func_rows(ps)
    rec, _ = size_run(n, q, delegate=delegate, structured=structured, seed=seed)

    print(f"\n--- ranked kernels (cProfile, total profiled = {total*1000:.0f} "
          f"ms) ---")
    print(f"{'kernel':>30} {'tottime_ms':>10} {'%tot':>6} {'cumtime_ms':>10} "
          f"{'%cum':>6} {'ncalls':>9} {'mean_sz':>8} {'max_sz':>7}")
    named = [(nm, agg[nm]) for nm in _NAMED if nm in agg]
    named.sort(key=lambda kv: -kv[1]["tottime"])
    cov = 0.0
    for nm, a in named:
        cov += a["tottime"]
        ms = (rec.size_sum[nm] / rec.n[nm]) if rec.n.get(nm) else None
        mx = rec.size_max.get(nm)
        ms_s = f"{ms:8.1f}" if ms is not None else f"{'-':>8}"
        mx_s = f"{mx:7d}" if mx else f"{'-':>7}"
        print(f"{nm:>30} {a['tottime']*1000:>10.1f} "
              f"{100*a['tottime']/total:>5.1f}% {a['cumtime']*1000:>10.1f} "
              f"{100*a['cumtime']/total:>5.1f}% {a['ncalls']:>9} {ms_s} {mx_s}")
    print(f"{'(named kernels tottime)':>30} {cov*1000:>10.1f} "
          f"{100*cov/total:>5.1f}% of profiled wall")

    print(f"\n--- where fmul is called (direct callers, ncalls) ---")
    for cn, nc in sorted(callers_of(ps, "fmul").items(), key=lambda kv: -kv[1]):
        print(f"  {cn:>34}: {nc:>9}")
    print(f"--- where sumcheck is called (direct callers, ncalls) ---")
    for cn, nc in sorted(callers_of(ps, "sumcheck").items(),
                         key=lambda kv: -kv[1]):
        print(f"  {cn:>34}: {nc:>9}")

    print(f"\n--- sumcheck calls bucketed by (nvars, degree, n_terms) ---")
    print(f"{'nvars':>6} {'deg':>4} {'nterms':>7} {'count':>7} {'mean_sz':>8}  "
          f"role")
    roles = {(nb, 3): "layer phase-A / trace zero-test (nb-cube)",
             (nb, 2): "affine/B1 region or small phase-B (nb-cube)",
             (nb, nb + 2): "trace B2 (nb-cube, deg nb+2)"}
    for key in sorted(rec.sc_buckets, key=lambda k: -rec.sc_buckets[k][1]):
        nvars, deg, nterms = key
        cnt, ssum = rec.sc_buckets[key]
        role = roles.get((nvars, deg), "")
        if nvars < nb and not role:
            role = "inner wiring cube (nvars~log2 p)"
        print(f"{nvars:>6} {deg:>4} {nterms:>7} {cnt:>7} {ssum/cnt:>8.1f}  "
              f"{role}")

    _verdict(agg, rec, total, nb, K)
    return ref, ps, rec


def _verdict(agg, rec, total, nb, K):
    fmul = agg.get("fmul", {})
    fmul_calls = rec.n.get("fmul", 0)
    fmul_mean = (rec.size_sum["fmul"] / fmul_calls) if fmul_calls else 0
    print(f"\n--- VERDICT ---")
    print(f"fmul: {fmul.get('tottime',0)*1000:.0f} ms "
          f"({100*fmul.get('tottime',0)/total:.0f}% of profiled wall) over "
          f"{fmul_calls} calls, mean array size {fmul_mean:.1f} elems.")
    print("This is the op-count-bound signature S500 predicted: the wall is a "
          "few-hundred-thousand-to-million tiny element-wise field multiplies, "
          "each paying full Python+numpy per-call dispatch on arrays of ~%.0f "
          "elements." % fmul_mean)
    print("HIGHEST-LEVERAGE BATCHING TARGET: the K per-layer sum-checks are "
          "independent in their challenge draws across the trace zero-tests / "
          "region checks (each draws its own tau/alpha) and share an identical "
          "shape per logical step.  Stacking the K layers' nb-cube tables into "
          "one (K*2^nb)-wide array and running ONE sum-check per logical step "
          "cuts the fmul COUNT ~K-fold while widening each fmul ~K-fold -- "
          "converting op-count-bound to width-bound, where numpy (and fmul's "
          "vectorised Mersenne path) finally wins.  verify_trace_region / "
          "verify_constraints is the fattest single target (most terms => most "
          "fmul per call).")


# ----------------------------------------------------------------------
# selftest
# ----------------------------------------------------------------------

def selftest():
    rng_seed = 5
    # 1. FAITHFULNESS: the size-instrumented run reproduces the reference
    #    run_chain exactly (claimed pi, comm, accept) -- the monkeypatch only
    #    counts, never alters arithmetic -- over both fields and both modes.
    for n in (8, 10):
        x = (1 << n) - 1
        truth = cl.sieve_pi(x)
        for q in (Q, BIG_Q):
            for (deleg, struct) in [(True, True), (False, False)]:
                ref = run(n, q, deleg, struct, rng_seed)
                rec, ires = size_run(n, q, delegate=deleg, structured=struct,
                                     seed=rng_seed)
                assert ref["claimed"] == truth, (n, q, deleg, ref["claimed"])
                assert ires["claimed"] == ref["claimed"], (n, q, deleg)
                assert ires["comm"] == ref["comm"], (n, q, deleg)
                assert ires["accepted"] == ref["accepted"], (n, q, deleg)

    # 2. INSTRUMENT AGREEMENT: cProfile ncalls == size-pass counts for the hot
    #    primitives (both instruments observe the same calls).  Different runs
    #    of the same (n,q,seed,mode) are deterministic, so counts must match.
    for n in (8, 10):
        for q in (Q, BIG_Q):
            ps, _ = cprofile_run(n, q)
            rec, _ = size_run(n, q)
            agg = func_rows(ps)
            for nm in ("fmul", "sumcheck", "eq_table", "mle_eval", "_asum"):
                assert agg[nm]["ncalls"] == rec.n[nm], \
                    (n, q, nm, agg[nm]["ncalls"], rec.n[nm])

    # 3. THE GAP FINDING: t_prover+t_verifier is a small fraction of wall
    #    (sum-check loop is untimed).  Robust, huge-margin bound.
    ref = run(12, BIG_Q)
    frac = (ref["t_prover"] + ref["t_verifier"]) / ref["wall"]
    assert frac < 0.30, ("accounted fraction unexpectedly high", frac)

    # 4. ATTRIBUTION: fmul is the dominant tottime and operates on SMALL arrays;
    #    the named kernels cover the bulk of the profiled wall.
    ps, _ = cprofile_run(12, BIG_Q)
    rec, _ = size_run(12, BIG_Q)
    agg = func_rows(ps)
    total = ps.total_tt
    named_tot = sum(agg[nm]["tottime"] for nm in _NAMED if nm in agg)
    assert named_tot / total > 0.85, ("named kernels miss too much", named_tot / total)
    # fmul is the single largest tottime contributor
    top = max((agg[nm]["tottime"], nm) for nm in agg)[1]
    assert top == "fmul", ("expected fmul dominant", top)
    fmul_mean = rec.size_sum["fmul"] / rec.n["fmul"]
    assert fmul_mean < 200, ("fmul arrays not small", fmul_mean)   # op-count bound

    # 5. sumcheck cumtime dominates the wall (the loop, not the verifier algebra)
    assert agg["sumcheck"]["cumtime"] / total > 0.5, "sumcheck not dominant"

    print("selftest OK")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=16)
    ap.add_argument("--field", choices=list(FIELDS), default="big")
    ap.add_argument("--rows", type=int, default=16,
                    help="(unused placeholder kept for symmetry)")
    ap.add_argument("--seed", type=int, default=5)
    ap.add_argument("--automaton", action="store_true",
                    help="profile the O(nb p) automaton-wiring chain instead "
                         "of delegate+structured")
    ap.add_argument("--dense", action="store_true",
                    help="delegate but dense (O(nb p^2)) inner prover")
    ap.add_argument("--fast-big", action="store_true",
                    help="route BIG_Q through the uint64 Mersenne path")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return
    if args.field == "big" and args.fast_big:
        cpmt.FAST_BIG = True
    q = FIELDS[args.field]
    delegate = not args.automaton
    structured = delegate and not args.dense
    report(args.n, q, rows=args.rows, delegate=delegate, structured=structured,
           seed=args.seed)


if __name__ == "__main__":
    main()
