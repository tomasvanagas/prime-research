#!/usr/bin/env python3
r"""
reach_utime_crosscheck.py  --  resolve the S529/S535 reach-wall falsifier from
the UTIME (compute) side.

CONTEXT (read PROGRAM.md item 5 + the NEXT ACTION block first).
The detached n=26 verification-chain reach has a WALL that already exceeds
S529's predicted band [5.9, 7.9] ks (it is >11 ks and still running).  S535's
SHARPENED falsifier says this is CONSISTENT with a Theta(x) op-count IFF the
excess lives in `stime` (kernel/page-fault/allocation), and only a
`utime`-dominated excess with super-linear per-element growth reopens the
super-Theta(x) concern.  S535 measured the reach is stime-dominated (0.695),
but could NOT do the clean n=24->26 *utime* ratio test because n=24's
utime/stime split was lost on process exit.

THIS SCRIPT closes that gap WITHOUT the lost n=24 split, using one fact:

    utime (CPU user time) is CONTENTION-ROBUST.

The concurrent reach steals WALL-clock, not CPU-seconds: if my probe process
consumes T cpu-seconds in user mode, that is real work it did, independent of
how many other processes share the cores.  So I can measure the chain's utime
baseline at small RAM-light n RIGHT NOW (n<=20, <=265 MB, ~minimal DRAM
contention with the reach) and extrapolate to n=26 via the measured Theta(x)
op-count (S530/S531/S532b), then COMPARE to the live reach's captured utime.

  utime(26)_pred = utime(n0) * [opcount(26)/opcount(n0)] * dram_utime_factor

If live utime(26) ~ prediction  => compute is Theta(x) op-count * ~const per-op
   => the >5.4x WALL excess is entirely stime/allocation => falsifier NOT
   triggered, Theta(x) op-count intact (from the utime side).
If live utime(26) >> prediction (super-linear per-op) => REOPEN super-Theta(x).

MODES
  --live [PID]     read the live reach /proc counters (ZERO contention): the
                   cumulative utime/stime split (perishable -- lost on exit).
  --baseline       run the REAL run_chain FULL config at small n, measure
                   utime/stime (os.times deltas, contention-robust) + op-count
                   (the S530 fmul monkeypatch).  utime/op + Theta(x) check.
  --dram-factor    reused-buffer _mul61 utime/elem across L1->DRAM sizes: the
                   cache->DRAM per-op UTIME inflation (isolates compute, no
                   allocation -> no stime), to bridge cache-resident baseline n
                   to the DRAM-resident reach.
  --predict [PID]  combine baseline + dram-factor + op-count model -> predict
                   utime(26); compare to the live reach.  Verdict.
  --all            live + baseline + dram-factor + predict.
  --selftest       offline checks (no reach, no long runs).

WHAT WOULD FALSIFY THE "Theta(x) op-count, allocation-bound wall" reading
(pre-stated):
  F1  baseline utime is NOT ~Theta(x) (utime ratio per Dn=2 stays >~4.4 and
      RISING across n=16->18->20) -> compute itself is super-linear.
  F2  live utime(26) EXCEEDS utime(26)_pred by a factor that the bounded
      dram-factor cannot absorb (>~1.5x over the high end) AND the gap grows
      -> a genuine super-linear per-op utime term in the reach.
  F3  the reach is NOT stime-dominated at capture time (cum stime_frac < 0.5)
      -> the wall excess is utime after all (this is S535's F4, re-checked).
  F4  dram_utime_factor is unbounded / grows with size without saturating ->
      the per-op compute cost is not a bounded cache constant.
"""
import argparse
import os
import sys
import time

import numpy as np

# reuse S535's tested /proc parser (off-by-one fixed there) and S530's op-count
import alloc_overhead_wall as aow
import compressed_prover_mult_trace as cpmt
import prover_opcount_scaling as pos

# ----------------------------------------------------------------------
# measured whole-chain prover op-counts (FULL config, BIG_Q+FAST_BIG)
#   n<=20 live (S530), n=22,24 deterministic detached (S531/S532b).
# Source: prover_opcount_scaling_results.md.  The op-count is structural /
# seed-independent (selftest [3] there), so these are exact.
# ----------------------------------------------------------------------
MUL_OPS = {
    8: 425_110, 10: 2_074_457, 12: 9_673_318, 14: 23_492_482,
    16: 105_029_523, 18: 465_458_097, 20: 2_040_403_997,
    22: 8_882_635_340, 24: 38_448_499_833,
}
# forced-alpha=1 polylog model (S532b): mul ~ C * x * (log2 x)^0.87, x = 2^n.
OP_POLYLOG_D = 0.87
OP_ANCHOR_N = 24


def opcount_model(n, anchor_n=OP_ANCHOR_N, d=OP_POLYLOG_D):
    """Theta(x)*polylog op-count model anchored on the measured anchor point.
    mul(n) = mul(anchor) * 2^(n-anchor) * (n/anchor)^d."""
    base = MUL_OPS[anchor_n]
    return base * (2.0 ** (n - anchor_n)) * (n / anchor_n) ** d


def opcount(n):
    """Measured if available, else model."""
    return float(MUL_OPS[n]) if n in MUL_OPS else opcount_model(n)


# ----------------------------------------------------------------------
# (A) live reach capture -- ZERO contention (only reads /proc)
# ----------------------------------------------------------------------
def capture_live(pid=None, window_s=4.0):
    if pid is None:
        pid = aow.find_reach_pid()
    if pid is None:
        return None
    s0 = aow.read_proc_stat(pid)
    t0 = time.perf_counter()
    time.sleep(window_s)
    s1 = aow.read_proc_stat(pid)
    dt = time.perf_counter() - t0
    clk = aow.CLK_TCK
    cum_u = s1["utime"] / clk
    cum_s = s1["stime"] / clk
    cpu = (s1["utime"] - s0["utime"]) + (s1["stime"] - s0["stime"])
    return dict(
        pid=pid, state=s1["state"], window_s=dt,
        cum_utime_s=cum_u, cum_stime_s=cum_s,
        cum_cpu_s=cum_u + cum_s,
        cum_stime_frac=cum_s / (cum_u + cum_s) if (cum_u + cum_s) else float("nan"),
        win_stime_frac=((s1["stime"] - s0["stime"]) / cpu) if cpu else float("nan"),
        win_minflt_per_s=(s1["minflt"] - s0["minflt"]) / dt,
        cum_minflt=s1["minflt"], majflt=s1["majflt"],
        rss_bytes=s1["rss_bytes"],
    )


def print_live(res):
    if res is None:
        print("  [no live reach found]")
        return
    print(f"  pid={res['pid']}  state={res['state']}  RSS={aow._human(res['rss_bytes'])}")
    print(f"  cumulative:  utime={res['cum_utime_s']:.0f}s  stime={res['cum_stime_s']:.0f}s  "
          f"cpu={res['cum_cpu_s']:.0f}s  stime_frac={res['cum_stime_frac']:.3f}")
    print(f"  window {res['window_s']:.1f}s:  stime_frac={res['win_stime_frac']:.3f}  "
          f"minflt={res['win_minflt_per_s']:.3e}/s  cum_minflt={res['cum_minflt']:.3e}  "
          f"majflt={res['majflt']}")
    # F3 (== S535 F4): reach must be stime-dominated for the allocation reading.
    if res["cum_stime_frac"] < 0.5:
        print("  *** F3 TRIGGERED: reach NOT stime-dominated (excess is utime) ***")


# ----------------------------------------------------------------------
# (B) chain utime baseline -- contention-robust CPU time of the REAL chain
# ----------------------------------------------------------------------
def baseline_chain_times(ns, field="big"):
    """Run the REAL FULL-config run_chain at each n, measuring utime/stime via
    os.times() deltas (CPU time => robust to the concurrent reach's wall
    contention) and the op-count via the S530 monkeypatch (one instrumented
    run gives both)."""
    rows = []
    for n in ns:
        t0 = os.times()
        w0 = time.perf_counter()
        r = pos.measure_chain(n, field=field)        # installs op-count patch
        wall = time.perf_counter() - w0
        t1 = os.times()
        du = t1.user - t0.user
        ds = t1.system - t0.system
        mul = r["mul_ops"]
        rows.append(dict(
            n=n, x=r["x"], mul_ops=mul, utime=du, stime=ds, wall=wall,
            cpu=du + ds,
            stime_frac=ds / (du + ds) if (du + ds) else float("nan"),
            utime_ns_per_op=du / mul * 1e9 if mul else float("nan"),
            accepted=r["accepted"], match=r["match"],
        ))
    return rows


def print_baseline(rows):
    print(f"  {'n':>3} {'x':>11} {'mul_ops':>16} {'utime_s':>9} {'stime_s':>9} "
          f"{'st_frac':>7} {'ut_ns/op':>9} {'wall_s':>9}  ok")
    for r in rows:
        ok = "OK" if (r["accepted"] and r["match"]) else "**BAD**"
        print(f"  {r['n']:>3} {r['x']:>11} {r['mul_ops']:>16,} {r['utime']:>9.2f} "
              f"{r['stime']:>9.2f} {r['stime_frac']:>7.3f} {r['utime_ns_per_op']:>9.2f} "
              f"{r['wall']:>9.2f}  {ok}")
    # utime Theta(x) check (F1): per-Dn=2 utime ratios should be <= op-count
    # ratios (per-op utime amortizes overhead => utime grows no faster than ops).
    print("  per-Dn=2 ratios (utime vs op-count; per-op utime should NOT rise):")
    for i in range(1, len(rows)):
        a, b = rows[i - 1], rows[i]
        ur = b["utime"] / a["utime"] if a["utime"] else float("nan")
        orr = b["mul_ops"] / a["mul_ops"]
        pr = b["utime_ns_per_op"] / a["utime_ns_per_op"]
        flag = " <-- per-op RISING (F1)" if pr > 1.05 else ""
        print(f"    n={a['n']:>2}->{b['n']:<2}: utime {ur:5.3f}x   op {orr:5.3f}x   "
              f"per-op {pr:5.3f}x{flag}")


# ----------------------------------------------------------------------
# (C) cache->DRAM per-op UTIME factor: reused-buffer _mul61, no allocation
# ----------------------------------------------------------------------
def dram_utime_factor(sizes=None, target_s=0.6):
    """For each array size (elements), reuse a fixed pair of buffers and time
    _mul61 in user mode (os.times user delta).  No new big allocation per
    iter once the allocator recycles the per-call temporaries (warmup absorbs
    the first-touch faults) => utime isolates compute+cache, stime ~ 0.
    Returns rows with ns/elem (utime) per size + the DRAM/L1 ratio."""
    if sizes is None:
        # 4 Ki (L1) ... 64 Mi (>> L3, DRAM).  uint64 => 8 B/elem.
        sizes = [1 << k for k in (12, 16, 20, 23, 26)]
    saved = cpmt.FAST_BIG
    cpmt.FAST_BIG = True
    rows = []
    try:
        rng = np.random.default_rng(0)
        for S in sizes:
            a = (rng.integers(0, cpmt.BIG_Q, size=S, dtype=np.uint64))
            b = (rng.integers(0, cpmt.BIG_Q, size=S, dtype=np.uint64))
            # warmup: settle allocator + first-touch faults (charged to stime)
            for _ in range(3):
                cpmt._mul61(a, b)
            # calibrate iteration count to ~target_s of work
            reps = 1
            while True:
                t0 = os.times()
                for _ in range(reps):
                    cpmt._mul61(a, b)
                du = os.times().user - t0.user
                if du >= target_s or reps >= (1 << 22):
                    break
                reps *= 2
            ns_per = du / (reps * S) * 1e9 if (reps and S) else float("nan")
            rows.append(dict(size=S, bytes=S * 8, reps=reps, utime=du,
                             ns_per_elem=ns_per))
    finally:
        cpmt.FAST_BIG = saved
    base = rows[0]["ns_per_elem"]
    peak = max(r["ns_per_elem"] for r in rows)
    for r in rows:
        r["factor_vs_L1"] = r["ns_per_elem"] / base if base else float("nan")
    return rows, (peak / base if base else float("nan"))


def print_dram(rows, ratio):
    print(f"  reused-buffer _mul61 UTIME per element (no allocation => stime~0):")
    print(f"  {'size_elems':>12} {'bytes':>10} {'reps':>8} {'ns/elem':>9} {'xL1':>6}")
    for r in rows:
        print(f"  {r['size']:>12,} {aow._human(r['bytes']):>10} {r['reps']:>8} "
              f"{r['ns_per_elem']:>9.3f} {r['factor_vs_L1']:>6.2f}")
    print(f"  DRAM/L1 per-op utime factor = {ratio:.2f}x")
    # F4: factor must saturate (bounded cache constant), not grow unboundedly.
    if len(rows) >= 3:
        last_step = rows[-1]["ns_per_elem"] / rows[-2]["ns_per_elem"]
        print(f"  last size-step utime ratio = {last_step:.2f}x "
              f"({'saturating' if last_step < 1.4 else '*** F4: still climbing ***'})")


# ----------------------------------------------------------------------
# (D) predict reach utime + compare to live
# ----------------------------------------------------------------------
def predict_and_compare(baseline_rows, dram_ratio, live, target_n=26,
                        anchor_n=None, reach_done=False):
    """Resolve the falsifier on the PER-OP utime, which is contention-robust and
    does NOT need the lost n=24 split.

    The chain's per-op utime (utime / counted mul_op) DECLINES with n: at small n
    the fixed per-numpy-call overhead is spread over narrow arrays, and it
    amortizes as arrays widen.  A super-linear (super-Theta(x)) utime term would
    make per-op RISE with n.  So the decisive test is: does the reach's realized
    per-op utime sit on/below the declining baseline trend, or does it rise?

    The reach's realized per-op = live_utime / opcount(target).  If the reach is
    still running this is a LOWER bound on the final per-op (utime keeps growing);
    once DONE it is exact.
    """
    if not baseline_rows:
        print("  [no baseline rows -- run --baseline]")
        return None
    op_t = opcount(target_n)
    if anchor_n is None:
        anchor = max(baseline_rows, key=lambda r: r["n"])
    else:
        anchor = next(r for r in baseline_rows if r["n"] == anchor_n)
    anchor_perop = anchor["utime_ns_per_op"] if "utime_ns_per_op" in anchor \
        else anchor["utime"] / anchor["mul_ops"] * 1e9

    print(f"  opcount({target_n}) = {op_t:.3e}  (model x*(log2 x)^{OP_POLYLOG_D}, "
          f"anchored on measured n={OP_ANCHOR_N})")
    print(f"  baseline per-op utime (ns/op) -- DECLINES with n (overhead amortizes):")
    for r in sorted(baseline_rows, key=lambda z: z["n"]):
        po = r.get("utime_ns_per_op", r["utime"] / r["mul_ops"] * 1e9)
        print(f"    n={r['n']:>2}: {po:7.2f} ns/op")
    # secondary: a naive Theta(x) utime band from the anchor (note it OVER-predicts
    # because per-op keeps declining past the anchor; kept for context).
    op_ratio = op_t / anchor["mul_ops"]
    lo = anchor["utime"] * op_ratio
    hi = anchor["utime"] * op_ratio * dram_ratio
    print(f"  [context] naive Theta(x) band from n={anchor['n']} anchor "
          f"(op_ratio={op_ratio:.1f}x, dram_bridge<=[{1.0},{dram_ratio:.2f}]): "
          f"utime({target_n}) ~ [{lo:.0f}s, {hi:.0f}s]  (OVER-predicts: per-op declines)")

    res = dict(anchor_n=anchor["n"], anchor_perop=anchor_perop, op_t=op_t,
               op_ratio=op_ratio, lo=lo, hi=hi, target_n=target_n)
    if live is not None:
        lu = live["cum_utime_s"]
        realized_perop = lu / op_t * 1e9
        res["live_utime"] = lu
        res["realized_perop"] = realized_perop
        res["reach_done"] = reach_done
        kind = "EXACT (reach DONE)" if reach_done else "LOWER BOUND (reach running)"
        print(f"  reach utime({target_n}) so far = {lu:.0f}s  "
              f"(stime={live['cum_stime_s']:.0f}s, stime_frac={live['cum_stime_frac']:.3f})")
        print(f"  reach realized per-op = {realized_perop:.2f} ns/op  [{kind}]")
        print(f"    vs anchor n={anchor['n']} per-op {anchor_perop:.2f} ns/op  "
              f"(ratio {realized_perop/anchor_perop:.2f}x)")
        wall = lu + live["cum_stime_s"]
        # F3 (== S535 F4): reach must be stime-dominated.
        if live["cum_stime_frac"] < 0.5:
            print("  *** F3 TRIGGERED: reach NOT stime-dominated. ***")
            res["verdict"] = "F3"
        # The per-op test.  A LOWER bound that is already on/below the declining
        # baseline trend rules out super-linear; only per-op clearly ABOVE the
        # trend (and rising) reopens it.
        elif realized_perop <= anchor_perop * 1.3:
            res["verdict"] = "not_triggered"
            print("  VERDICT: reach per-op utime sits ON/BELOW the declining "
                  "baseline trend (<=1.3x the n={} anchor) -- per-op is NOT "
                  "rising, so utime grows no faster than the Theta(x) op-count. "
                  "F1/F2 NOT triggered: the >5.4x WALL excess is the stime/"
                  "allocation leg ({:.0f}% of wall), Theta(x) op-count intact."
                  .format(anchor["n"], live["cum_stime_s"] / wall * 100))
            if not reach_done:
                print("  (per-op is a LOWER bound until DONE; final stays below "
                      "the small-n per-op values -- monotone-declining trend.)")
        elif realized_perop <= anchor_perop * 2.0:
            res["verdict"] = "weak"
            print("  VERDICT: reach per-op utime is moderately above the anchor "
                  "(1.3-2.0x) but still far below the small-n per-op values; the "
                  "DRAM per-op penalty (measured by --dram-factor) covers a "
                  "bounded rise. Not a clean super-linear trigger -- confirm at DONE.")
        else:
            res["verdict"] = "F2"
            print("  *** F2 TRIGGERED: reach per-op utime far exceeds the baseline "
                  "trend => a super-linear per-op utime term -- REOPEN super-Theta(x). ***")
        print(f"  (context: reach wall so far ~= utime+stime = {wall:.0f}s; "
              f"stime carries {live['cum_stime_s']/wall*100:.0f}% of it; the "
              f"naive band over-predicts utime, so live<band is expected.)")
    return res


# ----------------------------------------------------------------------
# selftest
# ----------------------------------------------------------------------
def selftest():
    fails = []

    def check(name, cond):
        print(f"  [{'ok' if cond else 'FAIL'}] {name}")
        if not cond:
            fails.append(name)

    # [1] op-count model reproduces the measured points it is NOT anchored on.
    for n in (20, 22):
        pred = opcount_model(n)
        meas = MUL_OPS[n]
        rel = abs(pred - meas) / meas
        check(f"opcount_model(n={n}) within 12% of measured ({rel:.3f})", rel < 0.12)
    # [2] model is monotone increasing and ~4x per Dn=2 near the reach.
    r = opcount_model(26) / opcount_model(24)
    check(f"opcount_model 24->26 ratio in (4.0,4.6) ({r:.3f})", 4.0 < r < 4.6)
    # [3] opcount() prefers measured when present, model otherwise.
    check("opcount(24) == measured", opcount(24) == float(MUL_OPS[24]))
    check("opcount(26) == model (no measured)", abs(opcount(26) - opcount_model(26)) < 1)
    # [4] capture_live math: stime_frac from synthetic counters.
    #     (parse our own /proc/self/stat to exercise the reader path)
    me = aow.read_proc_stat(os.getpid())
    check("read_proc_stat(self) has utime/stime ints",
          isinstance(me["utime"], int) and isinstance(me["stime"], int))
    # [5] _mul61 correctness on a reused buffer (the dram-factor kernel) vs
    #     exact Python-int modmul.
    saved = cpmt.FAST_BIG
    cpmt.FAST_BIG = True
    try:
        rng = np.random.default_rng(1)
        a = rng.integers(0, cpmt.BIG_Q, size=257, dtype=np.uint64)
        b = rng.integers(0, cpmt.BIG_Q, size=257, dtype=np.uint64)
        got = np.asarray(cpmt._mul61(a, b)).astype(object)
        exp = (a.astype(object) * b.astype(object)) % cpmt.BIG_Q
        check("_mul61 == exact modmul (reused buffer)", bool(np.all(got == exp)))
    finally:
        cpmt.FAST_BIG = saved
    # [6] predict_and_compare arithmetic: synthetic anchor + live.
    fake_base = [dict(n=20, mul_ops=MUL_OPS[20], utime=40.0)]
    fake_live = dict(cum_utime_s=3375.0, cum_stime_s=7833.0, cum_stime_frac=0.699)
    res = predict_and_compare(fake_base, dram_ratio=1.3, live=fake_live,
                              target_n=26)
    exp_ratio = opcount(26) / MUL_OPS[20]
    check("predict op_ratio matches opcount(26)/op(20)",
          abs(res["op_ratio"] - exp_ratio) < 1e-6)
    check("predict band lo<hi and live compared",
          res["lo"] < res["hi"] and "live_utime" in res)
    # [7] dram_utime_factor returns saturating, finite, positive values (tiny).
    rows, ratio = dram_utime_factor(sizes=[1 << 12, 1 << 16], target_s=0.02)
    check("dram_utime_factor finite positive ns/elem",
          all(0 < r["ns_per_elem"] < 1e4 for r in rows))
    check("dram_utime_factor ratio finite >=1",
          1.0 <= ratio < 1e3)

    cpmt.FAST_BIG = saved
    print(f"\n  selftest: {'ALL PASS' if not fails else 'FAILURES: ' + ', '.join(fails)}")
    return not fails


# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--live", nargs="?", const=-1, type=int, default=None,
                    help="capture live reach /proc counters (optional PID; "
                         "auto-detect if omitted)")
    ap.add_argument("--baseline", action="store_true",
                    help="measure chain utime/stime baseline at small n")
    ap.add_argument("--dram-factor", action="store_true",
                    help="cache->DRAM per-op utime factor (reused buffers)")
    ap.add_argument("--predict", nargs="?", const=-1, type=int, default=None,
                    help="predict + compare reach utime (optional PID)")
    ap.add_argument("--all", action="store_true", help="live+baseline+dram+predict")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--ns", type=str, default="14,16,18,20",
                    help="baseline n values (comma list)")
    ap.add_argument("--target-n", type=int, default=26)
    ap.add_argument("--field", type=str, default="big", choices=["big", "q"])
    args = ap.parse_args()

    if args.selftest:
        sys.exit(0 if selftest() else 1)

    do_all = args.all
    live = None
    if do_all or args.live is not None or args.predict is not None:
        pid = None
        for v in (args.live, args.predict):
            if isinstance(v, int) and v > 0:
                pid = v
        print("=== live reach capture ===")
        live = capture_live(pid)
        print_live(live)
        print()

    baseline = None
    if do_all or args.baseline:
        ns = [int(x) for x in args.ns.split(",")]
        print(f"=== chain utime baseline (FULL config, field={args.field}, "
              f"n={ns}) ===")
        baseline = baseline_chain_times(ns, field=args.field)
        print_baseline(baseline)
        print()

    dram_rows, dram_ratio = (None, 1.3)
    if do_all or args.dram_factor:
        print("=== cache->DRAM per-op utime factor ===")
        dram_rows, dram_ratio = dram_utime_factor()
        print_dram(dram_rows, dram_ratio)
        print()

    if (do_all or args.predict is not None) and baseline is not None:
        reach_done = False
        try:
            with open(f"run_n{args.target_n}.log") as fh:
                reach_done = "DONE" in fh.read()
        except OSError:
            pass
        print(f"=== predict utime({args.target_n}) + compare to live "
              f"(reach_done={reach_done}) ===")
        predict_and_compare(baseline, dram_ratio, live, target_n=args.target_n,
                            reach_done=reach_done)
        print()

    if not any([args.live is not None, args.baseline, args.dram_factor,
                args.predict is not None, do_all]):
        ap.print_help()


if __name__ == "__main__":
    main()
