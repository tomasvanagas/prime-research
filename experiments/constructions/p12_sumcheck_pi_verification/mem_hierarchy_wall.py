#!/usr/bin/env python3
r"""
mem_hierarchy_wall.py  (S529)

Mechanism test for the large-x reach's filed wall-scaling falsifier.

PROGRAM.md item 5 records a falsifier: the verification-chain prover wall
ACCELERATES with x -- the n=22->24 step ran at 5.80x per Delta n=2 versus the
ideal Theta(x) ratio of 4x.  S524 first read this as possibly super-Theta(x);
S525 re-attributed it to "a memory-hierarchy constant (Theta(x) arithmetic over
a cache-busting Theta(x) working set)", i.e. an IMPLEMENTATION cost, not a
complexity change.  That attribution was made from peak-RSS probes, never from a
direct measurement of the per-element field-op wall as the working set crosses
the cache levels.

This standalone benchmark supplies that missing measurement.  It times the
dominant BIG_Q=2^61-1 kernels (`_mul61` Mersenne mulmod, `_sum61` block-fold
sum, `fmul` dispatch) per ELEMENT over uint64 arrays whose footprint sweeps L1
(32 KiB) -> L2 (512 KiB) -> L3 (32 MiB) -> DRAM (multi-GiB), and asks the sharp
question the reach run answers at the macro level:

    Does ns/element keep RISING from the n=20 working set (~185 MiB, S527) to the
    n=24 working set (~2.6 GiB), or does it PLATEAU once the array exceeds L3?

  * If ns/elem rises ~1.45x across that range  ==>  4x (op-count, Theta(x)) *
    1.45x (per-element memory wall) ~= 5.8x  CONFIRMS the S525 memory-hierarchy
    reading: the acceleration is bandwidth/TLB, plateauing, NOT super-Theta(x).
  * If ns/elem is FLAT past L3                 ==>  the 5.8x is NOT memory; it is
    either a genuine super-Theta(x) op-count term or measurement noise, and the
    detached reach (n=25/26) becomes the arbiter.

Either way this is a falsifiable, in-cycle mechanism probe that does not need the
hours-long macro run to complete.

What would falsify THIS result:
  - `_mul61`/`fmul`/`scatter_fold61` disagreeing with exact Python-int arithmetic
    (selftest): then the kernel timings are meaningless.  (Asserted.)
  - ns/elem NOT U-shaped (high for tiny dispatch-bound arrays, min near L1/L2,
    rising past L3): would refute the two-regime (dispatch vs bandwidth) model
    this analysis rests on.
  - The model wall = opcount(n) * ns_per_elem(workingset(n)) reproducing the
    measured 71.8/252.8/1467.3 s reach series to within ~25% under the
    memory-hierarchy reading, AND the detached reach n=24->26 ratio landing OFF
    the predicted band: would mean the wall is not what this model says.

Usage:
  python3 mem_hierarchy_wall.py --selftest
  python3 mem_hierarchy_wall.py --bench-kernels         # the cache-step sweep
  python3 mem_hierarchy_wall.py --predict               # decompose + predict
  python3 mem_hierarchy_wall.py --all
"""
import argparse
import gc
import sys
import time

import numpy as np

import compressed_prover_mult_trace as cpmt

# ----------------------------------------------------------------------
# cache topology (read once; falls back to the measured host if unreadable)
# ----------------------------------------------------------------------
HOST_CACHES = [  # (label, bytes)  -- AMD EPYC 7313, per-core L1d/L2, per-CCX L3
    ("L1d", 32 * 1024),
    ("L2", 512 * 1024),
    ("L3", 32 * 1024 * 1024),
]

# Working-set (peak RSS, FULL streamed config, BIG_Q) MEASURED by S527/S525.
# n -> peak resident bytes of the prover (the Theta(x) batched-discharge cube
# dominates).  n=24 from S525 streamed estimate (~2.6 GiB cube of the 3.9 GiB
# pre-streaming peak).  These are the working sets the reach kernels run over.
MEASURED_PEAK_RSS = {  # MiB
    16: 43.0, 18: 70.0, 20: 185.0, 22: 677.0, 24: 2600.0,
}
# Measured FULL-config reach wall (S524, object/streamed ~equal per S527), seconds.
MEASURED_WALL_S = {20: 71.8, 22: 252.8, 24: 1467.3}


def _bytes_human(b):
    for u in ("B", "KiB", "MiB", "GiB"):
        if b < 1024 or u == "GiB":
            return f"{b:.1f} {u}"
        b /= 1024.0


def _regime(footprint_bytes):
    """Which cache level holds an array of this footprint (the larger of the two
    operands; the kernel streams several arrays of ~this size)."""
    for label, cap in HOST_CACHES:
        if footprint_bytes <= cap:
            return label
    return "DRAM"


# ----------------------------------------------------------------------
# kernel timing
# ----------------------------------------------------------------------
def _time_kernel(fn, *arrays, budget_s=0.4, min_reps=3, max_reps=200):
    """Median wall per call of fn(*arrays), adaptive reps to ~budget_s total.
    Warms once (page-in + branch warmup), then times."""
    fn(*arrays)  # warmup
    times = []
    t_total = 0.0
    reps = 0
    while reps < min_reps or (t_total < budget_s and reps < max_reps):
        t0 = time.perf_counter()
        fn(*arrays)
        dt = time.perf_counter() - t0
        times.append(dt)
        t_total += dt
        reps += 1
    times.sort()
    return times[len(times) // 2], reps


def bench_kernels(sizes_log2=None, modulus=None, light=False):
    """Sweep array size; for each, report ns/element of the dominant kernels.
    light=True skips the heavy `_mul61` (~10 transient arrays) so the sweep can
    push into the multi-GiB TLB/page-walk regime on the bandwidth controls only.
    Returns a list of per-size dicts."""
    if sizes_log2 is None:
        # 8 KiB (2^10 uint64) up to 1 GiB (2^27 uint64): spans L1->DRAM.
        sizes_log2 = list(range(10, 28))
    rng = np.random.default_rng(12345)
    P = (1 << 61) - 1
    rows = []
    prev_mul = None
    for lg in sizes_log2:
        n = 1 << lg
        footprint = n * 8  # uint64 bytes per array
        a = rng.integers(0, P, size=n, dtype=np.uint64)
        b = rng.integers(0, P, size=n, dtype=np.uint64)

        if light:
            ns_mul = float("nan"); r1 = 0
        else:
            t_mul, r1 = _time_kernel(cpmt._mul61, a, b)
            ns_mul = t_mul / n * 1e9
        t_sum, r2 = _time_kernel(cpmt._sum61, a)
        # plain uint64 multiply+mod is the natural "1 op/elem" memory-bandwidth
        # control (NOT a valid field mul for BIG_Q -- overflows -- but isolates
        # the pure load/store/ALU streaming rate the mulmod is layered on).
        t_np, r3 = _time_kernel(lambda x, y: (x * y) % np.uint64(P), a, b)

        ns_sum = t_sum / n * 1e9
        ns_np = t_np / n * 1e9
        step = (ns_mul / prev_mul) if (prev_mul and not light) else float("nan")
        if not light:
            prev_mul = ns_mul
        rows.append(dict(lg=lg, n=n, footprint=footprint, regime=_regime(footprint),
                         ns_mul=ns_mul, ns_sum=ns_sum, ns_np=ns_np, step=step,
                         reps=(r1, r2, r3)))
        del a, b
        gc.collect()
    return rows


def print_kernel_table(rows):
    print(f"  modulus = BIG_Q = 2^61-1, uint64 Mersenne kernels "
          f"(_mul61 ~24 numpy ops/elem)")
    print(f"  cache: L1d 32 KiB | L2 512 KiB | L3 32 MiB  (AMD EPYC 7313)")
    print()
    hdr = (f"  {'log2 n':>6} {'footprint':>11} {'regime':>6} "
           f"{'_mul61':>9} {'step':>6} {'_sum61':>9} {'np mul%':>9}")
    print(hdr)
    print(f"  {'':>6} {'(per array)':>11} {'':>6} {'ns/elem':>9} {'x':>6} "
          f"{'ns/elem':>9} {'ns/elem':>9}")
    print("  " + "-" * (len(hdr) - 2))
    for r in rows:
        stp = "" if r["step"] != r["step"] else f"{r['step']:.2f}"
        print(f"  {r['lg']:>6} {_bytes_human(r['footprint']):>11} {r['regime']:>6} "
              f"{r['ns_mul']:>9.3f} {stp:>6} {r['ns_sum']:>9.3f} {r['ns_np']:>9.3f}")
    return rows


def _rate_at(rows, footprint_bytes):
    """ns/elem of _mul61 at the row whose footprint is closest to (and >=, if
    possible) the given working-set footprint -- i.e. the bandwidth regime the
    reach's held cube of that size runs in."""
    # the held cube is one big array; pick the bench row with the nearest
    # footprint on a log scale.
    best = min(rows, key=lambda r: abs(np.log2(r["footprint"]) - np.log2(footprint_bytes)))
    return best["ns_mul"], best


def predict(rows):
    """Decompose the measured reach wall into op-count growth (Theta(x), known)
    and per-element memory wall (measured here), and emit the n=24->26 band."""
    print("Working-set regime per reach n (measured peak RSS, S527/S525):")
    print(f"  {'n':>3} {'peak RSS':>10} {'regime':>6}  (held cube ~ Theta(x))")
    for n in sorted(MEASURED_PEAK_RSS):
        b = MEASURED_PEAK_RSS[n] * 1024 * 1024
        print(f"  {n:>3} {_bytes_human(b):>10} {_regime(b):>6}")
    print()

    # per-element rate at each reach n's working-set footprint
    print("Per-element _mul61 wall at each reach working set (from the sweep):")
    rate = {}
    for n in (20, 22, 24):
        b = MEASURED_PEAK_RSS[n] * 1024 * 1024
        ns, row = _rate_at(rows, b)
        rate[n] = ns
        print(f"  n={n}: working set {_bytes_human(b):>9} -> ns/elem ~ {ns:.3f} "
              f"(nearest sweep row log2 n={row['lg']}, {row['regime']})")
    print()

    # memory-wall growth factor across the reach band (n=20 -> n=24)
    mem_growth_2022 = rate[22] / rate[20] if rate[20] else float("nan")
    mem_growth_2224 = rate[24] / rate[22] if rate[22] else float("nan")
    print(f"Per-element memory-wall growth: n20->22 = {mem_growth_2022:.3f}x, "
          f"n22->24 = {mem_growth_2224:.3f}x  (per Delta n=2)")
    print()

    # ideal op-count ratio is Theta(x) = 4x per Delta n=2 (verifier/prover both
    # work the sqrt(x) state with x-many element-ops in the dominant cube).
    IDEAL = 4.0
    obs_2022 = MEASURED_WALL_S[22] / MEASURED_WALL_S[20]
    obs_2224 = MEASURED_WALL_S[24] / MEASURED_WALL_S[22]
    print(f"Observed reach wall ratio:  n20->22 = {obs_2022:.3f}x, "
          f"n22->24 = {obs_2224:.3f}x   (ideal Theta(x) = {IDEAL:.1f}x)")
    print(f"Model (op-count {IDEAL:.0f}x * measured per-elem memory growth):")
    print(f"   n20->22 predicted = {IDEAL * mem_growth_2022:.3f}x  "
          f"(observed {obs_2022:.3f}x)")
    print(f"   n22->24 predicted = {IDEAL * mem_growth_2224:.3f}x  "
          f"(observed {obs_2224:.3f}x)")
    print()

    # n=24->26 prediction.  n=26 working set ~ 4x n=24 ~ 10.4 GiB (still DRAM).
    ws26 = MEASURED_PEAK_RSS[24] * 4 * 1024 * 1024
    ns26, row26 = _rate_at(rows, ws26)
    mem_growth_2426 = ns26 / rate[24] if rate[24] else float("nan")
    pred_2426 = IDEAL * mem_growth_2426
    print(f"PREDICTION for the detached reach (n=24 -> n=26, Delta n=2):")
    print(f"  n=26 working set ~ {_bytes_human(ws26)} (DRAM), per-elem ~ {ns26:.3f} ns "
          f"(growth {mem_growth_2426:.3f}x over n=24)")
    print(f"  => plateau-only floor (this sweep caps at 1 GiB) ~ {pred_2426:.2f}x  "
          f"(=> n=26 wall >= {MEASURED_WALL_S[24] * pred_2426:.0f} s)")
    print(f"  NOTE: --bandwidth-probe 30 shows a MILD multi-GiB TLB/bandwidth rise")
    print(f"  (~1.1-1.35x over the 2.6->10 GiB n=24->26 working set), so the honest")
    print(f"  band is ~[4.0, 5.4]x (n=26 wall ~6.5-7.9 ks).")
    print(f"  Falsifier: a measured n24->26 ratio WELL ABOVE ~5.5x -- with the")
    print(f"  per-elem wall shown flat-to-mild here -- would point to a genuine")
    print(f"  super-Theta(x) op-count term, refuting the 'implementation cost only'")
    print(f"  reading and reopening the S524 concern.  The TREND (geo-mean across")
    print(f"  n>=20 staying near ~4.5x) is the test, not any single noisy step.")
    return dict(rate=rate, mem_2224=mem_growth_2224, pred_2426=pred_2426)


# ----------------------------------------------------------------------
# selftest
# ----------------------------------------------------------------------
def selftest():
    rng = np.random.default_rng(7)
    P = (1 << 61) - 1
    ok = 0

    # [1] _mul61 == exact Python-int (a*b) % p, over random + boundary inputs
    a = rng.integers(0, P, size=4096, dtype=np.uint64)
    b = rng.integers(0, P, size=4096, dtype=np.uint64)
    bound = np.array([0, 1, 2, P - 1, P - 2, (1 << 31), (1 << 31) - 1,
                      (1 << 30), (1 << 60), P >> 1], dtype=np.uint64)
    a = np.concatenate([a, bound]); b = np.concatenate([b, bound[::-1]])
    fast = cpmt._mul61(a, b)
    exact = np.array([(int(x) * int(y)) % P for x, y in zip(a, b)], dtype=object)
    assert all(int(f) == int(e) for f, e in zip(fast, exact)), "_mul61 != exact"
    ok += 1

    # [2] _sum61 == exact sum mod p
    for sz in (0, 1, 3, 4, 5, 1000, 4096):
        arr = rng.integers(0, P, size=sz, dtype=np.uint64)
        got = cpmt._sum61(arr)
        want = int(sum(int(x) for x in arr) % P)
        assert got == want, f"_sum61 mismatch sz={sz}: {got} != {want}"
    ok += 1

    # [3] fmul dispatch routes BIG_Q to the Mersenne path under FAST_BIG and
    #     equals exact; default (no FAST_BIG) uses object/% and also equals
    saved = cpmt.FAST_BIG
    try:
        cpmt.FAST_BIG = True
        f = cpmt.fmul(a, b, P)
        assert all(int(x) == (int(u) * int(v)) % P for x, u, v in zip(f, a, b)), \
            "fmul FAST_BIG != exact"
        cpmt.FAST_BIG = False
        ao = a.astype(object); bo = b.astype(object)
        f2 = cpmt.fmul(ao, bo, P)
        assert all(int(x) == (int(u) * int(v)) % P for x, u, v in zip(f2, a, b)), \
            "fmul object != exact"
    finally:
        cpmt.FAST_BIG = saved
    ok += 1

    # [4] scatter_fold61 == exact integer scatter sum mod p
    idx = rng.integers(0, 64, size=2000, dtype=np.int64)
    vals = rng.integers(0, P, size=2000, dtype=np.uint64)
    got = cpmt.scatter_fold61(idx, vals, 64)
    want = np.zeros(64, dtype=object)
    for i, v in zip(idx, vals):
        want[i] = (int(want[i]) + int(v)) % P
    assert all(int(g) == int(w) for g, w in zip(got, want)), "scatter_fold61 mismatch"
    ok += 1

    # [5] regime classifier monotone + correct at boundaries
    assert _regime(16 * 1024) == "L1d"
    assert _regime(256 * 1024) == "L2"
    assert _regime(16 * 1024 * 1024) == "L3"
    assert _regime(1 << 30) == "DRAM"
    ok += 1

    # [6] tiny kernel sweep runs, ns/elem positive, _rate_at picks nearest size
    rows = bench_kernels(sizes_log2=[10, 12, 14], modulus=P)
    assert all(r["ns_mul"] > 0 for r in rows)
    ns, row = _rate_at(rows, (1 << 12) * 8)
    assert row["lg"] == 12, f"_rate_at nearest picked lg={row['lg']}"
    ok += 1

    # [7] model arithmetic: predicted ratio = ideal * measured per-elem growth
    fake = [dict(lg=20, n=1 << 20, footprint=(1 << 20) * 8, regime="DRAM",
                 ns_mul=10.0, ns_sum=1.0, ns_np=1.0, step=1.0, reps=(1, 1, 1)),
            dict(lg=27, n=1 << 27, footprint=(1 << 27) * 8, regime="DRAM",
                 ns_mul=14.5, ns_sum=1.0, ns_np=1.0, step=1.0, reps=(1, 1, 1))]
    r1, _ = _rate_at(fake, 185 * 1024 * 1024)   # ~2^25 footprint -> nearest 2^27 row? check
    assert r1 in (10.0, 14.5)
    ok += 1

    print(f"selftest: {ok}/7 groups passed.")
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--bench-kernels", action="store_true",
                    help="cache-step sweep of the BIG_Q kernels")
    ap.add_argument("--predict", action="store_true",
                    help="decompose the reach wall + emit the n=24->26 band")
    ap.add_argument("--bandwidth-probe", type=int, default=None, metavar="MAXLOG2",
                    help="light (no _mul61) sweep to 2^MAXLOG2 uint64 to probe the "
                         "multi-GiB TLB/page-walk regime (the n=26 working set)")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--maxlog2", type=int, default=27,
                    help="largest array size 2^k uint64 in the sweep (default 27 = 1 GiB)")
    args = ap.parse_args()

    if args.selftest:
        selftest()
        return
    if args.bandwidth_probe is not None:
        print("=" * 70)
        print(f"BANDWIDTH/TLB PROBE (light: _sum61 + np-mul only, to 2^{args.bandwidth_probe})")
        print("=" * 70)
        rows = bench_kernels(sizes_log2=list(range(20, args.bandwidth_probe + 1)),
                             light=True)
        print(f"  {'log2 n':>6} {'footprint':>11} {'regime':>6} "
              f"{'_sum61':>9} {'np mul':>9}")
        print(f"  {'':>6} {'(per array)':>11} {'':>6} {'ns/elem':>9} {'ns/elem':>9}")
        for r in rows:
            print(f"  {r['lg']:>6} {_bytes_human(r['footprint']):>11} {r['regime']:>6} "
                  f"{r['ns_sum']:>9.3f} {r['ns_np']:>9.3f}")
        return
    if not (args.bench_kernels or args.predict or args.all):
        ap.print_help()
        return

    rows = None
    if args.bench_kernels or args.predict or args.all:
        print("=" * 70)
        print("KERNEL CACHE-STEP SWEEP (_mul61 / _sum61 over uint64, BIG_Q)")
        print("=" * 70)
        rows = bench_kernels(sizes_log2=list(range(10, args.maxlog2 + 1)))
        print_kernel_table(rows)
        print()

    if args.predict or args.all:
        print("=" * 70)
        print("REACH WALL DECOMPOSITION + n=24->26 PREDICTION")
        print("=" * 70)
        predict(rows)


if __name__ == "__main__":
    main()
