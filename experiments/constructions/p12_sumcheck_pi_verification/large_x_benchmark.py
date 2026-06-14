#!/usr/bin/env python3
"""Large-x reach benchmark for the compressed succinct-verification chain
(open item 5).

This is the reproducible driver for the "large-x benchmark" headline: it runs
ONE honest verification of `pi(2^n - 1)` through the FULL succinct config of
`compressed_layer.run_chain` over the sound-characteristic field BIG_Q=2^61-1,
with the uint64 Mersenne fast path (FAST_BIG) enabled, and reports the
wall-clock plus the proof's structural numbers.

The full succinct config (identical to certificate_profile.py's FULL, S509):

    delegate=True       both wirings O(nb log p)            -> verifier Õ(sqrt x)
    structured=True     wiring prover O(nb p) (S496)        -> prover  Õ(x)
    pcs=True            real sum-check leaf openings (S505)
    batch_trace=True    K trace zero-tests -> 1 (S502)
    batch_ub=True       per-layer Ub opens -> 1 (S506)      -> per-layer leaf-free
    batch_wiring=True   K wiring chains -> 1 (S503/S504)     -> comm 5x down
    commit_base=True    S_0 base via tensor PCS (S508)       -> large-table ops Õ(x^{1/4})

PLUS, set on `compressed_prover_mult_trace.FAST_BIG`:

    FAST_BIG=True       uint64 Mersenne mulmod for 2^61-1 (S500/S504), a NET win
                        once BOTH big kernels are batched (the per-fmul array is
                        wide enough to amortise the 24-op mulmod).

This is NOT new protocol machinery -- every flag was built and cheat-tested in
S492-S512.  The contribution here is purely the REACH: running the bundle at
larger n than the S504 n=20 headline (66 s), recording the verified
`pi(2^n)==sieve_pi(2^n)` and the wall-clock, and confirming the proof's
size/verifier-op profile (comm dominated by the K sequential outer reductions
~Theta(sqrt x); per-layer verifier leaf-eval count exactly 0; large-table opens
only the one-time base commitment ~Theta(x^{1/4})) holds at the new reach.

Soundness at the headline n is spot-checked with ONE delta_pi liar (claim
pi(x)+1) -- it must be REJECTED.  The full cheat panel lives in
compressed_layer.py's selftest; here a single liar at the reach n is enough to
witness that the headline accept is not vacuous.

USAGE
  python3 large_x_benchmark.py --selftest          # fast correctness gate
  python3 large_x_benchmark.py --n 20              # reproduce the S504 headline
  python3 large_x_benchmark.py --n 22              # the reach push (~minutes)
  python3 large_x_benchmark.py --n 24 --no-cheat   # skip the liar check to halve wall

WHAT WOULD FALSIFY THE RESULT
  - claimed pi(x) != sieve_pi(x) at any n (the chain mis-counts);
  - the honest run REJECTS;
  - the delta_pi liar is ACCEPTED (soundness broken at the reach n);
  - FAST_BIG ON changing claimed/comm/accepted vs OFF (the fast path is not a
    verbatim drop-in -- selftest asserts bit-identical verdict);
  - any per-layer verifier large-table leaf-eval (vleaf_ops_pl > 0) reappearing
    (the per-layer verifier is supposed to be leaf-eval-free end-to-end).
"""
import argparse
import resource
import time

import numpy as np

import compressed_layer as cl
import compressed_prover_mult_trace as _cpmt


def _peak_rss_mb():
    """Process peak resident set (MB) from getrusage -- captures numpy buffer
    memory (OS-level RSS), unlike tracemalloc which sees only Python-object
    allocations.  ru_maxrss is in KiB on Linux."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0

# The full succinct config (matches certificate_profile.FULL, S509).
FULL = dict(delegate=True, structured=True, pcs=True, batch_trace=True,
            batch_ub=True, batch_wiring=True, commit_base=True)

# S525 memory-localization control: the FULL config WITHOUT the two batched
# discharges that build/hold the K-stacked witness cube (batch_trace builds
# `Ws = [build_witness(...) for p in primes]` and folds the stacked K*2^nb*Lv
# cube; batch_ub re-stacks it).  Everything else identical.  Comparing peak RSS
# of FULL vs NOBATCH isolates whether the chain's Theta(x)-ish working set comes
# from the batched-discharge prover cubes or from the per-layer reductions.
NOBATCH = dict(FULL, batch_trace=False, batch_ub=False)

# S527 LIST-streaming A/B control: the FULL config with stream_witnesses OFF, i.e.
# the chain materializes the held K-witness list Ws = [build_witness(...) for p in
# primes] before the batched discharges (the pre-S527 behaviour).  FULL streams it
# slice-by-slice (default stream_witnesses=True).  Comparing peak RSS of FULL vs
# FULL_NOSTREAM isolates the ~1.4 GB held-list term the streaming removes; both are
# bit-identical (same cube), so claimed/comm/accepted must match exactly.
FULL_NOSTREAM = dict(FULL, stream_witnesses=False)


def run_headline(n, seed=1, q=None, fast=True, check_cheat=True, verbose=True,
                 config=None):
    """Run ONE honest full-config verification of pi(2^n - 1) and (optionally)
    one delta_pi liar.  Returns a dict of the measured numbers.

    config overrides the chain flags (default FULL); used by the memory probe to
    run FULL vs NOBATCH.  fast toggles compressed_prover_mult_trace.FAST_BIG (the
    uint64 Mersenne path); it is restored on exit so callers are unaffected."""
    if q is None:
        q = cl.BIG_Q
    if config is None:
        config = FULL
    x = (1 << n) - 1
    V = cl.isqrt(x)
    nb = max(1, V.bit_length())
    K = len(cl.compressed_lucy(x)[0])
    truth = cl.sieve_pi(x)

    saved = _cpmt.FAST_BIG
    try:
        _cpmt.FAST_BIG = bool(fast) and (q == cl.BIG_Q)
        t0 = time.perf_counter()
        r = cl.run_chain(x, np.random.default_rng(seed), q=q, **config)
        wall = time.perf_counter() - t0
        ok = bool(r["accepted"]) and r["claimed"] == truth

        cheat_rejected = None
        cheat_wall = None
        if check_cheat:
            t1 = time.perf_counter()
            rc = cl.run_chain(x, np.random.default_rng(seed + 7),
                              cheat="delta_pi", q=q, **config)
            cheat_wall = time.perf_counter() - t1
            cheat_rejected = not rc["accepted"]
    finally:
        _cpmt.FAST_BIG = saved

    out = dict(n=n, x=x, V=V, nb=nb, K=K, truth=truth, peak_rss_mb=_peak_rss_mb(),
               claimed=r["claimed"], accepted=bool(r["accepted"]),
               match=(r["claimed"] == truth), ok=ok, wall=wall,
               comm=r["comm"], comm_outer=r["comm_outer"],
               comm_base=r["comm_base"], comm_bt=r["comm_bt"],
               comm_bw=r["comm_bw"], comm_ub=r["comm_ub"],
               vleaf_ops_pl=r["vleaf_ops_pl"], vleaf_ops_ot=r["vleaf_ops_ot"],
               vcommit_ops=r["vcommit_ops"],
               t_prover=r["t_prover"], t_verifier=r["t_verifier"],
               cheat_rejected=cheat_rejected, cheat_wall=cheat_wall,
               fast=_cpmt.FAST_BIG if False else (bool(fast) and q == cl.BIG_Q))

    if verbose:
        fld = "BIG_Q=2^61-1" + ("  [uint64 Mersenne FAST]" if out["fast"]
                                else "  [object dtype]") if q == cl.BIG_Q \
            else f"q={q}"
        print(f"n={n}: x = 2^{n}-1 = {x}, V = floor(sqrt x) = {V}, "
              f"nb = {nb} (cube 2^{nb} = {1 << nb}), layers K = {K}")
        print(f"  field = {fld}; config = FULL succinct "
              f"(delegate+structured+pcs+batch_trace+batch_ub+batch_wiring"
              f"+commit_base)")
        print(f"  HONEST: {'ACCEPTED' if out['accepted'] else 'REJECTED'}; "
              f"claimed pi(x) = {out['claimed']}, sieve = {truth}, "
              f"match = {out['match']}")
        print(f"  wall = {wall:.2f} s   "
              f"(t_prover {r['t_prover']*1000:.0f} ms, "
              f"t_verifier {r['t_verifier']*1000:.1f} ms reported region)")
        print(f"  peak RSS = {out['peak_rss_mb']:.0f} MB "
              f"({out['peak_rss_mb']/1024:.2f} GB, process-wide getrusage peak)")
        print(f"  certificate: comm = {out['comm']} field elems "
              f"(~{out['comm']*8/1024:.1f} KB @ 61-bit); "
              f"comm_outer = {out['comm_outer']} "
              f"({100*out['comm_outer']/max(1,out['comm']):.0f}% = the K "
              f"sequential reductions, Theta(sqrt x))")
        print(f"    comm_base={out['comm_base']} comm_bt={out['comm_bt']} "
              f"comm_bw={out['comm_bw']} comm_ub={out['comm_ub']}")
        print(f"  verifier large-table ops: per-layer leaf = {out['vleaf_ops_pl']} "
              f"(must be 0), one-time leaf = {out['vleaf_ops_ot']}, "
              f"base-commit opens = {out['vcommit_ops']} (Theta(x^1/4))")
        if check_cheat:
            print(f"  SOUNDNESS spot-check: delta_pi liar (claim pi+1) "
                  f"{'REJECTED' if out['cheat_rejected'] else 'ACCEPTED -- BUG!'}"
                  f"  ({cheat_wall:.2f} s)")
    return out


def selftest():
    """Fast correctness gate (small n).  Asserts, over BOTH Q and BIG_Q:
      1. the full config returns claimed pi(x) == sieve_pi(x), accepted;
      2. delta_pi and a self-consistent liar are both REJECTED;
      3. FAST_BIG ON gives a bit-identical verdict (claimed/comm/accepted) to
         FAST_BIG OFF over BIG_Q (the fast path is a verbatim drop-in);
      4. the per-layer verifier large-table leaf-eval count is EXACTLY 0
         (the Õ(sqrt x) end-to-end property), and comm_outer dominates comm.
    Boundary cases that surfaced during the build become cases here."""
    print("=== large_x_benchmark selftest ===")

    # 1+2: honest correctness + two cheat classes, both fields.
    for n in (8, 10, 12):
        x = (1 << n) - 1
        truth = cl.sieve_pi(x)
        K = len(cl.compressed_lucy(x)[0])
        for q in (cl.Q, cl.BIG_Q):
            r = cl.run_chain(x, np.random.default_rng(1), q=q, **FULL)
            assert r["accepted"] and r["claimed"] == truth, (n, q, r["claimed"], truth)
            rd = cl.run_chain(x, np.random.default_rng(2), cheat="delta_pi",
                              q=q, **FULL)
            assert not rd["accepted"], ("delta_pi accepted", n, q)
            i0 = max(1, K // 2)
            rl = cl.run_chain(x, np.random.default_rng(3), corrupt_layer=i0,
                              q=q, **FULL)
            assert not rl["accepted"], ("liar accepted", n, q, i0)
    print("  1. full config: claimed pi==sieve, accepted; delta_pi + "
          "self-consistent liar rejected (q & BIG_Q, n=8,10,12) OK")

    # 3: FAST_BIG bit-identical verdict (the headline runs with FAST_BIG=True;
    #    if it changed the transcript the headline would be untrustworthy).
    saved = _cpmt.FAST_BIG
    try:
        for n in (8, 10, 12):
            x = (1 << n) - 1
            _cpmt.FAST_BIG = False
            ro = cl.run_chain(x, np.random.default_rng(5), q=cl.BIG_Q, **FULL)
            _cpmt.FAST_BIG = True
            rf = cl.run_chain(x, np.random.default_rng(5), q=cl.BIG_Q, **FULL)
            assert (ro["claimed"] == rf["claimed"]
                    and ro["comm"] == rf["comm"]
                    and ro["accepted"] == rf["accepted"]), (
                "FAST_BIG not a verbatim drop-in", n,
                (ro["claimed"], ro["comm"], ro["accepted"]),
                (rf["claimed"], rf["comm"], rf["accepted"]))
            # and the headline driver's fast run must still reject a liar
            _cpmt.FAST_BIG = True
            rc = cl.run_chain(x, np.random.default_rng(6), cheat="delta_pi",
                              q=cl.BIG_Q, **FULL)
            assert not rc["accepted"], ("FAST delta_pi accepted", n)
    finally:
        _cpmt.FAST_BIG = saved
    print("  2. FAST_BIG ON == OFF (claimed/comm/accepted bit-identical over "
          "BIG_Q); FAST liar still rejected OK")

    # 4: the structural milestone numbers -- per-layer verifier leaf-eval count
    #    is 0, comm_outer dominates.
    for n in (10, 12):
        x = (1 << n) - 1
        r = cl.run_chain(x, np.random.default_rng(1), q=cl.BIG_Q, **FULL)
        assert r["vleaf_ops_pl"] == 0, ("per-layer leaf ops nonzero", n,
                                        r["vleaf_ops_pl"])
        assert r["comm_outer"] > 0.5 * r["comm"], ("outer not dominant", n,
                                                   r["comm_outer"], r["comm"])
    print("  3. per-layer verifier leaf-eval count == 0; comm_outer dominant "
          "(Õ(sqrt x) end-to-end) OK")

    # 5: the driver wrapper itself produces a correct, accepted, sound headline
    #    at a small n (exercises run_headline + the delta_pi spot-check path).
    h = run_headline(10, seed=1, fast=True, check_cheat=True, verbose=False)
    assert h["ok"] and h["match"] and h["cheat_rejected"], h
    assert h["vleaf_ops_pl"] == 0, h
    print("  4. run_headline wrapper: ok + match + cheat rejected + leaf-free OK")

    # 6 (S525): the --no-scatter-fold A/B control is a verbatim drop-in.  With
    #    FAST_BIG on over BIG_Q, USE_SCATTER_FOLD=False must take the ORIGINAL
    #    object scatter (NOT the uint64 _dt(q) path, which overflows and was the
    #    bug that produced a spurious REJECT), giving claimed/comm/accepted
    #    bit-identical to the scatter_fold61 path -- and still ACCEPT the honest
    #    run.  Guards the control so the A/B numbers are trustworthy.
    saved_sf = cl.USE_SCATTER_FOLD
    saved_fb = _cpmt.FAST_BIG
    try:
        _cpmt.FAST_BIG = True
        for n in (8, 10, 12):
            x = (1 << n) - 1
            truth = cl.sieve_pi(x)
            cl.USE_SCATTER_FOLD = True
            rs = cl.run_chain(x, np.random.default_rng(7), q=cl.BIG_Q, **FULL)
            cl.USE_SCATTER_FOLD = False
            ro = cl.run_chain(x, np.random.default_rng(7), q=cl.BIG_Q, **FULL)
            assert rs["accepted"] and rs["claimed"] == truth, ("fold reject", n)
            assert ro["accepted"] and ro["claimed"] == truth, ("object reject", n)
            assert (rs["claimed"] == ro["claimed"]
                    and rs["comm"] == ro["comm"]
                    and rs["accepted"] == ro["accepted"]), (
                "scatter_fold61 != object A/B control", n,
                (rs["claimed"], rs["comm"]), (ro["claimed"], ro["comm"]))
    finally:
        cl.USE_SCATTER_FOLD = saved_sf
        _cpmt.FAST_BIG = saved_fb
    print("  5. --no-scatter-fold control == scatter_fold61 (claimed/comm/"
          "accepted bit-identical over BIG_Q, both accept honest) OK")

    # 7 (S527): stream_witnesses (default True, FULL) is a verbatim drop-in --
    #    claimed/comm/accepted BIT-IDENTICAL to stream_witnesses=False
    #    (FULL_NOSTREAM), over BOTH Q and BIG_Q, with FAST_BIG on for BIG_Q.  This
    #    is what makes the landed n=24 artifact reproduce exactly with lower peak
    #    RSS; if streaming changed the transcript the reach numbers would shift.
    saved_fb = _cpmt.FAST_BIG
    try:
        for n in (8, 10, 12):
            x = (1 << n) - 1
            truth = cl.sieve_pi(x)
            for q in (cl.Q, cl.BIG_Q):
                _cpmt.FAST_BIG = (q == cl.BIG_Q)
                rs = cl.run_chain(x, np.random.default_rng(9), q=q, **FULL)
                rn = cl.run_chain(x, np.random.default_rng(9), q=q,
                                  **FULL_NOSTREAM)
                assert rs["accepted"] and rs["claimed"] == truth, ("stream", n, q)
                assert rn["accepted"] and rn["claimed"] == truth, ("nostream", n, q)
                assert (rs["claimed"] == rn["claimed"]
                        and rs["comm"] == rn["comm"]
                        and rs["accepted"] == rn["accepted"]
                        and rs["comm_bt"] == rn["comm_bt"]
                        and rs["comm_ub"] == rn["comm_ub"]), (
                    "stream_witnesses not a verbatim drop-in", n, q,
                    (rs["claimed"], rs["comm"]), (rn["claimed"], rn["comm"]))
                # the delta_pi liar is still rejected under streaming
                rd = cl.run_chain(x, np.random.default_rng(10), cheat="delta_pi",
                                  q=q, **FULL)
                assert not rd["accepted"], ("stream delta_pi accepted", n, q)
    finally:
        _cpmt.FAST_BIG = saved_fb
    print("  6. stream_witnesses ON == OFF (claimed/comm/comm_bt/comm_ub/accepted "
          "bit-identical over Q & BIG_Q; streamed liar still rejected) OK")

    print("ALL SELFTESTS PASSED")


CONFIGS = {"full": FULL, "nobatch": NOBATCH, "full_nostream": FULL_NOSTREAM}


def mem_probe(ns, seed=1, configs=("full", "nobatch")):
    """Localize the chain's peak-RSS source.  For each n, run each config in a
    FRESH subprocess (ru_maxrss is a process-wide monotonic peak, so configs must
    NOT share a process) and report peak RSS + wall.

    Default (full vs nobatch): if FULL peak >> NOBATCH peak, the Theta(x polylog)
    working set is the batched-discharge prover cubes (the held K-witness list +
    the stacked sum-check), NOT the per-layer reductions -- the corrected S524
    attribution.  S527 A/B (full vs full_nostream): the gap is exactly the held
    K-witness LIST that stream_witnesses removes (claimed/comm identical)."""
    import subprocess
    import sys
    print(f"=== memory-localization probe ({' vs '.join(configs)}, "
          f"fresh process each) ===")
    print(f"{'n':>3} {'K':>5} {'config':>13} {'peak_MB':>9} {'peak_GB':>8} "
          f"{'wall_s':>8} {'claimed':>10}")
    rows = []
    for n in ns:
        for cfg in configs:
            out = subprocess.run(
                [sys.executable, __file__, "--mem-one", str(n),
                 "--config", cfg, "--seed", str(seed)],
                capture_output=True, text=True)
            line = [l for l in out.stdout.splitlines()
                    if l.startswith("MEMONE ")]
            if not line:
                print(f"  n={n} {cfg}: FAILED\n{out.stdout}\n{out.stderr}")
                continue
            d = dict(kv.split("=") for kv in line[0].split()[1:])
            rows.append((n, cfg, d))
            print(f"{n:>3} {d['K']:>5} {cfg:>13} {float(d['peak_mb']):>9.0f} "
                  f"{float(d['peak_mb'])/1024:>8.2f} {float(d['wall']):>8.1f} "
                  f"{d['claimed']:>10}")
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=20,
                    help="run the full succinct config at x = 2^n - 1 (default 20)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--field", choices=["big", "q"], default="big",
                    help="big = BIG_Q=2^61-1 (sound to n~60); q = 2^31-1 "
                         "(demo, sound only for n<~30)")
    ap.add_argument("--no-fast", action="store_true",
                    help="disable the uint64 Mersenne path (object dtype)")
    ap.add_argument("--no-cheat", action="store_true",
                    help="skip the delta_pi soundness spot-check (halves wall)")
    ap.add_argument("--no-scatter-fold", action="store_true",
                    help="A/B control: force the object-dtype np.add.at for the "
                         "two per-layer outer-reduction scatter-sums even on the "
                         "fast path (isolates scatter_fold61's effect, S525)")
    ap.add_argument("--mem-probe", type=str, default=None,
                    help="comma-separated n list: localize peak RSS (FULL vs "
                         "NOBATCH, fresh process each), e.g. --mem-probe 16,18,20")
    ap.add_argument("--stream-probe", type=str, default=None,
                    help="comma-separated n list: S527 LIST-streaming A/B peak RSS "
                         "(FULL vs FULL_NOSTREAM), e.g. --stream-probe 18,20,22")
    ap.add_argument("--mem-one", type=int, default=None,
                    help="(internal) run ONE honest config and print MEMONE line")
    ap.add_argument("--config", choices=list(CONFIGS), default="full",
                    help="chain config for --mem-one (full, nobatch, full_nostream)")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    cl.USE_SCATTER_FOLD = not args.no_scatter_fold
    if args.selftest:
        selftest()
        return
    if args.mem_one is not None:
        cfg = CONFIGS[args.config]
        n = args.mem_one
        h = run_headline(n, seed=args.seed, q=cl.BIG_Q, fast=True,
                         check_cheat=False, verbose=False, config=cfg)
        print(f"MEMONE n={n} K={h['K']} config={args.config} "
              f"peak_mb={h['peak_rss_mb']:.1f} wall={h['wall']:.2f} "
              f"claimed={h['claimed']} accepted={int(h['accepted'])}")
        return
    if args.mem_probe is not None:
        ns = [int(t) for t in args.mem_probe.split(",") if t.strip()]
        mem_probe(ns, seed=args.seed)
        return
    if args.stream_probe is not None:
        ns = [int(t) for t in args.stream_probe.split(",") if t.strip()]
        mem_probe(ns, seed=args.seed, configs=("full", "full_nostream"))
        return
    q = cl.BIG_Q if args.field == "big" else cl.Q
    run_headline(args.n, seed=args.seed, q=q, fast=not args.no_fast,
                 check_cheat=not args.no_cheat, verbose=True)


if __name__ == "__main__":
    main()
