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
import time

import numpy as np

import compressed_layer as cl
import compressed_prover_mult_trace as _cpmt

# The full succinct config (matches certificate_profile.FULL, S509).
FULL = dict(delegate=True, structured=True, pcs=True, batch_trace=True,
            batch_ub=True, batch_wiring=True, commit_base=True)


def run_headline(n, seed=1, q=None, fast=True, check_cheat=True, verbose=True):
    """Run ONE honest full-config verification of pi(2^n - 1) and (optionally)
    one delta_pi liar.  Returns a dict of the measured numbers.

    fast toggles compressed_prover_mult_trace.FAST_BIG (the uint64 Mersenne path);
    it is restored on exit so callers/selftests are unaffected."""
    if q is None:
        q = cl.BIG_Q
    x = (1 << n) - 1
    V = cl.isqrt(x)
    nb = max(1, V.bit_length())
    K = len(cl.compressed_lucy(x)[0])
    truth = cl.sieve_pi(x)

    saved = _cpmt.FAST_BIG
    try:
        _cpmt.FAST_BIG = bool(fast) and (q == cl.BIG_Q)
        t0 = time.perf_counter()
        r = cl.run_chain(x, np.random.default_rng(seed), q=q, **FULL)
        wall = time.perf_counter() - t0
        ok = bool(r["accepted"]) and r["claimed"] == truth

        cheat_rejected = None
        cheat_wall = None
        if check_cheat:
            t1 = time.perf_counter()
            rc = cl.run_chain(x, np.random.default_rng(seed + 7),
                              cheat="delta_pi", q=q, **FULL)
            cheat_wall = time.perf_counter() - t1
            cheat_rejected = not rc["accepted"]
    finally:
        _cpmt.FAST_BIG = saved

    out = dict(n=n, x=x, V=V, nb=nb, K=K, truth=truth,
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

    print("ALL SELFTESTS PASSED")


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
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return
    q = cl.BIG_Q if args.field == "big" else cl.Q
    run_headline(args.n, seed=args.seed, q=q, fast=not args.no_fast,
                 check_cheat=not args.no_cheat, verbose=True)


if __name__ == "__main__":
    main()
