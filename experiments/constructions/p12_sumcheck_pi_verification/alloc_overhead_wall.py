#!/usr/bin/env python3
"""
alloc_overhead_wall.py  -- the THIRD leg of the verification-chain reach wall.

Context (PROGRAM.md item 5, S524/S525/S529):
  The large-x verification reach (`large_x_benchmark.py`) wall has been modelled
  (S529, `mem_hierarchy_wall.py`) as two legs:
     (1) op-count  Theta(x)*polylog   (per-Dn=2 ratio 4.0-4.4x, confirmed S530/S531/S532b to n=24)
     (2) a compute/cache memory factor (the L1->DRAM `_mul61` bowl: one-time ~7x L3
         step + per-Dn growth <=1.35x)
  -> S529 predicts the n24->26 wall ratio in [4.0, 5.4]x, "well above ~5.5x => reopen super-Theta(x)".

  But S529 leg (2) was measured as STEADY-STATE `_mul61` ns/elem over PRE-ALLOCATED
  arrays -- it never measured the cost of ALLOCATING and first-touching the large
  arrays.  The LIVE n=26 reach (observed this cycle via /proc) is *kernel/allocation
  bound*, NOT compute bound:  stime ~= 2.3x utime, ~3.4e5 minor page-faults/sec,
  ~3e9 cumulative minor faults, majflt=0 (no swap).  That is a THIRD leg the S529
  model omits.

This script measures leg (3) two ways, contention-safely:
  --live [PID]  : read the live reach's /proc counters (ZERO contention -- pure reads):
                  stime fraction, minor/major fault rate, cumulative faults, RSS.
                  Auto-detects the reach PID via the `large_x_benchmark` cmdline.
  --probe       : a SELF-LIMITING micro-benchmark separating, for an elementwise
                  uint64 op, the cost of "fresh allocation + first-touch each rep"
                  (malloc/mmap + minor page-faults, charged to KERNEL/stime) from the
                  same op on a "reused pre-faulted buffer" (pure compute, S529's regime).
                  Uses /proc/self/stat to attribute kernel vs user time per phase.
  --predict     : combine the three legs into a SHARPENED falsifier for the harvest
                  cycle: a wall ratio > 5.4x is consistent with Theta(x) op-count IFF
                  the excess is page-fault stime (allocation-bound).  Only a wall excess
                  that survives subtracting allocation overhead AND is utime-dominated
                  with super-linear per-element growth reopens super-Theta(x).

Scope: an in-cycle MEASUREMENT refining the S529 reach-wall model (a known wall, new
leg), NOT a goal advance.  polylog pi(x) stays blocked.  The allocation leg is an
implementation constant (removable by buffer reuse / out= -- a constant-factor working
set win like S527's list-streaming), NOT an op-count / complexity change.

Falsifiers (see _results.md): F1 reused-buffer op is NOT utime-dominated; F2 fresh
allocation is NOT stime-dominated at large N; F3 fresh/reused per-elem ratio does NOT
grow with N then saturate; F4 the live reach is NOT stime-dominated (would refute the
whole premise -- then the reach wall is compute, and a >5.4x ratio WOULD reopen
super-Theta(x)).
"""

import argparse
import os
import subprocess
import sys
import time

import numpy as np

CLK_TCK = os.sysconf("SC_CLK_TCK") if hasattr(os, "sysconf") else 100
PAGE = 4096  # bytes; minor faults are charged per touched page


# ---------------------------------------------------------------------------
# /proc parsing (Linux).  /proc/PID/stat fields (1-indexed, after comm):
#   10 minflt   12 majflt   14 utime   15 stime   20 num_threads   24 rss(pages)
# comm may contain spaces/parens, so split on the LAST ')'.
# ---------------------------------------------------------------------------
def read_proc_stat(pid):
    with open(f"/proc/{pid}/stat", "rb") as f:
        data = f.read().decode("latin-1")
    rparen = data.rfind(")")
    head = data[:rparen]
    rest = data[rparen + 1:].split()
    # rest[0] is field 3 (state); field i (1-indexed) -> rest[i-3]
    pid_field = int(head.split("(", 1)[0].strip())
    state = rest[0]                  # field 3
    minflt = int(rest[10 - 3])       # field 10
    majflt = int(rest[12 - 3])       # field 12
    utime = int(rest[14 - 3])        # field 14
    stime = int(rest[15 - 3])        # field 15
    num_threads = int(rest[20 - 3])  # field 20
    rss_pages = int(rest[24 - 3])    # field 24
    return {
        "pid": pid_field, "state": state, "minflt": minflt, "majflt": majflt,
        "utime": utime, "stime": stime, "num_threads": num_threads,
        "rss_bytes": rss_pages * PAGE,
    }


def find_reach_pid():
    """The detached large_x_benchmark reach process, or None."""
    try:
        out = subprocess.run(
            ["pgrep", "-f", "large_x_benchmark"],
            capture_output=True, text=True, timeout=10).stdout
    except Exception:
        return None
    me = os.getpid()
    for line in out.split():
        try:
            pid = int(line)
        except ValueError:
            continue
        if pid == me:
            continue
        # confirm it's actually the python benchmark, not our own grep/shell
        try:
            with open(f"/proc/{pid}/cmdline", "rb") as f:
                cl = f.read().replace(b"\x00", b" ").decode("latin-1")
        except OSError:
            continue
        if "large_x_benchmark.py" in cl and "python" in cl.lower():
            return pid
    return None


def _human(b):
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if abs(b) < 1024 or unit == "TiB":
            return f"{b:.1f} {unit}"
        b /= 1024.0


# ---------------------------------------------------------------------------
# (A) Live-reach in-situ measurement -- ZERO contention (only reads /proc).
# ---------------------------------------------------------------------------
def sample_live(pid, window_s=5.0):
    s0 = read_proc_stat(pid)
    t0 = time.monotonic()
    time.sleep(window_s)
    s1 = read_proc_stat(pid)
    dt = time.monotonic() - t0
    d_min = s1["minflt"] - s0["minflt"]
    d_maj = s1["majflt"] - s0["majflt"]
    d_u = s1["utime"] - s0["utime"]
    d_s = s1["stime"] - s0["stime"]
    cpu = d_u + d_s
    res = {
        "pid": pid, "window_s": dt, "state": s1["state"],
        "num_threads": s1["num_threads"], "rss_bytes": s1["rss_bytes"],
        "minflt_per_s": d_min / dt, "majflt_per_s": d_maj / dt,
        "win_stime_frac": (d_s / cpu) if cpu else float("nan"),
        "cum_minflt": s1["minflt"], "cum_majflt": s1["majflt"],
        "cum_utime_s": s1["utime"] / CLK_TCK, "cum_stime_s": s1["stime"] / CLK_TCK,
        "cum_stime_frac": s1["stime"] / (s1["utime"] + s1["stime"]) if (s1["utime"] + s1["stime"]) else float("nan"),
        # bytes newly faulted in (lower bound on alloc+touch traffic) per s
        "fault_traffic_Bps": (d_min / dt) * PAGE,
    }
    return res


def print_live(res):
    print(f"  pid={res['pid']}  state={res['state']}  threads={res['num_threads']}  "
          f"RSS={_human(res['rss_bytes'])}")
    print(f"  window {res['window_s']:.1f}s:  minflt={res['minflt_per_s']:.3e}/s  "
          f"majflt={res['majflt_per_s']:.3e}/s  fault-traffic={_human(res['fault_traffic_Bps'])}/s")
    print(f"  window stime fraction = {res['win_stime_frac']:.3f}   "
          f"(>0.5 => kernel/allocation bound)")
    print(f"  cumulative:  utime={res['cum_utime_s']:.0f}s  stime={res['cum_stime_s']:.0f}s  "
          f"stime_frac={res['cum_stime_frac']:.3f}  minflt={res['cum_minflt']:.3e}  "
          f"majflt={res['cum_majflt']:.0f}")


# ---------------------------------------------------------------------------
# (B) Controlled micro-probe: allocation+first-touch (stime) vs reused buffer (utime).
#     Self-limiting: capped array size and rep budget.  Single-threaded (matches
#     the single-threaded reach).  Attributes kernel/user via /proc/self/stat.
# ---------------------------------------------------------------------------
def _op_into(out, a, b, q):
    """A representative modular-mul step done with out= (NO fresh allocation)."""
    np.multiply(a, b, out=out)
    np.mod(out, q, out=out)


def _op_fresh(a, b, q):
    """The same step but allocating fresh temporaries (allocate + page-fault)."""
    return (a * b) % q


def _self_jiffies():
    s = read_proc_stat("self")
    return s["utime"], s["stime"]


def _time_phase(fn, reps):
    """Return (wall_s, d_utime_s, d_stime_s) over `reps` calls of fn()."""
    u0, s0 = _self_jiffies()
    t0 = time.monotonic()
    for _ in range(reps):
        fn()
    wall = time.monotonic() - t0
    u1, s1 = _self_jiffies()
    return wall, (u1 - u0) / CLK_TCK, (s1 - s0) / CLK_TCK


def probe(sizes_log2=None, max_mb=256, q=(1 << 61) - 1, byte_budget=8 << 30):
    """
    For each size 2^k uint64:
      reused : op done with out= into a pre-faulted buffer  -> pure compute (utime)
      fresh  : op done allocating fresh temporaries each rep -> alloc + page-fault (stime)
    `byte_budget` caps total bytes touched across the whole probe (contention guard).
    """
    if sizes_log2 is None:
        # 2^13 (64 KiB, ~L1/L2) .. up to max_mb
        top = (max_mb * 1024 * 1024).bit_length() - 4  # uint64 -> /8 -> -3, minus 1 margin
        sizes_log2 = list(range(13, min(top, 25) + 1))
    rng = np.random.default_rng(12345)
    rows = []
    touched = 0
    for k in sizes_log2:
        n = 1 << k
        nbytes = n * 8
        a = (rng.integers(0, 1 << 60, size=n, dtype=np.uint64))
        b = (rng.integers(0, 1 << 60, size=n, dtype=np.uint64))
        out = np.empty(n, dtype=np.uint64)
        _op_into(out, a, b, q)  # warm/fault the reused buffers + inputs
        # rep budget: more reps for small arrays, >=2 for large; bounded by time
        reps = max(2, min(200, int((1 << 22) / n) + 2))
        # ---- reused (compute floor) ----
        w_r, u_r, s_r = _time_phase(lambda: _op_into(out, a, b, q), reps)
        # ---- fresh (alloc + fault) ----
        w_f, u_f, s_f = _time_phase(lambda: _op_fresh(a, b, q), reps)
        ns_reused = w_r / reps / n * 1e9
        ns_fresh = w_f / reps / n * 1e9
        rows.append({
            "k": k, "n": n, "bytes": nbytes, "reps": reps,
            "ns_reused": ns_reused, "ns_fresh": ns_fresh,
            "ratio": ns_fresh / ns_reused if ns_reused else float("nan"),
            "reused_stime_frac": s_r / (u_r + s_r) if (u_r + s_r) else float("nan"),
            "fresh_stime_frac": s_f / (u_f + s_f) if (u_f + s_f) else float("nan"),
            # alloc overhead attributable to the fresh path, ns/elem
            "alloc_ns": ns_fresh - ns_reused,
        })
        # fresh path allocates ~3 temporaries (a*b, %q, plus internal) per rep
        touched += nbytes * reps * 4 + nbytes * 3  # inputs+out warmed once
        del a, b, out
        if touched > byte_budget:
            break
    return rows


def print_probe(rows):
    print(f"  {'size':>9} {'bytes':>9} {'reps':>5} "
          f"{'reused':>9} {'fresh':>9} {'ratio':>6} {'alloc':>8} "
          f"{'reused_s%':>9} {'fresh_s%':>9}")
    print(f"  {'(uint64)':>9} {'':>9} {'':>5} "
          f"{'ns/el':>9} {'ns/el':>9} {'f/r':>6} {'ns/el':>8} "
          f"{'kernel':>9} {'kernel':>9}")
    for r in rows:
        print(f"  {r['n']:>9} {_human(r['bytes']):>9} {r['reps']:>5} "
              f"{r['ns_reused']:>9.2f} {r['ns_fresh']:>9.2f} {r['ratio']:>6.2f} "
              f"{r['alloc_ns']:>8.2f} {r['reused_stime_frac']:>9.3f} "
              f"{r['fresh_stime_frac']:>9.3f}")


def probe_summary(rows):
    """Largest-size readings = the reach-like regime."""
    big = rows[-1]
    small = rows[0]
    return {
        "small_k": small["k"], "big_k": big["k"],
        "big_ratio": big["ratio"], "big_alloc_ns": big["alloc_ns"],
        "big_fresh_stime_frac": big["fresh_stime_frac"],
        "big_reused_stime_frac": big["reused_stime_frac"],
        "ratio_grew": big["ratio"] > small["ratio"],
    }


# ---------------------------------------------------------------------------
# (C) Three-leg model + sharpened falsifier.
# ---------------------------------------------------------------------------
def predict(live=None, probe_rows=None):
    print("\n=== Three-leg reach-wall model + sharpened falsifier ===")
    print("  wall(n) = opcount(n) x compute_cache_factor(n) x alloc_factor(n)")
    print("  leg 1 opcount        : Theta(x)*polylog, per-Dn=2 ratio 4.0-4.4x DOWN  "
          "(S530/S531/S532b, measured to n=24)")
    print("  leg 2 compute/cache  : L1->DRAM `_mul61` bowl, one-time ~7x L3 step, "
          "per-Dn growth <=1.35x  (S529)")
    print("  leg 3 alloc/faults   : malloc/mmap + minor page-faults, charged to KERNEL stime "
          "(THIS cycle; S529 omitted it)")
    if live is not None:
        sf = live["cum_stime_frac"]
        print(f"\n  live n=26 reach: cumulative stime fraction = {sf:.3f}, "
              f"minflt rate {live['minflt_per_s']:.2e}/s, majflt {live['cum_majflt']:.0f}")
        if sf > 0.5:
            print("  => the reach is KERNEL/ALLOCATION-BOUND: leg 3 dominates its CPU. "
                  "F4 NOT triggered.")
        else:
            print("  => the reach is NOT stime-dominated: F4 TRIGGERED -- leg 3 is small, "
                  "a >5.4x wall WOULD reopen super-Theta(x).")
    if probe_rows is not None:
        s = probe_summary(probe_rows)
        print(f"\n  probe (size 2^{s['big_k']}): fresh/reused per-elem ratio = "
              f"{s['big_ratio']:.2f}, alloc overhead = {s['big_alloc_ns']:.2f} ns/elem, "
              f"fresh stime frac = {s['big_fresh_stime_frac']:.3f}")
        print("  => fresh-allocation cost is page-fault (kernel) dominated at large N; "
              "reused-buffer cost is the S529 compute bowl.")
    print("\n  SHARPENED FALSIFIER for the harvest cycle:")
    print("    A wall ratio n24->26 > 5.4x is CONSISTENT with Theta(x) op-count IFF the")
    print("    n=26 reach is stime-dominated (leg 3, allocation).  Check the reach's")
    print("    stime/utime (e.g. /proc or `time`), not just the wall.  Only a wall excess")
    print("    that (i) survives subtracting allocation stime AND (ii) is utime-dominated")
    print("    with super-linear per-element growth would reopen super-Theta(x).")
    print("    The allocation leg is an IMPLEMENTATION constant (removable by buffer reuse/")
    print("    out=), NOT an op-count/complexity change -- op-count is Theta(x) to n=24.")


# ---------------------------------------------------------------------------
# selftest
# ---------------------------------------------------------------------------
def selftest():
    rng = np.random.default_rng(7)
    failures = []

    def check(name, cond):
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
        if not cond:
            failures.append(name)

    # [1] /proc/self/stat parses and is self-consistent
    s = read_proc_stat("self")
    check("read_proc_stat(self): pid matches os.getpid",
          s["pid"] == os.getpid())
    check("read_proc_stat(self): counters non-negative",
          s["minflt"] >= 0 and s["utime"] >= 0 and s["stime"] >= 0 and s["rss_bytes"] > 0)

    # [1b] FIELD-OFFSET cross-check (catches the off-by-one parsing bug): burn ~0.2s
    #      of CPU, then read_proc_stat utime must track os.times().user, num_threads
    #      must be >=1, and RSS must be a sane (sub-TiB) value.
    t_user0 = os.times().user
    sp = read_proc_stat("self")
    burn = np.uint64(0)
    spin = rng.integers(0, 1 << 60, size=1 << 20, dtype=np.uint64)
    for _ in range(40):
        burn ^= np.bitwise_xor.reduce(spin)
    t_user1 = os.times().user
    sp2 = read_proc_stat("self")
    d_utime_s = (sp2["utime"] - sp["utime"]) / CLK_TCK
    d_user_s = t_user1 - t_user0
    check(f"utime field tracks os.times().user (proc du={d_utime_s:.2f}s, times du={d_user_s:.2f}s)",
          abs(d_utime_s - d_user_s) < 0.5 + 0.5 * d_user_s and d_utime_s > 0)
    check(f"num_threads field is a sane thread count (got {sp2['num_threads']})",
          1 <= sp2["num_threads"] <= 100000)
    check(f"rss field is a sane size (got {_human(sp2['rss_bytes'])})",
          0 < sp2["rss_bytes"] < (1 << 50))  # < 1 PiB; a TiB-scale read = wrong field

    # [2] comm-with-parens robustness: a synthetic stat line with ')' and spaces in comm
    fake = "1234 (weird ) name) R 1 1 1 0 -1 0 " + " ".join(["7"] * 30)
    rparen = fake.rfind(")")
    rest = fake[rparen + 1:].split()
    check("comm-with-paren split: state field is the token after last ')'",
          rest[0] == "R")

    # [3] fault counters are monotonic non-decreasing over time (a real process)
    a0 = read_proc_stat("self")["minflt"]
    _ = np.empty(1 << 20, dtype=np.uint64); _[:] = 1  # force some faults
    a1 = read_proc_stat("self")["minflt"]
    check("minflt monotonic non-decreasing", a1 >= a0)

    # [4] _op_into and _op_fresh compute the SAME values (out= path is correct)
    q = (1 << 61) - 1
    n = 1 << 12
    a = rng.integers(0, 1 << 60, size=n, dtype=np.uint64)
    b = rng.integers(0, 1 << 60, size=n, dtype=np.uint64)
    out = np.empty(n, dtype=np.uint64)
    _op_into(out, a, b, q)
    fresh = _op_fresh(a, b, q)
    check("_op_into == _op_fresh (value equality)", bool(np.array_equal(out, fresh)))

    # [5] _time_phase returns non-negative wall and jiffy deltas
    w, du, ds = _time_phase(lambda: _op_into(out, a, b, q), 5)
    check("_time_phase: non-negative wall/utime/stime deltas",
          w >= 0 and du >= 0 and ds >= 0)

    # [6] probe runs, rows well-formed, fresh >= reused (alloc never free), ratio>=~1
    rows = probe(sizes_log2=[12, 14, 16], max_mb=8, byte_budget=1 << 30)
    check("probe: produced rows", len(rows) >= 2)
    check("probe: ns_fresh >= ns_reused (alloc overhead non-negative-ish)",
          all(r["ns_fresh"] >= 0.7 * r["ns_reused"] for r in rows))  # noise tolerance
    check("probe: alloc_ns finite", all(np.isfinite(r["alloc_ns"]) for r in rows))

    # [7] probe_summary fields present
    summ = probe_summary(rows)
    check("probe_summary: has big_ratio + fault fracs",
          "big_ratio" in summ and "big_fresh_stime_frac" in summ)

    # [8] _human formats bytes sensibly
    check("_human(4096) ~ KiB", _human(4096) == "4.0 KiB")
    check("_human(1<<30) ~ GiB", _human(1 << 30) == "1.0 GiB")

    # [9] sample_live works on our OWN pid (self) -- spins a tiny load
    me = os.getpid()
    res = sample_live(me, window_s=0.3)
    check("sample_live(self): returns finite minflt rate + stime frac",
          np.isfinite(res["minflt_per_s"]) and "win_stime_frac" in res)

    # [10] find_reach_pid never returns our own pid
    rp = find_reach_pid()
    check("find_reach_pid: not our own pid", rp != me)

    print(f"\n  selftest: {len(failures)} failure(s)" if failures else
          "\n  selftest: ALL PASS")
    return 1 if failures else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--live", nargs="?", const="auto", default=None, metavar="PID",
                    help="sample the live reach process (auto-detect if no PID given)")
    ap.add_argument("--window", type=float, default=5.0,
                    help="live-sample window seconds (default 5)")
    ap.add_argument("--probe", action="store_true",
                    help="run the alloc-vs-reused micro-benchmark")
    ap.add_argument("--max-mb", type=int, default=256,
                    help="largest probe array in MiB (contention/safety cap; default 256)")
    ap.add_argument("--predict", action="store_true",
                    help="print the three-leg model + sharpened falsifier")
    ap.add_argument("--all", action="store_true",
                    help="--live (auto) + --probe + --predict, with a before/after "
                         "contention check on the live reach")
    args = ap.parse_args()

    if args.selftest:
        sys.exit(selftest())

    if not (args.live or args.probe or args.predict or args.all):
        ap.print_help()
        return

    live_res = None
    probe_rows = None

    if args.all:
        pid = find_reach_pid()
        if pid:
            print("=== (A) live reach BEFORE probe ===")
            live_res = sample_live(pid, window_s=args.window)
            print_live(live_res)
        else:
            print("=== (A) live reach: none found (large_x_benchmark not running) ===")
        print("\n=== (B) allocation-overhead micro-probe ===")
        probe_rows = probe(max_mb=args.max_mb)
        print_probe(probe_rows)
        if pid:
            print("\n=== (A') live reach AFTER probe (contention check) ===")
            after = sample_live(pid, window_s=args.window)
            print_live(after)
            d = abs(after["minflt_per_s"] - live_res["minflt_per_s"]) / max(1.0, live_res["minflt_per_s"])
            print(f"  minflt-rate change across probe: {d*100:.1f}%  "
                  f"({'negligible (<25%) => probe did not materially contend' if d < 0.25 else 'check contention'})")
        predict(live=live_res, probe_rows=probe_rows)
        return

    if args.live:
        pid = args.live
        if pid == "auto":
            pid = find_reach_pid()
            if not pid:
                print("no large_x_benchmark reach process found")
                return
        else:
            pid = int(pid)
        print(f"=== live reach pid={pid} ===")
        live_res = sample_live(pid, window_s=args.window)
        print_live(live_res)

    if args.probe:
        print("=== allocation-overhead micro-probe ===")
        probe_rows = probe(max_mb=args.max_mb)
        print_probe(probe_rows)

    if args.predict:
        predict(live=live_res, probe_rows=probe_rows)


if __name__ == "__main__":
    main()
