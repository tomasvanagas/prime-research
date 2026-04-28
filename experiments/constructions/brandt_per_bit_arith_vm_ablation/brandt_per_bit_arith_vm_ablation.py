#!/usr/bin/env python3
"""
C3.a.iv -- Arithmetic-primitive ablation of the C3.a bounded-Kt VM.

Re-runs S150's per-bit Kt_b' scan under six ablation conditions:

    baseline         all 16 ops enabled (sanity check vs S150).
    drop-LOG2        op 0x8 disabled.
    drop-LI          op 0x9 disabled.
    drop-DIV_LOG     op 0xa disabled.
    drop-GEO_SUM     op 0xb disabled.
    only-LI          {LOG2, DIV_LOG, GEO_SUM} all disabled (LI alone retained).

Composes:
    E5.8  Brandt-class structural barrier (inherited from C3 / S105).
    E1.3  per-bit difficulty cut (the per-bit family is the target).
    C3.a  S150's extended VM and 7s-per-L_max=24 C inner loop.

Question: which of {LOG2, LI, DIV_LOG, GEO_SUM} is causally responsible
for the S150 cut shift from J*(N) toward ceil(N/2)?  S150 found
compressing programs whose disassembly involved LI, DIV_LOG, or GEO_SUM
in different cells; this experiment isolates which primitive is
necessary in each cell and answers whether LI alone suffices for the
easy-zone shift.

Pre-stated falsifiers (definition.md):
    F1  LI is solely necessary for the easy-zone shift, and any single
        non-LI ablation preserves N=4 J=2 + N=5 J=3 compression.
    F2  Multiple primitives are jointly necessary -- some single drop
        among {DIV_LOG, GEO_SUM} breaks an easy-zone cell.
    F3  LI alone does NOT suffice -- only-LI ablation leaves at least
        one easy-zone cell saturated that the full VM compressed.
    F4  No primitive is necessary -- every single-drop preserves every
        S150 compressed cell, with alternative programs.
"""

import argparse
import ctypes as ct
import math
import os
import sys
import time
from typing import Dict, List, Tuple

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
BRANDT_MKTP_DIR = os.path.normpath(os.path.join(THIS_DIR, "..", "brandt_mktp"))
sys.path.insert(0, BRANDT_MKTP_DIR)
from brandt_mktp import sieve_primes  # noqa: E402

_LIB_PATH = os.path.join(THIS_DIR, "sim_ablation.so")
if not os.path.exists(_LIB_PATH):
    raise SystemExit(
        f"C simulator missing: {_LIB_PATH}.  Build with:\n"
        f"  gcc -O2 -shared -fPIC -o sim_ablation.so sim_ablation.c -lm"
    )
_lib = ct.CDLL(_LIB_PATH)
_lib.simulate.restype = ct.c_int
_lib.simulate.argtypes = [
    ct.c_uint64, ct.c_int, ct.c_int,
    ct.POINTER(ct.c_uint64), ct.POINTER(ct.c_int),
]
_lib.batch_enumerate_at_length_filtered.restype = ct.c_int64
_lib.batch_enumerate_at_length_filtered.argtypes = [
    ct.c_int,
    ct.POINTER(ct.c_uint64), ct.POINTER(ct.c_int), ct.c_int,
    ct.POINTER(ct.c_int), ct.POINTER(ct.c_uint64), ct.POINTER(ct.c_int),
    ct.c_int,
    ct.c_uint,
]
_lib.emit_ops_in_cycle.restype = ct.c_int
_lib.emit_ops_in_cycle.argtypes = [ct.c_uint64, ct.c_int]


T_MAX = 4096
INT_CAP = 10**9

OP_NAMES = [
    "PUSH0", "PUSH1", "PUSH_N", "INC", "ADD", "SUB", "MUL", "SHR1",
    "LOG2", "LI", "DIV_LOG", "GEO_SUM", "DUP", "EMIT_S", "EMIT0", "HALT",
]
OP_INDEX = {n: i for i, n in enumerate(OP_NAMES)}

ABLATIONS: List[Tuple[str, List[str]]] = [
    ("baseline",      []),
    ("drop_LOG2",     ["LOG2"]),
    ("drop_LI",       ["LI"]),
    ("drop_DIV_LOG",  ["DIV_LOG"]),
    ("drop_GEO_SUM",  ["GEO_SUM"]),
    ("only_LI",       ["LOG2", "DIV_LOG", "GEO_SUM"]),
]


def disassemble(prog_bits: int, prog_len_nyb: int) -> str:
    parts = []
    for i in range(prog_len_nyb):
        op = (prog_bits >> (4 * i)) & 0xF
        parts.append(OP_NAMES[op])
    return ", ".join(parts)


def pi_table(N_bits: int) -> List[int]:
    M = (1 << N_bits) - 1
    sieve = sieve_primes(M)
    out = [0] * (M + 1)
    count = 0
    for x in range(M + 1):
        if sieve[x]:
            count += 1
        out[x] = count
    return out


def per_bit_string(pi_vals: List[int], J: int) -> str:
    return "".join("1" if ((v >> J) & 1) else "0" for v in pi_vals)


def bitstr_to_int(s: str) -> int:
    v = 0
    for i, c in enumerate(s):
        if c == '1':
            v |= (1 << i)
    return v


def J_star(pivals: List[int]) -> int:
    max_pi = pivals[-1]
    if max_pi == 0:
        return 0
    return int(math.ceil(math.log2(max_pi + 1)))


def batch_run_filtered(target_strings: List[str], L_max_nyb: int,
                       forbidden_ops: List[str], verbose: bool = False):
    """Run the C inner loop with the given forbidden-op list."""
    forbidden_mask = 0
    for op in forbidden_ops:
        forbidden_mask |= (1 << OP_INDEX[op])

    n_targets = len(target_strings)
    target_lens = [len(t) for t in target_strings]
    target_ints = [bitstr_to_int(t) for t in target_strings]
    max_target_len = max(target_lens)

    INF = 4 * L_max_nyb + int(math.ceil(math.log2(T_MAX + 1)))

    targets_buf = (ct.c_uint64 * n_targets)(*target_ints)
    target_lens_buf = (ct.c_int * n_targets)(*target_lens)
    bests_buf = (ct.c_int * n_targets)(*([INF] * n_targets))
    best_progs_buf = (ct.c_uint64 * n_targets)(*([0] * n_targets))
    best_steps_buf = (ct.c_int * n_targets)(*([-1] * n_targets))

    sim_total = 0
    t_total = 0.0
    for L_nyb in range(1, L_max_nyb + 1):
        L_bits = 4 * L_nyb
        if all(bests_buf[i] <= L_bits + 1 for i in range(n_targets)):
            break
        t0 = time.time()
        n_sim = _lib.batch_enumerate_at_length_filtered(
            L_nyb, targets_buf, target_lens_buf, n_targets,
            bests_buf, best_progs_buf, best_steps_buf,
            max_target_len, ct.c_uint(forbidden_mask),
        )
        dt = time.time() - t0
        t_total += dt
        n_total = 1 << (4 * L_nyb)
        if verbose:
            print(f"    L = {L_bits:2d} bits  ({n_total:11,} total, "
                  f"{n_sim:11,} simulated   t = {dt:6.2f} s)")
        sim_total += n_sim

    bests = list(bests_buf)
    best_progs = list(best_progs_buf)
    best_lens = []
    for i in range(n_targets):
        if best_progs_buf[i] == 0 and bests_buf[i] == INF:
            best_lens.append(-1)
        else:
            steps = best_steps_buf[i]
            log_s = 1 if steps <= 1 else max(1, (steps - 1).bit_length())
            best_lens.append(bests[i] - log_s)
    best_steps = list(best_steps_buf)
    return {
        "bests": bests,
        "best_progs": best_progs,
        "best_lens": best_lens,
        "best_steps": best_steps,
        "INF": INF,
        "sim_total": sim_total,
        "t_total": t_total,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--L_max", type=int, default=28,
                    help="VM program length cap in bits (multiple of 4)")
    ap.add_argument("--N_values", type=int, nargs="+", default=[3, 4, 5, 6])
    ap.add_argument("--ablations", type=str, nargs="+",
                    default=[name for name, _ in ABLATIONS],
                    help="subset of ablation names to run")
    args = ap.parse_args()

    L_max_bits = args.L_max
    assert L_max_bits % 4 == 0
    L_max_nyb = L_max_bits // 4
    INF = 4 * L_max_nyb + int(math.ceil(math.log2(T_MAX + 1)))

    print("=" * 78)
    print("C3.a.iv -- Arithmetic-primitive ABLATION of bounded-Kt VM")
    print("=" * 78)
    print()
    print(f"L_MAX_EXT = {L_max_bits} bits ({L_max_nyb} nibbles), "
          f"T_MAX = {T_MAX}, INT_CAP = {INT_CAP}, INF = {INF}")
    print()

    pi_tables = {N: pi_table(N) for N in args.N_values}

    cell_keys: List[Tuple[int, int]] = []
    cell_strings: List[str] = []
    cell_meta: List[dict] = []
    for N in args.N_values:
        pivals = pi_tables[N]
        n = 1 << N
        max_pi = pivals[-1]
        Js = J_star(pivals)
        boundary_E13 = math.ceil(N / 2)
        for J in range(N):
            s = per_bit_string(pivals, J)
            ones = s.count("1")
            cell_keys.append((N, J))
            cell_strings.append(s)
            cell_meta.append({
                "N": N, "J": J, "len": n, "ones": ones, "zeros": n - ones,
                "id_zero": (max_pi >> J) == 0,
                "zone_E13": "easy" if J >= boundary_E13 else "hard",
                "J_star": Js,
                "above_J_star": J >= Js,
            })

    print(f"Targets: {len(cell_strings)} per-bit cells over N = {args.N_values}")
    print()

    selected = [a for a in ABLATIONS if a[0] in args.ablations]
    runs: Dict[str, dict] = {}
    for name, forbidden in selected:
        print(f"--- Ablation: {name}  forbidden = {forbidden}")
        result = batch_run_filtered(cell_strings, L_max_nyb,
                                    forbidden, verbose=True)
        print(f"    {result['sim_total']:,} sims in {result['t_total']:.1f}s")
        runs[name] = result
        print()

    n_cells = len(cell_strings)

    headers = ["N", "J", "zone", "J*", "id0"]
    for name, _ in selected:
        headers.append(name)
    col_w = [3, 3, 6, 4, 4] + [12] * len(selected)
    fmt = "  " + "  ".join(f"{{:>{w}s}}" for w in col_w)
    print(fmt.format(*headers))
    print("  " + "-" * (sum(col_w) + 2 * (len(col_w) - 1)))

    rows = []
    for i, key in enumerate(cell_keys):
        meta = cell_meta[i]
        cells_str = [
            f"{meta['N']:d}",
            f"{meta['J']:d}",
            meta["zone_E13"],
            f"{meta['J_star']:d}",
            "Y" if meta["id_zero"] else " ",
        ]
        per_ablation = {}
        for name, _ in selected:
            r = runs[name]
            kt = r["bests"][i]
            sat = (kt == r["INF"])
            cells_str.append(f"{kt:3d}{'*' if sat else ' '}     ")
            per_ablation[name] = {
                "kt": kt,
                "sat": sat,
                "prog_bits": r["best_progs"][i],
                "prog_len": r["best_lens"][i],
                "steps": r["best_steps"][i],
                "disasm": (
                    disassemble(r["best_progs"][i], max(0, r["best_lens"][i]) // 4)
                    if not sat else None
                ),
            }
        print(fmt.format(*cells_str))
        rows.append({"meta": meta, "ablations": per_ablation})

    print()
    print("PROGRAMS FOUND in each ablation (compressed cells):")
    for name, _ in selected:
        print(f"\n  >>> {name}")
        for r in rows:
            ab = r["ablations"][name]
            if not ab["sat"]:
                m = r["meta"]
                lbl = f"N={m['N']}, J={m['J']} ({m['zone_E13']})"
                trv = " /trivial" if m["above_J_star"] else ""
                print(f"    {lbl + trv:24s}  Kt = {ab['kt']:3d}  "
                      f"L = {ab['prog_len']:2d} bits, steps = {ab['steps']:5d}"
                      f"   prog: {ab['disasm']}")

    print()
    print("ABLATION DELTA TABLE (Kt change vs baseline; '+' means worse):")
    base_run = runs.get("baseline")
    if base_run is None:
        print("  (no baseline run; skipping)")
    else:
        ab_names = [name for name, _ in selected if name != "baseline"]
        hdr = "  " + "N J  zone     baseline  " + "  ".join(
            f"{n:>14s}" for n in ab_names
        )
        print(hdr)
        print("  " + "-" * (len(hdr) - 2))
        for i, r in enumerate(rows):
            m = r["meta"]
            base = base_run["bests"][i]
            base_str = f"{base:3d}{'*' if base == base_run['INF'] else ' '}"
            cells = []
            for name in ab_names:
                kt = runs[name]["bests"][i]
                if kt == runs[name]["INF"]:
                    s = "INF*"
                elif kt == base and base != base_run["INF"]:
                    s = "  =  "
                else:
                    delta = kt - base if base != base_run["INF"] else None
                    if delta is None:
                        s = f"{kt:3d}"
                    elif delta > 0:
                        s = f"{kt:3d}(+{delta})"
                    elif delta < 0:
                        s = f"{kt:3d}({delta})"
                    else:
                        s = f"{kt:3d}  "
                cells.append(f"{s:>14s}")
            print(f"  {m['N']:1d} {m['J']:1d}  {m['zone_E13']:8s} {base_str}     "
                  + "  ".join(cells))

    print()
    print("VERDICT ANALYSIS:")
    if base_run is None:
        print("  baseline run missing")
    else:
        easy_cells_compressed_baseline = []
        for i, r in enumerate(rows):
            m = r["meta"]
            if (m["zone_E13"] == "easy" and not m["above_J_star"]
                    and base_run["bests"][i] != base_run["INF"]):
                easy_cells_compressed_baseline.append(i)

        print(f"  Easy-zone (J in [ceil(N/2), J*)) cells COMPRESSED in baseline: "
              f"{[(rows[i]['meta']['N'], rows[i]['meta']['J']) for i in easy_cells_compressed_baseline]}")
        print()

        for name, _ in selected:
            if name == "baseline":
                continue
            r = runs[name]
            broken = [
                (rows[i]['meta']['N'], rows[i]['meta']['J'])
                for i in easy_cells_compressed_baseline
                if r["bests"][i] == r["INF"]
            ]
            print(f"  {name:16s}: easy-zone cells now SATURATED = {broken}")

        only_li_breaks = []
        if "only_LI" in [n for n, _ in selected]:
            r = runs["only_LI"]
            for i in easy_cells_compressed_baseline:
                if r["bests"][i] == r["INF"]:
                    only_li_breaks.append((rows[i]['meta']['N'], rows[i]['meta']['J']))

        drop_li_breaks = []
        if "drop_LI" in [n for n, _ in selected]:
            r = runs["drop_LI"]
            for i in easy_cells_compressed_baseline:
                if r["bests"][i] == r["INF"]:
                    drop_li_breaks.append((rows[i]['meta']['N'], rows[i]['meta']['J']))

        non_li_drops_break = {}
        for name in ("drop_LOG2", "drop_DIV_LOG", "drop_GEO_SUM"):
            if name not in [n for n, _ in selected]:
                continue
            r = runs[name]
            non_li_drops_break[name] = [
                (rows[i]['meta']['N'], rows[i]['meta']['J'])
                for i in easy_cells_compressed_baseline
                if r["bests"][i] == r["INF"]
            ]

        print()
        if drop_li_breaks and not only_li_breaks and all(
            len(v) == 0 for v in non_li_drops_break.values()
        ):
            print("  --> F1: LI is solely NECESSARY and SUFFICIENT for "
                  "the easy-zone shift.")
        elif only_li_breaks:
            print("  --> F3: LI alone does NOT suffice; "
                  f"only-LI saturates {only_li_breaks}.")
        elif drop_li_breaks and any(non_li_drops_break.values()):
            print(f"  --> F2: multiple primitives are jointly necessary.  "
                  f"drop_LI breaks {drop_li_breaks}; "
                  f"non-LI drops also break: {non_li_drops_break}")
        elif not drop_li_breaks and all(
            len(v) == 0 for v in non_li_drops_break.values()
        ):
            print("  --> F4: no single primitive is strictly necessary "
                  "(all easy-zone cells survive every single-drop).")
        else:
            print(f"  --> mixed.  drop_LI breaks {drop_li_breaks}; "
                  f"only_LI breaks {only_li_breaks}; "
                  f"non-LI drops break {non_li_drops_break}")

        print()
        print("  Per-cell hard-zone / J* trivial impact (informational):")
        for i, r in enumerate(rows):
            m = r["meta"]
            base = base_run["bests"][i]
            if base == base_run["INF"]:
                continue
            damage = []
            for name, _ in selected:
                if name == "baseline":
                    continue
                kt = runs[name]["bests"][i]
                if kt == runs[name]["INF"]:
                    damage.append(f"{name}=SAT")
                elif kt > base + 4:
                    damage.append(f"{name}=+{kt - base}")
            if damage:
                print(f"    N={m['N']}, J={m['J']} ({m['zone_E13']}, "
                      f"baseline {base}): " + ", ".join(damage))

    print()
    print("=" * 78)


if __name__ == "__main__":
    main()
