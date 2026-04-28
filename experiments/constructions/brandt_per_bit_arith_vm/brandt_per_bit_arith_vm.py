#!/usr/bin/env python3
"""
C3.a -- Arithmetic-primitive bounded-Kt VM (S105 successor).

Composes E5.8 (Brandt-class barrier, structural) with E1.3 (per-bit
difficulty cut), via an EXTENDED bounded-Kt simulator that adds the
four arithmetic primitives proposed in S105's successor block:

    LOG2     -- TOS <- floor(log2 TOS)
    LI_APPROX -- TOS <- floor(TOS / log TOS)        (R^{-1} kernel)
    DIV_LOG  -- pop a, b; push floor(b / log max(a, 2))
    GEO_SUM  -- pop a; push 1 + a + a^2 + ... while sum <= T_MAX

Question: does the bounded-Kt cut on the per-bit family

    s_J^(N) := ( bit_J(pi(0)), ..., bit_J(pi(2^N - 1)) ) in {0,1}^{2^N}

shift from S105's J*(N) = ceil(log2(pi(2^N) + 1)) toward E1.3's
0.5*N once the VM can compute Li-style approximations?

The simulator's hot loop is implemented in C (sim.c -> sim.so) and
called via ctypes for ~50x speedup over pure Python.  Pre-check skips
programs whose pc cycle contains no EMIT op (would emit nothing in
T_MAX steps).

Pre-stated falsifiers (definition.md):
  F1 -- full shift to ceil(N/2)
  F2 -- no shift; easy zone saturates
  F3 -- intermediate hierarchy (some easy J compress, others saturate)
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
from brandt_mktp import sieve_primes, Kt_bounded as Kt_bounded_orig  # noqa: E402

# Load C simulator.
_LIB_PATH = os.path.join(THIS_DIR, "sim.so")
if not os.path.exists(_LIB_PATH):
    raise SystemExit(
        f"C simulator missing: {_LIB_PATH}.  Build with:\n"
        f"  gcc -O2 -shared -fPIC -o sim.so sim.c -lm"
    )
_lib = ct.CDLL(_LIB_PATH)
_lib.simulate.restype = ct.c_int
_lib.simulate.argtypes = [
    ct.c_uint64, ct.c_int, ct.c_int,
    ct.POINTER(ct.c_uint64), ct.POINTER(ct.c_int),
]
_lib.batch_enumerate_at_length.restype = ct.c_int64
_lib.batch_enumerate_at_length.argtypes = [
    ct.c_int,
    ct.POINTER(ct.c_uint64), ct.POINTER(ct.c_int), ct.c_int,
    ct.POINTER(ct.c_int), ct.POINTER(ct.c_uint64), ct.POINTER(ct.c_int),
    ct.c_int,
]
_lib.emit_ops_in_cycle.restype = ct.c_int
_lib.emit_ops_in_cycle.argtypes = [ct.c_uint64, ct.c_int]


# Op set
T_MAX = 4096
INT_CAP = 10**9

OP_NAMES = [
    "PUSH0", "PUSH1", "PUSH_N", "INC", "ADD", "SUB", "MUL", "SHR1",
    "LOG2", "LI", "DIV_LOG", "GEO_SUM", "DUP", "EMIT_S", "EMIT0", "HALT",
]


def disassemble(prog_bits: int, prog_len_nyb: int) -> str:
    parts = []
    for i in range(prog_len_nyb):
        op = (prog_bits >> (4 * i)) & 0xF
        parts.append(OP_NAMES[op])
    return ", ".join(parts)


def simulate_py(prog_bits: int, prog_len_nyb: int, target_len: int) -> Tuple[int, int, int]:
    """Python wrapper around the C simulator (for debugging / single-shot use).

    Returns (out_bits, out_count, steps).
    """
    out = ct.c_uint64(0)
    st = ct.c_int(0)
    oc = _lib.simulate(prog_bits, prog_len_nyb, target_len, ct.byref(out), ct.byref(st))
    return out.value, oc, st.value


# ---------------------------------------------------------------------
# Per-bit pi family
# ---------------------------------------------------------------------

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
    """Pack 0/1 string LSB-first to int (bit 0 = first char)."""
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


# ---------------------------------------------------------------------
# Batch Kt_b' evaluator (C-backed)
# ---------------------------------------------------------------------

def batch_kt_ext_c(target_strings: List[str], L_max_nyb: int):
    """Returns parallel lists: bests, best_prog_bits, best_prog_lens (in bits),
    best_steps."""
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
    for L_nyb in range(1, L_max_nyb + 1):
        L_bits = 4 * L_nyb
        # If every best is already <= L_bits + 1, no longer programs can beat them.
        if all(bests_buf[i] <= L_bits + 1 for i in range(n_targets)):
            break
        t0 = time.time()
        n_sim = _lib.batch_enumerate_at_length(
            L_nyb, targets_buf, target_lens_buf, n_targets,
            bests_buf, best_progs_buf, best_steps_buf,
            max_target_len,
        )
        dt = time.time() - t0
        n_total = 1 << (4 * L_nyb)
        print(
            f"  L = {L_bits:2d} bits  ({n_total:11,} total, {n_sim:11,} simulated"
            f"   t = {dt:6.2f} s)"
        )
        sim_total += n_sim
    print(f"  Total simulations: {sim_total:,}")

    bests = list(bests_buf)
    best_progs = list(best_progs_buf)
    best_lens = []
    for i in range(n_targets):
        if best_progs_buf[i] == 0 and bests_buf[i] == INF:
            best_lens.append(-1)
        else:
            # Reconstruct length from bests = L + log2_steps
            # We don't store the L directly; figure from bests and best_steps.
            steps = best_steps_buf[i]
            if steps <= 1:
                log_s = 1
            else:
                log_s = (steps - 1).bit_length()
                if log_s == 0:
                    log_s = 1
            best_lens.append(bests[i] - log_s)
    best_steps = list(best_steps_buf)
    return bests, best_progs, best_lens, best_steps, INF


# ---------------------------------------------------------------------
# Main scan
# ---------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--L_max", type=int, default=24,
                    help="extended-VM program length cap in bits (multiple of 4)")
    ap.add_argument("--N_values", type=int, nargs="+", default=[3, 4, 5, 6])
    ap.add_argument("--skip-orig", action="store_true",
                    help="skip slow original-VM Kt_b reference column")
    args = ap.parse_args()

    L_max_bits = args.L_max
    assert L_max_bits % 4 == 0, "L_max must be a multiple of 4"
    L_max_nyb = L_max_bits // 4
    INF_EXT_local = 4 * L_max_nyb + int(math.ceil(math.log2(T_MAX + 1)))

    print("=" * 78)
    print("C3.a -- Arithmetic-primitive bounded-Kt VM (S105 successor)")
    print("=" * 78)
    print()
    print(f"Extended VM: 16 ops, 4-bit encoding, L_MAX_EXT = {L_max_bits}")
    print(f"             T_MAX = {T_MAX}, INT_CAP = {INT_CAP}")
    print(f"             INF_EXT = {INF_EXT_local}")
    print()
    print("Op set: " + ", ".join(f"{i:x}={n}" for i, n in enumerate(OP_NAMES)))
    print()

    # Build pi tables once.
    N_values = args.N_values
    pi_tables = {N: pi_table(N) for N in N_values}
    Js_per_N = {N: J_star(pi_tables[N]) for N in N_values}

    cal_targets = {
        "0^16":     "0" * 16,
        "1^16":     "1" * 16,
        "(01)^8":   "01" * 8,
        "step8":    "0" * 8 + "1" * 8,
        "step7_9":  "0" * 7 + "1" * 9,
    }
    cal_names = list(cal_targets.keys())
    cal_strings = [cal_targets[k] for k in cal_names]

    cell_keys: List[Tuple[int, int]] = []
    cell_strings: List[str] = []
    cell_meta: List[dict] = []
    for N in N_values:
        pivals = pi_tables[N]
        n = 1 << N
        max_pi = pivals[-1]
        Js = Js_per_N[N]
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

    all_targets = cal_strings + cell_strings
    n_cal = len(cal_strings)
    n_cells = len(cell_strings)
    print(f"Targets: {n_cal} calibration + {n_cells} per-bit cells = {len(all_targets)}")
    print()

    print(f"Enumerating extended-VM programs up to L = {L_max_bits} bits "
          f"(C inner loop, batch evaluator):")
    t0 = time.time()
    bests, best_progs, best_lens, best_steps, INF = batch_kt_ext_c(all_targets, L_max_nyb)
    dt = time.time() - t0
    print(f"  Done in {dt:.1f} s.")
    print()

    print("CALIBRATION -- Kt_b' on canonical targets:")
    for i, name in enumerate(cal_names):
        kt = bests[i]
        if kt == INF:
            print(f"  {name:12s}  Kt_b' = {kt:3d}  [saturated]")
        else:
            asm = disassemble(best_progs[i], best_lens[i] // 4)
            print(f"  {name:12s}  Kt_b' = {kt:3d}  L={best_lens[i]}b, steps={best_steps[i]:5d}  prog: {asm}")
    print()

    if not args.skip_orig:
        print("PER-BIT SCAN  (computing original Kt_b column for reference)...")
        kt_orig_list: List[int] = []
        t0 = time.time()
        for s in cell_strings:
            kt_orig_list.append(Kt_bounded_orig(s))
        print(f"  original column done in {time.time() - t0:.1f} s")
    else:
        kt_orig_list = [-1] * n_cells
    print()

    rows = []
    for i, key in enumerate(cell_keys):
        idx = n_cal + i
        meta = cell_meta[i]
        row = dict(meta)
        row["kt_orig"] = kt_orig_list[i]
        row["kt_ext"] = bests[idx]
        row["ext_compressed"] = bests[idx] < INF
        if best_lens[idx] > 0:
            row["prog_bits"] = best_progs[idx]
            row["prog_len_bits"] = best_lens[idx]
            row["prog_steps"] = best_steps[idx]
            row["prog_disasm"] = disassemble(best_progs[idx], best_lens[idx] // 4)
        else:
            row["prog_bits"] = -1
            row["prog_len_bits"] = -1
            row["prog_steps"] = -1
            row["prog_disasm"] = None
        rows.append(row)

    header = (
        "  N   J  | len |  ones zeros | zone_E13 | J* | id-0 | "
        "Kt_b(orig) | Kt_b'(ext)"
    )
    print(header)
    print("  " + "-" * (len(header) - 2))
    cur_N = None
    for r in rows:
        if cur_N is not None and r["N"] != cur_N:
            print()
        cur_N = r["N"]
        sat_orig = "*" if r["kt_orig"] == 61 else " "
        sat_ext = "*" if r["kt_ext"] == INF else " "
        iz = "Y" if r["id_zero"] else " "
        kt_orig_str = "  -- " if r["kt_orig"] < 0 else f"  {r['kt_orig']:3d} {sat_orig}"
        print(
            f"  {r['N']:2d}  {r['J']:2d}  | {r['len']:3d} |  {r['ones']:4d}  {r['zeros']:4d} |"
            f" {r['zone_E13']:8s} | {r['J_star']:2d} |  {iz}   |"
            f"{kt_orig_str}    |"
            f"  {r['kt_ext']:3d} {sat_ext}"
        )
    print()

    print("PROGRAMS FOUND  (compressed cells in extended VM):")
    for r in rows:
        if r["ext_compressed"]:
            zone = r["zone_E13"]
            label = f"N={r['N']:1d}, J={r['J']:1d} ({zone}"
            if r["above_J_star"]:
                label += "/trivial)"
            else:
                label += ")"
            print(
                f"  {label:24s}  Kt_b' = {r['kt_ext']:3d}  "
                f"L = {r['prog_len_bits']:2d} bits, steps = {r['prog_steps']:5d}"
                f"   prog: {r['prog_disasm']}"
            )
    print()

    by_N: Dict[int, List[dict]] = {}
    for r in rows:
        by_N.setdefault(r["N"], []).append(r)
    f1_per_N = {}
    f2_per_N = {}
    diag = {}
    for N, group in by_N.items():
        if N < 4:
            continue
        boundary = math.ceil(N / 2)
        easy = [r for r in group if boundary <= r["J"] < r["J_star"]]
        hard = [r for r in group if r["J"] < boundary]
        triv = [r for r in group if r["J"] >= r["J_star"]]
        diag[N] = {
            "boundary": boundary,
            "J_star": group[0]["J_star"],
            "easy_J": [r["J"] for r in easy],
            "hard_J": [r["J"] for r in hard],
            "triv_J": [r["J"] for r in triv],
            "easy_compressed": [r["J"] for r in easy if r["ext_compressed"]],
            "easy_saturated": [r["J"] for r in easy if not r["ext_compressed"]],
            "hard_compressed": [r["J"] for r in hard if r["ext_compressed"]],
            "hard_saturated": [r["J"] for r in hard if not r["ext_compressed"]],
            "triv_compressed": [r["J"] for r in triv if r["ext_compressed"]],
        }
        f1_per_N[N] = (
            len(easy) > 0
            and len(diag[N]["easy_saturated"]) == 0
            and len(diag[N]["hard_compressed"]) == 0
        )
        f2_per_N[N] = (
            len(easy) > 0 and len(diag[N]["easy_compressed"]) == 0
        )

    if f1_per_N and all(f1_per_N.values()):
        verdict = "F1 (full shift to ceil(N/2))"
    elif f2_per_N and all(f2_per_N.values()):
        verdict = "F2 (no shift; easy zone saturates as in original VM)"
    else:
        verdict = "F3 (intermediate hierarchy)"

    print(f"VERDICT: {verdict}")
    print()
    print("PER-N DIAGNOSTIC:")
    for N in sorted(diag.keys()):
        d = diag[N]
        print(
            f"  N={N}: ceil(N/2)={d['boundary']}, J*(N)={d['J_star']}"
        )
        print(
            f"        easy zone J = {d['easy_J']};   hard J = {d['hard_J']};   "
            f"trivial J = {d['triv_J']}"
        )
        print(
            f"        easy compressed = {d['easy_compressed']};   "
            f"easy saturated = {d['easy_saturated']}"
        )
        print(
            f"        hard compressed = {d['hard_compressed']};   "
            f"trivial compressed = {d['triv_compressed']}"
        )
    print()
    print("=" * 78)


if __name__ == "__main__":
    main()
