"""
F6.a — Cross-modulus generalisation of S246's dyadic structural test.

For m in {3, 5, 6, 10, 30}, compute residuals
    r_R(k) = pi(m^k) - R(m^k),  c_R(k) = r_R(k) / m^{k/2}
and run the same three structural tests as S246 (dyadic_pi_structure):
  F1: Berlekamp-Massey GF(2) linear complexity of sign(r_R(k)) vs MC random-shuffle null.
  F2: Max |autocorr(c_R, lag)| over lags 1..min(10, K_m//3) vs MC iid-Gaussian null.
  F3: omitted (HKM cost test ran in S246 already and depends on x-magnitude not m).

Hypothesis (B-NEGATIVE shape carries over): for every m with log(m) / 2*pi
irrational (which holds for all m in {3, 5, 6, 10, 30} by
Lindemann-Weierstrass), Weyl equidistribution of {gamma_n * log(m) * k mod 2*pi}
forces F1+F2 to land at random-baseline values.
A-grade falsifier: any single m exhibits BM <= MC p05 OR max |ac| >= MC p999.

Pre-cached pi(m^k) values (verified vs sympy.primepi where feasible).
m^k computed live with sympy.primepi for k beyond the cache; cache keys
verified against sympy for the smallest k in each m.

CLI:
    python3 cross_modulus_pi_structure.py [--m 3,5,6,10,30] [--max-runtime 600]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import mpmath
from sympy import primepi

mpmath.mp.dps = 80


# OEIS A006880: pi(10^n) for n = 0..28. Used as ground truth (m=10).
PI_10K = {
    0: 0,
    1: 4,
    2: 25,
    3: 168,
    4: 1229,
    5: 9592,
    6: 78498,
    7: 664579,
    8: 5761455,
    9: 50847534,
    10: 455052511,
    11: 4118054813,
    12: 37607912018,
    13: 346065536839,
    14: 3204941750802,
    15: 29844570422669,
    16: 279238341033925,
    17: 2623557157654233,
    18: 24739954287740860,
    19: 234057667276344607,
    20: 2220819602560918840,
    21: 21127269486018731928,
    22: 201467286689315906290,
    23: 1925320391606803968923,
    24: 18435599767349200867866,
    25: 176846309399143769411680,
    26: 1699246750872437141327603,
    27: 16352460426841680446427399,
    28: 157589269275973410412739598,
}


# k-budgets per m. Largest m^k stays in sympy.primepi tractable range
# (~10^11 max for the heavy m's; m=10 served by OEIS up to k=28).
K_BUDGET = {
    3: 22,    # 3^22 = 3.1e10
    5: 15,    # 5^15 = 3.1e10
    6: 14,    # 6^14 = 7.8e10
    10: 28,   # served by OEIS A006880; sympy verifies k=1..10
    30: 8,    # 30^8 = 6.6e11
}


def compute_pi_m_k(m: int, max_k: int, oeis_table: dict | None = None,
                   verify_up_to: int = 9) -> dict:
    """Return {k: pi(m^k)} for k in 1..max_k.

    For m=10, prefer OEIS A006880 cache; verify k <= verify_up_to via sympy.
    For other m, compute via sympy.primepi for all k (no OEIS cache).
    """
    out = {}
    print(f"  computing pi({m}^k) for k=1..{max_k}", flush=True)
    for k in range(1, max_k + 1):
        x = m ** k
        if oeis_table is not None and k in oeis_table:
            cache_val = oeis_table[k]
            if k <= verify_up_to:
                t = time.perf_counter()
                sympy_val = int(primepi(x))
                dt = time.perf_counter() - t
                if sympy_val != cache_val:
                    raise RuntimeError(
                        f"OEIS-vs-sympy mismatch m={m} k={k}: "
                        f"oeis={cache_val} sympy={sympy_val}")
                out[k] = {"x": x, "pi_x": cache_val, "source": "oeis_verified",
                          "verify_dt_s": dt}
            else:
                out[k] = {"x": x, "pi_x": cache_val, "source": "oeis"}
        else:
            t = time.perf_counter()
            v = int(primepi(x))
            dt = time.perf_counter() - t
            out[k] = {"x": x, "pi_x": v, "source": "sympy", "compute_dt_s": dt}
            print(f"    k={k:2d}: pi({m}^{k}) = pi({x}) = {v}  [{dt:.2f}s]",
                  flush=True)
    return out


def explicit_residuals(m: int, pi_table: dict) -> dict:
    """Compute R(m^k), residual r_R(k) = pi(m^k) - R(m^k), and
    normalised c_R(k) = r_R(k) / m^{k/2}.
    """
    rows = {}
    for k, info in pi_table.items():
        x = mpmath.mpf(m) ** k
        Rx = mpmath.riemannr(x)
        pi_val = mpmath.mpf(info["pi_x"])
        r_R = pi_val - Rx
        c_R = r_R / mpmath.sqrt(x)
        rows[k] = {
            "pi": str(info["pi_x"]),
            "R": mpmath.nstr(Rx, 30),
            "r_R": mpmath.nstr(r_R, 20),
            "c_R": float(c_R),
            "sign_r_R": int(mpmath.sign(r_R)),
        }
    return rows


def autocorr(seq, lag):
    n = len(seq)
    if n - lag < 2:
        return float("nan")
    mn = sum(seq) / n
    num = sum((seq[i] - mn) * (seq[i + lag] - mn) for i in range(n - lag))
    den = sum((seq[i] - mn) ** 2 for i in range(n))
    if den == 0:
        return 0.0
    return num / den


def berlekamp_massey_gf2(s):
    """Linear complexity of binary sequence s over GF(2)."""
    n = len(s)
    b = [0] * n
    c = [0] * n
    b[0] = c[0] = 1
    L = 0
    m = 1
    for i in range(n):
        d = s[i]
        for j in range(1, L + 1):
            d ^= c[j] & s[i - j]
        if d == 1:
            t = c[:]
            for j in range(n - m):
                if j + m < n:
                    c[j + m] ^= b[j]
            if 2 * L <= i:
                L = i + 1 - L
                b = t
                m = 1
            else:
                m += 1
        else:
            m += 1
    return L


def monte_carlo_baselines(n, n_pos, n_trials=4000, seed=42):
    """MC: BM null on random binary with given pos/neg split; iid-Gaussian
    autocorr null at length n.
    """
    import random
    rng = random.Random(seed)
    bm_vals = []
    for _ in range(n_trials):
        s = [1] * n_pos + [0] * (n - n_pos)
        rng.shuffle(s)
        bm_vals.append(berlekamp_massey_gf2(s))
    bm_vals.sort()
    bm_mean = sum(bm_vals) / n_trials
    bm_std = (sum((v - bm_mean) ** 2 for v in bm_vals) / n_trials) ** 0.5

    ac_vals = []
    for _ in range(n_trials):
        seq = [rng.gauss(0, 1) for _ in range(n)]
        ac_vals.append(autocorr(seq, 1))
    ac_vals_abs = sorted(abs(v) for v in ac_vals)
    return {
        "n_trials": n_trials,
        "bm_mean": bm_mean,
        "bm_std": bm_std,
        "bm_p05": bm_vals[int(0.05 * n_trials)],
        "bm_p95": bm_vals[int(0.95 * n_trials)],
        "abs_ac_lag1_p95": ac_vals_abs[int(0.95 * n_trials)],
        "abs_ac_lag1_p99": ac_vals_abs[int(0.99 * n_trials)],
        "abs_ac_lag1_p999": ac_vals_abs[int(0.999 * n_trials)],
    }


def runs_test(s):
    if not s:
        return {"runs": 0, "n": 0}
    n = len(s)
    n1 = sum(s)
    n0 = n - n1
    runs = 1 + sum(1 for i in range(1, n) if s[i] != s[i - 1])
    if n1 == 0 or n0 == 0:
        return {"runs": runs, "n": n, "n1": n1, "n0": n0,
                "expected": float("nan"), "z": float("nan"),
                "note": "degenerate"}
    expected = 2.0 * n1 * n0 / n + 1.0
    var = (expected - 1.0) * (expected - 2.0) / (n - 1)
    if var <= 0:
        return {"runs": runs, "n": n, "n1": n1, "n0": n0,
                "expected": expected, "z": float("nan"), "note": "var<=0"}
    z = (runs - expected) / (var ** 0.5)
    return {"runs": runs, "n": n, "n1": n1, "n0": n0,
            "expected": expected, "z": z}


def run_one_m(m, max_k):
    print(f"\n=== m = {m} ===", flush=True)
    oeis = PI_10K if m == 10 else None
    pi_table = compute_pi_m_k(m, max_k, oeis_table=oeis)
    residuals = explicit_residuals(m, pi_table)

    ks = sorted(residuals)
    c_R_seq = [residuals[k]["c_R"] for k in ks]
    sign_R_seq = [residuals[k]["sign_r_R"] for k in ks]
    bit_R_seq = [0 if s < 0 else 1 for s in sign_R_seq]
    n = len(bit_R_seq)
    n_pos = sum(bit_R_seq)
    n_neg = n - n_pos

    # Lag set: 1..min(10, n//3) so each lag uses >= 2n/3 pairs.
    max_lag = min(10, max(1, n // 3))
    ac_R = {lag: autocorr(c_R_seq, lag) for lag in range(1, max_lag + 1)}
    bm = berlekamp_massey_gf2(bit_R_seq)
    rt = runs_test(bit_R_seq)

    print(f"  n={n}, sign split: pos={n_pos} neg={n_neg}", flush=True)
    print(f"  BM linear complexity = {bm}", flush=True)
    print(f"  autocorr_R lag 1..{max_lag}: "
          f"{[round(ac_R[l], 3) for l in range(1, max_lag + 1)]}", flush=True)

    print(f"  MC baselines (length n={n}, pos={n_pos})...", flush=True)
    mc = monte_carlo_baselines(n, n_pos, n_trials=4000)
    print(f"    BM null: mean={mc['bm_mean']:.2f} std={mc['bm_std']:.2f} "
          f"[p05={mc['bm_p05']}, p95={mc['bm_p95']}]", flush=True)
    print(f"    ac null: p95={mc['abs_ac_lag1_p95']:.3f} "
          f"p99={mc['abs_ac_lag1_p99']:.3f} "
          f"p999={mc['abs_ac_lag1_p999']:.3f}", flush=True)

    skewes_trivial = (n_pos == 0 or n_neg == 0)
    f1 = (not skewes_trivial) and (bm <= mc["bm_p05"])
    max_abs_ac = max(abs(ac_R[l]) for l in range(1, max_lag + 1))
    max_abs_ac_lag = max(range(1, max_lag + 1),
                        key=lambda l: abs(ac_R[l]))
    f2 = max_abs_ac >= mc["abs_ac_lag1_p999"]

    print(f"  F1 BM <= p05={mc['bm_p05']}: {'HOLDS' if f1 else 'FAILS'} "
          f"(bm={bm})", flush=True)
    print(f"  F2 max|ac| >= p999={mc['abs_ac_lag1_p999']:.3f}: "
          f"{'HOLDS' if f2 else 'FAILS'} "
          f"(max={max_abs_ac:.3f} at lag {max_abs_ac_lag})", flush=True)

    return {
        "m": m,
        "max_k": max_k,
        "n_anchors": n,
        "pi_table": pi_table,
        "residuals": residuals,
        "ac_R": ac_R,
        "max_abs_ac": max_abs_ac,
        "max_abs_ac_lag": max_abs_ac_lag,
        "bm_linear_complexity_sign_R": bm,
        "n_pos": n_pos,
        "n_neg": n_neg,
        "skewes_trivial": skewes_trivial,
        "runs_test": rt,
        "mc": mc,
        "f1_holds": bool(f1),
        "f2_holds": bool(f2),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", default="3,5,6,10,30",
                    help="comma list of moduli")
    ap.add_argument("--out", default=str(Path(__file__).parent),
                    help="output directory")
    args = ap.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    ms = [int(x) for x in args.m.split(",") if x.strip()]
    results = {}
    t0 = time.perf_counter()
    for m in ms:
        max_k = K_BUDGET[m]
        results[m] = run_one_m(m, max_k)
    elapsed = time.perf_counter() - t0
    print(f"\n  total wall time: {elapsed:.1f}s", flush=True)

    # Aggregate verdict
    any_f1 = any(r["f1_holds"] for r in results.values())
    any_f2 = any(r["f2_holds"] for r in results.values())
    if any_f1 or any_f2:
        verdict = "PARTIAL_POSITIVE_SURPRISE"
    else:
        verdict = "B_NEGATIVE_universal_cross_modulus"

    summary = {
        "mode": "F6.a cross-modulus structural test",
        "moduli_tested": ms,
        "verdict": verdict,
        "any_f1_holds": any_f1,
        "any_f2_holds": any_f2,
        "elapsed_s": elapsed,
        "per_m": {m: {
            "m": r["m"],
            "n_anchors": r["n_anchors"],
            "bm_value": r["bm_linear_complexity_sign_R"],
            "bm_p05": r["mc"]["bm_p05"],
            "bm_mean": r["mc"]["bm_mean"],
            "max_abs_ac": r["max_abs_ac"],
            "max_abs_ac_lag": r["max_abs_ac_lag"],
            "ac_p999": r["mc"]["abs_ac_lag1_p999"],
            "n_pos": r["n_pos"], "n_neg": r["n_neg"],
            "skewes_trivial": r["skewes_trivial"],
            "f1_holds": r["f1_holds"],
            "f2_holds": r["f2_holds"],
        } for m, r in results.items()},
    }

    raw_path = out_dir / "raw_data.json"
    with open(raw_path, "w") as f:
        # convert int-keyed dicts to str for JSON
        def convert(o):
            if isinstance(o, dict):
                return {str(k): convert(v) for k, v in o.items()}
            if isinstance(o, list):
                return [convert(v) for v in o]
            return o
        json.dump(convert(results), f, indent=2, default=str)
    print(f"\n  wrote {raw_path}", flush=True)

    sum_path = out_dir / "summary.json"
    with open(sum_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"  wrote {sum_path}", flush=True)
    print(f"\nVERDICT: {verdict}", flush=True)
    print(f"  any F1 holds: {any_f1}", flush=True)
    print(f"  any F2 holds: {any_f2}", flush=True)


if __name__ == "__main__":
    main()
