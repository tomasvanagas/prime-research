"""
F6 / P7 — Structural tests on the dyadic family π(2^k), k = 1..56.

Pre-registered F-claims (see results.md).

Outputs:
  raw_data.json       — π(2^k), Li(2^k), r(k), c(k), sign(r(k))
  hkm_cost.json       — wall-clock T(2^k), T(2^k ± 1) for k ∈ {15, 18, 20}
  summary.json        — verdicts on F-claims 1 / 2 / 3

CLI:
    python3 dyadic_pi_structure.py [--cost-anchors 15,18,20] [--cost-trials 5]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Iterable

import mpmath
from sympy import primepi

mpmath.mp.dps = 80  # high precision for Li(2^56) ≈ 2^56 / 39

# OEIS A007053: π(2^k) for k = 0..62.
PI_2K = {
    0: 0,
    1: 1,
    2: 2,
    3: 4,
    4: 6,
    5: 11,
    6: 18,
    7: 31,
    8: 54,
    9: 97,
    10: 172,
    11: 309,
    12: 564,
    13: 1028,
    14: 1900,
    15: 3512,
    16: 6542,
    17: 12251,
    18: 23000,
    19: 43390,
    20: 82025,
    21: 155611,
    22: 295947,
    23: 564163,
    24: 1077871,
    25: 2063689,
    26: 3957809,
    27: 7603553,
    28: 14630843,
    29: 28192750,
    30: 54400028,
    31: 105097565,
    32: 203280221,
    33: 393615806,
    34: 762939111,
    35: 1480206279,
    36: 2874398515,
    37: 5586502348,
    38: 10866266172,
    39: 21151907950,
    40: 41203088796,
    41: 80316571436,
    42: 156661034233,
    43: 305761713237,
    44: 597116381732,
    45: 1166746786182,
    46: 2280998753949,
    47: 4461632979717,
    48: 8731188863470,
    49: 17094432576778,
    50: 33483379603407,
    51: 65612899915304,
    52: 128625503610475,
    53: 252252704148404,
    54: 494890204904784,
    55: 971269945245201,
    56: 1906879381028850,
}


def verify_via_sympy(max_k: int) -> dict:
    """Verify hard-coded values against sympy.primepi for k <= max_k."""
    out = {}
    for k in range(1, max_k + 1):
        x = 1 << k
        v = int(primepi(x))
        ok = (v == PI_2K[k])
        out[k] = {"sympy_pi_2k": v, "table_pi_2k": PI_2K[k], "agree": ok}
        if not ok:
            print(f"  MISMATCH at k={k}: sympy={v}, table={PI_2K[k]}", file=sys.stderr)
    return out


def explicit_residuals() -> dict:
    """Compute Li(2^k), R(2^k), residuals against both, normalised by 2^{k/2}.

    The relevant baseline for "zero-driven oscillation" is R(x), not Li(x):
    Li(x) = R(x) + 0.5 Li(sqrt(x)) + (1/3) Li(x^{1/3}) + ... so
    π(x) - Li(x) carries a *smooth* -0.5 Li(sqrt(x)) trend, while
    π(x) - R(x) is the pure zero-sum part (cosines of γ_n log x with
    ~1/γ_n weights). Test pseudorandomness on c_R, not c_Li.
    """
    rows = {}
    for k in range(1, 57):
        x = mpmath.mpf(2) ** k
        Li = mpmath.li(x)
        Rx = mpmath.riemannr(x)
        pi_val = mpmath.mpf(PI_2K[k])
        r_li = pi_val - Li
        r_R = pi_val - Rx
        c_li = r_li / mpmath.sqrt(x)
        c_R = r_R / mpmath.sqrt(x)
        rows[k] = {
            "pi_2k": str(PI_2K[k]),
            "Li_2k": mpmath.nstr(Li, 30),
            "R_2k": mpmath.nstr(Rx, 30),
            "r_li": mpmath.nstr(r_li, 20),
            "r_R": mpmath.nstr(r_R, 20),
            "c_li": float(c_li),
            "c_R": float(c_R),
            "sign_r_li": int(mpmath.sign(r_li)),
            "sign_r_R": int(mpmath.sign(r_R)),
        }
    return rows


def autocorr(seq: list[float], lag: int) -> float:
    n = len(seq)
    if n - lag < 2:
        return float("nan")
    m = sum(seq) / n
    num = sum((seq[i] - m) * (seq[i + lag] - m) for i in range(n - lag))
    den = sum((seq[i] - m) ** 2 for i in range(n))
    if den == 0:
        return 0.0
    return num / den


def berlekamp_massey_gf2(s: list[int]) -> int:
    """Return linear complexity of binary sequence s over GF(2)."""
    n = len(s)
    b = [0] * n
    c = [0] * n
    b[0] = c[0] = 1
    L = 0
    m = 1
    bb = 1
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
                bb = 1
                m = 1
            else:
                m += 1
        else:
            m += 1
    return L


def monte_carlo_baselines(n: int, n_pos: int, n_trials: int = 4000,
                          seed: int = 42) -> dict:
    """Monte Carlo: BM complexity and lag-1 autocorr distributions for
    (a) uniform random binary sequences of length n with given marginal
        n_pos / n_neg split, and
    (b) standard-Gaussian iid sequences for autocorrelation.
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

    # Empirical CDF for autocorrelation under iid Gaussian
    ac_vals = []
    for _ in range(n_trials):
        seq = [rng.gauss(0, 1) for _ in range(n)]
        ac_vals.append(autocorr(seq, 1))
    ac_vals_abs = sorted(abs(v) for v in ac_vals)
    ac_p95 = ac_vals_abs[int(0.95 * n_trials)]
    ac_p99 = ac_vals_abs[int(0.99 * n_trials)]
    ac_p999 = ac_vals_abs[int(0.999 * n_trials)]

    return {
        "n_trials": n_trials,
        "bm_mean": bm_mean,
        "bm_std": bm_std,
        "bm_p05": bm_vals[int(0.05 * n_trials)],
        "bm_p95": bm_vals[int(0.95 * n_trials)],
        "bm_min": bm_vals[0],
        "bm_max": bm_vals[-1],
        "abs_ac_lag1_p95": ac_p95,
        "abs_ac_lag1_p99": ac_p99,
        "abs_ac_lag1_p999": ac_p999,
    }


def runs_test(s: list[int]) -> dict:
    """Wald-Wolfowitz runs test: count runs of identical bits."""
    if not s:
        return {"runs": 0, "n": 0, "expected": 0.0, "z": 0.0}
    n = len(s)
    n1 = sum(s)
    n0 = n - n1
    runs = 1 + sum(1 for i in range(1, n) if s[i] != s[i - 1])
    if n1 == 0 or n0 == 0:
        return {"runs": runs, "n": n, "n1": n1, "n0": n0,
                "expected": float("nan"), "z": float("nan"),
                "note": "degenerate (constant)"}
    expected = 2.0 * n1 * n0 / n + 1.0
    var = (expected - 1.0) * (expected - 2.0) / (n - 1)
    if var <= 0:
        return {"runs": runs, "n": n, "n1": n1, "n0": n0,
                "expected": expected, "z": float("nan"),
                "note": "var<=0"}
    z = (runs - expected) / (var ** 0.5)
    return {"runs": runs, "n": n, "n1": n1, "n0": n0,
            "expected": expected, "z": z}


def _subprocess_time_primepi(x: int) -> float:
    """Cold-start sympy.primepi(x) in a fresh subprocess; return wall time."""
    import subprocess
    code = (
        "import time\n"
        "from sympy import primepi\n"
        f"x = {x}\n"
        "t0 = time.perf_counter()\n"
        "v = int(primepi(x))\n"
        "t1 = time.perf_counter()\n"
        "print(f'{t1 - t0:.6f} {v}')\n"
    )
    proc = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True, text=True, timeout=300,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"subprocess failed: {proc.stderr}")
    parts = proc.stdout.strip().split()
    return float(parts[0]), int(parts[1])


def hkm_cost_test(anchors: Iterable[int], trials: int) -> dict:
    """Measure sympy.primepi cost at 2^k vs 2^k ± 1 in cold-start subprocesses.

    Sympy caches sieve tables across calls within a process, which would mask
    any genuine x-dependent cost difference. Each (k, x_label, trial) gets a
    fresh subprocess.
    """
    out = {}
    for k in anchors:
        x_pow = 1 << k
        triples = {"x_pow": x_pow, "x_pow_minus_1": x_pow - 1,
                   "x_pow_plus_1": x_pow + 1}
        timings = {}
        for label, x in triples.items():
            ts = []
            v_seen = None
            for _ in range(trials):
                t, v = _subprocess_time_primepi(x)
                ts.append(t)
                v_seen = v
            ts.sort()
            timings[label] = {
                "value": v_seen,
                "median_s": ts[len(ts) // 2],
                "min_s": ts[0],
                "max_s": ts[-1],
                "all_s": ts,
            }
        # ratio
        med_pow = timings["x_pow"]["median_s"]
        med_neighbour = 0.5 * (timings["x_pow_minus_1"]["median_s"]
                               + timings["x_pow_plus_1"]["median_s"])
        ratio = med_pow / med_neighbour if med_neighbour > 0 else float("nan")
        out[k] = {
            "timings": timings,
            "median_x_pow": med_pow,
            "median_neighbour": med_neighbour,
            "ratio_pow_over_neighbour": ratio,
        }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cost-anchors", default="15,18,20",
                    help="comma list of k values for HKM cost test")
    ap.add_argument("--cost-trials", type=int, default=5,
                    help="trials per (k, x_label) cell")
    ap.add_argument("--verify-up-to", type=int, default=20,
                    help="verify table vs sympy for k <= this (slow)")
    ap.add_argument("--out", default=str(Path(__file__).parent),
                    help="output directory")
    args = ap.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("[1/4] Verifying table vs sympy.primepi...")
    verify = verify_via_sympy(args.verify_up_to)
    n_ok = sum(1 for v in verify.values() if v["agree"])
    n_check = len(verify)
    print(f"  agreement: {n_ok}/{n_check}")

    print("[2/4] Computing Li(2^k), R(2^k), and residuals...")
    residuals = explicit_residuals()

    # Primary normalised signals
    c_R_seq = [residuals[k]["c_R"] for k in range(1, 57)]
    c_li_seq = [residuals[k]["c_li"] for k in range(1, 57)]
    sign_R_seq = [residuals[k]["sign_r_R"] for k in range(1, 57)]
    bit_R_seq = [0 if s < 0 else 1 for s in sign_R_seq]
    n_pos_R = sum(bit_R_seq)
    n_neg_R = len(bit_R_seq) - n_pos_R

    raw_path = out_dir / "raw_data.json"
    with open(raw_path, "w") as f:
        json.dump({"verify": verify, "residuals": residuals}, f, indent=2)
    print(f"  wrote {raw_path}")

    print("[3/4] Monte Carlo baselines (BM and autocorr distributions)...")
    n_pos_R_pre = sum(1 for k in range(1, 57) if residuals[k]["sign_r_R"] >= 0)
    mc = monte_carlo_baselines(56, n_pos_R_pre, n_trials=4000)
    print(f"  random-binary BM: mean={mc['bm_mean']:.2f} std={mc['bm_std']:.2f} "
          f"[p05={mc['bm_p05']}, p95={mc['bm_p95']}]")
    print(f"  Gaussian-iid |ac_lag1|: p95={mc['abs_ac_lag1_p95']:.3f} "
          f"p99={mc['abs_ac_lag1_p99']:.3f} p999={mc['abs_ac_lag1_p999']:.3f}")

    print("[3b/4] Pseudorandomness tests on c_R(k) and sign(r_R(k))...")
    ac_R = {lag: autocorr(c_R_seq, lag) for lag in range(1, 11)}
    ac_li = {lag: autocorr(c_li_seq, lag) for lag in range(1, 6)}
    bm = berlekamp_massey_gf2(bit_R_seq)
    rt = runs_test(bit_R_seq)

    mean_R = sum(c_R_seq) / len(c_R_seq)
    pr = {
        "n": len(c_R_seq),
        "c_R_min": min(c_R_seq),
        "c_R_max": max(c_R_seq),
        "c_R_mean": mean_R,
        "c_R_std": (sum((x - mean_R) ** 2 for x in c_R_seq)
                    / len(c_R_seq)) ** 0.5,
        "autocorr_R_lag_1": ac_R[1],
        "autocorr_R_lag_2": ac_R[2],
        "autocorr_R_lag_3": ac_R[3],
        "autocorr_R_lag_4": ac_R[4],
        "autocorr_R_lag_5": ac_R[5],
        "autocorr_R_lag_6": ac_R[6],
        "autocorr_R_lag_7": ac_R[7],
        "autocorr_R_lag_8": ac_R[8],
        "autocorr_R_lag_9": ac_R[9],
        "autocorr_R_lag_10": ac_R[10],
        "autocorr_li_lag_1": ac_li[1],  # for comparison; smooth-trend confounded
        "bm_linear_complexity_sign_R": bm,
        "n_positive_R": n_pos_R,
        "n_negative_R": n_neg_R,
        "runs_test_R": rt,
        "skewes_trivial_R": (n_pos_R == 0 or n_neg_R == 0),
    }
    print(f"  c_R stats: mean={pr['c_R_mean']:.5f} std={pr['c_R_std']:.5f}")
    print(f"  c_li (smooth-confounded): mean={sum(c_li_seq)/len(c_li_seq):.4f}")
    print(f"  autocorr_R lag1..10: {[round(ac_R[l], 3) for l in range(1, 11)]}")
    print(f"  autocorr_li lag1: {ac_li[1]:.4f} (expected high — smooth trend)")
    print(f"  BM linear complexity (sign_R seq, GF2): {bm}")
    print(f"  signs(r_R): pos={n_pos_R} neg={n_neg_R}  trivial={pr['skewes_trivial_R']}")

    print("[4/4] HKM cost test on dyadic vs neighbour...")
    anchors = [int(x) for x in args.cost_anchors.split(",") if x.strip()]
    cost = hkm_cost_test(anchors, args.cost_trials)
    cost_path = out_dir / "hkm_cost.json"
    with open(cost_path, "w") as f:
        json.dump(cost, f, indent=2)
    print(f"  wrote {cost_path}")
    for k, d in cost.items():
        print(f"  k={k}: T(2^k)={d['median_x_pow']:.4f}s, "
              f"T(neighbour)={d['median_neighbour']:.4f}s, "
              f"ratio={d['ratio_pow_over_neighbour']:.3f}")

    # Verdict on F-claims (refined: Riemann R baseline removes smooth trend).
    # F1: tightened post-MC. The naïve "BM <= log²(n) ≈ 33" threshold is
    # not informative at n=56 because random-baseline E[BM] ≈ 28 sits
    # inside [0, 33]. We re-bind F1 to "BM at most p05 of the
    # random-shuffle distribution, given the empirical sign split", i.e.
    # 5%-tail compressibility.
    f1 = ((not pr["skewes_trivial_R"]) and
          (pr["bm_linear_complexity_sign_R"] <= mc["bm_p05"]))
    # We pre-registered |lag-1 ac| >= 0.3 against a Li-baseline in v1; with
    # the corrected R-baseline we instead test the *strongest* lag from
    # 1..10. Rationale: zero-driven oscillation at γ_1 log 2 = 9.797 rad has
    # period ≈ 2π/9.797 = 0.641 in k, so γ_1 alone gives no integer-lag
    # peak; γ_1 + γ_2 mixing might. Use max-magnitude as an honest proxy
    # for "structure beyond noise".
    max_ac_R = max(abs(pr[f"autocorr_R_lag_{l}"]) for l in range(1, 11))
    max_ac_R_lag = max(range(1, 11),
                       key=lambda l: abs(pr[f"autocorr_R_lag_{l}"]))
    # F2 threshold post-MC: against MC p99 floor for a single lag.
    # We test 10 lags so we apply Bonferroni at p99 (10× single-lag p99).
    f2 = max_ac_R >= mc["abs_ac_lag1_p999"]
    ratios = [d["ratio_pow_over_neighbour"] for d in cost.values()]
    avg_ratio = sum(ratios) / len(ratios) if ratios else float("nan")
    f3 = avg_ratio <= 0.5

    summary = {
        "f_claim_1_BM_below_random_p05": {
            "holds": bool(f1),
            "bm_value": pr["bm_linear_complexity_sign_R"],
            "random_baseline_mean": mc["bm_mean"],
            "random_baseline_std": mc["bm_std"],
            "random_baseline_p05": mc["bm_p05"],
            "skewes_trivial": pr["skewes_trivial_R"],
            "n_pos": pr["n_positive_R"], "n_neg": pr["n_negative_R"],
        },
        "f_claim_2_max_abs_autocorr_above_MC_p999": {
            "holds": bool(f2),
            "max_abs_ac": max_ac_R,
            "max_abs_ac_lag": max_ac_R_lag,
            "MC_abs_ac_p99": mc["abs_ac_lag1_p99"],
            "MC_abs_ac_p999": mc["abs_ac_lag1_p999"],
            "all_lags": {l: pr[f"autocorr_R_lag_{l}"] for l in range(1, 11)},
        },
        "f_claim_3_dyadic_speedup_ratio_le_0_5": {
            "holds": bool(f3),
            "avg_ratio": avg_ratio,
            "per_k_ratios": dict(zip(anchors, ratios)),
        },
        "verdict": ("PARTIAL_POSITIVE_SURPRISE" if (f1 or f2 or f3)
                    else "B_NEGATIVE_no_dyadic_amortisation"),
        "pseudorandomness_tests": pr,
        "monte_carlo": mc,
    }
    sum_path = out_dir / "summary.json"
    with open(sum_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"  wrote {sum_path}")
    print(f"\nVERDICT: {summary['verdict']}")
    print(f"  F1 (BM at or below MC p05={mc['bm_p05']}): "
          f"{'HOLDS' if f1 else 'FAILS'}  "
          f"(value={pr['bm_linear_complexity_sign_R']}, "
          f"MC mean={mc['bm_mean']:.2f}±{mc['bm_std']:.2f})")
    print(f"  F2 (max_lag |ac_R| >= MC p999={mc['abs_ac_lag1_p999']:.3f}): "
          f"{'HOLDS' if f2 else 'FAILS'}  "
          f"(max={max_ac_R:.3f} at lag {max_ac_R_lag})")
    print(f"  F3 (dyadic speedup ratio <= 0.5): "
          f"{'HOLDS' if f3 else 'FAILS'}  (avg_ratio={avg_ratio:.3f})")


if __name__ == "__main__":
    main()
