"""
Thread 9 (P4 / k-tuple narrow-window count batched on x), slot 1: baseline.

Question: does Thread 5's correlation dichotomy (S224, 33x speedup at M=64
between correlated and uncorrelated batched pi(x) queries) transpose to
narrow-window k-tuple counts pi_H(x; w) = #{n in [x, x+w] : n+h prime for
all h in H}?

Hypothesis: NO — per-query cost is ALREADY polylog (window sieve over
w = polylog(x)), so the speedup ceiling is bounded constant (sieve-setup
sharing across batched-x), not 33x.

Slot 1 establishes:
  - per-x amortised cost vs M for M in {1, 4, 16, 64};
  - across THREE distributions: uncorrelated (Uniform[x_max/2, x_max]),
    correlated-W=M (consecutive integers), correlated-W=polylog;
  - across THREE k-tuples: k=2 (h=(0,2) twin), k=3 (h=(0,2,6) admissible),
    k=4 (h=(0,2,6,8) admissible);
  - at x_max in {10^6, 10^7}.
  - HL-approximation cost (singular series x li_2(x) — polylog by construction).

Empirical falsifiable hypotheses (recorded in slot1_baseline_results.md):
  H1: T_batched / T_naive --> 1 + O(1) as M --> infty (bounded constant
      saving), NOT --> 0 (Thread-5-shape).
  H2: speedup ratio is independent of x_max within the polylog-window
      regime (sieve-share scales as O(M / (M+1)), kernel polylog).
  H3: HL_H(x) cost is polylog and matches empirical k-tuple count up
      to typical error ~ sqrt(C(H) * w / log^k x) (Hardy-Littlewood).
"""

import csv
import math
import os
import random
import sys
import time

# ----------------------------------------------------------------------
# Sieve primitives
# ----------------------------------------------------------------------

def isqrt(n):
    return math.isqrt(n)


def primes_upto(N):
    """Sieve of Eratosthenes; returns list of primes <= N."""
    s = bytearray(b"\x01") * (N + 1)
    s[0] = 0
    if N >= 1:
        s[1] = 0
    for i in range(2, isqrt(N) + 1):
        if s[i]:
            step = i
            for j in range(i * i, N + 1, step):
                s[j] = 0
    return [i for i in range(2, N + 1) if s[i]]


def seg_sieve(lo, hi, primes):
    """Segmented sieve; returns bytearray of length (hi - lo + 1) where
    sieve[i] = 1 iff (lo + i) is prime. Assumes lo >= 2.
    """
    n = hi - lo + 1
    s = bytearray(b"\x01") * n
    for p in primes:
        if p * p > hi:
            break
        # Smallest multiple of p in [lo, hi] is max(p*p, ceil(lo/p)*p)
        start = max(p * p, ((lo + p - 1) // p) * p)
        for j in range(start, hi + 1, p):
            s[j - lo] = 0
    return s


# ----------------------------------------------------------------------
# k-tuple counting
# ----------------------------------------------------------------------

def count_ktuples_in(s, base_off, w, H):
    """Count #{i in [0, w] : s[base_off + i + h] = 1 for all h in H}."""
    cnt = 0
    H_local = H  # local to avoid attr lookup
    for i in range(w + 1):
        idx = base_off + i
        ok = True
        for h in H_local:
            if not s[idx + h]:
                ok = False
                break
        if ok:
            cnt += 1
    return cnt


def naive_ktuples(x_list, w, H, primes):
    """Naive: independent segmented sieve at each x_i."""
    h_max = max(H)
    counts = []
    for x in x_list:
        s = seg_sieve(x, x + w + h_max, primes)
        counts.append(count_ktuples_in(s, 0, w, H))
    return counts


def batched_ktuples(x_list, w, H, primes):
    """Batched: single segmented sieve over [min(x), max(x) + w + h_max]."""
    h_max = max(H)
    lo = min(x_list)
    hi = max(x_list) + w + h_max
    s = seg_sieve(lo, hi, primes)
    return [count_ktuples_in(s, x - lo, w, H) for x in x_list]


# ----------------------------------------------------------------------
# Hardy-Littlewood singular series for k-tuple H
# ----------------------------------------------------------------------

def singular_series(H, P_max=10000):
    """Approximate C(H) = product_p (1 - nu_H(p)/p) / (1 - 1/p)^k where
    nu_H(p) = #{distinct h mod p}. Truncates at P_max (sufficient for
    H with small entries).
    """
    H = list(H)
    k = len(H)
    primes = primes_upto(P_max)
    val = 1.0
    for p in primes:
        residues = set(h % p for h in H)
        nu = len(residues)
        val *= (1.0 - nu / p) / ((1.0 - 1.0 / p) ** k)
    return val


def li_2(x):
    """li_2(x) = integral_2^x dt / log^2(t). Numerical via Simpson, modest
    precision (sufficient for slot-1 comparison)."""
    if x < 2:
        return 0.0
    n = 4096
    h = (x - 2) / n
    f = lambda t: 1.0 / (math.log(t) ** 2)
    s = f(2) + f(x)
    for i in range(1, n, 2):
        s += 4 * f(2 + i * h)
    for i in range(2, n, 2):
        s += 2 * f(2 + i * h)
    return s * h / 3


def hl_estimate_window(x, w, H):
    """Hardy-Littlewood estimate for k-tuple count in [x, x+w]:
        pi_H(x; w) ~ C(H) * (li_k(x+w) - li_k(x))
                   ~ C(H) * w / log^k(x)   for w << x.
    Uses the leading term w / log^k(x).
    """
    k = len(H)
    return singular_series(H) * w / (math.log(x) ** k)


# ----------------------------------------------------------------------
# Verification helpers
# ----------------------------------------------------------------------

def _try_isprime():
    try:
        from sympy import isprime as _sympy_isprime
        return _sympy_isprime
    except Exception:
        return None


def verify_small(seed=0):
    """Sanity check: compare naive_ktuples to brute force at small x."""
    isp = _try_isprime()
    if isp is None:
        print("  (sympy unavailable; skipping verify)")
        return True
    H = (0, 2)
    x = 1000
    w = 50
    primes = primes_upto(40)
    s = seg_sieve(x, x + w + max(H), primes)
    sieve_cnt = count_ktuples_in(s, 0, w, H)
    brute_cnt = 0
    for n in range(x, x + w + 1):
        if all(isp(n + h) for h in H):
            brute_cnt += 1
    if sieve_cnt != brute_cnt:
        print(f"  VERIFY FAIL at x={x}, w={w}, H={H}: sieve={sieve_cnt}, brute={brute_cnt}")
        return False
    print(f"  verify OK at x={x}, w={w}, H={H}: count = {sieve_cnt}")
    return True


# ----------------------------------------------------------------------
# Experiment harness
# ----------------------------------------------------------------------

def time_call(fn, *args, n_repeat=2):
    best = float("inf")
    last = None
    for _ in range(n_repeat):
        t0 = time.perf_counter()
        last = fn(*args)
        dt = time.perf_counter() - t0
        if dt < best:
            best = dt
    return best, last


def run_block(x_max, w, H, primes_full, distributions, M_list, k_label):
    """Run one (x_max, H) block over the list of distributions and Ms.
    Returns rows: (k_label, x_max, w, H_str, distribution, M, T_naive_s,
    T_batched_s, T_amort_naive_s, T_amort_batched_s, ratio_batched_naive,
    counts_match).
    """
    rows = []
    H_str = "(" + ",".join(str(h) for h in H) + ")"
    log_x = math.log(x_max)
    w_polylog = max(2, int(log_x ** 2))  # polylog window size, e.g. 16^2 = 261
    for dist_label, dist_factory in distributions:
        for M in M_list:
            random.seed(hash((x_max, dist_label, M, k_label)) & 0xFFFFFFFF)
            x_list = dist_factory(x_max, M, w_polylog)
            x_list = sorted(x_list)
            # Naive
            T_naive, counts_naive = time_call(naive_ktuples, x_list, w, H, primes_full)
            # Batched
            T_batched, counts_batched = time_call(batched_ktuples, x_list, w, H, primes_full)
            # Compare
            match = (counts_naive == counts_batched)
            T_amort_naive = T_naive / M
            T_amort_batched = T_batched / M
            ratio = T_batched / T_naive
            rows.append({
                "k_label": k_label,
                "x_max": x_max,
                "w": w,
                "H": H_str,
                "distribution": dist_label,
                "M": M,
                "T_naive_s": T_naive,
                "T_batched_s": T_batched,
                "T_amort_naive_ms": T_amort_naive * 1000,
                "T_amort_batched_ms": T_amort_batched * 1000,
                "ratio_batched_naive": ratio,
                "counts_match": int(match),
            })
            tag = "OK" if match else "MISMATCH"
            print(f"  {dist_label:20s} M={M:>3} naive={T_naive*1000:>7.2f}ms "
                  f"batched={T_batched*1000:>7.2f}ms  amort_naive={T_amort_naive*1000:>6.3f}ms/x "
                  f"amort_batched={T_amort_batched*1000:>6.3f}ms/x  ratio={ratio:>5.3f}  [{tag}]")
    return rows


def dist_uncorrelated(x_max, M, w_polylog):
    return [random.randint(x_max // 2, x_max) for _ in range(M)]


def dist_correlated_consecutive(x_max, M, w_polylog):
    # Tight cluster of M consecutive offsets.
    return [x_max + i for i in range(M)]


def dist_correlated_polylog(x_max, M, w_polylog):
    # M offsets sampled from [0, w_polylog].
    if M > w_polylog + 1:
        # Can't sample M distinct in window of width w_polylog
        deltas = list(range(M))
    else:
        deltas = sorted(random.sample(range(w_polylog + 1), M))
    return [x_max + d for d in deltas]


def main():
    out_dir = os.path.dirname(os.path.abspath(__file__))
    out_csv = os.path.join(out_dir, "slot1_baseline.csv")
    out_hl_csv = os.path.join(out_dir, "slot1_hl_compare.csv")

    print("=" * 72)
    print("Thread 9 / P4 slot 1: k-tuple narrow-window batched-on-x baseline")
    print("=" * 72)
    print()

    # Sanity check
    print("[verify] sanity check on small case")
    verify_small()
    print()

    # H-tuples
    H_list = [
        ("k2_twin", (0, 2)),
        ("k3_admiss", (0, 2, 6)),
        ("k4_admiss", (0, 2, 6, 8)),
    ]

    # Distributions
    distributions = [
        ("uncorrelated", dist_uncorrelated),
        ("corr_w=M", dist_correlated_consecutive),
        ("corr_w=polylog", dist_correlated_polylog),
    ]
    M_list = [1, 4, 16, 64]

    rows = []
    hl_rows = []
    for x_max in [10**6, 10**7]:
        log_x = math.log(x_max)
        w = max(50, int(log_x ** 2))  # polylog window, ensure not too tiny
        h_max_global = 8  # max h across all H tested
        N_sqrt = isqrt(x_max + w + h_max_global) + 1
        print(f"--- x_max = 10^{int(round(math.log10(x_max)))}, w = {w}, "
              f"sqrt-bound primes <= {N_sqrt} ---")
        primes_full = primes_upto(N_sqrt)
        for k_label, H in H_list:
            print(f"  H = {H}  ({k_label})")
            block = run_block(x_max, w, H, primes_full, distributions, M_list, k_label)
            rows.extend(block)
            # HL-comparison: average count over a small uncorrelated batch
            random.seed(7777 + x_max + sum(H))
            x_check_list = [random.randint(x_max // 2, x_max) for _ in range(20)]
            counts = naive_ktuples(x_check_list, w, H, primes_full)
            mean_count = sum(counts) / len(counts)
            hl_pred = hl_estimate_window(x_max, w, H)
            err = mean_count - hl_pred
            sigma_pois = math.sqrt(max(hl_pred, 1e-12))
            hl_rows.append({
                "k_label": k_label,
                "x_max": x_max,
                "w": w,
                "H": str(H),
                "C_H": singular_series(H),
                "hl_pred": hl_pred,
                "mean_emp": mean_count,
                "err": err,
                "sigma_pois": sigma_pois,
                "err_over_sigma": err / sigma_pois,
                "n_samples": len(counts),
            })
            print(f"    HL pred = {hl_pred:.4f}, mean emp = {mean_count:.4f}, "
                  f"err/sigma_pois = {err/sigma_pois:.3f}")
        print()

    # Write CSVs
    with open(out_csv, "w", newline="") as f:
        if rows:
            w_csv = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w_csv.writeheader()
            w_csv.writerows(rows)
    print(f"wrote {out_csv}  ({len(rows)} rows)")

    with open(out_hl_csv, "w", newline="") as f:
        if hl_rows:
            w_csv = csv.DictWriter(f, fieldnames=list(hl_rows[0].keys()))
            w_csv.writeheader()
            w_csv.writerows(hl_rows)
    print(f"wrote {out_hl_csv}  ({len(hl_rows)} rows)")

    # Summary: speedup at M = 64 for each distribution / k
    print()
    print("=" * 72)
    print("SLOT 1 SUMMARY — speedup ratio T_batched / T_naive at M = 64")
    print("=" * 72)
    for x_max in [10**6, 10**7]:
        for k_label, _ in H_list:
            for dist_label, _ in distributions:
                cands = [r for r in rows
                         if r["x_max"] == x_max and r["k_label"] == k_label
                         and r["distribution"] == dist_label and r["M"] == 64]
                if cands:
                    r = cands[0]
                    print(f"  x_max=10^{int(round(math.log10(x_max)))} "
                          f"{k_label:>10s} {dist_label:>16s} M=64  "
                          f"ratio = {r['ratio_batched_naive']:.3f}  "
                          f"(amort batched = {r['T_amort_batched_ms']:.3f} ms/x)")
    print()
    print("Hypothesis check:")
    print("  H1 (Thread-5-shape FAILS to transpose): expect ratios > ~0.3 at M=64,")
    print("      i.e., bounded constant saving (compare vs Thread 5 ratio = 0.030 = 33x).")
    print()


if __name__ == "__main__":
    main()
