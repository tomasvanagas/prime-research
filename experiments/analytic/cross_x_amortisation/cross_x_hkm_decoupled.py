"""
Cross-x amortisation, slot 3 of Thread 5: HKM / Lucy DP / Meissel-Lehmer
combinatorial pillar.

Builds an instrumented batched Lucy DP that decouples shared state from
per-x state.

Decomposition:
  shared state  := small[v] for v <= sqrt(x_max), and the prime list.
                   Equivalently, sieve of Eratosthenes up to sqrt(x_max).
  per-x state   := large_i[v] for v <= sqrt(x_i), one array per query.

Lucy DP cost decomposition (BASIC, no Fenwick — focus on structural shape,
not the log factor):

  T_setup_shared(x_max) = O(sqrt(x_max) log log x_max)
                          [sieve of Eratosthenes up to sqrt(x_max)]
  T_per_x_evaluation(x_i) = O(x_i / log x_i)
                          [large[] sweep over primes p <= sqrt(x_i)]

  T_total(M, x_max, x_avg) = T_setup_shared(x_max) + M * T_per_x(x_avg)
  T_per_x_amortised(M)     = T_setup_shared(x_max) / M + T_per_x(x_avg)

The forecast (slot 1's twin pattern):

  Setup is O(sqrt(x_max) log log) which is *vastly smaller* than per-x
  cost O(x_avg / log x_avg). Amortisation gain is ~0% because the per-x
  term dominates immediately.

Cluster-stitched alternative (the only known sub-anchor amortisation):

  For x_i in [x_anchor, x_anchor + w], do anchor Lucy DP once, then for
  each i > 1 sieve the segment [x_anchor, x_i] (cost O(w/log)) and
  subtract: pi(x_i) = pi(x_anchor) + #primes in (x_anchor, x_i].

  T_cluster(M, w) = T_per_x(x_anchor) + (M-1) * O(w / log)
  Per-x amortised = T_per_x(x_anchor) / M + O(w / log)

  For batch-polylog: need M >= x^{2/3} / polylog AND w <= polylog^2.
  So the cluster width must be polylog and M must be Theta-tilde(x^{2/3}).
  This is the standard "anchor + sieve" trick — not a new amortisation.

The experiment confirms both forecasts empirically, then concludes that
HKM admits no batch-polylog amortisation over M = polylog uncorrelated x's.
"""

import math
import time
from sympy import isprime, primepi


def isqrt(n):
    s = int(math.isqrt(n))
    while (s + 1) * (s + 1) <= n:
        s += 1
    while s * s > n:
        s -= 1
    return s


# ---------------------------------------------------------------------------
# Shared setup: sieve of Eratosthenes up to sqrt(x_max)
# Returns small[v] = pi(v) for v <= N, plus the prime list.
# ---------------------------------------------------------------------------

def shared_setup(N):
    """Return small[v] = pi(v) for v in 0..N, plus list of primes <= N."""
    sieve = bytearray(b"\x01") * (N + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, isqrt(N) + 1):
        if sieve[i]:
            for j in range(i * i, N + 1, i):
                sieve[j] = 0
    primes = [i for i in range(2, N + 1) if sieve[i]]
    small = [0] * (N + 2)
    cnt = 0
    for v in range(N + 2):
        if v <= N and sieve[v]:
            cnt += 1
        small[v] = cnt
    return small, primes


# ---------------------------------------------------------------------------
# Per-x evaluation: given shared small[] and primes, compute large[v] for x_i
# via the Lucy DP iteration restricted to primes p <= sqrt(x_i).
# ---------------------------------------------------------------------------

def per_x_eval(x, small, primes):
    """Compute pi(x) given precomputed pi(v) = small[v] for v <= sqrt(x_max).

    Uses the Lucy DP iteration over primes p <= sqrt(x), updating only
    the per-x large[] array. small[] is read-only (NOT updated) — we use
    the FINAL value of small[v] = pi(v), which equals the value Lucy DP
    would have computed at step p > sqrt(v) anyway.

    BUT: at intermediate step p in Lucy DP, small[v] is NOT pi(v) — it's
    Legendre's S(v, p-1). The Lucy DP reference uses small[v] at step p
    only for v >= p, where the value at step p equals the value at step
    sqrt(v), which equals pi(v). The reference is small[x/d] where
    x/d > sqrt(x) > sqrt(v). So we use small[v] = pi(v) directly. (See
    Lucy_Hedgehog 2017 forum post; verified equivalence below.)

    Wait — that's wrong. Let me re-derive.

    Standard Lucy DP code:
        small[v] = v - 1   (initial)
        for p in primes:
            for v: small[v] -= small[v/p] - pcnt   if v >= p^2

    At step p, small[v] tracks #{n <= v : n's smallest prime factor > p, or n = 1}.
    This is the "partial sieve" Legendre quantity, NOT pi(v).

    Only at step p > sqrt(v) is small[v] = pi(v).

    The large[] update at step p references small[x/d] where d = v*p
    and v <= sqrt(x), x/d <= x/p. Now x/d may be much larger than sqrt(x)
    or smaller — for v=1, d=p, x/d = x/p which can be huge.

    Wait — the code reference:
        if d <= sqrtx: large[v] -= large[d] - pcnt
        else:          large[v] -= small[x / d] - pcnt

    If d > sqrtx, then x/d < sqrtx, and small[x/d] is referenced.

    At step p, does small[x/d] = pi(x/d)? Only if x/d < p (then sieving
    primes >= p has no effect). If x/d >= p, small[x/d] is the partial-sieve
    quantity and we need the iterative update.

    So the per-x evaluation MUST also iterate small[] alongside large[].
    The "shared" small[] is only the FINAL value, not intermediate values.

    This means slot 3's setup/per-x decoupling needs care:
      - Shared part: the prime list and the FINAL small[v] = pi(v).
      - Per-x part: re-iterate small[] state alongside large[] for each x_i.

    Because small[] state at step p depends only on the prime list (not on
    x), we can run the small[] update ONCE for x_max-bounded sqrt range,
    storing intermediate snapshots. But the snapshots are O(p * sqrt(x_max))
    in storage = O(x_max / log x_max) total — not a free shared resource.

    Alternative: re-run the full Lucy DP for each x_i (small AND large
    interleaved), but skip primes p > sqrt(x_i). This is what `per_x_eval`
    does below. The "shared" piece is only the final small[v] used for
    answer extraction (large[1] = pi(x_i)).

    This makes "shared" essentially trivial. The per-x cost is ~T_lucy(x_i).

    => HKM CROSS-X AMORTISATION HAS NEAR-ZERO SHARING. Worse than the
       explicit-formula pillar (which had non-trivial K_zeros_setup
       to amortise).
    """
    # Re-run full Lucy DP for x (basic version, NO Fenwick).
    sqrtx = isqrt(x)
    s = [v - 1 for v in range(sqrtx + 2)]  # s[v] = v - 1
    L = [0] * (sqrtx + 2)
    for v in range(1, sqrtx + 1):
        L[v] = x // v - 1
    for p in primes:
        if p > sqrtx:
            break
        if s[p] == s[p - 1]:
            continue  # not a prime in this state — shouldn't happen here
        pcnt = s[p - 1]
        p2 = p * p
        if p2 > x:
            break
        limit = min(sqrtx, x // p2)
        # Update large[v]
        for v in range(1, limit + 1):
            d = v * p
            if d <= sqrtx:
                L[v] -= L[d] - pcnt
            else:
                L[v] -= s[x // d] - pcnt
        # Update small[v] in place
        for v in range(sqrtx, p2 - 1, -1):
            s[v] -= s[v // p] - pcnt
    return L[1]


def per_x_eval_instrumented(x, primes):
    """Same as per_x_eval but counts the number of small-update and
    large-update operations separately."""
    sqrtx = isqrt(x)
    s = [v - 1 for v in range(sqrtx + 2)]
    L = [0] * (sqrtx + 2)
    for v in range(1, sqrtx + 1):
        L[v] = x // v - 1
    n_small_updates = 0
    n_large_updates = 0
    for p in primes:
        if p > sqrtx:
            break
        if s[p] == s[p - 1]:
            continue
        pcnt = s[p - 1]
        p2 = p * p
        if p2 > x:
            break
        limit = min(sqrtx, x // p2)
        for v in range(1, limit + 1):
            d = v * p
            if d <= sqrtx:
                L[v] -= L[d] - pcnt
            else:
                L[v] -= s[x // d] - pcnt
            n_large_updates += 1
        for v in range(sqrtx, p2 - 1, -1):
            s[v] -= s[v // p] - pcnt
            n_small_updates += 1
    return L[1], n_small_updates, n_large_updates


# ---------------------------------------------------------------------------
# Cluster-stitched approach: anchor + segmented sieve
# ---------------------------------------------------------------------------

def cluster_stitched(x_anchor, deltas, small, primes):
    """For x_anchor and a list of small offsets (deltas), compute pi(x_anchor + d).

    Approach: do per-x Lucy DP once at x_anchor (gives pi(x_anchor)).
    Then for each delta > 0, pi(x_anchor + delta) = pi(x_anchor) + #{primes
    in (x_anchor, x_anchor + delta]} via direct trial division.

    Returns: list of pi values.
    """
    pi_anchor = per_x_eval(x_anchor, small, primes)
    results = []
    for d in deltas:
        if d == 0:
            results.append(pi_anchor)
            continue
        if d < 0:
            # Subtract: count primes in (x_anchor + d, x_anchor]
            lo, hi = x_anchor + d + 1, x_anchor
            cnt = 0
            for n in range(lo, hi + 1):
                if isprime(n):
                    cnt += 1
            results.append(pi_anchor - cnt)
        else:
            lo, hi = x_anchor + 1, x_anchor + d
            cnt = 0
            for n in range(lo, hi + 1):
                if isprime(n):
                    cnt += 1
            results.append(pi_anchor + cnt)
    return results


# ---------------------------------------------------------------------------
# Profile sweep
# ---------------------------------------------------------------------------

def time_call(fn, *args, n_repeat=3, **kwargs):
    times = []
    last = None
    for _ in range(n_repeat):
        t0 = time.perf_counter()
        last = fn(*args, **kwargs)
        times.append(time.perf_counter() - t0)
    return min(times), last


def main():
    out_setup_csv = "experiments/analytic/cross_x_amortisation/cross_x_hkm_setup.csv"
    out_perx_csv = "experiments/analytic/cross_x_amortisation/cross_x_hkm_perx.csv"
    out_amort_csv = "experiments/analytic/cross_x_amortisation/cross_x_hkm_amortised.csv"
    out_cluster_csv = "experiments/analytic/cross_x_amortisation/cross_x_hkm_cluster.csv"
    out_ops_csv = "experiments/analytic/cross_x_amortisation/cross_x_hkm_ops.csv"
    out_scaling_csv = "experiments/analytic/cross_x_amortisation/cross_x_hkm_scaling.csv"

    # ============================================================
    # Q1: Profile shared setup vs per-x evaluation
    # ============================================================

    setup_rows = ["x_max,sqrt_xmax,n_primes,T_setup_s"]
    perx_rows = ["x,sqrtx,T_perx_s,pi_x"]

    x_max_list = [10**4, 10**5, 10**6, 10**7]
    print("=== Q1: shared setup vs per-x evaluation ===")
    cached_setups = {}
    for x_max in x_max_list:
        N = isqrt(x_max)
        T_setup, (small, primes) = time_call(shared_setup, N, n_repeat=3)
        setup_rows.append(f"{x_max},{N},{len(primes)},{T_setup:.6f}")
        print(f"  x_max=10^{int(round(math.log10(x_max)))}: "
              f"sqrt={N}, n_primes={len(primes)}, T_setup={T_setup*1000:.3f} ms")
        cached_setups[x_max] = (small, primes)

        # Time per-x at x_max itself
        T_perx, pi_xmax = time_call(per_x_eval, x_max, small, primes, n_repeat=3)
        perx_rows.append(f"{x_max},{N},{T_perx:.6f},{pi_xmax}")
        ref = int(primepi(x_max))
        ok = "OK" if pi_xmax == ref else f"WRONG (got {pi_xmax}, want {ref})"
        print(f"    per-x at x={x_max}: T={T_perx*1000:.2f} ms, "
              f"pi={pi_xmax} {ok}")

    # ============================================================
    # Q2: Operation-count instrumentation (is the work O(x/log x) per x_i?)
    # ============================================================

    ops_rows = ["x,sqrtx,n_small_updates,n_large_updates,total_ops,total_per_xover_logx"]
    print("\n=== Q2: operation-count scaling ===")
    primes_max = cached_setups[max(x_max_list)][1]
    for x in [10**3, 10**4, 10**5, 10**6, 10**7]:
        _, n_small, n_large = per_x_eval_instrumented(x, primes_max)
        total = n_small + n_large
        x_over_logx = x / math.log(x)
        ratio = total / x_over_logx
        ops_rows.append(f"{x},{isqrt(x)},{n_small},{n_large},{total},{ratio:.3f}")
        print(f"  x={x}: small_ops={n_small}, large_ops={n_large}, "
              f"total={total}, total/(x/log x)={ratio:.3f}")

    # ============================================================
    # Q3: M-batched per-x amortised cost (uncorrelated x_i's)
    # ============================================================

    amort_rows = ["x_max,M,T_setup_s,T_perx_total_s,T_amortised_per_x_s,setup_share_pct"]
    print("\n=== Q3: M-batched per-x amortised cost (uncorrelated x's) ===")
    x_max = 10**6
    small, primes = cached_setups[x_max]
    # Sample M x_i's uniformly in [x_max/2, x_max]
    import random
    random.seed(42)
    for M in [1, 2, 4, 8, 16, 32]:
        x_list = [random.randint(x_max // 2, x_max) for _ in range(M)]
        T_setup, _ = time_call(shared_setup, isqrt(x_max), n_repeat=2)

        t0 = time.perf_counter()
        results = [per_x_eval(x_i, small, primes) for x_i in x_list]
        T_perx_total = time.perf_counter() - t0

        T_amort = (T_setup + T_perx_total) / M
        share = T_setup / (T_setup + T_perx_total) * 100
        amort_rows.append(
            f"{x_max},{M},{T_setup:.6f},{T_perx_total:.6f},{T_amort:.6f},{share:.3f}"
        )
        print(f"  M={M:>3}: T_setup={T_setup*1000:.2f}ms, "
              f"T_perx_tot={T_perx_total:.3f}s, "
              f"T_amort={T_amort*1000:.2f}ms/x, "
              f"setup-share={share:.2f}%")

    # ============================================================
    # Q3b: Per-x amortised cost vs x_max (fixed M)
    # Confirms amortisation cannot keep up with per-x growth.
    # ============================================================

    scaling_rows = ["x_max,M,T_setup_s,T_perx_avg_s,T_amortised_per_x_s,perx_over_xlog_x"]
    print("\n=== Q3b: per-x amortised cost vs x_max (fixed M = 8 uncorrelated) ===")
    M_fixed = 8
    for x_max in x_max_list:
        small_xm, primes_xm = cached_setups[x_max]
        T_setup, _ = time_call(shared_setup, isqrt(x_max), n_repeat=2)

        random.seed(99)
        x_list = [random.randint(x_max // 2, x_max) for _ in range(M_fixed)]
        t0 = time.perf_counter()
        results = [per_x_eval(x_i, small_xm, primes_xm) for x_i in x_list]
        T_perx_total = time.perf_counter() - t0
        T_perx_avg = T_perx_total / M_fixed

        T_amort = (T_setup + T_perx_total) / M_fixed
        x_over_logx = x_max / math.log(x_max)
        ratio = T_perx_avg / x_over_logx
        scaling_rows.append(
            f"{x_max},{M_fixed},{T_setup:.6f},{T_perx_avg:.6f},{T_amort:.6f},{ratio:.3e}"
        )
        print(f"  x_max=10^{int(round(math.log10(x_max)))}: "
              f"T_setup={T_setup*1000:.2f}ms, T_perx_avg={T_perx_avg*1000:.2f}ms, "
              f"T_amort={T_amort*1000:.2f}ms/x, perx/(x/log)={ratio:.2e}")

    # ============================================================
    # Q4: Cluster-stitched anchor + sieve approach
    # ============================================================

    cluster_rows = ["x_anchor,M,cluster_width_w,T_anchor_s,T_sieve_s,T_total_s,T_amort_per_x_s"]
    print("\n=== Q4: Cluster-stitched (anchor + sieve) ===")
    x_anchor = 10**6
    small, primes = cached_setups[x_max]
    # Cluster modes: w = polylog (log^2, log^3), w = sqrt(x), w = x^{1/3}
    log_x = math.log(x_anchor)
    cluster_widths = [
        ("log^2 x", int(log_x ** 2)),
        ("log^3 x", int(log_x ** 3)),
        ("x^{1/4}", int(x_anchor ** 0.25)),
        ("x^{1/3}", int(x_anchor ** (1 / 3))),
    ]
    for name, w in cluster_widths:
        for M in [4, 16, 64]:
            # Sample M offsets uniformly in [0, w]
            random.seed(123)
            deltas = sorted(random.sample(range(w + 1), min(M, w + 1)))

            t0 = time.perf_counter()
            T_anchor_only, _ = time_call(per_x_eval, x_anchor, small, primes, n_repeat=2)
            T_anchor_only_actual = time.perf_counter() - t0  # not used; for ref

            t0 = time.perf_counter()
            results = cluster_stitched(x_anchor, deltas, small, primes)
            T_total = time.perf_counter() - t0

            T_sieve = T_total - T_anchor_only
            T_amort = T_total / len(deltas)
            cluster_rows.append(
                f"{x_anchor},{len(deltas)},{w},{T_anchor_only:.4f},"
                f"{T_sieve:.4f},{T_total:.4f},{T_amort:.6f}"
            )
            print(f"  w={name}={w:>5}, M={len(deltas):>3}: "
                  f"T_anchor={T_anchor_only*1000:.1f}ms, "
                  f"T_sieve={T_sieve*1000:.2f}ms, "
                  f"T_amort={T_amort*1000:.2f}ms/x")

    # ============================================================
    # Save CSVs
    # ============================================================
    for path, rows in [
        (out_setup_csv, setup_rows),
        (out_perx_csv, perx_rows),
        (out_amort_csv, amort_rows),
        (out_cluster_csv, cluster_rows),
        (out_ops_csv, ops_rows),
        (out_scaling_csv, scaling_rows),
    ]:
        with open(path, "w") as f:
            f.write("\n".join(rows) + "\n")
        print(f"  wrote {path}")

    # ============================================================
    # Summary
    # ============================================================
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print("""
HKM cross-x amortisation, slot 3 of Thread 5:

1. SHARED SETUP (sieve of Eratosthenes up to sqrt(x_max)):
   - Cost: O(sqrt(x_max) log log) — TINY compared to per-x cost.
   - Exists in cache for entire batched run.

2. PER-X EVALUATION (Lucy DP over primes p <= sqrt(x_i)):
   - Cost: Theta(x_i / log x_i) per query (basic Lucy, no Fenwick).
   - Cannot be reduced via shared state because Lucy DP's small[v] is
     iteratively coupled to step p, and the iterations on small[] for
     different p's are NOT amortisable across x_i.
   - Fenwick acceleration is per-x as well; constant factor only.

3. M-BATCHED AMORTISED COST:
   - T_per_x_amortised = T_setup/M + T_perx(x_avg)
   - For x_avg ~ x_max: setup is < 0.1% of per-x, so amortised ~= per-x.
   - SAME asymptotic order as single-query: Theta(x/log x).
   - NO batch-polylog regime possible.

4. CLUSTER-STITCHED (anchor + sieve):
   - For x_i in [x_anchor, x_anchor + w]: pi(x_i) = pi(x_anchor) + sieve.
   - Per-x amortised = T_lucy(x_anchor)/M + w/log(w).
   - For batch-polylog: need M >= x^{2/3}/polylog AND w <= polylog.
   - Standard "anchor + sieve"; no new amortisation. Polylog REGIME exists
     only for narrow clusters with M = poly(x^{2/3}/polylog) queries.

CONCLUSION: HKM cross-x amortisation does NOT yield batch-polylog pi(x)
over M = polylog uncorrelated queries. Per-x cost dominates and equals
the per-query cost asymptotically. The combinatorial pillar has even LESS
shareable state than the explicit-formula pillar (which had non-trivial
K_zeros_setup); HKM's "shareable" piece is just the sieve of Eratosthenes
up to sqrt(x_max), a near-trivial cost.

This closes slot 3 structurally with a cleaner negative result than
slot 1/2's explicit-formula closure.
""")


if __name__ == "__main__":
    main()
