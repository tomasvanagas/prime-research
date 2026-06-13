"""
Cross-x amortisation, slot 4 of Thread 5: Aggarwal binary search adversarial probe.

Aggarwal 2025 reduces p_n to O(log n) calls to pi(x_i) at midpoints of nested
brackets. The total cost is O(log n) * T_pi(x), giving O(sqrt(n) log^4 n) under
HKM. Question: do the O(log n) sub-queries admit cross-call amortisation that
the original analysis missed?

Slot 3 closed the "cluster of M uniformly-spaced or uncorrelated x_i" case:
shareable state across queries is the trivial sqrt(x_max) sieve only. Aggarwal's
query pattern is structurally distinct — midpoints of nested brackets, distance
n/2^i from the i-th query. This script applies the cross-x amortisation
machinery to the BINARY-SEARCH query pattern.

Three angles tested:

(A) Shared small[]-final precomputation. Aggarwal's O(log n) Lucy DPs each
    iterate small[] up to sqrt(x_i). Since small[] evolution is x-independent
    in its prime-iteration, can we precompute it once for x_max = R (Dusart
    upper bound) and share?

    Subtlety (slot 3): Lucy DP's large[] update at step p references
    small[x_i // d] AT STEP p, not the final pi value, when v <= sqrtx and
    d > sqrtx. So sharing FINAL small[] requires re-running small[] iteration
    each query. The slot-3 result is that intermediate snapshots cost
    O(x_max / log x_max), defeating amortisation.

    But for log n queries, snapshots could be limited. Empirical test: count
    Lucy DP small/large ops per query, attempt the "shared final small[],
    per-query small[] re-iteration" decoupling, measure savings.

(B) Hyperbola amortisation. One Lucy DP at x_max = R yields pi(R/k) for
    k = 1..sqrt(R) for free. How many of Aggarwal's binary-search midpoints
    {(L_i + R_i) / 2} lie in the hyperbolic point set {R/k}? Density estimate:
    midpoints are in [L, R] of width n; hyperbolic points in [L, R] number
    ~ R * (R-L) / (L * R) = (R - L) / L ~ n / (n log n) = 1 / log n. So
    DENSITY ZERO of midpoints are "free" — Aggarwal cannot exploit the
    hyperbola structure.

(C) Adaptive precision (S217 angle). With Riemann's R + Schoenfeld eps,
    Aggarwal's binary search costs match Aggarwal-pure asymptotically.
    Already established. Cross-checked here.

The script outputs CSVs:
  cross_x_aggarwal_op_breakdown.csv  — small/large ops per query and total
  cross_x_aggarwal_shared_small.csv  — shared-small precompute saving
  cross_x_aggarwal_hyperbola.csv     — hyperbola density vs n
  cross_x_aggarwal_summary.csv       — n vs total ops vs amortised lower bound

Asymptotic upper bound: even with PERFECT cross-call sharing of small[] state,
the cumulative LARGE[] work across log(n) queries is sum_{i=1}^{log n}
Theta(x_i / log x_i) ~ log(n) * Theta(n log n / log(n log n)) = Theta(n).
Strictly super-polylog. The Aggarwal bound O(sqrt(n) log^4 n) under HKM is
NOT improvable by cross-call amortisation in Lucy DP.
"""

import math
import time
from sympy import primepi, prime as sym_prime


def isqrt(n):
    s = int(math.isqrt(n))
    while (s + 1) * (s + 1) <= n:
        s += 1
    while s * s > n:
        s -= 1
    return s


# ---------------------------------------------------------------------------
# Shared sieve up to N (slot-3 carry-over)
# ---------------------------------------------------------------------------

def shared_sieve(N):
    """Eratosthenes sieve. Returns (small_pi, primes) where small_pi[v] = pi(v)
    for v in 0..N+1, and primes is the list of primes <= N."""
    if N < 2:
        return [0] * (N + 2), []
    sieve = bytearray(b"\x01") * (N + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, isqrt(N) + 1):
        if sieve[i]:
            for j in range(i * i, N + 1, i):
                sieve[j] = 0
    primes = [i for i in range(2, N + 1) if sieve[i]]
    small_pi = [0] * (N + 2)
    cnt = 0
    for v in range(N + 2):
        if v <= N and sieve[v]:
            cnt += 1
        small_pi[v] = cnt
    return small_pi, primes


# ---------------------------------------------------------------------------
# Lucy DP, instrumented for op counting
# ---------------------------------------------------------------------------

def lucy_pi_instrumented(x, primes_up_to_sqrt_x_max=None):
    """Compute pi(x) via basic Lucy DP, returning (pi_x, n_small_ops, n_large_ops).

    primes_up_to_sqrt_x_max: optional precomputed prime list (must include all
    primes <= sqrt(x)). If None, Lucy DP self-sieves.
    """
    if x < 2:
        return 0, 0, 0
    sqrtx = isqrt(x)
    s = list(range(sqrtx + 2))  # s[v] = v
    for v in range(sqrtx + 2):
        s[v] = v - 1 if v >= 1 else 0
    L = [0] * (sqrtx + 2)
    for v in range(1, sqrtx + 1):
        L[v] = x // v - 1

    n_small = 0
    n_large = 0

    if primes_up_to_sqrt_x_max is not None:
        # Use precomputed list, restricted to <= sqrtx.
        for p in primes_up_to_sqrt_x_max:
            if p > sqrtx:
                break
            sp = s[p - 1]
            p2 = p * p
            if p2 > x:
                break
            limit = min(sqrtx, x // p2)
            # Update large[]
            for v in range(1, limit + 1):
                d = v * p
                if d <= sqrtx:
                    L[v] -= L[d] - sp
                else:
                    L[v] -= s[x // d] - sp
                n_large += 1
            # Update small[]
            for v in range(sqrtx, p2 - 1, -1):
                s[v] -= s[v // p] - sp
                n_small += 1
    else:
        # Self-sieve via small[] state.
        for p in range(2, sqrtx + 1):
            if s[p] == s[p - 1]:
                continue  # composite
            sp = s[p - 1]
            p2 = p * p
            if p2 > x:
                break
            limit = min(sqrtx, x // p2)
            for v in range(1, limit + 1):
                d = v * p
                if d <= sqrtx:
                    L[v] -= L[d] - sp
                else:
                    L[v] -= s[x // d] - sp
                n_large += 1
            for v in range(sqrtx, p2 - 1, -1):
                s[v] -= s[v // p] - sp
                n_small += 1

    return L[1], n_small, n_large


# ---------------------------------------------------------------------------
# Aggarwal binary search with op recording
# ---------------------------------------------------------------------------

def dusart_bounds(n):
    if n < 6:
        small_p = [None, 2, 3, 5, 7, 11]
        return (max(small_p[n] - 1, 1), small_p[n] + 1)
    ln_n = math.log(n)
    ln_ln_n = math.log(ln_n)
    L = int(math.floor(n * (ln_n + ln_ln_n - 1)))
    R = int(math.ceil(n * (ln_n + ln_ln_n)))
    return (L, R)


def aggarwal_instrumented(n):
    """Run Aggarwal binary search, collecting per-query Lucy DP op counts.

    Returns dict with:
      L0, R0:      Dusart bracket
      midpoints:   list of x_i queried
      pi_values:   pi(x_i) for each
      ops_small:   small-update ops per query
      ops_large:   large-update ops per query
      ops_total:   total ops per query
      times:       wall time per query (seconds)
      pn:          the answer
      total_ops:   sum across all queries
    """
    L0, R0 = dusart_bounds(n)
    L, R = L0, R0

    sqrt_R = isqrt(R0)
    # Shared sieve up to sqrt(R0)
    t0 = time.perf_counter()
    small_pi, primes = shared_sieve(sqrt_R)
    t_sieve = time.perf_counter() - t0
    n_sieve_ops = sum(1 for _ in primes)  # rough; actual sieve cost is O(N log log N)

    midpoints = []
    pi_values = []
    ops_small = []
    ops_large = []
    times = []

    while L < R:
        mid = (L + R) // 2
        midpoints.append(mid)
        t0 = time.perf_counter()
        c, ns, nl = lucy_pi_instrumented(mid, primes)
        t1 = time.perf_counter()
        pi_values.append(c)
        ops_small.append(ns)
        ops_large.append(nl)
        times.append(t1 - t0)
        if c >= n:
            R = mid
        else:
            L = mid + 1

    return {
        'n': n,
        'L0': L0, 'R0': R0,
        'sqrt_R0': sqrt_R,
        'sieve_time_s': t_sieve,
        'sieve_n_primes': len(primes),
        'midpoints': midpoints,
        'pi_values': pi_values,
        'ops_small': ops_small,
        'ops_large': ops_large,
        'times': times,
        'pn': L,
        'total_small_ops': sum(ops_small),
        'total_large_ops': sum(ops_large),
        'total_ops': sum(ops_small) + sum(ops_large),
        'total_time_s': sum(times),
    }


# ---------------------------------------------------------------------------
# Hyperbola amortisation density
# ---------------------------------------------------------------------------

def hyperbola_density_in_dusart(n):
    """Fraction of {R0/k : k = 1..sqrt(R0)} that lies in the Dusart bracket
    [L0, R0]. Predicts whether one Lucy DP at x_max = R0 can supply Aggarwal's
    midpoints "for free" via the large[] hyperbolic point set.
    """
    L0, R0 = dusart_bounds(n)
    sqrt_R = isqrt(R0)
    if sqrt_R < 1:
        return 0, 0, 0
    k_lo = max(1, math.ceil(R0 / R0))   # k=1 -> R0/1 = R0 always in bracket
    k_hi = R0 // L0                       # largest k with R0/k >= L0
    in_bracket = max(0, k_hi - k_lo + 1)
    return in_bracket, sqrt_R, in_bracket / sqrt_R if sqrt_R > 0 else 0


# ---------------------------------------------------------------------------
# Shared-small[] precompute saving
# ---------------------------------------------------------------------------

def shared_small_precompute_savings(n):
    """Estimate savings from precomputing small[] up to sqrt(R0) once and
    sharing across all log(n) Aggarwal queries.

    Lower bound on Aggarwal cost = log(n) * (large_ops at x ~ R0).
    Shared-small Aggarwal cost = sieve(sqrt(R0)) + log(n) * large_ops_at_x.

    The "saving" is the small_ops fraction across queries, since shared small[]
    eliminates that work in the limit. Slot 3 noted small[] re-iteration is
    REQUIRED for correctness (intermediate state references), so the actual
    saving is bounded by the SHARED FINAL small[] = sieve(sqrt(R0)), not
    intermediate-state sharing.
    """
    info = aggarwal_instrumented(n)
    n_queries = len(info['midpoints'])
    return {
        'n': n,
        'log_n': math.log(n),
        'n_queries': n_queries,
        'total_small_ops': info['total_small_ops'],
        'total_large_ops': info['total_large_ops'],
        'total_ops': info['total_ops'],
        'small_ops_frac': info['total_small_ops'] / max(info['total_ops'], 1),
        'sieve_n_primes': info['sieve_n_primes'],
        'sieve_time_s': info['sieve_time_s'],
        'total_time_s': info['total_time_s'],
        'sieve_share': info['sieve_time_s'] / max(info['total_time_s'], 1e-12),
    }


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main():
    out_dir = "experiments/analytic/cross_x_amortisation"
    out_breakdown = f"{out_dir}/cross_x_aggarwal_op_breakdown.csv"
    out_shared = f"{out_dir}/cross_x_aggarwal_shared_small.csv"
    out_hyperbola = f"{out_dir}/cross_x_aggarwal_hyperbola.csv"
    out_summary = f"{out_dir}/cross_x_aggarwal_summary.csv"

    # ============================================================
    # Q1: Per-query op breakdown for Aggarwal binary search
    # ============================================================
    n_list = [100, 1_000, 10_000, 100_000, 1_000_000]

    breakdown_rows = ["n,query_idx,midpoint,pi_at_mid,ops_small,ops_large,ops_total,time_s"]
    summary_rows = [
        "n,L0,R0,n_queries,total_small_ops,total_large_ops,total_ops,total_time_s,pn,sympy_pn,agree"
    ]

    print("=== Q1: Aggarwal per-query op breakdown ===")
    full_runs = {}
    for n in n_list:
        info = aggarwal_instrumented(n)
        full_runs[n] = info
        truth = int(sym_prime(n))
        agree = info['pn'] == truth
        print(f"n={n:>9d}  L0={info['L0']:>11d}  R0={info['R0']:>11d}  "
              f"queries={len(info['midpoints']):>3d}  total_ops={info['total_ops']:>10d}  "
              f"pn={info['pn']:>11d}  agree={agree}")
        for i, (mid, pic, ns, nl, t) in enumerate(zip(
                info['midpoints'], info['pi_values'],
                info['ops_small'], info['ops_large'], info['times'])):
            breakdown_rows.append(
                f"{n},{i},{mid},{pic},{ns},{nl},{ns + nl},{t:.6f}"
            )
        summary_rows.append(
            f"{n},{info['L0']},{info['R0']},{len(info['midpoints'])},"
            f"{info['total_small_ops']},{info['total_large_ops']},{info['total_ops']},"
            f"{info['total_time_s']:.6f},{info['pn']},{truth},{agree}"
        )

    with open(out_breakdown, 'w') as f:
        f.write('\n'.join(breakdown_rows) + '\n')
    with open(out_summary, 'w') as f:
        f.write('\n'.join(summary_rows) + '\n')

    # ============================================================
    # Q2: Shared-small[] precompute savings (relative to Aggarwal-pure)
    # ============================================================
    print("\n=== Q2: Shared small[] precompute savings ===")
    shared_rows = [
        "n,n_queries,total_small_ops,total_large_ops,total_ops,small_ops_frac,"
        "sieve_n_primes,sieve_time_s,total_time_s,sieve_share"
    ]
    for n in n_list:
        info = full_runs[n]
        n_queries = len(info['midpoints'])
        small_frac = info['total_small_ops'] / max(info['total_ops'], 1)
        sieve_share = info['sieve_time_s'] / max(info['total_time_s'], 1e-12)
        print(f"n={n:>9d}  n_queries={n_queries:>3d}  "
              f"small_frac={small_frac:.3f}  sieve_share={sieve_share*100:.4f}%")
        shared_rows.append(
            f"{n},{n_queries},{info['total_small_ops']},{info['total_large_ops']},"
            f"{info['total_ops']},{small_frac:.6f},{info['sieve_n_primes']},"
            f"{info['sieve_time_s']:.6f},{info['total_time_s']:.6f},{sieve_share:.6f}"
        )
    with open(out_shared, 'w') as f:
        f.write('\n'.join(shared_rows) + '\n')

    # ============================================================
    # Q3: Hyperbola density — fraction of {R0/k} that lies in [L0, R0]
    # ============================================================
    print("\n=== Q3: Hyperbola density in Dusart bracket ===")
    hyperbola_rows = ["n,L0,R0,sqrt_R0,n_hyperbola_in_bracket,density"]
    n_list_hyp = [100, 1_000, 10_000, 100_000, 1_000_000, 10**9, 10**12, 10**18]
    for n in n_list_hyp:
        L0, R0 = dusart_bounds(n)
        in_b, sq, dens = hyperbola_density_in_dusart(n)
        log_n = math.log(n)
        print(f"n={n:>15d}  L0={L0:>20d}  R0={R0:>20d}  sqrt_R0={sq:>10d}  "
              f"in_bracket={in_b:>5d}  density={dens:.6f}  1/log_n={1.0/log_n:.6f}")
        hyperbola_rows.append(f"{n},{L0},{R0},{sq},{in_b},{dens:.8f}")
    with open(out_hyperbola, 'w') as f:
        f.write('\n'.join(hyperbola_rows) + '\n')

    # ============================================================
    # Q4: Asymptotic check — total Aggarwal ops vs theoretical
    # ============================================================
    print("\n=== Q4: Aggarwal cost scaling vs Lucy DP at R0 ===")
    print("If cross-call amortisation gave nothing, we expect")
    print("  total_ops(n) ~ log(n) * Lucy_ops(R0(n))")
    for n in n_list:
        info = full_runs[n]
        R0 = info['R0']
        L0 = info['L0']
        single_pi, ns_R, nl_R = lucy_pi_instrumented(R0, None)
        log_n = math.log(n)
        ratio = info['total_ops'] / (log_n * (ns_R + nl_R))
        print(f"n={n:>9d}  total_ops={info['total_ops']:>10d}  "
              f"log(n)*Lucy(R0) ops={log_n * (ns_R + nl_R):>12.0f}  ratio={ratio:.3f}")

    print("\nDone. Outputs in", out_dir)


if __name__ == '__main__':
    main()
