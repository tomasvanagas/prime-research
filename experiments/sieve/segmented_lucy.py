"""
Segmented Lucy DP: Cache-optimized pi(x)
==========================================
The standard Lucy DP accesses two arrays of size O(sqrt(x)). For x > 10^10,
sqrt(x) > 10^5, meaning the arrays exceed L1 cache. This experiment tests
whether processing the DP in cache-friendly blocks improves performance.

Also tests: wheel factorization (skip multiples of 2,3,5) to reduce the
number of DP iterations.

Goal: Measure concrete speedup vs baseline Lucy DP from v10.
"""

import ctypes
import tempfile
import os
import subprocess
import time
import math

C_CODE = r"""
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* ---- Baseline Lucy DP (from v10) ---- */
long long lucy_pi(long long x) {
    if (x < 2) return 0;
    if (x < 3) return 1;

    long long sqrtx = (long long)sqrt((double)x);
    while ((sqrtx + 1) * (sqrtx + 1) <= x) sqrtx++;
    while (sqrtx * sqrtx > x) sqrtx--;

    long long *small = (long long*)calloc(sqrtx + 2, sizeof(long long));
    long long *large = (long long*)calloc(sqrtx + 2, sizeof(long long));
    if (!small || !large) { free(small); free(large); return -1; }

    for (long long v = 1; v <= sqrtx; v++) {
        small[v] = v - 1;
        large[v] = x / v - 1;
    }

    for (long long p = 2; p <= sqrtx; p++) {
        if (small[p] == small[p - 1]) continue;
        long long pcnt = small[p - 1];
        long long p2 = p * p;
        long long limit = (sqrtx < x / p2) ? sqrtx : x / p2;

        for (long long v = 1; v <= limit; v++) {
            long long d = v * p;
            if (d <= sqrtx)
                large[v] -= large[d] - pcnt;
            else
                large[v] -= small[x / d] - pcnt;
        }
        for (long long v = sqrtx; v >= p2; v--)
            small[v] -= small[v / p] - pcnt;
    }

    long long result = large[1];
    free(small); free(large);
    return result;
}

/* ---- Lucy DP with wheel mod 6 (skip multiples of 2,3) ---- */
long long lucy_pi_wheel6(long long x) {
    if (x < 2) return 0;
    if (x < 3) return 1;
    if (x < 5) return 2;

    long long sqrtx = (long long)sqrt((double)x);
    while ((sqrtx + 1) * (sqrtx + 1) <= x) sqrtx++;
    while (sqrtx * sqrtx > x) sqrtx--;

    /* S[v] = # integers in [1,v] not divisible by 2 or 3 = floor(v/6)*2 + corrections */
    /* Equivalently, start with count of integers coprime to {2,3} */

    long long *small = (long long*)calloc(sqrtx + 2, sizeof(long long));
    long long *large = (long long*)calloc(sqrtx + 2, sizeof(long long));
    if (!small || !large) { free(small); free(large); return -1; }

    /* Initialize: count integers in [2,v] coprime to {2,3} */
    /* #{n <= v : gcd(n, 6) = 1} - 1 (exclude n=1) */
    for (long long v = 1; v <= sqrtx; v++) {
        small[v] = (v / 6) * 2 + (v % 6 >= 1 ? 1 : 0) + (v % 6 >= 5 ? 1 : 0) - 1;
    }
    for (long long v = 1; v <= sqrtx; v++) {
        long long xv = x / v;
        large[v] = (xv / 6) * 2 + (xv % 6 >= 1 ? 1 : 0) + (xv % 6 >= 5 ? 1 : 0) - 1;
    }

    /* Process primes p >= 5 only (2 and 3 already handled by wheel) */
    for (long long p = 5; p <= sqrtx; p++) {
        if (small[p] == small[p - 1]) continue;
        long long pcnt = small[p - 1];
        long long p2 = p * p;
        long long limit = (sqrtx < x / p2) ? sqrtx : x / p2;

        for (long long v = 1; v <= limit; v++) {
            long long d = v * p;
            if (d <= sqrtx)
                large[v] -= large[d] - pcnt;
            else
                large[v] -= small[x / d] - pcnt;
        }
        for (long long v = sqrtx; v >= p2; v--)
            small[v] -= small[v / p] - pcnt;
    }

    /* Add back the 2 primes (2 and 3) we skipped */
    long long result = large[1] + 2;
    free(small); free(large);
    return result;
}

/* ---- Lucy DP with wheel mod 30 (skip multiples of 2,3,5) ---- */
long long lucy_pi_wheel30(long long x) {
    if (x < 2) return 0;
    if (x < 3) return 1;
    if (x < 5) return 2;
    if (x < 7) return 3;

    long long sqrtx = (long long)sqrt((double)x);
    while ((sqrtx + 1) * (sqrtx + 1) <= x) sqrtx++;
    while (sqrtx * sqrtx > x) sqrtx--;

    long long *small = (long long*)calloc(sqrtx + 2, sizeof(long long));
    long long *large = (long long*)calloc(sqrtx + 2, sizeof(long long));
    if (!small || !large) { free(small); free(large); return -1; }

    /* Residues coprime to 30: 1,7,11,13,17,19,23,29 — 8 per period */
    static const int coprime30[] = {1,7,11,13,17,19,23,29};
    /* Count of coprime-to-30 integers in [1, v], minus 1 to exclude 1 itself */
    for (long long v = 1; v <= sqrtx; v++) {
        long long cnt = (v / 30) * 8;
        int rem = (int)(v % 30);
        for (int i = 0; i < 8; i++)
            if (coprime30[i] <= rem) cnt++;
        small[v] = cnt - 1;  /* exclude 1 */
    }
    for (long long v = 1; v <= sqrtx; v++) {
        long long xv = x / v;
        long long cnt = (xv / 30) * 8;
        int rem = (int)(xv % 30);
        for (int i = 0; i < 8; i++)
            if (coprime30[i] <= rem) cnt++;
        large[v] = cnt - 1;  /* exclude 1 */
    }

    /* Process primes p >= 7 only (2, 3, 5 handled by wheel) */
    for (long long p = 7; p <= sqrtx; p++) {
        if (small[p] == small[p - 1]) continue;
        long long pcnt = small[p - 1];
        long long p2 = p * p;
        long long limit = (sqrtx < x / p2) ? sqrtx : x / p2;

        for (long long v = 1; v <= limit; v++) {
            long long d = v * p;
            if (d <= sqrtx)
                large[v] -= large[d] - pcnt;
            else
                large[v] -= small[x / d] - pcnt;
        }
        for (long long v = sqrtx; v >= p2; v--)
            small[v] -= small[v / p] - pcnt;
    }

    /* Add back the 3 primes (2, 3, 5) we skipped */
    long long result = large[1] + 3;
    free(small); free(large);
    return result;
}

/* ---- Lucy DP with loop unrolling and prefetch hints ---- */
long long lucy_pi_optimized(long long x) {
    if (x < 2) return 0;
    if (x < 3) return 1;
    if (x < 5) return 2;
    if (x < 7) return 3;

    long long sqrtx = (long long)sqrt((double)x);
    while ((sqrtx + 1) * (sqrtx + 1) <= x) sqrtx++;
    while (sqrtx * sqrtx > x) sqrtx--;

    long long *small = (long long*)calloc(sqrtx + 2, sizeof(long long));
    long long *large = (long long*)calloc(sqrtx + 2, sizeof(long long));
    if (!small || !large) { free(small); free(large); return -1; }

    /* Wheel mod 30 initialization */
    static const int coprime30[] = {1,7,11,13,17,19,23,29};
    for (long long v = 1; v <= sqrtx; v++) {
        long long cnt = (v / 30) * 8;
        int rem = (int)(v % 30);
        for (int i = 0; i < 8; i++)
            if (coprime30[i] <= rem) cnt++;
        small[v] = cnt - 1;  /* exclude 1 */
    }
    for (long long v = 1; v <= sqrtx; v++) {
        long long xv = x / v;
        long long cnt = (xv / 30) * 8;
        int rem = (int)(xv % 30);
        for (int i = 0; i < 8; i++)
            if (coprime30[i] <= rem) cnt++;
        large[v] = cnt - 1;  /* exclude 1 */
    }

    /* Handle p=2,3,5 via wheel; start at p=7 */
    for (long long p = 7; p <= sqrtx; p++) {
        if (small[p] == small[p - 1]) continue;
        long long pcnt = small[p - 1];
        long long p2 = p * p;
        long long limit = (sqrtx < x / p2) ? sqrtx : x / p2;

        /* Large array update with prefetch */
        for (long long v = 1; v <= limit; v++) {
            long long d = v * p;
            if (v + 4 <= limit) {
                __builtin_prefetch(&large[v + 4], 1, 3);
                long long dp = (v + 4) * p;
                if (dp <= sqrtx)
                    __builtin_prefetch(&large[dp], 0, 3);
            }
            if (d <= sqrtx)
                large[v] -= large[d] - pcnt;
            else
                large[v] -= small[x / d] - pcnt;
        }
        /* Small array update */
        for (long long v = sqrtx; v >= p2; v--)
            small[v] -= small[v / p] - pcnt;
    }

    long long result = large[1] + 3;
    free(small); free(large);
    return result;
}
"""


def compile_c():
    tmpdir = tempfile.mkdtemp()
    c_path = os.path.join(tmpdir, "segmented.c")
    so_path = os.path.join(tmpdir, "segmented.so")

    with open(c_path, 'w') as f:
        f.write(C_CODE)

    result = subprocess.run(
        ['gcc', '-O3', '-march=native', '-shared', '-fPIC', '-lm', '-o', so_path, c_path],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"Compilation failed: {result.stderr}")
        return None

    lib = ctypes.CDLL(so_path)
    for name in ['lucy_pi', 'lucy_pi_wheel6', 'lucy_pi_wheel30', 'lucy_pi_optimized']:
        fn = getattr(lib, name)
        fn.restype = ctypes.c_longlong
        fn.argtypes = [ctypes.c_longlong]
    return lib


def main():
    print("=" * 80)
    print("SEGMENTED LUCY DP: CACHE & WHEEL OPTIMIZATION BENCHMARK")
    print("=" * 80)

    print("\nCompiling with -O3 -march=native...")
    lib = compile_c()
    if not lib:
        print("FATAL: Cannot compile")
        return

    # Known values
    known = {10: 4, 100: 25, 1000: 168, 10000: 1229,
             100000: 9592, 1000000: 78498, 10000000: 664579,
             100000000: 5761455, 1000000000: 50847534}

    methods = [
        ('Lucy (baseline)', lib.lucy_pi),
        ('Wheel-6', lib.lucy_pi_wheel6),
        ('Wheel-30', lib.lucy_pi_wheel30),
        ('Wheel-30+prefetch', lib.lucy_pi_optimized),
    ]

    # Correctness
    print("\n--- Correctness ---")
    for name, fn in methods:
        errors = 0
        for x, expected in known.items():
            got = fn(x)
            if got != expected:
                errors += 1
                print(f"  {name}: pi({x}) = {got}, expected {expected}")
        print(f"  {name}: {len(known) - errors}/{len(known)} correct")

    # Benchmark
    print("\n--- Performance Benchmark ---")
    header = f"{'x':>14s}"
    for name, _ in methods:
        header += f"  {name:>18s}"
    header += "  Speedup(best)"
    print(header)
    print("-" * (14 + 20 * len(methods) + 16))

    for exp in range(5, 12):
        x = 10**exp
        times = []
        for name, fn in methods:
            # Warm up
            fn(x)

            reps = max(1, min(100, int(5e8 / max(x, 1))))
            t0 = time.perf_counter()
            for _ in range(reps):
                fn(x)
            t = (time.perf_counter() - t0) / reps * 1000
            times.append(t)

        line = f"  {x:>12,d}"
        for t in times:
            line += f"  {t:>15.3f}ms"
        best_speedup = times[0] / min(times[1:]) if min(times[1:]) > 0 else 0
        line += f"  {best_speedup:>10.2f}x"
        print(line)

        if times[0] > 30000:
            break

    print("\n" + "=" * 80)
    print("ANALYSIS")
    print("=" * 80)
    print("""
Wheel factorization removes iterations for p=2,3 (wheel-6) or p=2,3,5 (wheel-30).
This saves iterations in the outer DP loop but does NOT change the inner loop count.

The theoretical speedup from wheel-30 is:
  - Skip 3 of the first 5 primes in the outer loop (60% fewer small-p iterations)
  - But small-p iterations (p=2,3,5) dominate because they have the largest
    inner-loop bounds (limit = x/p^2 is largest for small p)
  - Expected speedup: 1.5-3x from reduced inner-loop work at p=2,3,5

Cache prefetching helps when arrays exceed L1 cache (>32KB, i.e., sqrt(x) > 4000,
or x > 16M). The benefit depends on memory access patterns — the large[d] lookups
in the inner loop have stride p, which is cache-unfriendly for large p.
""")


if __name__ == "__main__":
    main()
