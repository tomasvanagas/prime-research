"""
Lehmer-method pi(x) vs Lucy DP benchmark
==========================================
Implements the Meissel-Lehmer formula for pi(x) in C and benchmarks
against the naive Lucy DP from v10. The Lehmer method decomposes:

  pi(x) = phi(x, a) + a - 1 - P2(x, a)

where a = pi(x^{1/3}), phi is the partial sieve function, and P2
counts semiprimes. This decomposition enables faster computation
by separating "easy" and "hard" parts.

Goal: Measure the constant-factor speedup from this decomposition.
"""

import ctypes
import tempfile
import os
import subprocess
import time
import math

# ============================================================
# C implementation: both Lucy DP and Lehmer method
# ============================================================

C_CODE = r"""
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* ---- Sieve of Eratosthenes ---- */
static int *sieve_primes(int limit, int *count) {
    char *is_p = (char*)calloc(limit + 1, 1);
    if (!is_p) return NULL;
    memset(is_p, 1, limit + 1);
    is_p[0] = is_p[1] = 0;
    for (int i = 2; (long long)i * i <= limit; i++)
        if (is_p[i])
            for (int j = i * i; j <= limit; j += i)
                is_p[j] = 0;
    *count = 0;
    for (int i = 2; i <= limit; i++)
        if (is_p[i]) (*count)++;
    int *primes = (int*)malloc((*count) * sizeof(int));
    int idx = 0;
    for (int i = 2; i <= limit; i++)
        if (is_p[i]) primes[idx++] = i;
    free(is_p);
    return primes;
}

/* ---- Lucy DP (baseline from v10) ---- */
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

/* ---- Lehmer method ---- */

/* pi_small: count primes <= x using a sieve, for x up to ~10^7 */
static long long *pi_table = NULL;
static int pi_table_size = 0;

static void build_pi_table(int limit) {
    if (pi_table && pi_table_size >= limit) return;
    if (pi_table) free(pi_table);
    pi_table_size = limit;
    pi_table = (long long*)calloc(limit + 1, sizeof(long long));
    char *is_p = (char*)calloc(limit + 1, 1);
    memset(is_p, 1, limit + 1);
    is_p[0] = is_p[1] = 0;
    for (int i = 2; (long long)i * i <= limit; i++)
        if (is_p[i])
            for (int j = i * i; j <= limit; j += i)
                is_p[j] = 0;
    pi_table[0] = 0;
    for (int i = 1; i <= limit; i++)
        pi_table[i] = pi_table[i-1] + is_p[i];
    free(is_p);
}

static long long pi_small(long long x) {
    if (x <= 0) return 0;
    if (x <= pi_table_size) return pi_table[x];
    /* Fallback to Lucy for values beyond table */
    return lucy_pi(x);
}

/*
 * phi(x, a) = count of integers <= x not divisible by any of first a primes.
 * Computed via recursive inclusion-exclusion with memoization.
 */
#define PHI_CACHE_SIZE 16384
typedef struct { long long x; int a; long long val; char used; } PhiEntry;
static PhiEntry phi_cache[PHI_CACHE_SIZE];

static void phi_cache_clear(void) {
    memset(phi_cache, 0, sizeof(phi_cache));
}

static long long phi_func(long long x, int a, int *primes) {
    if (a == 0) return x;
    if (a == 1) return x - x / 2;

    /* Check cache */
    unsigned int h = ((unsigned long long)x * 2654435761ULL + a * 40503ULL) % PHI_CACHE_SIZE;
    if (phi_cache[h].used && phi_cache[h].x == x && phi_cache[h].a == a)
        return phi_cache[h].val;

    long long result = phi_func(x, a - 1, primes) - phi_func(x / primes[a - 1], a - 1, primes);

    /* Store in cache */
    phi_cache[h].x = x;
    phi_cache[h].a = a;
    phi_cache[h].val = result;
    phi_cache[h].used = 1;

    return result;
}

/*
 * Lehmer's formula:
 *   pi(x) = phi(x, a) + a - 1 - P2(x, a)
 * where a = pi(x^{1/3})
 * P2(x, a) = sum_{a < i <= pi(sqrt(x))} [pi(x/p_i) - i + 1]
 */
long long lehmer_pi(long long x) {
    if (x < 2) return 0;
    if (x < 3) return 1;
    if (x < 5) return 2;

    long long cbrtx = (long long)cbrt((double)x);
    while ((cbrtx + 1) * (cbrtx + 1) * (cbrtx + 1) <= x) cbrtx++;
    while (cbrtx * cbrtx * cbrtx > x) cbrtx--;

    long long sqrtx = (long long)sqrt((double)x);
    while ((sqrtx + 1) * (sqrtx + 1) <= x) sqrtx++;
    while (sqrtx * sqrtx > x) sqrtx--;

    /* Build pi table up to sqrt(x) */
    int table_lim = (int)(sqrtx < 100000000LL ? sqrtx : 100000000LL);
    build_pi_table(table_lim);

    /* Sieve primes up to sqrt(x) */
    int nprimes;
    int *primes = sieve_primes((int)sqrtx, &nprimes);
    if (!primes) return -1;

    int a = (int)pi_small(cbrtx);  /* a = pi(x^{1/3}) */
    int b = (int)pi_small(sqrtx);  /* b = pi(x^{1/2}) */

    /* phi(x, a) */
    phi_cache_clear();
    long long phi_val = phi_func(x, a, primes);

    /* P2(x, a) = sum_{a < i <= b} [pi(x/p_i) - i + 1] */
    long long P2 = 0;
    for (int i = a; i < b; i++) {
        long long w = x / primes[i];
        long long pi_w;
        if (w <= table_lim)
            pi_w = pi_table[w];
        else
            pi_w = lucy_pi(w);  /* Fall back to Lucy for large values */
        P2 += pi_w - (i + 1) + 1;
    }

    long long result = phi_val + a - 1 - P2;
    free(primes);
    return result;
}

/* ---- Extended Meissel-Lehmer with P3 term ---- */
long long meissel_lehmer_pi(long long x) {
    if (x < 2) return 0;
    if (x < 3) return 1;
    if (x < 5) return 2;

    /* Use x^{1/4} as the sieve limit for better decomposition */
    long long x14 = (long long)pow((double)x, 0.25);
    while (x14 > 1 && x14 * x14 * x14 * x14 > x) x14--;

    long long cbrtx = (long long)cbrt((double)x);
    while ((cbrtx + 1) * (cbrtx + 1) * (cbrtx + 1) <= x) cbrtx++;
    while (cbrtx * cbrtx * cbrtx > x) cbrtx--;

    long long sqrtx = (long long)sqrt((double)x);
    while ((sqrtx + 1) * (sqrtx + 1) <= x) sqrtx++;
    while (sqrtx * sqrtx > x) sqrtx--;

    int table_lim = (int)(sqrtx < 100000000LL ? sqrtx : 100000000LL);
    build_pi_table(table_lim);

    int nprimes;
    int *primes = sieve_primes((int)sqrtx, &nprimes);
    if (!primes) return -1;

    int a = (int)pi_small(x14);    /* pi(x^{1/4}) */
    int b = (int)pi_small(sqrtx);  /* pi(x^{1/2}) */
    int c = (int)pi_small(cbrtx);  /* pi(x^{1/3}) */

    /* phi(x, a) with a = pi(x^{1/4}) — fewer recursions */
    phi_cache_clear();
    long long phi_val = phi_func(x, a, primes);

    /* P2(x, a) = sum_{a < i <= b} [pi(x/p_i) - i + 1] */
    long long P2 = 0;
    for (int i = a; i < b; i++) {
        long long w = x / primes[i];
        long long pi_w;
        if (w <= table_lim)
            pi_w = pi_table[w];
        else
            pi_w = lucy_pi(w);
        P2 += pi_w - (long long)(i + 1) + 1;
    }

    /* P3(x, a) = sum_{a < i <= c} sum_{i <= j <= pi(sqrt(x/p_i))}
     *            [pi(x/(p_i*p_j)) - j + 1]
     * Counts triprimes with smallest factor > p_a
     */
    long long P3 = 0;
    for (int i = a; i < c; i++) {
        long long xi = x / primes[i];
        long long sqrtxi = (long long)sqrt((double)xi);
        int bi = (int)pi_small(sqrtxi);

        for (int j = i; j < bi; j++) {
            long long w = xi / primes[j];
            long long pi_w;
            if (w <= table_lim)
                pi_w = pi_table[w];
            else
                pi_w = lucy_pi(w);
            P3 += pi_w - (long long)(j + 1) + 1;
        }
    }

    long long result = phi_val + a - 1 - P2 - P3;
    free(primes);
    return result;
}

/* ---- Cleanup ---- */
void cleanup(void) {
    if (pi_table) { free(pi_table); pi_table = NULL; pi_table_size = 0; }
}
"""


def compile_c():
    """Compile the C code to a shared library."""
    tmpdir = tempfile.mkdtemp()
    c_path = os.path.join(tmpdir, "pi_methods.c")
    so_path = os.path.join(tmpdir, "pi_methods.so")

    with open(c_path, 'w') as f:
        f.write(C_CODE)

    result = subprocess.run(
        ['gcc', '-O3', '-shared', '-fPIC', '-lm', '-o', so_path, c_path],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"Compilation failed: {result.stderr}")
        return None

    lib = ctypes.CDLL(so_path)
    lib.lucy_pi.restype = ctypes.c_longlong
    lib.lucy_pi.argtypes = [ctypes.c_longlong]
    lib.lehmer_pi.restype = ctypes.c_longlong
    lib.lehmer_pi.argtypes = [ctypes.c_longlong]
    lib.meissel_lehmer_pi.restype = ctypes.c_longlong
    lib.meissel_lehmer_pi.argtypes = [ctypes.c_longlong]
    lib.cleanup.restype = None
    lib.cleanup.argtypes = []
    return lib


def sieve_count(limit):
    """Reference pi(x) via sieve."""
    if limit < 2:
        return 0
    s = bytearray(b'\x01') * (limit + 1)
    s[0] = s[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if s[i]:
            s[i*i::i] = bytearray(len(s[i*i::i]))
    return sum(s)


def main():
    print("=" * 70)
    print("LEHMER vs LUCY DP BENCHMARK")
    print("=" * 70)

    print("\nCompiling C code...")
    lib = compile_c()
    if not lib:
        print("FATAL: Cannot compile C code")
        return

    # ---- Correctness verification ----
    print("\n--- Correctness Verification ---")
    errors = {'lucy': 0, 'lehmer': 0, 'meissel': 0}
    test_values = [10, 100, 1000, 10000, 100000, 1000000]

    # Known values of pi(x)
    known_pi = {
        10: 4, 100: 25, 1000: 168, 10000: 1229,
        100000: 9592, 1000000: 78498
    }

    for x in test_values:
        pi_lucy = lib.lucy_pi(x)
        pi_lehmer = lib.lehmer_pi(x)
        pi_meissel = lib.meissel_lehmer_pi(x)
        expected = known_pi[x]

        status_l = "OK" if pi_lucy == expected else "FAIL"
        status_lh = "OK" if pi_lehmer == expected else "FAIL"
        status_ml = "OK" if pi_meissel == expected else "FAIL"

        if pi_lucy != expected:
            errors['lucy'] += 1
        if pi_lehmer != expected:
            errors['lehmer'] += 1
        if pi_meissel != expected:
            errors['meissel'] += 1

        print(f"  pi({x:>10,d}) = {expected:>8,d}  "
              f"Lucy={pi_lucy:>8,d}[{status_l}]  "
              f"Lehmer={pi_lehmer:>8,d}[{status_lh}]  "
              f"M-L={pi_meissel:>8,d}[{status_ml}]")

    if any(errors.values()):
        print(f"\n  ERRORS: {errors}")
        print("  Stopping — correctness failures must be investigated.")

    # Exhaustive check at many values
    print("\n  Exhaustive check (pi(x) for x in [2, 10000])...")
    exhaustive_errors = {'lucy': 0, 'lehmer': 0, 'meissel': 0}
    for x in range(2, 10001):
        expected = sieve_count(x)
        pl = lib.lucy_pi(x)
        plh = lib.lehmer_pi(x)
        pml = lib.meissel_lehmer_pi(x)
        if pl != expected:
            exhaustive_errors['lucy'] += 1
        if plh != expected:
            exhaustive_errors['lehmer'] += 1
        if pml != expected:
            exhaustive_errors['meissel'] += 1

    print(f"  Lucy:   {10000-1-exhaustive_errors['lucy']}/{10000-1} correct")
    print(f"  Lehmer: {10000-1-exhaustive_errors['lehmer']}/{10000-1} correct")
    print(f"  M-L:    {10000-1-exhaustive_errors['meissel']}/{10000-1} correct")

    if any(exhaustive_errors.values()):
        print(f"  ERRORS in exhaustive check: {exhaustive_errors}")

    # Spot-check larger values
    print("\n  Spot-checking larger values...")
    large_tests = [10**7, 5*10**7, 10**8, 5*10**8, 10**9]
    for x in large_tests:
        pl = lib.lucy_pi(x)
        plh = lib.lehmer_pi(x)
        pml = lib.meissel_lehmer_pi(x)
        agree = (pl == plh == pml)
        print(f"  pi({x:>12,d}): Lucy={pl:>12,d}  Lehmer={plh:>12,d}  M-L={pml:>12,d}  {'AGREE' if agree else 'DISAGREE!'}")

    # ---- Benchmark ----
    print("\n--- Performance Benchmark ---")
    print(f"{'x':>14s}  {'Lucy (ms)':>12s}  {'Lehmer (ms)':>12s}  {'M-L (ms)':>12s}  {'Speedup L':>10s}  {'Speedup ML':>10s}")
    print("-" * 76)

    benchmark_values = []
    for exp in range(5, 11):
        benchmark_values.append(10**exp)
        if exp < 10:
            benchmark_values.append(5 * 10**exp)

    for x in sorted(benchmark_values):
        # Warm up
        lib.lucy_pi(x)
        lib.lehmer_pi(x)
        lib.meissel_lehmer_pi(x)

        # Time Lucy
        reps = max(1, min(100, int(1e9 / max(x, 1))))
        t0 = time.perf_counter()
        for _ in range(reps):
            lib.lucy_pi(x)
        t_lucy = (time.perf_counter() - t0) / reps * 1000

        # Time Lehmer
        t0 = time.perf_counter()
        for _ in range(reps):
            lib.lehmer_pi(x)
        t_lehmer = (time.perf_counter() - t0) / reps * 1000

        # Time Meissel-Lehmer
        t0 = time.perf_counter()
        for _ in range(reps):
            lib.meissel_lehmer_pi(x)
        t_meissel = (time.perf_counter() - t0) / reps * 1000

        speedup_l = t_lucy / t_lehmer if t_lehmer > 0 else float('inf')
        speedup_ml = t_lucy / t_meissel if t_meissel > 0 else float('inf')

        print(f"  {x:>12,d}  {t_lucy:>12.3f}  {t_lehmer:>12.3f}  {t_meissel:>12.3f}  {speedup_l:>10.2f}x  {speedup_ml:>10.2f}x")

        if t_lucy > 30000:  # 30 seconds
            print("  (stopping — too slow)")
            break

    lib.cleanup()

    print("\n" + "=" * 70)
    print("ANALYSIS")
    print("=" * 70)
    print("""
The Lehmer and Meissel-Lehmer methods decompose pi(x) computation into:
  1. phi(x, a) — partial sieve via inclusion-exclusion
  2. P2 — semiprime counting
  3. P3 — triprime counting (Meissel-Lehmer only)

vs Lucy DP which computes ALL floor-value pi(x/v) simultaneously.

Key differences:
- Lucy DP: O(x^{2/3}) time, O(sqrt(x)) space, very simple
- Lehmer: O(x^{2/3}/log x) time, more complex but can exploit
  pre-sieved prime tables to skip computation
- Meissel-Lehmer: Similar complexity, uses x^{1/4} splitting for
  shallower phi recursion at the cost of computing P3

Neither method breaks the O(x^{2/3}) barrier. The improvement is
in constants and cache behavior. The full Deleglise-Rivat/Gourdon
method (as in primecount) adds segmented sieve for the P2 term,
achieving O(x^{2/3}/log^2 x) with excellent constants.
""")


if __name__ == "__main__":
    main()
