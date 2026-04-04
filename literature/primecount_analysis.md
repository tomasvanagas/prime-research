# Primecount Library: Engineering Analysis

**Repository:** [kimwalisch/primecount](https://github.com/kimwalisch/primecount)
**Version analyzed:** v8.4 (April 2026)
**Purpose:** Understanding what makes primecount 100-1000x faster than naive Lucy DP

---

## 1. Algorithm: Gourdon's Variant vs Lucy DP vs Deleglise-Rivat

### 1.1 Lucy DP (Baseline)

Lucy_Hedgehog's algorithm (popularized on Project Euler) computes pi(x) via dynamic
programming on the ~2*sqrt(x) distinct values of floor(x/n). The recurrence is:

    S(v, p) = S(v, p-1) - [S(v/p, p-1) - S(p-1, p-1)]

for each prime p, iterating over all key values v >= p^2.

**Complexity:** O(x^{3/4}) time, O(sqrt(x)) space.
**With Fenwick tree enhancement:** O(x^{2/3} * (log x)^{1/3}) time, O(x^{2/3}) space.

Lucy DP is elegant and simple (~30 lines of Python) but fundamentally limited:
- No segmented computation (entire state must be in memory)
- No natural parallelism (each sieving step depends on the previous)
- No way to decompose into independent subproblems
- Cache-hostile access patterns for large x

### 1.2 Deleglise-Rivat Algorithm

The Deleglise-Rivat (1996) formula decomposes pi(x) as:

    pi(x) = pi(y) + S1(x, a) + S2(x, a) - 1 - P2(x, a)

where y = alpha * x^{1/3}, a = pi(y), and:

- **S1 (Ordinary leaves):** Sum over mu(n) * phi(x/n, c) for squarefree n <= y
  with lpf(n) > p_c. Computed via PhiTiny wheel cache for O(1) evaluation.
- **S2 (Special leaves):** Further decomposed into three parts:
  - **S2_trivial:** Direct calculation of trivial cases (small contribution)
  - **S2_easy:** Special leaves where pi(x/(n*p)) can be looked up in a pi-table
  - **S2_hard:** The bottleneck -- requires actual sieving to evaluate
- **P2:** Second prime power contribution

**Complexity:** O(x^{2/3} / log^2 x) time, O(x^{1/3} * log^2 x) space.

The key insight is that most of Lucy DP's work goes into computing values that
correspond to "special leaves" in a combinatorial tree, and these leaves can be
classified by difficulty and handled with different strategies.

### 1.3 Gourdon's Algorithm (2001) -- Default in primecount

Xavier Gourdon's variant introduces two independent tuning parameters instead of one:

    y = x^{1/3} * alpha_y
    z = y * alpha_z

with constraints 1 <= alpha_y, alpha_z <= x^{1/6}.

The formula becomes:

    pi(x) = A + B + C + D + Phi0 + Sigma

where A, B, C, D, Phi0, and Sigma (actually 7 sub-formulas) represent different
contributions. The two-parameter family allows independent optimization of the
easy-leaf and hard-leaf computational domains.

**Complexity:** O(x^{2/3} / log^2 x) time, O(x^{2/3} * log^3 x) memory.

### 1.4 Why Gourdon Beats Deleglise-Rivat in Practice

Benchmark on AMD EPYC Zen4 (32 cores, 3.7 GHz) for pi(10^22):
- **Deleglise-Rivat:** 1,690.78 seconds
- **Gourdon:** 436.77 seconds (**3.9x faster**)

The practical advantage comes from:
1. **Two alpha parameters** allow finer balancing of A/B/C/D formula work
2. **Independent chunk decomposition** (Gourdon 2002 modification): the hard special
   leaves computation can be split into independent chunks, enabling embarrassingly
   parallel execution without thread synchronization
3. **Better constant factors** in the formulas themselves

For pi(10^25), Gourdon needs only ~330 CPU core hours on EPYC Zen4, versus
Jan Buthe's ~40,000 hours using analytic methods (zeta zeros). That is a 120x
advantage for the combinatorial approach at this scale.

---

## 2. Key Engineering Optimizations

### 2.1 Hard Special Leaves: Counter Array + POPCNT (Most Important Innovation)

The hard special leaves (S2_hard in DR, D formula in Gourdon) are the computational
bottleneck. The standard academic approach uses a **binary indexed tree (Fenwick tree)**
to maintain prefix sums during sieving.

**Problem with Fenwick trees:** Each update/query touches O(log n) cache lines
scattered across memory. For large x, this is catastrophically cache-inefficient.

**Primecount's innovation:** Replace the Fenwick tree with a **linear counter array**
combined with hardware **POPCNT** (population count) instructions.

The approach works as follows:
1. Maintain a **bit array** representing the sieve state (1 bit per candidate)
2. When a count is needed at position i, use POPCNT to count set bits up to i
3. Use a **multi-level counter structure** to avoid scanning the entire array:
   - Level 0: the raw bit array
   - Higher levels: cached prefix sums over blocks of the bit array
   - This gives O(1) amortized counting with O(log log x) overhead per sieve event

**Why this is faster:**
- Bit array is ~8x more compact than byte array, fitting more in cache
- POPCNT on a 64-bit word counts 64 entries in a single instruction (~1 cycle)
- Linear memory access pattern is prefetch-friendly
- Multi-level counters reduce the amortized cost per counting operation

The v7.2 changelog (2021) states this approach has a proven O(log log x) factor
improvement over the Deleglise-Rivat binary indexed tree approach for counting.

As of v8.4 (2026), this is documented in "Hard-Special-Leaves-SIMD-Filtering.pdf"
which describes SIMD-accelerated filtering within the hard special leaves computation.

### 2.2 SIMD Acceleration (AVX512, ARM SVE)

Primecount uses SIMD at multiple levels:

1. **POPCNT vectorization:** AVX512 VPOPCNT (available on Ice Lake+) counts bits
   across 512-bit vectors, processing 512 sieve entries per instruction.
   ARM SVE provides equivalent functionality.

2. **Branchfree bitmask calculations** (v7.14+): SIMD enables branchfree computation
   of sieve masks, avoiding branch misprediction penalties in the inner loop.

3. **D formula SIMD kernels** (v8.4): Branchfree implementations for the D formula
   across AVX512, ARM SVE, and portable fallback code.

4. **Runtime dispatch** (v7.16+): CPUID check moved outside the hot loop; SIMD code
   inlined into the main algorithm to avoid function call overhead. The --status
   option prints which SIMD backend is active (AVX512, ARM SVE, etc.).

The SIMD optimizations provide roughly 4-8x speedup in the inner sieve loops
compared to scalar code, depending on the platform.

### 2.3 Segmented Computation and Cache Optimization

**SegmentedPiTable:** A per-segment pi(x) lookup table sized to fit in L1 cache.
Instead of a global pi-table of size O(x^{2/3}), primecount maintains a small
segment that is recomputed as the sieve advances. This keeps the working set
in L1 cache (~32-64 KB), avoiding expensive L2/L3/DRAM accesses.

**Mod-240 wheel encoding:** The sieve uses 8 bits per byte corresponding to
offsets {1, 7, 11, 13, 17, 19, 23, 29} mod 30, then packs 8 such residues
into a byte. This is a wheel-30 factorization that eliminates all multiples
of 2, 3, and 5, reducing storage by ~73%. The byte-level encoding with
offsets mod 240 (= 8 * 30) enables efficient POPCNT-based counting.

**FactorTable:** A compact table storing the smallest prime factor (or a related
quantity) for numbers up to x^{1/3}. Uses the "fast integer division trick"
(32-bit division when possible, since 64-bit division is 2-5x slower on x86).

### 2.4 Fast Integer Division (libdivide)

Integer division is one of the most expensive single instructions on modern CPUs
(~30-90 cycles for 64-bit division on x86 vs ~3-5 cycles for multiplication).

Primecount uses **libdivide** (by ridiculousfish) to replace division by
multiplication + bit shifts when dividing by the same divisor many times.
A custom **branchfree divider** was added to libdivide specifically for primecount.

Impact: **~40% speedup** for easy special leaves computation.

Exception: On Apple Silicon, native division is fast enough that libdivide
is actually slower, so it is disabled on that platform.

### 2.5 PhiTiny and Partial Sieve Function Caching

**PhiTiny:** For small values of a (the number of primes sieved), phi(x, a) can
be computed in O(1) using precomputed lookup tables based on wheel factorization.
This handles the ordinary leaves (S1) extremely efficiently.

**Phi cache:** Results of phi(x, a) for a <= 100 are cached in a ~16 MB hash table
(sized to slightly exceed L3 cache). This provides >10x speedup for the partial
sieve function. Using a larger cache degrades performance due to cache thrashing,
especially under multi-threading.

### 2.6 32-bit vs 64-bit Arithmetic

Christian Bau's insight: use 32-bit integers whenever the value fits, since:
- 32-bit division is 2-5x faster than 64-bit on x86
- 32-bit operations use narrower execution units
- More values fit in cache lines and SIMD registers

Primecount systematically uses 32-bit arithmetic for intermediate values below 2^32,
with 64-bit and 128-bit paths for larger inputs. Separate template instantiations
(not runtime checks) ensure zero overhead from the type dispatch.

---

## 3. OpenMP Parallelization Strategy

### 3.1 Gourdon's Independent Chunk Decomposition

The fundamental parallelization enabler is Gourdon's 2002 modification: the hard
special leaves computation is restructured so that it decomposes into **independent
chunks** that can be processed without inter-thread communication.

In the Deleglise-Rivat formulation, S2_hard requires maintaining a global sieve
state that is updated sequentially. Gourdon's reformulation partitions the work
into segments along one axis, where each segment can be processed independently.

### 3.2 Load Balancing

Primecount features a **novel adaptive load balancer** shared across all algorithm
implementations:

- **Dynamic chunk sizing:** Work chunks are sized dynamically based on observed
  throughput, not statically partitioned. This handles the irregular work
  distribution inherent in prime counting (some ranges have more special leaves).

- **LoadBalancerAC** (v7.0+): Complete rewrite achieving "100% CPU core usage"
  with minimal synchronization overhead. Uses dynamic/adaptive scheduling.

- **Scaling:** The load balancer scales to hundreds of CPU cores. This is critical
  for modern server CPUs (AMD EPYC with 64-128 cores).

- **Clang workaround** (v6.4): Discovered that Clang's OpenMP implementation had
  scaling issues with dynamic scheduling; implemented workarounds.

### 3.3 Formula-Level Parallelism

The Gourdon formula components (A, C, B, D, Phi0, Sigma) have different
computational costs and some can run concurrently:

- **A + C** are typically computed together
- **B** is independent and can run in parallel
- **D** (hard special leaves) is the bottleneck, gets most cores
- **Phi0** and **Sigma** are relatively cheap

Within each formula, the segmented computation enables further parallel decomposition
via OpenMP parallel for loops over segments.

### 3.4 Thread Reuse

Rather than spawning new threads for each formula, primecount reuses thread pools
across the computation phases. This avoids the overhead of thread creation/destruction
(v6.1+).

---

## 4. Complexity Analysis: Where the 100-1000x Comes From

### 4.1 Asymptotic Comparison

| Method              | Time Complexity          | Space              |
|---------------------|--------------------------|--------------------|
| Lucy DP (basic)     | O(x^{3/4})              | O(sqrt(x))         |
| Lucy + Fenwick      | O(x^{2/3} log^{1/3} x)  | O(x^{2/3})         |
| Meissel-Lehmer      | O(x^{2/3})              | O(x^{1/3})         |
| Deleglise-Rivat     | O(x^{2/3} / log^2 x)    | O(x^{1/3} log^2 x) |
| Gourdon (primecount)| O(x^{2/3} / log^2 x)    | O(x^{1/3} log^3 x) |

For x = 10^22:
- Lucy basic: ~10^{16.5} operations
- Gourdon: ~10^{14.7} / (log 10^22)^2 ~ 10^{14.7} / 2567 ~ 10^{11.3} operations

That is a **~10^5 asymptotic ratio** before constant factors.

### 4.2 Constant Factor Breakdown

On top of the asymptotic advantage, primecount gains from engineering:

| Optimization                    | Estimated Speedup Factor |
|---------------------------------|--------------------------|
| POPCNT counter vs Fenwick tree  | 5-10x (cache efficiency) |
| SIMD (AVX512 inner loops)       | 4-8x                     |
| libdivide (division elimination)| 1.4x (easy leaves)       |
| 32-bit arithmetic where possible| 2-3x (division-heavy code)|
| Segmented L1-cached PiTable     | 2-4x (memory latency)   |
| Mod-240 wheel sieve             | ~3.3x (storage reduction)|
| PhiTiny + phi cache             | 10x+ (partial sieve fn)  |
| OpenMP parallelism (32 cores)   | 20-30x (near-linear)     |

**Combined effect:** The multiplicative product of these factors easily reaches
100-1000x over a naive Lucy DP implementation, depending on input size and hardware.

### 4.3 Practical Benchmarks

| x        | Lucy DP (est.) | primecount (Gourdon) | Ratio     |
|----------|----------------|----------------------|-----------|
| 10^12    | ~2s            | ~0.01s               | ~200x     |
| 10^15    | ~77s (optimized)| <1s                  | ~100x     |
| 10^22    | infeasible     | 436s (32 cores)      | --        |
| 10^25    | infeasible     | 330 core-hours       | --        |

For x >= 10^18, Lucy DP becomes impractical while primecount remains feasible.

---

## 5. Relevance to Our Research

### 5.1 Implications for the O(x^{2/3}) Barrier

Primecount represents the **practical limit** of combinatorial prime counting:
- Asymptotically O(x^{2/3} / log^2 x) -- the best known combinatorial complexity
- Engineering optimizations squeeze out every constant factor
- 32-core parallelism provides another ~25x

For p(10^100): even with primecount's optimizations, the computation requires
~(10^100)^{2/3} / log^2(10^100) ~ 10^{66.7} / 53000 ~ 10^{62} operations.
At 10^15 ops/sec = 10^47 seconds. **Still completely infeasible.**

### 5.2 What We Already Knew, Confirmed

1. The combinatorial approach (Meissel-Lehmer-LMO-DR-Gourdon) is a family of
   increasingly refined decompositions of the SAME underlying computation.
   They all compute the Legendre sieve via different organizational strategies.

2. The log^2 x factor improvement over bare O(x^{2/3}) comes from careful
   decomposition into ordinary/easy/hard leaves and optimized handling of each.

3. No algorithmic breakthrough has improved the exponent 2/3 for combinatorial
   methods since Lagarias-Miller-Odlyzko (1985). All subsequent work (DR, Gourdon,
   Staple) improves constants and log factors only.

4. The hard special leaves are the irreducible bottleneck. All the SIMD/POPCNT/
   cache optimization in the world cannot change the fundamental O(x^{2/3}) scaling.

### 5.3 Our v10 Implementation Gap

Our best implementation (`algorithms/v10_c_accelerated.py`) achieves O(x^{2/3})
with p(10^9) in 0.175s. Compared to primecount:

- We lack the log^2 x denominator (Gourdon decomposition)
- We lack POPCNT-based counting (using Python/basic C)
- We lack SIMD acceleration
- We lack segmented cache-aware computation
- We lack multi-threaded parallelism

Closing this gap could give us ~100-1000x improvement at the same exponent,
but would NOT change the fundamental infeasibility for p(10^100).

### 5.4 The Only Path Forward

Primecount's engineering excellence confirms that the combinatorial approach is
fully optimized. The path to p(10^100) in <1 second requires either:
1. A fundamentally new algorithm with sublinear exponent (< 2/3), or
2. An analytic method avoiding the O(sqrt(x)) zeta zero barrier, or
3. A complexity-theoretic breakthrough (pi(x) in NC? #TC^0 in NC?)

None of these exist today. This is consistent with our project status.

---

## 6. Technical Architecture Summary

```
primecount v8.4 Architecture
=============================

Input x (up to 10^31)
  |
  v
Alpha optimization (auto-tune alpha_y, alpha_z)
  |
  v
+-- Phi0 + Sigma (cheap, O(x^{1/3}))
|
+-- A + C formulas (medium, parallelized)
|     |-- SegmentedPiTable (L1-cached)
|     |-- FactorTable (compact, 32-bit where possible)
|
+-- B formula (medium, parallelized)
|     |-- PhiTiny wheel cache
|     |-- libdivide for fast division
|
+-- D formula (BOTTLENECK, most cores allocated)
      |-- Gourdon independent chunk decomposition
      |-- Linear counter array + POPCNT (replaces Fenwick tree)
      |-- SIMD kernels: AVX512 / ARM SVE / portable
      |-- Branchfree bitmask calculations
      |-- Mod-240 wheel sieve encoding
      |-- Adaptive load balancer (scales to 100+ cores)
      |-- 32-bit arithmetic fast path for x < 2^64

All components: OpenMP parallelized
               128-bit int support for x > 2^64
               Double-check mode with alternative alpha values
```

---

## References

1. Deleglise, M.; Rivat, J. "Computing pi(x): The Meissel, Lehmer, Lagarias,
   Miller, Odlyzko Method." Mathematics of Computation 65:213 (1996), 235-245.
2. Gourdon, X. "Computation of pi(x): improvements to the Meissel, Lehmer,
   Lagarias, Miller, Odlyzko, Deleglise and Rivat method." (2001).
3. Staple, D. "The combinatorial algorithm for computing pi(x)." MSc Thesis,
   Dalhousie University (2015). arXiv:1503.01839.
4. Lagarias, J.; Miller, V.; Odlyzko, A. "Computing pi(x): The Meissel-Lehmer
   method." Mathematics of Computation 44 (1985), 537-560.
5. Oliveira e Silva, T. "Computing pi(x): the combinatorial method." Revista
   do DETUA 4:6 (2006), 759-768.
6. Walisch, K. primecount documentation: Hard-Special-Leaves-SIMD-Filtering.pdf,
   Easy-Special-Leaves.pdf, Partial-Sieve-Function.pdf. GitHub repository.
7. griff's math blog, "Lucy's Algorithm + Fenwick Trees" (2023).
   https://gbroxey.github.io/blog/2023/04/09/lucy-fenwick.html

---

*Analysis prepared April 2026 for the Prime Research project.*
*Sources: GitHub kimwalisch/primecount, arXiv, Codeforces, MathWorld, ManKier.*
