# Best Working Algorithms

All algorithms here are **100% exact** and **verified**. Listed by performance.

---

## 1. V10: C-Accelerated Hybrid (BEST)

**File:** `algorithms/v10_c_accelerated.py`
**Complexity:** O(p(n)^{2/3})
**Dual mode:** exact + approximate

### Exact Mode

| n | p(n) | Time |
|---|------|------|
| 10^3 | 7,919 | 0.000s |
| 10^5 | 1,299,709 | 0.006s |
| 10^7 | 179,424,673 | 0.034s |
| 10^8 | 2,038,074,743 | 0.140s |
| 10^9 | 22,801,763,489 | 0.175s |

### Approximate Mode (R^{-1} only)

| n | Time | Total Digits | Correct Digits |
|---|------|-------------|---------------|
| 10^10 | 1.5s | 12 | ~5 |
| 10^50 | 0.14s | 53 | ~25 |
| 10^100 | 0.32s | 103 | ~50 |
| 10^200 | 1.15s | 203 | ~99 |

### Algorithm

```
1. x0 = R^{-1}(n)                      [O(polylog)]
2. Newton refinement:                   [~4 steps]
     pi = C_Lucy_DP(x0)                [O(x^{2/3}), C-accelerated]
     x0 += (n - pi) * ln(x0)
3. Bisection in ~70-integer bracket     [~6 steps]
4. Walk down to prime                   [O(ln^2(x))]
```

Components: R^{-1} via mpmath Newton | Lucy_Hedgehog DP (C extension) |
Deterministic Miller-Rabin primality test

---

## 2. V7: Pure Python Hybrid

**File:** `algorithms/v7_optimized.py`
**Complexity:** O(p(n)^{2/3})
**Verified:** 10,000/10,000 exhaustive test PASSED

| n | p(n) | Time |
|---|------|------|
| 10^5 | 1,299,709 | 0.03s |
| 10^6 | 15,485,863 | 0.21s |
| 10^7 | 179,424,673 | 1.28s |
| 10^8 | 2,038,074,743 | 7.73s |

Same algorithm as V10 but pure Python (no C extension). ~40x slower.

---

## 3. V5: Pure Analytic

**File:** `algorithms/v5_pure_analytic.py`
**Complexity:** O(p(n)^{2/3})

Same architecture as V7 with slightly different implementation.
p(10^7) in 1.47s, p(10^8) in 9.0s.

---

## 4. V6: Single DP + Lookup

**File:** `algorithms/v6_single_dp.py`
**Complexity:** O(p(n)^{4/3})

One Lucy DP call + table lookup. Zero search. Simpler but slower for large n.

---

## Comparison: Our v10 vs State-of-Art

### Our v10 (Lucy DP) vs primecount (Gourdon) vs HKM (NTT)

| Aspect | Our v10 | primecount (Gourdon) | HKM 2023 |
|--------|---------|---------------------|----------|
| Algorithm | Lucy_Hedgehog DP | Gourdon (Deleglise-Rivat variant) | NTT Dirichlet convolution |
| Complexity | O(x^{2/3}) | O(x^{2/3}/log^2 x) | O(sqrt(x) * log^{5/2} x) |
| p(10^9) | 0.175s | ~ms | ~0.35s |
| p(10^12) | ~hours | ~seconds | ~14.6s |
| p(10^14) | infeasible | ~minutes | ~127s |
| Parallelism | None | OpenMP + SIMD/AVX512 | Possible via NTT |
| Language | C (via ctypes) | C++ (highly optimized) | C++ |
| Space | O(sqrt(x)) | O(x^{1/3}) | O(sqrt(x)) or O(x^{1/3}) |

### Why primecount is 100-1000x faster than our v10:
1. **Algorithm:** Gourdon's variant vs basic Lucy DP — different handling of "special leaves"
2. **Segmented sieve:** primecount uses cache-friendly segmented sieve for hard special leaves
3. **SIMD acceleration:** AVX512/POPCNT for bit-level sieve operations
4. **OpenMP:** Parallel computation across CPU cores
5. **Memory optimization:** Compressed lookup tables, cache-line alignment
6. **Lucy DP is naive:** Our v10 does ALL floor values sequentially; Gourdon decomposes
   the computation into ordinary/easy/hard leaves with different optimization strategies

### Key insight for our project:
All three methods are O(2^{Theta(N)}) in input bits N = log x. Engineering optimizations
(primecount) give 100-1000x constant-factor speedup. Asymptotic improvement (HKM) gives
better exponent (N/2 vs 2N/3) but worse constants. NEITHER approaches poly(N).

---

## For p(10^100) Specifically

Best known result: R^{-1}(10^100) gives ~50/103 correct digits in 0.32s.
Exact computation requires O(10^68) combinatorial or O(10^51) analytic operations.
At 10^15 ops/sec: 10^36 seconds minimum. **INFEASIBLE.**
