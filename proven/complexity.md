# Complexity Bounds for Prime Counting and nth Prime

---

## Upper Bounds: Algorithms for pi(x)

### Combinatorial Methods

| Algorithm | Time | Space | Year |
|-----------|------|-------|------|
| Sieve of Eratosthenes | O(x log log x) | O(x) | Ancient |
| Legendre | O(x) | O(sqrt(x)) | 1808 |
| Meissel | O(x/(log x)^3) | -- | 1870 |
| Meissel-Lehmer | O(x/(log x)^4) | O(x^{1/3}/log x) | 1959 |
| LMO (Lagarias-Miller-Odlyzko) | O(x^{2/3}/log x) | O(x^{1/3}) | 1985 |
| Deleglise-Rivat | O(x^{2/3}/log^2 x) | O(x^{1/3} log^3 x) | 1996 |
| Gourdon | O(x^{2/3}/log^2 x) | O(x^{2/3} log^3 x) | 2001 |
| Lucy Hedgehog DP | O(x^{3/4}) basic, O(x^{2/3}) with Fenwick | O(sqrt(x)) | 2013 |

### Analytic Methods

| Algorithm | Time | Space | Year |
|-----------|------|-------|------|
| Lagarias-Odlyzko | O(x^{1/2+epsilon}) | O(x^{1/4+epsilon}) | 1987 |
| Lagarias-Odlyzko (space-efficient) | O(x^{3/5+epsilon}) | O(x^epsilon) | 1987 |

### For p(n) Specifically (Aggarwal 2025)

| Method | Complexity |
|--------|-----------|
| Best sieve | O(n log n / log log n) |
| Binary search + combinatorial pi(x) | O(sqrt(n) * (log n)^4) |
| Binary search + analytic pi(x) | O(sqrt(n) * polylog) conditional on RH+Cramer |

**Best proven:** O(sqrt(n) * (log n)^4) unconditionally.

---

## Lower Bounds

| Bound | Applies To | Source |
|-------|-----------|--------|
| Omega(log x) | Any algorithm for pi(x) (trivial: read input) | - |
| Not in AC^0 | PRIMES decision | Allender-Saks-Shparlinski 2001 |
| Not in AC^0[p] | PRIMES decision | Allender-Saks-Shparlinski 2001 |

**CRITICAL: No super-logarithmic unconditional lower bound exists for pi(x).**

---

## Complexity Class Memberships

### PRIMES (decision: "Is n prime?")

| Class | Status |
|-------|--------|
| P | YES (AKS 2002) |
| coRP | YES (Miller-Rabin) |
| RNC | YES (randomized NC) |
| NC | UNKNOWN (major open problem) |
| TC^0 | UNKNOWN |
| L (logspace) | UNKNOWN |
| P-complete | UNKNOWN (would rule out NC) |
| AC^0 | NO (proven) |
| AC^0[p] | NO (proven) |

### pi(x) (counting function)

| Class | Status |
|-------|--------|
| P (binary input model) | UNKNOWN -- all known algorithms exponential in log(x) |
| FP (unary model) | YES -- O(x^{2/3}) is sublinear |
| NC | UNKNOWN |
| TC^0 | UNKNOWN |
| #P-hard | NOT KNOWN (but exact pi(x) is related) |

### p(n) (nth prime)

| Result | Source |
|--------|--------|
| p(n) in L (logspace) | Follows from AKS: iterate and count |
| p(n) in NC^2 | Follows from L subset NC^2 |
| Circuit size for p(n) | UNKNOWN -- polylog-size circuits not ruled out |

---

## The Input Model Distinction

**Binary input (n = log x bits):** All known pi(x) algorithms are EXPONENTIAL.
Best is 2^{Omega(n^{1/2})}. "Is pi(x) in P?" is wide open.

**Unary input (input size = x):** pi(x) is computable in SUBLINEAR time
O(x^{2/3}). Question: how far below x can we go?

---

## Circuit Complexity of pi(x) (Session 12)

### Best Known Circuits (in terms of N = log x)

| Method | Circuit Size | Circuit Depth |
|--------|-------------|---------------|
| Lucy DP (standard) | O(2^{N/2}) | O(2^{N/2}/N) [= pi(sqrt(x))] |
| Lucy DP (parallel) | O(2^{N/2}) | O(2^{N/3}/N) [= pi(x^{1/3})] |
| Meissel-Lehmer | O(2^{2N/3}) | O(log N) [= O(log log x)] |
| Lagarias-Odlyzko | O(2^{N/2+ε}) | O(poly(N)) |

All are EXPONENTIAL in N. No polynomial-size circuit is known.

### Structural Results (Session 12)

1. **Lucy DP DAG depth = pi(sqrt(x)) exactly**: The computation has an
   unavoidable sequential chain through S(x) at every prime step.
   Depth/pi(sqrt(x)) = 1.000 for all x tested (100 to 100000).

2. **Parallel round depth = pi(x^{1/3})**: With maximal parallelism,
   pi(x^{1/3}) sequential rounds suffice. Ratio → 1.00 as x grows.

3. **Floor-value mapping matrices are non-commutative**: Zero commuting
   pairs in all tests. Sieve steps cannot be reordered.

4. **Linear transformation is full-rank**: pi(x) depends on ~80% of
   the O(sqrt(x)) floor-value initial conditions with large coefficients.

5. **pi(x) mod m has invariant conditional entropy**: H = 0.537 bits
   for all moduli m (tested 2-30). Computing pi(x) mod m is as hard
   as pi(x) itself.

### Equivalence Result

**"Is pi(x) computable in polylog(x) time?" ⟺ "Is pi(x) in NC?"**

Both require polynomial-size (poly(N)) circuits. All known approaches
produce exponential-size (2^{Theta(N)}) circuits because they compute
O(sqrt(x)) = O(2^{N/2}) intermediate values (floor-values or zeta zeros).

### Session 13: TC^0 Primality Test Analysis

| Test | In TC^0? | Deterministic? | Condition |
|------|----------|---------------|-----------|
| MR(2) | YES | NO | - |
| MR(2,3,...,37) | YES | YES for n<3×10^24 | unconditional |
| MR(2,...,2ln²(n)) | YES | YES for all n | GRH |
| Strong Lucas (Selfridge) | YES | NO | - |
| BPSW (MR(2)+Strong Lucas) | YES | YES for n<2^64 | unconditional |
| BPSW | YES | YES for all n | BPSW conjecture |
| QFT (Grantham) | YES | NO (error<1/7710) | - |
| AKS | UNKNOWN | YES | unconditional |
| Wilson's | NO (needs 2^N mults) | YES | unconditional |

**Key:** PRIMES is in NONUNIFORM TC^0 unconditionally (for any fixed N,
a finite MR base set suffices). PRIMES is in UNIFORM TC^0 conditional
on GRH or BPSW correctness.

**Prime indicator GF(2) structure:** ANF degree = Theta(N), 50% sparsity.
Indistinguishable from random over GF(2). Confirms PRIMES not in AC^0[2].

### What Would Be Needed

A polylog(x) algorithm must avoid:
- Enumerating floor values {floor(x/k)} (O(sqrt(x)) values)
- Summing over zeta zeros (O(sqrt(x)) terms)
- Testing individual integers for primality (O(x) tests)
- Sieve-type inclusion-exclusion (exponential terms)

No mathematical framework is currently known that could achieve this.

---

## Related Summatory Functions

| Function | Best Combinatorial | Best Analytic |
|----------|-------------------|--------------|
| pi(x) = sum 1_{p prime} | O(x^{2/3}/log^2 x) | O(x^{1/2+eps}) |
| M(x) = sum mu(n) | O(x^{2/3}) | O(x^{1/2+eps}) |
| sum phi(n) | O(x^{2/3}) | O(x^{1/2+eps}) |
| sum d(n) | O(x^{1/3}) | O(x^{1/3}) |

The O(x^{2/3}) barrier appears for ALL multiplicative function sums involving
primes or Mobius. Lower bounds for any would likely transfer to others.

---

## Practical Records (primecount library, 2026)

- pi(10^25) in 330 CPU core hours (Gourdon, AMD EPYC Zen4)
- pi(10^29) computed (2022 record)
- Combinatorial is 100x faster than analytic in practice for all reachable x
- primecount v8.4 (April 2026): SIMD-accelerated Gourdon D formula with AVX512

---

## Key References

- Aggarwal (2025): arXiv:2510.16285
- Lagarias-Odlyzko (1987): Analytic method
- Deleglise-Rivat (1996): Combinatorial method
- AKS (2002): PRIMES is in P
- Allender-Saks-Shparlinski (2001): PRIMES not in AC^0
- Hesse-Allender-Barrington (2002): Division in TC^0
- Guth-Maynard (2024): arXiv:2405.20552 (zero density improvement)
- Kim Walisch: github.com/kimwalisch/primecount
