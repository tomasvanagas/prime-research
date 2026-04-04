# Aggarwal 2025: "A Note on Algorithms for Computing p_n"

## Analysis and Documentation

**Paper:** Ansh Aggarwal, "A Note on Algorithms for Computing p_n"
**Source:** arXiv:2510.16285, October 2025
**Affiliation:** University of Wisconsin-Madison
**Status:** Preprint (not yet peer-reviewed as of April 2026)

---

## 1. What the Paper Actually Does

**This is NOT a new algorithm.** Aggarwal provides a rigorous complexity analysis of
combining existing tools:

1. **Dusart's bounds** on p_n to establish a search interval
2. **Binary search** on that interval using a pi(x) oracle
3. **Hirsch-Kessler-Mendlovic (HKM) algorithm** as the pi(x) oracle

The main contribution is proving that this combination is provably superior to any
sieve-based approach for computing individual primes.

### Key Results

| Approach | Complexity | Assumptions |
|----------|-----------|-------------|
| Binary search + HKM | O(sqrt(n) * (log n)^4) | Unconditional |
| Improved hybrid (sieve) | O(sqrt(n) * (log n)^{7/2} * log log n) | RH + Cramer |
| Sublinear sieve | O(n * log n / log log n) | Unconditional |
| Segmented sieve | O(s * log log n) where s = interval | Unconditional |

### The Derivation of O(sqrt(n) * (log n)^4)

1. **Dusart bounds:** For n >= 6, p_n is in [n(log n + log log n - 1), n(log n + log log n)]
   - Interval size is O(n), so binary search needs O(log n) evaluations of pi(x)
2. **HKM oracle:** pi(x) computable in O(sqrt(x) * (log x)^{5/2})
3. **Substitution:** x ~ n log n (by PNT), so pi(x) costs O(sqrt(n log n) * (log(n log n))^{5/2})
   = O(sqrt(n) * (log n)^3)
4. **Total:** O(log n) evaluations * O(sqrt(n) * (log n)^3) each = O(sqrt(n) * (log n)^4)

### Proof that Sieves Cannot Compete

**Theorem 4.2:** The sublinear sieve has complexity O(n log n / log log n). The ratio
to binary search time is (n log n / log log n) / (sqrt(n) * (log n)^4) ~ sqrt(n) / (log n)^3,
which diverges. Hence the sieve is asymptotically worse by a factor of ~sqrt(n).

For interval sieves to compete, the interval size must be < sqrt(n) * (log n)^4 / log log n.
Current unconditional prime gap bounds (Dusart/Axler) give intervals of size O(n / log n),
which is much larger. Only under RH + Cramer can the interval be narrowed to
O(sqrt(n) * (log n)^{7/2}).

---

## 2. The Real Hero: Hirsch-Kessler-Mendlovic (HKM) Algorithm

The Aggarwal paper's complexity relies entirely on the HKM algorithm for pi(x).

**Paper:** Dean Hirsch, Ido Kessler, Uri Mendlovic, "Computing pi(N): An elementary
approach in O~(sqrt(N)) time"
**Source:** arXiv:2212.09857 (December 2022, revised August 2023)
**Published:** Mathematics of Computation (received March 2024, published electronically
November 2024). Peer-reviewed in a top journal.

### What HKM Achieves

| Variant | Time | Space |
|---------|------|-------|
| Time-optimal | O~(sqrt(N)) | O~(sqrt(N)) |
| Space-time tradeoff | O~(N^{8/15}) | O~(N^{1/3}) |

Where O~ hides polylogarithmic factors. The exact complexity is O(sqrt(N) * (log N)^{5/2}).

### Core Technique

HKM is an **elementary** (non-analytic) algorithm:
- Does NOT use complex analysis, zeta zeros, or Perron's formula
- Does NOT use arbitrary-precision complex arithmetic
- Uses **Dirichlet convolutions** and **Number Theoretic Transforms (NTT)**
- Extends the Dirichlet hyperbola method with fast convolution techniques
- Key innovation: applies NTT-based fast multiplication to Dirichlet series
  computations that previously required O(N^{2/3}) via Lucy_Hedgehog/Meissel-Lehmer

The algorithm components (from the C++ implementation):
1. NTT (Number Theoretic Transform)
2. Mobius function calculation
3. INTT (Inverse NTT)
4. Small prime convolution
5. Error correction

### Practical Performance (HKM implementation)

Benchmarks on Intel Core i3-1005G1, single core, memory tradeoff factor 5:

| N | Time |
|---|------|
| 10^7 | 0.06 sec |
| 10^10 | 1.23 sec |
| 10^14 | 127 sec |

### Comparison with Existing Methods

| Method | Time Complexity | Space | Practical Speed | Analytic? |
|--------|----------------|-------|-----------------|-----------|
| Lucy_Hedgehog DP | O(N^{2/3}) | O(N^{1/2}) | Very fast | No |
| Deleglise-Rivat/Gourdon | O(N^{2/3}/log^2 N) | O(N^{1/3}) | Fastest practical | No |
| Lagarias-Odlyzko | O(N^{1/2+eps}) | O(N^{1/2+eps}) | 100x slower than combinatorial | Yes |
| Platt (rigorous L-O) | O(N^{1/2+eps}) | O(N^{1/2+eps}) | 40,000 CPU hrs for pi(10^25) | Yes |
| **HKM** | O(N^{1/2} * log^{5/2} N) | O(N^{1/2}) | ~400x slower than primecount | No |

**Critical practical comparison:** For pi(10^14):
- primecount (Gourdon): ~0.3 seconds (estimated from pi(10^25) in 330 CPU hours)
- HKM: 127 seconds (from benchmark)

HKM is roughly **400x slower** in practice despite better asymptotics. The crossover
point where HKM becomes faster than combinatorial methods is likely around N ~ 10^{30}
or higher -- a regime that is currently computationally unreachable for either method.

---

## 3. Complexity Analysis for p(10^100)

### Using Aggarwal's binary search + HKM:

p(10^100) ~ 2.35 * 10^{102} (by PNT)

The nth prime with n = 10^100:
- sqrt(n) = 10^{50}
- (log n)^4 = (100 * ln 10)^4 = (230.26)^4 ~ 2.81 * 10^9
- **Total: ~2.81 * 10^{59} operations**

### Using combinatorial (Deleglise-Rivat):

- x ~ 10^{102}, so x^{2/3} ~ 10^{68}
- **Total: ~10^{68} operations**

### Using Lagarias-Odlyzko analytic:

- x^{1/2+eps} ~ 10^{51+eps}
- **Total: ~10^{51} operations** (but huge constants)

### Feasibility comparison for p(10^100):

| Method | Operations | At 10^{15} ops/sec | Years |
|--------|-----------|-------------------|-------|
| Combinatorial | 10^{68} | 10^{53} sec | 3 * 10^{45} |
| **Aggarwal/HKM** | **3 * 10^{59}** | **3 * 10^{44} sec** | **10^{37}** |
| Analytic (L-O) | 10^{51+eps} | 10^{36+eps} sec | 3 * 10^{28+eps} |

**Verdict:** HKM/Aggarwal saves ~9 orders of magnitude over combinatorial methods for
p(10^100), but is still ~10^{37} years -- utterly infeasible. The analytic method
(Lagarias-Odlyzko) remains theoretically better by ~8 more orders of magnitude,
but with much worse practical constants.

---

## 4. Implementation Feasibility

### HKM Implementation (exists, works)
- C++ and Python implementations available: https://github.com/PrimeCounting/PrimeCounting
- Published in Mathematics of Computation (peer-reviewed)
- No complex arithmetic needed (unlike Lagarias-Odlyzko)
- Memory: O(sqrt(N)) or O(N^{1/3}) with tradeoff
- Practical: works today, tested up to at least 10^{14}

### Obstacles for Large N
1. **Memory:** For N = 10^{102}, sqrt(N) ~ 10^{51} entries -- impossible
2. **Time:** 10^{59} operations -- impossible
3. **The space-time tradeoff** (N^{8/15} time, N^{1/3} space) at N = 10^{102}:
   - Time: 10^{102 * 8/15} ~ 10^{54.4} -- still impossible
   - Space: 10^{34} -- still impossible

### Advantages over Lagarias-Odlyzko
- No zeta zero computation needed
- No numerical stability issues with complex arithmetic
- No precision management for oscillatory integrals
- Simpler to implement and verify
- Same asymptotic class: O(N^{1/2+eps}) vs O(N^{1/2} * polylog(N))

### Disadvantages vs Lagarias-Odlyzko
- HKM's polylog factors may be larger in practice
- L-O can use precomputed zeta zeros (amortized cost)
- L-O has been used for record computations (pi(10^25))

---

## 5. Related 2025-2026 Developments

### Published/Accepted
- **Guth-Maynard (2024, Annals 2025):** Zero-density bound 3/5 -> 13/25. Does NOT
  improve algorithms but tightens theoretical guarantees. Could marginally reduce
  the number of zeta zeros needed for rigorous analytic pi(x).

- **Tao-Trudgian-Yang (Jan 2025):** Systematic exponent pair database (ANTEDB).
  New zero-density estimates via computational optimization. Infrastructure, not
  a new algorithm.

- **Gafni-Tao (May 2025):** Explicit exceptional intervals for PNT in short intervals.
  Tighter guarantees for binary search bracketing.

- **Pascadi (May 2025):** Equidistribution to x^{5/8} without Selberg eigenvalue
  conjecture. Creeping toward the x^{2/3} algorithmic threshold but not there yet.

### No New Algorithmic Breakthroughs
As of April 2026, no paper has improved upon:
- O(N^{2/3}/log^2 N) for practical combinatorial pi(x) [Deleglise-Rivat 1996]
- O~(sqrt(N)) for elementary pi(x) [HKM 2022/2024]
- O(N^{1/2+eps}) for analytic pi(x) [Lagarias-Odlyzko 1987]

### Deuring-Heilbronn Phenomenon
Ongoing work (Bristol seminar March 2025) on explicit versions, but no computational
implications for pi(x) algorithms.

---

## 6. Impact on Project Goals

### Does Aggarwal change the barrier?

**No.** The paper provides the tightest known complexity analysis for computing p(n)
but introduces no new algorithm. The complexity O(sqrt(n) * (log n)^4) was already
implicitly achievable by combining binary search with HKM.

### Does HKM change the barrier?

**Marginally.** HKM replaces the analytic method (Lagarias-Odlyzko) as the best
*elementary* algorithm for pi(x), matching the O(N^{1/2+eps}) analytic bound without
requiring complex analysis. However:

- The barrier for p(10^100) remains: ~10^{59} operations minimum (HKM) or ~10^{51}
  operations minimum (L-O), both utterly infeasible
- The gap between the trivial Omega(log x) lower bound and the O(x^{1/2+eps}) upper
  bound remains one of the largest in computational complexity
- No polylogarithmic algorithm is in sight

### What should be updated in the project?

1. **BEST_ALGORITHMS.md** should note HKM as the best elementary algorithm for pi(x),
   matching analytic methods asymptotically
2. **state_of_art_2026.md** already documents Aggarwal correctly
3. The project's v10_c_accelerated.py (O(x^{2/3})) remains the best *practical*
   implementation for reachable values of x

---

## 7. Summary Table

| Question | Answer |
|----------|--------|
| Is this a new algorithm? | No -- complexity analysis of existing tools |
| What pi(x) algorithm? | Hirsch-Kessler-Mendlovic (elementary, O~(sqrt(N))) |
| Based on Lagarias-Odlyzko? | No -- HKM is elementary, avoids complex analysis |
| Peer-reviewed? | HKM: yes (Math. Comp. 2024). Aggarwal: no (preprint) |
| Practical improvement? | Not for reachable N; primecount ~400x faster at 10^{14} |
| Changes the barrier? | No. Same O(x^{1/2+eps}) regime |
| Helps with p(10^100)? | Reduces from 10^{68} to 10^{59} ops. Still 10^{37} years |
| Implementation available? | Yes: github.com/PrimeCounting/PrimeCounting |

---

## References

- [Aggarwal 2025] A. Aggarwal, "A Note on Algorithms for Computing p_n," arXiv:2510.16285
- [HKM 2022/2024] D. Hirsch, I. Kessler, U. Mendlovic, "Computing pi(N): An elementary
  approach in O~(sqrt(N)) time," arXiv:2212.09857; Mathematics of Computation (2024)
- [Lagarias-Odlyzko 1987] J. Lagarias, A. Odlyzko, "Computing pi(x): An Analytic Method,"
  J. reine angew. Math. 387 (1987)
- [Deleglise-Rivat 1996] M. Deleglise, J. Rivat, "Computing pi(x): The Meissel, Lehmer,
  Lagarias, Miller, Odlyzko method," Math. Comp. 65 (1996)
- [Dusart 2010] P. Dusart, "Estimates of functions over primes without R.H.,"
  arXiv:1002.0442 (2010)

---

Generated: 2026-04-04
