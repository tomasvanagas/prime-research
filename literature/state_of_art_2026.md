# State of the Art: Prime Computation (2024-2026)

Consolidated survey of all developments relevant to computing p(n) exactly.
Deduplicated from 6 overlapping source files.

Last updated: 2026-04-04

---

## 1. Prime Counting Algorithms

### 1.1 Combinatorial Methods (Best Practical)

**Current best asymptotic:** O(x^{2/3} / log^2 x) -- Deleglise-Rivat / Gourdon (1996/2001)

**Best implementation:** Kim Walisch's primecount library
- Version 8.4 released 2026-04-01
- Implements Gourdon's algorithm with SIMD acceleration (AVX512, ARM SVE)
- primecount 8.3 (March 2026): bidirectional clustered easy leaves, faster compression
- primecount 8.1 (January 2026): 15% faster AC.cpp
- primecount 8.0 (December 2025): new 128-bit API
- pi(10^25) computed in 330 CPU core hours on AMD EPYC Zen4
- The SIMD-accelerated Gourdon D formula in v8.4 vectorizes the hard special leaves bottleneck, estimated 2-4x speedup on modern CPUs

**Relevance to project:** This is the engine behind our best practical approach (binary search on pi(x)). Constant-factor improvements continue but no asymptotic breakthrough.

**Does it change the barrier?** No. Still O(x^{2/3}).

### 1.2 Hirsch-Kessler-Mendlovic: Elementary O~(sqrt(N)) (2022/2024)

**Source:** arXiv:2212.09857; Mathematics of Computation (2024, peer-reviewed)
**Implementation:** https://github.com/PrimeCounting/PrimeCounting (C++ and Python)

**The breakthrough:** First ELEMENTARY (non-analytic) algorithm matching the analytic
O(x^{1/2+eps}) bound for pi(x). Uses Dirichlet convolutions via NTT (Number Theoretic
Transform) to speed up the combinatorial approach from O(x^{2/3}) to O(sqrt(x) * log^{5/2} x).

**Key features:**
- No complex analysis, no zeta zeros, no arbitrary-precision complex arithmetic
- Space-time tradeoff: O~(sqrt(N)) time/space, or O~(N^{1/3}) space with O~(N^{8/15}) time
- Also improves: Mertens function to O~(sqrt(N)) (from O~(N^{0.6})), Euler totient sum, etc.

**Benchmarks (single core, Intel i3-1005G1):**
| N | Time |
|---|------|
| 10^10 | 1.23s |
| 10^12 | 14.6s |
| 10^14 | 127s |

**Practical comparison:** ~400x slower than primecount at 10^14 despite better asymptotics.
Crossover point likely around N ~ 10^{30}+, computationally unreachable for either method.

**Does it change the barrier?** No. Still O(x^{1/2+eps}) = O(2^{N/2+eps}), exponential in input bits.

### 1.3 Analytic Methods (Best Theoretical)

**Current best asymptotic:** O(x^{1/2+eps}) -- Lagarias-Odlyzko (1987)

**Implementation status:**
- Platt's method: rigorous interval arithmetic, computed pi(10^25) unconditionally
- FKBJ method: Weil-Barner explicit formula with precomputed zeta zeros, 128-bit fixed-point
- In practice, 100x slower than combinatorial methods due to zeta zero computation overhead
- ~40,000 CPU core hours for pi(10^25) vs 330 for combinatorial

**Relevance to project:** Theoretically better asymptotics but impractical constants. The 100x gap persists for all currently reachable values of x.

**Does it change the barrier?** No. O(x^{1/2+eps}) has stood since 1987.

### 1.4 Aggarwal: Complexity of Computing p(n) (October 2025)

**Source:** arXiv:2510.16285

**Key results:**
- p(n) computable in O(sqrt(n) * (log n)^4) unconditionally via binary search on pi(x)
- Conditional (RH + Cramer): O(sqrt(n) * (log n)^{7/2} * log log n)
- Proves sieve-based methods CANNOT asymptotically beat binary-search-on-pi(x)
  - Best sieve for p(n): O(n log n / log log n)
  - Binary search + pi(x): O(sqrt(n) * polylog(n))
  - Binary search wins by factor ~sqrt(n)

**Relevance:** Most rigorous recent analysis. Confirms our implemented approach (Lucy_Hedgehog/Meissel-Lehmer + binary search) is provably superior to any sieve for computing a single p(n).

**Does it change the barrier?** No new algorithm, but provides the tightest known complexity analysis.

---

## 2. Riemann Hypothesis and Zero-Density Estimates

### 2.1 Guth-Maynard: Large Value Estimates (May 2024)

**Source:** arXiv:2405.20552; published in Annals of Mathematics, 2025

**The breakthrough:**
- First improvement to Ingham's 1940 zero-density bound in 80+ years
- Improved exponent from 3/5 = 0.6 to 13/25 = 0.52
- New zero-density estimate: N(sigma, T) <= T^{30(1-sigma)/13 + o(1)}

**Propagated results:**
- PNT in short intervals: [x, x + x^{17/30}], improving from x^{3/5}
- Almost all short intervals: theta improved from 1/6 to 2/15
- Tao: "a remarkable achievement"; Radziwill: "first new idea in the hunt for zeta zeros in 50 years"

**Relevance:** Does NOT directly yield a new algorithm, but tightens theoretical guarantees about prime distribution. Better zero-density bounds could reduce the number of zeta zeros needed for rigorous analytic pi(x) computation; practical impact unclear.

**Does it change the barrier?** No. Improves density estimates, not computational methods.

### 2.2 Gafni-Tao: Exceptional Intervals (May 2025)

**Source:** arXiv:2505.24017

Makes Guth-Maynard results explicit and quantitative for short intervals. Provides explicit upper bounds on exceptional set where PNT fails in [x, x+x^theta]. Computer-assisted optimization.

**Relevance:** Tighter guarantees for binary search bracketing. Marginal practical impact since Visser's exponential bounds already give tight brackets.

### 2.3 Connes: Letter Through Time (February 2026)

**Source:** arXiv:2602.04022

42-page survey plus novel "Letter to Riemann" approach:
- Using primes < 13, approximates first 50 zeta zeros to accuracies 2.6 x 10^{-55} to 10^{-3}
- Proves these approximating values lie EXACTLY on the critical line
- Connects Weil quadratic form to information theory
- Proposes proof strategy: convergence of zeros from finite to infinite Euler products

**Relevance:** Fascinating analytic tool but NOT computational. No pathway to polylog prime counting. Our Session 9 analyzed the Connes-Weil quadratic form and found: "finds zeros but can't sum them fast."

**Does it change the barrier?** No.

### 2.4 Overall RH Status (April 2026)

- 20 trillion zeros verified on critical line (Platt-Trudgian, 2021)
- De Bruijn-Newman constant Lambda = 0 (Rodgers-Tao, 2020)
- Lean 4 Mathlib: machine-checked proof of Prime Number Theorem
- RH remains OPEN after 167 years

---

## 3. Prime Distribution and Sieve Advances

### 3.1 Lichtman: Modified Linear Sieve (Published 2025)

**Source:** Algebra & Number Theory 19(1), 2025; arXiv:2109.02851

Equidistributes primes in arithmetic progressions to moduli up to x^{10/17}, surpassing the x^{4/7} (BFI classical) and x^{7/12} (Maynard) barriers. New upper bound on twin prime count -- largest improvement since BFI (1986).

**Relevance:** Fundamental sieve infrastructure. The x^{10/17} ~ x^{0.588} threshold is approaching the x^{2/3} ~ x^{0.667} used by Deleglise-Rivat but has not reached it.

**Does it change the barrier?** No, but closing the gap between equidistribution threshold and algorithmic threshold.

### 3.2 Pascadi: Exponents of Distribution to x^{5/8} (May 2025)

**Source:** arXiv:2505.00653

Primes and smooth numbers equidistributed in APs to moduli x^{5/8 - o(1)}. Eliminates dependency on Selberg eigenvalue conjecture. Uses triply-well-factorable weights and new large sieve inequalities for exceptional Maass forms.

**Relevance:** x^{5/8} = x^{0.625} is tantalizingly close to x^{2/3} = x^{0.667}. If pushed to x^{2/3}, might enable algorithmic improvements to combinatorial pi(x). Currently below threshold.

**Does it change the barrier?** Not yet, but worth monitoring.

### 3.3 Le Duc Hieu: AP of Primes in Short Intervals (September 2025)

**Source:** arXiv:2509.04883

Once theta > 17/30, every sufficiently long interval [x, x+x^theta] contains many k-term APs of primes. Combines Guth-Maynard uniform short-interval PNT with Green-Tao transference principle.

**Relevance:** Purely theoretical. Confirms richness of prime structure in computationally accessible intervals but provides no new algorithms.

### 3.4 Visser: Exponential Bounds on Primes and Gaps (2025)

**Sources:** Mathematics (MDPI) 13(11), 1844; Int. Math. Forum 2025; arXiv:2507.14410

Three papers providing:
- Exponentially tight bounds on p(n) location
- Explicit bounds on prime gaps g(n) with proven validity ranges
- Analysis of Chebyshev theta at primes

**Relevance:** Tighter bracketing for binary search when computing p(n). Reduces number of pi(x) evaluations needed.

**Does it change the barrier?** No. Constant-factor improvement to search phase, not asymptotic.

---

## 4. Exact Formulas and Characterizations

### 4.1 Ono-Craig-van Ittersum: Partition Detection of Primes (2024)

**Source:** PNAS September 2024; arXiv:2405.06451

n is prime iff certain equations in MacMahon partition functions M_k(n) hold. Infinitely many such characterizations exist.

**Follow-up (2025):**
- Structural conjecture PROVED by two independent teams (Kane-Krishnamoorthy-Lau; van Ittersum-Mauth-Ono-Singh) using analytic and l-adic methods
- Extended to detect cubes of primes and primes in arithmetic progressions
- AMS Notices feature (June 2025); Cozzarelli Prize finalist

**Relevance:** New theoretical lens. Computing M_k(n) is not competitive with sieving, but the extensions show the framework is more powerful than initially apparent.

**Does it change the barrier?** No. Structural, not computational.

### 4.2 Cloitre: Effective Analytic Recurrence (2025)

**Source:** arXiv:2508.02690

First proven effective analytic recurrence:
p(n+1) = ceil((-1 + zeta(2*p_n) * prod(1 - 1/p_j^{2*p_n}))^{-1/(2*p_n)})

Exact for all n >= 1. No sieving -- pure zeta evaluation. Connects Gandhi's arithmetic approach with Golomb's analytic framework.

**Relevance:** Theoretically significant as the first fully constructive analytic prime recurrence. Practically irrelevant: requires extreme-precision zeta evaluation.

**Does it change the barrier?** No. Trades sieving for equally expensive zeta computation.

### 4.3 Prunescu-Shunia: Arithmetic Terms (December 2024 - 2025)

**Source:** arXiv:2412.14594; J. Integer Seq. 28, 2025

First explicit fixed-length arithmetic terms for pi(n) and p(n). Uses only +, -, *, div, exp applied a fixed finite number of times.

**The catch:** Tower-of-exponentials intermediate values. Even p(1)=2 requires computing integers with more digits than atoms in the observable universe. This is provably inherent (Mazzanti 2002): any arithmetic term encoding variable-length computation into fixed operations MUST use exponentially large intermediates.

**Extended work (2025):** Arithmetic terms for Mersenne/Fermat/twin prime counting (arXiv:2512.01680); smallest/largest prime divisor (arXiv:2510.26939).

**Relevance:** Settles Hardy-Wright's question: yes, a formula for p(n) exists. But it is "utterly useless for computation" as Hardy-Wright anticipated.

**Does it change the barrier?** No. Theoretical existence result only.

### 4.4 Saito: Mills' Constant Is Irrational (2024-2025)

**Source:** arXiv:2404.19461; Mathematika 71(3), 2025

Resolves long-standing open problem. Partial transcendence results obtained. The irrationality means Mills' constant cannot be expressed as a ratio of integers, but the formula floor(A^{3^n}) remains circular (computing A requires primes).

**Does it change the barrier?** No.

### 4.5 Bouras: New Prime-Generating Recurrence (2025)

**Source:** arXiv:2509.09745

GCD-based recurrence producing only 1s and primes, improving on Rowland (2008). Represented as finite continued fraction.

**Relevance:** Does not produce ALL primes. Theoretical interest only.

**Does it change the barrier?** No.

---

## 5. AI/ML in Mathematics

### 5.1 AlphaEvolve (Google DeepMind, 2025-2026)

Applied to 50+ open problems. In 20% of cases, improved on best known solutions:
- Improved lower bounds for 5 classical Ramsey numbers (records stood 10+ years)
- Erdos minimum overlap problem: improved upper bound
- New matrix multiplication methods

**Assessment for primes:** Impressive for combinatorial optimization but NOT applicable to prime counting. The problem is exact computation, not optimization.

### 5.2 AI Solving Erdos Problems (2025-2026)

- ~100 Erdos problems moved to "solved" since October 2025 (per Tao's tracking)
- 15 solved since Christmas 2025; 11 credited AI involvement
- GPT-5.2 Pro solved Problem #397 (central binomial coefficients), verified by Tao
- Axiom AI solved 4 open problems including Chen-Gendron conjecture
- Ken Ono used AxiomProver for overnight proof derivation

**Tao's caveat:** "AI is still nowhere near being able to solve major open problems." Most AI assistance is "souped-up literature search" + theorem assembly.

**Does it change the barrier?** No. No AI has produced a new algorithm or formula for primes.

### 5.3 Gauss Agent: Formalization (Math, Inc., 2026)

Autoformalization agent for Lean proof verification. Formalized Strong Prime Number Theorem in 3 weeks (~25,000 lines of Lean, 1,000+ theorems). Previous human effort: 18+ months.

**Assessment:** Important for verification infrastructure. Formalizing PNT does not yield faster prime algorithms.

### 5.4 ML for Prime Classification

The withdrawn paper (arXiv:2402.03363, 2024) was the most notable attempt -- ResNet + Transformer achieving ~99% recall, never 100%. Authors retracted it, stating the approach "does not lead to meaningful or valid conclusions."

**Theoretical limitation:** Prime Coding Theorem (arXiv:2308.10817, 2023) proves information content of prime sequence exceeds any learnable pattern.

**Does it change the barrier?** No. ML cannot achieve exactness for primes.

---

## 6. Practical Implementations

### 6.1 Our Best Implementation

**v7 Optimized Hybrid (Session 3):**
- R^{-1} + Lucy_Hedgehog DP + Newton/bisection + Miller-Rabin
- p(10^7) in 1.28s, p(10^8) in 7.73s
- 10,000/10,000 verified correct

### 6.2 Theoretical Best for p(10^100)

p(10^100) ~ 2.35 x 10^102. Computing it exactly requires:
- pi(x) at x ~ 10^102, costing O(x^{2/3}) = O(10^68) minimum operations
- Even analytic: O(x^{1/2+eps}) = O(10^{51+eps}) operations
- At 10^15 ops/sec: 10^36 seconds (billion billion billion years)
- Quantum: best O(x^{1/3}) = 10^34 quantum ops -> 10^25 seconds

---

## Updated Assessment (April 2026)

### What has changed since the project began:

1. **Aggarwal (2025)** provided rigorous proof that binary-search-on-pi(x) is provably superior to sieving for individual p(n). Validates our approach.

2. **Guth-Maynard (2024)** broke an 80-year barrier in zero-density estimates. The most significant analytic number theory result in decades, but it improves distribution knowledge, not algorithms.

3. **Prunescu-Shunia (2024-2025)** proved a formula for p(n) EXISTS as a fixed-length expression. But tower-exponential intermediates make it computationally useless.

4. **Sieve equidistribution** (Lichtman x^{10/17}, Pascadi x^{5/8}) is creeping toward the x^{2/3} algorithmic threshold. Worth watching but no immediate impact.

5. **AI tools** are accelerating mathematical research but have produced no new algorithms or complexity results for primes.

6. **primecount** continues constant-factor improvements (SIMD, AVX512) but no asymptotic changes.

### What has NOT changed:

- **No new algorithm** for pi(x) or p(n) with better asymptotic complexity
- **No unconditional lower bound** beyond Omega(log x) for computing pi(x)
- **No polylogarithmic method** for pi(x) in any computational model
- **The O(x^{2/3}) combinatorial barrier** and **O(x^{1/2+eps}) analytic barrier** are intact
- **Quantum computing** still cannot break the O(x^{1/3}) barrier (Nayak-Wu lower bound)
- **The problem is genuinely OPEN** -- no proof exists that polylog is impossible

7. **Ono-Craig-van Ittersum (2024)** proved primes satisfy infinitely many Diophantine equations in partition functions. Beautiful but computationally worse (requires divisor computation). No algorithmic impact.

8. **Tao-Gafni (2025)** confirmed Erdős conjecture on rough numbers in prime gaps. Gap structure research, not counting.

9. **Chen-Tal-Wang (ECCC 2026)** proved n^{2.5-ε} lower bounds for depth-2 threshold circuits. Strongest THR∘THR lower bound, but hard function is in E^NP, not number-theoretic.

10. **Session 17** established EXACT communication complexity: rank(pi_N) = 2^{N/2-1}+2 for balanced bit partition. The √x barrier is universal across all algebraic measures.

### Bottom line:

The problem of computing p(n) in polylogarithmic time remains open in the strongest theoretical sense (no impossibility proof), but all 480+ approaches across 17 sessions confirm the barrier empirically. The 2024-2026 literature -- including the most significant results in decades (Guth-Maynard, Prunescu-Shunia, Aggarwal, Ono) -- provides deeper understanding without changing the computational landscape. Session 17's exact communication complexity formula (rank = 2^{N/2-1}+2) gives the most precise characterization yet of the √x barrier. The gap between the trivial Omega(log x) lower bound and the O(x^{1/2+eps}) upper bound for pi(x) remains one of the least-explored frontiers in computational complexity.

---

Generated: 2026-04-04
Consolidated from: latest_2026_breakthroughs.md, latest_research_2026.md, exact_formulas_2026.md, exact_formulas_research.md, april_2026_new_findings.md, session10_internet_search.md
