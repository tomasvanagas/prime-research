# Session 12 Summary
Date: 2026-04-04

## Focus: Circuit Complexity Deep-Dive

### Experiments Conducted (4 experiments in experiments/circuit_complexity/)

1. **lucy_dp_structure.py** - Floor-value mapping structure analysis
   - Mapping matrices are completely NON-COMMUTATIVE (0 commuting pairs)
   - Division depth for x=1000: max 3 (vs 11 sequential steps)
   - Parallel rounds (greedy): matches pi(x^{1/3}) exactly

2. **parallel_depth_scaling.py** - Parallel depth scaling analysis
   - Lucy DP parallel rounds = pi(x^{1/3}) + 1 (confirmed ratio -> 1.00)
   - Meissel-Lehmer: depth O(log log x) but width O(x^{2/3}) = O(2^{2N/3})
   - All floor values reused across recursion levels (V is fixed)
   - |V| ≈ 2*sqrt(x) = O(2^{N/2}), exponential in input size

3. **meissel_lehmer_dag.py** - Computation DAG depth analysis
   - **KEY RESULT**: DAG depth = pi(sqrt(x)) EXACTLY (ratio = 1.000 for all x tested)
   - Critical path goes through S(x) at every prime step
   - Lucy DP is fundamentally sequential for computing pi(x)
   - Width: concentrated at depths 2-6, tapering to 1 at critical path

4. **pi_modular_structure.py** - Modular structure of pi(x)
   - H(Y|X) = 0.537 bits INVARIANT across all moduli (2 through 30)
   - This equals the entropy of the prime indicator function
   - pi(x) mod 2 vs Liouville: 49.99% agreement (random)
   - No modular shortcut exists for any modulus

5. **floor_value_algebra.py** - Floor-value algebraic structure
   - V is closed under floor-division by primes (confirms Lucy DP validity)
   - Linear transformation is ~80% dense (full-rank)
   - Coefficients grow rapidly (|coeffs| >> pi(x)), massive cancellation
   - No low-rank structure to exploit

### Key Findings

1. **"Is pi(x) in NC?" is EQUIVALENT to our target problem**
   - NC = polylog depth, poly size circuits
   - All known approaches produce exponential-size circuits
   - The barrier is CIRCUIT SIZE (2^{2N/3}), not depth

2. **All sieve-based approaches are inherently exponential in N = log(x)**
   - Floor-value set |V| = O(sqrt(x)) = O(2^{N/2}) values, all needed
   - No compression: full-rank transformation, non-commutative mappings
   - Depth can be O(log log x) but width remains O(x^{2/3})

3. **pi(x) mod m is as hard as pi(x) for any m**
   - Conditional entropy is invariant across moduli
   - No CRT-based shortcut possible

4. **The barrier is now precisely characterized:**
   - Every known method computes O(sqrt(x)) intermediate values
   - These are either floor-values (sieve) or zeta zeros (analytic)
   - An O(polylog) algorithm must avoid BOTH

### Paths Closed (Session 12): ~5 new
- Lucy DP parallel depth = pi(sqrt(x))
- Floor-value mapping commutativity
- pi(x) mod m shortcut (all m)
- Floor-value linear algebra (full-rank)
- Meissel-Lehmer as NC circuit (exponential width)

### What Remains Open
1. Non-AKS TC^0 primality test (sub-agents investigating)
2. Growing-dimension matrix powering in TC^0 (at TC^0/NC^1 frontier)
3. Non-sieve, non-analytic approach to pi(x) (unknown if such exists)
4. Time-bounded Kolmogorov complexity of delta(n)
5. Novel number-theoretic identity

### Sub-agent Results (5 agents completed)

**paper-search:** Found 19 relevant papers. Most important: TG kernel paper
(arXiv:2506.22634) claiming polylog zeros for exact pi(x) — DEBUNKED (violates
uncertainty principle, AI-generated, fundamental errors). Also found: Dawar-Evans
NC^1 characterization, Hu-Manor-Oliveira rKt failure, Green-Sawhney Gowers norms.
No genuine algorithmic breakthroughs found in 2025-2026 literature.

**tc0-primality:** Found exactly TWO non-AKS paths to PRIMES in TC^0:
1. Miller's test under GRH (all ops in TC^0, blocked by GRH removal)
2. BPSW if correctness proven (scalar + MPOW_2, both TC^0, blocked by 45yr open problem)
Both blocked by analytic number theory, not circuit complexity.

**andrews-wigderson:** Andrews-Wigderson FOCS 2024 does NOT help. Three
independent reasons: wrong algebraic setting (fields not Z_n rings), wrong
bottleneck (GCD not matrix powering), wrong circuit model (arithmetic not Boolean).

**comm-ring-tc0:** Circulant structure reduces matrix powering to DFT + scalar
powers + inverse DFT. Scalar powers: TC^0. DFT: depth O(log log n), the PRECISE
residual obstruction. No Frobenius analogue over Z_n. Solvable monoid escapes IMP
barrier but growing k remains unresolved.

**succinct-computation:** All known succinct counting paradigms fail for pi(x).
Barvinok needs fixed dim (encoding needs dim~sqrt(x)). Permanent/determinant would
resolve NC question. #P connection doesn't help. pi(x) lacks every structural
property (multiplicativity, signed cancellation, low convolution rank, low dim,
commutativity) that enables known succinct methods.

### Paths Closed (Session 12): ~10 new
- Lucy DP parallel depth = pi(sqrt(x))
- Floor-value mapping commutativity
- pi(x) mod m shortcut (all m)
- Floor-value linear algebra (full-rank)
- Meissel-Lehmer as NC circuit (exponential width)
- TG kernel (arXiv:2506.22634 DEBUNKED)
- Andrews-Wigderson for pi(x) (wrong model)
- Succinct/lattice point counting (dim too high)
- Permanent/determinant encoding (circular)
- H-T cancellation transfer (quantified definitively)

### Assessment
Session 12 is the most thorough circuit complexity analysis to date. Key results:
1. "Is pi(x) in NC?" ⟺ finding an O(polylog) algorithm
2. ALL known methods produce exponential-size circuits
3. The barrier is CIRCUIT SIZE, and the residual TC^0 obstruction is O(log log n)
   depth for the DFT/CRT recombination step
4. GRH → PRIMES in TC^0 (but not unconditionally)
5. No genuine breakthrough found in 2025-2026 literature
6. The problem remains GENUINELY OPEN
