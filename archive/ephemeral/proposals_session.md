# Proposal Session: Fresh Approaches to O(polylog) p(n)
## Session 29 — April 5, 2026

**Methodology:** 7 proposals developed, each with runnable Python code + computational tests on n < 10000.
**Internet search:** 3 parallel agents searched arxiv, Google Scholar for 2024-2026 papers.

---

## NEW LITERATURE FINDING (URGENT)

**TG Kernel Method (ArXiv 2506.22634, June 2025)**
- Authors: Kilictas & Alpay
- Claims: Using a truncated Gaussian kernel in the explicit formula, only ~1200 zeta zeros suffice for x with 10^8 decimal digits. Error rigorously bounded below 1/2.
- If true: pi(x) via O(polylog) zeros with O(polylog) per-zero computation = O(polylog) total.
- Status: Preprint, unverified. Critical question: how does zero count scale with x?
- **ACTION NEEDED: Verify this paper's claims immediately.**

**Ono Partition-Prime Detection (PNAS, Sept 2024, ArXiv 2405.06451)**
- Authors: Craig, van Ittersum, Ono
- n >= 2 is prime iff certain polynomial equations in MacMahon partition functions M_k(n) = 0
- Example: (3n^3-13n^2+18n-8)M_1(n) + (12n^2-120n+212)M_2(n) - 960*M_3(n) = 0
- Computing M_k(n) requires O(n) work — not faster than sieving
- But the quasi-modular forms connection (Eisenstein series) is deep and might inspire analytic shortcuts
- Follow-ups: ArXiv 2412.19180, 2511.04030 (higher-level modular forms)

**Pair Correlation Advances (Goldston et al. 2024-2025)**
- PCC implies 100% of zeta zeros simple and on critical line (ArXiv 2503.15449)
- Unconditionally: >= 61.7% of zeros are simple (ArXiv 2306.04799)

Also found: Lamzouri effective LI conjecture (2311.04860), Brandt MKtP lower bounds (2024).
**No dequantization-to-number-theory bridge exists yet.** Gap in literature.
**No compressed-sensing approach to prime counting exists yet.** Another gap worth exploring.

---

## PROPOSAL 1: Sparse Fourier Transform on the Explicit Formula

### Mathematical Idea
The Riemann explicit formula pi(x) = R(x) - sum_rho R(x^rho) is a sum over zeta zero contributions. Each contribution oscillates with frequency gamma_k * ln(x). If this signal is "effectively sparse," the Sparse Fourier Transform (Hassanieh et al. 2012) can recover dominant components in O(k * polylog(N)) time.

### Pseudocode
```python
def sparse_fourier_pi(x, num_zeros=200):
    R_val = R_function(x)         # O(polylog)
    correction = 0
    for gamma in top_k_zeros(x):  # k dominant zeros
        correction -= zero_contribution(x, gamma)
    return round(R_val - correction)
```

### Complexity
O(k * polylog(x)) where k = effective sparsity.

### Experimental Results
- **Energy distribution**: 90% of energy in top ~40-55 of 200 zeros (NOT sparse)
- **Convergence**: 200 zeros insufficient for n=5000 (error 6.7)
- For small n (< 500), even 1 zero suffices
- **VERDICT: Sparsity NOT observed.** GUE statistics mean all zeros contribute comparably.

### Code: `experiments/proposals/proposal9_sparse_fourier_explicit.py`

---

## PROPOSAL 2: CRT Modular Reconstruction of pi(x)

### Mathematical Idea
Compute pi(x) mod m for many small moduli m, then reconstruct via CRT. Key sub-observation: floor(x/d) mod m has only m distinct values out of sqrt(x) terms.

### Key Finding: Floor Value Collapse
```
x=100000, mod 7: only 7 distinct values out of 316 terms (ratio 0.022)
x=1000000, mod 7: only 7 distinct values out of 1000 terms (ratio 0.007)
```

### Experimental Results
- CRT works: 7 prime moduli suffice for n=5000
- Modular sieve is correct for small n (verified up to n=500)
- **VERDICT: Correct but not faster.** The Legendre sieve has 2^pi(sqrt(x)) terms regardless of modulus. The modular reduction doesn't bypass sieve combinatorics.

### Complexity: O(x^{2/3} * log(x)/loglog(x)) — worse than direct.
### Code: `experiments/proposals/proposal10_crt_modular_pi.py`

---

## PROPOSAL 3: Compressed Sensing on Delta(n) = p(n) - R^{-1}(n)

### Mathematical Idea
delta(n) has only O(log n) bits but costs O(x^{2/3}) to compute. Analyze its compressibility in FFT/Haar bases.

### Key Experimental Findings

**Structure discovered in delta(n):**
| Property | Value |
|----------|-------|
| Autocorrelation at lag 1 | r = 0.79 |
| FFT sparsity (90% energy) | 53-72 / 256 (~25%) |
| Haar sparsity (90% energy) | 50-57 / 255 (~22%) |
| Oscillatory period | ~10-13 |

**Binary search with approximate oracle (R(x) guiding):**
| n | Steps | Oracle accuracy |
|---|-------|-----------------|
| 100 | 9 | 89% |
| 500 | 11 | 100% |
| 1000 | 12 | 100% |
| 2000 | 13 | 85% |
| 5000 | 15 | 73% |

### Novel Observation
**Autocorrelation r(1) = 0.79 is striking.** Consecutive deltas share most of their information. This is because consecutive primes see similar zeta-zero landscapes. An incremental approach exploiting this structure might achieve sublinear amortized cost.

### Complexity: Not polylog — sparsity is ~25%, not polylog.
### Code: `experiments/proposals/proposal11_delta_compressed_sensing.py`

---

## PROPOSAL 4: Recursive Interval Refinement

### Mathematical Idea
1. R^{-1}(n) locates p(n) within O(sqrt(x)/log(x))
2. Computing pi(lo) at interval start gives rank
3. pi(lo) is a SMALLER instance of the same problem
4. RECURSE

### Key Experimental Finding: Dramatic Recursive Reduction
```
n=1000:  7919 → 9.9 (one level!)
n=10000: 104729 → 28 → 1.6 (two levels!)
```

### Why It Doesn't Help
Each level still requires exact pi() computation. Level 0: O(x^{2/3}). Level 1: O(x^{1/3}). Total dominated by level 0.

### Asymptotic Expansion Analysis
| n | p(n) | R^{-1} error | 6-term expansion error |
|---|------|-------------|----------------------|
| 100 | 541 | 4.7 | 26.3 |
| 1000 | 7919 | 3.4 | 53.4 |
| 5000 | 48611 | 56.2 | 135.6 |
| 10000 | 104729 | 38.6 | 29.8 |

### Complexity: O(x^{2/3}) — circularity barrier.
### Code: `experiments/proposals/proposal12_recursive_interval_refinement.py`

---

## PROPOSAL 5: Adelic Interpolation / Multi-Residue Collapse

### Mathematical Idea
Combine real (R^{-1}) and p-adic (mod q) information. Use Lemke Oliver-Soundararajan bias in prime residue transitions to predict p(n) mod q.

### Key Experimental Findings

**Markov structure in primes mod 3:**
```
     →1    →2
1:  0.382  0.618  (strong anti-repetition!)
2:  0.613  0.386
```
Diagonal ≈ 0.38 vs expected 0.50 for uniform. Bias confirmed for all moduli.

**CRT reconstruction efficiency:**
| n | Moduli needed | Product |
|---|--------------|---------|
| 10 | 3 | 30 |
| 100 | 5 | 2310 |
| 5000 | 7 | 510510 |

**Prediction accuracy:** Bias too weak to predict individual residues (42% vs 50% baseline for mod 3).

### Complexity: Cannot compute p(n) mod q without knowing p(n). Circular.
### Code: `experiments/proposals/proposal13_adelic_interpolation.py`

---

## PROPOSAL 6: Galois Cohomology / Schoof Analogue

### Mathematical Idea
Schoof computes |E(F_p)| mod l via Frobenius action on l-torsion. By analogy, compute pi(x) mod q using Dirichlet characters.

### Key Experimental Findings

**Character sum scaling (GRH verification):**
| x | q | |S_nonprincipal| / (sqrt(x)*ln(x)) |
|---|---|--------------------------------------|
| 100 | 3 | 0.043 |
| 1000 | 3 | 0.032 |
| 10000 | 3 | 0.007 |

Ratio is DECREASING — non-principal sums may be O(sqrt(x)) or even smaller.

### The Deep Problem
S(chi_0) = pi(x) - (primes dividing q), so computing pi(x) mod q from characters still requires pi(x). Schoof works because E[l] is finite-dimensional. For primes, the "torsion" involves infinitely many L-function zeros.

### Complexity: Cannot break O(x^{2/3}).
### Code: `experiments/proposals/proposal14_galois_cohomology_count.py`

---

## PROPOSAL 7: Neural Delta Oracle + Certification

### Mathematical Idea
Train f(n) ≈ delta(n) using zeta-zero-derived features. If accurate to within 1/2, get p(n) = round(R^{-1}(n) + f(n)). Certify via primality testing.

### Experimental Results

**Linear regression (46 features):**
- Test RMSE: 11.85 (need < 0.5 for correctness)
- Round accuracy: 4% — useless

**Verification is O(polylog):**
| n | Gap | Verification ops |
|---|-----|-----------------|
| 100 | 18 | 713 |
| 1000 | 12 | 967 |
| 10000 | 6 | 802 |

**Feature importance:** Zeta zero cos/sin features rank 7th-15th (coefficients ~1-3). Dominant features are n, n*ln(ln(n)), sqrt(x).

**Gap prediction:** RMSE 7.06 (only 2.3% better than mean prediction). Gaps unpredictable.

### Key Insight
**Verification is CHEAP (O(log^4 x)). The entire difficulty is prediction.** If ANY oracle could predict delta(n) to within the prime gap, certification is polylog.

### Complexity: O(polylog) IF prediction works. It doesn't with linear models.
### Code: `experiments/proposals/proposal15_neural_delta_oracle.py`

---

## SYNTHESIS

### Universal Barriers Confirmed (Again)
1. **Information barrier**: delta(n) encodes ~log(n)/2 bits from ~sqrt(x)/log(x) zeta zeros
2. **Circularity barrier**: Most approaches reduce to computing pi(x) exactly
3. **Sparsity barrier**: Zero contributions are NOT sparse (GUE statistics)

### Top 3 Promising Directions

**1. TG Kernel paper (HIGHEST PRIORITY)**
If ~1200 zeros suffice for 10^8-digit numbers and the scaling is polylog, this IS the solution. Must read and verify the paper's proofs.

**2. Autocorrelation exploitation**
delta(n) has r(1) = 0.79. An INCREMENTAL algorithm computing delta(n+1) from delta(n) might bypass the full O(x^{2/3}) cost per query, achieving O(polylog) amortized cost for sequential queries. This doesn't give random-access O(polylog) but might give batch-mode speedup.

**3. Verification-prediction separation**
Since verification is O(polylog), the problem reduces ENTIRELY to prediction. Any new prediction oracle — quantum, ML, algebraic — that achieves delta accuracy < 0.5 would solve the problem. Focus future efforts on the prediction problem alone, not end-to-end computation.

---

## ADDITIONAL LITERATURE FINDINGS (from parallel research agents)

**Connes-Consani-Moscovici: Zeta Spectral Triples (ArXiv 2511.22755, Nov 2025)**
- Primes up to 13 reproduce first 50 zeta zeros to 10^{-55} precision via rank-one perturbation operators
- Beautiful but circular: requires primes as input to find zeros

**Brandt: MKtP Lower Bounds (TCC 2024, ePrint 2024/687)**
- Proves MKtP not in DTIME[O(n^2)] unconditionally via novel diagonalization
- **Bypasses BOTH the algebrization barrier AND the natural proofs barrier**
- Most promising route to eventually proving impossibility of polylog p(n)
- If extended to polynomial time, implies strong circuit lower bounds

**Kolpakov & Rocke: "Impossibility of discovering a formula for primes using AI" (ArXiv 2308.10817, PLOS ONE 2024)**
- Prime Coding Theorem: information content of prime sequence exceeds any learnable pattern
- Formally confirms ML cannot achieve exact prime prediction
- Closes Proposal 7 (neural delta oracle) in principle

**Lamzouri: Effective LI Conjecture (ArXiv 2311.04860, 2024)**
- Effective version of Q-linear independence of zeta zero ordinates
- If true: no finite linear combination of zero contributions vanishes → no sparse shortcut
- Formally confirms the sparsity barrier from Proposal 1

**Ono, Craig, van Ittersum: "Integer partitions detect the primes" (PNAS Sept 2024, ArXiv 2405.06451)**
- n >= 2 prime iff polynomial equations in MacMahon partition functions M_k(n) = 0
- Deep connection to quasi-modular forms (Eisenstein series)
- Computing M_k(n) requires O(n) work — not faster, but opens analytic avenues
- Follow-ups: ArXiv 2412.19180, 2511.04030, Archiv der Mathematik 2025

**Automated Conjecture Tools (2024-2026)**
- Conservative Matrix Field (PNAS 2024): structured PSLQ for continued fractions
- Ramanujan Library (ICLR 2025): 75 new constant relations via hypergraph + PSLQ
- IntSeqBERT (ArXiv 2603.05556, 2026): CRT-based integer sequence Transformer, 7.4x improvement
- Best tools for searching for delta(n) identities (Open Problem #5)
