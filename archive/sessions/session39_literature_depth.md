# Session 39: Literature Monitoring + Threshold Depth Profile

**Date:** 2026-04-05
**Focus:** Literature monitoring for TC^0/NC^1 results; threshold circuit depth experiment
**Sub-agents used:** 2 (codebase exploration, duplicate script search)

---

## Summary

This session focused on two goals: (1) monitoring for new 2026 publications
that could change the project's barrier landscape, and (2) running a new
experiment measuring the threshold circuit depth profile of primality functions.

**Bottom line:** No new results change the picture. The problem remains open
but all experimental directions are exhausted. Future work is limited to
literature monitoring and engineering improvements.

---

## 1. Literature Monitoring

### ECCC 2026 papers reviewed:
- **TR26-002** (Behera, Hansen, Limaye, Srinivasan): Multilinear IPS separations.
  Separates multilinear NC^1 from NC^2 in proof complexity. Not relevant to circuits.
- **TR26-008**: Circuit lower bounds for linear functions. Not relevant.
- **TR26-014**: Switching lemmas over finite fields. Tangentially related but no impact.
- **TR26-024**: Hilbert's Nullstellensatz in counting hierarchy. No computational relevance.
- **TR26-025**: Tensor circuit complexity. No direct connection.
- **TR26-033**: Not reviewed in detail.
- **TR26-038**: XOR lemma for F_p sums. Graph counting, not primes.
- **TR26-039** (Chen-Tal-Wang): n^{2.5-ε} THR∘THR lower bounds for E^NP functions.
  Already tracked. Does not apply to number-theoretic functions.

### arXiv 2026 papers reviewed:
- **2601.04072** (Gurumukhani et al.): Optimal monotone depth-3 lower bounds for MAJORITY.
  Uses local enumeration technique. Applies to MONOTONE circuits, not threshold circuits.
  MAJORITY is trivially in TC^0 (single gate), so this is about a different model.

### Brandt MKtP follow-ups:
- Liu-Pass 2025 connected MKtP to one-way functions. No 2026 follow-up found that
  extends the diagonalization technique to relevant complexity classes.

### Verdict:
No 2026 paper proposes or achieves progress on TC^0/NC^1 separation, prime counting
algorithms, or growing-dimension matrix powering. The frontier has not moved.

---

## 2. Threshold Depth Profile Experiment

**Script:** `experiments/circuit_complexity/threshold_depth_profile.py`

### Method:
For N = 4, 6, 8, 10 and depths 1-5:
- Generate k random LTFs per hidden layer (k = N, 2N, N^2, 2N^2)
- Find optimal final-layer LTF via LP (margin maximization)
- Report whether exact computation is achieved

### Results:

| N | is_prime | pi_mod2 | random |
|---|----------|---------|--------|
| 4 | depth=2, k=8 | depth=2, k=16 | depth=2, k=8 |
| 6 | depth=2, k=36 | depth=2, k=36 | depth=2, k=36 |
| 8 | depth=2, k=128 | depth=2, k=128 | depth=2, k=128 |
| 10 | >5 (0.832) | >5 (0.541) | >5 (0.539) |

At N=10 with k=1024 and 100 trials: still no exact depth-2 circuit found.
But S35 PTF results show depth-2 IS possible with 638 STRUCTURED gates (monomials).

### Interpretation:
- All three functions behave identically under random threshold features
- Random LTFs don't correlate with primality (extends pseudorandomness evidence)
- Depth alone (without structured features) doesn't help
- The experiment cannot distinguish TC^0 from non-TC^0

---

## 3. Growing-Dimension Matrix Powering Analysis

Analyzed the precise barrier for the one remaining open sub-question:

**The question:** Can M^n mod m be computed in TC^0 when M is r×r with r = O(polylog(n))?

**Why it's hard:**
1. Fixed-r MPOW IS in TC^0 (Mereghetti-Palano 2000): reduce to r scalar powerings
2. Growing-r MPOW requires combining polylog eigenvalue powers at depth O(log r) = O(log log n)
3. O(log log n) depth is between TC^0 (constant) and NC^1 (log n) — undefined complexity class

**Why known techniques fail for the AKS case:**
- HAB scalar powering uses CRT + discrete log mod small primes → requires factoring the modulus (circular for primality)
- Healy-Viola F_{2^n} uses Frobenius endomorphism (x→x^p) → no analog over Z_n (unknown characteristic)
- Commutative ring structure: Z_n[x]/(x^r-1) decomposes via CRT into cyclotomic factors, but decomposition depends on n's factorization (circular)

**Positive observations:**
- Powers of a single matrix generate a CYCLIC (hence solvable) subgroup — avoids the IMP barrier
- The BPSW path bypasses matrix powering entirely (only needs scalar pow + 2×2 MPOW)
- But BPSW correctness is unproven and appears very hard to resolve

---

## 4. Housekeeping

Flagged 7 groups of duplicate/versioned scripts in TODO.md:
- weil_optimized (3 versions), ht_transfer (2), k_party_nof (2),
  approx_degree (5 variants), information_shortcut (2),
  breakthrough_attempt_v2 (no v1), self_correcting_v2 (no v1)
- All have companion _results.md files. Human review needed before cleanup.

---

## 5. Assessment

**The project has reached a natural plateau.** The only remaining open direction
(circuit complexity of pi(x)) is a major open problem in complexity theory
(#TC^0 ⊆ NC?) with no known experimental attack path. The 2026 literature
confirms no progress on this front from the research community either.

**Productive future work:**
- Literature monitoring (quarterly): new ECCC/CCC/STOC/FOCS papers on TC^0/NC^1
- Engineering: integrate primecount v8.4 SIMD benchmarks when available
- Theoretical: tighten existing novel results (extend pseudorandomness measures)
- Cleanup: resolve duplicate scripts (human review needed)

**What would change the picture:**
- A TC^0 ≠ NC^1 proof (would likely close the problem)
- A new primality test in TC^0 that's unconditionally correct (bypasses BPSW)
- An unconditional lower bound for pi(x) beyond Omega(log x)
- A breakthrough in zeta zero summation methods
