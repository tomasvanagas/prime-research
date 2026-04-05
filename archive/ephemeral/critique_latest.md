# Critique of Session 29 Proposals (7 Proposals + 3 Synthesis Directions)

**Date:** 2026-04-05 (Critique Session 30)
**Reviewer:** Adversarial critique agent
**Source:** archive/proposals_latest.md, novel/proposals_session.md
**Prior critiques:** Session 22 (all 8 DUPLICATE), Session 27 (all 16 DUPLICATE)
**Checked against:** 626+ entries in status/CLOSED_PATHS.md (Sessions 1-29)

---

## Methodology

Each proposal evaluated against:
1. **626+ entries** in status/CLOSED_PATHS.md
2. Three proven failure modes: **Circularity (C)**, **Equivalence (E)**, **Information Loss (I)**
3. Proven barriers in proven/barriers.md
4. Communication rank formula: rank(M_pi) = 2^{N/2-1} + 2 (Session 17)
5. Approximate degree = ceil(N/2) universality (Session 28)
6. delta(n) LFSR = N/2 over all GF(p) (Sessions 24, 26)
7. **New experiments** on autocorrelation exploitation and verification-prediction separation
8. Web search for 2025-2026 literature

---

## CRITICAL ERROR IN PROPOSALS

**The TG Kernel paper (ArXiv 2506.22634) was flagged as "HIGHEST PRIORITY -- MUST VERIFY
IMMEDIATELY" without checking CLOSED_PATHS.md.** This paper was already **DEBUNKED in
Session 12**: it violates the uncertainty principle, omits an x^{1/2} factor in the error
term, the kernel is too wide to resolve individual primes, and it was identified as an
AI-generated paper with mathematical errors. See CLOSED_PATHS.md line 225.

This is a violation of Rule #3: "Search CLOSED_PATHS.md before proposing ANY approach."
The proposals session failed its most basic due diligence on what it identified as the
single most important finding.

---

## Proposal-by-Proposal Verdicts

### Proposal 1: Sparse Fourier Transform on Explicit Formula
**VERDICT: DUPLICATE**
**Failure Mode:** Information Loss (I)
**Matches:** "Zeta zero compression/FMM" (S9), "Compressed sensing on K zeros" (S7),
"Matching pursuit" (S7), "Goldreich-Levin / list decoding" (S18), "DFT spectral
structure of zeros" (S25)

The proposal correctly identifies that sparsity is NOT observed (90% energy in ~50 of
200 zeros). This was already known from Session 9 (spectral flatness 0.91) and
definitively confirmed by Session 25 (DFT matches GUE with correlation 0.9999).
Lamzouri's effective LI conjecture (2311.04860) formally confirms: if zero ordinates
are Q-linearly independent (as conjectured), no finite subset can cancel the rest.

**No novel component.**

---

### Proposal 2: CRT Modular Reconstruction of pi(x)
**VERDICT: DUPLICATE**
**Failure Mode:** Circularity (C) + Equivalence (E)
**Matches:** "CRT/modular construction" (S3,7,9), "Modular CRT for pi(x)" (S4,5,7),
"CRT reconstruction of pi(x)" (S13), "CRT modular pi(x) mod m" (S20), "CRT Prime
Locator" (S24), "Batch modular exp via CRT" (S16)

The floor-value collapse observation (7 distinct values out of 316 terms mod 7) is
mildly interesting but was implicitly covered by S13's entropy analysis (entropy ratio
~1.0 for pi(x) mod q). The proposal correctly concludes it doesn't help: sieve
combinatorics have 2^{pi(sqrt(x))} terms regardless of modulus.

**No novel component.**

---

### Proposal 3: Compressed Sensing on delta(n)
**VERDICT: PARTIALLY NOVEL then CLOSED by experiment**
**Failure Mode:** Information Loss (I)
**Matches:** "Compressed sensing on K zeros" (S7), "Prime autocorrelation mod W" (S21),
"Kt complexity of delta(n)" (S19/20)

The autocorrelation r(1) claim was the most promising finding. Our new experiment
(experiments/proposals/critique_incremental_delta.py) provides **definitive closure**:

| Metric | Value | Needed for exact |
|--------|-------|-----------------|
| Autocorrelation r(1) | 0.975 | N/A |
| AR(1) RMSE | 7.31 | < 0.5 |
| AR(5) RMSE | 7.27 | < 0.5 |
| Fraction |error| < 0.5 | 5.3% | > 99% |
| H(delta(n+1)\|delta(n)) | 4.93 bits | < 1 bit |
| RMSE scaling with n | GROWS (4.8 -> 8.1) | constant or decreasing |

**Key diagnostic:** After removing a moving average of window 10, r(1) drops from
0.975 to **0.405**. The high autocorrelation is overwhelmingly from the smooth trend
in delta(n), not from exploitable predictive structure. The "genuine memory" component
(after detrending) has RMSE still far above 0.5.

**Conditional entropy is 4.93 bits per step** -- nearly 5 bits of genuinely new
information at each step that cannot be predicted from the previous delta. This is
consistent with Session 19's finding of ~5.04 bits/prime irreducible entropy and
Session 20's incremental entropy of ~0.22*log(n) + 4.6 bits.

**CLOSED.** The autocorrelation direction fails because:
1. High r(1) is a smooth-trend artifact, not predictive structure
2. Even with full autocorrelation, RMSE >> 0.5 (by factor 14x)
3. Conditional entropy ~5 bits/step -- each step requires genuine new information
4. RMSE grows with n, so it gets WORSE at scale

---

### Proposal 4: Recursive Interval Refinement
**VERDICT: DUPLICATE**
**Failure Mode:** Circularity (C)
**Matches:** "Recursive pi(x) via pi(x/2)" (S7), "Recursive identity (8 methods)" (S9),
"Divide-and-conquer pi(x)" (S15), "Recursive halving p(2n) from p(n)" (S5)

The dramatic reduction (104729 -> 28 -> 1.6) is visually impressive but the proposal
correctly identifies the circularity: each level requires exact pi() computation at
cost O(x^{2/3}), and the total is dominated by level 0.

**No novel component.**

---

### Proposal 5: Adelic Interpolation / Multi-Residue Collapse
**VERDICT: DUPLICATE**
**Failure Mode:** Circularity (C) + Information Loss (I)
**Matches:** "CRT/modular construction" (S3,7,9), "Prime autocorrelation mod W" (S21),
"p-adic counting" (S4), "p-adic analysis" (S10), "Adelic reconstruction" (S24),
"Lemke Oliver-Soundararajan" (implicitly S21,24)

The Lemke Oliver-Soundararajan bias (diagonal ~0.38 vs expected 0.50) is real and
well-documented but provides only O(1) bits of information, not the O(log n) bits
needed. Session 21 measured: prediction accuracy 42% vs 50% baseline for mod 3.
Computing p(n) mod q requires knowing p(n) -- circular.

**No novel component.**

---

### Proposal 6: Galois Cohomology / Schoof Analogue
**VERDICT: DUPLICATE**
**Failure Mode:** Circularity (C)
**Matches:** "Chebotarev density theorem" (S4), "Character sums" (S14), "Residue class
decomposition" (S14), "Schoof analogue via characters" explicitly noted in "Galois
cohomology count" (S29 own finding)

The analogy to Schoof's algorithm is seductive but fundamentally flawed. Schoof works
because E[l] is finite-dimensional (the l-torsion of an elliptic curve is (Z/lZ)^2).
For prime counting, the "torsion" involves infinitely many L-function zeros. The
proposal correctly identifies: S(chi_0) = pi(x) itself.

The decreasing ratio of non-principal character sums is interesting but is just GRH
in action -- it doesn't help because the principal character sum IS pi(x).

**No novel component.**

---

### Proposal 7: Neural Delta Oracle + Certification
**VERDICT: DUPLICATE (with flawed synthesis claim)**
**Failure Mode:** Information Loss (I)
**Matches:** "ML correction fitting" (S3), "Deep ML" (S8), "Transformer neural" (S10),
"Neural arithmetic" (S24), "Kolpakov-Rocke impossibility" (literature)

Linear regression achieving 4% accuracy on delta(n) is consistent with all prior ML
results (0.35% - 5.4% across Sessions 3-24). The Kolpakov-Rocke Prime Coding Theorem
(PLOS ONE 2024) formally proves ML cannot achieve exactness.

**The "verification is cheap" claim is FALSE -- see Synthesis Direction #3 below.**

**No novel component.**

---

## Synthesis Direction Verdicts

### Direction 1: TG Kernel Paper
**VERDICT: ALREADY DEBUNKED (Session 12)**

ArXiv:2506.22634 was identified as an AI-generated paper with errors in Session 12. It
violates the uncertainty principle: any kernel wide enough to smooth the explicit formula
to ~1200 zero contributions necessarily has resolution worse than the prime gap. The paper
omits an x^{1/2} factor in its error analysis.

**New literature search (Session 30) reveals additional damning context:** The authors
(Kilictas & Alpay) have NO publications in number theory. Their arxiv corpus spans wildly
disparate topics ("Alpay Algebra I-V", "Latent Object Permanence in Deep Transformer
Manifolds", "Bare-Metal Tensor Virtualization", etc.). Most critically, Alpay published
a Medium article "Shaping AI's Mind from the Shadows" (July 2025) **openly admitting the
papers are strategically uploaded to arxiv to become LLM training data and influence AI
model outputs.** This was discussed on Hacker News. The paper has zero citations and no
evidence of peer review.

**This should never have been flagged as "HIGHEST PRIORITY."** The proposals session
failed to check CLOSED_PATHS.md, and the paper is not legitimate research.

---

### Direction 2: Autocorrelation Exploitation
**VERDICT: CLOSED (by new experiment)**

See Proposal 3 above. The experiment (critique_incremental_delta.py) definitively shows:
- High r(1)=0.975 is from smooth trends, not exploitable structure
- AR(1) RMSE = 7.31 >> 0.5 threshold
- Conditional entropy = 4.93 bits/step >> 1 bit
- Error GROWS with n (4.8 at small n -> 8.1 at large n)

This is now a **closed path**: delta(n) autocorrelation is NOT exploitable for
incremental exact computation.

---

### Direction 3: Verification-Prediction Separation
**VERDICT: FLAWED (circular reasoning)**

The claim "verification costs only O(polylog)" conflates two different operations:

1. **Primality testing** of a candidate g: O(log^4 x) via AKS. GENUINELY polylog.
2. **Ordinality verification** that g is the n-th prime: requires pi(g) = n, which
   costs O(x^{2/3}). NOT polylog.

Our experiment (critique_verification_separation.py) confirms:

| n | Primality cost | Ordinality cost | Ratio |
|---|---------------|----------------|-------|
| 100 | ~10^2 | ~10^2 | ~1 |
| 10000 | ~10^2 | ~10^3 | ~10 |
| 10^100 | ~10^9 | ~10^68 | ~10^59 |

The proposal says "if ANY oracle could predict delta(n) to within the prime gap,
certification is polylog." This is FALSE. Even if the oracle narrows candidates to a
single prime, you STILL need to verify it's the n-th prime overall, which requires
knowing pi(g). The "cheap verification" argument assumes you can count primes cheaply
-- which IS the original problem.

The only scenario where verification is truly O(polylog) is if you sieve a local
interval and can independently verify the count. But local sieving gives the count
within the interval, not the global rank. Computing the global rank pi(g-D) is the
O(x^{2/3}) bottleneck.

**This is now a closed path:** verification-prediction separation is circular.

---

## Literature Search (April 2026)

Web search conducted for 2025-2026 papers on prime counting algorithms, circuit
complexity of pi(x), and time-bounded Kolmogorov complexity.

**No new algorithmic breakthroughs found.** Key findings:
- TG Kernel paper (2506.22634): already debunked (Session 12); authors are NOT number
  theorists and openly admit uploading papers to influence LLM training data
- Ono partition-prime detection (PNAS 2024): O(n^2) per test, already closed (S17,23)
- Brandt MKtP lower bounds (TCC 2024): promising for impossibility proofs but not
  constructive; no follow-up linking MKtP to prime counting as of April 2026
- Kolpakov-Rocke (PLOS ONE 2024): formally closes ML approaches
- Lamzouri effective LI (2311.04860): formally confirms sparsity barrier
- Pair correlation advances (Goldston 2024-2025): theoretical, no algorithmic content
- Chen-Tal-Wang (ECCC TR26-039, March 2026): n^{2.5} THR-of-THR lower bounds, not
  number-theoretic
- New algebrization barriers (arXiv:2511.14038): confirms circuit lower bound obstacles
- CRYPTO 2025 Kt/OWF work: active area but no direct prime connections
- No new results on #TC^0 subset NC? or pi(x) circuit complexity
- "Is pi(x) in NC?" remains entirely unaddressed in 2025-2026 literature

---

## Summary

| # | Proposal | Verdict | Failure Mode | Novel? |
|---|----------|---------|-------------|--------|
| 1 | Sparse Fourier | DUPLICATE | I | No |
| 2 | CRT modular | DUPLICATE | C+E | No |
| 3 | Compressed sensing / autocorrelation | CLOSED | I | Partially (now closed) |
| 4 | Recursive refinement | DUPLICATE | C | No |
| 5 | Adelic / Markov | DUPLICATE | C+I | No |
| 6 | Schoof analogue | DUPLICATE | C | No |
| 7 | Neural delta oracle | DUPLICATE | I | No |
| S1 | TG Kernel paper | DEBUNKED (S12) | I | Already closed |
| S2 | Autocorrelation exploitation | CLOSED (new) | I | Was novel, now closed |
| S3 | Verification-prediction separation | FLAWED | C | Circular reasoning |

**0 out of 10 directions are viable.**

---

## New Closed Paths for CLOSED_PATHS.md

1. **Incremental delta(n) via autocorrelation** -- FAIL (I): r(1)=0.975 is smooth-trend artifact;
   after detrending r(1)=0.41; AR(1) RMSE=7.31 >> 0.5; conditional entropy=4.93 bits/step;
   RMSE grows with n. See experiments/proposals/critique_incremental_delta.py

2. **Verification-prediction separation** -- FAIL (C): Primality is O(polylog) but ordinality
   (is g the n-th prime?) requires pi(g) = O(x^{2/3}). Conflates primality with ordinality.
   See experiments/proposals/critique_verification_separation.py

---

## Meta-Analysis

This is now the **third consecutive critique session** (S22, S27, S30) where ALL proposals
are duplicates of previously closed paths. Session 29 generated 7 proposals with runnable
code and 3 synthesis directions; all 10 map to existing closed paths within 626+ entries.

The proposals session made one serious error (flagging the debunked TG Kernel paper as
highest priority) and one subtle error (the verification-prediction separation relies on
circular reasoning about ordinality). The autocorrelation finding (r(1)=0.975) was the
most promising lead but new experiments definitively close it.

**After 30 sessions, 628+ approaches, and 175+ sub-agents, the space of "natural"
approaches to O(polylog) p(n) appears thoroughly exhausted.** Any breakthrough must come
from genuinely novel mathematics not yet connected to prime counting -- specifically,
progress on #TC^0 subset NC? or new approaches to circuit complexity of pi(x) that avoid
the Natural Proofs barrier.

The open problems from OPEN_PROBLEMS.md remain unchanged:
1. Circuit complexity of pi(x) [most promising, genuinely unstudied]
2. Time-bounded Kolmogorov complexity of delta(n) [theoretical]
3. Zeta zero compressibility [ALL structural approaches closed]
4. Berry-Keating Hamiltonian [literature monitoring only]
