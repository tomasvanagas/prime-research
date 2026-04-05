# Critique of Session 33 Proposals
## April 5, 2026

**Reviewer:** Critique agent (Session 34)
**Source:** `archive/ephemeral/proposals_session.md` (4 proposals, numbered 16–20)

---

## Proposal 16: Spectral Truncation with Adaptive Zero Selection

### Verdict: DUPLICATE

**CLOSED_PATHS matches:**
- "Explicit formula + few zeros" (S5): K^{−0.01} convergence — established that few zeros don't suffice
- "Convergence acceleration of zero sum" (S32): Richardson, Shanks, etc. — errors GROW as N^{0.8–1.0}
- "Hybrid analytic + local sieve" (S32): Optimal hybrid = O(x^{1/2+ε}), matching Lagarias-Odlyzko
- "GRH explicit formula optimal T" (S33): Under GRH, need T = O(√x · log²x) zeros for error < 0.5
- "Riemann-Siegel analog" (S6): = Lagarias-Odlyzko method
- "Zeta zero compression/FMM" (S9): Spectral flatness 0.91, phase catastrophe

**Failure mode:** EQUIVALENCE. Selecting "resonant" zeros is just reordering the explicit formula summation. The error after K zeros is O(x^{1/2}/K) regardless of selection order. The proposal rediscovers the Lagarias-Odlyzko bound from a different angle.

**Specific mathematical obstacle:** The explicit formula error term after K zeros is:
  |π(x) − R(x) − Σ_{k=1}^{K} R(x^{ρ_k})| = O(x^{1/2} · T^{−1} · log²T)
where T ~ γ_K. For error < 0.5, need T ~ √x, i.e., K ~ √x · log(√x)/(2π). No reordering changes the total information content needed.

**What the experiment actually showed:** For n ≥ 500, 1000 zeros insufficient. This is consistent with 6+ prior experiments. The "1–2 zeros suffice for n < 100" observation is a known small-number artifact (R^{−1}(n) is already within 1 of p(n) for small n).

**Novel content:** None. The "adaptive selection" framing is cosmetically different from "convergence acceleration" (S32) but mathematically identical — both attempt to reorder/select zeros to reduce the count below O(√x).

---

## Proposal 17: p-adic Lifting via CRT

### Verdict: DUPLICATE

**CLOSED_PATHS matches (extensive):**
- "CRT/modular construction" (S3,7,9): Fundamentally circular
- "Modular CRT for pi(x)" (S4,5,7): pi(x) mod m costs same as pi(x)
- "CRT reconstruction of pi(x)" (S13): Each pi(x) mod q costs O(x^{2/3}); k moduli → k× worse
- "Prime race pi(x;q,a) shortcut for CRT" (S22): L-function zeros same barrier
- "CRT Prime Locator" (S24): Only 4–5 CRT moduli needed but each requires full prime counting
- "CRT from modular periods of p(n)" (S29): p(n) mod m NOT periodic for m ≥ 5
- "Adelic local-global reconstruction" (S24): 2 moduli suffice but each needs L-function zeros
- "p-adic interpolation of pi(x)" (S29): pi(x) NOT p-adically continuous
- "CRT novel encoding (178 bits)" (S7): Identifies info but no way to obtain it
- "CRT modular π(x) mod m reconstruction" (S20): random walk, no structure
- "Batch modular exp via CRT" (S16): Exponent-modulus coupling prevents batch
- "Class field tower prime reconstruction" (S25): Circular + more zeros needed

This is the **single most tested approach** in the entire project, with 12+ closed entries spanning sessions 3–29. The proposal's "novel insight" (information-theoretic cost is low but computational cost is high) was explicitly noted in S7 ("CRT novel encoding: identifies info but no way to obtain it") and is the core finding of `novel/info_computation_gap.md`.

**Failure mode:** CIRCULARITY + EQUIVALENCE. floor(x/d) mod q ≠ floor((x mod q)/(d mod q)). The Legendre recursion cannot be reduced modulo q because floor division is not a ring homomorphism. Computing π(x) mod p requires the full Meissel-Lehmer computation.

**Novel content:** None. The floor-division-is-not-a-homomorphism observation was established in S13.

---

## Proposal 18: Dequantized Grover Counting

### Verdict: PARTIALLY NOVEL (framing only)

**CLOSED_PATHS matches:**
- "Grover on Deleglise-Rivat sieve" (S7): O(x^{1/3}) quantum, still 10^34 ops
- "Entanglement-assisted counting" (S10): No advantage over classical
- "Batched primality for pi(x)" (S14): No shared structure in modular exponentiation

**Failure mode:** INFORMATION LOSS. The dequantization framework (Tang 2024, Chia et al. 2025) applies when the input has low-rank structure. The Möbius function μ(d) and the prime indicator are full-rank (confirmed: ANF degree = Θ(N), approximate degree = N/2, communication rank = 2^{N/2−1}+2). Dequantization cannot help for full-rank problems.

**What's partially novel:** The *framing* of applying dequantization results to prime counting is new — no prior CLOSED_PATHS entry specifically references Tang 2024 or Chia et al. 2025 in connection with π(x). However, the underlying reason for failure (full-rank / no low-rank shortcut) has been established from multiple angles (S13 GF(2) analysis, S17 communication rank, S23 approximate degree, S28 tensor rank).

**Specific mathematical obstacle:** For dequantization to give a speedup, one needs the input matrix (encoding the sieve/counting problem) to have rank r ≪ N. The prime indicator's approximate degree is N/2 and its communication matrix rank is 2^{N/2−1}+2 — both are maximal up to constant factors. There is no low-rank structure to exploit.

**Experiment assessment:** The experiment correctly identifies the obstacle but the test (truncated Möbius sums, sampling-based estimates) doesn't actually implement dequantization — it tests simpler proxies. A proper test would construct the quantum circuit for Grover counting over the sieve, then apply Tang's sample-and-query framework. This is infeasible for interesting sizes but the conclusion is correct regardless: full-rank blocks dequantization.

---

## Proposal 19: Ramanujan Library / PSLQ on delta(n)

### Verdict: DUPLICATE

**CLOSED_PATHS matches:**
- "Symbolic regression/PSLQ" (S5): No pattern in delta(n)
- "PSLQ/LLL exhaustive identity search" (S18,19): All 6 relation types tested, all spurious
- "Kt complexity of delta(n) empirical" (S19): AR(1) R²=0.996 but RMSE=10.5
- "Cipolla residual autoregression" (S29): AR(1) 91.4% reduction but irreducible O(log n) error
- "Kolmogorov complexity of δ(n)" (S24): 5× more compressible than random but insufficient for exactness
- "CF / Stern-Brocot structure of p(n)/n" (S29): CF partial quotients follow Gauss-Kuzmin (random)
- "Gap-based predictability" (S19): AR(1..50) gives NO improvement over baseline

**Failure mode:** INFORMATION LOSS. δ(n) is determined by O(√p(n)) zeta zeros with GUE-random phases. Any finite basis of computable functions misses the oscillatory component. PSLQ finds different relations per n because there IS no universal formula.

**What the experiment confirms:** Strong lag-1 autocorrelation (0.95) — known from S19 (AR(1) R²=0.996). No modular patterns — confirmed S20, S29. PSLQ finds different relations per n — confirmed S18, S19. RMSE = 31 from basis fit — comparable to S19's RMSE = 10.5 with AR(1). The "Ramanujan Library" inspiration is new branding but the method (PSLQ/LLL on mathematical constants) is exactly what was done in S18/S19.

**Novel content:** None. The ICLR 2025 Ramanujan Library paper is for discovering *new* identities among constants, not for finding structure in pseudorandom sequences. Applying it to δ(n) was already tried (S18/S19 used PSLQ which is the same core method).

---

## Proposal 20: Étale Cohomology Point-Counting

### Verdict: DUPLICATE

**CLOSED_PATHS matches:**
- "Etale cohomology of Spec(Z)" (S7): Recovers Euler product = knowing all primes
- "Algebraic geometry point counting for pi(x)" (S13): Genus must be Ω(x/ln x); H¹_ét(Spec(Z)) infinite-dim; Weil duality = explicit formula
- "Algebraic variety F_q point count" (S14): Low-dim: too few points. High-dim: Kedlaya O(N³ · 2^N) = sieve
- "Elliptic curve point counts" (S4): Can't invert a_p to p
- "Motivic cohomology" (S7): H¹_mot(Spec(Z)) = Z^∞ — isomorphic to primes

**Failure mode:** CIRCULARITY + EQUIVALENCE. The Grothendieck-Lefschetz trace formula counts points on a FIXED variety over F_q in polylog(q) time. But encoding "count primes up to x" as a point-counting problem requires a variety whose dimension/genus grows with x. Specifically:
- Fixed variety: polylog, but encodes the WRONG function
- Variable variety encoding π(x): correct function, but construction is circular (requires knowing primes) and dimension grows as Ω(x/ln x)

The Frobenius eigenvalues of the variety ARE the zeta zeros (for Spec(Z)), completing the circle: étale cohomology → Frobenius eigenvalues → zeta zeros → explicit formula.

**What the experiment shows:** Character sums work but cost O(x). Hasse-Weil provides enough information bits but each costs O(x^{2/3}). Elliptic curve traces have correlation 0.39 with π(x) — too weak. All consistent with prior findings (S7, S13, S14).

**Novel content:** None.

---

## Summary

| # | Proposal | Verdict | CLOSED_PATHS matches | Failure mode |
|---|----------|---------|---------------------|--------------|
| 16 | Spectral Truncation | **DUPLICATE** | S5, S6, S9, S32, S33 (6+) | Equivalence |
| 17 | p-adic Lifting / CRT | **DUPLICATE** | S3–S29 (12+) | Circularity + Equivalence |
| 18 | Dequantized Grover | **PARTIALLY NOVEL** (framing) | S7, S10, S14 | Information Loss |
| 19 | PSLQ on δ(n) | **DUPLICATE** | S5, S18, S19, S24, S29 (7+) | Information Loss |
| 20 | Étale Cohomology | **DUPLICATE** | S4, S7, S13, S14 (5+) | Circularity + Equivalence |

### Cross-cutting assessment

The proposals' "cross-cutting insight" (all hit the same barrier: O(log n) bits encoded in O(x^α) independent pieces) is correct but was already articulated in `novel/info_computation_gap.md` and confirmed across 30+ sessions. No new barrier understanding emerged.

**No experiments warranted for NOVEL testing** — all core hypotheses have been tested in prior sessions with matching conclusions. Proposal 18's dequantization framing is the only genuinely new angle, but since it fails for the same reason (full-rank), there is nothing to test that hasn't been tested.

### Recommendations

1. **Update CLOSED_PATHS.md** with the five new entries from Session 33 (proposals 16–20)
2. **Do not add to OPEN_PROBLEMS.md** — no new viable direction found
3. **Dequantization literature** (Tang 2024, Chia et al. 2025) should be noted in `literature/references.md` as surveyed and found inapplicable due to full-rank barrier
4. The high autocorrelation of δ(n) at lag 1 (0.95) is well-known but could be cross-referenced with S19's AR(1) finding for completeness
