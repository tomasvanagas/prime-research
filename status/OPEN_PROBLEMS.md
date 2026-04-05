# Open Problems: Viable Research Directions

Last updated: 2026-04-05 (Session 37)

These are the ONLY directions not yet proven closed. Everything else has been
tested (525+ approaches across 37 sessions) and confirmed to hit one of three
failure modes: Circularity, Equivalence, or Information Loss.

**Key reframing (S36):** E(x) = pi(x) - Li(x) has only O(log x) bits of information,
well within polylog bounds. The barrier is COMPUTATIONAL (extracting O(log x) bits from
a sum of ~x^{1/2} cancelling terms), not information-theoretic. This reframes the problem:
find algebraic structure in the zero sum enabling bulk cancellation.

**Pseudorandomness (S37 synthesis):** pi(x) mod 2 is indistinguishable from random under
21 independent structural measures. See `novel/pseudorandomness_of_pi.md`.

---

## 1. Circuit Complexity of pi(x) [ONLY REMAINING DIRECTION]

**Question:** What is the circuit complexity of the prime-counting function?

**The gap:** Omega(log x) proven lower bound vs O(x^{1/2+epsilon}) best upper bound.
This is the largest unexplored gap in complexity theory for a natural problem.

**Equivalence:** "Is pi(x) in NC?" is EQUIVALENT to our polylog target (S12).

**What's known:**
- PRIMES is in P (AKS), NOT in AC^0 or AC^0[p] (Allender-Saks-Shparlinski 2001)
- PRIMES is in NONUNIFORM TC^0 unconditionally (S13)
- BPSW correct => PRIMES in uniform TC^0 (S13, novel: `novel/bpsw_tc0_reduction.md`)
- If BPSW in TC^0, then pi(x) in NC iff #TC^0 ⊆ NC (S15)
- AKS path to TC^0 is BLOCKED (matrix powering k>2 implies TC^0=NC^1) (S11)
- All sieve-based TC^0/NC paths closed (S11-14)
- pi(x) mod 2 is random-like under 21+ measures (S28-37, `novel/pseudorandomness_of_pi.md`)

**What's been tried and failed:**
- Sieve/floor-value methods: all produce exponential circuits (S12-14)
- TC^0 batch counting: 5 routes closed (S16)
- Meta-complexity (MKtP): reformulation, not tool (S35)
- SAT-based minimization: matches random function growth (S35)
- GF(2) algebraic structure: fully random (S35)
- Determinantal complexity: dc >= 2^{N/2-1}+2 exponential (S17)
- Communication complexity: bounded by N, cannot give super-polylog (S23)
- Novel intermediate quantities: 15 families closed (S16)

**What might still work:**
- Non-AKS TC^0 primality test using only scalar operations
- Growing-dimension matrix powering in TC^0 (genuine frontier)
- Exploiting commutativity of polynomial ring multiplication
- A fundamentally new intermediate quantity not based on floor values or zeta zeros

---

## 2. Time-Bounded Kolmogorov Complexity of delta(n) [THEORETICAL ONLY]

**Status:** Attack path CLOSED (S35), but the theoretical question remains open.

Kt framework reformulates but doesn't solve. Kt(T_N) = O(2^N*N) by sieve regardless
of circuit size. Brandt framework too generic. Empirically: Kt(delta) ~ 5.58*N (linear,
incompressible). Transfer entropy n->delta is 0.013 bits (n adds no info).

Resolving Kt(delta(n)|n) would connect to Problem 1 but no viable attack path exists.

---

## 3. Zeta Zero Compressibility — CLOSED (S25, S36)

**CLOSED.** Zeros are GUE-random in every sense tested. 82% of Fourier modes needed,
Kt ~ 5.58*N linear, no algebraic structure. See `experiments/information_theory/kt_complexity/SYNTHESIS.md`.

---

## 4. Berry-Keating / Hilbert-Polya Hamiltonian [LITERATURE MONITORING ONLY]

No concrete self-adjoint operator with zeta zero eigenvalues is known. Even if found,
QPE requires 10^51 zeros for p(10^100). Connes' 2026 paper advances the Weil quadratic
form approach but doesn't resolve it. **No experimental work possible — monitor only.**

---

## 5. Novel Number-Theoretic Identity — CLOSED (S29)

**CLOSED.** 7 experiments: f(x) = pi(x) - R(x) is algebraically independent of all
tested bases. No computable shortcut identity exists.

---

## Priority Assessment (updated Session 37)

| Direction | Status | Action |
|-----------|--------|--------|
| Circuit complexity (#TC^0 ⊆ NC?) | OPEN, Low feasibility | THE only viable question |
| Berry-Keating | OPEN, Very Low feasibility | Literature monitoring only |
| Kt complexity | Theoretical question OPEN, attack path CLOSED | No experiments possible |
| Zero compressibility | **CLOSED** (S36) | — |
| Novel identity | **CLOSED** (S29) | — |
| Determinantal complexity | **CLOSED** (S17) | — |
| Space-time tradeoff | **CLOSED** (S23) | — |
| Non-sieve pi(x) | **CLOSED** (S16) | — |
| GF(2) algebraic | **CLOSED** (S35) | — |
| Proving impossibility | **BLOCKED** (Natural Proofs) | — |

**For detailed per-session findings, see `status/SESSION_INSIGHTS.md`.**
**For the full list of 525+ closed paths, see `status/CLOSED_PATHS.md`.**
