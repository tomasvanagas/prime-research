# Open Problems: Viable Research Directions

Last updated: 2026-04-26 (Session 67 + S51 FOCUS-3 closure + S49 FOCUS-4 closure)

> **Note for agents:** this file maintains the canonical record of which
> direction is *technically* still open in the polylog-π(x) frame. As of
> S67 the answer is "circuit complexity, but with all known technique
> families closed."
>
> Active research now lives in `NOVELTY_CHALLENGES.md` (composition,
> frame-shift, Lean, synthesis targets) and `RESEARCH_AGENDA.md`
> (multi-session arcs). Sessions should pick targets from those files,
> not from this one — this file's content has been mined exhaustively.

These are the ONLY directions not yet proven closed. Everything else has been
tested (538+ approaches across 66 sessions) and confirmed to hit one of three
failure modes: Circularity, Equivalence, or Information Loss.

**S66 status update — Chain E is computationally cornered.** All three
FOCUS-1 sub-attacks have now closed (S47 cyclotomic-CRT line 690; S61
non-cyclotomic ring line 714; S64 Frobenius transplant line 719; S66
Bernstein 2003 strengthened-gcd this session). Every modulus-twist and
gcd-strengthening of the AKS family reduces to growing-dim r×r MPOW over
the same r ~ log²n.

**S51 status update — Chain E now closed for both known technique
families.** FOCUS-3 (Brandt MKtP) closed: Brandt 2024 (TCC, IACR ePrint
2024/687) diagonalisation is structurally welded to MKtP itself and does
not extend to fixed natural functions like π(x) mod 2 (E5.8). The hard
string is an oracle-dependent Kt-random prefix, the contradiction is
self-referential on Kt, the Chaitin-Ω density argument has no analog
for fixed totals, and the bound is on uniform time not circuits —
Brandt explicitly avoids the Williams/Hirahara algorithmic-method route
that *would* yield circuit bounds, precisely because that route faces T4
on stronger classes. With E7.10 (AKS) and E5.8 (Brandt), the only
remaining levers on E5.3 are non-AKS TC⁰ primality tests or entirely-new
lower-bound techniques. Construction-flavoured work in the current
taxonomy is exhausted.

**S67 status update — FOCUS-2 (concrete fourth-encoding sweep) closed.**
The E2.11 finite-difference pre-test ruled out 9 candidate intermediate
quantities (8 NEW + φ-summatory re-checked) as fourth encodings of π(x):
all produced GUE-white residuals indistinguishable from f(x)=π(x)−R(x),
or were entirely smooth (T_2=ΣH_n, T_3=ΣH_n²). Cumulative fourth-
encoding routes empirically closed: ~78. The "concrete candidate sweep"
sub-task is exhausted; only the abstract "fundamentally new intermediate
quantity" remains as a theoretical aspiration. See CLOSED_PATHS.md last
line and `experiments/algebraic/fourth_encoding_search/e211_pretest_focus2_results.md`.

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
- Growing-dimension matrix powering in TC^0 (genuine frontier; Chain E open primitive E5.3)
- Exploiting commutativity of polynomial ring multiplication
- A fundamentally new intermediate quantity not based on floor values or zeta zeros
- Brandt MKtP-style diagonalization for natural functions — **CLOSED S51 (E5.8)**: technique structurally welded to MKtP, does not extend to π(x) mod 2

**S66/S51 — Chain E status:** closed for both known technique families.
All three AKS-family sub-attacks closed (S61/S64/S66, E7.10) and Brandt
2024 MKtP-diagonalisation closed (S51, E5.8). The Chain-E open frontier
E5.3 is unchanged, but no construction-style attack in either the AKS
or diagonalisation-via-meta-complexity family remains. Pure new-technique
work (non-AKS TC⁰ primality tests, fundamentally new lower-bound
techniques) is the only remaining lever.

---

## 2. Time-Bounded Kolmogorov Complexity of delta(n) [THEORETICAL ONLY]

**Status:** Attack path CLOSED (S35), but the theoretical question remains open.

Kt framework reformulates but doesn't solve. Kt(T_N) = O(2^N*N) by sieve regardless
of circuit size. Brandt framework too generic. Empirically: Kt(delta) ~ 5.58*N (linear,
incompressible). Transfer entropy n->delta is 0.013 bits (n adds no info).

Resolving Kt(delta(n)|n) would connect to Problem 1 but no viable attack path exists.

---

## 3. Zeta Zero Compressibility — CLOSED (S25, S36, S45, S57, S49)

**CLOSED.** Zeros are GUE-random in every sense tested. 82% of Fourier modes needed,
Kt ~ 5.58*N linear, no algebraic structure. See `experiments/information_theory/kt_complexity/SYNTHESIS.md`.

**S49 update — FOCUS-4 large-N + Bogomolny-Keating arithmetic-correction probe also closed.**
Extended the S25/S45/S57 correlation battery to N=8000 and tested whether
the empirical pair-correlation residual D(s) = R_2_emp(s) − R_2_GUE(s)
carries the BK prime-arithmetic shape D_BK(s; T) ≈ -(2/L²) ∑_{p,k}
(log p)²/p^k · cos(2π s · k log p / L). Calibrated against a gap-shuffled
null (the only valid null per new edge E1.10): zeta gives Pearson z =
−10.85σ, phase coherence z = −2.20σ — both BELOW the null. Zeta is
*more* GUE-like than gap-shuffled controls; no detectable BK arithmetic
correction at N≤8000, T≤8148. See CLOSED_PATHS.md last line and
`experiments/analytic/zeta_structure/large_n_correlations/large_n_battery_results.md`.

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
