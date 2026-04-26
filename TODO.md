# TODO

> **For agents starting a session:** read this file first, then `EDGES.md`,
> then `status/OPEN_PROBLEMS.md`. EDGES catalogues every real mathematical
> edge across 67+ sessions, tagged with IDs (E1.x .. E7.x). Cite edge IDs
> by name in CLOSED_PATHS entries and session syntheses (CLAUDE.md step 10).
>
> # === Active priority order ===
>
> **FOCUS-3 (Brandt MKtP) is the new critical path.** With FOCUS-1 (AKS
> sub-attacks) closed S61/S64/S66 — all three hitting the same E5.3 wall
> per E7.10 — and FOCUS-2 (4th-encoding sweep) closed S67 — 9/9 candidates
> mode I via the E2.11 pre-test — Brandt is the only construction-flavoured
> critical-path option still untouched. It has been logged for 30+ sessions
> without engagement. **Pick this up in the next focused-mode session.**
>
> Other active items (lower priority): FOCUS-4 (3-point correlation at
> N≥10⁴), FOCUS-5 (literature watch). FOCUS-1 and FOCUS-2 retained below
> as archived references with closure pointers.
>
> # === Recently closed ===
>
> - **FOCUS-1 entire (S61, S64, S66):** all three AKS sub-attacks closed
>   FAIL with mode E or E+I or E+C; modulus / coefficient-ring / gcd-
>   strengthening are orthogonal to depth (E7.10).
> - **FOCUS-2 entire (S67):** all 9 candidate intermediate-quantity
>   summatories closed mode I via E2.11 pre-test (Δᵏ ratio → 2.0 white-
>   noise signature). Six WHITE-A (GUE-equivalent to π−R), two WHITE-B
>   (entirely smooth, no prime info), one mixed.
> - Connes operator scaling (S53), π(x) mod q for q ∈ {2..13} via
>   Liouville (S55+S56), Borel+Fejér hybrid (S52), 3-point GUE correlation
>   (S57), N=2000 zeta structure (S45), Helfgott-Thompson M(x) benchmark
>   (S60, sieve route tight via E6.7+E7.6).
>
> # === Recent EDGES extensions ===
>
> - **S60 batch (10 edges + 1 caveat):** E1.7, E1.8, E2.11, E2.12, E3.12,
>   E6.7, E6.8, E6.9, E7.8, E7.9, plus E6.3 finite-size DCT caveat.
> - **S66 batch (2 edges + 2 refinements):** E1.9 (φ 2D rank, 22nd
>   pseudorandomness measure), E7.10 (AKS modulus-twist orthogonality
>   theorem), K_min non-monotonicity caveat to E3.11, AKS-gcd-extraction
>   note under E5.3.
> - **Critique-46 batch (1 edge):** E7.11 (convergence-acceleration family
>   systematically exhausted).

---

# === CRITICAL PATH ===

## [FOCUS-3] Brandt MKtP framework deep dive  [NEW critical path, S67+]

**Why this is now THE critical item.** With FOCUS-1 (AKS) and FOCUS-2
(4th-encoding sweep) both fully closed, Brandt is the only remaining
*construction-flavoured* attack on the only open problem (circuit
complexity of π(x), per `status/OPEN_PROBLEMS.md`).

Brandt 2024 (TCC) proved `MKtP ∉ DTIME[O(n²)]` via a diagonalisation
that **bypasses Natural Proofs**. S30 flagged this as "the only known
technique that could lead to unconditional superpolynomial lower bounds
without hitting the barriers." S39 confirmed no follow-up papers
through April 2026. The project has logged Brandt for 30+ sessions
without anyone engaging with the technique.

### Concrete actions

1. **Read Brandt 2024 carefully.** TCC proceedings or arXiv listing.
   Identify the diagonalisation skeleton: what is being diagonalised
   against, what natural-function class is the result for, what
   exact ingredient bypasses Razborov-Rudich.
2. **Adapt to π(x) mod 2.** This function is total, computable in
   sub-exponential time O(x^{2/3}), conjectured outside TC⁰, and
   pseudorandom in 22+ measures (`novel/pseudorandomness_of_pi.md`).
   Does Brandt's hypothesis class admit π mod 2 as a natural target,
   or does the construction require an artificially-defined function
   (like MKtP itself) that doesn't extend to natural NT?
3. **Even a *conditional* Brandt-style lower bound on π(x) mod 2
   would be the first non-trivial circuit lower bound the project has
   produced.** If Brandt's hypothesis only applies to artificial
   classes, document this rigorously as the closure mode and elevate
   to `proven/circuit_size_barrier.md`.
4. Save analysis at `experiments/circuit_complexity/brandt_mktp/`
   (theory dump + any concrete computations on small N).
5. Cite EDGES IDs in the closure: this is a Chain E attempt against
   E5.3 via a non-AKS technique; T4 (Natural Proofs) is the constraint
   to thread.

This is exploratory — most likely closure mode is "Brandt's hypothesis
doesn't extend to natural functions like π mod 2", but the upside is
unique among the remaining options. **It is also the only critical-path
work that satisfies CLAUDE.md's "Construction is encouraged" rebalance**:
even a careful theoretical construction (formal definitions + small-N
simulation) qualifies. Save under `experiments/constructions/brandt_mktp/`
if the work produces a concrete object/circuit, otherwise the
`brandt_mktp/` path above.

---

# === ARCHIVED / CLOSED ===

## [FOCUS-1] AKS growing-dim MPOW in TC^0  [CLOSED S61/S64/S66 — archive]

> **All three sub-attacks closed FAIL.** Per E7.10 (AKS modulus-twist
> orthogonality theorem), every modulus / coefficient-ring / gcd-
> strengthening twist of the AKS test reduces to growing-dim r×r MPOW
> at the same `r ~ log²n`, leaving E5.3 untouched.

S47 closed cyclotomic-CRT splitting (AKS-prescribed r is prime in 21/22
sampled n). Three uncovered sub-attacks remain — each tractable for a
focused-mode session:

### Sub-attack 1: Smaller-r AKS variant (Bernstein 2003) — CLOSED (S66)

**Status:** CLOSED, FAIL (E with weak C). See CLOSED_PATHS.md line ~722,
`experiments/circuit_complexity/aks_alternative/bernstein_smaller_r/`,
and `archive/sessions/session66_bernstein_smaller_r.md`.

**Headline:** Empirical r already satisfies r/log²n mean = 1.207 on
the S47 22-sample grid (mean r − log²n = 25.05, r prime in 21/22).
"Smaller r" is empirical reality already; the question is *deterministic
correctness* at this r. The Bernstein-style gcd certificate works
empirically (gcd extracts factor of n in 13/13 composite failures
including all 7 Carmichaels) but the integer-gcd subroutine on
O(log n)-bit integers is in NC¹ (Hesse-Allender-Barrington 2002), NOT
known in TC⁰ — same NC¹/TC⁰ frontier as growing-dim MPOW (E5.3) it
purportedly replaces. Closure mode (E): the strengthening trades one
frontier problem for another at the same frontier.

### Sub-attack 2: Non-cyclotomic ring decomposition

AKS works in `Z_n[x]/(x^r - 1)`. Replace the modulus by an irreducible
`f(x)` of degree exactly polylog(n) over `Z_n`. The cyclotomic-CRT
split (S47) closes because `x^r - 1` factors trivially when r is prime;
choosing `f` irreducible bypasses that.

**Concrete next step:** prototype with `f(x) = x^d + a` (Eisenstein-style)
for d = ceil(log^2 n), check (a) irreducibility certification cost in
TC^0, (b) whether the AKS congruence still characterises primality
when the ring is `Z_n[x]/f(x)`.

Save under `experiments/circuit_complexity/aks_alternative/non_cyclotomic_ring/`.

### Sub-attack 3: Healy-Viola Frobenius transplant — CLOSED (S64)

**Status:** CLOSED, FAIL (E + I). See CLOSED_PATHS.md line 719,
`experiments/circuit_complexity/aks_alternative/frobenius_transplant/`,
and `archive/sessions/session64_frobenius_transplant.md`.

**Headline:** primes do NOT satisfy mod-q AKS for non-trivial a (0/399
across 19 primes; trivial a≡0 mod q always passes but carries zero
info). The Frobenius identity `(a+x)^n ≡ a^n + x^n` is mod-n specific;
reducing mod q≠n loses it (residual polynomial fills 8-102 of r
coefficients, by Lucas's theorem on binom(n,k) mod q). Even granted
the test worked, the q-twist gives only a `log(q)/log(2)` constant-factor
saving — both schemes reduce to length-O(log n) growing-dim r×r MPOW
over F_q (E5.3), unchanged.

### Backup objective: ALL THREE NOW CLOSED — Chain E computationally cornered (S66)

**S66 milestone:** with sub-attacks 1 (S66), 2 (S61), 3 (S64) all
closed, the AKS-family sub-attack space is exhausted. Every modulus-
twist, ring-replacement, and gcd-strengthening of the AKS test
reduces to growing-dim r×r MPOW over the same r ~ log²n. The
"computationally cornered" marker is now active in
`status/OPEN_PROBLEMS.md`.

**Remaining levers on Chain E:**
1. **Brandt MKtP** (FOCUS-3 below) — un-engaged framework with a
   possible diagonalization argument for natural functions.
2. **Fundamentally new lower-bound technique** — open problem class,
   no concrete proposal.
3. **A non-AKS TC⁰ primality test** using only scalar operations
   (long-standing aspiration since S15; no active attack line).

Future Chain-E construction work should look outside the AKS
family — within it, no further leverage point is known.

---

# === SECONDARY (sharpening) ===

## [FOCUS-2] 4th encoding of pi(x) — concrete candidate sweep — CLOSED (S67)

**Status:** CLOSED, FAIL (mode I). See CLOSED_PATHS.md last line,
`experiments/algebraic/fourth_encoding_search/e211_pretest_focus2_results.md`,
and `archive/sessions/session67_focus2_e211_pretest.md`.

**Headline:** All 9 listed candidates (T_1, T_2, T_3, Ψ × 3 regimes, σ_2,
σ_3, Q, T_6) close at the E2.11 pre-test stage *without reaching PSLQ*.
T_2 = ΣH_n and T_3 = ΣH_n² are entirely smooth (closed-form, residual at
f64 precision — they carry no prime information at all). The other 7
produce GUE-white residuals indistinguishable from f(x) = π(x) − R(x) —
they encode the *same* zeta-zero oscillation, hence a polylog evaluator
for any of them would yield polylog π(x). None has a polylog evaluator.
Cumulative fourth-encoding routes empirically closed: ~78. The "concrete
candidate sweep" sub-task of FOCUS-6 is now exhausted; only abstract
"fundamentally new intermediate quantity not based on floor values or
zeta zeros" remains as theoretical aspiration.

> **Body of FOCUS-2 archived (S67):** the original 6-candidate sweep
> definition + E2.11 pre-test methodology lived here. With all 9
> candidates closed mode I, the active body has been collapsed to this
> archive note. See `experiments/algebraic/fourth_encoding_search/
> e211_pretest_focus2_results.md` for the 9-candidate verdict table
> and `archive/sessions/session67_focus2_e211_pretest.md` for synthesis.
> Methodology preserved in `EDGES.md` E2.11.

---

# === SECONDARY (sharpening) ===

## [FOCUS-4] Higher-order zeta-zero correlation at N >= 10^4

S57 closed 3-point correlation at N=2000 to GUE noise floor, but
explicitly noted: *"arithmetic corrections [Conrey-Snaith] become
visible only at large lag in the Fourier-conjugate variable, not here."*
The S57 windows-disjoint test at L=32 gave c_3 ~ 10^-3, but the
disjoint-window count was only ~62 windows.

### Concrete action

1. Use Odlyzko's tabulated zeros (online) to extend N from 2000 to
   10^5 (factor 50 in sample).
2. Re-run S57's pair / triple / cumulant battery with proper disjoint
   windows.
3. Specifically test the **Conrey-Snaith arithmetic correction** for
   pair correlation:
   `R_2_arith(s) = R_2_GUE(s) + (1/2) (log(gamma_n / 2 pi))^{-2} A_HL(s)`
   where `A_HL` is a Hardy-Littlewood-style correction.
4. If the arithmetic correction is empirically confirmed at this scale,
   it is genuinely novel structural info (not in current literature
   for pi(x)) and could feed back into Chain D-style attacks if anyone
   ever revisits them.

Save at `experiments/analytic/zeta_structure/large_n_correlations/`.

## [FOCUS-5] Literature watch (recurring, lightweight)

Last run: S58 (window 2026-04-05 -> 2026-04-26, NO-DELTA).
See `archive/sessions/session58_literature_watch.md`.

Re-run every 2-3 weeks. Update `literature/state_of_art_2026.md` only if
the asymptotic landscape changes.

Sources to scan:
- arXiv math.NT recent submissions for "pi(x)", "p_n", "zeta zero"
- ECCC TR2026 for TC^0/NC^1 separations and **Brandt MKtP follow-ups**
- Connes / Yakaboylu / van Ittersum / Ono author streams
- primecount / kim-walisch GitHub for new releases
- **Unconditional level-of-distribution past x^{5/8} (NEW from S60 EDGES
  extension).** Pascadi 2025 (E3.12) reached x^{5/8} = 0.625 unconditionally,
  within 0.04 of the algorithmic Meissel-Lehmer threshold x^{2/3} ≈ 0.667.
  Any unconditional improvement past x^{2/3+epsilon} would convert several
  conditional pi(x) algorithms (Lagarias-Odlyzko under GRH, etc.) into
  unconditional ones — this is a publishable status-change for our project.
  Watch for: Pascadi follow-ups, Maynard streams, Bombieri-Vinogradov
  improvements, Friedlander-Iwaniec extensions.

Note: most months produce zero algorithmic deltas. Absence of news is
itself information.

---

# === HOUSEKEEPING & BENCHMARKING (non-blocking) ===

## [HOUSEKEEPING] Flag duplicate scripts (human review required)

Rule 11 says no duplicate scripts. Suspects flagged S39, awaiting human review:

- [ ] `experiments/analytic/weil_optimized.py` vs `weil_optimized_v2.py` vs `weil_optimized_v3.py` (v3 same size as v1, possible revert)
- [ ] `experiments/sieve/ht_transfer_attempt.py` vs `ht_signed_transfer_v2.py`
- [ ] `experiments/circuit_complexity/k_party_nof.py` vs `k_party_nof_v2.py`
- [ ] `experiments/circuit_complexity/approx_degree.py` + `approx_degree_quick.py` + `approx_degree_small.py` + `approx_degree_prime.py` + `approx_degree_counting.py`
- [ ] `experiments/information_theory/information_shortcut.py` vs `information_shortcut_v2.py`
- [ ] `experiments/other/breakthrough_attempt_v2.py` (no v1 exists)
- [ ] `experiments/other/self_correcting_v2.py` (no v1 exists)

Each pair has a companion `_results.md`. Do NOT delete without human approval.

> **Removed S60:** `[BENCHMARK] Helfgott-Thompson O(x^{3/5}) for M(x)`
> dropped from this file. E6.7 (HKM time-space curve O~(N^{8/15}) /
> O~(N^{1/3})) combined with E7.6 (Lucy-DAG pebbling lower bound
> Omega(x^{5/6}/ln x)) shows the sieve route is asymptotically tight to
> within x^{0.034}. H-T sits inside this Pareto frontier and benchmarking
> it would not change project status. The script `experiments/sieve/
> mertens_speedup.py` remains in place for archival reference.
