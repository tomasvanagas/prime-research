# TODO

> **For agents starting a session:** read this file first, then `EDGES.md`,
> then `status/OPEN_PROBLEMS.md`. EDGES catalogues every real mathematical
> edge across 67+ sessions, tagged with IDs (E1.x .. E7.x). Cite edge IDs
> by name in CLOSED_PATHS entries and session syntheses (CLAUDE.md step 10).
>
> **Active critical path is FOCUS-3 (Brandt MKtP).** It is the only
> construction-flavoured attack still untouched on the only open problem.
> Pick this up in the next focused-mode session.
>
> Recurring lightweight tasks: FOCUS-5 (literature watch).
>
> Recently closed and removed from active work: FOCUS-1 (S61, S64, S66 —
> AKS sub-attacks → E7.10), FOCUS-2 (S67 — fourth-encoding sweep, E2.11
> pre-test), FOCUS-4 (S49 — large-N zeta correlation + BK probe → E1.10,
> E3.13). Closure details are in `status/CLOSED_PATHS.md` and the named
> session syntheses; the EDGES.md footer logs all recent edge additions.

---

# === CRITICAL PATH ===

## [FOCUS-3] Brandt MKtP framework deep dive

**Why this is THE critical item.** With FOCUS-1 (AKS) and FOCUS-2
(4th-encoding sweep) and FOCUS-4 (large-N zeta correlations) all closed,
Brandt is the only remaining *construction-flavoured* attack on the only
open problem (circuit complexity of π(x), per `status/OPEN_PROBLEMS.md`).

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

# === RECURRING ===

## [FOCUS-5] Literature watch (lightweight, every 2-3 weeks)

Last run: S58 (window 2026-04-05 → 2026-04-26, NO-DELTA).
See `archive/sessions/session58_literature_watch.md`.

Update `literature/state_of_art_2026.md` only if the asymptotic landscape
changes.

Sources to scan:
- arXiv math.NT recent submissions for "pi(x)", "p_n", "zeta zero"
- ECCC TR2026 for TC⁰/NC¹ separations and **Brandt MKtP follow-ups**
- Connes / Yakaboylu / van Ittersum / Ono author streams
- primecount / kim-walisch GitHub for new releases
- **Unconditional level-of-distribution past x^{5/8}.** Pascadi 2025
  (E3.12) reached x^{5/8} = 0.625 unconditionally, within 0.04 of the
  algorithmic Meissel-Lehmer threshold x^{2/3} ≈ 0.667. Any unconditional
  improvement past x^{2/3+ε} would convert several conditional π(x)
  algorithms (Lagarias-Odlyzko under GRH, etc.) into unconditional ones.
  Watch for: Pascadi follow-ups, Maynard streams, Bombieri-Vinogradov
  improvements, Friedlander-Iwaniec extensions.

Note: most months produce zero algorithmic deltas. Absence of news is
itself information.

---

# === HOUSEKEEPING (non-blocking, human review required) ===

## Flag duplicate scripts

Rule 11 says no duplicate scripts. Suspects flagged S39, awaiting human
review (do NOT delete without human approval — each pair has a companion
`_results.md`):

- [ ] `experiments/analytic/weil_optimized.py` vs `weil_optimized_v2.py` vs `weil_optimized_v3.py` (v3 same size as v1, possible revert)
- [ ] `experiments/sieve/ht_transfer_attempt.py` vs `ht_signed_transfer_v2.py`
- [ ] `experiments/circuit_complexity/k_party_nof.py` vs `k_party_nof_v2.py`
- [ ] `experiments/circuit_complexity/approx_degree.py` + `approx_degree_quick.py` + `approx_degree_small.py` + `approx_degree_prime.py` + `approx_degree_counting.py`
- [ ] `experiments/information_theory/information_shortcut.py` vs `information_shortcut_v2.py`
- [ ] `experiments/other/breakthrough_attempt_v2.py` (no v1 exists)
- [ ] `experiments/other/self_correcting_v2.py` (no v1 exists)
