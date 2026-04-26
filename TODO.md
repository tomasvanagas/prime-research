# TODO

> **For agents starting a session:** read this file first, then `EDGES.md`,
> then `status/OPEN_PROBLEMS.md`. EDGES catalogues every real mathematical
> edge across 67+ sessions, tagged with IDs (E1.x .. E7.x). Cite edge IDs
> by name in CLOSED_PATHS entries and session syntheses (CLAUDE.md step 10).
>
> **No active critical-path item.** All four FOCUS-N construction lines
> are now closed. Brandt MKtP (FOCUS-3) was closed in S51 — see E5.8
> in EDGES.md and the S51 row in `status/CLOSED_PATHS.md`.
>
> Recurring lightweight tasks: FOCUS-5 (literature watch).
>
> Recently closed and removed from active work: FOCUS-1 (S61, S64, S66 —
> AKS sub-attacks → E7.10), FOCUS-2 (S67 — fourth-encoding sweep, E2.11
> pre-test), FOCUS-3 (S51 — Brandt 2024 MKtP-diagonalisation does not
> extend to natural functions like π(x) mod 2 → E5.8), FOCUS-4 (S49 —
> large-N zeta correlation + BK probe → E1.10, E3.13). Closure details
> are in `status/CLOSED_PATHS.md` and the named session syntheses; the
> EDGES.md footer logs all recent edge additions.
>
> With E7.10 + E5.8 closing both known technique families on Chain E
> (AKS-style and diagonalisation-via-meta-complexity), the open problem
> on E5.3 now requires either a non-AKS TC⁰ primality test or an
> entirely-new circuit lower-bound technique. Construction-flavoured
> work is therefore exhausted in the current taxonomy.

---

# === CRITICAL PATH ===

## [FOCUS-3] Brandt MKtP framework deep dive — CLOSED in S51 (E5.8)

**Verdict (S51):** FAIL/E. Brandt 2024 (TCC, IACR ePrint 2024/687)
proves MKtP ∉ DTIME[O(n)] (Thm 1, unconditional) plus three sharper
conditional results, via length-monotonic depth-first traversal that
uses 1-Kt-randomness of Chaitin Ω prefixes to bypass the black-box
barrier (page 5).  The technique relativizes and **does** thread T4.

**However it does NOT extend to π(x) mod 2** for four orthogonal
reasons documented in `experiments/constructions/brandt_mktp/`: (O1)
the hard string is an oracle-dependent Kt-random prefix, not a fixed
function; (O2) the contradiction uses self-referential Kt on both
sides; (O3) the Chaitin-Ω density argument has no analog for fixed
naturals; (O4) the bound is on uniform time, not on circuits — and
Brandt explicitly avoids the Williams/Hirahara algorithmic-method
route precisely because that route IS subject to Natural Proofs on
stronger circuit classes.

Files: `experiments/constructions/brandt_mktp/{brandt_mktp.py,
brandt_mktp_results.md, definition.md}`. Edge: E5.8. Closure row in
`status/CLOSED_PATHS.md` (S51).

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
