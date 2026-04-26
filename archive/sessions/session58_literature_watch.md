# Session 58 — FOCUS-7 Literature Watch (2026-04-05 → 2026-04-26)

**Date:** 2026-04-26
**Mode:** normal (FOCUS-7 from `TODO.md`, last ran S47)
**Verdict:** **NO-DELTA.** No publication in the window changes the
asymptotic landscape, opens a new path, or invalidates any closed path.

## Why this watch

Project is in mature steady state.  Per `CLAUDE.md` and `TODO.md`,
literature watch is a recurring, lightweight FOCUS-7 task — the
appropriate work when no concrete experiment is open.  Last full
watch was S47 (window 2026-04-05 → 2026-04-25).  Three weeks of
window since then; project guidelines call for a new pass.

## Method

Single agent did the search.  Coverage:

| Source                            | Span scanned                  |
|-----------------------------------|-------------------------------|
| arXiv math.NT recent submissions  | All 245 entries for April 2026 |
| ECCC TR2026                       | TR26-040 .. TR26-061           |
| Author streams                    | Connes, Yakaboylu, van Ittersum, Ono, Granville, Tao |
| GitHub                            | kimwalisch/primecount, PrimeCounting/PrimeCounting |

Targeted arXiv search terms: `pi(x)`, `nth prime`, `prime counting
function`, `zeta zeros`, `explicit formula`, `Riemann hypothesis`,
`sieve`, `prime distribution`, `Hardy-Littlewood`, `Selberg`.

## Findings (what was added)

Four minor catalog entries appended to
`literature/state_of_art_2026.md` (S58 update block).  None changes
the asymptotic barrier.

### 1. Inoue, arXiv:2604.05733 (April 2026) — INCREMENTAL

*"Small gaps between consecutive zeros of the Riemann zeta-function."*
Resonance-correlation method combining Montgomery pair-correlation
and Montgomery-Odlyzko approach.  Proves μ < 0.50895 under RH (was
0.515).  Statistic on zero gaps; no algorithmic implication.  Marginal
addition to GUE-statistics line in §2.4 if anything.

### 2. primecount post-v8.4 (master, April 15-25, 2026) — INCREMENTAL

~30 commits on GitHub master since v8.4 release (April 6).  Notable:
- "New lockfree thread load balancer (#114)" (Apr 24)
- "Improve load balancing on many core systems" (Apr 20)
- "Use more accurate zeta zeros" (Apr 16)
- "Use Brun–Titchmarsh theorem" (Apr 15)
- "Improve dist_approx formula" (Apr 15)

v8.5 release imminent.  All constant-factor work; touches the same
nth_prime bracketing path our V10 implementation uses.

### 3. Kakkar, arXiv:2604.02383 (April 2026) — NONE

*"Neural Prime Sieves: Density-Driven Generalisation and Empirical
Evidence for Hardy–Littlewood Asymptotics."*  PrimeFamilyNet ML model.
~60-90% search-space reduction with >95% recall.  Confirms §5.4
barrier (no exactness).  No complexity improvement.

### 4. Li, arXiv:2604.14596 (April 2026) — NONE (flagged)

*"Prime-Zero Duality: Fractal Geometry, Renormalization-Group Flow,
and an Information-Ontological Framework for Number Theory."*
103-page speculative single-author preprint with explicit
"speculative spirit" disclaimer.  Author acknowledges core claim is
unproven.  Same profile as Kilictas-Alpay TG kernel (debunked
S12/S30) — no peer review, no institutional affiliation, no
computational pathway.  Flagged for completeness so future watches
recognise the pattern.

## Other in-window items checked, no project relevance

- arXiv:2604.15396 (González-Negrín, Salem integral RH equivalence note)
- arXiv:2604.03051 (CUE moment computation)
- arXiv:2604.16530 (zeta deficiency identity)
- ECCC TR26-061 (Raz, non-commutative ABP lower bound)
- ECCC TR26-056 (Frick et al., Z₂-topological sign-rank framework — adjacent
  to χ_P sign-rank but no number-theoretic application)
- ECCC TR26-051 (Liu-Mazor-Pass, time-bounded Kolmogorov / Pessiland —
  conceptually adjacent to Brandt MKtP line, no direct π(x) connection)
- ECCC TR26-052 (Hirahara-Nanashima, Pessiland characterisation)
- ECCC TR26-047 (Ko, dynamic Boolean cell-probe lower bound)
- ECCC TR26-043 (Kush, multilinear ABP lower-bound barrier)

No TC^0/NC^1 separation.  No matrix-powering circuit lower bound.
No number-theoretic hardness in the window.  No Brandt MKtP follow-ups.

## Author-stream specifics

- **Connes:** No new April 2026 submission post-arXiv:2602.04022 (Feb 2026).
- **Yakaboylu (arXiv:2408.15135):** Still v15 (March 11, 2026); no v16.
  Note: S47 entry references "v10," now stale-by-naming; content unchanged.
  Updated implicitly by absence-of-v16.
- **van Ittersum:** No new April 2026 NT submission.  Latest related is
  arXiv:2412.19180 (quasi-modularity, Dec 2024) and arXiv:2507.00147 (July 2025).
- **Ono:** No April 2026 prime-detection partition follow-ups.
- **Granville / Tao:** No April 2026 prime/zeta papers found.  Tao's most
  recent NT-adjacent posting is the March 2026 Erdős-Graham paper on
  "products of consecutive integers with unusual anatomy" — irrelevant.

## Implication for the project

The mature-state hypothesis from S47 holds.  Most weeks produce no
relevant news.  Brandt MKtP (S30) remains the only identified
theoretical technique that could in principle bypass Natural Proofs
en route to circuit lower bounds for π(x), and no follow-up appeared
this window.

The only genuinely open research direction
(`status/OPEN_PROBLEMS.md`: circuit complexity of π(x), specifically
the Ω(log x) vs O(x^{1/2+ε}) gap and growing-dim MPOW in TC^0)
remains unchanged.  Aggarwal (October 2025) is still the most recent
direct statement on n-th-prime complexity.

## Files written

- `literature/state_of_art_2026.md` (S58 block appended; header date bumped)
- `status/SESSION_INSIGHTS.md` (Session 58 entry appended)
- `archive/sessions/session58_literature_watch.md` (this file)

## Files NOT written

- No new `experiments/` (literature watch is a non-experimental task).
- No `novel/` updates (no novel finding produced).
- No `status/CLOSED_PATHS.md` updates (no path closed).
- No `status/OPEN_PROBLEMS.md` updates (none resolved).

## Next steps

Per the FOCUS-7 cadence note, the next watch should fire in ~3 weeks
(estimated session 70).  Earlier triggers:
- A van Ittersum or Connes announcement on partition-detection or
  spectral-triple identities.
- Any post-v8.4 primecount asymptotic claim (currently constant-factor only).
- A Brandt MKtP follow-up extending to P or to natural functions —
  the strongest external watch signal for the project.

## State of project

No breakthrough.  No closures.  No new attack surfaces.  The polylog
gap remains: exact O(x^{2/3}) (`algorithms/v10_c_accelerated.py`)
vs O(polylog) ~50% digits (`R^{-1}(n)`).  Steady state preserved.
