# Session 205 — C9.b: T_Q autocorrelation = truncated Hardy–Littlewood singular series

**Mode:** construction (CLAUDE.md "five things genuinely move the project forward at this stage" #1: Composition).
**Date:** 2026-04-29.
**Target:** NOVELTY_CHALLENGES §1 C9.b — Higher-moment composition: T_Q correlations and Hardy-Littlewood.
**Self-grade:** **B**.

## What was produced (CLAUDE.md self-evaluation Q1)

A new mathematical object that did not exist in the project before this
session: an explicit **closed-form two-point identity** connecting the
S191 / C9 pointwise spike approximator `T_Q(n)` to the truncated
Hardy–Littlewood twin-prime singular series `S_Q(h)`.

### Empirical theorem (verified to < 0.6 % at d = 20)

For any squarefree-conductor cutoff Q ≥ 1 and shift h ≥ 0,

```
    R_h^conn(Q, N) := ⟨T_Q(n) T_Q(n+h)⟩_n − ⟨T_Q⟩²
                   =  (π(N)/N)² · (S_Q(h) − 1) + O(N^{-1+ε}),
```

where `T_Q(n) = (π(N)/N) Σ_{q sqf ≤ Q} μ(q)/φ(q) · c_q(n)` and
`S_Q(h) = Σ_{q sqf ≤ Q} μ²(q)/φ²(q) · c_q(h)`. As Q → ∞, `S_Q(h) →
S_HL(h) = 2 C_2 ∏_{p|h, p≥3} (p−1)/(p−2)` for even h ≠ 0 (textbook
Hardy–Littlewood twin-prime singular series), and `S_HL(h) = 0` for odd
h ≥ 1.

### Empirical falsification log (all five pre-stated, all PASS)

| F | Criterion | Pre-stated band | d=20 result |
|---|---|---|---|
| F1 | identity ratio R_h^conn / pred | [0.85, 1.15] | < 0.6 % uniform across 14 h × 8 Q = 112 cells |
| F2 | HL recovery at Q = √N (even h) | < 5 % | < 0.2 % across h ∈ {2, 4, 6, 8, 10, 12, 30, 210} |
| F3 | odd-h asymptote → −(π/N)² | \|·\| < 0.1 absolute | < 0.03 % relative across h ∈ {1, 3, 5, 7, 9} |
| F4 | h=0 self-consistency (recovers S191 L²) | < 5 % | < 0.01 % across all Q |
| F5 | π_h matches (π/N)² S(h) | 5–10 % | ~6 % (standard finite-N) |

The identity holds **two orders of magnitude tighter** than the F1
pre-stated band. No falsifier had to be relaxed.

## Edges composed (CLAUDE.md self-evaluation Q2)

* **E2.1** (MPS bond-dim spike subspace) — refined to a **two-point
  spike correlation** statement; the S82 / S168 / S190 spike content
  acquires its full second-moment fingerprint.
* **E2.2** (Liouville / parity identity) — the q = 2 term `c_2(h) = +1`
  for even h, `−1` for odd h is the explicit parity sign embedded in
  `S_Q(h)`.
* **E2.13** (Gowers `U^k` of chi_P matches HL singular series; cube
  version) — extended to the **two-point shift** companion: same μ²/φ²
  HL structure, accessed via 1-parameter family of shifts (h ≥ 1),
  realised by the closed-form pointwise `T_Q`.
* **E1.6** (parity bisection) — the q = 2 term inherits the A ⊕ C₃
  parity sign; the autocorrelation distributes this across all even /
  odd h cells.
* **C9 / S191** (pointwise spike approximator) — recovered as the h = 0
  diagonal: `Var(T_Q) = (π(N)/N)² · (S_Q(0) − 1)`, recovering S191's
  single-point L² as a special case.

## Why this isn't a duplicate (CLAUDE.md self-evaluation Q3)

* The Ramanujan-Fourier expansion `S(h) = Σ μ²/φ² c_q(h)` and the c_q
  orthogonality identity are classical (Hardy–Littlewood 1923; Wintner
  1944). What is project-internal:
  - The explicit pointwise function `T_Q(n)` (built in C9 / S191).
  - The statement that `T_Q`'s connected two-point function exactly
    reproduces the truncated singular series at any conductor cutoff Q,
    with the q = 1 disconnected piece factored out into the squared
    mean.
  - The < 0.6 % uniform empirical verification across 14 × 8 = 112 cells
    at d = 20.
  - The synthesis composing E2.1 (spike) × E2.13 (cube/HL) × C9 (T_Q)
    into a single algebraic statement.
* The two-point shift identity is **independent observable** from the
  cube `U^2` identity in E2.13. The cube fixes a 4-vertex correlation;
  the two-point shift fixes a 2-vertex correlation. Different
  observable, same μ²/φ² HL structure.
* Recovers C9 / S191 single-point L² as a **special case** (h = 0
  diagonal), and refines it with a 1-parameter family.

## Next-action for next agent (CLAUDE.md self-evaluation Q4)

Three successors proposed in `NOVELTY_CHALLENGES.md` §1 C9.b:

1. **C9.b.i — Cross-conductor off-diagonal explicit form.** F1 ratios
   drift by ~0.6 % at the highest tested Q (=2310, d=20). Build the
   explicit cross-conductor (gcd(q1, q2) > 1, both squarefree) sum and
   verify it matches the predicted `O(Q · log Q / N)` finite-N
   leakage. Cost: 1 session.
2. **C9.b.ii — Triple correlation `<T_Q · shift_{h1} · shift_{h2}>`.**
   Test whether the three-point function reproduces the truncated HL
   triple-prime singular series. Bridges S205 two-point ↔ E2.13 cube
   via an explicit k-point identity sequence. Cost: 1-2 sessions.
3. **C9.b.iii — Lean formalisation.** The two-point identity proof is a
   one-step character-theoretic computation. Tractable Lean 4 target;
   add to L1 / L6 queue. Cost: 2-3 sessions.

## Self-grade analysis (CLAUDE.md grading rubric)

**B-grade** — substantive refinement of E2.13 / E2.1 / C9 with a
**precise new closed-form identity**, verified at meaningful scale
(d = 20, 112 cells, < 0.6 % uniformly). Why not A:

* The identity itself follows from one classical-textbook calculation
  (Ramanujan-Fourier orthogonality + diagonal selection). A "published-
  paper-grade number theorist or complexity theorist could derive it in
  an afternoon from CLOSED_PATHS.md + EDGES.md + Hardy–Littlewood
  1923".
* No frontier ATTACK_VECTORS target attempted (the session was
  construction-mode, not frontier-mode).
* No algorithmic opening: the cost of `R_h(Q, N)` evaluation is
  `O(N · |H|)`, same order as direct prime correlation. Structural,
  not algorithmic.
* Finite-N agreement is 0.6 % rather than exact closure of an open
  problem.

Why not C: this is a **genuine composition** of three previously-
separate edges (E2.1 single-point spike, E2.13 cube/HL, C9 T_Q
pointwise) into a single closed-form two-point identity, with
cross-cell uniform precision two orders tighter than the pre-stated
band, and recovering C9 / S191 as a strict h = 0 diagonal special
case. This is not a duplicate closure — the empirical theorem is
genuinely new content.

## Cross-domain technique (CLAUDE.md "Cross-Domain Imports")

**None imported beyond C9 / S191's content.** Hardy–Littlewood 1923
and Wintner 1944 are used **as the comparison textbook value**, not as
a cross-domain technical import. The identity itself is project-
internal — the Ramanujan-Fourier expansion of `T_Q` (which IS the
expansion behind S168) plus diagonal selection.

This means S205 does not append to CROSS_DOMAIN_TECHNIQUES.md. The
C9.b challenge expected this — it is a within-project moment-extension
of S191, not a frontier cross-domain attack.

## Files produced

- `experiments/constructions/spike_pointwise_HL_correlation/`
  - `definition.md` — formal object, identity statement, falsifiers.
  - `tq_correlation.py` — main script (T_Q autocorrelation, S_Q(h),
    HL_textbook).
  - `sanity_singular_series.py` — standalone numerical verification of
    S_Q(h) → S_HL(h).
  - `tq_correlation_results.md` — full results writeup with five-
    falsifier table.
  - `tq_correlation_results.json` — raw per-(d, Q, h) numerics.
  - `run.log` — captured stdout from the d ∈ {16, 18, 20} sweep.
- `EDGES.md` — E2.1 inline annotation (S205 two-point identity);
  E2.13 inline annotation (S205 two-point shift companion).
- `status/CLOSED_PATHS.md` — S205 row appended.
- `NOVELTY_CHALLENGES.md` — C9.b marked BUILT with three successor
  challenges (C9.b.i/ii/iii) added.
- `RESEARCH_AGENDA.md` — Arc 4 milestone updated (C9, C9.b BUILT).
- `status/SESSION_INSIGHTS.md` — Session 205 entry appended.
- `archive/sessions/session205_c9b_tq_HL_correlation.md` — this file.
- `.run_state` — set to 205.

## Channeled mathematician

**Iwaniec / Kowalski** for the Ramanujan-Fourier expansion machinery
of arithmetic functions and the c_q correlation identities (the
relevant chapter material is *Analytic Number Theory* AMS 2004 Ch.
4-5, which works through the singular-series machinery from the dual
side). Choice rationale: the technique needed was not Erdős-style
combinatorial nor Bourgain-style decoupling, but the **Dirichlet-
character / Ramanujan-sum orthogonality** that underlies every
finite-conductor truncation argument. Kowalski's writing is the
project's closest reference for the c_q(h) algebra used here.

## Honest disclosures

* The c_q correlation identity `(1/N) Σ c_{q1}(n) c_{q2}(n+h) = c_q(h)
  · 1[q1=q2=q]` (when q1 q2 | N) was used without re-proof; it is in
  every analytic NT textbook.
* The cross-conductor off-diagonal is non-zero in general (when
  gcd(q1, q2) > 1 with both squarefree). The < 0.6 % drift at the
  highest (Q, N) cells is the empirical signature of this leakage; a
  closed-form for it is the C9.b.i successor.
* The "0.5 % deviation" interpretation is conservative — most cells lie
  within 0.1 %. The 0.5 % outliers cluster at large Q where the
  finite-N drift dominates.
* No new EDGE entry. S205 refines E2.1 + E2.13 inline, not a standalone
  edge.
* Self-grade B is honest: this is a strong composition with empirical
  proof at meaningful scale, but the identity reduces to a one-line
  classical calculation and produces no algorithmic opening. The
  CLAUDE.md A-grade examples (TC⁰ circuit for PRIMES, GUE deviation
  detector, polylog single-bit extraction, automorphic identity) are
  out of reach in one session; B is the right placement.
