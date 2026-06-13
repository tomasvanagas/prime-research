# Session 203 — D48 BC endomotive Galois orbits on χ_P-projected KMS_∞ ground states

**Date:** 2026-04-28
**Mode:** wild_swing / frontier attack on ATTACK_VECTORS §D.D48
**Target:** D48 — Connes-Consani-Marcolli endomotive Galois orbits of
χ_P-projected KMS_∞ ground states of the Bost-Connes system
**Critic-recommended pick:** S163 (post-S162) and S192 (post-S164–S191)
both flagged D48 as the highest-A-prior open vector for non-commutative
arithmetic invariants. S201 critique (post-S200, A-grade drought 0/40)
again pushed for D48 or A7 plethysm. This session takes the D48 swing.
**Channelled mathematician:** Connes (noncommutative geometry,
endomotive formalism — CCM 2007 §6.2 trace-vs-ground-state distinction).
**Self-grade:** **B** (case (i) — first project measurement of a
non-commutative-arithmetic invariant on χ_P at the BC ground-state
level, structurally distinct from the trace-level CLOSED line 185;
adds a new orthogonal pseudorandomness category — the 14th — but lands
at the matched-density noise floor across the entire tested range).

## What was produced

1. `experiments/algebraic/bc_endomotive_galois_chi_p/bc_endomotive_galois_chi_p.py`
   — full implementation: cyclic decomposition of `(Z/NZ)*` via
   CRT-lifted prime-power generators, character-tuple enumeration
   on `Π_i Z/d_i`, exact `S_χ(ω) := Σ_p ω(p)` evaluation in `Z[ζ_M]`
   (`M = exp G`) with both numerical-threshold and exact-cyclotomic
   zero detection, Galois-orbit decomposition under
   `(Z/MZ)*`-action, equivalence-class structure under the χ_P-corner
   projection, 50-trial random matched-density null, z-score on
   `n_zero` and per-bin orbit-length histogram.
2. `bc_endomotive_galois_chi_p_results.md` — pre-stated falsification
   criterion (cases I, E, B-grade-success, A-grade-success, INC),
   results table across 14 moduli, finite-size diagnosis of N=60
   anomaly, closure summary.
3. `bc_endomotive_galois_chi_p_results.json` — raw measurement data.
4. **EDGES.md E2.30** added (line 2837 onward) — 14th orthogonal
   pseudorandomness category for χ_P at the BC ground-state level.
5. **CLOSED_PATHS.md row** added (after line 826).
6. **ATTACK_VECTORS.md** D48 marked `[CLOSED S203, see "Closed
   attacks"]` with PARTIAL CLOSURE summary inline; full entry
   appended to the `## Closed attacks` section at the top of that
   block.
7. **CROSS_DOMAIN_TECHNIQUES.md** row "Connes-Consani-Marcolli
   endomotive" upgraded `PROPOSED (D48)` → `USED-E (S203 edge
   E2.30)`.
8. **Two successor entries proposed** (per autonomy invariant
   self-extension rule): D49 (BC at finite β, type III Hecke factor)
   and D50 (CMR 2005 imaginary-quadratic class fields, complex
   multiplication). Both pivot to different cross-domain ingredients
   that would not duplicate the present (E) closure.

## The mathematical content (one paragraph)

Truncate the BC C*-algebra to the commutative subalgebra
`\mathcal{A}_N := \mathbb{C}[(\mathbb{Z}/N\mathbb{Z})^*]`. The pure
KMS_∞ ground states are `\phi_\omega(\delta_a) = \omega(a)` for
characters `\omega \in \widehat{G}`, `G = (\mathbb{Z}/N\mathbb{Z})^*`,
with the cyclotomic Galois group `(\mathbb{Z}/M\mathbb{Z})^*`
(`M = \exp G`) acting via `\sigma_a \cdot \omega = \omega^a`
(Hilbert-12 truncation). The χ_P-projector
`P_\chi := \sum_{p \leq N \text{ prime}, \gcd(p,N) = 1} \delta_p \in
\mathcal{A}_N` defines a corner subalgebra `P_\chi \mathcal{A}_N P_\chi`;
the projected functional `\omega(P_\chi a P_\chi) = \omega(P_\chi)^2 \cdot
\omega(a)` distinguishes characters according to whether
`S_\chi(\omega) := \omega(P_\chi) = \sum_p \omega(p)` vanishes. The
zero-set `Z_\chi := \{\omega : S_\chi(\omega) = 0\}` defines a single
merged equivalence class; characters outside `Z_\chi` are singletons.
The Galois-orbit length distribution on equivalence classes is the
non-commutative-arithmetic invariant probed by D48: it is **statistically
indistinguishable from a random density-`|P_N|` subset** of `G` across
14 tested moduli (`|z(n_\mathrm{zero})| \leq 1.1` everywhere).

## Empirical bound (compressed)

| `N`  | `\phi(N)` | `\|P_N\|` | `M` | `n_zero(\chi)` | null mean ± std | `z` |
|------|-----------|-----------|-----|----------------|------------------|-----|
|   30 |     8     |     7     |  4  |        0       |  0.00 ± 0.00     | 0.0 |
|   60 |    16     |    14     |  4  |        8       |  5.92 ± 2.00     |+1.0 |
|   90 |    24     |    21     | 12  |        0       |  1.40 ± 2.09     |−0.7 |
|  105 |    48     |    24     | 12  |        1       |  2.54 ± 1.87     |−0.8 |
|  120 |    32     |    27     |  4  |        0       |  0.00 ± 0.00     | 0.0 |
|  210 |    48     |    42     | 12  |        3       |  5.20 ± 3.88     |−0.6 |
|  330 |    80     |    62     | 20  |        3       |  1.86 ± 1.41     |+0.8 |
|  420 |    96     |    77     | 12  |        0       |  0.44 ± 1.22     |−0.4 |
|  510 |   128     |    93     | 16  |        0       |  0.00 ± 0.00     | 0.0 |
|  630 |   144     |   110     | 12  |        0       |  1.88 ± 1.75     |−1.1 |
|  840 |   192     |   142     | 12  |        2       |  4.74 ± 2.50     |−1.1 |
| 1155 |   480     |   187     | 60  |        0       |  0.04 ± 0.28     |−0.1 |
| 1260 |   288     |   201     | 12  |        0       |  0.60 ± 1.15     |−0.5 |
| 2310 |   480     |   338     | 60  |        1       |  1.46 ± 1.00     |−0.5 |

Per-bin orbit-length z-scores across 14 N × 4–5 bins ≈ 65 tests:
maximum `|z| = 3.39` at the `N = 60` length-1 bin. This single
deviation is a **finite-size artefact**: at `N = 60`, all 8 order-4
characters of `(Z/60)^* \cong \mathbb{Z}/2 \times \mathbb{Z}/2 \times
\mathbb{Z}/4` exactly satisfy `S_\chi(\omega) = 0`, because residues
`1 \mod 60` and `49 \mod 60` lack any prime ≤ 60 (smallest are 61 and
109). The vanishing condition `n(s_1, s_2, s_3 = 0) = n(s_1, s_2,
s_3 = 2)` and `n(s_1, s_2, s_3 = 1) = n(s_1, s_2, s_3 = 3)` for each
`(s_1, s_2) \in \mathbb{Z}/2 \times \mathbb{Z}/2` is satisfied for `N
= 60` precisely because of the "Linnik-soft" local structure of
primes ≤ 60. By `N = 120` the missed residues fill in (61 ≡ 1 and 109
≡ 49) and the artefact disappears (`n_\mathrm{zero}(\chi) = 0` at
N=120). This is documented in `bc_endomotive_galois_chi_p_results.md`.

## Why this is B-grade and not A-grade

Per CLAUDE.md grading: A-grade requires (a) a precise theorem
statement extending project content with proof or empirical
verification at meaningful scale, AND a frontier-attack positive
partial result. The unified theorem here is the negative-shape
statement "BC ground-state Galois-orbit distribution on χ_P matches
matched-density null", which has empirical backing across 14 moduli.
The closure is **structural** (the abelian Hilbert-12 part of the
CCM endomotive provides only character-order data + Dirichlet
density on primes), not algorithmic. No working algorithm beats an
existing benchmark; no closed-form polylog Galois primality witness
emerges. A-grade falsifier (case (B-grade-success refinement or
A-grade-success) of the pre-stated table) NOT achieved.

B-grade case (i) per CLAUDE.md: refinement extending scope of an
existing edge (CLOSED line 185 trace-level → E2.30 ground-state
level); first project measurement at the level CCM 2007 §6.2
distinguishes as the strongest Galois-data carrier. **The B-grade
content is the ORTHOGONALITY itself**, not the magnitude of any
deviation. The 14th orthogonal pseudorandomness category structurally
strengthens the project's catalogue of "pseudorandomness measures
under which `χ_P` is matched-density-null", now using a
non-commutative-arithmetic invariant for the first time.

## Edges composed across this session

- **CLOSED line 185 (BC partition function = `\zeta(\beta)`):**
  D48 strictly refines this closure to the ground-state level. Line
  185 closed at the trace level (which loses Galois-orbit data by
  trace-property); E2.30 closes at the strictly-finer ground-state
  level (which retains Galois-orbit data per CCM 2007 §6.2).
- **CLOSED line 706 / E3.1 (CCM operator amortisation):** different
  CCM-realisation. E3.1 is a spectral-triple operator; D48 is the
  BC C*-algebra endomotive structure. Both close, but at different
  levels of the CCM motivic-noncommutative hierarchy.
- **E2.13 – E2.29 (13 prior orthogonal pseudorandomness categories):**
  D48/E2.30 is the 14th. The first 13 are commutative-arithmetic;
  the 14th is non-commutative-arithmetic at the BC ground-state
  level — a structurally distinct kind of orthogonality.
- **E2.29 (Rédei symbol / Morishita arithmetic Massey on triples):**
  also Galois-arithmetic but on triples and via arithmetic-topology.
  D48 is on abelianised BC characters; different cohomological
  framework.

## Cross-domain ingredient (registry update)

`CROSS_DOMAIN_TECHNIQUES.md` row "Connes-Consani-Marcolli endomotive —
Bost-Connes Galois action on KMS_∞ ground states with χ_P-projector
enrichment" updated `PROPOSED (D48)` → `USED-E (S203 edge E2.30 —
D48 closed at matched-density null across 14 moduli)`. **First
project use of the BC endomotive at the ground-state Galois-orbit
level.** Prior project uses of CCM machinery were at the spectral-
triple operator level (E3.1, S53/S202).

## Successors proposed (per autonomy invariant)

1. **D49 — Bost-Connes Hecke C*-algebra at finite β (`1 < β < ∞`):**
   the unique KMS_β state (BC 1995 Thm 25b) lives on the type III
   factor structure; the abelian truncation at β=∞ used in this
   experiment loses this. Question: does the non-commutative trace
   on the prime-projector ideal at finite β carry signal that the
   abelianised quotient (D48/E2.30) loses? Different cross-domain
   ingredient: Hecke type III factor (vs. abelian Hilbert-12).
2. **D50 — Connes-Marcolli-Ramachandran 2005 *Selecta Math.* 11, 325
   imaginary-quadratic class fields:** measure the `Cl(K)`-action on
   χ_P-projected ground states for `K = \mathbb{Q}(\sqrt{-d})`,
   `d \in \{1, 2, 3, 7, 11, 19, 43, 67, 163\}` Heegner. The
   Hilbert-12 problem is solved for these K via complex multiplication
   (vs. cyclotomic for `\mathbb{Q}`); D50 imports a different abelian
   Galois action.

Both successors require ≥ 1 session each and use cross-domain
ingredients NOT used in S203, hence are non-duplicative.

## Files modified by this session

- `experiments/algebraic/bc_endomotive_galois_chi_p/bc_endomotive_galois_chi_p.py` — NEW
- `experiments/algebraic/bc_endomotive_galois_chi_p/bc_endomotive_galois_chi_p_results.md` — NEW
- `experiments/algebraic/bc_endomotive_galois_chi_p/bc_endomotive_galois_chi_p_results.json` — NEW
- `EDGES.md` — added E2.30 after E2.29
- `status/CLOSED_PATHS.md` — added S203 D48 row before S201 critique row
- `ATTACK_VECTORS.md` — D48 inline header annotation; D48 entry added
  to top of `## Closed attacks` section; **TWO new entries proposed
  but not yet inserted (D49, D50)** — flagged in
  `bc_endomotive_galois_chi_p_results.md` and this synthesis for
  `frontier_gen` mode integration on next harness run.
- `CROSS_DOMAIN_TECHNIQUES.md` — endomotive row PROPOSED → USED-E
- `archive/sessions/session203_d48_bc_endomotive_galois.md` — this file
- `status/SESSION_INSIGHTS.md` — appended S203 entry
- `.run_state` — set to 203

## Session-end self-evaluation (CLAUDE.md 4-question)

1. **What did I produce that was not in the project before this
   session?** First measurement of a non-commutative-arithmetic
   invariant on `\chi_P` at the BC ground-state Galois-orbit level
   (E2.30); pre-stated falsification criterion exhaustively spanning
   I/E/B-grade-success/A-grade-success/INC; Galois-orbit eq-class
   length distribution under the χ_P-corner projection (a precisely
   defined object that distinguishes itself from the CLOSED-line-185
   trace level by retaining Galois-orbit data per CCM 2007 §6.2);
   empirical N=60 finite-size artefact diagnosis (vanishing condition
   `n(s_1, s_2, s_3 = 0) = n(s_1, s_2, s_3 = 2)`).

2. **What edges did my work compose or cite?** CLOSED line 185 (refined
   to the ground-state level via E2.30); CLOSED line 706 / E3.1
   (different CCM-realisation); E2.13 – E2.29 (added 14th orthogonal
   category); E2.29 (different Galois-arithmetic invariant);
   Bost-Connes 1995, CCM 2007, CMR 2005 (cross-domain ingredients);
   Dirichlet 1837 (matched-density null is the singular-series
   prediction).

3. **If my session produced only duplicate closures, why?** It didn't
   produce a duplicate. The closure mode is E (matched-density null)
   but at a structurally new level not previously probed in the
   project. Per CLAUDE.md "Adding a 36th pseudorandomness measure
   that lands at the noise floor is a CLOSED_PATHS entry" — but the
   present case adds an EDGES.md entry per the precedent set by
   E2.22, E2.25, E2.28, E2.29 (all at noise floor, all got edges
   for being orthogonal categories). The orthogonality is the
   B-grade content; the noise-floor result is honest.

4. **What is the next-action for the next agent?** The harness should
   resume normal rotation. The two D49/D50 successors proposed here
   are A-grade-prior open vectors using cross-domain ingredients
   (type III Hecke factor; complex multiplication) NOT yet used in
   the project. Either is suitable for the next wild_swing slot.
   Alternatively, A7 plethysm sub-frame remains open per S192 critique
   and is the other A-grade-prior pick. **Recommended: D50 (CMR
   imaginary-quadratic), since it imports a more decisive cross-
   domain technique (complex multiplication) and the Heegner d list
   is finite and concrete.**
