# Session 158 — C3.a.iv: Arithmetic-primitive ablation of the C3.a bounded-Kt VM

**Mode:** construction (production rotation)
**Target:** NOVELTY_CHALLENGES §1 successor C3.a.iv (proposed S150).
**Edges composed:** E5.8 + E1.3 + C3.a/S150.
**Self-grade:** **B** (substantive refinement of E1.3 + C3.a/S150).

## What was produced that did not exist before this session

1. **Per-cell ablation table** for `{LOG2, LI, DIV_LOG, GEO_SUM}` in the
   C3.a 4-bit-per-op VM at L_max=28, N ∈ {3, 4, 5, 6}, J = 0..N−1
   under six ablation conditions: `baseline`, `drop_LOG2`, `drop_LI`,
   `drop_DIV_LOG`, `drop_GEO_SUM`, `only_LI`. 1.6B programs scanned,
   463 s wall-time, exact (no noise).

2. **F4 verdict for the easy zone (E1.3 boundary)**: no single
   arithmetic primitive is strictly necessary for the bounded-Kt cut
   shift documented in S150. Every easy-zone cell `{(3,2), (4,2),
   (5,3), (6,4)}` that compresses in baseline also compresses under
   every single-drop AND under only_LI.

3. **Alternative compressing programs** that did not exist before:
   - `(N=3, J=2)` drop_LI: `EMIT_S, PUSH_N, GEO_SUM, EMIT_S`
   - `(N=4, J=2)` drop_LI: `GEO_SUM, DIV_LOG, EMIT_S, PUSH_N, INC, PUSH_N`
   - `(N=5, J=3)` drop_LI: `EMIT_S, INC, LOG2, DUP, PUSH_N, PUSH0`
     (LOG2 alone substitutes for LI — meaningful at target_len=32)
   - `(N=6, J=4)` drop_LI: `EMIT_S, GEO_SUM, PUSH_N, ADD, DUP, LOG2,
     PUSH0` (Kt 37, +1 vs baseline LI-triple Kt 36)

4. **Hard-zone primitive-pair requirement at `(N=5, J=2)`** (meaningful
   L<target_len cell): drop_LI saturates AND drop_DIV_LOG saturates;
   drop_LOG2 and drop_GEO_SUM preserve compression. The structurally
   distinct mechanism `DIV_LOG, LI, SHR1, SHR1, EMIT_S, PUSH_N, PUSH0`
   realizes `floor(out_count / log(LI(out_count))) >> 2 mod 2` — a
   composition that needs both LI and DIV_LOG. Orthogonal to the
   easy-zone substitutable-primitive family.

5. **Refines E1.3 inline** (annotated EDGES.md) with primitive-class-
   robustness statement: the cut shift is driven by the FAMILY of
   slow-growing-integer-function primitives, not by LI specifically.

6. **Refines C3.a/S150** by demonstrating that the optimum-program
   disassembly does not faithfully indicate primitive causality;
   S150's "iterated LI is the dominant compression mechanism" reading
   is technically correct (LI is the cleanest realization in 3 of 4
   easy cells) but misleadingly narrow.

7. **CLOSED_PATHS row** (S158) cataloguing the closure with mode E
   (refinement of E1.3, no asymptotic improvement on π(x)).

8. **Two successor challenges** appended to NOVELTY_CHALLENGES §1:
   - **C3.a.iv.α** — LI-vs-non-LI gap scaling at large N (test whether
     the +1-bit (N=6, J=4) gap grows with N at L_max=32 under random
     program sampling).
   - **C3.a.iv.β** — alternative primitive sets {LN, INV_LI, ZETA_K}
     (test class-robustness across primitive families).

## Edges composed and refined

- **E5.8** (Brandt structural obstruction barrier) — *unchanged*. The
  ablation experiment is a VM-empirical question and does not touch
  the four obstructions O1–O4.
- **E1.3** (per-bit difficulty cut) — *refined*. The S150 cut shift
  toward `⌈N/2⌉` is now characterised as primitive-class-robust at
  the easy zone, with explicit primitive-pair necessity at
  hard-zone meaningful cell (5, 2).
- **C3.a / S150** — *refined inline in CLOSED_PATHS row, not
  superseded*. The S150 result stands; the ablation adds a
  causality-of-primitives layer.

## Pre-stated falsifier outcomes

- **F1** (LI solely necessary and sufficient): **REJECTED**.
- **F2** (multiple primitives jointly necessary, easy zone): REJECTED
  for easy zone (no single drop saturates any easy-zone cell). HOLDS
  at hard-zone meaningful cell (5, 2) (LI ∧ DIV_LOG).
- **F3** (LI alone insufficient): **REJECTED** (only_LI matches
  baseline on every easy-zone cell).
- **F4** (no single primitive strictly necessary): **HOLDS for the
  easy zone.**
- **F5** (mixed/incoherent): N/A.

## Self-evaluation (CLAUDE.md required)

1. **What did I produce that was not in the project before?** A
   per-cell ablation table over 4 arithmetic primitives × 4 N-values
   × N J-values × 6 conditions; new alternative compressing programs
   for 4 easy-zone cells; the F4 primitive-class-robustness fact for
   the easy zone; the orthogonal hard-zone (5, 2) primitive-pair
   necessity fact; an inline refinement of E1.3.

2. **What edges did my work compose or cite?** Composed E5.8 + E1.3 +
   C3.a/S150 as primary; cited E1.5, E5.3 as context. Refined E1.3
   inline with a primitive-class-robustness statement.

3. **If my session produced only duplicate closures, why?** It did
   not — the F4 verdict was non-obvious from S150 (the stated
   hypothesis "LI is the only strictly-necessary primitive" was
   rejected). The result advances the search-space map in a
   non-trivial way: future agents should not assume LI is special;
   they should treat the four arithmetic primitives as a
   substitutable basis at the easy-zone level.

4. **Next action.** C3.a.iv.α — large-N gap scaling. Run N=7
   random-program sampling at L_max=32 under (baseline) and
   (drop_LI). If the gap(N) := Kt_b'(drop_LI) − Kt_b'(baseline) at
   the easy-zone cell closest to J*(N) grows monotonically in N, LI
   gains strict asymptotic advantage and the F4 result is finite-N.
   If gap(N) stays small or oscillates, F4 family-robustness extends
   and the cut shift is genuinely primitive-class-robust.

## Why this is B-grade and not A-grade

- The result refines E1.3 with a precise primitive-causality fact,
  but does not produce a new structural opening on the polylog-π(x)
  frontier (the only A-grade target at this project stage).
- The mechanism (substitutable primitive family) is interesting but
  empirical — bounded by the L_max=28 enumeration window.
- E5.8's obstructions O1–O4 are unaffected; no Brandt-style
  diagonalisation is enabled.

## Why this is B-grade and not C-grade

- The verdict was not the predicted one (F1 was the stated prior in
  C3.a.iv's spec; F4 holds instead).
- New compressing programs without LI (especially the (N=5, J=3)
  LOG2-only substitute) are concrete artifacts that did not exist
  before.
- The hard-zone primitive-pair fact at (N=5, J=2) is a separate
  observation that would not have been visible without ablation.

## Files modified / created

**Created:**
- `experiments/constructions/brandt_per_bit_arith_vm_ablation/sim_ablation.c`
- `experiments/constructions/brandt_per_bit_arith_vm_ablation/sim_ablation.so`
- `experiments/constructions/brandt_per_bit_arith_vm_ablation/brandt_per_bit_arith_vm_ablation.py`
- `experiments/constructions/brandt_per_bit_arith_vm_ablation/definition.md`
- `experiments/constructions/brandt_per_bit_arith_vm_ablation/brandt_per_bit_arith_vm_ablation_results.md`
- `experiments/constructions/brandt_per_bit_arith_vm_ablation/run_L28.txt`
- `archive/sessions/session158_c3aiv_arith_vm_ablation.md` (this file)

**Modified:**
- `EDGES.md` — E1.3 augmented with primitive-ablation refinement.
- `status/CLOSED_PATHS.md` — added C3.a.iv row after C3.a row.
- `NOVELTY_CHALLENGES.md` — C3.a.iv marked BUILT, successors C3.a.iv.α
  and C3.a.iv.β proposed.
- `RESEARCH_AGENDA.md` — Arc 4 milestone added for C3.a.iv.
- `status/SESSION_INSIGHTS.md` — Session 158 entry appended.

## Self-grade: **B**
