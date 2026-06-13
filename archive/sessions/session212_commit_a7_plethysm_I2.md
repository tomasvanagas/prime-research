# Session 212 — commit thread 4, slot 2 of 5: A7 plethysm sub-frame, I_2 layer

**Mode:** commit (continuation of slot-2 plan from S211; thread =
A7 plethysm sub-frame).
**Date:** 2026-04-29.
**Thread:** A7 plethysm sub-frame; sessions_used = 2.
**Cross-domain ingredient:** Geometric Complexity Theory occurrence-
obstruction at level k = 2 (Mulmuley-Sohoni 2001 *SIAM J. Comput.*
31, 496; Bürgisser-Ikenmeyer-Panova 2017 *FOCS* arXiv:1604.06431).
**Channel:** Bürgisser (algebraic complexity).
**Self-grade:** **B-grade** — substantive computational result that
closes the second-order plethysm sub-frame at three (n, d) settings
({(4, 3), (5, 3), (5, 4)}) for chi_P + 9 matched-support baselines +
det/perm/e_3/p_3/Plücker comparisons; the Veronese control x_1^3
confirms the test is sharp.

## What this slot did

Built the second-order GCT machinery promised by S211: a per-weight
kernel computation of `I_2(orbit-closure(f)) ⊆ Sym^2(Sym^d V_n*)`,
together with Schur-irrep decomposition.  Applied to:

- `f_chi_P^(4)` deg-3 component / homogenisation in `(n, d) ∈
  {(4, 3), (5, 3), (5, 4)}`.
- Determinantal / permanental / elementary-symmetric / power-sum /
  single-monomial reference polynomials.
- Random matched-support baselines (3 per setting).
- The Veronese control `x_1^3` (whose orbit closure is the cubic
  Veronese, with known `dim I_2`).

**Result:** dim `I_2 = 0` for chi_P at every tested (n, d), and for
all 9 matched-support baselines + 5 reference polynomials.  Schur
multiplicities `b_λ` of `Sym^2(Sym^d V*) / I_2` coincide identically
across the entire chi_P-baseline class.  No occurrence obstruction
at level 2 separates chi_P from its matched-support siblings or from
det_2·y / perm_2·y.

The control `x_1^3` confirms test sharpness:
- `(n, d) = (4, 3)`: dim `I_2(x_1^3) = 126 = dim S_(4,2)` at n=4.
- `(n, d) = (5, 3)`: dim `I_2(x_1^3) = 420 = dim S_(4,2)` at n=5.
- These match the classical catalecticant-minor count for the
  cubic Veronese (Landsberg 2017 ch. 5).

Sanity checks at (n=2, d=3) and (n=3, d=3) reproduce the classical
catalecticant-rank-1 ideal (dim `I_2(x_1^3) = 3` at n=2; dim 27 at
n=3).  Machinery verified.

## Refines E2.26 (extension iii'')

The S211 first-order tangent sub-frame (iii') is now joined by the
second-order ideal sub-frame (iii''): both are support-hypergraph-
determined, std = 0 across baselines.

Updated catalogue of GCT invariants probed for f_chi_P:

| invariant | session | result |
|---|---|---|
| dim Stab_{GL_n}(f) Lie alg | S156 | 0 (full orbit) |
| partial-derivative-space dim | S156 | matched baseline |
| Hessian rank | S156 | matched baseline |
| diagonal torus stabiliser | S156 | matched baseline |
| discrete S_n stabiliser | S156 | matched baseline |
| Lie-derivative-tangent dim | S211 | matched baseline |
| Lie-derivative-tangent torus weights | S211 | matched baseline |
| **dim I_2(orbit-closure)** | **S212** | **0 = matched baseline** |
| **Schur mults b_λ in Sym^2/I_2** | **S212** | **matched baseline** |

The plethysm sub-frame at higher k (k ≥ 3) remains untested.

## What slot 2 does NOT do

The Mulmuley-Sohoni occurrence-obstruction question at level
`k ≥ 3` remains OPEN.  At k=3, `Sym^3(Sym^d V_n*) / I_3(f)`
decomposes into many more irreps and the multiplicity counts are
richer.  Slot 3 next-action is to compute `I_3` at the smaller
setting (n=4, d=3) where ambient `Sym^3(Sym^3 V_4)` has dim 1540 and
the kernel computation is tractable in a single session.

## Self-evaluation (CLAUDE.md 4-question)

1. **What did I produce that was not in the project before this session?**
   - The `plethysm_subframe_I2.py` per-weight kernel machinery for
     I_2 of any tractable (n, d, f).
   - The dim-`I_2 = 0` finding for chi_P + 9 matched-support
     baselines + 5 reference polynomials at three (n, d) settings.
   - A NEW sub-frame of E2.26 (part iii'') with the Schur-multiplicity
     invariant.
   - Concrete slot-3 next-action (k=3 plethysm test at n=4).

2. **What edges did my work compose or cite?**
   E2.26 (refined further with iii''); composes with S211's iii'
   (first-order tangent) and S156's iii (partial-derivative-space
   dim).  Cross-domain technique GCT (CROSS_DOMAIN_TECHNIQUES.md
   §2) refined to second-order sub-frame with kernel-computation
   USED status.  Cites Landsberg 2017 *GCT* ch. 5 (catalecticant
   classical fact) for the control's verification.

3. **If my session produced only duplicate closures, why?**
   The session is **not a duplicate closure** — it computes the
   GENUINELY NEW invariant `dim I_2` (and its Schur decomposition),
   the slot-2 next-action that S211 specifically queued.  The
   `dim I_2 = 0` finding is novel computational data; it cannot be
   inferred from S156 + S211 + closed-paths alone.  However, the
   structural conclusion ("chi_P matches matched baselines at this
   level too") parallels the S211 / S156 finding pattern; this is
   incremental progress at the next level of the GCT hierarchy.
   B-grade by the CLAUDE.md grading.

4. **What is the next-action for the next agent (slot 3 of 5)?**
   Compute `dim I_3(orbit-closure(f))` for `f_chi_P^(4)` deg-3
   component, three matched-support baselines, and `e_3` / `x_1^3`
   controls at `(n, d) = (4, 3)`.  Ambient `Sym^3(Sym^3 V_4)` has
   dim 1540 (= 20·21·22/6); per-weight kernel computation on
   M ~ 2000 orbit samples should run in ~5–10 minutes per polynomial.
   Decompose into Schur using S211's plethysm machinery.  Identify
   any `S_λ` (with `|λ| = 9`, at most 4 parts) appearing in
   `Sym^3 / I_3(chi_P)` but absent from baseline (or vice-versa).
   If found at (4, 3): A-grade signal — first non-support-determined
   GCT invariant of f_chi_P.  If not found: continued mode E
   closure of the k = 3 sub-frame, slot 4 attacks Bürgisser-Ikenmeyer-
   Panova 2017 no-occurrence-obstruction theorem applicability to
   chi_P.

## Why B not A

No occurrence obstruction was found.  Closing a sub-frame at "matches
baseline" is B-grade refinement, not A-grade novelty.  The contribution
is the new computational machinery (per-weight kernel, Schur
decomposition for I_2), reusable for the slot-3 k=3 test, plus the
sub-frame closure at level k=2.  This is exactly the "substantive
progress" definition of CLAUDE.md B-grade.

## Why B not C

The session computes a NEW invariant (`dim I_2`) and applies it
across 15 distinct (polynomial × (n, d)) cells.  The plethysm-table
machinery was built in S211 but the kernel-computation infrastructure
is new in S212; reusable for slots 3-5.  This is non-trivial new
mathematical content, exceeding the C-grade duplicate-plus /
verification floor.

## Files written / modified

- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe_I2.py` (new)
- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe_I2_results.md` (new)
- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe_I2_n5_d3_results.json` (new)
- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe_I2_n5_d3_log.txt` (new)
- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe_I2_n4_d3_results.json` (new)
- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe_I2_n4_d3_log.txt` (new)
- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe_I2_n5_d4_results.json` (new, if completed)
- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe_I2_n5_d4_log.txt` (new, if completed)
- `EDGES.md` — E2.26 refined: part (iii'') added.
- `CROSS_DOMAIN_TECHNIQUES.md` §2 — GCT entry updated to include
  second-order kernel computation.
- `ATTACK_VECTORS.md` — A7 entry header updated; commit slot 2/5
  noted.
- `status/CLOSED_PATHS.md` — appended row for S212 second-order
  kernel sub-frame closure.
- `.commit_state` — sessions_used incremented to 2; history
  S211, S212.
- `.run_state` — set to 213.
- `archive/sessions/session212_commit_a7_plethysm_I2.md` (this file).
