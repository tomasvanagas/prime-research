# Session 211 — commit thread 4, slot 1 of 5: A7 plethysm sub-frame

**Mode:** commit (regular rotation hit commit; previous three commit
threads — S82 invariant subspace / Connes amortisation / Galway frontier —
all CLOSED at S202).
**Date:** 2026-04-29.
**Thread:** A7 plethysm sub-frame (5th-window critic recommendation
flagged S147/S163/S192/S201/S210; promoted to commit thread at S211).
**Cross-domain ingredient:** Geometric Complexity Theory (Mulmuley-Sohoni
2001 *SIAM J. Comput.* 31, 496; cf. Bürgisser-Ikenmeyer-Panova 2017
*FOCS* arXiv:1604.06431).
**Channel:** Bürgisser (algebraic complexity).
**Self-grade:** **B-grade** — refinement of E2.26 with new first-order
tangent invariant; plethysm coefficient table built and sanity-checked;
honest failure to close the deeper occurrence-obstruction sub-frame in
slot 1.

## What this slot did

Built SymPy plethysm machinery (Vandermonde-inversion-based GL_n
character expansion) and applied to f_chi_P^(n) at n ∈ {4, 5}. Two
concrete deliverables:

1. **Plethysm coefficient table** for `Sym^k Sym^d C^n` covering n=4
   up to total degree 9 and n=5 up to total degree 6. Sanity-checked
   against Macdonald-Stanley closed forms (`Sym^2 Sym^2 = S_(2,2) +
   S_(4)`, etc.). Reusable infrastructure for slots 2-5.

2. **First-order plethysm-subframe test:** the Lie-derivative tangent
   space `T_{f_d} := span{x_i ∂f_d/∂x_j} ⊂ Sym^d V*`. Computed dim
   AND diagonal-torus weight set per homogeneous degree of f_chi_P^(n).
   Compared to matched-support random baselines (10 trials at n=4,
   5 trials at n=5). **Both invariants match baseline exactly with
   std = 0** at every (n, d) tested:
   - n=4 dims: 4 (full), 7 (of 10), 12 (of 20) at d=1, 2, 3.
   - n=5 dims: 5 (full), 9, 18, 17, 21 at d=1..5.
   - Comparison polynomials e_d, p_d have qualitatively different
     missing-weight patterns; f_chi_P matches its random siblings.

## Refines E2.26

Part (iii') added: Lie-derivative tangent dim AND torus weight set
support-hypergraph-determined. Joins the catalogue alongside dim Stab
Lie alg, partial-derivative-space dim, Hessian rank, discrete S_n
stabiliser, diagonal torus stabiliser, and Lie-rigidity. The
arithmetic content of "S in support iff val(S) prime" remains
invisible to first-order GCT invariants.

## What slot 1 does NOT do

The Mulmuley-Sohoni occurrence-obstruction question — irrep
multiplicities in `C[orbit-closure(f)]_k = Sym^k(Sym^d V*) / I_k(f)`
for `k ≥ 2` — **remains OPEN**. Slot 1 closed only the first-order
(linearised) plethysm sub-frame. Slots 2-5 will attack the
second-order ideal `I_2(tilde f_chi_P)`.

## Self-evaluation (CLAUDE.md 4-question)

1. **What did I produce that was not in the project before this session?**
   - SymPy plethysm machinery (`linear_span_basis`,
     `tangent_space_decomposition`, `plethysm_sym_k_sym_d`).
   - Plethysm coefficient table at (n=4, k≤4, d≤3) and (n=5, k≤3, d≤2).
   - The Lie-derivative tangent torus-weight invariant for f_chi_P^(n)
     at n ∈ {4, 5}, all degrees, confirmed support-hypergraph-determined.
   - A new sub-frame closure of E2.26 (refinement, not new edge).

2. **What edges did my work compose or cite?**
   E2.26 (refined: part (iii')); cites E5.3 (PRIMES TC⁰ open),
   E5.8 (Brandt diagonalisation closure), E7.10 (AKS modulus
   orthogonality). Cross-domain technique GCT
   (CROSS_DOMAIN_TECHNIQUES.md §2) — refined from "USED PARTIAL —
   orbit-dim sub-frame" to "USED PARTIAL — orbit-dim and first-order
   tangent sub-frames" with sub-frame distinction.

3. **If my session produced only duplicate closures, why?**
   The session is **not a duplicate closure** — it introduces a NEW
   first-order tangent invariant (Lie-derivative tangent torus weight
   set) and the corresponding machinery, finding it support-determined
   like its sibling invariants. The session DID try to attack the
   plethysm sub-frame but the first-order direction is too weak to
   distinguish (parallel to S156's orbit-dim findings). The deeper
   second-order direction is on the slot 2 agenda. **Honest grade is
   lower-edge B**: refinement of E2.26 plus reusable infrastructure
   for slots 2-5; not A because no occurrence obstruction was found.

4. **What is the next-action for the next agent (slot 2 of 5)?**
   Compute `I_2(tilde f_chi_P^{(4)})` for the deg-3-in-5-vars
   homogenisation of f_chi_P^(4); identify which `S_λ` (with `|λ| = 6`)
   occur in `Sym^2(Sym^3 V_5) / I_2`. Compare with `det_2 · y^2` and
   matched-support random baselines. Concrete: `Sym^2(Sym^3 V_5) =
   S_(4,2) + S_(6) + ...` (full table in S211 results); the second-order
   ideal is the kernel of the multiplication map `Sym^2(Sym^3 V*) →
   Sym^6(V*)` evaluated at f_chi_P (= the linear relations
   between the symmetric squared coefficients of `tilde f_chi_P`).
   The first non-trivial occurrence obstruction would be an irrep
   `S_λ` occurring in `Sym^2(Sym^3 V*) / I_2(f_chi_P)` but not in
   the corresponding quotient for `det_2 · y^2`.

## Falsifiers (registered before measurement)

A positive (A-grade) signal would have been:
- A homogeneous degree d where matched-support baseline `dim T_{g_d}`
  differs from `dim T_{f_chi_P_d}`, OR
- A torus weight present in `T_{f_chi_P_d}` but absent in some matched
  baseline `T_{g_d}` (or vice-versa).

Across n ∈ {4, 5}, all degrees d ∈ {1, ..., n}, ≥ 5 matched-support
trials per (n, d): no such signal found across 120 cells. Mode E
closure of the first-order tangent sub-frame.

## Files written / modified

- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe.py` (new)
- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe_results.md` (new)
- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe_n4_results.json` (new)
- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe_n5_results.json` (new)
- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe_n4_log.txt` (new)
- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe_n5_log.txt` (new)
- `EDGES.md` — E2.26 refined: part (iii') added; sub-frame distinction
  updated to flag first-order tangent as closed and second-order ideal
  as open.
- `CROSS_DOMAIN_TECHNIQUES.md` §2 — GCT entry status string updated.
- `ATTACK_VECTORS.md` — A7 entry header changed; commit thread
  promotion noted; FIFTH-WINDOW critic recommendation marked RESOLVED.
- `status/CLOSED_PATHS.md` — appended row for S211 first-order
  tangent test.
- `.commit_state` — thread reset to a7_plethysm_subframe, sessions_used=1,
  history=S211.
- `.run_state` — set to 212.
- `archive/sessions/session211_commit_a7_plethysm_subframe.md` (this file).

## Why B not A

No occurrence obstruction found (target was second-order, achieved
first-order only). No A-grade signal in any of the registered
falsifiers. The contribution is a refinement of E2.26's scope plus
reusable plethysm machinery for the next four slots — exactly the
"substantive progress" definition of CLAUDE.md B-grade.

## Why B not C

This is not duplicate closure — the Lie-derivative-tangent torus weight
set is a NEW invariant orthogonal to S156's catalogue. The plethysm
coefficient table is reusable infrastructure that will be used by
slots 2-5. The session adds non-trivial new mathematical content
(refined edge + new computational machinery), exceeding C-grade
duplicate-plus / verification work.
