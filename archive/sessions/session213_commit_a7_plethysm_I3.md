# Session 213 — commit thread 4, slot 3 of 5: A7 plethysm sub-frame, I_3 layer

**Mode:** commit (continuation of slot-3 plan from S212; thread =
A7 plethysm sub-frame).
**Date:** 2026-04-29.
**Thread:** A7 plethysm sub-frame; sessions_used = 3.
**Cross-domain ingredient:** Geometric Complexity Theory occurrence-
obstruction at level k = 3 (Mulmuley-Sohoni 2001 *SIAM J. Comput.*
31, 496; Bürgisser-Ikenmeyer-Panova 2017 *FOCS* arXiv:1604.06431;
Iarrobino-Kanev 1999 *Power Sums, Gorenstein Algebras and
Determinantal Loci* LNM 1721 — cubic-Veronese Hilbert-function
benchmark for sanity verification).
**Channel:** Bürgisser (algebraic complexity).
**Self-grade:** **B-grade** — substantive computational result that
closes the third-order plethysm sub-frame at (n, d) = (4, 3) for
chi_P + 3 matched-support baselines + e_3/p_3/Plücker references,
with the Veronese control x_1^3 confirming sharpness (dim I_3 =
1320 = entire complement of S_(9) in Sym^3(Sym^3 V_4), the
catalecticant rank-1 ideal).

## What this slot did

Built the third-order GCT machinery promised by S212: a per-weight
kernel computation of `I_3(orbit-closure(f)) ⊆ Sym^3(Sym^d V_n*)`,
together with full Schur-irrep decomposition.  Applied to:

- `f_chi_P^(4)` deg-3 component = `x_1 x_2 x_3 + x_1 x_2 x_4 + x_1 x_3 x_4`.
- `e_3`, `p_3`, single-monomial `x_1 x_2 x_3` reference polynomials.
- 3 random matched-support baselines.
- The Veronese control `x_1^3` (whose orbit closure is the cubic
  Veronese ν_3(P^3); known I_3 from Iarrobino-Kanev 1999).

**Result:** dim `I_3 = 0` for chi_P + 6 other chi_P-class
polynomials at (n, d) = (4, 3); the Schur multiplicities `b_λ` of
`Sym^3(Sym^3 V_4) / I_3` coincide identically across all 7 cells
with the full plethysm decomposition `Sym^3(Sym^3 V_4) = S_(9) +
S_(7,2) + S_(6,3) + S_(5,2,2) + S_(4,4,1)` (each multiplicity 1).
No occurrence obstruction at level 3 separates chi_P from its
matched-support siblings or from the natural references.

The control `x_1^3` confirms test sharpness: dim I_3(x_1^3) = 1320
= S_(7,2) + S_(6,3) + S_(5,2,2) + S_(4,4,1) (the catalecticant rank-1
ideal at degree 3 in 4 vars; quotient = S_(9) of dim 220 = the cubic-
Veronese Hilbert function value `dim Sym^9 V_4`).

Sanity verifications at (n, d) = (2, 3) and (3, 3): dim I_3(x_1^3) =
10 and 165 respectively, both reproducing classical cubic-Veronese
Hilbert-function values from Iarrobino-Kanev 1999 textbook.

Verification at M = 3000 with separate seed: largest weight block
(γ = (3, 2, 2, 2), block size 24) has all 24 singular values
exceeding SVD tolerance by 10^10 (top SV = 3.16 × 10^5, smallest =
8.92 × 10^3, tol = 2.10 × 10^-7).  The kernel is empty by an
enormous margin, not a numerical artefact.

## Refines E2.26 (extension iii''')

The S211/S212 sub-frames (iii'/iii'') are now joined by the third-
order ideal sub-frame (iii'''): all three are support-hypergraph-
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
| dim I_2(orbit-closure) | S212 | 0 = matched baseline |
| Schur mults b_λ in Sym^2/I_2 | S212 | matched baseline |
| **dim I_3(orbit-closure)** | **S213** | **0 = matched baseline** |
| **Schur mults b_λ in Sym^3/I_3** | **S213** | **matched baseline (full plethysm)** |

The plethysm sub-frame at higher k (k ≥ 4) remains untested.

## What slot 3 does NOT do

The Mulmuley-Sohoni occurrence-obstruction question at level k ≥ 4
remains OPEN.  At k = 4, `Sym^4(Sym^3 V_4) / I_4(f)` has ambient
dim 8855 and decomposes into ~10 Schur components; per-weight kernel
computation requires ~10× the cost of k = 3 (~30 minutes per
polynomial with M ~ 12000 samples).

The slot-4 next-action will instead PRIORITISE testing the Bürgisser-
Ikenmeyer-Panova 2017 no-occurrence-obstruction theorem applicability
to chi_P.  S211/S212/S213 cumulatively show chi_P matches matched-
support baselines at levels k = 1, 2, 3 with full plethysm
multiplicities; this parallels the BIP no-OCB result for det vs
padded-perm.  If chi_P can be shown to inherit BIP no-OCB structurally
(via a containment chi_P ∈ closure(GL · g) for natural g, or via a
representation-theoretic theorem), this is a B-grade structural result.

## Self-evaluation (CLAUDE.md 4-question)

1. **What did I produce that was not in the project before this session?**
   - The `plethysm_subframe_I3.py` per-weight kernel machinery for
     I_3 of any tractable (n, d, f) with full Schur-irrep
     decomposition.
   - The dim-`I_3 = 0` finding for chi_P + 6 matched baselines / refs
     at (n, d) = (4, 3); verified at M = 3000 with 10^10 SV gap.
   - The new computational fact `dim I_3(x_1^3) = 1320 = S_(7,2) +
     S_(6,3) + S_(5,2,2) + S_(4,4,1)` at n = 4, d = 3 (catalecticant
     ideal at k = 3, complement of S_(9) in Sym^3 plethysm).
   - A NEW sub-frame of E2.26 (part iii''') with the third-order
     Schur-multiplicity invariant.
   - Sanity verification reproducing Iarrobino-Kanev 1999 textbook
     values at (n, d) = (2, 3) and (3, 3).
   - Concrete slot-4 next-action (BIP 2017 no-OCB applicability for
     chi_P; fallback I_4 at (4, 3)).

2. **What edges did my work compose or cite?**
   E2.26 (refined further with iii'''); composes with S211's iii'
   (first-order tangent), S212's iii'' (second-order ideal), and
   S156's iii (partial-derivative-space dim).  Cross-domain technique
   GCT (CROSS_DOMAIN_TECHNIQUES.md §2) refined to third-order
   sub-frame with kernel-computation USED status.  Cites Iarrobino-
   Kanev 1999 *Power Sums, Gorenstein Algebras and Determinantal Loci*
   LNM 1721 ch. 1 (cubic-Veronese Hilbert function in degree k for the
   sanity verification).  Cites S211 plethysm decomposition table for
   Sym^3(Sym^3 V_4) = S_(9) + S_(7,2) + S_(6,3) + S_(5,2,2) +
   S_(4,4,1).

3. **If my session produced only duplicate closures, why?**
   The session is **not a duplicate closure** — it computes a
   GENUINELY NEW GCT invariant (`dim I_3` and its Schur decomposition)
   at one structural level deeper than S212; the slot-3 next-action
   that S212 specifically queued.  The `dim I_3 = 0` finding is
   novel computational data; it cannot be inferred from S211 / S212
   / closed-paths alone.  However, the structural conclusion
   ("chi_P matches matched baselines at this level too") parallels
   the S211 / S212 finding pattern; this is incremental progress at
   the next level of the GCT hierarchy.  B-grade by the CLAUDE.md
   grading.

4. **What is the next-action for the next agent (slot 4 of 5)?**

   PRIORITY: Test the Bürgisser-Ikenmeyer-Panova 2017 no-occurrence-
   obstruction theorem applicability to chi_P at (n, d) = (4, 3).
   BIP proved that for det_m vs perm_n^{m²-n} (padded permanent), no
   S_λ occurs in Sym^d / I_d(perm_n) but absent from Sym^d / I_d(det_m).
   Our S211 + S212 + S213 cumulatively show chi_P exhibits the same
   pattern (no occurrence obstruction at k = 1, 2, 3 vs matched
   baselines).  If chi_P is shown to inherit BIP no-OCB barrier
   structurally — e.g., by exhibiting a containment chi_P ∈ closure(GL ·
   g) for some BIP-style natural g, or by a representation-theoretic
   structural argument — this is a B-grade result: chi_P would be the
   FIRST natural-NT polynomial known to inherit the BIP no-OCB
   barrier.  An A-grade signal would be a proof that the no-OCB
   barrier extends to chi_P at all k ≥ 1, closing the entire plethysm
   sub-frame in one move.

   FALLBACK: compute `dim I_4` at (n, d) = (4, 3); ambient
   `Sym^4(Sym^3 V_4)` has dim 8855 (= C(20+3, 4)).  Per-weight kernel
   with M ~ 12000 orbit samples: ~30 minutes per polynomial (~10× cost
   of k=3).  Largest weight blocks ~200, so SVD step becomes the main
   cost.  Decompose via S211 plethysm table (10 Schur components for
   k = 4 plethysm at n = 4).

   Slot 5 (WRAP) will synthesise the 5-session arc.  If slot 4 yields
   the BIP-applicability result, the synthesis becomes a first-of-its-
   kind structural claim.  If only mode E closures stack, slot 5
   declares the entire plethysm sub-frame closed at this (n, d) with a
   unified theorem plus a CLOSED_PATHS row tying chi_P to the BIP
   no-OCB pattern.

## Why B not A

No occurrence obstruction was found at k = 3 either.  Closing a
sub-frame at "matches baseline" is B-grade refinement, not A-grade
novelty.  The contribution is the new computational machinery (per-
weight kernel + Schur decomposition for I_3, reusable for slot 4 if
the fallback I_4 path is taken), plus the sub-frame closure at level
k = 3, plus the Iarrobino-Kanev verification at (n, d) = (4, 3).
This is exactly the "substantive progress" definition of CLAUDE.md
B-grade.

## Why B not C

The session computes a NEW invariant (`dim I_3` plus its Schur
decomposition) and applies it across 8 distinct polynomial cells at
(n, d) = (4, 3) plus 2 sanity-verification cells at (n, d) = (2, 3)
and (3, 3).  The kernel-computation infrastructure is genuinely new
in S213 (the S212 I_2 code only handles pairs); reusable for slot 4
if k = 4 is attempted.  This is non-trivial new mathematical content,
exceeding the C-grade duplicate-plus / verification floor.

## Falsifiers (registered before measurement)

A positive (A-grade) signal would have been:
- A partition `λ ⊢ 9` (≤ 4 parts) with `b_λ(chi_P) ≠ b_λ(g)` for
  any matched-baseline `g`, OR
- `dim I_3(chi_P) ≠ dim I_3(e_3)`, OR
- `dim I_3(chi_P) ≠ dim I_3(plucker)`, OR
- A non-zero kernel `dim (I_3)_γ` for chi_P at any of the 220 weight
  blocks.

Across `5` partitions × `7` polynomials × `220` weight blocks at
M = 2000 + M = 3000 verification: no such signal found.  Mode E
closure of the third-order plethysm sub-frame.

## Files written / modified

- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe_I3.py` (new)
- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe_I3_results.md` (new)
- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe_I3_n4_d3_results.json` (new)
- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe_I3_n4_d3_log.txt` (new)
- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe_I3_n2_d3_results.json` (new, sanity)
- `experiments/algebraic/gct_chi_p_orbit/plethysm_subframe_I3_n3_d3_results.json` (new, sanity)
- `EDGES.md` — E2.26 refined: part (iii''') added; sub-frame distinction
  updated.
- `CROSS_DOMAIN_TECHNIQUES.md` §2 — GCT entry status updated to include
  third-order ideal sub-frame.
- `ATTACK_VECTORS.md` — A7 entry header updated; commit slot 3/5 noted.
- `status/CLOSED_PATHS.md` — appended row for S213 third-order kernel
  sub-frame closure.
- `status/SESSION_INSIGHTS.md` — updated with S213 entry.
- `.commit_state` — sessions_used incremented to 3; history S211, S212,
  S213; slot3_summary recorded.
- `.run_state` — set to 214.
- `archive/sessions/session213_commit_a7_plethysm_I3.md` (this file).
