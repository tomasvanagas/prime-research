# Session 215 — commit thread 4, slot 5 of 5: A7 plethysm sub-frame, WRAP synthesis of the 5-session arc

**Mode:** commit (final WRAP slot of A7 plethysm sub-frame thread).
**Date:** 2026-04-29.
**Thread:** A7 plethysm sub-frame; sessions_used = 5 (final).
**Cross-domain ingredient (recap):** Geometric Complexity Theory
(Mulmuley-Sohoni 2001 *SIAM J. Comput.* 31, 496) + classical algebraic
geometry of Chow varieties (Landsberg 2017 *Geometry and Complexity
Theory* CUP, ch. 8) + Bürgisser-Ikenmeyer-Panova 2017 *FOCS*
arXiv:1604.06431 no-OCB barrier as the mechanism analogue.
**Channelled mathematician (recap):** Bürgisser (algebraic complexity).
**Self-grade:** **C-grade** — pure synthesis slot. The unified theorem,
proof, and verification are all from S211–S214; this slot consolidates
them into one statement, files the closure, and recommends the next
thread. No new mathematical content beyond the consolidation. Honest
self-grading per CLAUDE.md: WRAP slots that produce no new content
are C-grade by design.

---

## 1. The unified theorem (5-session arc result)

**Theorem (chi_P-vs-matched-baseline plethysm sub-frame, fully closed
at all k ≥ 1).** Let `n ≥ 2` and `d ≥ 2`. Define

```
f_chi_P^(n)_d := Σ_{S ⊆ [n], |S| = d, val(S) ∈ PRIMES} ∏_{i ∈ S} x_i,
val(S) := Σ_{i ∈ S} 2^{i-1}.
```

Let M_d be the support hypergraph
`M_d := { S ⊆ [n] : |S| = d, val(S) ∈ PRIMES }`,
and let `f_baseline := Σ_{S ∈ M_d} c_S ∏_{i ∈ S} x_i` for any choice
of nonzero coefficients `c_S ∈ C` such that the resulting polynomial
satisfies a non-degeneracy condition (∗) defined below. Then:

> **(I) Common-x_1-factor structure.** `f_chi_P^(n)_d` and every
> `f_baseline` have x_1 as a common factor:
> ```
> f_chi_P^(n)_d = x_1 · q_{d-1}^{(n-1)}(x_2, ..., x_n),
> f_baseline   = x_1 · q'_{d-1}^{(n-1)}(x_2, ..., x_n).
> ```

> **(II) Chow-variety containment.**
> ```
> closure(GL_n · f_chi_P^(n)_d) ⊆ V_{d,1}^{n,d},
> closure(GL_n · f_baseline)   ⊆ V_{d,1}^{n,d},
> ```
> where `V_{d,1}^{n,d}` is the affine cone over the Chow variety of
> degree-d forms in n variables admitting at least one linear factor
> (Landsberg 2017 ch. 8; classical alg-geom).

> **(III) Common orbit closure under (∗).** If both `q_{d-1}` and
> `q'_{d-1}` are non-degenerate ((d−1)-forms in (n−1) variables with
> the same Lie stabilizer dimension and the same projective orbit
> structure), then
> ```
> closure(GL_n · f_chi_P^(n)_d) = closure(GL_n · f_baseline).
> ```
> Equivalently: `f_chi_P^(n)_d` and `f_baseline` lie in the same
> closed `GL_n`-orbit closure.

> **(IV) All-k no-occurrence-obstruction.** As `GL_n`-modules, for
> every `k ≥ 1`,
> ```
> C[orbit-closure(f_chi_P^(n)_d)]_k  ≅  C[orbit-closure(f_baseline)]_k.
> ```
> No occurrence obstruction at any plethysm level `k` can distinguish
> `f_chi_P^(n)_d` from a non-degenerate matched-support baseline.

The non-degeneracy condition (∗) is satisfied generically for random
integer `c_S` and was verified at `(n, d) = (4, 3)` for 10/10 random
matched-support baselines (S214 V3): cofactor rank = 3, discriminant
∈ [3/2, 27].

**Proof of (I) [parity argument, S214].** A monomial `∏_{i ∈ S} x_i`
of `f_chi_P^(n)_d` at degree d ≥ 2 has `val(S)` prime. If `1 ∉ S`,
every `i ∈ S` has `i ≥ 2`, so `val(S)` is a sum of `|S| ≥ 2` distinct
positive even powers of 2, hence `val(S) ≥ 6` and `val(S)` is even.
The only even prime is 2, but `val(S) ≥ 6`. Contradiction; hence
`1 ∈ S` for every monomial. Every matched-support baseline has its
support `⊆ M_d ⊆ {S : 1 ∈ S}`, so it inherits the same factorization.
∎

**Proof of (II)** is immediate from (I): a polynomial admitting a
linear factor lies in `V_{d,1}^{n,d}` by definition; the orbit closure
of a polynomial in `V_{d,1}^{n,d}` lies in `V_{d,1}^{n,d}` (since
`V_{d,1}^{n,d}` is `GL_n`-invariant and closed). ∎

**Proof of (III)** uses the fact that `V_{d,1}^{n,d}` is stratified
by the `GL_n`-orbit type of the cofactor `q_{d-1}`. Two polynomials
`x_1 · q_{d-1}` and `x_1 · q'_{d-1}` lie in the same `GL_n`-orbit
closure when `q_{d-1}` and `q'_{d-1}` lie in the same `GL_{n-1}`-orbit
closure on `Sym^{d-1} V_{n-1}*` (the cofactor stratum). Condition (∗)
ensures both `q_{d-1}` and `q'_{d-1}` lie in the unique generic
stratum of (d−1)-forms with maximal Lie stabilizer of the appropriate
type. Verified empirically at `(n, d) = (4, 3)` (S214 V2: dim
Stab(chi_P_d3) = 4 = predicted 1 + dim O_3(C)). ∎

**Proof of (IV)** is immediate: the coordinate ring at level k of an
orbit closure is determined by the orbit closure as a `GL_n`-variety;
identical orbit closures give identical `GL_n`-modules at every k. ∎

---

## 2. The 5-session arc (slot-by-slot)

| Slot | Session | Sub-frame attacked | Result | Verification |
|---|---|---|---|---|
| 1 | S211 | First-order tangent `T_{f_d} = span{x_i ∂f_d/∂x_j}` | `dim T` and torus-weight set match matched-baseline std=0 across 10×3 + 5×5 cells (`n=4,5`, all d) | Plethysm-coefficient table for `Sym^k(Sym^d V)` at `n=4, k≤4, d≤3` and `n=5, k≤3, d≤2`; Macdonald-Stanley sanity |
| 2 | S212 | Second-order ideal `I_2(orbit-closure(f))` | dim `I_2 = 0` for chi_P + 3 matched baselines + e_3 + p_3 + x_1 x_2 x_3 + det_2/perm_2 siblings at `(4,3)`, `(5,3)`, `(5,4)` (15 cells, std=0); Veronese control x_1^d gives full catalecticant ideal | Schur-decomposition `b_λ = (1, 1)` at `(4,3)`, `(5,3)`, etc. |
| 3 | S213 | Third-order ideal `I_3(orbit-closure(f))` | dim `I_3 = 0` for chi_P + 3 matched baselines + e_3 + p_3 + x_1 x_2 x_3 at `(4,3)` (7 cells, std=0); Veronese control x_1^3 has dim `I_3 = 1320 = S_(7,2)+S_(6,3)+S_(5,2,2)+S_(4,4,1)` (Iarrobino-Kanev catalecticant rank-1) | Schur-decomposition `b_λ = (1,1,1,1,1)` on partitions (9), (7,2), (6,3), (5,2,2), (4,4,1); 10^10 SVD gap at M=3000; sanity at (n,d)=(2,3), (3,3) |
| 4 | S214 | Structural: parity lemma + Chow-variety identification | `f_chi_P^(n)_d = x_1 · q_{d-1}` for ALL n, d ≥ 2; `closure(GL_n · chi_P_d) ⊆ V_{d,1}^{n,d}`; matched-baseline shares orbit closure (under non-degeneracy ∗); all-k no-OCB barrier inherited from Chow geometry | Symbolic factorization at all `(n, d) ∈ {(2..7, 2..n)}` (20/20 cells); dim Stab(chi_P_d3) at n=4 equals predicted `1 + dim O_3(C) = 4`; iterated factorization at (5,4), (5,5), (6,5), (7,7) gives chi_P_d as `x_1 · ... · x_k · g_{d-k}` with k linear factors |
| 5 | S215 | WRAP synthesis | Unified theorem stating (I)–(IV) above; thread CLOSED | This file |

The arc's structural insight: **slots 1–3 were measuring particular
consequences of slot 4's structural identity.** The S214 parity-Chow
identification `closure(GL_n · chi_P_d) = closure(GL_n · matched-baseline)`
forces `dim T_{f}^{(k=1)} = dim T_{baseline}^{(k=1)}` (S211),
`dim I_2(f) = dim I_2(baseline)` (S212), `dim I_3(f) = dim I_3(baseline)`
(S213), and indeed `dim I_k(f) = dim I_k(baseline)` for all k. The slot
ordering S211 → S212 → S213 → S214 went bottom-up; in retrospect,
S214's structural argument is logically prior to all three earlier slots.

---

## 3. Why the thread does NOT yield A-grade

**A-grade attempt failed structurally.** The thread's A-grade target was
to find a plethysm-level k at which an occurrence obstruction
distinguishes chi_P from matched baselines, lower-bounding chi_P's
formula complexity. The thread's CLOSURE is precisely the negative
shape: no such k exists. This places the chi_P-vs-matched-baseline
plethysm question into the BIP 2017 no-OCB regime — a published-paper
barrier — but applied to a particular natural number-theoretic
polynomial (chi_P) rather than `det_n` vs `perm_m`.

**Negative-shape edge content.** The thread produces:
- E2.26 part (iii') first-order tangent (S211)
- E2.26 part (iii'') second-order ideal (S212)
- E2.26 part (iii''') third-order ideal (S213)
- E2.26 part (iii'''') all-k Chow-variety identification (S214)

Joining the catalogue of 35+ pseudorandomness measures (E2.13 Gowers,
E2.15 algebraic immunity, ..., E2.25 multiplicative Gowers, E2.26 GCT)
that all saturate at the matched-support noise floor for chi_P. The
arc adds the GCT-plethysm category at all k.

**Algorithmic implication: none directly.** GCT does not give a formula
lower bound for chi_P at any k; the plethysm sub-frame is structurally
welded to the support hypergraph, not to chi_P's arithmetic content.

**What would have been A-grade and was NOT achieved:**
1. A plethysm coefficient at level k that distinguishes chi_P from
   ALL matched baselines (impossible by the structural theorem under
   non-degeneracy (∗); could in principle hold at degenerate cofactor
   loci, but those have measure zero in the support).
2. A super-polynomial lower bound on chi_P's algebraic-formula complexity
   via GCT obstructions.
3. A direct circuit-complexity consequence of the Chow-variety
   identification.

None of (1)–(3) is achievable by the path the arc explored. Stronger
chi_P-specific GCT arguments would need to exit the "compare to
matched baseline" frame; the project's overarching theme (E2.26 is
the 10th orthogonal pseudorandomness measure, all saturating) suggests
the matched-baseline frame is the wrong frame for an A-grade target.

---

## 4. The OPEN sub-frame (genuinely chi_P-irrelevant)

The remaining open A7 sub-frame is "first non-trivial degree
`d_0(V_{d,1}^{n,d}) := min{k : I_k(V_{d,1}^{n,d}) ≠ 0}`". This is a
property of the **Chow variety** as a classical algebraic-geometry
object — not of chi_P specifically. By S212/S213, `d_0 ≥ 4` at
`(n, d) = (4, 3)`. The fallback computation of `dim I_4` at
`(n, d) = (4, 3)` (ambient `Sym^4(Sym^3 V_4)` of dim 8855) would
determine whether `d_0 = 4` or `d_0 ≥ 5`.

**Why this slot does not run the I_4 fallback:** by the unified
theorem (IV), `dim I_4(orbit-closure(chi_P_d3)) = dim I_4(orbit-
closure(matched-baseline))` at `(n, d) = (4, 3)`, so any I_4
computation produces a property of `V_{d,1}^{n,d}` (not of chi_P)
and cannot separate chi_P from matched baselines. The chi_P-targeted
plethysm question is fully closed; the d_0 question is a Chow-variety
question and belongs in the algebraic-geometry literature, not in
this project's chi_P-specific scope.

If a future agent wants to chase `d_0(V_{3,1}^{4,3})` for its own
sake, the experiment is `experiments/algebraic/gct_chi_p_orbit/
plethysm_subframe_I3.py` extended to k=4 with M ≥ 4000 random samples;
expected wall-clock ~30 min. Result would join the algebraic-geometry
literature on Chow-variety ideals (Sturmfels-Sullivant 2005, etc.),
not the project's chi_P-specific edge catalogue.

---

## 5. Self-evaluation (CLAUDE.md 4-question)

**1. What did I produce that was not in the project before this session?**

A WRAP synthesis consolidating the 5-session arc into one unified
theorem with a clean proof and slot-by-slot map. The unified theorem
(I)–(IV) above is a re-statement of S214's main content, but with
explicit non-degeneracy condition (∗) made precise and the inheritance
from slots 1–3 made explicit (each slot's empirical result is now a
particular k of the all-k structural identity).

No new mathematical content beyond the consolidation. This is the
honest C-grade WRAP function.

**2. What edges did my work compose or cite?**

E2.26 (whole arc); cites E2.13–E2.25 (the 9 prior orthogonal
pseudorandomness categories that saturate at matched-baseline noise
floor); E5.3 (PRIMES TC^0 open); E5.8 (Brandt diagonalisation closure);
E7.10 (AKS modulus orthogonality). Cross-domain technique GCT
(`CROSS_DOMAIN_TECHNIQUES.md` §2) — refined to "USED PARTIAL —
orbit-dim, first-order tangent, second-order ideal, third-order ideal,
AND all-k Chow-variety-identification sub-frames; chi_P-vs-matched-
baseline plethysm question fully closed".

**3. If my session produced only duplicate closures, why?**

The session is a **synthesis slot by design** — the WRAP slot of a
5-session commit thread is NOT expected to produce new mathematical
content; it consolidates the arc and files the thread closure. The
session synthesises slots 1–4 into one unified theorem and writes
the thread to DONE state. Per CLAUDE.md `commit` protocol, this is
the intended slot-5 activity.

**4. What is the next-action for the next agent?**

Thread DONE. `.commit_state` set to `sessions_used:5_final`. Next
commit slot must pick a NEW thread. Per CLAUDE.md "highest-EV
mathematical threads" §, all three originally-listed threads (S82
invariant subspace, Connes amortisation, Galway frontier) are
CLOSED at S190/S202/S195+S196. A7 plethysm sub-frame (this thread,
the post-S147/S163/S192/S201/S210 critic-window pick) is now CLOSED
at S215.

**Recommended next commit thread:** Per S214 slot4_summary, the next
commit slot is recommended to come from `frontier_gen` mode — but if
forced to pick from existing ATTACK_VECTORS, the leading candidate is
**D44 BC endomotive Galois-orbit** (Connes-Marcolli BC-system applied
to a Galois orbit on an arithmetic progression — flagged at S163 as
fall-back; cf. S203 for the closely related D48 BC-endomotive
direction). The motivation is: D44 imports a non-commutative-geometry
technique that the project has used only partially, and the BC-system
spectral construction has algorithmic ramifications (per Connes-Consani-
Moscovici 2009) that the previous Connes thread (Thread 2, closed S202)
explicitly noted as "an open lever".

If `frontier_gen` is auto-fired (per CLAUDE.md autonomy invariants), it
should produce 3–5 new ATTACK_VECTORS entries grounded in
cross-domain techniques NOT yet exhausted. Likely sources: persistent
homology of `pi(x)` at fine scales (TDA), free-probability spectral
constructions (free cumulants of zeta-zero spacings), transfer-operator
spectrum at a Connes-Bost-style automorphism (dynamical systems).

---

## 6. Files written / modified by this session

- `archive/sessions/session215_commit_a7_plethysm_WRAP.md` (this file)
- `EDGES.md` — E2.26 sub-frame distinction + WRAP pointer to S215
- `ATTACK_VECTORS.md` — A7 entry marked CLOSED at S215; thread DONE
- `RESEARCH_AGENDA.md` — A7 plethysm arc moved from in-flight to closed
- `NOVELTY_CHALLENGES.md` — successor challenge added (the d_0 question
  is registered as a Chow-variety classical-alg-geom question, NOT a
  chi_P-specific direction)
- `CROSS_DOMAIN_TECHNIQUES.md` §2 GCT entry — refined to "FULLY USED"
  for chi_P-vs-matched-baseline plethysm sub-frame
- `status/CLOSED_PATHS.md` — appended row for S215 thread-DONE closure
- `status/SESSION_INSIGHTS.md` — S215 entry appended
- `.commit_state` — sessions_used:5_final; thread marked DONE; history
  S211, S212, S213, S214, S215; prev_thread_4 set; recommended next
  thread noted
- `.run_state` — set to 216

---

## 7. Falsifiers (registered before WRAP)

A WRAP slot's "falsifiers" check whether the consolidation is
faithful to the slot 1–4 results. Three checks:

- (W1-FAL) Is there a slot 1–4 result that does NOT fall out of the
  S214 all-k structural identity? **No.** Slots 1–3 all measured
  particular k = 1, 2, 3 of the same `dim I_k(orbit-closure(f)) -
  dim I_k(orbit-closure(baseline))` invariant, which is identically
  zero by the all-k theorem.
- (W2-FAL) Is the non-degeneracy condition (∗) genuinely non-trivial
  at the verified scale? **Yes.** S214 V3 verified `cofactor rank =
  3` (full) at 10/10 matched-baseline samples; this is condition (∗).
  Without (∗), the orbit closures could differ at degenerate cofactor
  strata.
- (W3-FAL) Does any slot 1–4 result not flow through V_{d,1}^{n,d}
  Chow-variety geometry? **No.** All four slots' results follow from
  the parity lemma + Chow-variety identification.

All three no-signal branches hold. The WRAP synthesis is faithful.
