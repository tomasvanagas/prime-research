# Session 427 — Möbius bisection of π(x) (paradigm-shift mode)

**Mode.** Paradigm-shift — no cross-domain technique imported, no
WebFetch/WebSearch, no new ATTACK_VECTORS entries. Pure recombination
of edges already in the project plus elementary Dirichlet-series
identity manipulation.

**Edges composed.** E2.2 (Liouville bisection
`π = (x − L)/2 − C_3`), E1.6 / E2.10 (parity bisection,
`L(x) mod 2 = x mod 2` trap), folklore Möbius decomposition
`M = S_e − S_o, Q = S_e + S_o`.

## Self-grade — B

(Refinement of E2.2 with a precise new statement that extends its
scope. Honest demotion to C-grade DUPLICATE-PLUS in a verify session
is defensible — see "what this is not" below.)

## Q1: What did this session produce that was not in the project before?

Four pieces of project content that did not appear in `EDGES.md`,
`CLOSED_PATHS.md`, `NOVELTY_CHALLENGES.md` or `archive/sessions/`
before this session:

1. **The Möbius bisection of π(x)**: for `Q(x) = #{n ≤ x : sqfree}`,
   `M(x) = Σμ(n)`,
   `C_3*(x) = #{n ≤ x : sqfree, composite, ω odd, ω ≥ 3}`,
   ```
       π(x) = (Q(x) − M(x))/2 − C_3*(x)
   ```
   integer-exact for all `x ∈ [1, 10⁶]` (paradigm-shift discipline:
   only x ≤ 10⁶ verified; sieve scales as `Õ(N log log N)` so larger
   `N` is rerunnable).

2. **The bridge identity** linking the Liouville and Möbius bisections:
   ```
       NS_o(x) := #{n ≤ x : not sqfree, Ω odd}
                = (x − Q(x) − L(x) + M(x)) / 2.
   ```
   One-line analytic proof from `Σ_{n ≤ x} μ²(n) λ(n) = Σ_{n ≤ x} μ(n)
   = M(x)`, because `λ = μ` on the squarefree subset and `μ² = 0` off
   it. Verified integer-exact at every `x ∈ [1, 10⁶]`.

3. **The 4-cell decomposition** `(squarefree × Ω-parity)` with explicit
   closed forms for every cell:
   ```
       S_e = (Q + M)/2,  S_o = (Q − M)/2,
       NS_e = ((x − Q) + (L − M))/2,  NS_o = (x − Q − L + M)/2.
   ```

4. **Empirical Möbius-side parity independence**:
   `I(B(x) mod 2 ; C_3*(x) mod 2) ≈ 1.15 × 10⁻⁵` bits at `N = 10⁶`,
   the same order as E1.6's `I(A ; C_3) ≈ 2 × 10⁻⁵` bits. The Möbius
   parity bisection inherits E1.6's near-independence.

## Q2: What edges did the work compose or cite?

* **E2.2** (Liouville bisection) — directly extended with the parallel
  Möbius identity and the bridge.
* **E1.6** (parity bisection of π mod 2 with `A`/`C_3` near-independent)
  — Möbius-side analogue measured.
* **E2.10** (`L(x) mod 2 = x mod 2` trap) — used for the `(x − L)/2`
  integrality and parity convention.
* The folklore identity `Q(x) − M(x) = 2·#{sqfree, ω odd}` (Mertens
  1874-style splitting).

Inline annotations added to E2.2 in `EDGES.md` (S427 block), and a row
appended to CLOSED_PATHS.md "Sieve / Combinatorial / Counting" section.

## Q3: If only duplicate closures, why?

This session is not a duplicate-only closure. It is a **structural
refinement of E2.2** that adds a parallel bisection plus a closed-form
bridge identity, none of which were recorded in the project. However,
honesty requires noting that the construction is a textbook-level
identity rearrangement — a working analytic NT theorist would derive
this in well under an hour. The "novelty" is project-catalogue,
not discipline-of-mathematics. The B-grade self-claim is the highest
defensible grade and may be demoted in verify.

The session **did succeed at the paradigm-shift mode discipline**: no
external technique imported, no new ATTACK_VECTORS, no WebFetch. The
output is a recombination of existing project edges with elementary
Dirichlet-convolution identity moves.

## Q4: Next-action

Two paradigm-shift-mode successors are recorded in
`mobius_bisection_results.md` §"Successor proposals" and offered to
NOVELTY_CHALLENGES (not added in this session, to keep the
self-extension pressure clean):

* **C13.a** Wheel-W-coprime generalisation:
  `π_W(x) = (Q_W(x) − M_W(x))/2 − C_{3,W}*(x)` and the corresponding
  bridge `NS_o,W(x) = (x_W − Q_W − L_W + M_W)/2`. Composes with the
  S219 wheel-graded Liouville bisection lift.
* **C13.b** Test prime-q joint distribution `(A mod q, B mod q)` for
  `q ∈ {3, 5, 7, 11, 13}` (parallel of S70 / C1 on the Möbius side).

Both are single-session B-target experiments. Neither requires an
external technique.

## Empirical numbers at N = 10⁶ (wall 6.27 s)

| symbol     | value      | predicted lag-1 ac (1−2·incr-density) | empirical |
|------------|-----------:|--------------------------------------:|----------:|
| `π(N)`     |     78498  |                                       |           |
| `Q(N)`     |    607926  | 0 (every step) [N/A]                  |           |
| `M(N)`     |       212  |                                       |           |
| `L(N)`     |      −530  |                                       |           |
| `C_3*(N)`  |    225359  | 0.550                                 |   0.549   |
| `C_3(N)`   |    421767  | 0.156                                 |   0.156   |
| `NS_o(N)`  |    196408  | 0.608                                 |   0.607   |
| `S_o(N)`   |    303857  |                                       |           |
| **`A` par**  |          | 0.000                                 |  −0.001  |
| **`B` par**  |          | 0.392                                 |   0.392   |
| **`π` par**  |          | 0.843                                 |   0.843   |

All lag-1 autocorrelations match the cumulative-counter density-only
prediction to within 1 % — the parities are **density-structured but
not predictability-structured**, exactly the state of E1.6 transposed
to the Möbius side.

Mutual-information matrix at `N = 10⁶`:

| pair                   | bits         |
|------------------------|--------------|
| `I(A par ; B par)`     | 1.09 × 10⁻⁶  |
| `I(C₃ par ; C₃* par)`  | 1.09 × 10⁻⁶  |
| `I(A par ; C₃ par)`    | 1.15 × 10⁻⁵  |
| `I(B par ; C₃* par)`   | 1.15 × 10⁻⁵  |
| `I(π par ; NS_o par)`  | 2.89 × 10⁻⁷  |

## What this is NOT

* **Not algorithmic**. `(Q − M)/2 ~ 0.304·x` and `C_3* ~ 0.304·x` with
  `π ~ x/log x` as the needle — same C-circular obstacle as E2.2.
* **Not deep mathematics**. Textbook Möbius/squarefree partition,
  one-line proof of the bridge identity. The derivation is well below
  the published-paper threshold.
* **Not falsifying any open conjecture**. The four-cell decomposition
  is a definitional consequence; the bridge identity follows by direct
  Dirichlet calculation; the Möbius bisection is a re-arrangement of
  these.

## Verify-mode hooks

* **Verify F1–F6 at larger N** (the script supports any `N`; at
  N = 10⁷ wall ≈ 60 s estimated linear-extrapolating from 6.27 s).
* **Re-derive the bridge identity from the indicator
  `(1 − μ²)(1 − λ)/2`** as cross-check (matches the script output
  by definition; verifier should attempt an INDEPENDENT derivation
  via Dirichlet generating functions).
* **Demote to C-grade DUPLICATE-PLUS** if the verifier judges the
  Möbius bisection to be a textbook re-arrangement rather than a
  genuine refinement of E2.2's scope. That demotion would be
  defensible. The B-grade claim rests on the bridge identity and the
  4-cell decomposition being explicit closed forms not previously
  recorded.

## Cleanup

* `experiments/constructions/mobius_bisection_of_pi/` —
  `mobius_bisection.py`, `mobius_bisection_results.md`,
  `definition.md`, `mobius_bisection_results.json`.
* `EDGES.md` — S427 inline block under E2.2.
* `status/CLOSED_PATHS.md` — S427 row under "Sieve / Combinatorial /
  Counting".
* `archive/sessions/session427_mobius_bisection.md` — this file.
