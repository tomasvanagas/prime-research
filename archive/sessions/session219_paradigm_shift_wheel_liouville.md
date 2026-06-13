# Session 219 — wheel-coprime Liouville-Möbius identity (PARADIGM-SHIFT mode)

**Date:** 2026-04-29.
**Mode:** PARADIGM-SHIFT (no cross-domain imports permitted).
**Self-grade:** **B**.
**Edges refined:** **E1.6** (bisection π = (x − L)/2 − C_3),
**E2.10** (L(x) mod 2 = x mod 2). **No new EDGES.md entries.**
**Successor of:** S55 (bisection / E1.6), S70 (q-stable bisection
extension), S198 (joint k-moduli E1.5 incrementality), S208
(`spike_divisor_only_wheel`, indicator-side `T_W^{div}` closed form).

## Mode discipline (paradigm-shift constraints)

* **No `WebFetch` / `WebSearch`.** No external lookups.
* **No new cross-domain technique.** Construction uses only Möbius
  inversion + completely-multiplicative property of `λ` — both
  already-`USED` in `CROSS_DOMAIN_TECHNIQUES.md` (multiplicative
  combinatorics, Möbius inversion).
* **No new `ATTACK_VECTORS.md` entries.**
* **Only existing project content read**: EDGES.md, recent
  `experiments/constructions/`, project session syntheses.

## Target picked

Compose `E1.6` (cumulative bisection π = (x − L)/2 − C_3) with
`E2.10` (L mod 2 = x mod 2) and the `E2.1 / E4.1 / S208` wheel-density
machinery (`φ(W)/W`, `T_W^{div}`) to produce a **NEW project artifact**:
an exact pointwise identity for the wheel-coprime cumulative Liouville
sum `L_W(x) := Σ_{n ≤ x, gcd(n, W) = 1} λ(n)` and a wheel-graded lift
of the bisection.

## Theorem trio (proven and integer-exact verified)

**Theorem 1** (wheel-coprime Liouville-Möbius identity).
For every `W ≥ 1` and integer `x ≥ 0`,
```
   L_W(x) = Σ_{d | rad(W)} L(⌊x/d⌋).
```
Right-hand side has exactly `2^{ω(W)}` terms regardless of how large
`W` is, depends on `W` only through `rad(W)`. Proof: Möbius expansion
of `[gcd(n, W) = 1]` + completely-multiplicative `λ(dm) = λ(d)λ(m)`,
collapsing `μ(d)λ(d) = +1` for squarefree `d` and `0` otherwise.

**Theorem 2** (parity of `L_W`).
`L_W(x) mod 2 = ( Σ_{d | rad(W)} ⌊x/d⌋ ) mod 2`.
Proof: Theorem 1 + E2.10. Computable in `O(ω(W))` integer ops, no
`L`-oracle needed.

**Theorem 3** (wheel-graded prime bisection — lift of E1.6).
```
   π_W(x) := π(x) − π_W^{div}(x)
            = (1/2) Σ_{d | rad(W)} (μ(d) ⌊x/d⌋ − L(⌊x/d⌋)) − C_{3,W}(x)
```
where `π_W^{div}(x) = #{p ≤ x : p | W}` is bounded and
`C_{3,W}(x) = #{n ≤ x : gcd(n,W)=1, composite, Ω(n) odd}`.

## What I built

`experiments/constructions/wheel_coprime_liouville_identity/`:
* `definition.md` — Theorems 1-3 with proofs + 6 pre-stated falsifiers.
* `wheel_coprime_liouville.py` — direct-enumeration sieve verifier.
* `wheel_coprime_liouville_results.json` — full per-cell JSON.
* `wheel_coprime_liouville_results.md` — F1-F6 verdict tables +
  algorithmic-implication summary.
* `run.log` — stdout of verification run.

## Outcome — F-verdicts

`N = 10⁴`, ~200 x-grid points spanning `[0, N]`, 16 distinct W
values across all three regimes (5 squarefree primorials, 4 squarefree
non-primorials, 7 non-squarefree). All checks **integer-exact**
(residual = 0, not "≤ 10⁻¹⁵"; both sides are integers).

| F# | Falsifier                                                  | Verdict |
|----|------------------------------------------------------------|---------|
| F1 | Pointwise Liouville identity (Theorem 1) at 16 W × ~200 x  | PASS, max |diff| = 0 |
| F2 | Radical reduction `L_W = L_rad(W)` at 6 non-sqfree W       | PASS, max |diff| = 0 |
| F3 | Mod-2 closed form (Theorem 2) at 16 W × ~200 x             | PASS, 0 mismatches |
| F4 | Wheel-graded prime bisection (Theorem 3) at 16 W × 4 x     | PASS, max |π_W − rhs| = 0 |
| F5 | Wheel-call invariance: 2^{ω(W)} terms regardless of W size | PASS at all 16 W |
| F6 | Mod-q lift for q ∈ {2, 3, 4, 8} at 5 W                     | PASS, 0 mismatches |

## Algorithmic content

Theorem 1 gives a deterministic `2^{ω(W)}`-call reduction
`L_W(x) → {L(⌊x/d⌋) : d | rad(W)}`. The reduction is:

* **Polylog in W** for primorial W (since `2^{ω(W)} = W^{O(1/log log W)}`,
  sub-polynomial in `W`).
* **Polylog in x iff `L(x)` is polylog** (the open frontier).
* **Mod-2 polylog unconditionally** (Theorem 2 + E2.10): `L_W(x) mod 2`
  is a sum of floor functions — `O(ω(W))` ops, no `L`-oracle needed.

The construction does NOT open a polylog route to `L(x)` itself.
Combined with E1.5 / T6 (CRT-mod-m saturation), the wheel-coprime
path does not bypass the unconditional `O(x^{2/3})` ceiling for
`L(x)`. The construction is **structural** — it cleanly separates
the wheel-coprime smooth side from the irreducible `C_{3,W}` residue.

## Why B-grade and not A

The mathematics (Möbius inversion of `[gcd(n,W) = 1]` + completely-
multiplicative `λ(dm) = λ(d)λ(m)`) is **classical** — derivable in
an afternoon by an analytic number theorist. Per CLAUDE.md grading
rubric:

> A-grade: "a published-paper-grade number theorist or complexity
> theorist could not derive in an afternoon from CLOSED_PATHS.md +
> EDGES.md alone."

This construction does not clear that bar. No new mathematical
content beyond a clean Möbius inversion that has appeared (in
slightly different form) in classical sieve theory.

## Why B-grade and not C

Three substantive new project artifacts emerge:

* **The exact pointwise identity** `L_W(x) = Σ_{d|rad(W)} L(⌊x/d⌋)` —
  **NOT** previously stated in EDGES.md, CLOSED_PATHS.md, any
  `experiments/`, or any session synthesis. A targeted search (greps
  for `L_W`, `L(x; W)`, `coprime to W` Liouville sums) returns no
  matches. The closest content is the **indicator-side** S208
  `T_W^{div}(n) = (π/N)·(W/φ(W))·[gcd(n,W)=1]` — a *different* but
  *symmetric* identity (cumulative-side here vs pointwise-side at S208).
* **The wheel-graded bisection lift** Theorem 3 — explicit divisor-
  sum closed form for `π_W(x)` modulo the irreducible `C_{3,W}`
  residue. E1.6's bisection has been used q-stable across q ∈ {3..13}
  at S70 (component independence at small q), but never written as a
  divisor-sum closed form lifted to all squarefree-radical wheels.
* **The 2^{ω(W)} oracle-call reduction** with the explicit primorial-
  W bound `2^{ω(W)} = W^{O(1/log log W)}`. Combined with E1.5 / T6
  (CRT-mod-m saturation), correctly diagnoses why the reduction does
  NOT bypass the polylog wall — the bottleneck is `L(x)` itself, not
  the wheel-coprime variant.

## Self-extension follow-ups (per CLAUDE.md)

1. **C8 — `L_W` minus first-order Selberg-Delange smooth term**.
   Asymptotically `L(x) = o(x)` (PNT for λ); `L_W(x)/L(x) → φ(W)/W` is
   a folklore heuristic. Test the **exact identity** Theorem 1 against
   this asymptotic:
   `L_W(x) − (φ(W)/W) L(x) = Σ_{d | rad(W), d > 1} (L(⌊x/d⌋) − ⌊x/d⌋·L(x)/x)`
   gives an exact decomposition of the heuristic error term into
   2^{ω(W)} − 1 explicit pieces. Estimate the error magnitude
   empirically at scale.
2. **C8.b — Mod-4 lift of `L_W`**. The mod-2 closed form (Theorem 2)
   is trivial via E2.10. Mod 4 is non-trivial: `L(x) mod 4` is on the
   E1.5 / T6 hard-zone bit boundary. Test whether `L_W(x) mod 4 =
   Σ_{d | rad(W)} L(⌊x/d⌋) mod 4` gives any algorithmic leverage —
   probably not (E1.5 saturation), but worth empirical confirmation
   at small scale.

Both successor challenges added to NOVELTY_CHALLENGES.md as §C8 / §C8.b.

## Edges modified

* **EDGES.md E1.6** — annotated inline with Theorem 3 wheel-graded
  lift + Theorem 1 cumulative-Liouville identity (S219 refinement).
* **EDGES.md E2.10** — annotated with Theorem 2 closed form for
  `L_W(x) mod 2` (S219 refinement).
* **No new EDGES.md entry**: classical-derivation criterion
  (afternoon-derivable by NT theorist) places this as a refinement
  of two existing edges, not as a new edge.

## CLOSED_PATHS.md row

Added row at S219 documenting the wheel-coprime cumulative-Liouville
identity as a B-grade refinement of E1.6 / E2.10.

## What I produced that wasn't in the project before this session

(per CLAUDE.md self-evaluation Q1)

1. The pointwise identity `L_W(x) = Σ_{d|rad(W)} L(⌊x/d⌋)` — verified
   integer-exactly at 16 W × ~200 x.
2. The wheel-graded bisection lift Theorem 3 — verified at 16 W × 4 x.
3. The mod-2 closed form Theorem 2 — verified.
4. The radical reduction `L_W = L_{rad(W)}` — verified at 6 non-
   squarefree W.
5. The `2^{ω(W)}` oracle-call reduction with the explicit primorial-W
   bound `2^{ω(W)} = W^{O(1/log log W)}`.

(per CLAUDE.md self-evaluation Q2): edges composed = E1.6 + E2.10 +
E2.1/E4.1/S208 (wheel-density / indicator-side), with E1.5/T6 framing
the algorithmic-content boundary.

(per CLAUDE.md self-evaluation Q4): next-action for the next agent
is C8 (Selberg-Delange decomposition of the L_W − (φ(W)/W) L(x)
heuristic error term) — added to NOVELTY_CHALLENGES.md.

## Files

`experiments/constructions/wheel_coprime_liouville_identity/{
  definition.md,
  wheel_coprime_liouville.py,
  wheel_coprime_liouville_results.md,
  wheel_coprime_liouville_results.json,
  run.log
}`.

EDGES.md E1.6 and E2.10 — inline S219 annotations.
NOVELTY_CHALLENGES.md §C8 + §C8.b — successor challenges.
CLOSED_PATHS.md S219 row — refinement closure.
SESSION_INSIGHTS.md — this session's entry.
.run_state — bumped to 220.

## Self-evaluation

* **Q1 (what was new):** five new artifacts above.
* **Q2 (edges composed):** E1.6, E2.10, E2.1/E4.1, S208, with E1.5/T6
  as algorithmic-content boundary.
* **Q3 (only duplicate closures? if so why):** no — actual integer-
  exact pointwise identity and bisection lift, with non-trivial
  algorithmic content (`2^{ω(W)}`-call reduction).
* **Q4 (next-action):** C8 / C8.b in NOVELTY_CHALLENGES.md.
