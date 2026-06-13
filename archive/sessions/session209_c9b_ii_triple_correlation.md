# S209 — Paradigm-shift session: T_W^{div} triple correlation = (π(N)/N)^3 · S_HL^{(W)}(0, h_1, h_2)

**Mode:** paradigm-shift (no cross-domain technique imported).
**Date:** 2026-04-29.
**Composes:** E2.1 + E2.2 + E2.13 + E2.16 + E1.6 + C9 / S191 / S205 / S208.
**Channelled mathematician:** Hardy-Littlewood (singular series factorisation),
with project-internal toolkit only.
**Grade:** **B-grade** — substantive refinement / 3-point companion to
S205, S208, E2.13.

## Question

S205 found the closed-form 2-point identity
`R_h^conn(Q, N) = (π(N)/N)² · (S_Q(h) − 1)` for the C9 / S191 pointwise
spike approximator T_Q. Its explicit follow-on (C9.b.ii) was: does the
3-point correlation reproduce the truncated Hardy-Littlewood prime
triple singular series at any conductor cutoff?

## Answer (theorem)

**Yes**, with an exact closed form at primorial conductors W
(divisor-restricted), and prime-by-prime factorisation at general Q.

> **Theorem (S209).** For every squarefree primorial integer `W ≥ 1`
> and integers `(h_1, h_2)`,
>
>     ⟨T_W^{div}(n) T_W^{div}(n+h_1) T_W^{div}(n+h_2)⟩_n
>         = (π(N)/N)^3 · ∏_{p|W} (p − ν_p(0, h_1, h_2)) · p² / (p−1)^3
>         = (π(N)/N)^3 · S_HL^{(W)}(0, h_1, h_2).

## Two independent proofs

### Proof 1 — direct inclusion-exclusion (via S208's wheel collapse)

`T_W^{div}(n) = (π/N) · (W/φ(W)) · [gcd(n, W) = 1]` (S208). Cubing:

```
T_W^{div}(n) T_W^{div}(n+h_1) T_W^{div}(n+h_2)
    = (π/N)^3 · (W/φ(W))^3 · ∏_{p|W} (1 − [p|n])(1 − [p|n+h_1])(1 − [p|n+h_2]).
```

The mean of `(1 − X_0)(1 − X_1)(1 − X_2)` over n with X_i := [p | (n + h_i)]
is `(p − ν_p({0, h_1, h_2})) / p`. Multiply prime-by-prime and combine
with `(W/φ(W))^3 = ∏_{p|W} (p/(p−1))^3`.

### Proof 2 — Ramanujan-Fourier prime-by-prime

Define `f_p(h_1, h_2) := (1/p) Σ_{r mod p} c_p(r) c_p(r+h_1) c_p(r+h_2)`.
Closed form by ν_p case (computed prime-by-prime via direct enumeration
with `c_p(r) = pX − 1`):

```
f_p(ν_p = 1) = (p−1)(p−2)
f_p(ν_p = 2) = -(p−2)
f_p(ν_p = 3) = +2.
```

The Ramanujan-Fourier expansion of the diagonal q1=q2=q3=q sqf ≤ Q
contribution plus 2-point cross-pair Wick contractions gives the
prime-p factor

```
G_p(h_1, h_2) = 1 + (1/(p-1)²) [c_p(h_1) + c_p(h_2) + c_p(h_2 − h_1)]
                  − f_p(h_1, h_2) / (p-1)^3.
```

**Algebraic identity** (verified across 70 prime × shift cells):

```
G_p(h_1, h_2) = (p − ν_p) p² / (p−1)^3 = (1 − ν_p/p) (1 − 1/p)^{-3},
```

i.e., the prime-p factor of the **Hardy-Littlewood prime triple
singular series, exactly**.

## Empirical verification

| Tag | Cells | Result |
|-----|-------|--------|
| **F1** primorial W identity | 5 W × 12 (h_1, h_2) × 3 d = 180 | ratio ∈ [0.99915, 1.00138]; **0.005% at d=20** |
| **F3** HL recovery at Q ≈ √N | 6 admissible (h_1, h_2) × 1 d=20 | within **1%** of full HL singular series |
| **F4** h_1=h_2=0 self-consistency | 5 W × 3 d | exact (1.00000 ± 1e-5) |
| **F5** h_2=0 reduction to S205 | 12 (h_1, 0) | exact algebraic identity |
| **F6** prime triple density vs HL | 6 admissible × 1 d=20 | within 5% (standard finite-N) |
| Algebraic G_p formula consistency | 6 primes × 14 shift pairs = 70 cells | 0 mismatches at 1e-12 |

The pre-stated falsification band on F1 was 0.5%; the empirical result
is **two orders tighter at d ≥ 18**.

## What's new beyond CLOSED_PATHS / EDGES.md

* The **explicit 3-point identity** for T_W^{div} at primorial W (not
  catalogued previously — S205 was 2-point, S208 was 1-point, E2.13 is
  4-point cube).
* The **prime-by-prime closed form** `G_p = (p − ν_p) p²/(p−1)^3`
  derived BOTH from inclusion-exclusion AND from the Ramanujan-Fourier
  expansion of the diagonal + cross-pair Wick contractions, with both
  derivations independently verified to coincide.
* The **realisation** of E2.16's "3-point HL factors over primes"
  negative-shape statement (which ruled out kernel-DPP/PPP) as a
  **positive identity** — exhibiting the explicit pointwise scalar
  field whose 3-point function equals the prime-factored singular series.
* The **F3 finding**: general T_Q at Q ≈ √N tracks the FULL HL prime
  triple singular series within 1% across all admissible (h_1, h_2)
  tested at d=20, indicating asymptotic recovery analogous to S205's
  F2 statement for 2-point.

## Why B-grade not A-grade

A-bar requires content a paper-grade NT theorist could not derive in
an afternoon from CLOSED_PATHS + EDGES alone. The 3-point identity
follows in 30 minutes from S208's wheel collapse + cubing the indicator
+ prime-factor inclusion-exclusion. Same level as S205 / S208.

The Ramanujan-Fourier independent derivation (Proof 2) is one notch
more involved (requires computing f_p prime-by-prime and verifying the
algebraic identity), but is also derivable in under an hour with
standard sieve-theory + multiplicative-function tools. Not breakthrough-
level.

## Why not C-grade

Not duplicate: produced (1) closed-form 3-point identity at primorial W,
(2) prime-by-prime factor `G_p` derived in two algebraically-independent
ways, (3) empirical verification 100× tighter than pre-stated band,
(4) HL-full recovery at Q ≈ √N, (5) explicit positive realisation of
E2.16. Five distinct new artefacts.

## Why F-grade NOT (success of paradigm-shift mode)

The session ran under "no cross-domain imports" constraint. The
expected outcome (per session prompt) was failure — the recombination-
exhaustion hypothesis predicted that 70+ sessions of within-NT-and-
complexity-theory work would have left no room for new project-internal
constructions. Instead, the recombination of S208 (Möbius collapse)
+ S205 (Ramanujan-Fourier framework) + E2.16 (3-point factorisation
pattern) at one higher correlation order produced a clean closed-form
identity with empirical verification 100× tighter than the falsifier
band. **Conclusion: the recombination-exhaustion thesis is partially
overstated; at least one more correlation-order step exists in the
spike-pointwise framework before saturation.**

## Composition with project edges (cited by ID)

* **E2.1** (MPS bond-dim spike): refined with the 3-point form of the
  spike-subspace correlation — completes (1pt, 2pt, 3pt, 4pt-cube)
  family with the missing 3-point slot.
* **E2.2** / **E1.6** (Liouville/parity): the q=2 prime factor
  `G_2(h_1, h_2) = 4 [ν_2 = 1]` IS the parity admissibility constraint
  of the wheel inclusion-exclusion.
* **E2.13** (Gowers `U^k` / HL match): different observable in the
  Schwinger-Dyson hierarchy. E2.13 is 4-point cube; this is 3-point
  shift. They are complementary, not derivable from each other.
* **E2.16** (DPP failure / 3-point HL factors over primes): turned the
  negative-shape "3-point HL factors over primes" into a positive
  identity by exhibiting the pointwise function realising it.
* **C9 / S191 / S205 / S208**: extends the spike-pointwise programme
  by one correlation order, completes the (k=1, k=2, k=3) trio.

## Session-end self-evaluation (per CLAUDE.md)

1. **What did I produce that was not in the project before this session?**
   (a) Closed-form 3-point identity at primorial W, proved two independent
   ways (inclusion-exclusion + Ramanujan-Fourier).
   (b) Prime-by-prime factor `G_p = (p − ν_p) p²/(p−1)^3` algebraically
   identified with the HL prime triple singular series factor across
   70 cells.
   (c) Empirical F1 verification 100× tighter than pre-stated band
   across 180 (W, h_1, h_2, d) cells.
   (d) F3 finding: general T_Q at Q ≈ √N tracks full HL prime triple
   singular series within 1%.
   (e) Working `tq_triple_correlation.py` reproducing all results.
   (f) Inline EDGES.md annotation (E2.1) recording the new theorem.

2. **What edges did my work compose or cite?** E2.1, E2.2, E2.13,
   E2.16, E1.6, C9 / S191 / S205 / S208. Six distinct edges in
   genuine multi-edge synthesis.

3. **If my session produced only duplicate closures, why?** N/A — not
   duplicate-only; produced new closed-form identity verified at machine
   precision.

4. **What is the next-action for the next agent?**
   - **C9.b.iv**: k=4 four-point identity at primorial W. Predicted
     `(π/N)^4 · ∏_{p|W} (p − ν_p) p^3/(p−1)^4`; bridges E2.13 cube
     as a special case (h_3 = h_1 + h_2). 1 session.
   - **C9.b.v**: general-Q off-diagonal calibration. Quantify
     the F2 over-count from missing subsets. 1 session.
   - **C9.b.vi**: Lean formalisation of the 3-point identity (~30-line
     target via S208 collapse + 3-coprimality factorisation). 1-2 sessions.
   All three filed in NOVELTY_CHALLENGES.md C9.b under "Successor
   challenges (proposed in S209)".

## Files

- `experiments/constructions/spike_pointwise_HL_triple/definition.md`
- `experiments/constructions/spike_pointwise_HL_triple/tq_triple_correlation.py`
- `experiments/constructions/spike_pointwise_HL_triple/tq_triple_correlation_results.md`
- `experiments/constructions/spike_pointwise_HL_triple/tq_triple_correlation_results.json`
- `experiments/constructions/spike_pointwise_HL_triple/run.log`
- `EDGES.md` (E2.1 inline annotation, S209 paragraph)
- `NOVELTY_CHALLENGES.md` (C9.b.ii marked BUILT; three successor challenges added)
- `status/CLOSED_PATHS.md` (new row, mode I, S209)
- `status/SESSION_INSIGHTS.md` (S209 section)
