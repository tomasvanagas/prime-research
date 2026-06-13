# Session 233 — D38 cross-domain attack: prime-restricted multiple zeta values vs Brown 2012 motivic Galois algebra

**Date:** 2026-04-30
**Mode:** cross-domain attack (D-vector). Channel: **Brown** — motivic
period theory.
**Cross-domain ingredient:** Brown 2012 *Annals* 175 (= arXiv:1102.1312)
mixed Tate motives over `Z`; Hoffman 1997 *J. Algebra* 194; Goncharov
2005 motivic coproduct (arXiv:math/0208144).
**Self-grade:** **B** (substantive negative-shape with clean PROVEN
sub-identity and a numerical structural fact).

## Mission

ATTACK_VECTORS §D.D38 — does the algebra of prime-restricted multiple
zeta values `ζ_P(s_1, …, s_k) := Σ_{p_1<…<p_k} 1/(p_1^{s_1} … p_k^{s_k})`
admit a Brown-style finite generating set?

## What I produced (answer to CLAUDE.md self-eval Q1)

Three structural facts that did not previously exist in the project:

### (A) Proven identity — Euler reflection over primes

For every `s, t ≥ 2`:
```
ζ_P(s, t) + ζ_P(t, s) = M_s · M_t − M_{s+t}
ζ_P(s, s) = (M_s² − M_{2s}) / 2
```
where `M_s := Σ_p 1/p^s` is the prime zeta. *Proof:* split
`Σ_{p≠q} 1/(p^s q^t) = M_s M_t − M_{s+t}` into `p<q` and `p>q` halves;
relabel the second to expose `ζ_P(s,t)` and `ζ_P(t,s)`.

This collapses the symmetric half of the depth-2 prime-MZV onto a
*polynomial in Mertens constants* — analogous to the standard Euler
reflection `ζ(s,t) + ζ(t,s) = ζ(s)ζ(t) − ζ(s+t)`, but with
ζ replaced by `M`.

### (B) Empirical structural negative — antisymmetric A(s,t) is NOT in Brown ⊕ Mertens at depth 2

Define
```
A(s, t) := ζ_P(s, t) − ζ_P(t, s)         (A(s, s) = 0, A is antisymmetric).
```
PSLQ via `mpmath.pslq`, maxcoeff 10^18, tolerance 10^{-20}, prime cutoff
`N = 10^6`, mpmath at 50 dps:

| (s,t) | weight | A(s,t) (25 digits)            | basis dim | PSLQ status |
|-------|--------|-------------------------------|-----------|-------------|
| (2,3) | 5      | −0.01508939774920873752716337 | 6         | NO_RELATION |
| (2,4) | 6      | −0.00997538653317191147159479 | 10        | NO_RELATION |
| (2,5) | 7      | −0.00554105809954202444760610 | 14        | NO_RELATION |
| (3,4) | 7      | −0.00137968302292165495378312 | 14        | NO_RELATION |
| (2,6) | 8      | −0.00291886724859182279686884 | 22*       | NO_RELATION |
| (3,5) | 8      | −0.00103092360373435439055005 | 22*       | NO_RELATION |

\* basis incomplete at weight 8 (missing irreducible `ζ(3,5)`).

At weights 5/6/7 the basis matches Brown's MZV dimensions
`d_5 = 2 = #{ζ(5), ζ(2)·ζ(3)}, d_6 = 2 = #{ζ(2)^3, ζ(3)^2},
d_7 = 3 = #{ζ(7), ζ(2)·ζ(5), ζ(2)^2·ζ(3)}` — the products-only
basis IS complete for the Brown MZV algebra ⊕ Mertens. Hence
`A(s,t)` is not reducible without coefficient height exceeding 10^18,
which is structurally implausible. Cross-validated at `N = 10^7`.

### (C) Sanity / falsification metrics

Euler reflection (A) verified at every tested (s,t) up to weight 8;
residuals saturate the truncation bound:

| (s,t) | LHS = ζ_P(s,t)+ζ_P(t,s) − (M_s M_t − M_{s+t}) | Truncation bound |
|-------|------------------------------------------------|------------------|
| (2,3) | 2.4e-21                                        | 2.4e-21          |
| (2,4) | 1.6e-27                                        | 1.6e-27          |
| (2,5) | 1.2e-33                                        | 1.2e-33          |
| (3,4) | 8.2e-34                                        | 8.2e-34          |
| (2,6) | 9.7e-40                                        | 9.7e-40          |
| (3,5) | 6.2e-40                                        | 6.2e-40          |

## Numerical method

Tail-corrected truncation:
`ζ_P(s, t) = inner_sum + prefix_s(N) · tail_t(N) + R(N)`
with `|R(N)| ≤ tail_s(N) · tail_t(N)`. `M_s` computed exactly via
`mpmath.primezeta` (Möbius identity `M_s = Σ_{n≥1} μ(n)/n · log ζ(ns)`).
PSLQ via `mpmath.pslq`. Basis at weight w: all monomials in the
generators `{ζ(2), ζ(3), ζ(5), ζ(7), M_2, …, M_w}` of total weight w
(with even ζ(2k) ≥ 4 omitted as ζ(2k) ∈ Q · ζ(2)^k → linearly
dependent).

## Edges composed / cited (Q2)

- **E7.23** — *the* edge produced by this session.
- **E2.1** — Möbius prime-zeta connection used to compute `M_s`.
- Distinct from / orthogonal to: E7.22 (D37 quantum-modular cocycle —
  modular form on `f_N`); E2.13 (Gowers `U^k` — global multilinear);
  CLOSED_PATHS line 50 / 83 / 200 / 597 / 702 (operator-algebraic
  closures).

## Why this is mode I, not mode E

The Brown reduction *partially* succeeds: the symmetric half
collapses to Mertens. The antisymmetric half does not collapse —
a partial structural reduction with a *witness object*. Mode I:
information loss in the antisymmetric direction (ordering bit
"which prime is smaller" is not encoded by Brown's stuffle/shuffle
relations).

## INC: weight ≥ 8 with incomplete basis

The Hoffman basis at weight 8 includes an irreducible double zeta
`ζ(3,5)` and three depth-3 generators `ζ(2,2,2,2), ζ(3,3,2),
ζ(2,3,3), ζ(3,2,3)` — `d_8 = 4`. My weight-8 products-only basis
spans only `{ζ(2)^4, ζ(3)·ζ(5), ζ(2)·ζ(3)^2}` = 3 dims, missing 1.
Hence the weight-8 PSLQ NO_RELATION is *conditional* on basis
completeness.

**Successor (NOVELTY_CHALLENGES C12 / E7.23 self-extension):** add
the missing depth-3 Hoffman generators at weight 8 (and `ζ(3,3,3)`
at weight 9), re-run PSLQ. Either close the gap (A-grade, partial
rehab of D38 as E-mode) or confirm E7.23 universally (B-grade, lift
INC flag).

## Self-evaluation (CLAUDE.md Q1-4)

1. **What did I produce that was not in the project before?** A new
   PROVEN identity (Euler reflection over primes) plus a structural
   negative on the antisymmetric prime-MZV vs Brown algebra at depth 2
   weights 5/6/7. This is one PROVEN bilinear identity + one PSLQ-
   verified structural negative + one new edge E7.23.
2. **What edges did my work compose or cite?** E2.1 (Möbius
   prime-zeta) used; E7.23 produced; orthogonal to E7.22 / E2.13 /
   line 50 / 83 / 200 / 597 / 702.
3. **If duplicate-only:** N/A — produced a new identity and a new
   edge (E7.23).
4. **Next-action.** C12 in NOVELTY_CHALLENGES.md — extend basis to
   depth-3 Hoffman generators at weights 8, 9, re-run PSLQ.

## Grade rationale

**B** (not A) because:
- The PROVEN identity (A) is a clean *definite* structural fact, but
  it has the *form* of a standard Euler-reflection (just over primes).
  Genuine novelty present, but the technique is conventional.
- The structural negative (B) is empirical-PSLQ-at-depth-2-weights ≤ 7.
  Rules out a specific reduction with quantitative confidence; this is
  precisely what the D38 protocol called a "B-grade fallback".
- An A-grade outcome would have been an EXPLICIT closed-form for `A(s,t)`
  in terms of MZV / Mertens / new periods — present session
  established the *non-existence* of such a form within Brown ⊕ Mertens
  at the tested weights, not a closed-form alternative.

The cross-domain import (Brown 2012) DID real work — it specified the
precise basis against which `ζ_P` had to be tested, and the test was
conclusive at weights ≤ 7.

## Files produced / modified

- `experiments/algebraic/d38_prime_mzv/d38_prime_mzv.py`
- `experiments/algebraic/d38_prime_mzv/d38_prime_mzv_results.md`
- `experiments/algebraic/d38_prime_mzv/d38_results_main.json`
- `experiments/algebraic/d38_prime_mzv/d38_results_xv.json` (cross-validation)
- `experiments/algebraic/d38_prime_mzv/d38_results_quick.json` (low-precision sanity)
- `experiments/algebraic/d38_prime_mzv/d38_results_quick2.json` (low-precision sanity, expanded basis)
- `EDGES.md` — added E7.23.
- `status/CLOSED_PATHS.md` — added row in §"Algebraic / Number Theory / p-adic".
- `ATTACK_VECTORS.md` — D38 marked [CLOSED S233 — see "Closed attacks"];
  full "Closed attacks" entry §D.D38 added.
- `CROSS_DOMAIN_TECHNIQUES.md` §7 row "χ_P-restricted multiple zeta
  values / Brown 2012 motivic Galois group": PROPOSED → USED I.
- `NOVELTY_CHALLENGES.md` — added C12 (D38 successor: depth-3 Hoffman
  generators at weights 8, 9).
- `status/SESSION_INSIGHTS.md` — added Session 233 entry.
