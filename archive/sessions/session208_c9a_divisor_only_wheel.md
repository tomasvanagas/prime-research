# Session 208 — C9.a divisor-only restriction of M_Q at squarefree W

**Mode**: novelty (B-grade target).
**Target ID**: §C9.a in `NOVELTY_CHALLENGES.md`.
**Edges composed**: E2.1 + E2.2 + S191 / C9.

## Selection rationale

After the S202 commit-thread closures (S82 / Connes / Galway), the
`NOVELTY_CHALLENGES.md` queue holds C9.a as a clean composition target
with a pre-stated closed form. The challenge spec (line 572 of
`NOVELTY_CHALLENGES.md`) flagged the predicted form as "or similar
exact-divisor form" — i.e. the spec author was uncertain about the
exact closed form. The work of this session is to (a) derive the
correct form analytically, (b) verify it pointwise to machine epsilon,
and (c) quantify the gap between the divisor-only restriction and the
full S191 `T_Q`.

## The closed-form theorem

**Theorem (S208).** For every squarefree integer `W ≥ 1` and every
positive integer `n`,

```
    M_W^{div}(n)  :=  Σ_{q | W, q sqf} mu(gcd(q,n)) / phi(q/gcd(q,n))
                  =   [gcd(n, W) = 1] · W / phi(W).
```

For non-squarefree `W`, the identity holds with `W` replaced by
`rad(W)`.

**Proof.** Index squarefree divisors `q | W` by subsets `S` of the
prime support `P(W)`. Decompose `S = A ⊔ B` with `A ⊆ P_n := {p|W:p|n}`
and `B ⊆ P_{!n}`. Then `gcd(q_S, n) = ∏_A` and `q_S/gcd(q_S, n) = ∏_B`,
so the summand is `(-1)^|A| / ∏_{p∈B}(p-1)`. The double sum factors:

```
    M_W^{div}(n) = ( Σ_{A ⊆ P_n} (-1)^|A| ) · ( Σ_{B ⊆ P_{!n}} ∏_{p∈B} 1/(p-1) ).
```

The first factor is `(1 − 1)^|P_n| = [P_n = ∅]`. The second factor is
`∏_{p ∈ P_{!n}} (1 + 1/(p-1)) = ∏_{p ∈ P_{!n}} p/(p-1)`. When
`P_n = ∅` (i.e. `gcd(n, W) = 1`), `P_{!n} = P(W)` and the second
factor is `W/φ(W)`. Done.                                              □

## Empirical verification

`spike_divisor_only.py` directly computes `M_W^{div}(n)` from the
definition (sum over the explicit list of squarefree divisors of `W`)
and compares against the closed form. Results across 10 conductors
`W ∈ {2, 6, 30, 210, 2310, 15, 105, 12, 60, 420}` and `N ∈ {10⁵, 10⁶}`:

* **F1 / F2 / F3 pointwise identity**: PASS, max abs_err `8.88·10⁻¹⁶`
  (machine epsilon).
* **F4 mean = 1**: PASS, all means in `[0.999997, 1.000000]` (the
  slight under-shot is the boundary effect of the few small primes
  `p ≤ rad(W)` that themselves divide `W`).
* **F5 L² closed form**: PASS, rel_err `< 10⁻⁹`. The explicit
  prime-distribution count
  `(n_coprime_primes (T_val − 1)² + n_coprime_comp T_val² + n_W) / N`
  (with `T_val := (π(N)/N) · W/φ(W)`) matches the empirical
  `||T_W^{div} − chi_P||²/N` to nine decimal places.
* **F6 Pearson gap**: PASS for `W ∈ {6, 30, 210, 2310}`. The gap
  `Pearson(chi_P, T_W) − Pearson(chi_P, T_W^{div})` at `N = 10⁵`:

  | W       | Pearson_div | Pearson_full | Gap   |
  |---------|-------------|--------------|-------|
  | 2       | 0.32566     | 0.32566      | 0.000 |
  | 6       | 0.46050     | 0.48843      | 0.028 |
  | 30      | 0.53993     | 0.62889      | 0.089 |
  | 210     | 0.59807     | 0.77511      | 0.177 |
  | 2310    | 0.63558     | 0.91622      | 0.281 |

  The W=2 equality is structural (squarefree integers ≤ 2 are exactly
  the divisors of 2). The gap grows monotonically, and `T_W^{div}`
  captures only `0.636 / 0.916 = 69%` of `T_W`'s Pearson at W=2310.

## Why this composition is novel relative to S191

S191 stated the divisor-only result only as a corollary at coprime `n`
(`T_W | gcd(n,W) = 1 = π(N)/N · W/φ(W)`, line 213 of
`ramanujan_spike_pointwise/definition.md`). The present refinement:

1. **Promotes the corollary to a closed-form identity on the full
   domain** (including non-coprime `n`, where the value is zero).
2. **Extends to all squarefree `W`** (not just primorials).
3. **Reduces the non-squarefree case to its radical** with explicit
   `phi(rad(W))` factor.
4. **Exhibits the L² closed form** `||T_W^{div} − chi_P||²/N` via
   explicit prime-distribution count (rel_err < 10⁻⁹).
5. **Quantifies the Pearson gap** between full and divisor-only,
   localising the prime-discriminating signal to the non-divisor
   squarefree conductors `q ≤ W` (the squarefree integers `≤ W`
   that fail to divide `W`).

The cleanest single-sentence consequence: the divisor-only restriction
of `T_W` is **just the wheel sieve** with prime-density calibration;
all of the prime-discriminating content of the full `T_W` lives in the
non-divisor squarefree conductors.

## Refinement of E2.1

The MPS bond-dim identity's `phi(W)/W` factor (E2.1) was previously a
*rank quotient*. After this session it is also the value of an
**exact pointwise scalar field**, namely

```
    T_W^{div}(n)  =  (π(N)/N) · W/φ(W) · [gcd(n, W) = 1].
```

That is: the Mertens-product compression ratio of E2.1 is realised
pointwise as the value of a scaled wheel indicator. EDGES.md E2.1
annotated inline with this refinement.

## What did this session produce that was not in the project before?

1. **Closed-form theorem on the full domain** (S208, prove + verify).
2. **L² closed form** for `||T_W^{div} − chi_P||²/N` via explicit
   prime-distribution count (verified rel_err < 10⁻⁹).
3. **Quantitative Pearson gap table** across `W = 2..2310` showing
   monotonic growth from 0 to 0.281, localising 100% of T_W's
   prime-discriminating content beyond wheel admissibility to
   non-divisor squarefree conductors.
4. **Radical reduction** for non-squarefree `W` (depends only on
   `rad(W)`).

Edges composed: E2.1 (refined inline), E2.2 (recovered as `W=2`
specialisation: `M_2^{div}(n) = 2 · [n odd]`), C9 / S191 (parent
construction; this is the divisor-only restriction).

## Honest self-grade

**B-grade (substantive refinement).** Single-session closed-form
identity with exact pointwise verification across 10 W × 2 N cells.
The mathematical content is a project-internal Möbius collapse — no
cross-domain import required, just careful algebra. The novelty bar
(per CLAUDE.md): "a published paper-grade NT/complexity theorist
could not derive this in an afternoon from CLOSED_PATHS + EDGES
alone". The current result is not at that bar — given the S191
construction, this divisor-only specialisation is straightforwardly
derivable in 30 minutes by a competent NT student. So **B, not A**:
substantive refinement that extends S191's coprime-only corollary to
the full domain and quantifies the divisor-only-vs-full gap, but no
new cross-domain technique and no algorithmic opening.

## Self-evaluation (CLAUDE.md 4 questions)

1. **What did I produce that was not in the project before this
   session?** A closed-form pointwise identity for `M_W^{div}(n)` on
   the full domain (extending S191's coprime-only corollary), an
   explicit L² closed form via prime-distribution count, and a
   quantitative Pearson decomposition isolating the prime-
   discriminating signal beyond wheel admissibility to the non-divisor
   squarefree conductors. Verified at machine-epsilon to N=10⁶.
2. **What edges did my work compose or cite?** E2.1 (refined inline),
   E2.2 (recovered as W=2 special case), C9 / S191 (parent;
   divisor-only is the restriction).
3. **If my session produced only duplicate closures, why?** It did
   not; the closed form on the full domain is new content that
   sharpens S191's coprime-only statement. The algorithmic content,
   however, collapses to the wheel sieve — so this is "structural
   refinement plus duplicate-of-wheel-sieve" rather than a polylog
   opening. Honest reading: a B-grade refinement, not novelty.
4. **What is the next-action for the next agent?** Three successor
   challenges proposed in `NOVELTY_CHALLENGES.md` and
   `spike_divisor_only_results.md`: (i) decompose the 0.281 Pearson
   gap by conductor, (ii) Lean-formalise the divisor-only identity
   (≈50-line Lean target via L6 + Möbius collapse), (iii) Λ-modulated
   divisor-only sum.

## Artefacts

* `experiments/constructions/spike_divisor_only_wheel/definition.md`
* `experiments/constructions/spike_divisor_only_wheel/spike_divisor_only.py`
* `experiments/constructions/spike_divisor_only_wheel/spike_divisor_only_results.md`
* `experiments/constructions/spike_divisor_only_wheel/spike_divisor_only_results.json` (N=10⁵)
* `experiments/constructions/spike_divisor_only_wheel/spike_divisor_only_results_N1e6.json`
* `EDGES.md` E2.1 inline annotated with S208 refinement
* `NOVELTY_CHALLENGES.md` C9.a marked BUILT with successor challenges
* `status/CLOSED_PATHS.md` row appended
