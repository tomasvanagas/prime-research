# Session 92 — Algebraic Immunity / Polynomial-Method Invariants of chi_P

**Mode:** wild swing — §B1 of ATTACK_VECTORS.md
**Date:** 2026-04-26
**Self-grade:** **B-grade** (ambitious frontier attack, structural
informative failure; adds new EDGE; closes §B1)
**Channel:** Tao (slice rank, additive combinatorics) + Croot-Lev-Pach
(polynomial method) + Courtois-Meier (algebraic immunity from
cryptanalysis).

## What I tried

§B1 of ATTACK_VECTORS.md asks: does the cap-set polynomial method
(Croot-Lev-Pach / slice rank / Tao) yield a polylog primality witness
via low-degree annihilators of chi_P? Three structurally distinct
invariants attacked in one session:

1. **Algebraic immunity AI(chi_P) over F_2** for N=4..13 — direct LP
   over F_2 (Gauss elimination on the `2^N x sum_i C(N,i)`
   monomial-evaluation matrix). Plus annihilator extraction via
   nullspace recovery.
2. **F_p multilinear ANF degree** of chi_P viewed as `(F_p)^k -> F_p`
   via base-p expansion, for `(p, k)` in
   `{(3,2..5), (5,2..3), (7,2..3), (11,2)}`.
3. **Slice-rank brackets** via Tao 2016 inequality
   `slice_rank(T) >= max_axis rank(flatten_axis(T))` (LB) plus
   greedy heaviest-slice peeling (UB).

Plus a W-trick correction test: AI of `chi_P_{W,b}(n) := chi_P(W*n + b)`
for `W in {1, 2, 6, 30}`, b coprime, to determine whether the AI
deviation persists after sieving.

## What I produced

`experiments/algebraic/algebraic_immunity_chi_p/`:
- `algebraic_immunity_chi_p.py` — three-invariant driver
- `extract_annihilator.py` — annihilator extraction with verification
- `wtrick_AI.py` — W-trick correction test
- `algebraic_immunity_chi_p_results.md` — full writeup with falsifiers
- JSON data files for AI, F_p ANF, slice rank, and W-trick

Plus updates to: CLOSED_PATHS.md (S92 row), EDGES.md (E2.15),
ATTACK_VECTORS.md (B1 closed + closure entry with self-extension).

## Main empirical finding

**`AI(chi_P, N) = 2` for ALL N in [4, 13]**, while
`AI(random matched-density, N) = Theta(log_2(1/rho))` (reaches 4 at
N=11/12 with zero std). AI deviation `AI_chi_P - AI_random` reaches
`-2` by N=11.

Liouville+ AI grows like AI_random; Mobius!=0 also has constant AI=2.

## Mechanism (the structural identification)

Annihilator extraction returns the **identical polynomial g for every
N from 5 through 13**:

```
g(x) = 1 + x_0 + x_1 + x_{0,1} = (1 + x_0)(1 + x_1)
```

`g(n) = 1` iff `n ≡ 0 mod 4`. The annihilation works because no prime
> 2 is divisible by 4, and the only prime with bit_0 = 0 is n = 2,
which is annihilated by `x_1(2) = 1`. Mobius!=0 inherits the same
annihilator (n ≡ 0 mod 4 implies non-squarefree). Liouville+ does NOT
inherit (`lambda(4) = +1`).

**AI(chi_P) = 2 is the polynomial-method encoding of the trivial
mod-4 sieve fact.**

## W-trick correction — confirms the picture

`chi_P_{W,b}(n) := chi_P(W*n + b)` for gcd(b, W) = 1:

| W  | b      | N=8 | N=11 | random_mean N=11 |
|----|--------|-----|------|------------------|
| 1  | 0      | 2   | 2    | 4.00 (std 0)     |
| 2  | 1      | 3   | 4    | 4.50             |
| 6  | 1      | 4   | 5    | 5.00 (std 0)     |
| 6  | 5      | 4   | 5    | 5.00 (std 0)     |
| 30 | 1      | 4   | 5    | 5.00 (std 0)     |
| 30 | 7      | 4   | 5    | 5.00 (std 0)     |
| 30 | 11     | 4   | 5    | 5.00 (std 0)     |

**The AI deviation is fully removed by W >= 6**, with chi_P_{W,b}
tracking random matched-density within ±1 step at all tested N
(exact match with zero std at most cells).

## F_p multilinear ANF degree

chi_P near-saturates max possible degree `(p-1)*k` at all tested
`(p, k)` other than (3, 3) and (5, 3) where chi_P drops by 1 below
the random max — within "almost all coefficients populated" noise.

## Slice-rank brackets

p=2 non-informative (2-row flattenings cap rank at 2 trivially, both
chi_P and random); p=3 chi_P matches random for k>=3.

## What edges I composed or cited

- **Cited**: E2.13 (Gowers U^k of chi_P -> HL singular series, S85),
  E2.14 (Anderson Lyapunov W-cascade, S88), E1.10 (35-measure
  pseudorandomness), E3.13 (BK detection limits).
- **Added**: E2.15 — Algebraic immunity AI_F_2(chi_P) = 2 with
  explicit annihilator (1+x_0)(1+x_1), removed by W-trick at W >= 6.
- **Composed**: triple of independent confirmations of "chi_P
  structure = HL equidistribution mod q" in three orthogonal
  mathematical categories:
  1. additive combinatorics (Gowers, E2.13)
  2. spectral / transfer-matrix Lyapunov (Anderson, E2.14)
  3. Boolean polynomial method / algebraic cryptanalysis
     (algebraic immunity, E2.15)

## Why this is B-grade and not A-grade

The session attempted a **frontier target** (§B1, listed at high
A-grade probability in the prompt's priority order) and produced a
**structural informative failure** with a new EDGE.

A-grade was reserved for "polynomial-method shortcut bypassing the
sqrt(x) wall" or "AI deviation persistent past the W-trick" — the
W-trick fully explained the deviation, so A-grade was not reached.

B-grade because:
- the failure is STRUCTURAL: AI(chi_P) = 2 mechanism is fully
  explained by an explicit polynomial encoding the mod-4 sieve fact,
  with W-trick verification;
- adds a NEW edge (E2.15) with closed-form annihilator;
- provides a THIRD independent confirmation of the W-trick picture in
  a category orthogonal to E2.13 and E2.14;
- session produced 4 concrete deliverables (script, results.md,
  extracted annihilator, W-trick verdict).

## What was novel beyond S17 (tensor rank closure)

S17 closed tensor rank of chi_P (full / matches random) using a
3-way bit-partition approach. **This session's invariants are
structurally distinct:**

- **Algebraic immunity** is the smallest-DEGREE polynomial annihilating
  chi_P; tensor rank is the smallest-RANK CP decomposition. They
  measure different things (degree vs rank in different bases).
- **F_p multilinear ANF degree** is the polynomial degree over F_p
  via base-p expansion — not addressed by S17's binary tensor work.
- **Slice rank brackets via Tao's flattening lower bound** — S17's
  greedy slice rank gave only upper bounds.

The mod-4 annihilator (1+x_0)(1+x_1) was explicitly extracted with
mathematical interpretation, not done in any prior session.

## Self-extension follow-ups (per CLAUDE.md)

Filed in ATTACK_VECTORS.md "Closed attacks" entry for §B.B1:

- **(a)** AI of multiplicative-regime functions (Liouville centered,
  Mobius), per §G of ATTACK_VECTORS — fourth-leg confirmation in the
  multiplicative regime.
- **(b)** Higher-degree structured annihilators of chi_P_{W,b}:
  beyond the mod-4 fact, are there unexpectedly-low-degree degree-O(1)
  annihilators of W-tricked chi_P? Concrete sub-falsification:
  AI(chi_P_{6,1}) < AI(random_{6,1}) by >= 1 at N >= 10.

## Self-evaluation (CLAUDE.md questions)

1. **What did I produce that was not in the project before?**
   - The **explicit closed-form algebraic immunity** AI_F_2(chi_P) = 2,
     with the specific annihilator polynomial (1+x_0)(1+x_1) extracted
     and verified across N=5..13.
   - The W-trick correction at AI level, confirming W >= 6 fully removes
     the deviation (zero std at multiple cells at scale N=11).
   - F_p multilinear ANF degree of chi_P at 9 (p, k) cases — never
     computed in this project before.
   - **EDGE E2.15** as a third leg of the W-trick triple.

2. **Edges composed/cited:** E2.13 + E2.14 (parallel mod-q
   confirmations); E1.10 / E3.13 (35-measure pseudorandomness).

3. **Did I produce only duplicate closures?** No — produced a new
   edge with explicit closed-form annihilator + cross-domain
   confirmation. The closure mode E was structural, not "duplicate
   of existing CLOSED_PATHS row." The S17 tensor rank closure is
   structurally distinct (rank vs degree).

4. **Next-action for next agent:** §G2 (Liouville Gowers U^k) or
   §G1 (Liouville Anderson Lyapunov) — the multiplicative-regime
   parallels of E2.13/E2.14, currently §G of ATTACK_VECTORS. The
   AI follow-up (a) is even simpler — single script change + one
   N range run, ~ 1 hour of work.

## Honest failure clause

This session DID produce a novel artefact (E2.15 + the explicit
annihilator). This is NOT a session that "produced only DUPLICATE
closures of fresh-perspective brainstorms with no structural reason
added" — the structural reason (mod-4 sieve fact -> (1+x_0)(1+x_1)
annihilator) IS the structural content.

The session is graded B (not C) because:
- The attack vector (§B1) was a frontier target, not a refinement.
- The deliverable is a new edge with closed-form structure.
- The triple (E2.13 + E2.14 + E2.15) constitutes a meaningful
  scientific picture: chi_P deviates from random in three orthogonal
  categories, all via the same mod-q sieve mechanism.

The session is NOT graded A because:
- The deviation was fully captured by the W-trick — no NEW
  arithmetic content beyond mod-q.
- Polynomial method gave no compression / polylog opening.

This grade should NOT be inflated by future verify-mode sessions.
