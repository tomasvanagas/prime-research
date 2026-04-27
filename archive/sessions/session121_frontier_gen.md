# Session 121 — Frontier generation: 5 new ATTACK_VECTORS entries (D20–D24)

**Mode:** frontier_gen (auto-fired).
**Date:** 2026-04-27.

## Why frontier_gen fired

The harness fires `frontier_gen` when ATTACK_VECTORS open count drops below
4, OR no A-grade in last 20 sessions, OR 2 consecutive F-grade sessions.
Open vectors are not low (30+ remain), so this is most likely an
A-grade-scarcity fire. Either way, the job here is to extend the
frontier by importing UNUSED cross-domain techniques.

## Method

1. Read ATTACK_VECTORS.md (3784 lines) and CROSS_DOMAIN_TECHNIQUES.md
   (220 lines) to enumerate UNUSED rows of the cross-domain registry.
2. Cross-checked candidates against status/CLOSED_PATHS.md (773 lines)
   to avoid silently duplicating closures: critical hits at lines 356
   (`Cay(Z/xZ, primes)`), 387 (GCD spectrum = Ramanujan sums), 423
   (transfer operator = ζ-zeros), 425 / 654 (tropical), 671 (étale),
   598 / 680 (adelic). All ruled out matching candidates.
3. WebFetched survey/foundational refs for each picked technique. The
   restriction-theory candidate was dropped after 4 URL attempts failed
   to surface the Tao 2003 / Stein-Tomas survey directly; replaced
   with Reeb-graph (Edelsbrunner-Harer textbook ref, no URL needed).
4. Picked 5 techniques with concrete single-session protocols and
   pre-stated falsification criteria.
5. Appended D20–D24 to ATTACK_VECTORS.md under §D before §E. Updated
   CROSS_DOMAIN_TECHNIQUES.md status (UNUSED → PROPOSED) and Priority
   hints section.

## The 5 new vectors

### D20 — Friedman spectral gap of `Cay(Z/NZ, primes < N^c)`

Test if `λ_2(G_N) / 2√(d-1) → 1` (Ramanujan-typical), `< 1` (super-
Ramanujan, arithmetic gain), or `> 1` (sub-Ramanujan). Computable as a
prime exponential sum via FFT in `O(N log N)`. Distinct from CLOSED line
356 (which closed `λ_0 = π(N^c)`, the *trivial* eigenvalue) and from
CLOSED §A.A3 / E7.12 (which closed *fixed-generator* `Cay((Z/nZ)*,
{2,3,5})` as a *pointwise primality test* probing `ω(n)`). The
quantitative `λ_2` saturation question is OPEN in the literature.

### D21 — Reeb graph of arithmetic height functions on Z

Exact Reeb graph of `f : Z_{≤ N} → R` for arithmetic `f` (e.g.
`ω(n)`, the number of distinct prime factors). No metric (vs PH on
Takens embeddings, CLOSED S96), no lens (vs Mapper, D19 proposed), no
poset (vs discrete Morse, D17 proposed). Tests whether the loop count
`β_1(R_N)` of the Reeb graph deviates from an Erdős-Kac-Gaussian
matched control by `> 5σ`.

### D22 — Higher-order Hodge Laplacian on coprimality FLAG complex

`L_k = ∂_{k+1} ∂_{k+1}^* + ∂_k^* ∂_k` for `k ≥ 1` on the simplicial
flag complex of pairwise coprime `(k+1)`-tuples in `[2, N]`. CLOSED
lines 356, 387 closed only `L_0` (graph Laplacian = Ramanujan sums =
Meissel-Lehmer cost). Higher-order `L_k` encodes triple-coprime /
Möbius cohomological data orthogonal to bilinear forms. Tests
spectral gap and Betti-spectrum deviation from Erdős-Rényi flag complex.

### D23 — Density Hales-Jewett combinatorial-line density of primes

For `n ∈ [0, k^d)` viewed as base-`k` digit string in `[k]^d`, measure
the density of primes ON COMBINATORIAL LINES (parametrised by wildcard
patterns `w ∈ {0, ..., k-1, *}^d`). Combinatorial lines are STRUCTURED
APs of common difference `Σ_{i ∈ S} k^i` over `S ⊂ [d]`. Compares
empirical line-density to the Hardy-Littlewood prediction `1/log(n_0)
⋅ S(d)` for the specific singular series. Distinct from D6 (Gowers
`U^k`, CLOSED S85, E2.13) which measured the unrestricted-AP norm; D23
measures structured-AP density.

### D24 — Eynard-Orantin topological recursion of prime correlation moments

Define `W^{(n)}(z_1, ..., z_n; N) = Σ_{p_1 ≠ ... ≠ p_n ≤ N} ∏ 1/(z_i -
p_i)` and test if the hierarchy satisfies the Eynard-Orantin
topological recursion with an arithmetic spectral curve `y^2 = R(x)`.
A-grade success: prime correlations are polylog-determined from `R(x)`
via residue calculus. Failure: explicit obstruction at the loop-equation
level. Distinct from CLOSED §C.C2 (Conrey-Snaith) which gave
fixed-`n` arithmetic corrections; D24 asks for the FULL hierarchical
recursion.

## Cross-domain literature consulted

| Vector | Survey ref | Status of fetch |
|--------|-----------|-----------------|
| D20 | Friedman 2008 Memoirs AMS = arXiv:cs/0405020 | Confirmed: λ_2 ≤ 2√(d-1) + ε w.h.p. for random d-regular |
| D21 | Edelsbrunner-Harer 2010 *Computational Topology* AMS, ch. VI.3 | Reference; no URL needed (textbook) |
| D22 | Lim 2020 SIAM Review = arXiv:1507.05379 | Confirmed: elementary intro to combinatorial Hodge Laplacian on simplicial complex via boundary maps |
| D23 | Polymath1 2010 *Annals* 175 = arXiv:0910.3926 | Confirmed: density-δ subset of `[3]^n` contains combinatorial line for `n ≥ tower(O(1/δ³))` |
| D24 | Eynard-Orantin 2007 *Comm. Number Theory Phys.* 1 = arXiv:math-ph/0702045 | Confirmed: topological recursion / loop-equation framework with spectral-curve input; matches matrix-model genus expansion |

## Self-grade

**B.** All five vectors are grounded in genuinely UNUSED cross-domain
techniques (Friedman, Reeb, higher-order Hodge, DHJ, Eynard-Orantin).
Each has a concrete single-session first step and pre-stated
falsification criteria. Each is structurally distinguished from the
nearest existing closure, with the distinction documented in the
"Distinction from..." paragraph.

I withhold A because:
- Fetching surveys hit several arXiv-ID misfires (got wrong papers
  for cs/0405020, math/0306172, 1507.05379, 1010.1955, math/0311025,
  math/0301070); only ~50% of fetches surfaced the intended content.
  The correct papers were eventually accessible (cs/0405020, 1507.05379,
  math-ph/0702045) and I cited them faithfully, but the fetch
  latency suggests I have less hands-on familiarity with each
  technique than ideal for evaluating A-grade potential.
- Two of five (D21 Reeb, D22 Hodge) are in the topological /
  spectral category and may produce only B-grade outcomes — the
  HL-saturation pattern of recent closures (E2.13 / E2.14 / E2.16 /
  E2.17 / E2.18 / E2.19) suggests new "topological pseudorandomness
  measures" land at the noise floor by default.
- D23 (DHJ) overlaps in spirit with Gowers `U^k` (E2.13); the
  distinction (structured-AP vs unrestricted-AP) is real but might
  reduce to a Hardy-Littlewood saturation in the same family.

The A-grade candidates among the 5 are most plausibly D20 (Friedman
super-Ramanujan would be a genuine arithmetic gain over Vinogradov)
and D24 (matrix-model topological recursion of prime correlations
would, if it holds, falsify the information-theoretic floor). Both
are speculative.

## Next-action

Per CLAUDE.md autonomy invariants: future BUILD or CLOSE sessions
on D20–D24 should propose 1-2 follow-on challenges in
NOVELTY_CHALLENGES.md upon completion. The §1 spectral category now
has 2 fresh vectors (D20, D22); the §4 topological category has 1
fresh (D21); the §7 additive category has 1 fresh (D23); the §3
random-matrix category has 1 fresh (D24). The frontier should now be
adequate for at least the next 5 production sessions.

## Files modified

- `ATTACK_VECTORS.md` — appended D20–D24 (~5 entries × ~120 lines each)
- `CROSS_DOMAIN_TECHNIQUES.md` — 5 rows updated UNUSED → PROPOSED with
  full survey references; Priority hints section extended with 5 new
  bullets.
- `archive/sessions/session121_frontier_gen.md` — this file.
