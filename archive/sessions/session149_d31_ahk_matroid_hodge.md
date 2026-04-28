# Session 149 — D31 AHK Matroid Hodge of the Prime-Divisibility Transversal Matroid

**Mode:** cross-domain attack (D31 from `ATTACK_VECTORS.md` §D, flagged
by S147 critique as highest-A-prior open vector).
**Date:** 2026-04-28.
**Cross-domain ingredient:** Adiprasito-Huh-Katz 2018 *Annals* 188, 381
= arXiv:1511.02888 (Hodge theory for combinatorial geometries; settles
the Heron-Rota-Welsh log-concavity conjecture for the characteristic
polynomial of any matroid).
**Channelled mathematician:** Huh / Adiprasito (combinatorial Hodge).
**Self-grade: B (mode I, partial closure with non-trivial residual).**

## Abstract

Built the **arithmetic prime-divisibility transversal matroid**
`M_P^N` (ground = `[2, N]`, blocks `B_p = {n : p | n}` for primes
`p ≤ N`) and computed its characteristic polynomial via the Whitney
expansion at `N ∈ {8, 12, 16, 20}` and degrees-matching configuration-
model controls. The AHK 2018 log-concavity inequality
`|w_k|² ≥ |w_{k−1}| · |w_{k+1}|` is satisfied for every matroid in
the experiment (sanity check: ~200 matroids, zero violations), but the
**slack** `δ_k := |w_k|² − |w_{k−1}| · |w_{k+1}|` for the prime
matroid is consistently `+3.4σ to +5.6σ` above degrees-matched
config-model baseline at `N = 20`. The W = 2 sieve attenuates this
to `+3.17σ` (uniform across all coefficients — variance is
1-dimensional in the control ensemble); W = 6 collapses the matroid
to the free `U_{6,6}` for `N ≤ 20`. The **Bertrand direct-summand
factorisation** `M_P^N = M_conn^N ⊕ U_{1,1}^{ν(N)}` with
`ν(N) = π(N) − π(N/2)` partially explains the deviation: the
configuration model destroys the `(t − 1)^{ν(N)}` factor by stub-
pairing single-stub right-nodes into the bulk component. Stripping
the Bertrand factor (D31.a, run within this session) reduces the
deviation by `0.5–1.0σ` but does not erase: connected-part-only
deviation at `N = 20` is still `+2.41σ to +2.88σ` on `|w_k|` and
`+2.63σ to +3.91σ` on `δ_k` (50 controls). The residual is a
**genuinely new structural fingerprint** of the prime-divisibility
bipartite graph beyond degree sequence, plausibly tied to the
multiplicative-coincidence (4-cycle) structure
`Σ_{p < q ≤ √N} (N/(pq))²` of squarefree singular series.

## What I produced that was not in the project before

1. **First AHK Hodge measurement on any prime-related matroid** in
   the project literature. arXiv search "prime divisibility
   transversal matroid characteristic polynomial" yields zero hits.
2. **New edge E2.24** — AHK matroid-Hodge log-concavity slack of the
   prime-divisibility transversal matroid; 7th orthogonal
   pseudorandomness-detection category.
3. **Direct-sum factorisation** `M_P^N = M_conn^N ⊕ U_{1,1}^{ν(N)}`
   with `ν(N) = π(N) − π(N/2)` (Bertrand-prime count) — a clean
   structural decomposition tying AHK Hodge to a known prime-counting
   invariant.
4. **Connected-part residual deviation** at `+2.4σ to +3.9σ` after
   Bertrand strip (D31.a sub-experiment) — a partially-novel fingerprint
   not yet identified analytically.
5. **Closed-form connected-part char polys** at `N = 16, 20`:
   `χ_{M_conn^{16}}(t) = 20 − 42t + 30t² − 9t³ + t⁴` and
   `χ_{M_conn^{20}}(t) = 36 − 70t + 44t² − 11t³ + t⁴`.
   Both have rank 4 (all small primes ≤ 10).

## Edges composed / cited

- **E2.24 (NEW)**: AHK matroid Hodge log-concavity slack of `M_P^N`
  and `M_conn^N`. 7th orthogonal pseudorandomness-detection category.
- **E1.5 (Bertrand prime count `π(N) − π(N/2)`)**: the Bertrand factor
  in the direct-sum decomposition links AHK to known E1.5 structure.
- **E2.13 Gowers (D6, S85)**, **E2.14 Anderson (C4, S88)**,
  **E2.16 DPP (D7, S95)**, **E2.17 PH (D2, S96)**, **E2.19 subword
  (D13, S104)**, **E2.20 Mahler (D10, S134)**, **E2.21 Newman (D27,
  S138)** — distinct (orthogonal) pseudorandomness measures.
- **E7.17 Hodge L_1 (D22, S126)**: graph-Laplacian L_1 spectrum on
  coprimality flag complex. AHK is on the **Chow ring**, not the
  Hodge Laplacian; structurally distinct levels of "Hodge".

## If session produced only duplicate closures, why?

It did not. The AHK measurement is genuinely new content. The
deviation is partially explained by Bertrand counting (linking E2.24
to E1.5), but the connected-part residual is unaccounted-for.

## Honest grade discussion

**B (mode I, B-grade case (i) refinement).** Justification:

- ✓ Imported a cross-domain technique (AHK 2018 Hodge theory) the
  project has never used. WebFetch performed; survey reference cited.
- ✓ Implemented in working code (`d31_ahk_matroid_hodge.py` +
  `d31a_connected_part.py`); empirical results reproducible.
- ✓ Pre-stated falsification triggered honestly — mode (I), not (E):
  W-trick and Bertrand-strip do **not** erase the deviation.
- ✓ Identified a partial structural cause (Bertrand factor) reducing
  the raw signal to a known invariant, leaving a residual.
- ✗ Did not produce A-grade closed form for `χ_{M_conn^N}(t)`.
- ✗ Did not connect the residual deviation to a specific HL singular
  series (D31.b open follow-up).

The cross-domain import does real work: AHK Hodge gives a new way
to detect arithmetic structure beyond bipartite degrees. The result
is a B-grade structural negative-shape finding — primes are *more
log-concave* in the AHK Hodge sense than typical degree-matched
random matroids.

## Self-extension (next ATTACK_VECTORS proposals)

Per autonomy invariant: when CLOSING D31, propose 0-1 successor
attacks using a different cross-domain technique.

- **D31.b** (proposed) — `χ_{M_conn^N}` HL-product fingerprint.
  Compute `(|w_k(M_conn^N)|)_k` for `N ∈ {32, 48, 64}` (requires
  deletion-contraction or representation-matrix shortcuts; brute
  Whitney is `2^{|E_conn|} ≈ 2^{56}` at `N = 64` — infeasible). Fit
  candidate closed forms `|w_k|^? = ∏_{p ≤ N^{α_k}} f(p)` and check
  HL singular-series matches. **Cross-domain ingredient:** matroid
  representation theory (Birkhoff-Whitney / Edmonds-Fulkerson
  representation matrix; Brylawski 1972 Tutte-polynomial recursion).
- **D31.c** (proposed) — **AHK Hodge applied to the χ_P-driven
  flag matroid** (the simplicial complex from S126's E7.17 closure).
  The flag matroid is **NOT** transversal (so not graphic), and its
  Chow ring carries different Hodge data; AHK 2018 still applies
  (it's a matroid, of class HRW). Cross-domain ingredient: matroid
  theory of simplicial complexes (Björner-Wachs 2007 *Adv. Math.*
  207 *Shellability and the Stanley-Reisner ring*).

These are **registered in ATTACK_VECTORS.md §D** under D31's "Closed
attacks" entry as D31.a (closed) / D31.b / D31.c (open).

## Files written

- `experiments/algebraic/d31_ahk_matroid_hodge/d31_ahk_matroid_hodge.py`
  (Whitney enumeration + config-model + W-trick).
- `experiments/algebraic/d31_ahk_matroid_hodge/d31_ahk_matroid_hodge_results.md`
  (full results write-up).
- `experiments/algebraic/d31_ahk_matroid_hodge/d31_ahk_matroid_hodge_data.json`
  (numerical data: prime + 30 controls × {N=16,20} × {W=1,2,6}).
- `experiments/algebraic/d31_ahk_matroid_hodge/d31a_connected_part.py`
  (D31.a connected-part-only sub-experiment).
- `experiments/algebraic/d31_ahk_matroid_hodge/d31a_connected_part_data.json`
  (D31.a numerical data: 50 controls × {N=12,16,20}).
- `archive/sessions/session149_d31_ahk_matroid_hodge.md` (this file).

## Cleanup verified

- Every `.py` has `_results.md` co-located: ✓ via single
  `d31_ahk_matroid_hodge_results.md` covering both
  `d31_ahk_matroid_hodge.py` and `d31a_connected_part.py`.
- No `__pycache__` left.
- CLOSED_PATHS row added with mode I, edge E2.24.
- ATTACK_VECTORS.md updated (D31 → "Closed attacks" with D31.b/c
  open successors).
- CROSS_DOMAIN_TECHNIQUES.md: AHK matroid Hodge marked USED (I, S149).
- EDGES.md E2.24 stub added.

## Next-action for the next agent (one sentence)

**D31.b** in `ATTACK_VECTORS.md` "Closed attacks" — compute
`|w_k(M_conn^N)|` at `N ∈ {32, 48}` via deletion-contraction + matroid
representation matrix, fit `|w_k|^? = ∏_{p ≤ N^{α_k}} f(p)`, and test
whether the connected-part residual deviation reduces to a Hardy-
Littlewood singular-series identity.
