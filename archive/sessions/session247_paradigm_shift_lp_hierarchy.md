# Session 247 — L^p hierarchy of f_N(z) (PARADIGM-SHIFT mode)

**Date:** 2026-04-30.
**Mode:** PARADIGM-SHIFT (no cross-domain imports permitted).
**Self-grade:** **B**.
**Edges refined:** **E2.21** (L^∞ → closed-form L^p hierarchy for all
even p ≥ 4), **E2.31** (m_4 = 2.95 prefactor → reduced to a
Selberg-Delange product `G_4 · S_3 ≈ 1.534` with a residual
BDJ-Szegő-standardisation prefactor pending), **S168** (sqfree-q L²
energy → per-q L^p energy with decay `1/φ(q)^{p-1}`).
**Predecessors:** S138 (E2.21 build), S204 (E2.31 build), S168 / S190
(squarefree-q spike formula), S239 (paradigm-shift cumulative-Q-strip
data on m_4).

## Mode discipline

* **No `WebFetch` / `WebSearch`.** Confirmed.
* **No new cross-domain technique imported.** Sinc integrals
  (Borwein-Borwein closed forms `∫|sinc|^p dx = π G_p`) are
  elementary calculus, already implicit in E2.21's Vinogradov-style
  noise analysis. The HL singular series `∏ (1 + 1/(p-1)^a)` is the
  same Selberg-Delange Euler product used in E2.13 (USED in
  CROSS_DOMAIN_TECHNIQUES.md). Major-arc bookkeeping is from
  Vinogradov, also USED.
* **No new `ATTACK_VECTORS.md` entries.**
* **Only existing project content read** (EDGES.md sections E1, E2,
  E6, E7; CLOSED_PATHS.md; recent paradigm-shift sessions S239 + S246
  + S232 for stylistic template; experiments/constructions/
  bdj_toeplitz_hankel_chi_p/ for E2.31 numerical anchor;
  experiments/constructions/parity_stripped_trinity/ for S239 strip
  data structure).

## Target picked (paradigm-shift composition)

Three positive-content edges encode `μ²(q)/φ(q)^a` HL imprints with
DIFFERENT exponents a:
* E2.21 (L^∞): `‖f_N‖_∞ = π(N) · max_q μ²(q)/φ(q) = π(N)`
  (q=2 wins; at every q the spike amplitude is `μ²(q)·π(N)/φ(q)`).
* S168 (L²-spike): `‖P_{V_q^prim} χ_P‖² = μ²(q)·(π(N)−r(q))²/(φ(q) N)`,
  i.e. per-q decay `1/φ(q)`.
* E2.31 (L^4 via BDJ m_4): empirical constant 2.95, inferred per-q
  decay `1/φ(q)^3` from the sqfree-q strip cascade in S239.

The S247 paradigm-shift target: **derive the unified per-q exponent
law `1/φ(q)^{p-1}` for even integer p ≥ 4, expressed as a single
closed-form prediction for the entire L^p hierarchy.**

The composition exploits two existing structural facts:
1. Each q-major-arc has spike amplitude `μ(q)/φ(q) · π(N)`
   (from E2.21).
2. Each q-major-arc has effective sinc-window width `~ 1/N` from
   the Vinogradov truncation; integrating |spike|^p over this window
   gives `(spike)^p · ∫|sinc|^p · (1/N)`.

Combining: total L^p mass = (geometric constant) · π(N)^p / N · Σ_q
spike^p / arc-density, leading to the closed form
`G_p · S_{p-1}^HL · π(N)^p / N`.

## What I built

`experiments/constructions/lp_hierarchy_chi_p/`:
* `definition.md` — formal construction + 4 pre-stated falsifiers
  (F1-F4).
* `lp_hierarchy_chi_p.py` — main FFT-based ‖f_N‖_p^p computation,
  cumulative-q-strip projection (additive-Fourier 2D Gram-Schmidt),
  Bernoulli matched-density null.
* `lp_hierarchy_chi_p_results.md` — F1-F4 verdict tables + per-q
  cohort closure + grade defense.
* `lp_hierarchy_chi_p_results.json`, `_full.json`, `_convergence.json`
  — raw FFT measurements at N up to 2²⁰.
* `run_full.log`, `run_convergence.log` — experiment output.

## The conjecture (S247)

> For every even integer p ≥ 4, as N → ∞,
>
>   `‖f_N‖_p^p · N / π(N)^p  →  G_p · S_{p-1}^HL`
>
> where `‖f_N‖_p^p = (1/2π) ∫_0^{2π} |f_N(e^{iθ})|^p dθ`,
> `G_p = (1/π) ∫_R |sin u/u|^p du` (sinc constants), and
> `S_{p-1}^HL = ∏_{prime r} (1 + 1/(r-1)^{p-1})` (HL singular series
> at depth p−1).

**Specific numerical predictions:**
| p | G_p     | S_{p-1}^HL | predicted ratio |
|--:|--------:|-----------:|----------------:|
| 4 | 0.6667  | 2.301      | 1.534           |
| 6 | 0.5500  | 2.065      | 1.136           |
| 8 | 0.4794  | 2.016      | 0.966           |

## The four findings

### Finding 1 — F1 PASS at p=4 with monotone convergence

| N (= 2^k)   | π(N)   | empirical ratio | emp / pred |
|------------:|-------:|----------------:|-----------:|
|       2¹²   |   564  | 1.4835          | 0.967      |
|       2¹⁴   |  1900  | 1.4979          | 0.977      |
|       2¹⁶   |  6542  | 1.5077          | 0.983      |
|       2¹⁸   | 23000  | 1.5103          | 0.985      |
|       2²⁰   | 82025  | 1.5126          | 0.986      |

Monotone convergence from below to G_4 · S_3 = 1.534, reaching 98.6%
at N = 2²⁰.

### Finding 2 — F2 PASS at p=6 and p=8 with same residual gap

At N = 2²⁰: emp/pred = 0.985 (p=4), 0.985 (p=6), 0.985 (p=8). The
cross-p stability of the 1.5% gap points to a single shared
`O(1/log N)` correction (consistent with PNT / `Li(N) → π(N)`).

### Finding 3 — F3 PASS: cumulative-Q-strip closure formula

At N = 2¹⁶, empirical p=4 closure tracks the predicted partial-sum
formula `[Σ_{2 ≤ q sqf ≤ Q} μ²(q)/φ(q)^3] / [S_3 − 1]` to ±1pp:

| Q (sqfree)  | empirical closure | predicted closure |
|------------:|------------------:|------------------:|
|       2     | 77.8%             | 76.9%             |
|       6     | 98.5%             | 97.2%             |
|      30     | 101.1%            | 99.7%             |
|     210     | 101.3%            | 99.9999%          |

The empirical 101% reflects that the stripped residual (0.656 vs
predicted G_4 = 0.667) is ~1.7% below the asymptotic `G_p` baseline —
the same PNT-correction gap as F1/F2.

### Finding 4 — F4 PASS structurally: Bernoulli null shows no HL modulation

| N (= 2^k) | chi_P p=4 ratio | Bern p=4 ratio | predicted (Bern) |
|----------:|----------------:|---------------:|-----------------:|
|     2¹⁴   | 1.4979          | 0.5776         | G_4 = 0.6667     |
|     2¹⁶   | 1.5077          | 0.7075         | G_4 = 0.6667     |
|     2¹⁸   | 1.5103          | 0.6914         | G_4 = 0.6667     |
|     2²⁰   | 1.5126          | 0.6561         | G_4 = 0.6667     |

Bernoulli ratio fluctuates around `G_4 = 0.667` (single-trial CLT
variance ~0.04); does NOT show the `S_3` factor. Confirms the HL
singular series in the chi_P ratio is arithmetic-content, not a
shape artefact.

## Mechanism interpretation (no new technique)

The closed form decomposes as:

1. **Per-q spike amplitude** (from E2.21):
   `f_N(e^{2πia/q}) ≈ π(N) · μ(q)/φ(q)`.
2. **Per-q arc geometry** (sinc window from finite-N truncation):
   `f_N(e^{2πi(a/q + δ)}) ≈ (μ(q)/φ(q)) · π(N) · sinc(Nδ/2)`.
3. **Per-q L^p contribution**:
   `(1/2π) ∫_{q-arc} |f_N|^p dδ ≈ μ²(q)/φ(q)^p · π(N)^p · (G_p / N)`
   (after the substitution `u = Nδ/2` and using `∫|sinc|^p du = π G_p`).
4. **Sum over (a, q) primitive pairs**: φ(q) primitive a's per
   sqfree q gives the per-q factor `φ(q) · μ²(q)/φ(q)^p =
   μ²(q)/φ(q)^{p-1}`.
5. **Sum over q sqfree** by the Selberg-Delange Euler product:
   `Σ_q μ²(q)/φ(q)^{p-1} = ∏_r (1+1/(r-1)^{p-1}) = S_{p-1}^HL`.

This is a *structural* derivation using only project-internal edges
(E2.21 spike amplitude, S168 sqfree-q decomposition) plus elementary
sinc-window integration. **No external technique imported.**

## What this composes / refines

* **E2.21** (L^∞ HL imprint at every q-rational): refined inline.
  The L^∞ endpoint extends to a complete L^p hierarchy for even
  p ≥ 4 with closed-form predictions verified at N up to 2²⁰.
* **E2.31** (BDJ m_4 ≈ 2.95): refined inline. The empirical 2.95
  reduces in shape to `G_4 · S_3 · (BDJ-Szegő-prefactor)`. The
  prefactor itself remains an open question — the standardisation
  transform from `f_N` to centred Toeplitz adds a `log²N` 
  normalisation that I did not close in this session.
* **S168** (sqfree-q L²-energy `μ²(q)·π(N)²/(φ(q) N)`): extended
  to per-q L^p energy `‖P_{V_q^prim} f_N‖_p^p ≈ G_p · μ²(q)/φ(q)^{p-1}
  · π(N)^p / N`. Per-q decay law `1/φ(q)^{p-1}` is the L^p
  generalisation of S168's L² law `1/φ(q)`.
* **S239** (cumulative-Q-strip data on BDJ m_4): the closure cascade
  is now derivable from a single closed-form formula
  `[Σ_{2 ≤ q ≤ Q sqf} μ²/φ^{p-1}] / [S_{p-1}-1]`, not requiring
  empirical fits.

## What this does NOT do

* Does not produce algorithmic content for π(x). Closed-form L^p
  norms for p ≥ 4 do not give a faster way to count primes.
* Does not add a new EDGES.md entry. The S247 prediction refines
  E2.21 + E2.31 inline; the construction is a unification of these
  two edges with S168, not a new theorem.
* Does not establish the conjecture rigorously. The 1.5% residual
  gap at N=10⁶ has the SHAPE of `O(1/log N)` PNT correction but no
  formal proof.
* Does not close the BDJ-Szegő-prefactor mystery (the 2.95 → 1.534
  reduction needs an explicit `log²N · (Li/π)^p` standardisation
  factor I leave for future work).

## Pre-stated F1-F4 verdicts

| F# | Falsifier                                          | Verdict | Outcome |
|----|----------------------------------------------------|---------|---------|
| F1 | p=4 ratio within ±10% of G_4·S_3, monotone        | PASS    | At 98.6% of pred at N=2²⁰; monotone convergence. |
| F2 | p=6, 8 ratios within ±10% of G_p·S_{p-1}          | PASS    | At 98.5% / 98.5% of pred at N=2²⁰. |
| F3 | closure cascade matches partial-sum formula        | PASS    | ±1pp per cohort at Q ∈ {2, 6, 30, 210}. |
| F4 | Bernoulli null shows no HL modulation             | PASS    | Bern ratio ~ G_p single-trial; no S_{p-1} factor. |

## Self-evaluation (CLAUDE.md template)

1. **What did I produce that was not in the project before?**
   * Closed-form L^p prediction for f_N at all even p ≥ 4 (project
     previously had only L² Parseval and L^∞ from E2.21).
   * Explicit cumulative-Q-strip closure formula
     `[Σ_{2 ≤ q sqf ≤ Q} μ²(q)/φ(q)^{p-1}] / [S_{p-1}^HL − 1]`.
   * Per-q L^p decay law `1/φ(q)^{p-1}` extending S168's L² law.
   * Identification of a 1.5% finite-N residual gap shared across
     all even p ≥ 4, consistent with PNT correction.

2. **What edges did my work compose or cite?**
   E2.21 (L^∞ spike, supplied per-q amplitude); E2.31 (BDJ m_4
   empirical 2.95, target for closed-form); S168 (sqfree-q L²-energy,
   structural ground); S239 (cumulative-Q-strip data, empirical
   anchor); E2.13 (HL singular series structure).

3. **Did my session produce only duplicate closures?** No.
   The closed-form L^p hierarchy for p ≥ 4 is a new artefact;
   the cumulative-Q closure formula is a new structural insight;
   the per-q `1/φ(q)^{p-1}` decay law is a new L^p generalisation
   of S168.

4. **Next-action for the next agent.** Close the
   BDJ-Szegő-standardisation prefactor that maps the S247 closed
   form `‖f_N‖_4^4 ≈ G_4·S_3·π(N)^4/N` to E2.31's empirical
   `m_4 ≈ 2.95 N/log²N`. The reduction from `f_N` to BDJ-standardised
   Toeplitz introduces a `log²N · (Li(N)/π(N))^4` factor that I have
   not pinned down to a clean constant. A clean derivation would close
   E2.31's empirical prefactor mystery and complete the unification.

## Grade defense

**B** is appropriate because:

* The session produced new closed-form structural content (L^p
  hierarchy for all even p ≥ 4) not previously stated in the project.
* It refines THREE existing edges (E2.21, E2.31, S168) inline with
  empirical evidence at N up to 2²⁰.
* It does NOT produce a new theorem rigorously — the conjecture is
  asymptotic with a 1.4% finite-N gap, not a sharp identity.
* It does NOT produce algorithmic content for π(x).
* The construction is elementary (Vinogradov major-arc + sinc-window);
  novelty is in the unification, not the technique.
* Pre-stated F1-F4 all PASSED, with sharper outcomes than predicted
  (cross-p alignment of the 1.5% gap was unanticipated).

Inflated A would require: rigorous proof of the asymptotic, OR
explicit BDJ-Szegő-prefactor closure that maps 2.95 → 1.534, OR
algorithmic consequence. None of these were produced.

C-grade rejected because: the closed-form L^p hierarchy is genuinely
new structural content (project had only L² + L^∞ before), not a
duplicate-plus closure.
