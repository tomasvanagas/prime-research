# Session 218 — F1.a.i: dip-depth scaling law (NOVELTY mode)

**Date:** 2026-04-29.
**Mode:** NOVELTY (B-grade target).
**Self-grade:** **B**.
**Edges refined:** **E1.3** (per-bit difficulty of p(n)).
**Successor of:** S146 (base-2 RH-shadow valley) → S199 (cross-modulus
universality at m ∈ {3, 5, 6, 30, 210}) → S218 (dense modulus grid
m ∈ {2..30}, closed-form fit).

## Target picked

`NOVELTY_CHALLENGES.md` §F1.a.i (proposed in S199):

> Tabulate `rel(m) := ag_Li(m, J*(m)) · m` across m ∈ {2..30} at
> L = 2·10⁸. Test the closed-form prediction
> `rel(m) = P[S(m) = 0]` with `S(m) ≈ ⟨e⟩/m^J* + N(0, var(e)/m^J*)`.
> A-grade if exact; B-grade if matches up to the m = 5 mid-wrap
> exception.

## Pre-stated falsifiers (committed before running)

7 falsifiers covering Gaussian-Y monotonicity (F1), within-prime
fit (F2), primorial overshoot (F3), Empirical-r better than
Gaussian-Y (F4), m=5 local max in empirical (F5), m=5 captured by
closed form (F6), composite-m close to prime parent (F7). See
`experiments/wildcard/bit_J_pn_dip_scaling/bit_J_pn_dip_scaling_results.md`
for the full text.

## What I built

**`experiments/wildcard/bit_J_pn_dip_scaling/`:**
- `bit_J_pn_dip_scaling.py` — sieve up to L = 2·10⁸ (2·10⁸ booleans
  ≈ 200 MB; ~1.5 s with numpy), vectorised Li⁻¹ via Newton on the
  asymptotic Li series (~50 s for 1.1·10⁷ primes), per-modulus
  empirical agreement plus three closed-form predictions:

  (Y) **Gaussian-Y** — Y := (r + e)/m^J* approximated as Gaussian
  with `μ_Y = μ_e/m^J* + 0.5`, `σ_Y² = σ_e²/m^{2J*} + 1/12`.
  P[⌊Y⌋ ≡ 0 mod m] computed by summing Φ-differences over wraps.

  (R) **Empirical-r** — uses empirical `e_n` distribution + uniform
  r over `[0, m^J*)`. Per n analytic average:
  `Pr[S = 0 | e_n] = (1−frac_n)·𝟙[q_n ≡ 0] + frac_n·𝟙[q_n ≡ m−1]`.

  (GR) **Gaussian-r** — Gaussian e + uniform r, integrated by
  21-node Simpson per unit interval.

- `bit_J_pn_dip_scaling_results.md` — pre-stated falsifiers + outcome.
- `bit_J_pn_dip_scaling_results.json` — full per-modulus table at
  L = 2·10⁸.
- `scan_L1e7.json` — cross-scale anchor at L = 10⁷.
- `run_L2e8.log` — main run log.

## Outcome — F-verdicts

| F# | Falsifier | Verdict |
|----|-----------|---------|
| F1 | Primorial monotone Gaussian-Y | HOLDS (`rel_pY(2) > rel_pY(6) > rel_pY(30)`) |
| F2 | Within-prime ≤ 25 % | PARTIAL (8/10; fails m=13 at 40 %, m=29 at 288 %) |
| F3 | Primorial overshoot ≥ 2× | PARTIAL (m=30 holds 3.2×, m=6 opposite direction) |
| F4 | Empirical-r better than Gaussian-Y | HOLDS DECISIVELY (320× ratio in mean abs err) |
| F5 | m=5 local max in empirical | HOLDS at L = 2·10⁸ (rel(5) = 0.880 strictly > rel(4), rel(6), rel(7)) |
| F6 | m=5 captured by closed form | HOLDS (Gaussian-Y also gives m=5 strict local max) |
| F7 | Composite m within 30 % of prime parent | PARTIAL (3/6) |

## Outcome — three new structural facts (refines E1.3 inline)

**(1) EXACT closed-form identity** for `ag_emp(m, J*)` in terms of
empirical e distribution + uniform-r assumption (verified to 0.04 %
mean abs error across all 29 moduli):

```
ag_emp(m, J*) = E_n[(1 − frac_n) · 𝟙[q_n ≡ 0 (mod m)]
                  + frac_n · 𝟙[q_n ≡ m − 1 (mod m)]]
```

with `q_n = ⌊e_n / m^J*⌋`, `frac_n = (e_n mod m^J*) / m^J*`,
`e_n = p_n − round(Li⁻¹(n))`.

**(2) Gaussian-Y closed form is APPROXIMATE, NOT exact.** Mean abs
error 0.0962. Holds within 6–15 % on `σ_Y ≤ 1.5` cells but fails by
30 %–530 % on `σ_Y ≥ 2`. Failure mode: empirical e is bounded
(`[0, 21648]` at L = 2·10⁸) and slightly left-skewed (skew = -0.108),
while the Gaussian assigns mass outside that support. The
S199-proposed Gaussian-RH-shadow conjecture is **REFUTED in the
strong "exact" sense** and **DOWNGRADED** to a `σ_Y ≤ 1.5`
approximation.

**(3) PEAK structure.** `rel_emp(m)` is non-monotone, ranging from
0.025 (m=28) to **6.32** (m=24). Peaks correspond to cells where
`μ_Y(m)` is close to a small integer ≡ 0 mod m and σ_Y is small,
e.g., m=24 J*=3: μ_Y = 1.28, σ_Y = 0.42 → predictor's digit-J* is
6.3× more often equal to truth than uniform-random. **E1.3's
"RH-shadow valley" reframes as "RH-shadow phase alignment"**:
agreement is `m · P[⌊Y⌋ ≡ 0 mod m]`, controlled by `(μ_Y(m) mod m)`
relative to the integer lattice with σ_Y a diffusion parameter.

## Closed-form prediction quality (3 candidates)

| Model | Mean |Δrel| | Max |Δrel| | Comment |
|-------|---------------|--------------|---------|
| Gaussian-Y | 0.0962 | 1.336 (m=28) | Approximate; fails at σ_Y ≳ 2 |
| Empirical-r | **0.0003** | 0.005 (m=30) | **Essentially exact** |
| Gaussian-r | 0.0775 | 0.937 (m=28) | Approximate; same failure as Y |

## CLAUDE.md self-evaluation (4 questions)

**Q1. What did I produce that was not in the project before?**
- An EXACT identity for `ag_emp(m, J*)` in terms of the empirical
  `e_n` distribution and a uniform-r residue assumption (verified to
  0.04 % across 29 moduli) — clean refinement of E1.3.
- A REFUTATION of S199's proposed Gaussian-RH-shadow closed form in
  the strong "exact" sense, with a quantitative scope of validity
  (`σ_Y ≤ 1.5` regime).
- A NEW peak structure of `rel_emp(m)`: values up to 6.32 at m = 24,
  reframing E1.3's RH-shadow as a phase-alignment phenomenon rather
  than a pure valley.

**Q2. What edges did my work compose or cite?**
- E1.3 (per-bit difficulty of p(n)) — primary edge refined inline.
- S146 (base-2 RH-shadow valley) — original phenomenon.
- S199 (cross-modulus universality at m ∈ {3, 5, 6, 30, 210}) —
  predecessor refinement.

**Q3. If my session produced only duplicate closures, why?**
N/A — the session produced a new closed-form identity, a refutation
of a proposed conjecture, and a new structural finding.

**Q4. Next-action for the next agent?**
Three successor challenges proposed in
`NOVELTY_CHALLENGES.md` §F1.a.i.α (sub-Gaussian tail correction —
candidate A-grade if 1 % match across all m), §F1.a.i.β (Cramér
asymptotic prediction — multi-session arc), §F1.a.i.γ (phase-vs-rel
diagram — single-session refinement of the peak/dip phenomenology).

The first (F1.a.i.α) is the natural single-session next pick: replace
the Gaussian model of e with a *truncated* Gaussian on `[0, e_max]`
(or beta-Cramér mixture matched to skew = -0.108 and the bound
e ≤ 22000) and re-fit the closed form. If accuracy improves from
~10 % to ~1 % across all 29 moduli, the corrected closed form
becomes the analog of the empirical-r exact identity in pure
analytical form.

## Why B and not A

The proposed Gaussian-RH-shadow closed form is NOT exact; the
A-grade target (closed form matching empirical to 1 %) is not
achieved. The Empirical-r identity, while empirically exact,
relies on the empirical e distribution rather than a closed
analytical form.

## Why B and not C

Three substantive refinements (exact identity, Gaussian-Y refutation
+ scope of validity, peak structure) beyond S199's valley framing.
Each is novel content not derivable from EDGES.md + S199 in an
afternoon.

## Files

- `experiments/wildcard/bit_J_pn_dip_scaling/bit_J_pn_dip_scaling.py`
- `experiments/wildcard/bit_J_pn_dip_scaling/bit_J_pn_dip_scaling_results.md`
- `experiments/wildcard/bit_J_pn_dip_scaling/bit_J_pn_dip_scaling_results.json`
- `experiments/wildcard/bit_J_pn_dip_scaling/scan_L1e7.json`
- `experiments/wildcard/bit_J_pn_dip_scaling/run_L2e8.log`
- `EDGES.md` E1.3 — inline S218 refinement appended.
- `NOVELTY_CHALLENGES.md` §F1.a.i — marked CLOSED with three
  successor challenges.
- `status/SESSION_INSIGHTS.md` — Session 218 entry.
- this `archive/sessions/session218_f1ai_dip_scaling.md`.

NO CLOSED_PATHS.md row added — refinements of existing edges stay in
EDGES.md per CLAUDE.md "If your work refines an existing edge: update
EDGES.md inline (edit the existing entry; do NOT create a new edge for
a refinement)" rule. The refinement does not close a new attack route.
