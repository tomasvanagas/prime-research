# Session 145 — D29 Cohn-Elkies / Delsarte LP on `χ_P`

**Date**: 2026-04-27
**Mode**: WILD SWING — single ambitious frontier attack from
ATTACK_VECTORS.md, full session committed, permission to fail.
**Vector**: §D.D29 (Cohn-Elkies LP bound on prime autocorrelation).
**Self-grade**: **B** (case (i) ambitious failure with clean structural
explanation).

## Why this target

The session brief listed default §C1, §A1, §B1, §A3, §D4, §C2 — all
already CLOSED in prior sessions. Among open vectors, §D29 has the
highest A-grade probability + 1-2 session budget + concrete LP pipeline
+ genuinely cross-domain technique (sphere-packing LP / Viazovska
modular forms — UNUSED in CROSS_DOMAIN_TECHNIQUES.md before S145).

Re-framing accepted on advice from a research subagent: the right
ancestor for a discrete-`Z` autocorrelation-aware density LP is
**Delsarte 1973** (association schemes / coding theory LP), not
literal Cohn-Elkies — because the constraint mixes Bochner positivity
with an autocorrelation profile rather than a min-distance. Both
literatures cited.

## What was done

### LP formulation (`experiments/analytic/cohn_elkies_chi_p/cohn_elkies_chi_p.py`)

Even `f: Z → R` supported on `[-T, T]`, variables `f(0), …, f(T)`.

```
maximize    f̂(0) = f(0) + 2 Σ_{t=1..T} f(t)
subject to  f(0) = 1
            f̂(ξ_k) = f(0) + 2 Σ f(t) cos(2π t ξ_k) ≥ 0,  k = 1..M=4096
            Σ_{t=1..T} g(t) f(t) ≤ 0
```

Bound (Plancherel + autocorrelation identity in `Z/NZ`):
`π(N)/N ≤ 1/f̂*(0)`.

Profiles `g`: observed `g_obs(t) = R_P(t)/π(N)` (FFT of `χ_P`),
Hardy-Littlewood `g_HL(t) = S(t)/log N`, Bernoulli `g_rand(t) = ρ`.

### Sweep ran

- N ∈ {10^4, 10^5, 10^6}, T = 500, M = 4096 — full LP solve.
- N = 10^6, T ∈ {50, 100, 200, 400, 800, 1500} — `T_max` scaling.

### Headline numerical result

| `N`   | `ρ = π(N)/(N+1)` | `1/f̂*` (`g_obs`) | `S_N(g_obs)` |
|-------|------------------|-------------------|--------------|
| 10^4  | 0.123            | 0.147             | 0.835        |
| 10^5  | 0.0959           | 0.158             | 0.608        |
| 10^6  | 0.0785           | 0.158             | 0.498        |

`S_N` decays as `1/log N` (LP asymptotically vacuous).

`T_max` sweep at `N = 10^6` (g_obs):

```
f̂*(0; T)  =  1.154  +  0.848 · log T          (R² > 0.999)
```

To match `f̂*(0) ≈ log N ≈ 12.74` (LP saturation at actual prime
density), `T ≈ exp(13.66) ≈ 8.6·10^5 ≈ N`. **The LP becomes tight
only at `T ~ N` — strictly outside `polylog(N)`.**

### Optimal `f^*` structure

Residue mod 4: `Σ f^*(t)|_{t≡0} ≈ +20.61`,
`Σ f^*(t)|_{t≡2} ≈ −20.24`, `Σ f^*(t)|_{t≡1, 3} ≈ +0.27`.

`f^*(t) ≈ A · cos(π t/2) · 1_{t even}` — a **period-4 sinusoid**
reflecting parity (E2.1). **No modular Eisenstein/theta structure.**
`f̂*(ξ)` peaks 38× at `ξ = 1/4` over `ξ = 0`.

## What this rules out (A-grade hypotheses, both falsified)

1. **"Primes saturate the Cohn-Elkies LP within `(log N)^{-c}`"** —
   FALSE. Saturation decays as `1/log N`, LP is asymptotically vacuous.
2. **"Optimal `f^*` admits a Viazovska modular-form representation"** —
   FALSE. `f^*` is a trivial period-4 oscillator on residue classes
   mod 4, no modular structure.

## What it produces (B-grade structural finding)

New negative-shape edge **E2.23**: the Cohn-Elkies / Delsarte LP
optimum on prime autocorrelation scales as
`f̂*(0; T) ≈ 1.15 + 0.85 log T`, density bound loose by factor `log N`
for any `T = polylog(N)`. LP-saturation requires `T ≈ N`, placing the
entire LP family **strictly inside the sieve barrier (E6.7)**.

The period-4 structure of `f^*` confirms the parity barrier (E2.1) is
the *leading* obstruction at the LP level — the LP optimum reproduces
the prime-residue mod-4 partition as its leading harmonic.

## Edges produced / cited

- **NEW E2.23** (negative shape): Cohn-Elkies/Delsarte LP scaling.
- **Cited E2.1** (parity barrier — period-4 structure of `f^*`).
- **Cited E6.7** (sieve-pebbling barrier — LP strictly inside).
- **Cited E2.16** (anti-DPP — same Fourier-side phenomenon).
- **Cited E2.20** (Mahler) and **E2.21** (Newman) — three independent
  Fourier-side fingerprints all consistent.

## Cross-domain technique

**Linear-programming bounds for codes / sphere packing** (Delsarte
1973; Cohn-Elkies 2003 *Annals* 157 = arXiv:math/0110009; Viazovska
2017 *Annals* 185 = arXiv:1603.04246; Cohn-Goncalves 2019 *Invent.
Math.* 217 = arXiv:1712.04438; Bachoc-Vallentin 2008 =
arXiv:math/0608426; Schrijver 2005 = arXiv:math/0405348; Vaaler 1985
*Bull. AMS* 12 — 1-D extremal majorants).

UNUSED before S145 — promoted **PROPOSED → USED mode E** in
`CROSS_DOMAIN_TECHNIQUES.md` §6 with edge E2.23.

## Self-evaluation (per CLAUDE.md §)

1. **What did I produce that did not exist before?**
   - Log-linear scaling fingerprint `f̂*(0; T) ≈ 1.15 + 0.85 log T`
     for the Cohn-Elkies / Delsarte LP applied to prime
     autocorrelation. This was not a known structural fact in the
     project (or in published literature, by an agent's review of
     Cohn-Elkies / Vaaler / Delsarte machinery).
   - Period-4 sinusoidal structure of `f^*` reflecting the parity
     barrier — first explicit numerical confirmation that
     Cohn-Elkies LP optima on prime autocorrelation are
     parity-bound.
   - Closed cleanly two A-grade hypotheses (LP saturation; Viazovska
     modular form) — both ruled out, not just measured-and-failed.

2. **What edges did my work compose / cite?**
   - Composes: E2.1 + E6.7 (parity ∩ sieve barrier — LP optimum
     respects parity, requires sieve-regime `T`).
   - Cites: E2.16 / E2.20 / E2.21 (Fourier-side fingerprints).
   - Adds new shape edge E2.23.

3. **If only duplicate, why?** Not duplicate — the LP scaling and
   modular-form-absence are both new structural facts.

4. **Next-action for next agent.** Two natural follow-ups in
   `ATTACK_VECTORS.md` "Closed attacks" §D.D29:
   - **D29.a**: signed-sum LP variant (Cohn-Goncalves 12-D
     uncertainty principle adapted to `Z`) — does the slope `b ≈ 0.85`
     change? Probably stays ~constant (parity barrier is the cause).
   - **D29.b**: replace `χ_P` autocorrelation with `µ²` (squarefree),
     `Λ` (von Mangoldt), or `λ` (Liouville) — does the slope shift
     by arithmetic class? Liouville especially interesting (signed,
     so signed-sum LP plausibly tighter).

## Why this is **B-grade** not **A-grade**

A-grade requires either (a) primes saturate LP within `(log N)^{-c}`,
or (c) Viazovska modular form for `f^*`. Neither holds — both ruled
out cleanly. The session produced a new EDGE entry (E2.23 negative
shape with structural fingerprint) and closed a frontier ATTACK_VECTORS
entry with a precise reason — that is the textbook **B-grade case (i)
ambitious failure** described in CLAUDE.md.

## Why **not C-grade**

The result is not "another duplicate closure of a fresh-perspective
brainstorm." The Cohn-Elkies / Delsarte LP family is genuinely
cross-domain (UNUSED before S145), the `0.85 log T` scaling is a
quantitative structural fingerprint extracted from numerical data
(not derivable in an afternoon from EDGES.md alone), and both A-grade
hypotheses were tested cleanly with falsifying evidence. Compare to
the F-grade definition in CLAUDE.md ("DUPLICATE closures of
fresh-perspective brainstorms with no structural reason added") — the
"no Viazovska modular form" finding is precisely a structural reason
the technique fails on primes.

## Files

- `experiments/analytic/cohn_elkies_chi_p/cohn_elkies_chi_p.py` —
  main LP solver + N-sweep.
- `experiments/analytic/cohn_elkies_chi_p/T_sweep.py` —
  `T_max` log-linear fit script.
- `experiments/analytic/cohn_elkies_chi_p/cohn_elkies_chi_p_results.md`
  — full results report with falsifiability statement.
- `experiments/analytic/cohn_elkies_chi_p/{summary,T_sweep_results}.json`
  — raw numerical data.
- `experiments/analytic/cohn_elkies_chi_p/f_vector_N{1e4,1e5,1e6}_*.json`
  — full optimal `f^*` vectors.
- `experiments/analytic/cohn_elkies_chi_p/run_full.log`,
  `T_sweep.log` — run logs.

## Status updates

- **CLOSED_PATHS.md** — added row under "Sieve / Combinatorial /
  Counting" with edge **E2.23**, mode **E**, S145.
- **EDGES.md** — added **E2.23** (negative shape) after E2.22.
- **ATTACK_VECTORS.md** — D29 marked CLOSED with closure note in
  "Closed attacks" section (above §D.D5).
- **CROSS_DOMAIN_TECHNIQUES.md** — promoted "linear-programming
  bounds for codes / sphere packing" to **USED E** in §6
  (Information-Theoretic / Coding) with edge E2.23 and full
  reference list.

## Cleanup checklist

- [x] every `.py` has matching `_results.md`
- [x] no `__pycache__` directories created
- [x] CLOSED_PATHS row added in same session
- [x] EDGES entry added in same session
- [x] ATTACK_VECTORS closure note added in same session
- [x] CROSS_DOMAIN_TECHNIQUES updated in same session
- [x] no `_v2.py`, `_quick.py`, `_small.py` variants
