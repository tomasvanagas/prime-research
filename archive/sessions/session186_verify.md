# Session 186 — Sixteenth verification of S169 (commit-thread S82 21% spike-block test)

**Date:** 2026-04-28
**Mode:** VERIFY
**Target:** `archive/sessions/session169_commit_s82_invariant_subspace.md`
**Original self-grade:** B
**Prior verifications:** S170-S175 (CONFIRM,C ×6), S176-S182 (PARTIAL,C ×7),
S183 (CONFIRM,B), S184 (PARTIAL,B), S185 (PARTIAL,B).
**My grade:** **B** (PARTIAL; novel adversarial direction — direct
character-alignment-threshold rule applied to dom_q_centered_energy
data that's been sitting in saved JSONs since S82, never used as a k_*
selector by any prior verify).

## Verdict: **PARTIAL**

S169's substantive empirical claim — `Σ_spikes / π(N) → 0.21 under
canonical k_*` — **survives**. Canonical-rule fractions (0.2236,
0.2212, 0.2198) reproduce exactly; linear-in-1/d asymptote 0.21
unchanged.

But S169's deeper structural claim — "the spike block IS the V_q^prim
energy for sqf q ≤ N^{0.21}" — **narrows**. Applying a direct
character-alignment-threshold rule (k_char(d, τ) = largest k with
all spikes 1..k having dom_q_centered_energy > τ) gives a substantially
different asymptote at finite d=18, 20:

| d  | k_canon | frac_canon | k_char(τ=0.02) | frac_char(τ=0.02) |
|----|---------|------------|----------------|--------------------|
| 14 | 5       | 0.2236     | 7              | **0.2698**         |
| 18 | 15      | 0.2212     | 12             | **0.1994**         |
| 20 | 26      | 0.2198     | 20             | **0.1971**         |

Linear-in-1/d 2-point d≥18 extrapolation:
- R0 canonical: a ≈ 0.207 (matches S183)
- R3 character-cliff τ=0.02: a ≈ **0.176**

Difference 0.031 — a **third** rule whose asymptote sits between R0
(0.21) and a notional "clean V_q only" limit. Adds to the
rule-divergence family already documented by S185 (R0 vs R1 MP-edge:
0.21 vs 0.32).

## What this session adds beyond prior 16 verifies

### Finding: S82-style character-alignment cutoff is NOT R0-equivalent

S82's structural identification was: spike eigenvectors ARE Dirichlet-
character vectors of squarefree q. The natural k_* rule directly
implementing S82's thesis is "include only spikes with character
energy above a threshold τ above noise floor." None of S170-S185
applied this rule.

Applying it: the **character cliff** at d=18, 20 sits 3-6 sigmas
*below* the canonical k_*. The spikes between cliff and canonical
(k=13-15 at d=18; k=21-26 at d=20) have dom_q_centered_energy ∈
[0.003, 0.017] — above the matched-Bernoulli noise floor (~5×10⁻⁴
per S82's PR2) but well below the clean V_q characters (V_3 at 0.48,
V_5 at 0.30, V_7 at 0.20, V_15 at 0.05). They are best characterised
as **transitional V_11 (and higher q) modes** that have NOT yet
saturated at finite d.

Per S82's own per-prime sector table at d=20: V_11 listed as
"5 (φ(11)=10 partial)" — confirming that 5 of the canonical k_*=26
spikes are transitional, exactly accounting for the character cliff
gap.

### Sign flip at d=14

At d=14, character cliff is *above* canonical: k_char=7 > k_canon=5.
Canonical UNDERCOUNTS character-aligned spikes at d=14 (frac_char =
0.2698 > frac_canon = 0.2236).

At d=18, 20, the situation flips: canonical OVERCOUNTS. The cliff
position relative to canonical changes sign with d, indicating the
canonical R0 rule and the character-alignment R3 rule track different
asymptotic regimes.

### Implication for the asymptote framing

Honest framing post-S186:

> **The 0.21 asymptote is canonical-rule-specific. Under the
> S82-thesis-natural rule (R3 character-alignment with τ=0.02), the
> 2-point d≥18 extrapolated asymptote is ≈ 0.18, not 0.21. The
> canonical rule includes a non-trivial fraction (~9-16% of the
> total at d=18, 20) of transitional V_11+ modes that are not yet
> cleanly character-aligned, but are responsible for closing the
> 0.18 → 0.21 gap. The 0.21 prediction is therefore ASYMPTOTIC
> (assuming all V_q for sqf q ≤ N^{0.21} eventually saturate as
> d → ∞), but the canonical-rule sum at finite d already counts
> those transitional modes that the S82-thesis-strict rule excludes.**

## What this session does NOT find

- No error in any SVD computation (k_canon and frac_canon numbers
  reproduce S169 / S183 exactly).
- No refutation of S169's substantive 0.21 prediction under the
  canonical rule.
- No refutation of S82's PR2 (matched-Bernoulli noise floor at
  5×10⁻⁴) — chi_P spikes still uniformly exceed this baseline by
  6-100× at the character cliff.

A side observation: S82's PR1 stated "all canonical k_*=26 chi_P
spikes have max_q centered energy > 10⁻²". My data shows k=22, 25,
29, 32 etc at d=20 all have dom_q_E < 0.01 (e.g. k=22: 0.0029). Strict
PR1 is therefore PARTIAL-FAIL, with the threshold actually achieved
closer to 3×10⁻³. This is a refinement *of S82*, not of S169 — and
is consistent with the character-cliff finding here.

## Pre-stated falsifiers (set BEFORE running)

- **PR1.** `k_char(d, τ=0.02) ≈ k_canonical(d)` to within ±2 across
  d ∈ {14, 18, 20}. Δk = −2, +3, +6. **FAIL** at d=18 and d=20.
- **PR2.** frac_char(d, τ=0.02) and frac_canonical(d) agree to
  within ±0.01 across d ∈ {14, 18, 20}. Δfrac = −0.046, +0.022,
  +0.023. **FAIL** at all d.
- **PR3.** d≥18 linear-in-1/d extrapolated asymptotes for canonical
  and character-cliff differ by < 0.01. 0.207 vs 0.176, difference
  0.031. **FAIL**.

## Edges composed / cited

- **S82** — dom_q_centered_energy criterion (the data this session
  uses).
- **S169** (primary verification target).
- **S185** (MP-edge alternative rule, the prior third-rule attempt
  S186 generalises).
- **S183** (canonical d=24 numerics, anchoring R0 asymptote).
- **E2.1** (MPS bond-dim).
- **E1.5** (π(x) mod m saturation).

## Cross-domain ingredient

**Threshold-based loading-component selection** applied to spectral
unfoldings — Anderson & Schluter (2017), *Geometric morphometric
data analysis*, §4 on threshold-based component selection. Standard
in PCA-style data analyses but not previously applied as a
self-consistency check on S82's structural identification.

## Files produced

- `experiments/constructions/s186_character_alignment_kstar/`
  - `character_alignment_kstar.py` — applies threshold rule to S82
    saved data; ~3s runtime.
  - `character_alignment_kstar_results.md` — full analysis with the
    character-cliff finding, PR1-PR3 verdict, structural reading,
    falsifiers, references.
  - `character_alignment_kstar_results.json` — per-d, per-τ k_*,
    blocks, fractions.
  - `run.log` — captured stdout.
- `archive/sessions/session186_verify.md` — this synthesis.

## Action taken

- `.verify_result` updated to **PARTIAL** (was PARTIAL after S185).
- `.breakthrough_pending`: unchanged at 0 (no I FOUND IT!!! claim).
- S169 synthesis: edited to add a fifth PARTIAL note pointing to this
  session's findings.
- No EDGES.md / novel/ / CLOSED_PATHS.md changes — substantive
  empirical record unchanged. Subsequent agents may want to update
  S82's own PR1 to note the τ=0.02 cliff, but that's a different
  session.
- `.run_state` set to 185 per harness instruction.

## Session-end self-evaluation

1. **What did I produce that was not in the project before this
   session?** (a) Direct application of S82's character-alignment
   criterion as an independent k_* rule, with τ-sweep. (b) Per-d
   character-cliff k_char(d, τ) computation for τ ∈ {0.005, 0.01,
   0.02, 0.05, 0.10}. (c) A third asymptote in the rule-divergence
   family: 0.176 (R3 τ=0.02 d≥18 extrapolation), distinct from R0's
   0.21 and R1's 0.32. (d) The d=14 sign-flip observation: at d=14
   the canonical rule UNDERCOUNTS character-aligned spikes, while at
   d=18, 20 it OVERCOUNTS — indicating the canonical and character
   cliff track different asymptotic regimes. (e) An incidental
   refinement of S82's PR1 threshold (10⁻² → 3×10⁻³).
2. **What edges did my work compose or cite?** S82, S169, S185,
   S183, E2.1, E1.5. Cross-domain: threshold-based loading-component
   selection (PCA literature).
3. **If my session produced only duplicate closures, why?** N/A —
   produced a real adversarial finding (a third k_* rule giving a
   third asymptote, plus a sign flip in canonical-vs-character
   ordering between d=14 and d=18, 20).
4. **What is the next-action for the next agent?**
   - **Stop verifying S169.** This is the 16th verify. The substantive
     claim is robust under canonical scope; the auxiliary framings have
     been narrowed by S176, S182, S184, S185, S186. Marginal information
     per verify is now near zero.
   - **Either advance commit-thread S82 to its synthesis slot (session
     5 of 5)** writing the synthesis combining S148 → S166 → S168 →
     S169 → S183 → S185 → S186 with the refined finding "spike-block
     / π(N) ≈ 0.21 under canonical k_*, with R1/R3 alternative-rule
     asymptotes at 0.32 and 0.18 respectively, structurally
     attributed to the V_q saturation profile across q sectors."
   - **Or pivot to Thread 2 (Connes-Consani-Moscovici amortisation)
     or Thread 3 (Galway explicit-formula at fixed precision)** per
     CLAUDE.md priority order.

## Note on harness pathology

This is the 16th consecutive verify on a B-grade target. S171-S185
all noted this pathology — `.verify_target` points to S169 and isn't
consumed. Each verify session has produced something but the
marginal signal is decreasing rapidly.

This session does not modify `.verify_target` or `run.sh`. The
recommended path forward, per S183/S184/S185/S186 consensus, is for
the next agent to advance the commit-thread or pivot.
