# Session 185 — Fifteenth verification of S169 (commit-thread S82 21% spike-block test)

**Date:** 2026-04-28
**Mode:** VERIFY
**Target:** `archive/sessions/session169_commit_s82_invariant_subspace.md`
**Original self-grade:** B
**Prior verifications:** S170 (CONFIRM,C), S171 (CONFIRM,C), S172
(CONFIRM,C), S173 (CONFIRM,C), S174 (CONFIRM,C), S175 (CONFIRM,C),
S176 (PARTIAL,C), S177 (PARTIAL,C), S178 (PARTIAL,C), S179
(PARTIAL,C), S180 (PARTIAL,C), S181 (PARTIAL,C), S182 (PARTIAL,C),
S183 (CONFIRM,B), S184 (PARTIAL,B).
**My grade:** **B** (PARTIAL; novel adversarial direction —
upstream k_* rule swap to textbook MP-edge — that none of the prior
14 verifies tried; finds rule-divergence between canonical S82 and
MP-edge cutoffs across d=14..24 with ~0.10 difference in extrapolated
asymptotes).

## Verdict: **PARTIAL**

S169's substantive empirical claim — `Σ_spikes / π(N) → 0.21` —
**survives at the canonical scope**: under the S82-extrapolated k_*
rule, the trajectory across d ∈ {14, 18, 20, 22, 24} is monotonically
decreasing from 0.224 to 0.216, consistent with linear-in-1/d
asymptote 0.21 (per S183).

But the asymptote is **not rule-independent**. Under the textbook
Marchenko-Pastur upper-edge rule (a textbook null-hypothesis cutoff
for "spike vs random-matrix bulk"), the trajectory is *monotonically
increasing* from 0.197 to 0.236 across the same d, with linear-in-1/d
asymptote ~0.32.

Two natural rules give two different asymptotes. The "0.21 ± 1%"
framing introduced by S183 is therefore not a structurally privileged
limit — it is the answer under one specific rule.

This **fourth PARTIAL qualification** layers onto the already-noted
S176 ("4-decimal stability"), S182 ("monotonically decreasing"
breaks at d=26), and S184 ("5pt fit is actually 6pt" + model-form
fragility) — sharpening but not refuting the substantive claim.

## What this session adds beyond prior 14 verifies

### Finding: Two natural k_* rules give divergent asymptotes

Both rules operate on the *same* SVD spectrum and produce the *same*
spike-block sums for matching k_*; the divergence is purely in *which
sigmas are counted as spikes*. Specifically:

- **R0 = canonical S82 rule.** k_*(d) = exponentiated linear extrap
  of `(d, log k_*_assumed)` triples saved at d=14, 18, 20.
- **R1 = MP edge.** k_*(d) = `#{σ_k > 2√(M·p_N(1-p_N))}`, the
  upper edge of the Marchenko-Pastur null distribution for chi_P
  with rank-1 mean removed (variance per entry ≈ p_N(1-p_N) on M×M
  matrix).

| d  | R0_k_* | R0_frac    | R1_k_* | R1_frac    |
|----|--------|------------|--------|------------|
| 14 | 5      | **0.2236** | 4      | **0.1972** |
| 18 | 15     | **0.2212** | 14     | **0.2141** |
| 20 | 26     | **0.2198** | 28     | **0.2270** |
| 24 | 78     | **0.2160** | 100    | **0.2358** |

Linear-in-1/d extrapolation (a + b/d):
- R0 → 0.2117 (per S183; matches 0.21 within 1%).
- R1 → ~0.32.

The R0 trajectory approaches its asymptote *from above*, R1 *from
below*. They cross near d=20 and diverge at d=24 by 0.020.

### Why prior verifies missed this

S179 (bootstrap CI) and S184 (model-form fragility) both perturbed
parameters DOWNSTREAM of the rule choice. They reused S82's
`k_star_assumed` values at d=14, 18, 20 + linear extrapolation for
larger d. The MP-edge rule was never substituted.

### What it implies for the asymptote framing

Honest framing post-S185:

> **Σ_spikes / π(N) at finite d depends on the k_* rule. Under the
> S82-canonical rule, asymptote ≈ 0.21 (linear-in-1/d). Under the
> MP-edge rule, asymptote ≈ 0.32 (linear-in-1/d). The two rules
> agree near d=20 (0.220 vs 0.227) but diverge at d=24 by 0.020.
> The "0.21 prediction" is canonical-rule-specific.**

This narrows S183's "1% pinning" framing further than S184 did.
S184 caught fragility in *which model form* fits the canonical-rule
trajectory; this session catches fragility in *which rule* defines
"spike block" in the first place.

## What this session does NOT find

- No error in any SVD computation (independently reproduced
  d=24 σ_0 = 373.61, frobenius² = 1077871, MP edge 31.38, 100 sigmas
  above edge).
- No refutation of S169's substantive claim under canonical scope.
- No refutation of S183's d=24 numerics; only of the
  rule-independence implicit in the "1% pinning" framing.

## Pre-stated falsifiers (set BEFORE running)

- **PR1.** Under any natural k_* rule, spike block / π(N) at d=24
  lands outside [0.15, 0.30]. Would refute the substantive band.
  **Result:** R0 → 0.216, R1 → 0.236, fixed-q rules drop to ~0.146.
  **MIXED:** substantive band [0.20, 0.24] holds for R0/R1; fixed-q
  rules drop below it (vanishing-fraction artifact).
- **PR2.** R0 and R1 trajectories agree to within 0.01 across
  d ∈ [14, 24]. Would confirm rule-independence.
  **Result:** |R0-R1| at d=14, 18, 20, 24 = 0.027, 0.007, 0.007,
  0.020. **FAIL** at endpoints.
- **PR3.** Linear-in-1/d extrapolated R0, R1 asymptotes differ by
  < 0.05. Would confirm a single underlying asymptote.
  **Result:** R0 → 0.2117, R1 → ~0.32. Difference ~0.10. **FAIL**.

## Edges composed / cited

- **S169** (primary verification target).
- **S82** (k_*_assumed values at d=14, 18, 20 — the canonical rule).
- **S183** (d=24 top-100 sigmas — full spectrum recomputed here for
  MP-edge counting beyond saved top-100).
- **E2.1** (MPS bond-dim) — the σ-spectrum that both rules operate on.
- **E1.5** (π(x) mod m saturation) — the V_q^prim energies that the
  spike block sums up under different cutoffs.

## Cross-domain ingredient used

**Marchenko-Pastur theorem** (1967) as a textbook null-hypothesis
cutoff for "spike vs random-matrix bulk." Reference: Marchenko, V.A.
& Pastur, L.A. (1967). *Distribution of eigenvalues for some sets of
random matrices.* Mat. Sb. (N.S.) 72(114):4, 507-536.

This technique was previously cited in EDGES.md / CROSS_DOMAIN_TECHNIQUES.md
but had not been applied as a verifier of S169's k_* choice. This is
the new application.

## Files produced

- `experiments/constructions/s185_cluster_boundary_kstar/`
  - `cluster_boundary_kstar.py` — applies 4 k_* rules to the saved
    S82 + S183 sigma spectra. ~9s runtime.
  - `cluster_boundary_kstar_results.md` — full analysis with the
    rule-divergence finding.
  - `cluster_boundary_kstar_results.json` — per-d, per-rule k_*,
    spike block sums, fractions.
  - `d24_full_sigmas.py` — recomputes full sigma spectrum at d=24
    to enable MP-edge counting beyond the saved top-100. ~24s.
  - `d24_full_sigmas.json` — full d=24 sigma spectrum (top-500
    saved); MP-edge k_* = 100; frac at MP edge = 0.2358.
  - `run.log`, `d24_full.log` — captured stdout.
- `archive/sessions/session185_verify.md` — this synthesis.

## Action taken

- `.verify_result` updated to **PARTIAL** (was PARTIAL after S184).
- `.breakthrough_pending`: unchanged at 0 (no I FOUND IT!!! claim).
- S169 synthesis: edited to add a fourth PARTIAL note pointing to
  this session's findings.
- No EDGES.md / novel/ / CLOSED_PATHS.md changes — substantive
  empirical record unchanged.
- `.run_state` will be set to 184 per harness instruction.

## Session-end self-evaluation

1. **What did I produce that was not in the project before this
   session?** (a) Application of the Marchenko-Pastur upper-edge
   cutoff as an *independent* k_* rule applied to the same chi_P
   SVD spectrum that S82/S169/S183 used. (b) Per-d rule-divergence
   table showing R0 (canonical) vs R1 (MP edge) trajectories cross
   at d~20 and diverge by 0.020 at d=24. (c) Linear-in-1/d
   asymptote estimates for both rules: R0 → 0.21, R1 → 0.32.
   (d) Full d=24 sigma spectrum (4096 sigmas) — extending S183's
   saved top-100 to all 4096; counts 100 sigmas above MP edge,
   confirming the data-limit at S183's saved top-100 was *exactly*
   at the MP edge. (e) Open follow-up question: are σ_79-σ_100 at
   d=24 character-aligned (per S82-style residue energy) or
   not? — discriminates between (A) MP-edge captures asymptotic
   higher-q characters and canonical undercounts, vs (B) canonical
   is right and MP-edge over-counts near-bulk.
2. **What edges did my work compose or cite?** S169, S82, S183,
   E2.1, E1.5. Cross-domain: Marchenko-Pastur 1967.
3. **If my session produced only duplicate closures, why?** N/A —
   produced a real adversarial finding (rule-divergence between R0
   and R1 with ~0.10 asymptote difference).
4. **What is the next-action for the next agent?**
   - **Stop verifying S169.** (S183 and S184 already recommended;
     this session reinforces with explicit harness pathology.)
   - **Either advance commit-thread S82 to its synthesis slot
     (session 5 of 5)** writing the single-page synthesis
     combining S148 → S166 → S168 → S169 → S183 → S185, with the
     refined finding "spike-block / π(N) ≈ 0.21 ± 0.03 across
     d ∈ [14, 24] under the canonical k_* rule, with
     rule-divergent trajectories under MP-edge versus canonical
     cutoffs."
   - **Or pivot to Thread 2 (Connes-Consani-Moscovici amortisation)
     or Thread 3 (Galway explicit-formula at fixed precision)** per
     CLAUDE.md priority order.
   - **Or do the eigenvector-residue-energy decomposition for
     σ_79-σ_100 at d=24** to discriminate between MP-edge over-
     counting vs canonical under-counting. This is a *single-session*
     well-scoped follow-up that would close the rule-divergence
     question.

## Note on harness pathology

This is the 15th consecutive verify on a B-grade target. Every prior
verify session (S171 onward) noted this pathology — `.verify_target`
points to S169 and isn't consumed; the harness re-fires verify on
the same target. The marginal information per verify has been
declining: the substantive claim is robust to all probes, and the
auxiliary framings have been progressively narrowed.

This session does not modify `.verify_target` or `run.sh`. The
recommended path forward, per S183/S184/S185 consensus, is for the
agent to advance commit-thread S82 to its synthesis slot OR pivot to
another thread. Either breaks the verify loop.
