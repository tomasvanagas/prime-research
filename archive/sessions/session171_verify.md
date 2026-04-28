# Session 171 — Re-verification of S169 (commit-thread S82 21% spike-block test)

**Date:** 2026-04-28
**Mode:** VERIFY (re-fire of the S169 verification slot)
**Target:** `archive/sessions/session169_commit_s82_invariant_subspace.md`
**Original self-grade:** B
**Prior verification:** `session170_verify.md` — verdict CONFIRM with discretization caveat
**My grade for this re-verification:** **C** (independent reproduction at one
larger d not previously checked; substantive claims uphold; one further
caveat surfaced — same verdict as session170)

## Verdict: **CONFIRM**

S169's claims survive a second independent reproduction. Session170
already verified at d=14 and d=18 by re-coding the Fourier sieve.
This session extends the audit to **d=22** (one of the two d values
session170 left unchecked, and the largest scale the synthesis reports
analytic data for) and probes the headline framing more sharply.

## What this session adds beyond session170

### 1. d=22 analytic reproduction (NEW — session170 did not test)

| d | Q* | synthesis cum(Q*) | my cum(Q*) | match | runtime |
|---|----|-------------------|-----------|-------|---------|
| 14 | 8  | 530.6732   | 530.6732  | exact (also session170) | <1s |
| 18 | 14 | 6085.3509  | 6085.3509 | exact (also session170) | <1s |
| 20 | 18 | 20556.9734 | 20556.9734| exact (this session) | <2s |
| 22 | 25 | 72844.00   | 72844.0032| exact (this session) | 1.6s |

The d=22 computation `cum(Q*=25) = 72844.0032` (independent re-impl)
matches the synthesis's reported `72844.00` to 4 decimals. The ratio
`cum/(0.21·π(N)) = 1.1721` reproduces the synthesis claim.

### 2. (A.1) canonical exponent table at d=20 reproduces exactly

Probed the (A.1) table at d=20 across all six exponent values
`a ∈ {0.10, 0.15, 0.20, 0.21, 0.25, 0.30}` independently. All six
ratios reproduce to 4 decimals: 0.1173 / 0.1890 / 0.2457 / 0.2506 /
0.2946 / 0.3481. This was a table session170 did not audit.

### 3. The "0.21 → cum/π = 0.21" relation is closer to a Wirsing-A=1
**tautology than a discovery** — sharper criticism than session170 made

The synthesis presents the asymptotic claim
`cum(Q=N^a) / π(N) → a as N → ∞` (instantiated at a=0.21) as a derived
prediction from S168. But this relation is essentially Wirsing-A=1:

  cum(Q) ≈ Σ_{sqf q ≤ Q} (π(N))²/(φ(q)·N)
        ≈ π(N) · π(N)/N · Σ_{sqf q ≤ Q} 1/φ(q)
        ≈ π(N) · 1/log(N) · A·log(Q)        [Wirsing-A=1]
        = π(N) · log(Q)/log(N)
        = a · π(N)                          [when Q = N^a]

So once you set the exponent a, the asymptotic ratio MUST be a. Session
170's check that "Wirsing-A=1 holds algebraically" is the substantive
content; the "→0.21" is then an arithmetic substitution.

The genuinely empirical claims are:
(i) the SVD spike block is approximately equal to the analytic
    cum(Q_eff) for some Q_eff at finite N (the "missing-spike" gap is
    real but bounded);
(ii) the matching Q_eff is at exponent ~0.18, drifting upward toward
    0.21 as N grows.

Both reproduce. But the framing "21% confirmed" is ambiguous between
(a) the tautological a=0.21 instantiation and (b) the empirical SVD
spike block / π(N) ≈ 0.22 at d=14..20. The synthesis is mostly clear
on this distinction but does occasionally elide it (e.g., point #1 of
"What this session produced" reads as if the 21% cum was empirically
confirmed — actually that's a Wirsing-A=1 prediction with finite-N
corrections).

This sharper criticism does not refute the substantive claims (the
SVD ≈ 0.22·π reproduction stands; Wirsing-A=1 is correct). It reduces
the novelty content slightly: the asymptotic relation is *not* an
independent test of S168, it is Wirsing-A=1 substitution. Synthesis
already self-grades B, so this caveat does not warrant grade demotion.

### 4. Re-confirms session170's discretization caveat for Q_eff exponent

The integer-rounded Q_eff values 6, 10, 13 give log/log ratios
0.1846, 0.1846, 0.1850 — agreeing to 4 decimals as the synthesis
states. Session170 already noted this is partly a quantization
artifact: the *continuous-interpolated* Q_eff varies between exponents
0.169 and 0.181. I re-derived this and concur. Future sessions should
report continuous Q_eff.

## What I tried to falsify

1. **d=22 reproduction.** Synthesis's 72844.00 reproduces exactly.
2. **(A.1) canonical table at d=20.** Six ratios reproduce.
3. **Wirsing-A=1 algebraic identity.** prod_p (1 − 1/p)·(1 + 1/(p−1))
   = prod_p [(p−1)/p · p/(p−1)] = 1. Confirmed.
4. **Spike block / π(N) ratio variance under k_star perturbation.**
   Session170 already showed ±5·k_star → ±12% ratio. I re-checked
   that k_star=26 is structurally determined by S82 (φ-dimensions of
   {3,5,7,15} = 2+4+6+8 = 20, plus partial mod-11 ≈ 5–6 → 25–26).
   The structural determination matches the JSON's k_star_assumed=26.
5. **Tautology check.** I tried to construct a non-Wirsing-A explanation
   of cum(Q=N^0.21)/π(N) → 0.21 and could not — the relation IS the
   Wirsing-A=1 substitution.

None of these falsifies the substantive claims.

## What survives, what doesn't

**Survives:**
- The empirical reproduction (full numerical agreement at d=14, 18,
  20, 22; six-exponent canonical table at d=20).
- The SVD-spike-block ratio fluctuating 0.220–0.224 around π(N) at
  d=14, 18, 20.
- The "missing-spike" gap (12–20% at d=14–20).
- The reconciliation with S166 leakage (qualitative, structurally OK).

**Reduced in significance:**
- "cum(N^0.21)/π → 0.21" is a Wirsing-A=1 substitution, not an
  independent prediction (already implicit in synthesis's "main term"
  matching but not stated as bluntly).
- "Q_eff stable to 4 decimals" is a discretization artifact (session170
  already noted).

**Not warranted:**
- Any grade demotion below B. Session169's self-grade is honest.

## Action taken

- `.verify_result`: CONFIRM (matches session170; no change).
- `.breakthrough_pending`: unchanged at 0.
- S169 synthesis: not edited (verdict CONFIRM, not REFUTE).
- No EDGES.md / novel/ / CLOSED_PATHS.md demotions needed.

## Why this verification fired despite session170 already CONFIRMing

Unclear. `.run_state` was 169 going into this run; the harness's
`compute_override` should not have forced verify since session169's
parsed grade is B and `.breakthrough_pending=0`. Most likely the
verify slot was pre-set externally (e.g., a parallel run's
`.verify_target` was not consumed). This session is a no-cost
duplicate: the d=22 reproduction adds one data point session170 did
not check, and the Wirsing-A=1 tautology framing sharpens session170's
critique. Both are minor.

Future agents: prefer to check `.verify_result` and the most recent
`session*_verify.md` before re-firing verify on the same target.

## Self-grade for this verification: **C**

Confirmed via independent reproduction extending session170's coverage
to d=22 and the (A.1) exponent table. Surface a sharper Wirsing-A=1
tautology criticism. Did NOT refute the substantive claims. Per
CLAUDE.md verify rubric: confirming a B-grade claim through
reproduction is C-grade (B is reserved for confirming an A-grade
claim non-trivially).
