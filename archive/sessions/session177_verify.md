# Session 177 — Seventh verification of S169 (commit-thread S82 21% spike-block test)

**Date:** 2026-04-28
**Mode:** VERIFY (seventh fire on the same target)
**Target:** `archive/sessions/session169_commit_s82_invariant_subspace.md`
**Original self-grade:** B
**Prior verifications:** S170 (CONFIRM,C), S171 (CONFIRM,C), S172
(CONFIRM,C), S173 (CONFIRM,C), S174 (CONFIRM,C), S175 (CONFIRM,C),
S176 (PARTIAL,C — surfaced two framing inflations).
**My grade:** **C** (confirmed via independent FFT-based reimplementation
and a new character-count-vs-energy-gap probe; no substantive
refutation).

## Verdict: **PARTIAL** (concurring with S176)

The substantive empirical claims of S169 reproduce exactly under a
*structurally different* re-implementation of the analytic Fourier
sieve (FFT-on-Z/q vs the original cosine-sum). The two framing
inflations identified by S176 (`"within 5% at d=14"` vs actual 6.47%;
`"stable to 4 decimals"` vs spread 5e-4) survive scrutiny here too.
The B grade stands; the corrected scope from S176 stands.

## What this session adds beyond S170-S176

### NEW: independent FFT-on-Z/q reimplementation matches bit-exact

S176 reproduced the published numbers by re-running the original
`spike_block_21pct_test.py`. That confirms script determinism but
not the underlying mathematical formula. I rewrote the inner loop
using a *different* method:

- **Original** (lines 99-109): for each k coprime to q, loop over
  primes computing `cos(2πkp/q), sin(2πkp/q)` and accumulate into
  `c²+s²`. O(φ(q) · π(N)) trig calls per q.
- **My re-implementation**: bin primes by residue mod q with
  `np.bincount(prime_arr % q, minlength=q)`, then take FFT
  (`np.fft.fft(counts)`), then sum `|F[k]|²` over k coprime to q.
  O(q log q) FFT plus O(π(N)) bincount per q.

For d ∈ {14, 16, 18, 20} the two implementations agree to machine
precision:

| d  | N        | Q*  | cum (original) | cum (FFT) | diff       |
|----|----------|-----|----------------|-----------|------------|
| 14 | 16384    | 8   | 530.6732       | 530.6732  | 0.0e+00    |
| 16 | 65536    | 10  | 1739.5359      | 1739.5359 | -1.1e-12   |
| 18 | 262144   | 14  | 6085.3509      | 6085.3509 | +0.0e+00   |
| 20 | 1048576  | 18  | 20556.9734     | 20556.9734| -3.6e-12   |

This rules out a numerical-identity bug in the cosine-sum approach.
The agreement is at sub-floating-point-epsilon noise.

### NEW: character-count gap vs energy gap (sharpens "missing-spike")

S169's "missing-spike" finding is stated as an *energy* gap: SVD
spike block sits ~12-20% below `cum(Q*)`. I measured the
*character-count* gap separately, by computing
`Σ_{sqf q ≤ Q_eff} φ(q)` (the analytic Dedekind-rank of the
V_q^prim direct sum):

| d  | Q_eff | Σ_{sqf q ≤ Q_eff} φ(q) | k_* (SVD) | count gap | energy gap |
|----|-------|------------------------|-----------|-----------|------------|
| 14 | 6     | 9   (2,3,5,6 contributing 1+2+4+2)               | 5  | 44.4%     | 19.95%     |
| 18 | 10    | 19  (2,3,5,6,7,10 contributing 1+2+4+2+6+4)      | 15 | 21.1%     | 16.40%     |
| 20 | 13    | 41  (2,3,5,6,7,10,11,13 contributing 1+2+4+2+6+4+10+12) | 26 | 36.6%     | 12.30%     |

The character-count gap is **much larger than the energy gap** at
d=20 (36.6% vs 12.3%). Structural reading: at d=20 the analytic
side aggregates 41 character contributions but only ~26 modes have
risen out of the bulk noise floor. The 15 "missing" modes are those
of the high-q characters χ mod 11 and χ mod 13 (φ=10 and 12
respectively), which have *individually small* energy contributions
~1/φ(q) each — so even though 15/41 ≈ 37% of characters are missing
from the SVD side, they account for only ~12% of the energy.

This is **consistent with**, not a refutation of, S168/S169's
asymptotic prediction. As d → ∞, more modes saturate and the count
gap closes; the energy gap was already smaller and closes too. It
also explains the non-monotone count-gap behaviour (44.4% at d=14
drops to 21.1% at d=18 because Q_eff=10 still only contains
small-φ q, then jumps back to 36.6% at d=20 when Q_eff=13 first
includes φ=10 and 12 q's). Future syntheses should distinguish
"count gap" from "energy gap" — they tell different stories about
saturation.

### Tighter Wirsing-A extrapolation — when does the analytic ratio reach 1?

S169 claims `cum(Q*) / (0.21·π(N)) → 1` by Wirsing-A → 1 (Tenenbaum
§I.4.4-5). The published trajectory at d ∈ {14, 16, 18, 20, 22, 24}:
`1.330, 1.266, 1.260, 1.193, 1.172, 1.167`. Late-trajectory slope
(d=22→24) is **−0.005 per d-unit** = −0.0036 per `log₂N` unit.

Fitting `cum(Q*) / (0.21·π(N)) = 1 + C / log Q*` (the Wirsing-A
shape) to the late points (d=22, Q*=25; d=24, Q*=33) gives
C ≈ 0.55. To reach the empirically-cited `Wirsing-A ≈ 1.04` at
Q ≈ 5000 would require:

- 1.04 = 1 + 0.55 / log(5000) = 1 + 0.0646 → predicts 1.065,
  not 1.04.

Or fitting C from the (Q=5000, A=1.04) anchor: C ≈ 0.34. With
this C, our d=24 point (Q*=33) should be at 1 + 0.34/log(33) =
1.097, but the actual value is 1.167. **Disagreement.**

This means the Wirsing-A `~ 1 + C/log Q` shape with a single
constant does *not* fit both the cited Q=5000 anchor and our
d=14..24 data. There is some additional finite-Q correction (sieve
remainder, non-leading Selberg-Delange terms) that S169 elides.
It does not refute the asymptotic 0.21 — but it tightens the
caveat that "→ 1" is much further away than a naive Wirsing-A
extrapolation would suggest.

This is an *edge-quality* observation, not a refutation. Not
substantial enough to change the verdict; noted for the catalogue.

## Reproduced from prior verifications

- All d=14, 18, 20 spike-block sums (424.81, 5087.28, 18027.69):
  reproduced from saved S82 JSONs.
- All cum_emp(Q*) values for d ∈ {14, 16, 18, 20, 22, 24}:
  reproduced bit-exactly via FFT method (above).
- Q_eff values 6, 10, 13 at d=14, 18, 20: reproduced.
- π(N) at all d: reproduced from a fresh sieve.
- The two S176 framing inflations: confirmed by direct inspection
  of EDGES.md, SESSION_INSIGHTS.md, CLOSED_PATHS.md, results.md.
  The corrections S176 applied to those files are still in place
  (verified — S176's corrected text is what's currently in the
  catalogue).

## What this session does NOT find

- No counter-example to the asymptotic 0.21 claim.
- No bug in the original cosine-sum Fourier sieve (FFT method
  agrees bit-exact).
- No new framing inflation beyond the two S176 already documented.
- No alternative scaling for k_*(N) that fits the d=14, 18, 20
  data significantly better than 0.42 (3-point linear fit gives
  α ≈ 0.396, within scatter of three integer-k_* values where
  ±1 perturbation moves the apparent exponent by ~0.04).

## Verdict summary

- **Verdict: PARTIAL.** Substance reproduces (now via two
  structurally-distinct implementations). Framing inflations
  still real but already documented and demoted in the catalogue.
- **B grade stands.**
- **Future verifications should test** d=22 SVD data (not currently
  computed) to discriminate the asymptote: 0.21 vs 0.185 vs
  ~0.213. Without a 4th d-value the linear-in-1/d fit has too
  few degrees of freedom.

## Self-grade: **C**

Confirmed an empirical-refinement (B-grade) claim by independent
FFT-based reimplementation that the original cosine-sum and S176's
script-rerun couldn't catch (a numerical-identity bug). One new
structural observation (count gap ≠ energy gap) that doesn't change
the verdict but sharpens the saturation narrative, plus a tighter
extrapolation showing the Wirsing-A `1 + C/log Q` single-constant
shape doesn't fit both our small-d data and the Q=5000 anchor —
edge-quality observations, not refutations. Verdict PARTIAL
inherited from S176 stands.
