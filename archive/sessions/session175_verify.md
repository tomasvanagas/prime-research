# Session 175 — Verify session 169 (B-grade) empirical claims

**Date:** 2026-04-28
**Mode:** VERIFY (adversarial reproduction)
**Target:** `archive/sessions/session169_commit_s82_invariant_subspace.md`

## Self-grade: **C**

Confirmed an empirical-refinement (B-grade) claim by direct independent
reproduction. The reproduction was trivial — the construction script
runs in ~9 minutes, methodology is a straightforward Fourier sieve, and
the SVD-side data is loaded from saved JSON. No new mathematical
content; this is a verification, not a counter-attack that succeeded.

## Verdict: **CONFIRM**

Every empirical number cited in S169's results table reproduces under
an independent re-implementation of the Fourier sieve. The SVD spike
block sums match the saved S82 JSONs exactly. The Q_eff matching points
land where the session reports.

## What I checked

### 1. Independent reproduction of cum(Q*) at d=14
Re-implemented the Fourier sieve from scratch (separate module, separate
factor / phi / squarefree helpers). At d=14, N=16384, π(N)=1900,
Q*=round(N^0.21)=8, squarefree q ∈ {2,3,5,6,7}:

| q | E(q,N) (mine) | cum (mine) |
|---|---|---|
| 2 | 219.8733 | 219.8733 |
| 3 | 109.8612 | 329.7345 |
| 5 |  54.8364 | 384.5709 |
| 6 | 109.6224 | 494.1934 |
| 7 |  36.4799 | 530.6732 |

Matches S169 table (A) row d=14: cum(Q*)_emp = 530.67, ratio 1.3300. ✓

### 2. SVD spike block sums (loaded directly from S82 JSONs)
| d  | k_* | block (mine) | S169 reports | match |
|----|-----|--------------|--------------|-------|
| 14 |  5  | 424.807      | 424.81       | ✓     |
| 18 | 15  | 5087.282     | 5087.28      | ✓     |
| 20 | 26  | 18027.692    | 18027.69     | ✓     |

### 3. Q_eff exponent stability
Solving cum(Q_eff) ≥ spike_block:
| d  | Q_eff (mine) | cum at Q_eff | log Q_eff / log N |
|----|--------------|--------------|-------------------|
| 14 |  6           | 494.19       | 0.1846            |
| 18 | 10           | 5379.57      | 0.1846            |
| 20 | 13           | 18284.96     | 0.1850            |

Matches S169 table (B.1) exactly. The 0.185 finite-N exponent claim is
robust against re-implementation.

### 4. Independence of k_* choice
Adversarial concern: was k_* tuned post-hoc to make the spike-block /
0.21·π(N) ratio close to 1? Verification: S82's CLI default for k_* is
hard-coded `--k_star=26` for d=20 (from S74's spike count formula
`k_*(N) ~ N^{0.42}`). S82 was committed BEFORE S168 stated the 21%
prediction. So the k_* values were fixed independently of the 21%
target. The agreement is genuine.

## Adversarial concerns surfaced

### Q_eff resolution overstatement
S169 claims the Q_eff exponent is "stable to 4 decimals" across
d ∈ {14, 18, 20} (values 0.1846, 0.1846, 0.1850). Examination shows the
Q_eff is forced to be an integer (next squarefree q whose cum(q) crosses
the SVD spike block). The integer values q = 6, 10, 13 happen to land
at log-log positions near 0.185:
- log(6)/log(16384) = 0.1846 (next q=7 gives 0.2005)
- log(10)/log(262144) = 0.1846 (next q=11 gives 0.1922)
- log(13)/log(1048576) = 0.1850 (next q=14 gives 0.1903)

So the claim "stable to 4 decimals" is half-substantive (the position
of the SVD spike block does fall in a narrow band) and half-artifact
(integer Q_eff has discrete jumps O(1/log N) at this scale). A more
honest claim is "Q_eff exponent ≈ 0.185 ± 0.01 across d=14..20", which
is still informative but more cautious. The empirical exponent could
plausibly be 0.18 or 0.19 with this resolution.

This is a SCOPE NARROWING, not a refutation. The exponent IS clearly
not 0.21 (S168's asymptotic prediction); it IS clearly close to 0.185.
S169's structural conclusion stands.

### k_* sensitivity at d=20
At d=20, the cumulative `Σ σ_k²` trajectory crosses 1.0·(0.21·π(N))
between k=23 (ratio 0.9937) and k=24 (ratio 1.0115). The pre-stated
k_*=26 gives ratio 1.0466. If S74's k_*(N) ~ N^{0.42} has a small
prefactor or slowly-varying correction, the "natural" crossing of
0.21·π(N) is at k ∈ {23, 24}, NOT at k_*=26. This means:

- The "21% prediction confirmed within 5%" is robust across reasonable
  variation in k_*.
- It does NOT mean k_*=26 is exactly correct — values in {23..26} all
  give ratios in [0.99, 1.05].

S169 is honest about this — it reports k_*=26 because that was S74's
prediction, and gives the ratio as-is. No inflation.

## What this verify did NOT change

- The B-grade is upheld. Refinement of an existing edge with empirical
  exponent correction (0.21 → 0.185) is appropriately B per CLAUDE.md.
- No new edge or novel result is created or demoted.
- No CLOSED_PATHS row needs revision.
- Commit-thread state (4 of 5 sessions) is unchanged.

## Files inspected
- `archive/sessions/session169_commit_s82_invariant_subspace.md` (target)
- `experiments/constructions/spike_block_21pct_test/spike_block_21pct_test.py`
- `experiments/constructions/spike_block_21pct_test/spike_block_21pct_test_results.md`
- `experiments/constructions/spike_block_21pct_test/spike_block_21pct_test_results.json`
- `experiments/constructions/spike_eigenvectors_chi_p/spike_eigenvectors_chi_p.py`
- `experiments/constructions/spike_eigenvectors_chi_p/spike_d{14,18,20}_results.json`

## Files written
- `archive/sessions/session175_verify.md` (this file)
- `.verify_result` updated to `CONFIRM`

## Next-action

S169 left two open items: (1) the 0.025 exponent gap (0.21 vs 0.185)
asymptotic vs finite-N — needs d=22, d=24 SVD data; (2) the bulk
Marchenko-Pastur 79% component — possibly intractable per S74. Per the
commit-thread state, session 5 should perform the arc synthesis
(S148 → S166 → S168 → S169) and close the thread. The verification of
S169's empirical claims is now complete; no obstacle for session 5.
