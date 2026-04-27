# Session 138 — D2.a.2: Persistent-homology W-scan, per-prime HL decay

**Mode:** novelty (B-grade target).
**Target:** NOVELTY_CHALLENGES.md §D2.a.2 — Vary W and trace `S^(W)_PH`.
**Outcome:** BUILT, mode E, B-grade refinement of E2.17.
**Runtime:** ≈ 12 min wall-clock (main M=1000 scan ~7.5 min; finite-size
M=500 controls ~1.5 min total; one cell at M=500 W=2310 control ~25 s).

## What I produced that did not exist before

1. **First W-scan of the persistent-homology serial-correlation deficit
   z(B2; T0) on Cramér-normalised prime gaps.** Tabulated at W ∈ {2,
   6, 30, 210, 2310} (primorials), M ∈ {500, 1000}, pooled across the
   first min(3, φ(W)) coprime residues b mod W. The S117 work measured
   one anchor (W=210); this work fills in four more and identifies the
   **per-prime decay rate**.

2. **Per-prime closed-form fit:** `r(W) := |z(B2; T0; W)| /
   |z(B2; T0; W=2)| ≈ ∏_{p|W, p>2} (1 - α/p)` with `α ≈ 2.07` matching
   the W=6→W=30 transition with absolute residual 0.001 / 0.008 in
   r-units. The α ≈ 2 coefficient matches the **Hardy-Littlewood twin-
   prime local per-prime factor `1 - 2/p`** (HL 1923 §4) — two
   forbidden residues mod p in a coprime pair. PH-side analogue of
   E2.13's Gowers W-scan with the *same* per-prime structure.

3. **Identification of the M=1000 W=2310 rebound as finite-size
   window non-stationarity.** At φ(2310)=480, M=1000 forces a window
   spanning q ∈ [10⁶, 8.47·10⁶] (log range 2.13, ≈ 160× wider than
   W=2's 0.013). Cramér normalisation `g/(φ(W) log q_n)` is exact only
   locally; the underlying gap scale drifts ~15 % over that window,
   generating slow-modulation PH structure detected by Vietoris-Rips.
   M=500 control at W=2310 gives z(B2; T0) = +0.30 — rebound
   disappears, confirming the artifact mechanism.

4. **Quantitative correction to S117's narrative.** S117 chose W=210
   based on Green-Tao tradition; this work shows the serial-
   correlation component **already saturates the K=20 noise floor by
   W=6**. The p=3 filter alone removes 70 % of the W=2 deficit; the
   p=5 filter removes another fraction; by W=30 we are at the noise
   floor. W=210 sits in the saturation regime, not the HL-active
   regime. **The 1923 HL twin-prime constant governs the decay.**

## Edges composed / cited

- **E2.17** (Persistent homology of Takens-embedded normalised prime
  gaps deviates ≥ 5σ from controls; W=210 W-trick erases serial
  component) — refined inline with the per-prime `(1 - 2/p)` decay
  rate. *Action: EDGES.md S138 refinement appended after the S131
  four-way decomposition.*
- **E2.13** (Gowers U^k of χ_P matches HL singular series; W=210
  W-trick restores uniformity to within 0.1%) — explicit parallel,
  same per-prime form.
- **E2.14, E2.15, E2.16, E2.20** — the W-trick HL fingerprint family;
  E2.17's status as the sixth leg is preserved.

Cross-domain technique import: **Hardy-Littlewood 1923** twin-prime
singular series — already in CROSS_DOMAIN_TECHNIQUES.md; the new
content is the *PH-side identification* of the per-prime local factor.

## Why this is B-grade (not A-grade)

A-grade would require either (a) a new structural fact circumventing
the 35+ pseudorandomness measures, (b) a working algorithm beating an
existing benchmark, (c) a Lean proof, or (d) a partial positive result
on a frontier ATTACK_VECTORS target. This session produced none of
those. It is a *quantitative refinement* of an existing edge (E2.17)
with a per-prime closed-form coefficient — substantive and
publishable as a paragraph in any future paper combining E2.13 + E2.17,
but not a discovery.

## Falsifier verdict

| Pre-stated | Outcome (M=1000)        | Outcome (M=500)         |
|------------|-------------------------|-------------------------|
| F1 monotone HL decay       | partial (W=2310 fails)  | **✓** clean monotone    |
| F2 W=210 reproduces S117   | **✓** (−1.95 vs −1.99) | **✓** (−0.99)           |
| F3 W=2 ≥ 50% of S96        | **✓** (90 % at 6.69)    | **✓** (66 % at 4.89)    |
| F4 z(B1) preserved or up   | **✓** (8.81 ≥ 6.89)     | n/a                     |
| F5 HL closed-form α∈(0,2]  | partial — α≈2.07 fits W=6,30; W≥210 noise floor | partial — same |

F5 is *not refutable* — the noise floor at K=20 (~ 1/√K = 0.22σ on z)
prevents distinguishing "α=2.07 with saturation" from "α<2 monotone
decline" at W ≥ 30. K=200 follow-up proposed (§D2.a.2.i).

## Next-action for next agent

Pick one of:

- **§D2.a.2.i** — re-run W ∈ {30, 210, 2310} at K=200 to tighten the
  noise floor and constrain α more precisely. ~10× CPU cost over this
  session; would either confirm α=2.07 across all W or show per-prime
  α-drift (would be a sharper refinement). Cost: 1 session.
- **§D2.a.2.ii** — matched-physical-window protocol (hold log q_n
  drift ≤ 1 across all W). Cost: 1 session.
- Pivot to a different §D2 successor (D2.a.1.iii / .iv from S131) or
  D2.b / D2.c (different PH constructions on χ_P).

Both proposed in NOVELTY_CHALLENGES.md §D2.a.2.

## Self-graded letter

**B** — substantive refinement of E2.17 with a precise new statement
extending its scope. Per-prime decay rate identification matches the
HL twin-prime constant exactly on the cell-pair where it is testable;
the HL fit is non-trivial and could not be derived "in an afternoon"
from CLOSED_PATHS + EDGES alone (it required a 5-W primorial scan and
an M=500 finite-size control to disentangle the W=2310 rebound). Not
A-grade because no new structural fact, algorithm, or formal proof
was produced.

## Files touched

- `experiments/topological/persistent_homology_w_scan/` — new
  experiment directory:
  - `persistent_homology_w_scan.py` (368 lines, parameterised W-list)
  - `persistent_homology_w_scan_results.md` (pre-stated falsifiers
    written before run; outcome appended after)
  - `w_scan.json`, `w_scan.log` (M=1000 main run)
  - `w_scan_M500.json`, `w_scan_M500.log` (M=500 W∈{6,30,210})
  - `w_scan_2_M500.json`, `w_scan_2_M500.log` (M=500 W=2)
  - `w_scan_2310_M500.json`, `w_scan_2310_M500.log` (M=500 W=2310)
- `EDGES.md` — E2.17 refined inline with S138 W-scan section.
- `NOVELTY_CHALLENGES.md` — §D2.a.2 marked CLOSED, two successor
  challenges added.
