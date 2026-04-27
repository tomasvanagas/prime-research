# D2.a — Persistent homology of W=210 W-tricked normalised prime gaps

**Status:** PRE-REGISTERED (in-progress).
**Edge candidate:** test of E2.17 (PH deficit) reduction to E2.13 (Gowers
HL singular series, W-trick uniformity).
**Cross-domain technique:** persistent homology (Carlsson 2009 BAMS 46;
Bauer 2021 ripser) + W-trick (Green-Tao 2008 *Annals* 167).

## Question

NOVELTY_CHALLENGES.md §D2.a: the W=210 W-trick restores Gowers `U^2`
uniformity of `chi_P` to within 0.1% (E2.13). Does it also restore
PH-Poisson behaviour of the gap sequence (E2.17 result, S96)?

If yes — E2.17 is an alternate observable for the *same physics*
(HL singular series at W=210). If no — PH detects a structural
obstruction beyond singular-series cancellation.

## Protocol

Same Takens-embed + ripser PH pipeline as `persistent_homology_chi_p.py`
(S96), but the input is filtered to a single residue class mod W=210:

1. Sieve primes up to `N_max = 10^7` (664,579 primes).
2. For each `b ∈ {1, 11, 13}` (gcd(b, 210) = 1), filter primes to
   `q_n ≡ b (mod 210)`.
3. Take `M = 2000` consecutive primes in the residue class starting at
   the first `q ≥ 10^6`. Compute gaps `g_n = q_{n+1} - q_n`.
4. Cramér-normalise per residue:
   `x_n = g_n / (φ(W) log q_n)` with `φ(210) = 48`.
   Under Cramér + Dirichlet equidistribution, `x_n → Exp(1)`.
5. Takens-embed at `d = 3`, `τ = 1`. Run ripser with `thresh = 4.0`.
6. Compute `T0`, `T1`, `L0`, `L1`, `N1` (same as S96).
7. Baselines (K=20 each, fresh per residue):
   - `B1` IID `Exp(1)` of length `M`.
   - `B2` random permutation of the empirical W-tricked window.
8. Pool z-scores across the three residues; report range and mean.

## Falsifiers (pre-stated, before running the experiment)

The S96 unconditioned anchor at `M = 2000, d = 3, x ≈ 10^6`:

| Statistic | PRIMES | z(P;B1) | z(P;B2) |
|-----------|--------|---------|---------|
| T0        | 243.34 | **−10.31** | **−8.70** |
| T1        |  37.24 |  −4.20 | **−11.99** |

### F_a — W-trick erases the PH deviation

Pooled z (mean across `b ∈ {1, 11, 13}`) lies within `±3` for
**both** `T0` and `T1` against **both** B1 and B2.

- **Outcome:** PH deficit IS the HL singular-series obstruction.
- **Edge action:** annotate E2.17 as "reduces to E2.13 under W-trick";
  add closure row to CLOSED_PATHS marked DUPLICATE-PLUS of E2.13.
- **Grade:** C (verification of an existing structural reading).

### F_b — W-trick partially erases the deviation

At least one of `|z_pool(T0; B1)|`, `|z_pool(T0; B2)|`,
`|z_pool(T1; B1)|`, `|z_pool(T1; B2)|` is reduced by ≥ 50% relative
to the S96 anchor, but at least one `|z|` remains > 3.

- **Outcome:** PH detects HL singular-series content + a structural
  residual. The residual is the new content; localise via per-residue
  `T0` vs `T1` reduction (e.g., if `T1` z survives but `T0` does not,
  the loop deficit is residual; if `T0` survives but `T1` does not,
  the cluster deficit is residual).
- **Edge action:** **refine E2.17** with quantitative W-trick reduction
  factors per summary; flag the surviving observable as a candidate
  *cross-W invariant*.
- **Grade:** B (substantive refinement of E2.17 with new quantitative
  content).

### F_c — W-trick does NOT erase the deviation

All four S96 z-scores reproduced within `25%` relative on the pooled
W-tricked sequence — `z_pool(T0; B1) ∈ [−12.9, −7.7]`,
`z_pool(T0; B2) ∈ [−10.9, −6.5]`, `z_pool(T1; B1) ∈ [−5.3, −3.2]`,
`z_pool(T1; B2) ∈ [−15.0, −9.0]`.

- **Outcome:** PH deficit is NOT the HL singular-series obstruction.
  Genuine new structural content beyond E2.13 — a separate metric-
  topological observable on the gap sequence.
- **Edge action:** elevate E2.17 to a *negative-shape edge* — it
  detects a structure NOT killed by W-trick, separating it from the
  E2.13/E2.14/E2.15/E2.16 chain.
- **Grade:** B (negative-shape edge candidate).

(Outcomes F_a and F_c are mutually exclusive; F_b covers the middle.)

## Files

- `persistent_homology_w_trick.py` — driver.
- `w_trick_main.json` — main run results (M=2000, d=3, x≈10⁶,
  residues 1/11/13 mod 210).
- (post-run) `summary_table.md` — pooled z-score table vs S96.

## Edges cited / composed

- **Composes** E2.13 (Gowers `U^2` W-trick uniformity) with
  E2.17 (PH deficit on gap sequence).
- **Tests reduction** of E2.17 to E2.13 (F_a), partial reduction
  (F_b), or independence (F_c).
- **Cross-references** E2.14 (Anderson localisation × W-trick) for
  the "W-trick erases the HL fingerprint" pattern.
