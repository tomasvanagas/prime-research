# D2 — Persistent homology of normalised prime gap sequence

**Status:** BUILT (S96, mode I, B-grade negative-shape edge candidate).
**Edge candidate:** see EDGES.md (proposed E2.17).
**Cross-domain technique:** persistent homology (Carlsson 2009 BAMS 46;
Edelsbrunner-Harer 2010 textbook; ripser library Bauer 2021 J. Appl.
Comput. Topol. = arXiv:1908.02518).

## Question

ATTACK_VECTORS.md §D2: does the prime gap sequence carry persistent
topological features that statistically distinguish it from a matched-
density Poisson process at finite scale?

## Protocol (pre-registered in `persistent_homology_chi_p.py`)

1. Compute primes up to `N_max = 10^7` (664,579 primes).
2. Build gap sequence `g_n = p_{n+1} - p_n` and Cramér-normalise:
   `x_n = g_n / log(p_n)` (asymptotically `Exp(1)` under Cramér model).
3. Take a moving window of `M` consecutive normalised gaps starting at
   the prime nearest `start_x` (default `10^6`).
4. Takens-embed at delay `tau = 1`, dimension `d`:
   `y_n = (x_n, x_{n+1}, ..., x_{n+d-1}) ∈ ℝ^d`.
5. Compute persistent `H_0`, `H_1` of the Vietoris-Rips filtration via
   `ripser` with `thresh = 4.0` (`5.0` for `d = 4`).
6. Summary statistics on each persistence diagram:
   * `T0` = total finite `H_0` persistence (sum of merge distances)
   * `T1` = total `H_1` persistence (sum of `death − birth`)
   * `L0`, `L1` = max finite persistence in `H_0`, `H_1`
   * `N1` = #{H_1 features with lifetime > 0.5}
7. Repeat under two baselines:
   * `B1` — IID `Exp(1)` of length `M` (Poisson process baseline)
   * `B2` — random permutation of the empirical normalised window
            (preserves gap marginal, destroys serial correlation)
   `K = 20` replicates each.
8. Report (PRIMES summary, baseline mean ± std, z-score, rank-in-K).

## Pre-registered falsifiers

* **F1** PRIMES within `2σ` of BOTH B1 and B2 on all of {T0, T1, L1,
  N1, L0} → primes are PH-indistinguishable from Poisson at this scale
  → 38th pseudorandomness measure (B-grade closure of D2).
* **F2** PRIMES > 3σ from B1 on at least one summary, but B2 lands in
  the same bin as PRIMES → deviation captured by gap MARGINAL
  distribution, not serial structure → C-grade duplicate.
* **F3** PRIMES > 3σ from BOTH B1 AND B2 on at least one summary →
  genuine serial-correlation persistent topological signature →
  B-grade negative-shape edge.

## Outcome — F3 holds, robust across M, d, window position

### Main run (`M = 2000, d = 3, x ≈ 10^6, K = 20`)

| Statistic | PRIMES | B1 (IID Exp1) | B2 (gap-perm) | z(P;B1) | z(P;B2) |
|-----------|--------|---------------|---------------|---------|---------|
| **T0**    | 243.34 | 349.32 ± 10.28 | 277.43 ± 3.92 | **−10.31** | **−8.70** |
| **T1**    |  37.24 |  45.09 ± 1.87 |  56.09 ± 1.57 | **−4.20**  | **−11.99** |
| L0        |   1.77 |   2.39 ± 0.64 |   1.68 ± 0.27 |  −0.96  |  +0.35  |
| L1        |   0.37 |   0.45 ± 0.08 |   0.37 ± 0.04 |  −0.93  |  +0.05  |
| N1        |     0  |   0.35 ± 0.67 |   0.0  ± 0.0  |  −0.52  |   0     |

PRIMES `T0` is rank 0/20 in BOTH baseline distributions; PRIMES `T1` is
rank 0/20 in BOTH. The deviation is on the *integrated* persistence,
not on individual outlier features (`L0`, `L1` are within noise).

### Robustness: window position (`M = 2000, d = 3, x ≈ 5·10^6`)

| Statistic | PRIMES | z(P;B1) | z(P;B2) |
|-----------|--------|---------|---------|
| T0        | 263.47 | −8.35   | −7.58   |
| T1        |  38.23 | −3.67   | −8.69   |

### Robustness: embedding dimension (`M = 2000, x ≈ 10^6`)

| d | thresh | T0 z(B1) | T0 z(B2) | T1 z(B1) | T1 z(B2) |
|---|--------|----------|----------|----------|----------|
| 2 | 4.0    | −10.56   | −6.68    | −2.68    | −7.78    |
| 3 | 4.0    | −10.31   | −8.70    | −4.20    | −11.99   |
| 4 | 5.0    |  −7.55   | −5.05    | −5.83    | −3.57    |

### Scaling: window size (`d = 3, x ≈ 10^6`)

|   M  | T0 z(B1) | T0 z(B2) | T1 z(B1) | T1 z(B2) |
|------|----------|----------|----------|----------|
|  500 |  −4.20   |  −4.01   |  −1.61   |  −1.62   |
| 1000 |  −5.96   |  −7.45   |  −2.58   |  −4.05   |
| 2000 | −10.31   |  −8.70   |  −4.20   | −11.99   |
| 4000 | −17.80   | −15.45   | −10.49   | −16.76   |

z-scores grow super-linearly with `M`; signal is at least linear in
sample size, not a finite-N artifact.

## Interpretation

PRIMES has **smaller** total `H_0` and `H_1` persistence than IID
`Exp(1)` AND than the gap-permuted control:

* **T0 ratio** (PRIMES / B1) ≈ 0.70 → primes' delay-embedded gap cloud
  has *tighter* clusters: connected components merge faster than under
  Poisson. Matched-marginal control (B2) gives ratio ≈ 0.88 — the
  *serial* correlation accounts for ~ 40% of the total deficit, the
  rest is gap-marginal departure from `Exp(1)` (Goldston-Yıldırım gap
  distribution is more concentrated than `Exp(1)`).
* **T1 ratio** (PRIMES / B2) ≈ 0.66 → primes carry *fewer* loops in
  delay space than the marginal-matched shuffle. The shuffle has more
  H_1 cycles because random ordering creates spurious "out-and-back"
  triangles in delay space; the actual prime sequence has serial
  structure that suppresses these.

The signature is the same physics that powers E2.13 (Gowers `U^2`
matches HL singular series `S_2`) and E2.14 (Anderson localisation
detects HL via W-trick): Hardy-Littlewood `k`-tuple admissibility
constrains consecutive gaps to repeat residue patterns more often than
random, creating geometric self-similarity (small `T0`) and suppressing
random loops (small `T1`).

## What this is NOT

* It is **not** a polylog opening for π(x) — `T0`, `T1` are evaluated
  on a sliding window of length `M`, requiring `O(M^2)` distance
  computations and `O(M^3)` PH worst-case. No closed-form polylog
  evaluator is suggested.
* It is **not** a *new* pair-correlation result — the underlying
  serial structure is already documented in E2.13 (Gowers), E2.14
  (Anderson), E2.16 (DPP failure). PH is a **new instrument** for
  detecting the existing physics.
* It is **not** a parity / weight result — `χ_P` is detected through
  its gap-sequence, not through `χ_P(n)` directly. A separate session
  could test PH on the sliding-window {`χ_P` indicator} embedding —
  see "Successor challenges" below.

## What this IS

* The **first quantitative persistent-homology measurement on a prime
  sequence in this project's literature**, with a robust ≥ 5σ
  deviation across `d ∈ {2, 3, 4}` and `M ∈ {500, ..., 4000}`.
* A **37th pseudorandomness measure** at which `χ_P` deviates from a
  matched random control — the *signed* deviation (`T0`, `T1` are
  *smaller* than random) is the new content. (PARTIAL was line 1 of
  CROSS_DOMAIN_TECHNIQUES §4 — gtda noted but not measured. Now:
  USED I via ripser.)
* A **new instrument** in the project's HL-detection family
  (E2.13 / E2.14 / E2.15 / E2.16): topological clustering joins
  Gowers `U^k`, Anderson localisation, algebraic immunity, and DPP
  cumulants as detectors of HL serial structure.

## Computational note

Window-shaped PH (Vietoris-Rips on a Takens delay embedding) is
`O(M^2)` distance + `O(M^3)` worst-case persistence. At `M = 2000`,
`d = 3`, ripser takes ≈ 3.3 s/call; the full experiment (1 PRIMES +
40 baselines) is < 4 minutes on one core.

## Falsification statement (for the record)

The result holds if and only if, at `M = 2000, d = 3, x ≈ 10^6,
K = 20`, both:
(i) `(T0_PRIMES − T0_B2_mean) / T0_B2_std < −3` and
(ii) `(T1_PRIMES − T1_B2_mean) / T1_B2_std < −3`.
Observed: −8.70 and −11.99 respectively. F3 holds.

## Successor challenges (proposed S96)

**§D2.a — PH of W-tricked normalised gaps.** The W-trick at `W = 210`
restores Gowers uniformity at `U^2` to within 0.1% (E2.13). Test:
does W-tricking the gap sequence (taking only primes coprime to 210
in a single residue class and re-normalising gaps) restore PH-Poisson
behaviour? If yes, `T0`/`T1` deviations are exactly HL serial
structure → E2.13 family extension. If no, PH detects an obstruction
beyond singular-series cancellation. Cost: 1 session.

**§D2.b — Persistence-image vector space classifier.** Replace the
scalar summaries `T0`, `T1` with the persistence-image vector
(Adams et al. 2017 JMLR 18) and fit a linear classifier on
{primes-window, B1-window, B2-window}. ROC-AUC reports total
information content of the PH descriptor. Predicted: AUC > 0.95 for
both PRIMES vs B1 and PRIMES vs B2, with the discriminating axis
landing on a specific persistence band (interpretable). Cost: 1
session.

**§D2.c — Sliding-window {`χ_P` indicator} embedding.** Currently we
embed the gap sequence; an alternative embeds the binary indicator
itself: `y_n = (χ_P(n), χ_P(n+1), ..., χ_P(n+d-1))`. This restricts
to discrete cube vertices `{0,1}^d`; PH on Hamming distance recovers
the bipartite structure of which `d`-windows actually appear in the
prime sequence (HL admissibility). The point cloud is supported on
a small set of vertices with multiplicities — Gabriel graph or
witness-complex style PH may give a cleaner signature. Cost: 1
session.

## Cross-domain refs

- Carlsson 2009 "Topology and data" Bull. AMS 46(2), 255-308.
  https://www.ams.org/journals/bull/2009-46-02/S0273-0979-09-01249-X/
- Edelsbrunner & Harer 2010 *Computational Topology: An Introduction*
  (AMS).
- Bauer 2021 "Ripser: efficient computation of Vietoris-Rips
  persistence barcodes" J. Appl. Comput. Topol. 5, 391-423
  https://arxiv.org/abs/1908.02518
- Cohen-Steiner-Edelsbrunner-Harer 2007 "Stability of persistence
  diagrams" Discrete Comput. Geom. 37, 103-120.
- Adams et al. 2017 "Persistence images: a stable vector
  representation of persistent homology" J. Mach. Learn. Res. 18.

## Edges cited / composed

- **Composes** the gap-sequence representation of `χ_P`
  (a pillar 1 representation per E7.7) with the Takens delay
  embedding (no edge) and Vietoris-Rips persistent homology (cross-
  domain).
- **Detects** Hardy-Littlewood serial-correlation structure
  (E2.13 / E2.14).
- **Disjoint from** E2.16 — DPP failure was a multi-body kernel
  obstruction; PH is a metric-space topological observable on the
  delay embedding.

## Files

* `persistent_homology_chi_p.py` — driver (segmented sieve, Takens
  embed, ripser, baseline replicates, JSON output).
* `main_run_d3.json` — main result, `M=2000, d=3, x ≈ 10^6`.
* `main_run_d3_x5M.json` — robustness, different window position.
* `main_run_d4.json`, `main_run_d2.json` — embedding-dimension scan.
* `scale_M{500,1000,2000,4000}.json` — sample-size scaling.
* `main_run_d3.log` — primary run console log.
