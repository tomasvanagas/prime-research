# Session 216 — D8 wild swing: DARTS for TC^0 PRIMES at N=12

**Date:** 2026-04-29
**Mode:** WILD SWING (frontier attack, ATTACK_VECTORS.md §D8).
**Target:** §D8 — "Differentiable architecture search (DARTS) for π(x)
circuits". 1-2 session attack with high A-grade probability if PRIMES
admits a small TC^0 circuit; B-grade structural closure if not.
**Cross-domain technique imported:** Differentiable Architecture Search
(Liu-Simonyan-Yang ICLR 2019) — `CROSS_DOMAIN_TECHNIQUES.md §8` row
PROPOSED → USED-I.
**Self-grade:** **B** — substantive cross-domain refinement of S84
oddness mechanism. The attack imported DARTS cleanly, ran 9 model fits
(3 conditions × 3 seeds) at N=12 / depth-3 / G_1=G_2=12, and produced
informative structural failure: PRIMES converges to the parity floor,
calibrated-1-bit random matches PRIMES within seed noise (p=0.14),
argmax discretisation collapses to constant zero. The cross-domain
technique survives import but is not load-bearing for any new
conclusion — depth-2 SAT (S84) and depth-3 DARTS (S216) find the
same predictor.

## Why this attack was chosen

Wild Swing prompt — pick highest-A-grade-probability vector not yet
attempted. The default order:

- §C1 Odlyzko zeros — attempted S71 (closed).
- §A1 SAT-search TC^0 — attempted S84 (partial closure).
- §B1 polynomial method — closed S92.
- §A3 Cayley spectral — attempted S79.
- §D4 Szegedy walk — S80; closed S141.
- §C2 orders 4-6 zero correlations — closed S123.

All defaults were exhausted. §D8 was named in `session207_d9_bgk_sum_product.md`
synthesis as "high A-grade probability" alongside D11, untried, and
imports a cross-domain technique (DARTS) at `CROSS_DOMAIN_TECHNIQUES.md
§8` PROPOSED status. Tractable in 1 session with PyTorch. Complementary
to S84 SAT/ILP at N=8.

## What was built

`experiments/ml/darts_primes/`:

- `darts_primes.py` — DARTS-style depth-3 differentiable circuit search
  with gate library `{AND, OR, XOR, MAJ_3, MAJ_5, MAJ_7, ID, NOT}`,
  soft input selection via sigmoid logits β, gate temperature anneal.
- `baselines.py` — reference predictor losses (constant, 1-bit oddness,
  mod-6, mod-30, mod-210 wheels).
- `analyze.py` — Welch t / Mann-Whitney comparisons, architecture
  inspection.
- `analyze_full.py` — combined PRIMES vs density-random vs calibrated
  pairwise tests.
- `extrap_eval.py` — discrete-circuit eval on held-out windows
  ([2^N, 2^N+1000), [2^(N+1), 2^(N+1)+1000), [10000, 11000),
  [10^5, 10^5+1000)).
- `calibrated_control.py` — parity-matched random control runner.
- `darts_primes_results.md` — results write-up with falsifiers.

## Reference baselines (N=12, density 0.1377)

```
constant predictor BCE  = 0.4008  (entropy floor)
1-bit oddness BCE       = 0.2961
mod-6 wheel BCE         = 0.2298
mod-30 wheel BCE        = 0.1905
mod-210 wheel BCE       = 0.1614
```

## Results

### Final BCE loss (100 epochs, 3 seeds each)

| Condition | mean ± std | range |
|-----------|------------|-------|
| PRIMES (chi_P)                      | **0.2969 ± 0.0000** | [0.2969, 0.2969] |
| matched-density random              | **0.3967 ± 0.0113** | [0.3873, 0.4092] |
| calibrated-1-bit random             | **0.2963 ± 0.0004** | [0.2959, 0.2967] |

PRIMES converges deterministically to the parity-conditional entropy
floor (zero variance to 4 decimals). Random matched-density converges
to the unconditional entropy floor (0.40). The calibrated control —
random function with the same parity-conditional density as PRIMES
(563 primes among 2048 odd n's, 1 prime among 2048 even n's, randomly
placed) — converges to ≈ 0.2963, **statistically indistinguishable
from PRIMES**.

### Pairwise Welch t-tests on final loss

| Comparison | t | p | mean gap |
|------------|---|---|----------|
| PRIMES vs RANDOM_DENS  | −15.286 | 0.0043 | −0.0998 nats |
| PRIMES vs CALIBRATED   | +2.343  | 0.144  | +0.0006 nats |
| RANDOM_DENS vs CALIB.  | +15.362 | 0.0042 | +0.1004 nats |

The PRIMES-vs-density-random gap of 0.0998 nats is **exactly reproduced**
by the calibrated-vs-density-random gap of 0.1004 nats. The
PRIMES-vs-calibrated gap is 0.0006 nats — within seed noise. **Oddness
is sufficient AND exhaustive** for the DARTS PRIMES advantage at this
search-space scale.

### Discrete-circuit accuracy

The argmax architecture is hard-evaluated on test windows:

| n-window | discrete_acc | majority baseline | gap |
|----------|--------------|-------------------|-----|
| train [0, 4096)            | 0.8623 | 0.8623 | +0.0000 |
| [2^N, 2^N+1000)            | 0.8840 | 0.8840 | +0.0000 |
| [2^(N+1), 2^(N+1)+1000)    | 0.8890 | 0.8890 | +0.0000 |
| [10000, 11000)             | 0.8940 | 0.8940 | +0.0000 |
| [10^5, 10^5+1000)          | 0.9190 | 0.9190 | +0.0000 |

**The discrete circuit is the constant-zero function at every window.**
Argmax on the best PRIMES architecture yields:

```
Layer 1: 12 nodes, all ID gates, each picking 1 input wire
Layer 2: 12 nodes, all MAJ_5 gates, with empty selected-input set
Layer 3: 1 MAJ_7 gate with empty selected-input set
```

With empty selected input sets, MAJ_k(empty) = sigmoid(τ · (−k/2)) → 0
for large τ. The continuous-relaxation soft training advantage doesn't
survive discretisation — known DARTS failure mode (Zela et al. ICLR
2020 https://arxiv.org/abs/1909.09656).

## Verdict

**B-grade closure of §D8** (mode I — Information loss in argmax
discretisation, AND signal reduction to known oddness mechanism).

1. DARTS at depth 3 with gate library `{AND, OR, XOR, MAJ_3/5/7, ID,
   NOT}` and G_1 = G_2 = 12 **converges to the 1-bit oddness predictor
   on PRIMES** (loss 0.297, matching the parity-conditional entropy
   floor `H(prime | bit_0)` exactly).
2. The PRIMES advantage over matched-density random (0.10 nats) is
   **fully explained by oddness**: calibrated-1-bit random achieves
   the same loss within seed noise.
3. DARTS does NOT find mod-6 (loss 0.23) or mod-30 (loss 0.19)
   baselines at this depth/width. No super-oddness modular structure
   is extracted.
4. Argmax discretisation collapses to constant-zero. The continuous-
   relaxation prime signal is information-lost — this is the DARTS
   "argmax-collapse" failure mode (Zela 2020).
5. Adds **36th orthogonal pseudorandomness measure** (DARTS continuous-
   relaxation reachable BCE) with the same elementary mechanism as
   S84 (oddness mod 2). `novel/pseudorandomness_of_pi.md` updated.
6. **Adds new edge E2.33** (DARTS BCE saturation at parity floor;
   argmax-collapse to constant zero).

## Why B and not A

The session did not find a non-trivial generalising depth-3 PRIMES
circuit (which would have been the A-grade outcome). PRIMES advantage
saturates at oddness, which S84 already documented at N=6 / depth 2 /
SAT. The cross-domain DARTS technique imported cleanly but is not
load-bearing — it reproduces S84's mechanism in a different paradigm
without extending it.

## Why not C

The session produced 4+ artefacts not previously in the project:

- First DARTS application to PRIMES at any N. Cross-domain technique
  Liu-Simonyan-Yang ICLR 2019 import status PROPOSED → USED-I.
- First **calibrated-1-bit random control** for PRIMES circuit
  search. The S84 follow-up question "is there a residual PRIMES-vs-
  random gap beyond the oddness effect?" is now empirically answered
  at depth 3 / N = 12: NO, calibrated control matches PRIMES at p=0.14.
- DARTS argmax-collapse documentation specific to chi_P (Zela 2020
  generic phenomenon, applied here to a number-theoretic target).
- New EDGE E2.33; refined `novel/pseudorandomness_of_pi.md` to 36
  measures; new CLOSED_PATHS row 853.

## Why not F

Committed to ONE wild-swing target (§D8) and spent the full session on
it. No mid-session pivot. Honest negative result with structural
mechanism. Two new ATTACK_VECTORS sub-vectors proposed (D8.a, D8.b)
in the closure section.

## Q1 — What was produced not in the project before this session

A: DARTS continuous-relaxation BCE measurement on chi_P with three
matched controls (matched-density random, parity-calibrated-1-bit
random) at N=12, depth 3. The calibrated-1-bit-random control was the
S84 follow-up question previously sitting open in NOVELTY_CHALLENGES.
Empirical answer: parity is exhaustive at this depth/width.

## Q2 — Edges cited / composed

- E1.10 / E3.13 / `novel/pseudorandomness_of_pi.md` — pseudorandomness
  battery; this is a 36th measure.
- E5.3 / E7.10 — PRIMES ∈ TC^0 status; depth-vs-modulus orthogonal.
- S84 (`session84_a1_sat_tc0_primes.md`) — depth-2 SAT W=1 oddness
  mechanism, now confirmed at depth 3 / N = 12 / DARTS paradigm.
- E2.33 (NEW) — DARTS BCE saturation at parity floor.

## Q3 — Was this duplicate-closure?

A: No. The session imported a NEW cross-domain technique (DARTS,
PROPOSED → USED-I), ran the calibrated-1-bit control that was a
PRE-EXISTING NOVELTY_CHALLENGES follow-up from S84, and added a new
edge E2.33 with falsifiers and successors. The mechanism (oddness)
duplicates S84, but the measurement, the calibrated control, and the
edge are new.

## Q4 — Next-action

(a) Two successor vectors added to ATTACK_VECTORS.md §D.D8 closure:
    - **§D.D8.a** — DARTS at depth ≥ 4, G ≥ 32 with XOR-favoured init.
    - **§D.D8.b** — Boolean STE (Bengio 2013) + Gumbel-Softmax
      discretisation paradigm.
(b) `frontier_gen` mode should fire if open ATTACK_VECTORS count drops
    below 4 (currently several open: A2, A4, A6, B3, B4, B5, C3, C6,
    D1, D6, D11, D12, D14-19, D21, D23-26, D28, D33-47, G3) — D8.a/b
    are added but not all "open" vectors are easily tractable.
(c) `commit_state` recommendation per S215 was already to draw from
    new sources (TDA on π(x) fine scales; free-probability spectral;
    transfer-operator at Connes-Bost). S216 closure does not change
    that recommendation.

## Files produced

- `experiments/ml/darts_primes/darts_primes.py`
- `experiments/ml/darts_primes/baselines.py`
- `experiments/ml/darts_primes/calibrated_control.py`
- `experiments/ml/darts_primes/analyze.py`
- `experiments/ml/darts_primes/analyze_full.py`
- `experiments/ml/darts_primes/extrap_eval.py`
- `experiments/ml/darts_primes/darts_primes_results.md`
- `experiments/ml/darts_primes/run/results.json`
- `experiments/ml/darts_primes/run/calibrated_results.json`
- `experiments/ml/darts_primes/run/darts_run.log`
- `experiments/ml/darts_primes/run/darts_calibrated.log`
- `experiments/ml/darts_primes/run/analyze_full_output.txt`
- `experiments/ml/darts_primes/run/extrap_eval_output.txt`
- `archive/sessions/session216_d8_darts_primes.md` (this file)

Updated: `ATTACK_VECTORS.md` (D8 closed), `EDGES.md` (E2.33 added),
`CROSS_DOMAIN_TECHNIQUES.md` (§8 DARTS PROPOSED → USED-I),
`novel/pseudorandomness_of_pi.md` (35 → 36 measures),
`status/CLOSED_PATHS.md` (row 853).

## Channelled mathematician

**Williams** (algorithmic-method circuit-lower-bound thinking — circuit
search as a probe rather than a witness). Williams 2014 "Nonuniform
ACC circuit lower bounds" recasts circuit-search as algorithmic
probing of the limits of small classes; the present session's
finding that DARTS at depth 3 cannot extract anything beyond
parity at N=12 fits this frame: the absence of a small generalising
circuit is itself information about the structure (or its absence)
in the target function.
