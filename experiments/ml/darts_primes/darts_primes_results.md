# DARTS for PRIMES TC^0 Circuits — Results

**Attack vector:** `ATTACK_VECTORS.md §D8` — Differentiable Architecture
Search (DARTS) for π(x) circuits.
**Cross-domain technique:** Liu-Simonyan-Yang ICLR 2019
(`https://arxiv.org/abs/1806.09055`); Bender et al. ICML 2018; Baydin
et al. JMLR 18 (2018) — automatic differentiation.
**Session:** S216 (wild swing, 2026-04-29).
**Self-grade:** **B** — substantive cross-domain refinement of S84
oddness mechanism. Continuous-relaxation paradigm reproduces the S84
PRIMES-vs-random gap at N=12, but does NOT extend it: DARTS plateaus
exactly at the 1-bit oddness floor and finds no super-oddness structure
at depth 3 with the chosen gate library.

## What was tested

DARTS-style depth-3 differentiable circuit at N=12, gate library
`{AND, OR, XOR, MAJ_3, MAJ_5, MAJ_7, ID, NOT}`, soft input selection
via sigmoid logits β, gate-temperature anneal τ ≈ 2 → 10. Trained on
(n, χ_P(n)) for `n ∈ [0, 2^12)` (4096 inputs, 564 primes, density
0.1377). 3 seeds for each of three conditions (PRIMES, matched-density
random control, calibrated-1-bit random control).

`G_1 = G_2 = 12`, output 1, 100 epochs, Adam lr = 1e-2, BCE loss.
Architecture parameters (α, β, γ, log_τ) jointly optimised.

### Reference baselines (computed in `baselines.py`)

```
prime density           = 0.1377
constant predictor BCE  = 0.4008  (entropy floor)
1-bit oddness BCE       = 0.2961
mod-6 wheel BCE         = 0.2298
mod-30 wheel BCE        = 0.1905
mod-210 wheel BCE       = 0.1614
```

A predictor at loss `≈ 0.30` has identified the parity bit and is
applying conditional density to each parity class. A predictor at
loss `≈ 0.23` has additionally identified the multiples-of-3 residue
class. Anything below `≈ 0.20` requires constructing modular-arithmetic
features.

## Results

### Final BCE loss (100 epochs)

| Condition | seed 0 | seed 1 | seed 2 | mean ± std |
|-----------|--------|--------|--------|------------|
| PRIMES               | 0.2969 | 0.2969 | 0.2969 | **0.2969 ± 0.0000** |
| Random (matched ρ)   | 0.3934 | 0.4092 | 0.3873 | **0.3966 ± 0.0113** |
| Calibrated-1-bit (parity-matched ρ) | 0.2963 | 0.2959 | 0.2967 | **0.2963 ± 0.0004** |

PRIMES variance across seeds is zero to 4 decimals — the optimiser
converges deterministically to the parity floor. Random control
variance is ~0.01, hovering just below the entropy floor 0.4008.
The calibrated-1-bit control (matching PRIMES's parity-conditional
density: 563 primes among 2048 odd n's, 1 prime among 2048 even n's,
randomly placed) converges to **0.2963 ± 0.0004**, statistically
INDISTINGUISHABLE from PRIMES (gap = 0.0006 ± 0.0004 nats, well within
seed noise).

**This isolates the entire DARTS PRIMES-vs-random gap as a parity
artefact.** No residual prime-specific structure is detectable at
depth 3 beyond bit_0. This is the "calibrated-1-bit random" falsifier
proposed in `session84_a1_sat_tc0_primes.md` (Open question for next
session). Result: at this depth/width, the residual is ZERO within
seed noise.

### Welch t-tests on final loss (pairwise)

| Comparison | t | p | mean gap |
|------------|---|---|----------|
| PRIMES vs RANDOM_DENS  | −15.286 | 0.0043 | −0.0998 nats |
| PRIMES vs CALIBRATED   | +2.343  | 0.144  | +0.0006 nats |
| RANDOM_DENS vs CALIB.  | +15.362 | 0.0042 | +0.1004 nats |

The PRIMES-vs-density-random gap is **exactly reproduced** by the
calibrated-1-bit-vs-density-random gap (within seed noise). The
PRIMES-vs-calibrated gap is statistically zero. **Oddness is
sufficient AND exhaustive** for the DARTS PRIMES advantage at this
search-space scale.

### Discrete-circuit accuracy

All conditions achieve the majority-class baseline (0.8623 = 1 − ρ)
on the discretised circuit. NO condition discovers a circuit that
beats the constant-zero predictor on the binary level.

```
PRIMES       discrete_acc = 0.8623 (= majority baseline)
Random ctrl  discrete_acc = 0.8631 ± 0.0059
```

### Extrapolation [2^N, 2^N + 1000)

```
prime density on extrap window: 0.1160
majority baseline:              0.8840
PRIMES best discrete circuit:   0.8840  (= majority baseline)
```

The discrete circuit discovered for PRIMES at N=12 produces "always 0"
on inputs outside the training range — i.e., the architecture
*implements the constant-zero function* after argmax discretisation.
The continuous training advantage at PRIMES is fully buried in soft
mixture weights that fail to survive discretisation.

### Discovered architecture (best PRIMES seed, after argmax)

```
Layer 1: 12 nodes, all ID gates, each picking 1 input wire from {0..11}.
Layer 2: 12 nodes, all MAJ_5 gates, with empty selected-input set.
Layer 3: 1 MAJ_7 gate, with empty selected-input set.
```

With empty input sets, MAJ_k(empty) = sigmoid(τ · (-k/2)) → 0 for
large τ. So the discretised circuit is the constant-zero function.

The continuous architecture fits well to the 1-bit oddness predictor
during training (loss 0.297 ≈ entropy of `prime | bit_0 = 1`), but
the win lives entirely in the soft-mixture weights, NOT in the
argmax-discretised gate selection. This is a known failure mode of
DARTS (Zela et al. 2019 "Understanding and Robustifying Differentiable
Architecture Search").

## Verdict

**B-grade closure of §D8 (mode I — Information loss in argmax
discretisation).**

DARTS at depth 3 with the standard gate library:

1. **Reproduces** the S84 PRIMES-vs-random oddness gap at N=12. The
   PRIMES-vs-random loss gap of 0.10 nats matches `H(prime) − H(prime
   | bit_0)`, exactly the 1-bit oddness reduction.
2. **Does not extend** the gap. PRIMES converges to the oddness floor
   (0.297 ± 0.000), not to mod-6 (0.230) or mod-30 (0.190) baselines.
   DARTS does not synthesise modular-arithmetic features at this
   depth/width.
3. **Argmax discretisation is lossy.** The continuous-relaxation
   advantage doesn't survive: discrete accuracy hits the majority-
   class baseline. This is a known DARTS failure mode (Zela 2019).
4. **No A-grade signal.** The "discovered circuit at depth 3" hypothesis
   is falsified at this scale: the argmax circuit is degenerate
   (constant zero) regardless of training loss. To find a non-trivial
   PRIMES circuit, DARTS would need (a) different discretisation
   strategy (top-k > 0, gate set explicitly Boolean), (b) substantially
   larger depth/width than tested here, or (c) curriculum / structural
   priors (e.g., bias toward XOR for parity).

**Mode I (information loss):** the DARTS continuous relaxation contains
prime-specific signal during training (the 0.10-nat oddness gap) but
the signal is information-lost in the argmax-then-evaluate step. This
is structurally a 36th orthogonal pseudorandomness measure, with the
same "oddness only" mechanism documented at S84.

## What would falsify this

1. **Larger search space.** Re-run with `G_1, G_2 ≥ 32` and depth ≥ 4
   with explicitly XOR-favoured initialisation. If PRIMES converges
   to mod-6 (loss 0.23) or mod-30 (0.19), the depth-3 hardness claim
   is refuted.
2. **Different discretisation.** Replace argmax with top-k threshold
   + soft-to-hard rounding (Bender 2018 one-shot). If the PRIMES
   discrete circuit beats majority-class baseline (≥ 0.870), the
   "argmax loss is the bottleneck" claim is refuted.
3. **Calibrated-1-bit control showing PRIMES-residual gap.** This was
   tested directly: parity-matched random converges to 0.2963 ± 0.0004,
   indistinguishable from PRIMES at 0.2969. **Negative result —
   the DARTS gap is fully oddness.** Falsification would require a
   future re-run with a depth ≥ 4 / G ≥ 32 search showing PRIMES
   < 0.297 while parity-matched random plateaus at 0.297.

## Files produced

- `darts_primes.py` — DARTS-style training (3 conditions × 3 seeds)
- `baselines.py` — reference predictor losses
- `analyze.py` — Welch t / Mann-Whitney comparisons
- `extrap_eval.py` — discrete-circuit eval on held-out windows
- `calibrated_control.py` — parity-matched random control runner
- `run/results.json` — main experiment data
- `run/calibrated_results.json` — calibrated control data
- `run/darts_run.log` — training log

## Edges cited / composed

- **E1.10 / E3.13** — pseudorandomness battery (PRIMES indistinguishable
  from random under 35+ measures); DARTS = 36th measure with mode I.
- **`novel/pseudorandomness_of_pi.md`** — to update with this measure.
- **E5.3 / E7.10** — PRIMES ∈ TC^0 status; depth-vs-modulus orthogonality.
- **S84 (`session84_a1_sat_tc0_primes.md`)** — the SAT-based depth-2
  W=1 PRIMES-vs-random gap. THIS SESSION confirms the oddness mechanism
  at higher N (12 vs 6) and depth (3 vs 2) with a different paradigm.

## Cross-domain references

- Liu, Simonyan, Yang (ICLR 2019). "DARTS: Differentiable Architecture
  Search." https://arxiv.org/abs/1806.09055
- Bender, Kindermans, Zoph, Vasudevan, Le (ICML 2018). "Understanding
  and Simplifying One-Shot Architecture Search."
- Zela, Elsken, Saikia, Marrakchi, Brox, Hutter (ICLR 2020). "Understanding
  and Robustifying Differentiable Architecture Search."
- Baydin, Pearlmutter, Radul, Siskind (JMLR 18, 2018). "Automatic
  Differentiation in Machine Learning: a Survey."
