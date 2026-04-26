# Session 48 — Causal-state / excess-entropy of prime parity stream

**Date:** 2026-04-25
**Mode:** fresh-perspective wildcard

## Goal

Find a wildcard angle not represented among the 57 existing wildcard scripts.
Test whether the prime parity stream q(n) = 1{p(n) ≡ 5 mod 6} has finite
excess entropy E (a *stochastic* memory measure, strictly more general than
the linear complexity over GF(2) which was already known to be N/2).

If E plateaus → there exists a finite-state stochastic generator → potential
polylog model. If E grows linearly → no such generator.

## Brainstorm (5 angles)

`archive/ephemeral/session48_brainstorm.md` lists five fresh angles:
1. Causal-state / excess-entropy complexity (selected, this session).
2. Cumulant expansion of δ(n) past order 2.
3. Phase retrieval of π(x) from |ζ(½+it)|² magnitude only.
4. Quasi-modular Eichler / mock-modular completions of log ζ.
5. Multipole expansion on the prime side (Riemann–Weil dual).

## Experiment

`experiments/wildcard/causal_state_complexity.py`. N = 148 931 prime parity
bits (p ∈ (3, 2·10⁶]). Block entropy H_L computed empirically for L = 1..18.
Entropy rate h_μ estimated from Δ_L = H_L − H_{L−1} averaged over L = 8..10
(moderate-L window — finite-sample bias dominates from L ≈ 12 since 2^L
approaches the sample size). Excess entropy E_L = H_L − L·h_μ. Compared
against an i.i.d. Bernoulli(p₁) null with same single-bit bias.

## Findings

| metric | prime parity | Bernoulli null |
|---|---|---|
| h_μ (Δ_L mean, L = 8..10) | **0.97686 nats / step** | 0.99708 |
| E_L peak (L = 10) | +0.04973 nats | +0.02389 |
| E_L − E_L(null), L = 10 | **+0.026 nats** | — |

**Per-step deficit** vs i.i.d.: 0.02022 nats ≈ 0.029 bits/step. Consistent
with Hardy–Littlewood pair correlations between consecutive primes' residues
mod 6 (e.g. twin primes p, p+2 force the (5, 1) residue pair).

**Plateau** of E_L: the excess entropy is bounded by ≈ 0.026 nats; it does
*not* grow with L within the measurement window. So an ε-machine of size
exp(0.026) ≈ 1.026 (essentially zero extra states beyond Markov-1) suffices.

## Verdict

**CLOSED.** Failure mode **I (Information Loss).** The genuine excess
entropy is real (positive, plateaus at ≈ 0.026 nats above the Bernoulli
null), but the per-step compressibility is only ~3.3 % below the entropy
ceiling. Far below what a polylog encoding of π(x) (which needs ≈ log₂ x
bits/step) would require.

## Files

- `experiments/wildcard/causal_state_complexity.py` (+ `_results.md`)
- `archive/ephemeral/session48_brainstorm.md` (5 angles)
- `status/CLOSED_PATHS.md` line appended (entry, S48)
- `novel/pseudorandomness_of_pi.md` table extended to **measure #24**:
  excess entropy of prime parity stream

## Why this entry belongs in the pseudorandomness table

It is the first *stochastic-memory* measure (the 23 earlier are mostly
deterministic / spectral / algebraic). The prime parity stream is the
first measure where the project finds a tiny but **non-zero** structural
deviation from a random null (h_μ < 1 by ≈ 0.020 nats). Earlier measures
all rounded to "indistinguishable from random". The deviation is small
enough that the random-like narrative still holds, but precise enough to
record.

## What's next

Steady-state continues. The remaining brainstorm angles 2–5 are testable
in future sessions; angle 4 (quasi-modular Eichler) needs a literature
sweep before code.
