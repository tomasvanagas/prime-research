# Session 65 — Fresh-perspective audit of arXiv:2506.22634 + phi 2D rank

**Date:** 2026-04-26
**Mode:** Fresh-perspective (instructed not to read CLAUDE.md / CLOSED_PATHS)

## TL;DR

Two experiments, both closing paths:

1. **phi(x,a) 2D low-rank** — Apparent exp-decay of SVD spectrum is
   illusory for integer-precision recovery; linear rank required.
   Pseudorandomness wall confirmed in a 22nd measure.

2. **arXiv:2506.22634 audit** — Code-backed empirical falsification of
   the Kılıçtaş–Alpay TG-kernel "rigorous bound" claim. Multiple
   independent flaws nailed at small x (where π(x) is known). The paper
   was already loosely "DEBUNKED S12/S30" in the literature monitor; this
   session sharpens the verdict to a verified-false code-backed audit
   and replaces the stale "NEEDS VERIFICATION" annotation.

No new viable direction opened. Confirms the S47 mature-state hypothesis.

## Context — what the literature scan found

A sub-agent literature search across 2025–2026 arXiv preprints, ECCC,
and conference accepted-paper lists returned three categories of work:

- **Marginal**: Aggarwal 2510.16285 (binary-search wrapper for p_n,
  O(√n · log⁴ n) — uses existing π(x) primitives, does not change the
  asymptotic).
- **Subconvexity / zero-density**: Guth-Maynard 2405.20552 with v2
  April 2026 (first improvement on Ingham's 1940 zero-density estimate),
  arXiv:2411.13791, arXiv:2508.02041 — sharpen π(x) error terms but
  yield no algorithmic break.
- **Headline candidate**: Kılıçtaş-Alpay arXiv:2506.22634 — claims
  exact π(x) via truncated-Gaussian explicit formula with ~1200 zeros
  for x with 10⁸ digits. This is the only 2025 paper whose scaling
  would, if true, be near-polylog. Sub-agent flagged for full audit.

No quantum, ML, tensor-network, lattice, or topological algorithmic
result for π(x) appeared in window. AlphaProof / Gemini Deep Think
worked on olympiad problems, not analytic NT.

## Result 1 — phi 2D low-rank

**Question:** Does the Meissel φ function, viewed as a 2D matrix
M[i,j]=φ(x_i, a_j), admit a polylog-rank reconstruction useful for
attacking the O(x^{2/3}) Meissel-Lehmer-Deléglise-Rivat barrier?

**Method:** Build K×K matrices under four framings:
- A: φ(x_i, a_j) raw
- B: φ - x·∏(1−1/p) (residual from Mertens smooth)
- C: φ/x normalised
- D: column-difference φ(x,a)−φ(x,a−1)

φ computed by inclusion-exclusion (sum over smooth squarefree d ≤ x).
Tested K ∈ {18, 40, 60}.

**Findings:**

| Scale | Frame | b_exp  | alpha  |
|-------|-------|--------|--------|
| 18    | B     | -0.376 | -2.33  |
| 40    | B     | -0.302 | -2.42  |
| 60    | B     | -0.327 | -2.59  |

Relative singular value decay rate stable around exp(-0.33·k) across
scales — looks compressible. **But ‖M‖_F grows ∝ x.** For integer-
precision recovery (±0.5) we need σ_{r+1}/σ_1 < 0.5/σ_1 = O(1/x), so

  required rank ≈ -log(0.5/σ_1) / 0.33

For K=60, σ_1≈6.6×10⁴ ⇒ required rank ≈ 35 — linear in K, not polylog.
The relative compressibility is exactly cancelled by the absolute
scale growth. Same pseudorandomness wall.

**Closed**, mode I (information loss). Adds 22nd
structural-pseudorandomness measure.

## Result 2 — TG-kernel audit

**Paper claim:** With Φ_TG(t) = e^{-t²} on [0, α], cubic taper to 0 by
α+Δ, plugged into the Riesz-Weil explicit formula Σ Λ(n) Φ_TG(n/x) =
F(1) − Σ_ρ F(ρ) + E_triv, only ~1200 nontrivial zeta zeros suffice
to recover π(x) by rounding for x with 10⁸ digits.

**Audit (`experiments/wildcard/tg_kernel_audit.py`):**

### Flaw 1 — 0th moment fails by 0.886
Paper requires ∫Φ_TG = 0 (so F(1)=0 cancels main term). But Φ_TG ≥ 0
on its support [0, α+Δ]; we measured directly:

  ∫₀⁴ Φ_TG(t) dt = 0.8862  (α=3, Δ=1, cubic taper)

So F_TG(1) = 0.886, **not** 0. Premise self-contradictory.

### Flaw 2 — LHS off by 8 orders of magnitude
Paper's IBP derivation (p.7) yields S(x) := Σ Λ(n) Φ_TG(n/x) ≈
αe^{-α²} = 3.70×10⁻⁴.

Direct computation:

| x      | π(x) | S(x)      | S(x)/x  |
|--------|------|-----------|---------|
| 100    | 25   | 86.79     | 0.8679  |
| 1,000  | 168  | 884.39    | 0.8844  |
| 10,000 | 1,229| 8,860.41  | 0.8860  |
| 30,000 | 3,245|26,584.91  | 0.8862  |

S(x) ≈ 0.886·x, **eight orders of magnitude larger** than the paper's
prediction.

**Trace of error**: in the IBP step, paper substitutes Ψ(t) ≈ t (PNT
main term) over [xα, x(α+Δ)] and drops the (Ψ(t) − t) deviation.
But that deviation IS the oscillatory zero-sum signal that encodes
π(x). Removing it gives a tautological identity.

### Flaw 3 — Σ F_TG(ρ) is x-independent
F_TG depends only on Φ_TG, not on x. So the entire RHS is a fixed
constant. We computed |Σ F_TG(ρ)| over first 20 positive-imaginary
zeros: ≈ 2.5×10⁻². A fixed number. Rounding it cannot produce
the x-dependent integer π(x).

### Flaw 4 — Wrong zero-density estimate
Lemma 2 invokes N(σ,T) ≤ A T^{1−1/σ} (log T)^B. At σ=1/2:
N(1/2, T) ≤ A T^{−1} (log T)^B — decreasing in T, contradicting
N(T) ≈ (T/2π) log(T/2π).

### Flaw 5 — Crank Appendix B
"EmbedS(Faruk Alpay) := Φ_∞ … canonical identity fold," citing
"F. Alpay, Formal Proof: Faruk Alpay ≡ Φ_∞, Preprints, 2025." Not
mathematics.

### What's actually true about smoothed kernels
Smoothing on log scale with width h gives Mellin decay e^{-h²γ²/2}.
Zero-truncation error from |γ|>T is ~ T log(T) e^{-h²T²/2}; for ε
accuracy, T·h ~ √log(1/ε). Integer-precision π(x) requires
smoothing finer than the prime gap, h ≲ log(x)/x, so

  T ≳ √log(2/ε) · x / log(x) ~ x

The same Ω(x) zero-count barrier. No polylog escape.

## Methodological note

The audit shows a fast falsification recipe for "smoothed kernel beats
explicit formula" preprints:

1. Compute the LHS Σ Λ(n) Φ(n/x) directly at small x where π(x) is
   known.
2. Compare with the paper's claimed RHS (main term + zero sum).
3. If they don't match, the IBP / explicit-formula application is
   wrong.

Cost: < 100 lines of mpmath + Python; falsified this paper in
minutes. Worth keeping in the literature-monitor playbook.

## Files

- `experiments/wildcard/phi_2d_lowrank.py` + `_results.md`
- `experiments/wildcard/tg_kernel_audit.py` + `_results.md`
- Updated `status/CLOSED_PATHS.md` (2 new entries)
- Updated `literature/state_of_art_2026.md` (S29 line replaced
  "NEEDS VERIFICATION" → "VERIFIED FALSE" with empirical evidence)
- `status/SESSION_INSIGHTS.md` S65 block

## Outcome

No breakthrough. Two paths closed cleanly. The audit is the most
useful product: it converts a vague "DEBUNKED" annotation into
a code-backed verdict and provides a generalisable method for
future Alpay-style preprints.
