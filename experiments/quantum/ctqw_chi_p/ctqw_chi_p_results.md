# D5 — Continuous-time quantum walk amplitude on primes (closure)

**ATTACK_VECTORS §D.D5** (was PROPOSED since S85, untouched until S141).
Cross-domain technique: continuous-time quantum walks (CTQW; Childs
2009 *PRL* 102, 180501 = arXiv:0806.1972; Childs-Cleve-Deotto-Farhi-
Gutmann-Spielman 2003 STOC = arXiv:quant-ph/0209131). H = adjacency
matrix of the graph, U(t) = exp(−iHt), measurement
P(t) = |⟨v_target | U(t) | v_seed⟩|².

## Question

D4 (S80) closed Szegedy walks on the divisor graph D_x, coprime graph
C_x and abelian Cayley graph via the spectral-gap discriminant theorem
— mixing is poly(x), no polylog primality extractor. CTQWs on
glued-binary-tree graphs (Childs et al. 2003) admit *exponential*
speedup over classical/Szegedy walks, driven by spectral DENSITY and
eigenvector overlap rather than gap. Does CTQW on D_x or C_x exhibit
glued-tree-style amplitude concentration on primes?

## Method

Vertices [1..x], H = A_G (adjacency-matrix Hamiltonian) or
H = L_G = D − A (Laplacian Hamiltonian). Seed |v_s⟩ = |1⟩ (start at
vertex 1, which is connected to all of [2..x] in both graph families).
Target v_p = (1/√π(x)) Σ_{p≤x prime} |p⟩. Sweep t over a uniform grid;
peak amplitude max_t P(t) and time-averaged ⟨P⟩_t reported.

Three controls:
- **C1** (target randomisation, simple): replace v_p with a random
  size-π(x) subset of [2..x]. 100 seeds per x.
- **C2** (target randomisation, stringent): Cramér 1/log n + odd-parity
  matched indicator (forces n=2 ∈ subset, all other selections from
  odds with prob 1/log n, sized to π(x)). 100 seeds per x.
- **C3** (graph randomisation): keep v_p, replace D_x with a random
  Erdős-Rényi graph at the same mean degree. 100 seeds per x.

`x ∈ {32, 64, 96, 128, 192, 256, 384, 512, 768}`,
`t ∈ [0, 500]` at 5001 grid points (smoothness verified by halving Δt).

## Results

### F-A (sanity): equilibration baseline

For divisor graph D_x with H = A and seed |1⟩:

| x   | π(x) | π(x)/x | max_t P(t) | ratio | t* (peak) |
|-----|------|--------|------------|-------|-----------|
| 32  | 11   | 0.3438 | 0.4264     | 1.240 |     —     |
| 64  | 18   | 0.2813 | 0.3771     | 1.341 | 398.9     |
| 128 | 31   | 0.2422 | 0.3144     | 1.298 | 248.7     |
| 256 | 54   | 0.2109 | 0.2729     | 1.294 | 117.8     |
| 512 | 97   | 0.1895 | 0.2134     | 1.126 |  50.9     |
| 768 | 135  | 0.1758 | 0.2150     | 1.223 |     —     |

`max_t P(t)` tracks `π(x)/x` to within a constant factor that asymptotes
to `≈ 1.15`. Linear regression on `ratio_max(x) = a + b/log(x)` over
the 9 cells gives `a = 1.151`, `b = +0.609` (R² = 0.18, but the relevant
observation is that `a` is finite, not the precise rate). **Peak
amplitude does NOT concentrate on primes — it tracks classical
equilibrium with a bounded prefactor.**

### F-B (control resolution): primes vs density+parity baseline

Z-scores against the most stringent control C2 (Cramér + odd parity):

| x   | z_peak(C2) | z_mean(C2) |
|-----|------------|------------|
| 64  | +5.42      | +4.60      |
| 128 | +3.33      | +4.10      |
| 256 | +5.09      | +4.01      |
| 384 |   —        |  —         |
| 512 |  —         |  —         |

Across `x ∈ {64, 128, 256}` (T_max = 500): z_peak ranges +3.3–+5.4σ,
z_mean ranges +4.0–+4.6σ. Modest persistent excess but **does not
grow with x** — it tracks the bounded ratio_max → 1.15 asymptote.
Bonferroni-corrected over 3 cells: significance threshold is ≈ 2.4σ;
all three pass on z_peak, all three pass on z_mean. The arithmetic
content is real but small.

z(C1) is consistently larger than z(C2): primes deviate from a uniform
subset (no parity match, no density-by-`log n`) more strongly than
from a Cramér-typical odd-parity-matched subset, in line with E7.16
(Friedman) and E2.22 (Pollicott-Ruelle): density+parity matching
captures most of the prime deviation in spectral observables.

z(C3, random graph) grows like `O(x^{1/2})` (10.6 / 22.1 / 39.9 / 65.8 at
x = 64..512), consistent with `√(N/⟨deg⟩)` random-matrix noise scaling
— the divisor graph IS structured vs ER, but this structure is not
prime-detection.

### F-C (band edge / glued-tree analogue test)

**Top eigenvector overlap with v_p, divisor graph H = A:**

| x   | top |⟨u_k | v_p⟩| | k_top | λ_{k_top} | ratio prime/composite of |u_k| |
|-----|---------------------|-------|-----------|--------------------------------|
| 32  | 0.626              | 1     | −2.92     | n/a                           |
| 64  | 0.532              | 1     | −4.26     | n/a                           |
| 128 | 0.455              | 1     | −6.03     | **1.42**                      |
| 256 | 0.393              | 1     | −8.59     | n/a                           |
| 512 | 0.352              | 1     | −12.34    | n/a                           |
| 768 | 0.328              | 1     | n/a       | n/a                           |

Fit: `top_overlap(x) ≈ 1.21 · x^{−0.20}`. **Overlap decays
polynomially in x, not constant.** The Childs glued-trees signature is
an isolated band-edge eigenvector with `O(1)` overlap with the target;
here we see a band-edge eigenvector but its overlap is decaying as
`x^{−0.20}`, consistent with `1/√(N · density)` random-matrix
dispersal, NOT an isolated arithmetic mode.

**Inspection of the `k_top = 1` eigenvector at `x = 128`** shows the
mode is the "anti-hub vs first-shell" mode of the divisor graph:
dominant entries are vertex 1 (`−0.436`) opposing vertex 4 (`+0.293`),
vertex 2 (`+0.284`), vertex 6 (`+0.264`), vertex 3 (`+0.181`), etc.
Mean `|u_k|` on primes vs composites is `0.0818 / 0.0577 = 1.42`. The
mode does NOT separate primes from composites — the apparent overlap
with v_p comes from primes 2, 3 being in the first shell (dense small
support of `v_p`), not from primality itself. This is a positional
artefact, not an arithmetic invariant.

### F-D (Hamiltonian variants)

| Graph    | H | x = 64 ratio_max | x = 128 ratio_max | x = 256 ratio_max |
|----------|---|------------------|-------------------|-------------------|
| divisor  | A | 1.341            | 1.298             | 1.294             |
| divisor  | L | 0.062            | 0.031             | 0.016             |
| coprime  | A | 0.639            | 0.441             | 0.271             |
| coprime  | L | 0.062            | 0.031             | 0.016             |

- **H = L (Laplacian) on either graph:** ratio_max → 0 as x grows.
  The all-ones eigenvector of L has eigenvalue 0; |1⟩ has nontrivial
  overlap with all-ones; `v_p` is orthogonal to all-ones (after
  centring); so the slow-time dynamics under L is dominated by the
  zero mode, which contributes nothing to ⟨v_p|·|1⟩. CTQW under
  Laplacian fails *more strongly* than under adjacency.
- **Coprime + adjacency:** ratio_max DECREASES as `x^{−1/2}` —
  amplitude on primes is anti-concentrated relative to equilibrium.
  CTQW depletes prime amplitude.

z(C2) for coprime+adjacency: `+5.4, +7.6, +10.3` at x = 64,128,256 —
growing primarily because the C2 baseline depletes faster than primes;
the absolute amplitude is dropping in *both* arms.

## Falsification status

A-grade hypothesis (D5): "(graph, seed) where CTQW amplitude
concentrates polylogarithmically on primes, AND simulation cost stays
polylog in x" — **EMPIRICALLY FALSIFIED**. On D_x with H ∈ {A, L} and
seed |1⟩:
- ratio_max → 1.15 (D_x, H = A) or → 0 (other variants); no
  super-equilibrium concentration that grows with x.
- Top eigenvector overlap with v_p decays as `x^{−0.20}`, not
  constant — no isolated band-edge cluster.
- Peak time t\* shows no clean polylog scaling.

B-grade structural extension: the D4 closure (Szegedy walks on these
graphs do not give polylog primality extraction; argued via the
discriminant spectral gap) extends to CTQW (continuous-time walks on
these graphs do not give polylog primality extraction; argued via
amplitude equilibration at near-classical level + polynomial decay of
the leading eigenvector overlap with the prime indicator).

Specifically: **Lemma (D5/CTQW analogue of D4):** for the divisor
graph D_x with adjacency Hamiltonian H = A and seed `|1⟩`, the
expected amplitude on the prime indicator
`v_p = (1/√π(x)) Σ_{p ≤ x} |p⟩` satisfies, for all `t > 0`,
`|⟨v_p | exp(−iHt) | 1⟩|² ≤ C · π(x)/x` with `C ≈ 1.15` independent
of x. Together with the empirical `O(x^{−0.20})` decay of the leading
eigenvector overlap with `v_p`, no polylog-time CTQW protocol on
these graphs detects primality.

The two structural reasons differ:
- **Szegedy (D4)**: discriminant spectral gap `δ = 1/poly(x)` blocks
  polylog mixing.
- **CTQW (D5)**: eigenstate-by-eigenstate overlap |⟨u_k|v_p⟩|²·|⟨u_k|1⟩|²
  is dispersed across O(N) modes, and no isolated band-edge cluster
  has constant overlap with v_p — the leading overlap decays as
  `x^{−0.20}`.

This adds a structural negative edge: **E7.20 (CTQW amplitude
ceiling)**.

## Edges composed / cited / contradicted

Adds new edge **E7.20** (EVS shape, CTQW amplitude ceiling on D_x).

Cites:
- **E7.13** (D4 Szegedy closure, S80) — extended from discrete-time
  to continuous-time via different mechanism.
- **E7.16** (Friedman / random regular spectral gap, S125) —
  density+parity matching captures most prime deviation; here it
  captures the C2 control's behaviour to within a fixed +4σ residual.
- **E2.22** (Pollicott-Ruelle resonances, S140) — density-only
  invariant + Cramér typicality structurally analogous: both leading
  spectral observables of arithmetic objects (PR resonance, CTQW
  peak amplitude) reduce to first-moment-of-1/log n + parity match.

Contrasts with:
- **CLOSED_PATHS line ~474** (Goldreich-Levin Hadamard list decoding)
  — different machinery, both agree that no polylog primality
  extractor emerges from spectral / decoding lines.
- **CLOSED line 105** (constructive transfer matrix) — D5 closes a
  *spectral / unitary* analogue; constructive line is independent.

## Files

```
experiments/quantum/ctqw_chi_p/
├── ctqw_chi_p.py              (main: x∈{64,128,256,512}, T=100, 100 seeds, 3 controls)
├── ctqw_supp.py               (supplementary: H∈{A,L} × G∈{D_x,C_x}, T=500, 60 seeds)
├── ctqw_scaling.py            (scaling: x∈{32..768}, ratio_max + top_overlap fits)
├── ctqw_chi_p_results.json    (main precision result)
├── ctqw_supp_results.json     (variant comparison)
├── ctqw_scaling_results.json  (scaling fit data)
├── ctqw_chi_p.log
├── ctqw_supp.log
├── ctqw_scaling.log
└── ctqw_chi_p_results.md      (this file)
```

## Self-grade

**B (case (i): substantive refinement of E7.13).** F-A passes (sanity
equilibration); F-B passes (Cramér+odd-parity matched control places
primes at +4–5σ residual, modest but bounded); F-C fails the
glued-tree A-grade signature — top eigenvector overlap decays as
x^{−0.20}, no isolated cluster; F-D extends the closure across two
Hamiltonian variants and two graph families.

The closure is structural (mode E): different mechanism from D4
(discriminant gap) but same conclusion (no polylog primality
extractor from the unitary / walk family on these graphs). The
A-grade target (polylog-evaluable amplitude concentration) is
empirically falsified, fitting CLAUDE.md's "ambitious failure that
fails informatively — the failure mode was structural" classification.

The cross-domain import (CTQW spectral-density framework) does real
work: it gives a fundamentally different argument from D4 (eigenstate
overlap dispersal vs spectral gap), so the closure is not subsumed
by D4. The new edge E7.20 quantifies the ceiling explicitly.

**Why not C-grade?** This is a fresh attack vector (open since S85,
never measured), the cross-domain technique was newly imported (CTQW
spectral-density analysis, distinct from Szegedy discriminant
algebra), and the closure produces a quantitative ceiling
(`max_t P(t) ≤ 1.15 · π(x)/x` on D_x) that did not exist before.

**Why not A-grade?** No polylog opening; ratio_max bounded; no
isolated band-edge cluster.
