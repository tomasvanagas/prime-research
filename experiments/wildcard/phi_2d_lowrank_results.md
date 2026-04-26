# phi(x, a) 2D Low-Rank / SVD Test

## Question

Viewed as a 2D matrix M[i, j] = phi(x_i, a_j), does the Meissel phi function
have low numerical rank or low effective sparsity? If yes (rank polylog in
input size), phi(x, a) — the bottleneck of Meissel-Lehmer-Deléglise-Rivat —
could be reconstructed from polylog samples, attacking the O(x^{2/3}) wall.

## Method

Built K×K matrices under four framings:
- **A**: phi(x_i, a_j) raw
- **B**: phi - x·∏_{p≤p_a}(1 - 1/p) (residual from Mertens smooth)
- **C**: phi/x (normalised)
- **D**: phi(x, a) - phi(x, a-1) (column-difference)

phi computed via direct inclusion-exclusion (sum over smooth squarefree d ≤ x
with all factors ≤ p_a).

Tested at three scales: K=18 (x ≤ 2e4), K=40 (x ≤ 2e6), K=60 (x ≤ 1.7e7).

## Results

| Scale       | Frame | b_exp  | alpha  | Verdict (relative-decay)              |
|-------------|-------|--------|--------|---------------------------------------|
| small (18)  | B     | -0.376 | -2.33  | EXP-DECAY (relative)                  |
| medium (40) | B     | -0.302 | -2.42  | EXP-DECAY (relative)                  |
| large (60)  | B     | -0.327 | -2.59  | EXP-DECAY (relative)                  |

The exponential-decay rate b_exp ≈ -0.32 is **stable across scales** for the
residual framing B. Naive reading: low effective rank.

## Verdict — but with the integer-precision correction

**FAIL — closed.** Failure mode I (information loss / pseudorandomness wall).

The relative singular value decay is exponential, **but** ‖M‖_F ∝ x grows
proportionally. To recover phi(x, a) to **integer precision** (±0.5, required
for exact π(x)) we need σ_{r+1}/σ_1 < 0.5/σ_1.

Quick calculation at large scale (K=60, framing B):
- σ_1 ≈ 6.6×10^4, so we need σ_r/σ_1 < 7.6×10^{-6}
- With decay rate ~0.33/k (natural log), required r ≈ -log(7.6e-6) / 0.33 ≈ **35**

Required rank as a function of K:
- K=18 → r ≈ 12-15
- K=60 → r ≈ 35

This is **linear** in K, not polylog. The exponential decay in *relative*
terms is exactly cancelled by the *absolute* growth of ‖M‖, which is the
classic signature of pseudorandomness when one needs integer-precision
reconstruction.

This matches the project's prior finding: pi-related sequences are
pseudorandom-like in 20+ structural measures. The phi 2D structure is a 21st.

## Key numbers

- Framing B large (K=60): σ_1 = 6.6e4, σ_8/σ_1 ≈ 1.5e-3, full rank > 1e-10.
- All four framings show similar qualitative behavior at all three scales.
- Total runtime: 60s.

## Closed-Path entry

```
phi(x,a) 2D low-rank reconstruction        FAIL (I)
Framing: SVD of phi(x_i, a_j) matrix at K=18,40,60. Relative singular
spectrum decays as exp(-0.33 k), but ||M||_F grows linearly in x, so
required rank for integer-precision recovery is linear in K, not polylog.
The relative compressibility is illusory — same pseudorandomness wall.
Session 65 fresh-perspective.
```

## What this rules out vs. what it does NOT rule out

- **Rules out**: a generic SVD-based "interpolate phi from polylog samples"
  approach using these four framings.
- **Does NOT rule out**: structured low-rank recovery using a basis adapted
  to the multiplicative structure of phi — e.g., decomposition by smooth
  squarefree generating functions, or hierarchical block-low-rank
  (HSS / HODLR) structure aligned with the recursion tree. Those would be
  separate experiments.
- **Does NOT rule out**: integer-recovery with side information (e.g., known
  π(x_anchor) at sparse anchors; reconstruct via low-rank + sparse
  correction). Untested.

## Files
- `phi_2d_lowrank.py` — implementation.
