# Session 53 — FOCUS-1 closure: Connes operator scaling

**Date:** 2026-04-25
**Mode:** normal
**Target:** TODO.md FOCUS-1 / EDGES.md E3.1 / Chain A step 4
**Outcome:** FOCUS-1 closed; Chain A no longer a polylog candidate.

## Mission

Measure the prime-budget → zero-count scaling law for the
Connes-Consani-Moscovici Nov-2025 spectral-triple operator
(arXiv:2511.22755, `literature/state_of_art_2026.md` §2.5b).

Anchor data point: B = 6 primes (2,3,5,7,11,13) yields the first 50
nontrivial Riemann zeros with errors ranging from 2.5×10^{-55} (zero #1)
to ~10^{-3} (zero #50). The crucial scaling law K_accurate(B) had not
been measured. If polylog: Chain A is a real polylog architecture
(huge). If linear / sub-linear: Chain A closes via Equivalence.

## Approach

Numerical proxy (NOT a faithful CCM reproduction):
discretize scaling operator D = -i d/du on L^2([-L, L]) with
L = log(λ), x_cutoff = λ² = 10^4, in Fourier basis with N=1200 modes
(matrix size 2401). Add self-adjoint rank-one perturbation V = c|v⟩⟨v|
with v encoding primes ≤ p_B via von-Mangoldt-weighted delta-comb in
log space:

    v(u) = Σ_{p^k ≤ e^L} (log p)/√(p^k) · δ(u − log p^k),
    Fourier-projected onto the basis.

Tested two perturbation kernels (delta-comb vs ψ-step Mertens pairing);
delta-comb selected by tuning. Coupling c tuned at B=6 to minimize
median matched error (best: c = −2.0). Swept B ∈ {1..9}. Diagonalized
with `numpy.linalg.eigvalsh` (double precision; mpmath multi-precision
not needed since errors are O(0.1)).

## Key findings

1. **Median matched error flat** at 0.13–0.20 across B = 1..9. Adding 9×
   more primes does not improve eigenvalue accuracy.
2. **K_accurate(<0.5, matched)** saturates at 50/50 for ALL B including
   the B=0 control — pigeonhole / comb-density artefact, not signal.
3. **Monotone test** (μ_K vs γ_K, not nearest): K=0 for all B at the
   0.5 threshold. The architecturally correct test shows zero
   information signal.
4. **Linear regression**: K_accurate(B) = −0.000·B + 50.000
   (slope 0, R² = 1.0).
5. **Reproduction fidelity to CCM is poor**: at B=6 the proxy gives
   err[1] ≈ 9.1×10^{-2}, vs CCM's published 2.5×10^{-55} — off by 53
   orders of magnitude. We could not faithfully replicate CCM's
   specific kernel (likely tied to Mellin transform of Riemann xi)
   without the paper's detailed construction.

## Three independent closure arguments (all → mode E)

1. **Rank-one parameter count.** A self-adjoint rank-one perturbation
   has ≤ B parameters (entries of v indexed by primes); by Cauchy
   interlacing it can shift at most ~B eigenvalues substantially while
   the rest interlace with the unperturbed comb. To encode ≫ B zeta
   zeros to high accuracy is information-theoretically constrained.
2. **Diagonalization cost (kernel-independent).** Even granting CCM's
   published per-zero accuracy at face value, spectrum extraction from
   an N×N truncation costs O(N³) = O(K³) for K eigenvalues. For
   K = √x, this is O(x^{3/2}) — strictly worse than the existing
   O(x^{1/2+ε}) zero-summation barrier. Iterative methods (Lanczos)
   give O(K²) per eigenvalue, still O(K³) for K eigenvalues. This
   argument is independent of how compactly the operator is
   parameterized.
3. **Geometric per-zero error growth.** CCM's published B=6 → K=50
   data range 10^{−55..−3} implies err(k) ≈ A · r^k with r ≈ 11.3,
   yielding K_accurate(B=6) ≈ 53 even at face value. CCM does not
   extrapolate K_accurate(B); the literature contains no multi-B data
   point that establishes super-linear scaling. The burden of proof
   rests on anyone claiming polylog.

## Honesty notes

- The proxy operator is NOT a faithful CCM reproduction — off by
  ~53 orders on err[1]. We tried two natural prime-encoded kernels
  (delta-comb and ψ-step) and tuned coupling. The CCM construction
  must use a more specialized perturbation (likely tied to the explicit
  Mellin transform of the Riemann xi function) not captured by Fourier
  projection of a delta-comb.
- However, the **scaling law** result is robust: any rank-one
  perturbation has B parameters, and the diagonalization-cost argument
  is independent of kernel choice.
- A future faithful CCM reproduction with multi-B data points would
  strengthen the closure further but is not load-bearing — argument 2
  alone suffices.

## Files

- `experiments/analytic/connes_operator/connes_operator_scaling.py`
- `experiments/analytic/connes_operator/connes_operator_scaling_results.md`
- `experiments/analytic/connes_operator/connes_operator_scaling_data.csv`
- `status/CLOSED_PATHS.md` (1 entry, line 696)
- `status/SESSION_INSIGHTS.md` (Session 53 entry)
- `TODO.md` FOCUS-1 marked DONE
- `EDGES.md` E3.1 / Chain A status updated

## Bigger picture

Chain A — the highest-EVS chain in the EDGES.md catalogue — collapses.
The remaining open frontiers in the project are:

- **Chain B**: π(x) mod q polylog for fixed q ≤ 13 (FOCUS-2). Cleanest
  single-target objective; no known polylog algorithm but no closure
  either.
- **Chain C**: Liouville parity polylog (FOCUS-3). Identity
  π(x) = (x − L(x))/2 − C₃(x) is exact; bottleneck is L(x) mod 2
  polylog feasibility.
- **Chain E**: AKS growing-dim MPOW alternative attacks (FOCUS-5).
  Highest stakes / lowest feasibility.

After 533+ closed paths and 53 sessions, no polylog architecture for
π(x) is credibly visible. The project remains in steady state.
