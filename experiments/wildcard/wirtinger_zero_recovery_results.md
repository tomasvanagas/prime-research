# Wirtinger-flow recovery of γ_k from π(x) — results

**Question:** Can gradient-descent on the explicit-formula loss recover the
first K zeta-zero ordinates from a few values of π(x)? If so (and recovery
basin is large), the inverse problem is easier than the forward problem.

**Setup:** Generate π(x) values at x ∈ {20, 30, 50, 70, 100, 150, 200, 300, 500,
700, 1000} using the true first K zeros (oracle). Then run gradient descent on
`L(γ̂) = Σ_x (ψ_explicit(x; γ̂) − ψ_oracle(x))²` from various initializations.

## Convergence-basin diagnostic (perturbed init γ̂₀ = γ_true + ε·N(0,1))

| K  | ε=0.01           | ε=0.1            | ε=0.5            | ε=1.0            | ε=2.0            |
|----|------------------|------------------|------------------|------------------|------------------|
| 3  | dist=0.000 ✓     | dist=0.000 ✓     | dist=0.288       | dist=1.20        | dist=1.60        |
| 5  | dist=0.001 ✓     | dist=0.014 ✓     | dist=1.02        | dist=1.93        | dist=4.12        |
| 10 | dist=0.019 ✓     | dist=0.160       | dist=1.17        | dist=2.70        | dist=5.63        |

Recovery succeeds only when init is within **ε ≲ 0.1** of truth.

## Random initialization: 0/5 successes for every K ∈ {3, 5, 10}

| K  | trials | converged | best final dist |
|----|--------|-----------|-----------------|
| 3  | 5      | 0         | 1.50            |
| 5  | 5      | 0         | 5.17            |
| 10 | 5      | 0         | 14.11           |

## Verdict

**Loss landscape is exponentially non-convex.** The recovery basin shrinks fast
in K — radius ~0.1 around the truth at K=10. With unique γ_k spacing ~ 5–10
in this range, even *informed* random init lies outside every basin.

Why: the loss is a sum of M squared terms, each oscillating in (γ̂ log x). Each
spurious near-resonance creates a local minimum. A counting argument: at scale
γ_K ≈ 50 with M=11 evaluation points logarithmically spread, the loss has on
the order of `(γ_K log x_max)^K ≈ (50·log 1000)^K ≈ 350^K` competing minima.
For K=10 that's 10²⁵ — far beyond global-search reach.

## Failure mode

**E (Equivalence).** Inferring γ_k from π(x) values is at least as hard as the
forward problem. Information-theoretically reasonable (the explicit formula is
a continuous bijection), but optimization-theoretically *worse* than direct
computation. By "duality" no shortcut emerges.

## Comparison to phase retrieval (the original Wirtinger flow setting)
Candès-Li-Soltanolkotabi proved Wirtinger flow recovers signals from |Ax|² in
O(n log n) measurements with O(log n) iterations. The success there relies on
the *quadratic* (not oscillatory) loss landscape having a single basin once
spectral init is used. Our loss is a *sum of irrationally-related sinusoids*,
no analogous spectral init exists, and the K! permutation symmetry alone gives
a forest of local minima.

CLOSED. Cross-reference: this same conclusion via different framing was likely
hit by `iterative_zero_refinement` and `zero_convergence_proper` — but the
explicit basin-radius vs K table is a useful new diagnostic.
