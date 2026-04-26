# Session 45 Wildcard Brainstorm — 5 Genuinely Fresh Angles

Date: 2026-04-25
Scope: Each idea must (a) not appear by name in existing wildcard list, (b) be testable on n < 10000.

## Why these analogies?

| Breakthrough | Mechanism | Analog for π(x)? |
|---|---|---|
| Shor (factoring) | QFT exposes period | Spectral structure in oscillatory zero sum |
| Compressed sensing | Sparsity in transform basis | Wavelet/Padé sparsity of residual |
| AlphaFold | Learned energy landscape | NN learns δ(n) = p(n) − R⁻¹(n)? |
| Fast multipole | Hierarchical truncation of long-range kernel | The zero sum IS a long-range kernel — multipole moments? |
| Borel resummation | Divergent → fast-converging | Explicit formula is conditionally convergent |

## Selected experiments (each gets a script)

### 1. Borel resummation of the explicit-formula zero sum (`borel_explicit_formula.py`)
The truncated zero sum Σ_{|γ|<T} x^ρ/ρ converges only as T → ∞ at rate ~1/log T. Apply
Borel transform to the *coefficient sequence in T* and integrate. If the Borel
transform analytically continues to a tractable function on ℝ⁺, we collapse infinite
zeros into one integral evaluation. Free lunch test: compute partial sums S_K(x),
form Σ_k S_k z^k / k!, integrate against e^{−z}. Does it converge faster than S_K?

### 2. Padé approximants of the residual δ(n) = p(n) − R⁻¹(n) (`pade_residual.py`)
δ(n) is "random-like" in 21 measures, but Padé approximants exploit pole structure not
captured by Fourier/wavelet. If δ(n) has hidden meromorphic structure when interpolated
as a continuous function, even a low-order Padé should give super-Gaussian compression.
Compare RMSE of degree-(M, N) Padé vs polynomial of degree M+N.

### 3. Multipole expansion of the zero sum (`multipole_zero_sum.py`)
Treat γ_k as charges on the line {Re=1/2}, treat x as a probe. The "potential"
Σ x^{1/2+iγ}/(1/2+iγ) is a Coulomb-like interaction. Group γ's into clusters, expand
in spherical harmonics around cluster centers. If the multipole truncation gives
RMSE ε with O((log 1/ε)^d) terms, we have a fast-multipole prime counter.

### 4. Continued-fraction expansion of π(x)/li(x) (`cf_pi_over_li.py`)
The ratio approaches 1 with corrections of size ~1/√x. Khinchin's theorem says random
ratios have geometric mean of partial quotients ≈ 2.6854. If π(x)/li(x) has *bounded*
or *eventually periodic* CF, it's badly approximable in a structural way → hidden
algebraic identity. Compute CF for x = 10, 100, ..., 10⁵ to 100 digits.

### 5. Wirtinger-flow recovery of pi(x) from few zero ordinates (`wirtinger_zero_recovery.py`)
Pose: given oracle access to π(x) at small x, can we *recover the first few γ_k* via
nonconvex optimization (Wirtinger flow / phase retrieval)? Inverse direction of usual
question. If we can recover γ₁, ..., γ_K from O(K log K) values of π(x), then by
duality we can compute π at large x from O(K log K) zeros — gives effective polylog
*if* recovery is sub-exponential in K.

## What success looks like

Failure mode is probably C (circularity) or I (information loss) for all five.
But specific things to look for:
- Idea 1: Borel sum converges to wrong limit (essential singularity blocks it)
- Idea 2: Padé poles are random — confirms δ(n) genuinely irreducible
- Idea 3: Multipole order grows linearly with x — Coulomb analogy fails
- Idea 4: Partial quotients look like exponential distribution (Khinchin) — random
- Idea 5: Wirtinger landscape has exponentially many local minima

Even null results refine the barrier picture.
