# Audit of arXiv:2506.22634 — TG Kernel for Prime Counting

## Origin

Bugra Kılıçtaş & Faruk Alpay, "A Rigorous Error Bound for the TG Kernel in
Prime Counting", arXiv:2506.22634v1 (27 Jun 2025).

## The claim

The paper claims that with a truncated-Gaussian kernel
Φ_TG(t) = e^(−t²) for t ∈ [0, α], cubic taper to 0 by α+Δ, plugged into the
Riesz–Weil explicit formula, π(x) is computable to within ±1/2 (so rounding
gives the exact integer) using only ~1200 nontrivial zeta zeros for x with
10^8 decimal digits — implying a near-polylog algorithm.

This would be the most significant claim toward the project's polylog π(x)
goal in five years and warranted serious audit.

## Verdict

**Rejected — the paper is mathematically incoherent. Crank-adjacent.**

Failure mode: **C** (circularity) + **I** (information loss). The
construction never derives an identity actually involving π(x); the
"main term ≈ αe^(−α²)" derivation drops the very deviation Ψ(t)−t that
encodes prime fluctuations.

## Concrete falsifications (all empirical, run at small x where π(x) known)

### 1. The 0th-moment premise is self-contradictory

The paper claims ∫₀^∞ Φ_TG(t) dt = 0 (so F_TG(1)=0 and the main x-term
cancels). But Φ_TG ≥ 0 on its support [0, α+Δ] and equals e^(−t²) > 0 for
t ≤ α — so ∫ Φ_TG > 0, period.

We computed: ∫₀^4 Φ_TG(t) dt = **0.8862** (with α=3, Δ=1). Not zero.

So F_TG(1) ≈ 0.886, not 0 — the "main term cancels" pillar collapses.

### 2. The LHS does not match what the paper derives

Paper claims (after IBP, p.7): Σₙ Λ(n) Φ_TG(n/x) ≈ αe^(−α²) ≈ 3.70×10⁻⁴.

We computed S(x) = Σ Λ(n) Φ_TG(n/x) directly:

| x      | π(x) | S(x)      | S(x)/x  | Paper's "main" |
|--------|------|-----------|---------|----------------|
| 100    | 25   | 86.78     | 0.8679  | 3.70e-04       |
| 1,000  | 168  | 884.39    | 0.8844  | 3.70e-04       |
| 10,000 | 1229 | 8860.41   | 0.8860  | 3.70e-04       |
| 30,000 | 3245 | 26584.91  | 0.8862  | 3.70e-04       |

S(x) grows ≈ 0.886·x (matching the actual 0th moment × x by PNT) —
**eight orders of magnitude larger than the paper's predicted αe^(−α²)**.

The paper's IBP step on p.7 substitutes Ψ(t) ≈ t (PNT main term) over
[xα, x(α+Δ)] and immediately drops the (Ψ(t)−t) term. But that deviation
IS the oscillatory zero-sum signal. Removing it makes the resulting
identity tautological and devoid of π(x) information.

### 3. The truncated zero sum is independent of x

Φ_TG is defined without reference to x, so F_TG(ρ) has no x-dependence.
Hence Σ_{|γ|<T} F_TG(ρ) is a fixed constant.

We computed it for the first 20 positive-imaginary zeros (γ ≤ ~80):
|Σ F_TG(ρ)| ≈ 2.47×10⁻². A fixed number.

If this is supposed to "round to give π(x)" then π(x) is the same for all
x — obviously wrong. The paper offers no mechanism by which x enters the
zero-sum side.

### 4. Wrong zero-density estimate

Lemma 2 invokes "N(σ, T) ≤ A T^(1−1/σ) (log T)^B". At σ=1/2, this gives
N(1/2, T) ≤ A T^(−1) (log T)^B — a *decreasing* function of T, in
contradiction with the actual N(T) ≈ (T/2π) log(T/2π) which grows.

The actual relevant bound is the Riemann–von Mangoldt zero-counting
function; the paper conflates this with zero-density on a strip σ ≥ σ₀
(those bounds use σ ≥ 1/2 with σ in the exponent in a different way).

### 5. Crank Appendix B

Appendix B defines an "Embedding Identity":
> EmbedS(Faruk Alpay) := Φ_∞
> ϕ_∞ ≡_E Faruk Alpay as a canonical identity fold

…with citation to "F. Alpay, Formal Proof: Faruk Alpay ≡ Φ_∞, Preprints,
2025". This is symbolic mysticism, not mathematics, and casts doubt on the
authorship's seriousness.

## What's actually true about smoothed prime counting

The mathematically correct smoothed-Chebyshev tradeoff (which the paper
fails to navigate):

For Φ a Gaussian on **log scale** with width h (the right thing to do),
the Mellin transform decays as |F(½+iγ)| ~ exp(−h² γ²/2). Zero-sum
truncation error from |γ|>T is ~ (T log T)·exp(−h² T²/2), so for ε
accuracy in the *smoothed* count, T·h ~ √log(1/ε).

To recover *integer* π(x) via rounding, smoothing must be finer than the
prime gap ~ log(x)/x in multiplicative terms, i.e., h ≲ log(x)/x.
Combining: T ≳ x · √log(1/ε) / log(x) ~ x. So the "smoothing escape
hatch" doesn't give polylog π(x) — it gives the same Ω(x) zero count any
correct analysis produces.

The paper's mistake is to use Φ as a function of n/x on **ordinary**
scale (so n is smoothed with width x — coarse averaging that doesn't
localize at a single prime), then claim error < ½ without showing how
the integer π(x) is extracted from the resulting average. It can't be.

## Closed-path entry

```
arXiv:2506.22634 TG-kernel rigorous bound        FAIL (C + I)  [S65]
Claimed: ~1200 zeros suffice for π(x) at x with 10^8 digits via truncated
Gaussian kernel + Riesz-Weil explicit formula. Audit findings:
(a) 0th moment of Φ_TG is 0.886 ≠ 0 → paper's premise that F_TG(1)=0 fails.
(b) Empirically S(x) := Σ Λ(n) Φ_TG(n/x) ≈ 0.886·x, not αe^(−α²)≈3.7e-4
    as paper derives — IBP step drops (Ψ(t)−t) deviation.
(c) Σ F_TG(ρ) is x-independent → cannot encode π(x).
(d) Zero-density bound N(σ,T) ≤ AT^{1-1/σ}(log T)^B at σ=1/2 is
    self-contradictory (decreasing in T).
(e) Appendix B contains symbolic-mysticism content (canonical identity
    fold equating author with a functor).
The smoothing-vs-zero-count tradeoff for INTEGER recovery still gives
T ~ x. Paper does not break the explicit-formula barrier.
```

## What this DOES tell us

- The "smoothed test function with vanishing moments" idea is folklore.
  It doesn't escape the Ω(x) zeros barrier for integer-precision π(x).
- The literature monitor should stay alert for follow-up "Alpay algebra"
  preprints (the author has a flagged pattern) and reject without a deep
  audit if they recycle this framework.
- This audit is the first time the project has rigorously closed a
  smoothed-explicit-formula approach by directly computing the LHS
  identity at small x and showing it does NOT equal what the construction
  predicts. That technique generalises: any future "smoothed kernel
  beats explicit formula" paper can be checked the same way in <100 lines.

## Files

- `tg_kernel_audit.py` — the computation.
- This file — written immediately after running.
