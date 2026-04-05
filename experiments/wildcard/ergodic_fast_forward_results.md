# Ergodic Fast-Forwarding of Prime Gaps — Results

**Experiment:** Can prime gaps be modeled as a fast-forwardable dynamical system?

**Verdict:** CLOSED — FAILS. Positive Lyapunov exponent (chaotic), no linear recurrence structure.

## Key Findings

### 1. Conditional Entropy Decay
| Markov order | Cond. entropy | Reduction |
|---|---|---|
| 1 | 3.33 | 10.2% |
| 2 | 2.99 | 19.4% |
| 3 | 2.34 | 37.1% |
| 4 | 1.31 | 64.8% |
| 5 | 0.47 | 87.3% |
| 8 | 0.005 | 99.9% |
| 10 | 0.000 | 100.0% |

**CAUTION:** The 100% reduction at order 10 is a FINITE-SAMPLE ARTIFACT. With N=20000 gaps and order-10, most 10-tuples are unique (state space ~ 50^10 >> 20000), so each history maps to exactly one observed next gap. This does NOT mean gaps are order-10 deterministic.

### 2. Linear Recurrence — Essentially Zero
R^2 < 0.006 for all orders 1-10. Prime gaps have NO linear recurrence structure.

### 3. Lyapunov Exponent — CHAOTIC
| Steps | lambda |
|---|---|
| 1 | 1.06 |
| 2 | 0.75 |
| 5 | 0.38 |
| 10 | 0.16 |

All positive → the gap sequence is **chaotic**. Nearby trajectories diverge exponentially. Fast-forwarding a chaotic system requires computing all intermediate states.

### 4. False Nearest Neighbors
FNN drops to 0 at dim ~ 9, suggesting the gap sequence lives in a ~9-dimensional space. However, this is NOT useful for fast-forwarding because the dynamics within that space are chaotic (positive Lyapunov).

### 5. Comparison with Cramér Random Model
Surprisingly, **real gaps have HIGHER conditional entropy** than Cramér random gaps at order 3:
- Real: H = 2.34
- Cramér: H = 1.09

This means prime gaps are MORE unpredictable than the Cramér random model at short range! The Cramér model (exponential gaps scaled by ln(p)) has less diversity in gap values because it's a smooth model, while real gaps have more varied small factors.

## Implications
- No linear recurrence → can't use matrix exponentiation
- Chaotic dynamics → can't skip steps without computing each one
- Higher entropy than random model → prime gaps are maximally unpredictable at these scales
- No fast-forwarding possible for the gap dynamical system

## Failure Mode
**INFORMATION**: The gap sequence is chaotic (positive Lyapunov) with near-zero linear predictability. The apparent high-order "determinism" is a sample-size artifact. No dynamical shortcut exists.

**Session:** 40 (Fresh Perspective)
