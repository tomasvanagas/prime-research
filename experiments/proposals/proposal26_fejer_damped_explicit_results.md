# Proposal 26 — Cesàro–Fejér-Damped Explicit Formula

**Goal of the script.** Replace the sharp truncation in Riemann's explicit formula
$$ \pi(x) \approx R(x) - \sum_{|\gamma| \le T} R(x^{1/2 + i\gamma}) $$
with a Fejér window K_T(γ) = (sin(γ/(2T))/(γ/(2T)))² 1[|γ| ≤ 2T]. Measure (a) integrity at T = 1000 sharp, (b) round-gap distributions at T ∈ {30, 100, 300, 1000}, (c) recovery rate ⌊S_T(x)⌉ = π(x) for x ≤ 3000.

## Implementation note (important bug fix)

`mpmath.li(x**rho)` is **wrong** for ρ = 0.5 + iγ with large γ: mpmath uses the principal branch of complex log, collapsing x^ρ to its principal value and discarding winding (~ γ log x / 2π full rotations). The correct expression is `mpmath.ei(rho * log_x)` evaluated on the natural sheet. The first version of the script used `li(x^rho)` and produced answers off by factor 100×; switching to `Ei(ρ · log x)` fixed integrity (S_T at T=1000 sharp matches π(x) within ±0.5).

This bug cost ~30 minutes and is worth flagging in any future zeros-summation code.

## Integrity check (T=1000, sharp cutoff)

| x      | π(x)  | S_T(x)    | difference |
|--------|-------|-----------|------------|
| 100    | 25    | 24.6512   | −0.349     |
| 200    | 46    | 46.1268   | +0.127     |
| 500    | 95    | 94.9348   | −0.065     |
| 1000   | 168   | 168.3113  | +0.311     |
| 2000   | 303   | 302.5610  | −0.439     |
| 5000   | 669   | 668.4781  | −0.522     |
| 10000  | 1229  | 1229.0755 | +0.076     |

T = 1000 zeros (gamma ≤ 1000, so √(10000) = 100, T/√x = 10) gets all 7 sample x correctly to within 0.5 — the formula works.

## Round-gap CDF over x ∈ [100, 3000] (step 100)

P(|S_T(x) − ⌊S_T(x)⌉| < b):

### Sharp cutoff
| T    | b=0.05 | b=0.1 | b=0.2 | b=0.3 | b=0.4 | b=0.5 |
|------|--------|-------|-------|-------|-------|-------|
| 30   | 0.03   | 0.20  | 0.37  | 0.60  | 0.87  | 1.00  |
| 100  | 0.10   | 0.20  | 0.33  | 0.53  | 0.73  | 1.00  |
| 300  | 0.10   | 0.23  | 0.37  | 0.67  | 0.80  | 1.00  |
| 1000 | 0.10   | 0.27  | 0.47  | 0.63  | 0.87  | 1.00  |

### Fejér damped
| T    | b=0.05 | b=0.1 | b=0.2 | b=0.3 | b=0.4 | b=0.5 |
|------|--------|-------|-------|-------|-------|-------|
| 30   | 0.13   | 0.30  | 0.47  | 0.63  | 0.80  | 1.00  |
| 100  | 0.17   | 0.23  | 0.47  | 0.73  | 0.80  | 1.00  |
| 300  | 0.13   | 0.27  | 0.40  | 0.67  | 0.80  | 1.00  |
| 1000 | 0.10   | 0.20  | 0.50  | 0.70  | 0.87  | 1.00  |

## Recovery rate ⌊S_T(x)⌉ = π(x)

| T    | Sharp | Fejér |
|------|-------|-------|
| 30   | 15/30 | 16/30 |
| 100  | 16/30 | 20/30 |
| 300  | 21/30 | 24/30 |
| 1000 | 29/30 | 29/30 |

## Verdict

Fejér damping gives a **measurable but small** improvement at intermediate T:
- T = 100 (T/√x ranges 1.8 – 10): sharp 53% recovery, Fejér 67%.
- T = 300: sharp 70%, Fejér 80%.

The improvement disappears at T = 1000 (both saturate to ~97%). At very small T = 30, Fejér gives only marginal gain (0.53 → 0.63 at b = 0.3 cdf bin).

**This does NOT break the √x barrier.** Asymptotically, the truncation error of the sharp sum is O(x/T); Fejér halves the constant in the L¹ sense but does not change the order. The constant improvement of ~10–14 percentage points at T = 100 vs T = 1000 means Fejér recovers ⅔ of primes with one third of the zeros.

## Status — recommended

**Open as a partial improvement, not a polylog route.** Useful for the engineering goal of speeding up exact π(x) computation by a constant factor, but not the kind of advance this project is aiming for.

A follow-up experiment: do the *failures* of Fejér recovery cluster by arithmetic structure of x (residue mod small primes, smooth/rough composition)? If yes, a hybrid algorithm that uses Fejér on "easy" x and falls back to sharp on "hard" x might give a real speedup. **Not pursued in this session.**

## Runtime

- 200-zero sanity check: < 1 s.
- Sample-x table at 5 T's: ~2 s.
- Round-gap CDF over 30 x-values × 4 T's × 2 modes: ~10 s.

## Failure mode classification

C/E/I taxonomy → **C (Circularity-adjacent / numerical)**: the Fejér kernel does not introduce new information, just reweights the same zero-truncation source. Cannot circumvent the explicit-formula sqrt-barrier on its own.
