# Proposal 25 — Liouville-Parity Triangulation

**Goal of the script.** Verify the closed-form identity linking π(x) parity to the Liouville summatory L(x) and the count C₃(x) of integers ≤ x with Ω(n) odd and ≥ 3, and probe whether a truncated/cheap version of L(x) (or C₃(x)) suffices for the parity bit.

## Verdict

CLOSED — identity is correct, but bottleneck simply shifts from π(x) to C₃(x), which is empirically random-like in the same sense.

## Verified identity

For all x ∈ [2, 10000]:
```
pi(x) = (x - L(x))/2 - C_3(x)
```
**Integer match: 9999/9999.  Parity match: 9999/9999.** The identity is correct.

The original formulation in the proposal mis-identified C₃(x) with (Q(x) − 1) (Q = squarefree count). That was wrong; corrected here. The actual relationship comes from
- #{n ≤ x : Ω(n) odd} = (x − L(x))/2  (standard Liouville-summatory split)
- #{Ω odd} = π(x) + #{Ω odd ≥ 3} = π(x) + C₃(x)

## Growth of the residual term

| x      | π(x)   | C₃(x)  | (x−L)/2 |
|--------|--------|--------|---------|
| 10     | 4      | 1      | 5       |
| 100    | 25     | 26     | 51      |
| 1000   | 168    | 339    | 507     |
| 10000  | 1229   | 3818   | 5047    |

C₃(x) grows roughly *3× faster than π(x)*. Asymptotically C₃(x) ≈ x · Σ_{k odd, k ≥ 3} (log log x)^{k−1}/((k−1)! log x), dominated by k = 3 giving ≈ x (log log x)² / (2 log x). Computing C₃(x) is at least as hard as π(x).

## Parity of C₃(x) — empirical agreement with simple proxies

| Proxy            | P(parity matches C₃) on x ∈ [2, 10000] |
|------------------|----------------------------------------|
| π(x^{1/3})       | 0.504                                  |
| π(x^{1/2})       | 0.494                                  |
| x ÷ 8            | 0.505                                  |
| L(x)             | 0.500                                  |

All proxies are at chance level. C₃(x) parity is empirically as random as π(x) parity. **No simple-arithmetic shortcut found.**

## Truncated-L probe

Define L_K(x) = sum_{n ≤ x, Ω(n) ≤ K} λ(n) (only count integers with at most K prime factors). Parity match with full L(x):

| K  | Agree fraction |
|----|----------------|
| 1  | 0.500          |
| 2  | 0.501          |
| 3  | 0.507          |
| 4  | 0.492          |
| 5  | 0.496          |
| 6  | 0.526          |
| 8  | 0.518          |
| 10 | 0.589          |

Capping Ω(n) at small K does not preserve L parity. Even K = 10 gives only 59% agreement on x ≤ 10000.

## Conclusion

The identity π(x) = (x − L(x))/2 − C₃(x) is exact. But:
1. C₃(x) appears to be at least as hard as π(x).
2. C₃(x) parity is uncorrelated with simple cheap proxies.
3. Truncating L by Ω-cap does not preserve parity.

**Path closed.** Add to status/CLOSED_PATHS as: "Liouville-parity triangulation — identity correct but C₃(x) parity is random-like; bottleneck not eliminated, only renamed." Failure mode: **E (Equivalence)** — the parity question reduces to an equally hard parity question on a different summatory.

## What WOULD make it work
A polylog-time algorithm for L(x) mod 2 *or* C₃(x) mod 2 specifically. Neither is in the literature; this experiment provides empirical evidence both look pseudorandom on the parity bit.

## Runtime
~1 second on x ≤ 10000.
