# wirsing_check.py — independent verification of A_W = 1

Companion script to `asymptote_extrapolation.py`. Computes the partial
Wirsing sum

    Σ_{sqf q ≤ Q} 1/φ(q) = A_W · log(Q) + B_W + o(1)

via direct sieve up to Q = 5×10⁶ and fits A_W, B_W from the high-Q tail.

## TL;DR

Linear fit on Q ∈ [50000, 5000000]:

    Σ = 0.999972 · log(Q) + 1.332969.

So **A_W = 1.0000 (4 sf), B_W = 1.3330 (4 sf)**. This independently
confirms S168's foundational assumption that Wirsing-A = 1.

## Why this is needed

S168 line 68 stated `A = 1` by Selberg-Delange and cited the empirical
value `1.04 at Q=5000`. This script:
- Verifies the asymptote A=1 directly (PASS).
- Discovers that the cited Q=5000 value is wrong: actual is 1.157,
  not 1.04 (10% error in S168's intermediate citation).

## Tables

| Q       | Σ 1/φ(q) over sqf | A_W partial = sum/log(Q) | sum − log Q |
|---------|-------------------|--------------------------|-------------|
| 100     | 5.911             | 1.283                    | 1.305       |
| 500     | 7.550             | 1.215                    | 1.335       |
| 1000    | 8.240             | 1.193                    | 1.332       |
| 5000    | 9.851             | 1.157                    | 1.334       |
| 10000   | 10.543            | 1.145                    | 1.333       |
| 50000   | 12.153            | 1.123                    | 1.333       |
| 100000  | 12.846            | 1.116                    | 1.333       |
| 500000  | 14.455            | 1.102                    | 1.333       |
| 1000000 | 15.148            | 1.096                    | 1.333       |
| 5000000 | 16.758            | 1.086                    | 1.333       |

The (sum − log Q) column converges to a constant 1.3326 by Q = 100K
(stable to 4 sf afterwards). Successive deltas: 1e-5 by Q=5M, so the
asymptote B_W = 1.333 is reached to ~5 decimals at this scale.

The A_W partial column converges from above. At Q = 5×10⁶, A_W
partial = 1.086 (still 8.6% above asymptote). To reach within 1% of
the asymptote, one would need Q ~ 10^14 (extrapolating the 1/log Q
correction).

## Falsifier

If A_W partial at Q = 5×10⁶ were materially > 1.10 (say, > 1.15), the
asymptote claim A = 1 would be in question. **Result:** A_W(5×10⁶) =
1.086, well below 1.10. PASS.

## What this DOES NOT verify

- The empirical SVD spike-block sum / 0.21·π(N) at d=14, 18, 20
  (untouched by this script).
- Whether k_*(N) is N^{0.42} (S74 territory).
- Whether the squarefree restriction is the right one for the chi_P
  decomposition (S166 / S168 territory).

## Code

Linear smallest-prime sieve to Q_max, then Euler totient and squarefree
flag computed by trial-division using smallest_prime[]. Bug fix from
first iteration: the inner `while n % p == 0` loop must COUNT the
multiplicity to detect non-squarefree q, not just divide them all out
in one step. (First iteration had this bug; the squarefree filter
returned True for all q. Caught by manual sanity-check at Q=30:
{4, 8, 9, 12, ...} were all flagged squarefree.)

## Cross-reference

- S168 §"Corollary — additivity over squarefree q" (line 67).
- S182 verify session: `archive/sessions/session182_verify.md`.
