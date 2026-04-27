# C4 — Aggarwal × Dusart × BPSW unified `p_n` library

## Composition signature

A function `p_n(n: int) -> int` with three implementations and one
composition, each combining a subset of three EDGES.md edges:

| Mode      | Uses E6.6 | Uses E6.8 | Uses E5.1 | Description                                  |
|-----------|:---------:|:---------:|:---------:|----------------------------------------------|
| `agg`     | ✓         | ✓         | —         | Aggarwal binary search on `pi_lucy(x)`       |
| `bpsw`    | —         | —         | ✓         | BPSW-walk from 2                             |
| `hybrid`  | ✓         | ✓         | ✓         | Aggarwal narrow to width K, then BPSW-walk   |

`hybrid` is the C4 composition proper; `agg` and `bpsw` are subset
baselines that isolate the contribution of each component.

## Edge IDs cited

- **E6.6 (Aggarwal 2025 binary search optimality):** `p_n` reducible to
  `O(log n)` calls to `pi(x)` via binary search on the Dusart bracket.
  `literature/aggarwal_2025_analysis.md` §1.
- **E6.8 (Dusart bracket width = n):** for `n >= 6`,
  `n(log n + log log n - 1) <= p_n <= n(log n + log log n)`.
  Reference: Dusart 2010, 2018.
- **E5.1 (BPSW correctness ⇒ PRIMES in TC^0):** Strong-MR(2) ∧ Strong-Lucas
  with Selfridge `D` parameter. Verified deterministic to `2^64`; conditional
  above. `novel/bpsw_tc0_reduction.md`.

## Intended relationship to π(x)

The composition exploits the fact that `pi(x)` is the asymptotic bottleneck
in Aggarwal's bound `O(sqrt(n) log^4 n) = O(log n) · pi(x)_cost`. Once the
binary search has narrowed `[L, R]` to width `K`, every additional
narrowing step costs one full `pi(x)` evaluation, which dominates at large
`x`. Replacing the residual narrowing with BPSW-walk over `O(K)`
candidates avoids the last `O(log K)` `pi(x)` calls — substantively
cheaper when `pi(x)` is computed via Lucy DP / HKM (which scale as
`x^{2/3}` / `x^{1/2+o(1)}`).

The composition does NOT improve Aggarwal's worst-case asymptotic bound.
It improves the practical constant by replacing the trailing
`pi(x)`-call regime with cheap BPSW tests.

## Conditional propagation through Aggarwal's wrapper

The original `agg` mode is unconditional: `pi_lucy(x)` is exact (no
primality testing involved — Lucy DP only uses divisibility and
small-prime sieving). The `hybrid` mode is conditional on BPSW
correctness (E5.1) for any `p_n >= 2^64`. The conditionality propagates
**1-to-1**: a single BPSW pseudoprime in the final bracket `[L, R]`
shifts the answer by at most one prime. The wrapper does not amplify
the conditional — Aggarwal's binary search runs on `pi_lucy`, not on
BPSW, so the BPSW conditional enters only at the final walking step.

## Files

- `aggarwal_dusart_bpsw.py` — runnable code that builds + benchmarks all
  three modes.
- `aggarwal_dusart_bpsw_results.md` — falsification criterion and
  empirical outcomes.
- `definition.md` — this file.
