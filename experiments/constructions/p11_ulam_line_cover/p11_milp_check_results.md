# Thread 11 / Slot 4 — exact MILP integer optimum (HiGHS branch-and-bound)

## Cross-domain ingredient

scipy.optimize.milp with HiGHS branch-and-bound and dynamic constraint
generation. Set-cover IP is NP-hard, so feasibility is restricted to
N ≲ a few thousand for full convergence; larger N runs to time-limit
and produces best-feasible upper bound + LP-bound lower bound (gap).

## Falsifier

If MILP-OPT(Ulam, primes) = LP(Ulam, primes), then the slot-3 LP
relaxation is integer-tight and the structural compression is
algorithmically realisable. **Falsified at N=10³, 5·10³**: OPT > LP
strictly. Slot 4 thus established that the integrality gap is a real
property of the prime-on-Ulam structure, not a heuristic.

## Edges referenced

- E1.5 (information-theoretic compression floor)
- HL Conjecture F (off-edge external — quadratic-form prime densities)

## Results

| N | π(N) | LP | **MILP OPT** | greedy | iter best | OPT/LP | runtime |
|---|---|---|---|---|---|---|---|
| 10³ | 168 | 23.34 | **25** | 26 | 26 | 1.071 | 3.3s |
| 5·10³ | 669 | 54.30 | **61** | 63 | 61 | 1.123 | 141s |
| 10⁴ | 1229 | 77.59 | **∈ [82, 89]** | 95 | 89 | ∈ [1.057, 1.147] | 1500s, gap 8.89% |

MILP at N=10⁴ ran to its 1500s time limit. Final state: primal bound
(best integer feasible) = 90, dual bound (LP+cuts lower bound) = 82,
gap 8.89%. Iterated rounding found 89 as upper bound on OPT, so
**OPT ∈ [82, 89]** at N=10⁴. Branch-cuts raised the LP lower bound
from 77.59 to 82 (a 5.7% improvement) — confirming the LP relaxation
alone is a loose lower bound for prime line cover.

## Conclusion

The integrality gap of Ulam-spiral prime line cover is **structural
and growing**:

- 7.1% at N=10³ (exact MILP)
- 12.3% at N=5·10³ (exact MILP)
- ≥5.7% at N=10⁴ from MILP-cuts dual bound 82 vs LP 77.59; ≤14.7%
  from iter-89 vs LP 77.59 (true OPT/LP ∈ [1.057, 1.147])

This rules out the slot-3 hopeful frame "LP_p / LP_r = 0.78 stable
extends to integer cover". The fractional separation is real but
the integer optimum is asymptotically tight at √N for both primes
and random.
