# Proposal 20: Étale Cohomology Point-Counting Shortcut — Results

## Idea
Use algebraic geometry point-counting (Grothendieck-Lefschetz trace formula)
to compute pi(x) by encoding prime counting as counting F_p-points on varieties.

## Results

### Character sum approach
- Correctly computes pi(x, q, a) = #{p <= x : p ≡ a mod q} for all tested
  cases (q = 3, 5, 7, 11; x = 100, 1000)
- BUT: the computation iterates over primes up to x, same cost as direct counting
- Character sums just REORGANIZE the computation via Fourier analysis on (Z/qZ)*

### Hasse-Weil constraint analysis
| n | p(n) | Bits from 119 constraints | Bits needed |
|---|------|--------------------------|-------------|
| 100 | 541 | 480.4 | 9.1 |
| 1000 | 7919 | 759.0 | 13.0 |
| 5000 | 48611 | 807.3 | 15.6 |

The constraints provide MORE than enough bits of information.
**Catch**: computing each constraint requires pi(x, q, a), which is O(x^{2/3}).

### Elliptic curve traces (y^2 = x^3 - x)
- 167 a_p values computed, range [-62, 62], mean 0.168 (near 0 as expected)
- Cumulative sum follows sqrt(pi(x)) (Sato-Tate)
- Correlation with pi(x) residual: 0.39 — weak, not useful for prediction
- EC point counts encode CURVE information, not prime distribution information

## Why algebraic geometry doesn't shortcut prime counting
1. Point-counting on a FIXED variety is O(polylog) — but that counts points on
   the variety, not primes among integers
2. "Is k prime?" requires checking divisibility by ALL primes up to sqrt(k),
   which cannot be encoded as a fixed-dimension variety
3. The number of varieties needed to sieve grows with x
4. Character sums (L-functions) DO encode prime distribution but evaluating
   them requires O(x^{1/2+eps}) zeros — the same barrier as the explicit formula
5. The Hasse-Weil bound gives info about each constraint, but the constraints
   themselves are as expensive to compute as pi(x)

## Verdict: CLOSED
The gap is fundamental: algebraic geometry counts points on FIXED algebraic objects
in polylog time, but the "primality predicate" over [1,x] has unbounded algebraic
complexity as x grows. No fixed-dimensional variety encodes pi(x).
