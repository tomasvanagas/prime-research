# Proposal 4: Finite Field Lifting / Function Field Analogy -- Results

**Script:** proposal4_finite_field_lift.py

## What Was Tested
Exploit the exact closed-form formula for irreducible polynomial counts over F_q[x] (I_q(n) = (1/n) sum_{d|n} mu(n/d) * q^d) and take a q -> 1 limit ("field with one element" F_1) to try to extract pi(x) or p(n) for the integers.

## Key Findings
- I_q(n) for fixed q is exact and efficient: O(d(n)) divisor operations. No sqrt(x) barrier in function fields because RH is proved (Weil 1948).
- As q -> 1: I_q(n) -> (1/n) sum_{d|n} mu(n/d) = [n=1]/n. This degenerates to 0 for n > 1 -- not useful.
- Regularized limit I_q(n)/(q-1) as q -> 1 via L'Hopital gives sum_{d|n} d * mu(n/d) / n = phi(n)/n -- the totient ratio, not pi(x).
- The F_1 analogy (Connes-Consani, Borger, Deninger) is a deep open research program. No concrete algorithm emerges because the "arithmetic site" has infinite-dimensional cohomology (all zeta zeros), unlike F_q curves which have finite-dimensional H^1.
- The function field advantage is precisely that H^1 is finite-dimensional (2g); for Z, it is infinite-dimensional.

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (E) -- the q -> 1 limit degenerates; the function field advantage (finite H^1) has no analogue over Z.

## One-Line Summary
F_1 / function field analogy: q->1 limit of I_q(n) degenerates; infinite-dimensional cohomology of Spec(Z) blocks transfer.
