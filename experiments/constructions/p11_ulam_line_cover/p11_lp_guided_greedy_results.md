# Thread 11 / Slot 4 — LP-guided greedy variant (negative result)

## Idea

After solving the LP relaxation, restrict greedy to only LP-active
lines (x_l > 0). Hypothesis: LP support concentrates on structurally
important lines, so restricted greedy might match or beat full greedy.

## Result (Ulam, N=10⁴, K=5)

- bounded greedy (full): L = 95
- LP-guided greedy (restricted to 462 active LP lines): L = 97 (91
  multi-line + 6 singleton)

**LP-guided greedy is WORSE than pure greedy.** The LP support is
too narrow — the greedy needs the freedom to pick lines outside the
LP active set to cover certain primes.

## Conclusion

LP-guided greedy is not a useful heuristic for this problem. Iterated
rounding (slot 4) and pure greedy remain the better upper bounds;
true MILP optimum is between.

## Files

- `p11_lp_guided_greedy.py` (the script)

This is filed for completeness (a heuristic tested, found inferior).
