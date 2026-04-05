# Nonlinear Bitwise / TC0: Results

**Date:** 2026-04-04 (Session 14)
**Script:** nonlinear_bitwise_tc0.py

## What Was Tested
Bit decomposition of floor values, polynomial degree vs floor values needed, and mutual information analysis of floor value pairs.

## Key Findings
- Bit features + AND pairs: max_err=1.53, 428/498 exact, but ~15% error rate with no generalization guarantee
- Increasing K or polynomial degree improves training, worsens testing (structural overfit)
- I(floor(x/2); pi(x)) = 7.65 bits (95% of info in x/2 alone); pairs barely add info beyond individuals
- Critical insight: diversity of information requires k up to O(sqrt(x))

## Verdict
**CLOSED** -- see nonlinear_sieve_summary.md for detailed results.
**Failure Mode:** E (Equivalence)

## One-Line Summary
Bitwise/TC0 operations on floor values overfit training data; mutual information confirms O(sqrt(x)) floor values needed.
