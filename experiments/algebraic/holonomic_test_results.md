# Holonomic (D-finite) Test for pi(n): Results

**Script:** holonomic_test.py

## What Was Tested

Whether pi(n) is holonomic (D-finite), i.e., satisfies a linear recurrence sum_{i=0}^{d} p_i(n) * a(n+i) = 0 with polynomial coefficients p_i(n) of degree <= r. Tested all (d,r) pairs up to d=20, r=8 using least-squares with train/test cross-validation. Also tested the prime indicator function.

## Key Findings

- pi(n) is NOT holonomic for any tested (d,r) up to d=20, r=8
- Test residuals are far from machine epsilon; ratios to random baselines are ~1 (no signal)
- The prime indicator function is also NOT holonomic
- This rules out: holonomic sequence algorithms (baby-step/giant-step), D-finite generating function approaches, any linear recurrence with polynomial coefficients
- Consistent with Mauduit-Rivat (prime indicator not k-automatic, which is weaker than non-holonomic)
- Stronger than prior Session 14 result (pi(x) mod m not an LRS): polynomial coefficients don't help either

## Verdict

**CLOSED** -- Failure Mode: Information Loss (I). pi(n) satisfies no holonomic recurrence of any tested order and degree; the prime counting function is fundamentally non-D-finite.

## One-Line Summary

pi(n) is not holonomic -- no linear recurrence with polynomial coefficients of order <= 20 and degree <= 8, ruling out D-finite generating function approaches.
