# Ramanujan Library-style PSLQ on delta(n) -- Results

**Script:** ramanujan_pslq_delta.py

## What Was Tested
Search for algebraic/recurrence/modular integer relations among delta(n) = p(n) - R^{-1}(n) values using PSLQ algorithm (inspired by ICLR 2025 Ramanujan Library paper). Tests: (a) linear recurrences up to order 20; (b) polynomial relations of degree 2-3 among consecutive deltas; (c) mixed relations with simple functions of n; (d) modular recurrences.

## Key Findings
- Linear recurrences: PSLQ finds no exact integer relations among delta(n), delta(n-1), ..., delta(n-k) for k up to 20. Residuals are O(1), not improving with higher order.
- Polynomial relations: no degree-2 or degree-3 polynomial relations among consecutive delta values. PSLQ returns null or relations with huge coefficients (numerical artifacts).
- Mixed relations (delta vs n, log(n), sqrt(n), etc.): no relations found. delta(n) is algebraically independent of simple functions of n.
- Modular recurrences: delta(n) mod m for small m (2,3,5,7) shows no recurrence structure. Sequence appears uniformly distributed mod m.
- Training on n=1..5000, validation on n=5001..6000: any patterns found in training fail completely on validation -- consistent with overfitting noise.
- This confirms delta(n) has no low-complexity algebraic structure detectable by PSLQ.

## Verdict
**CLOSED**
**Failure Mode:** Information Loss (I) -- delta(n) has no detectable algebraic/recurrence structure; PSLQ finds no relations, consistent with GUE-random behavior.

## One-Line Summary
PSLQ on delta(n): no linear recurrences, polynomial relations, or modular structure found; delta is algebraically unstructured.
