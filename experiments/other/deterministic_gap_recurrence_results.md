# Deterministic Gap Recurrence: Results

**Script:** deterministic_gap_recurrence.py

## What Was Tested
Whether a deterministic nonlinear recurrence g(n) = F(g(n-1), ..., g(n-k), n) exists for prime gaps, enabling p(n) = 2 + sum g(k). Four approaches: Rowland's sequence, Cloitre's recurrence, the greedy prime construction, and Gandhi's formula.

## Key Findings
- Rowland's sequence a(n) = a(n-1) + gcd(n, a(n-1)) generates primes but is O(n^2) in total work and produces primes out of order with many repeats
- Cloitre's recurrence and similar analytic recurrences produce primes but at O(n) cost per prime, no better than sequential enumeration
- The greedy prime construction (smallest coprime to all previous) IS the primes by definition -- circular
- Gandhi's formula and recursive formulas all require O(n) work per prime at minimum
- No recurrence with polylog evaluation cost was found

## Verdict
**CLOSED** -- Failure Mode: C (Circularity)

## One-Line Summary
Known prime-generating recurrences (Rowland, Cloitre, Gandhi) all require at least O(n) work per prime.
