# Botkin-Dawsey-Hemmer-Just-Schneider Partition-Theoretic Prime Model: Results

**Script:** botkin_partition.py

## What Was Tested
Partition-theoretic model from arXiv:2501.00580v3 where p_n = 1 + 2*sum ceil(d(j)/2) + epsilon(n). Three models tested: epsilon=0 (Model 1), epsilon from semiprime counting (Model 2), and epsilon from Landau asymptotic (Model 2*). Accuracy measured up to n=10000.

## Key Findings
- The formula reduces to p_n = 1 + sum_{k=1}^{n-1} floor((n-1)/k) + floor(sqrt(n-1)), which requires O(n) divisor function evaluations
- Model 1 (epsilon=0) gives an approximation that diverges from true p(n) with growing error
- Model 2 (semiprime correction) improves accuracy but computing pi_2(x) is itself O(x^{2/3})
- The core computation sum_{k=1}^{N} floor(N/k) is the divisor summatory function, computable in O(sqrt(N)) but N = p(n) which is the unknown
- No polylog shortcut: the formula fundamentally requires enumerating up to p(n)

## Verdict
**CLOSED** -- Failure Mode: C (Circularity)

## One-Line Summary
Partition-theoretic model gives a valid formula but its evaluation cost is O(p(n)), worse than existing methods.
