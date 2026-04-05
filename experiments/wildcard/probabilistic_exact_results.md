# Probabilistic Exact: Results

**Script:** probabilistic_exact.py

## What Was Tested
Five probabilistic/randomized approaches: (1) randomized rounding of smooth approximation with CRT on random moduli, (2) Miller-Rabin style binary search, (3) sieve + verification with randomized sieves, (4) algebraic number field counting with "nicer" Dedekind zeta zeros, (5) trace formula shortcut via Selberg.

## Key Findings
- Randomized rounding + CRT: computing delta(n) mod M for random M still requires pi(x) information -- circular
- Binary search with local sieve: works but step 3 (counting pi(x_approx - Delta)) requires O(x^{2/3}) -- the hard problem returns
- Randomized sieve: concentration bounds narrow the window but still need exact count in the window, which costs O(window * sieve_cost)
- Number field counting: Dedekind zeta zeros for quadratic fields are NOT simpler than Riemann zeta zeros; they add Dirichlet L-function zeros
- Selberg trace formula: deep connection between geodesic lengths and eigenvalues, but computing the nth eigenvalue IS computing a zeta zero

## Verdict
**CLOSED**
**Failure Mode:** Circularity (C) for approaches 1,2,3; Equivalence (E) for approaches 4,5 -- all roads lead back to pi(x) or zeta zeros.

## One-Line Summary
Probabilistic exact approaches (random CRT, binary search, random sieve, number fields, trace formula): all circular or equivalent to known methods.
