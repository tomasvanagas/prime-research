# Session 15 Summary

**Date:** 2026-04-04
**Focus:** Determinantal complexity, randomized algorithms, #TC^0 counting, monotone complexity, novel intermediate quantities, lower bounds, literature search

## Key Results

### 1. Determinantal Complexity Connection (NOVEL)
pi(x) viewed as degree-N multilinear polynomial in bits of x. Found N×N determinantal
representations for N=2,3,4 via numerical optimization. "Is pi(x) in GapL?" ⟺
"Does pi(x) have polynomial determinantal complexity?" For N≥10, generic polynomials
DON'T have N×N det representations — pi(x) would need special structure from number theory.
Documented in novel/determinantal_complexity.md.

### 2. #TC^0 ⊆ NC? is THE Complexity Question
If BPSW ∈ TC^0 (conditional on BPSW correctness), then pi(x) ∈ NC iff #TC^0 ⊆ NC.
The Fermat residue coupling (2^{n-1} mod n where exponent and modulus both depend on n)
prevents batch counting. The #TC^0 counting problem is NOT known to be in NC.
Pure complexity theory question, independent of number theory.

### 3. Uniformity is the True Barrier (NOVEL)
Nonuniform poly(N)-size circuits for pi(x) exist trivially (hardcode primes as advice).
The hard part is UNIFORM construction — generating primes without knowing them.
Natural proofs barrier blocks proving super-TC^0 lower bounds.
Documented in novel/uniformity_barrier.md.

### 4. ALL Randomized Approaches Fail
- Zeta zero sampling: must use 100% of zeros (variance prohibitive)
- Probabilistic sieve: 10^6-10^9x WORSE than exhaustive (I-E cancellation)
- Hash-based counting: exact sketch needs O(x^2) bits
- Quantum counting: gap between sqrt(ln x) lower bound and sqrt(x) upper bound

### 5. Divide-and-Conquer Fails
Error at each level O(sqrt(x_level)), accumulates to O(sqrt(x)) total.
Buchstab tree has O(x^{2/3+eps}) nodes. No recursion gives sublinear error.

### 6. Monotone Complexity Inapplicable
[pi(x) >= k] = [x >= p(k)] — just integer comparison with a constant.
Monotone complexity: trivially O(N). Individual bits of pi(x) NOT monotone.
Razborov-style lower bounds don't apply.

### 7. Arithmetic Circuit Complexity Doesn't Help
VP = VNC^2 (VSBR): depth is free in algebraic model. But pi(x) is NOT a polynomial
over fields, making arithmetic complexity theory hard to apply directly.
Tau conjecture orthogonal. No new viable path.

### 8. Systematic Analysis of 8 Intermediate Quantity Families
All 8 candidates (residues, polynomial evals, matrix eigenvalues, topology,
representation theory, entropy, recursive/self-similar, physical/dynamical) either:
- Route through floor values or zeta zeros (Equivalence)
- Carry insufficient information (Information Loss)
- Require solving an open complexity theory problem

### 9. Literature Search (March-April 2026)
- **No new algorithms for pi(x)** — barriers intact
- **Chen-Tal-Wang STOC 2026**: Superquadratic depth-2 threshold lower bounds (near TC^0)
- **ECCC TR26-035**: Learning Read-Once Determinants — relevant to GapL
- **Raz ECCC TR26-008**: Natural proofs barrier extends to linear functions
- **arXiv:2507.23563**: Complexity of log-space bounded counting classes

## Experiments Created
- experiments/circuit_complexity/divide_and_conquer_pi.py
- experiments/circuit_complexity/det_perm_encoding.py
- experiments/circuit_complexity/randomized_pi.py
- experiments/circuit_complexity/det_complexity_search.py
- experiments/circuit_complexity/det_complexity_extract.py
- experiments/circuit_complexity/tc0_counting.py
- experiments/circuit_complexity/monotone_complexity.py

## Novel Findings
- novel/determinantal_complexity.md — dc(pi_N) = N for N=2,3,4; GapL ⟺ poly dc
- novel/uniformity_barrier.md — Uniformity, not monotonicity, is the true barrier

## Updated Files
- CLAUDE.md (Session 15 status, closed paths, viable directions)
- status/CLOSED_PATHS.md (~12 new entries, count 445+)
- status/OPEN_PROBLEMS.md (Session 15 refinements)

## Paths Still Open
The problem remains GENUINELY OPEN. Key remaining directions:
1. **Determinantal complexity**: Does pi_N have poly(N) dc? Needs algebraic geometry
2. **#TC^0 ⊆ NC?**: Pure complexity theory, independent of primes
3. **BPSW unconditional proof**: Would give PRIMES in TC^0
4. **Kt complexity**: Time-bounded Kolmogorov complexity of delta(n)
5. **Berry-Keating Hamiltonian**: Quantum, far from concrete
6. **Novel identity**: Serendipity only — all systematic searches exhausted
